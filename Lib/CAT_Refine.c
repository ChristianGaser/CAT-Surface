/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 *
 */

#include <bicpl.h>

#include "CAT_Curvature.h"
#include "CAT_Smooth.h"

private  BOOLEAN  lookup_edge_midpoint(
    hash2_table_struct    *edge_lookup,
    int                   p0,
    int                   p1,
    int                   *midpoint )
{
    int     k0, k1;

    k0 = MIN( p0, p1 );
    k1 = MAX( p0, p1 );

    return( lookup_in_hash2_table( edge_lookup, k0, k1, (void *) midpoint ) );
}

private  void  subdivide_edge(
    hash2_table_struct *edge_lookup,
    double               normalized_length,
    int                p1,
    int                p2,
    polygons_struct    *new_polygons,
    Point              *length_points[] )
{
    int    midpoint;
    Point  mid;

    midpoint = new_polygons->n_points;
    INTERPOLATE_POINTS( mid, new_polygons->points[p1],
                             new_polygons->points[p2], 0.5 );
    ADD_ELEMENT_TO_ARRAY( new_polygons->points, new_polygons->n_points,
                          mid, DEFAULT_CHUNK_SIZE );

    INTERPOLATE_POINTS( mid, (*length_points)[p1], (*length_points)[p2], 0.5 );
    --new_polygons->n_points;
    ADD_ELEMENT_TO_ARRAY( *length_points, new_polygons->n_points,
                          mid, DEFAULT_CHUNK_SIZE );

    insert_in_hash2_table( edge_lookup, MIN(p1,p2), MAX(p1,p2),
                           (void *) &midpoint );
}

private  void  add_polygons(
    int                p0,
    int                p1,
    int                p2,
    Point              length_points[],
    hash2_table_struct *edge_lookup,
    polygons_struct    *new_polygons,
    int                *n_indices )
{
    int      tri, edge, mid[3], n_present, indices[3], i0, i1, i2;
    int      n_new_triangles, new_indices[4][3], offset;
    BOOLEAN  mid_present[3];

    indices[0] = p0;
    indices[1] = p1;
    indices[2] = p2;

    n_present = 0;
    for_less( edge, 0, 3 )
    {
        mid_present[edge] = lookup_edge_midpoint( edge_lookup, indices[edge],
                                             indices[(edge+1)%3], &mid[edge] );

        if( mid_present[edge] )
            ++n_present;
    }

    if( n_present == 3 )
    {
        n_new_triangles = 4;
        new_indices[0][0] = p0;
        new_indices[0][1] = mid[0];
        new_indices[0][2] = mid[2];
        new_indices[1][0] = mid[0];
        new_indices[1][1] = p1;
        new_indices[1][2] = mid[1];
        new_indices[2][0] = mid[0];
        new_indices[2][1] = mid[1];
        new_indices[2][2] = mid[2];
        new_indices[3][0] = mid[2];
        new_indices[3][1] = mid[1];
        new_indices[3][2] = p2;
    }
    else if( n_present == 2 )
    {
        offset = 0;
        while( mid_present[offset] )
            ++offset;

        n_new_triangles = 3;

        i0 = offset;
        i1 = (offset+1) % 3;
        i2 = (offset+2) % 3;

        if( sq_distance_between_points( &length_points[indices[i0]],
                                        &length_points[mid[i1]] ) <
            sq_distance_between_points( &length_points[indices[i1]],
                                        &length_points[mid[i2]] ) )
        {
            new_indices[0][0] = indices[i0];
            new_indices[0][1] = indices[i1];
            new_indices[0][2] = mid[i1];
            new_indices[1][0] = indices[i0];
            new_indices[1][1] = mid[i1];
            new_indices[1][2] = mid[i2];
        }
        else
        {
            new_indices[0][0] = indices[i0];
            new_indices[0][1] = indices[i1];
            new_indices[0][2] = mid[i2];
            new_indices[1][0] = indices[i1];
            new_indices[1][1] = mid[i1];
            new_indices[1][2] = mid[i2];
        }
        new_indices[2][0] = mid[i2];
        new_indices[2][1] = mid[i1];
        new_indices[2][2] = indices[i2];
    }
    else if( n_present == 1 )
    {
        offset = 0;
        while( !mid_present[offset] )
            ++offset;

        n_new_triangles = 2;

        i0 = offset;
        i1 = (offset+1) % 3;
        i2 = (offset+2) % 3;

        new_indices[0][0] = indices[i0];
        new_indices[0][1] = mid[i0];
        new_indices[0][2] = indices[i2];
        new_indices[1][0] = mid[i0];
        new_indices[1][1] = indices[i1];
        new_indices[1][2] = indices[i2];
    }
    else if( n_present == 0 )
    {
        n_new_triangles = 0;

        ADD_ELEMENT_TO_ARRAY( new_polygons->indices, *n_indices,
                              p0, DEFAULT_CHUNK_SIZE );
        ADD_ELEMENT_TO_ARRAY( new_polygons->indices, *n_indices,
                              p1, DEFAULT_CHUNK_SIZE );
        ADD_ELEMENT_TO_ARRAY( new_polygons->indices, *n_indices,
                              p2, DEFAULT_CHUNK_SIZE );
        ADD_ELEMENT_TO_ARRAY( new_polygons->end_indices, new_polygons->n_items,
                              *n_indices, DEFAULT_CHUNK_SIZE );
    }

    for_less( tri, 0, n_new_triangles )
    {
        add_polygons( new_indices[tri][0],
                      new_indices[tri][1],
                      new_indices[tri][2],
                      length_points, edge_lookup, new_polygons, n_indices );
    }
}

int refine_mesh(
    Point              *length_points[],
    polygons_struct    *polygons,
    double             max_length,
    polygons_struct    *new_polygons,
    double             weight_curvature )
{
    int                  n_indices, i, p1, p2, size, point, edge, midpoint, poly;
    int                  *n_neighbours, **neighbours;
    double               normalized_length, *curvatures, max_length_weighted;
    double               mean_curv, max_curv, sum_curv;
    hash2_table_struct   edge_lookup;

    /* estimate (scaled) absolute mean curvature if weighting is defined */
    if (weight_curvature > 0.0) {
            ALLOC(curvatures, polygons->n_points);
            get_all_polygon_point_neighbours(polygons, &n_neighbours, &neighbours);
            get_polygon_vertex_curvatures_cg(polygons, n_neighbours, neighbours,
                                         3.0, 0, curvatures);
                                         
            /* get absolute value, transform data to more normally distributed data 
               with x^0.4 and finally get maximum */
            max_curv = 0.0; 
            for (i = 0; i < polygons->n_points; i++) {
                    curvatures[i] = pow(fabs(curvatures[i]), 0.4);
                    max_curv = MAX(curvatures[i], max_curv);
            }
            
            /* scale curvature to range (0..1) */
            sum_curv = 0.0;
            for (i = 0; i < polygons->n_points; i++) {
                    curvatures[i] = curvatures[i]/max_curv;
                    sum_curv += curvatures[i];
            }
                    
            /* force mean of 0 and add 1 to ensure that mean curvature vales are weighted with 1 */
            mean_curv = sum_curv/(double)polygons->n_points;
            for (i = 0; i < polygons->n_points; i++) 
                    curvatures[i] = curvatures[i] - mean_curv + 1.0;
    }
    
    initialize_polygons( new_polygons, WHITE, NULL );

    SET_ARRAY_SIZE( new_polygons->points, 0, polygons->n_points,
                    DEFAULT_CHUNK_SIZE );

    new_polygons->n_points = polygons->n_points;

    for_less( point, 0, polygons->n_points )
        new_polygons->points[point] = polygons->points[point];

    initialize_hash2_table( &edge_lookup, 10 * polygons->n_points,
                            sizeof(int), 0.25, 0.125 );

    for_less( poly, 0, polygons->n_items )
    {
        size = GET_OBJECT_SIZE( *polygons, poly );

        for_less( edge, 0, size )
        {
            p1 = polygons->indices[POINT_INDEX(polygons->end_indices,poly,edge)];
            p2 = polygons->indices[POINT_INDEX(polygons->end_indices,poly,
                                 (edge+1)%size)];
            
            if (weight_curvature > 0.0) {
                    max_length_weighted = max_length/pow(curvatures[p1],weight_curvature);
            } else  max_length_weighted = max_length;
            
            normalized_length = distance_between_points( &(*length_points)[p1],
                                &(*length_points)[p2])/max_length_weighted;

            if( normalized_length > 1.0 &&
                !lookup_edge_midpoint( &edge_lookup, p1, p2, &midpoint ) )
            {
                subdivide_edge( &edge_lookup, normalized_length,
                                p1, p2, new_polygons, length_points );
            }
        }
    }

    n_indices = 0;
    for_less( poly, 0, polygons->n_items )
    {
        add_polygons( polygons->indices[
                               POINT_INDEX(polygons->end_indices,poly,0)],
                      polygons->indices[
                               POINT_INDEX(polygons->end_indices,poly,1)],
                      polygons->indices[
                               POINT_INDEX(polygons->end_indices,poly,2)],
                      *length_points, &edge_lookup, new_polygons, &n_indices );
    }

    ALLOC( new_polygons->normals, new_polygons->n_points );

    compute_polygon_normals( new_polygons );

    delete_hash2_table( &edge_lookup );
    if (weight_curvature > 0.0) FREE(curvatures);

    return( new_polygons->n_items - polygons->n_items );
}

