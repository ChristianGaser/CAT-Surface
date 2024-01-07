#include "bicpl_internal.h"

static  void  coalesce_lines(
    lines_struct   *lines )
{
    int   p, n, l, *n_neighbours, **neighbours, total_neighbours;
    int   size, p1, p2, n_items, n_indices, *indices, *end_indices, edge;

    if( lines->n_points <= 0 ||
        lines->n_items <= 0 )
    {
        return;
    }

    ALLOC( n_neighbours, lines->n_points  );

    for_less( p, 0, lines->n_points )
        n_neighbours[p] = 0;

    for_less( l, 0, lines->n_items )
    {
        size = GET_OBJECT_SIZE( *lines, l );

        for_less( edge, 0, size-1 )
        {
            p1 = lines->indices[POINT_INDEX(lines->end_indices,l,edge)];
            p2 = lines->indices[POINT_INDEX(lines->end_indices,l,edge+1)];

            ++n_neighbours[p1];
            ++n_neighbours[p2];
        }
    }

    total_neighbours = 0;
    for_less( p, 0, lines->n_points )
    {
        if( n_neighbours[p] != 0 && n_neighbours[p] != 2 )
        {
            return;
        }
        total_neighbours += n_neighbours[p];
    }

    ALLOC( neighbours, lines->n_points );
    ALLOC( neighbours[0], total_neighbours );
    for_less( p, 1, lines->n_points )
    {
        neighbours[p] = &neighbours[p-1][n_neighbours[p-1]];
    }

    for_less( p, 0, lines->n_points )
        n_neighbours[p] = 0;

    for_less( l, 0, lines->n_items )
    {
        size = GET_OBJECT_SIZE( *lines, l );

        for_less( edge, 0, size-1 )
        {
            p1 = lines->indices[POINT_INDEX(lines->end_indices,l,edge)];
            p2 = lines->indices[POINT_INDEX(lines->end_indices,l,edge+1)];

            neighbours[p1][n_neighbours[p1]] = p2;
            ++n_neighbours[p1];
            neighbours[p2][n_neighbours[p2]] = p1;
            ++n_neighbours[p2];
        }
    }

    n_items = 0;
    end_indices = NULL;
    n_indices = 0;
    indices = NULL;

    for_less( l, 0, lines->n_items )
    {
        size = GET_OBJECT_SIZE( *lines, l );

        for_less( edge, 0, size-1 )
        {
            p1 = lines->indices[POINT_INDEX(lines->end_indices,l,edge)];
            p2 = lines->indices[POINT_INDEX(lines->end_indices,l,edge+1)];

            for_less( n, 0, n_neighbours[p1] )
            {
                if( neighbours[p1][n] == p2 )
                    break;
            }

            if( n >= n_neighbours[p1] )
                continue;

            ADD_ELEMENT_TO_ARRAY( indices, n_indices, p1, DEFAULT_CHUNK_SIZE );

            while( n < n_neighbours[p1] )
            {
                ADD_ELEMENT_TO_ARRAY( indices, n_indices, p2,
                                      DEFAULT_CHUNK_SIZE );

                neighbours[p1][n] = -1;

                for_less( n, 0, n_neighbours[p2] )
                {
                    if( neighbours[p2][n] == p1 )
                        break;
                }

                if( n >= n_neighbours[p1] )
                {
                    handle_internal_error( "coalesce_lines" );
                    return;
                }

                neighbours[p2][n] = -1;

                p1 = p2;

                for_less( n, 0, n_neighbours[p1] )
                {
                    if( neighbours[p1][n] >= 0 )
                        break;
                }

                p2 = neighbours[p1][n];
            }

            ADD_ELEMENT_TO_ARRAY( end_indices, n_items, n_indices,
                                  DEFAULT_CHUNK_SIZE );
        }
    }

    FREE( lines->end_indices );
    FREE( lines->indices );

    lines->n_items = n_items;
    lines->indices = indices;
    lines->end_indices = end_indices;
}

BICAPI  void   intersect_planes_with_polygons(
    polygons_struct   *polygons,
    Point             *plane_origin,
    Vector            *plane_normal,
    lines_struct      *lines )
{
    int                 n_points, p, index;
    int                 poly, edge, size;
    int                 point_index1, point_index2;
    Vector              v1, v2;
    Point               point;
    Real                t1, t2, ratios[2];
    int                 p1s[2], p2s[2];
    hash2_table_struct  hash;

    initialize_lines( lines, WHITE );

    initialize_hash2_table( &hash, polygons->n_items, (int) sizeof(int),
                            0.5, .25 );

    index = 0;

    for_less( poly, 0, polygons->n_items )
    {
        n_points = 0;

        size = GET_OBJECT_SIZE( *polygons, poly );

        for_less( edge, 0, size )
        {
            point_index1 = polygons->indices[
                     POINT_INDEX(polygons->end_indices,poly,edge)];
            point_index2 = polygons->indices[
                     POINT_INDEX(polygons->end_indices,poly,(edge+1)%size)];

            SUB_POINTS( v1, polygons->points[point_index1], *plane_origin );
            SUB_POINTS( v2, polygons->points[point_index2], *plane_origin );

            t1 = DOT_VECTORS( v1, *plane_normal );
            t2 = DOT_VECTORS( v2, *plane_normal );

            /*--- check if it intersects the line segment [p1 .. p2)  */

            if( t1 == 0.0 || (t1 > 0.0 && t2 < 0.0) || (t1 < 0.0 && t2 > 0.0) )
            {
                if( n_points < 2 )
                {
                    if( t1 == 0.0 )
                        ratios[n_points] = 0.0;
                    else
                        ratios[n_points] = t1 / (t1 - t2);

                    p1s[n_points] = MIN( point_index1, point_index2 );
                    p2s[n_points] = MAX( point_index1, point_index2 );
                    if( p1s[n_points] != point_index1 )
                        ratios[n_points] = 1.0 - ratios[n_points];

                    if( ratios[n_points] == 0.0 )
                    {
                        p2s[n_points] = n_points;
                        ratios[n_points] = 0.0;
                    }
                    else if( ratios[n_points] == 1.0 )
                    {
                        p1s[n_points] = p2s[n_points];
                        p2s[n_points] = n_points;
                        ratios[n_points] = 0.0;
                    }
                }

                ++n_points;
            }
        }

        if( n_points == 2 &&
            (ratios[0] != 0.0 || ratios[1] != 0.0 || p1s[0] != p1s[1]) )
        {
            start_new_line( lines );

            for_less( p, 0, 2 )
            {
                if( !lookup_in_hash2_table( &hash, p1s[p], p2s[p],
                                            (void *) &index ) )
                {
                    if( ratios[p] == 0.0 )
                        point = polygons->points[p1s[p]];
                    else
                    {
                        INTERPOLATE_POINTS( point,
                                            polygons->points[p1s[p]],
                                            polygons->points[p2s[p]],
                                            ratios[p] );
                    }

                    index = lines->n_points;
                    insert_in_hash2_table( &hash, p1s[p], p2s[p],
                                           (void *) &index );

                    ADD_ELEMENT_TO_ARRAY( lines->points, lines->n_points,
                                          point, DEFAULT_CHUNK_SIZE );
                }

                ADD_ELEMENT_TO_ARRAY( lines->indices,
                                      lines->end_indices[lines->n_items-1],
                                      index, DEFAULT_CHUNK_SIZE );
            }
        }
    }

    delete_hash2_table( &hash );

    coalesce_lines( lines );
}

BICAPI  void   intersect_planes_with_quadmesh(
    quadmesh_struct   *quadmesh,
    Point             *plane_origin,
    Vector            *plane_normal,
    lines_struct      *lines )
{
    int                 n_points, p, index, m, n, x_index, y_index;
    int                 edge;
    int                 point_index1, point_index2, indices[4];
    Vector              v1, v2;
    Point               point, *points;
    Real                t1, t2, ratios[2];
    int                 p1s[2], p2s[2];
    hash2_table_struct  hash;

    initialize_lines( lines, WHITE );

    get_quadmesh_n_objects( quadmesh, &m, &n );

    initialize_hash2_table( &hash, m * n, (int) sizeof(int),
                            0.5, .25 );

    points = quadmesh->points;

    index = 0;

    for_less( x_index, 0, m )
    {
        for_less( y_index, 0, n )
        {
            n_points = 0;

            indices[0] = IJ(x_index,y_index,quadmesh->n);
            indices[1] = IJ((x_index+1)%quadmesh->m,y_index,quadmesh->n);
            indices[2] = IJ((x_index+1)%quadmesh->m,(y_index+1)%quadmesh->n,
                         quadmesh->n);
            indices[3] = IJ(x_index,(y_index+1)%quadmesh->n,
                         quadmesh->n);

            for_less( edge, 0, 4 )
            {
                point_index1 = indices[edge];
                point_index2 = indices[(edge+1) % 4];

                SUB_POINTS( v1, points[point_index1], *plane_origin );
                SUB_POINTS( v2, points[point_index2], *plane_origin );

                t1 = DOT_VECTORS( v1, *plane_normal );
                t2 = DOT_VECTORS( v2, *plane_normal );

                /*--- check if it intersects the line segment [p1 .. p2)  */

                if( t1 == 0.0 || (t1 > 0.0 && t2 < 0.0) || (t1 < 0.0 && t2 > 0.0) )
                {
                    if( n_points < 2 )
                    {
                        if( t1 == 0.0 )
                            ratios[n_points] = 0.0;
                        else
                            ratios[n_points] = t1 / (t1 - t2);

                        p1s[n_points] = MIN( point_index1, point_index2 );
                        p2s[n_points] = MAX( point_index1, point_index2 );
                        if( p1s[n_points] != point_index1 )
                            ratios[n_points] = 1.0 - ratios[n_points];

                        if( ratios[n_points] == 0.0 )
                        {
                            p2s[n_points] = n_points;
                            ratios[n_points] = 0.0;
                        }
                        else if( ratios[n_points] == 1.0 )
                        {
                            p1s[n_points] = p2s[n_points];
                            p2s[n_points] = n_points;
                            ratios[n_points] = 0.0;
                        }
                    }

                    ++n_points;
                }
            }

            if( n_points == 2 &&
                (ratios[0] != 0.0 || ratios[1] != 0.0 || p1s[0] != p1s[1]) )
            {
                start_new_line( lines );

                for_less( p, 0, 2 )
                {
                    if( !lookup_in_hash2_table( &hash, p1s[p], p2s[p],
                                                (void *) &index ) )
                    {
                        if( ratios[p] == 0.0 )
                            point = points[p1s[p]];
                        else
                        {
                            INTERPOLATE_POINTS( point,
                                                points[p1s[p]], points[p2s[p]],
                                                ratios[p] );
                        }

                        index = lines->n_points;
                        insert_in_hash2_table( &hash, p1s[p], p2s[p],
                                               (void *) &index );

                        ADD_ELEMENT_TO_ARRAY( lines->points, lines->n_points,
                                              point, DEFAULT_CHUNK_SIZE );
                    }

                    ADD_ELEMENT_TO_ARRAY( lines->indices,
                                          lines->end_indices[lines->n_items-1],
                                          index, DEFAULT_CHUNK_SIZE );
                }
            }
        }
    }

    delete_hash2_table( &hash );

    coalesce_lines( lines );
}
