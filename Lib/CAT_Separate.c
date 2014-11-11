/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 *
 */

#include  "CAT_Separate.h"

private  BOOLEAN  recursive_triangulate_one_polygon(
    int     size,
    int     poly[],
    int     n_neighbours[],
    int     *neighbours[],
    int     indices[] )
{
    int       *left, *right, n_left, n_right, p1, p2;
    int       start_index, end_index, count, i, n;
    BOOLEAN   found;

    if( size < 3 )
        handle_internal_error( "recursive_triangulate_one_polygon" );
    
    if( size == 3 )
    {
        indices[0] = poly[0];
        indices[1] = poly[1];
        indices[2] = poly[2];
        return( TRUE );
    }

    found = FALSE;

    ALLOC( left, size );
    ALLOC( right, size );

    for_less( start_index, 0, size-2 )
    {
        p1 = poly[start_index];

        for_less( end_index, start_index+2, size )
        {
            if( start_index == 0 && end_index == size-1 )
                continue;

            p2 = poly[end_index];

            count = 0;
            for_less( n, 0, n_neighbours[p1] )
            {
                if( neighbours[p1][n] == p2 )
                    ++count;
            }

            if( count != 1 )
                continue;

            n_left = 0;
            for_inclusive( i, 0, start_index )
                left[n_left++] = poly[i];
            for_less( i, end_index, size )
                left[n_left++] = poly[i];

            n_right = 0;
            for_inclusive( i, start_index, end_index )
                right[n_right++] = poly[i];

            if( n_left + n_right != size + 2 )
            {
                handle_internal_error( "n_left" );
            }

            if( recursive_triangulate_one_polygon( n_left, left,
                                       n_neighbours, neighbours, indices ) &&
                recursive_triangulate_one_polygon( n_right, right,
                                       n_neighbours, neighbours,
                                       &indices[3*(n_left-2)] ) )
            {
                found = TRUE;
                break;
            }
        }

        if( found )
            break;
    }

    FREE( left );
    FREE( right );

    return( found );
}

private  BOOLEAN  triangulate_one_polygon(
    int     size,
    int     poly[],
    int     n_neighbours[],
    int     *neighbours[],
    int     indices[] )
{
    return( recursive_triangulate_one_polygon( size, poly,
                                n_neighbours, neighbours, indices ) );
}

void  triangulate_polygons(
    polygons_struct  *polygons,
    polygons_struct  *triangles )
{
    int                poly, size, index, ind, n_matches, n;
    int                *n_neighbours, **neighbours, *indices, max_size;
    progress_struct    progress;
    BOOLEAN            done;

    create_polygon_point_neighbours( polygons, TRUE,
                                     &n_neighbours, &neighbours,
                                     NULL, NULL );

    *triangles = *polygons;

    triangles->colour_flag = ONE_COLOUR;
    ALLOC( triangles->colours, 1 );
    triangles->colours[0] = polygons->colours[0];

    triangles->points = polygons->points;
    triangles->normals = polygons->normals;

    triangles->n_items = 0;
    for_less( poly, 0, polygons->n_items )
        triangles->n_items += GET_OBJECT_SIZE( *polygons, poly ) - 2;

    ALLOC( triangles->indices, 3 * triangles->n_items );

    max_size = 0;
    for_less( poly, 0, polygons->n_items )
        max_size = MAX( max_size, GET_OBJECT_SIZE( *polygons, poly ) );

    ALLOC( indices, max_size );

    initialize_progress_report( &progress, FALSE, polygons->n_items,
                                "Triangulating" );

    ind = 0;
    for_less( poly, 0, polygons->n_items )
    {
        size = GET_OBJECT_SIZE( *polygons, poly );
        for_less( index, 0, size )
        {
            indices[index] = polygons->indices[POINT_INDEX(
                                            polygons->end_indices,poly,index)];
        }

        for_less( index, 2, size-1 )
        {
            n_matches = 0;
            for_less( n, 0, n_neighbours[indices[0]] )
            {
                if( neighbours[indices[0]][n] == indices[index] )
                    ++n_matches;
            }

            if( n_matches != 1 )
                break;
        }

        if( size > 3 && index < size-1 )
        {
            done = triangulate_one_polygon( size, indices,
                                              n_neighbours, neighbours,
                                              &triangles->indices[ind] );
            if( !done )
                print( "Could not find good triangulation: %d\n", poly );
            else
                ind += 3 * (size-2);
        }
        else
            done = FALSE;

        if( !done )
        {
            for_less( index, 1, size-1 )
            {
                triangles->indices[ind] = indices[0];
                ++ind;
                triangles->indices[ind] = indices[index];
                ++ind;
                triangles->indices[ind] = indices[index+1];
                ++ind;
            }
        }

        update_progress_report( &progress, poly+1 );
    }

    terminate_progress_report( &progress );

    if( ind != 3 * triangles->n_items )
        handle_internal_error( "Summation of ind" );

    FREE( polygons->end_indices );
    FREE( polygons->indices );
    FREE( indices );

    delete_polygon_point_neighbours( polygons, n_neighbours, neighbours,
                                     NULL, NULL );

    ALLOC( triangles->end_indices, triangles->n_items );
    for_less( poly, 0, triangles->n_items )
        triangles->end_indices[poly] = 3 * (poly + 1);
}

private  int   make_connected_components(
    polygons_struct    *polygons,
    Smallest_int       polygon_classes[],
    int                n_in_class[] )
{
    int                poly, current_poly, edge, size;
    int                neigh;
    int                n_components;
    Smallest_int       not_done;
    QUEUE_STRUCT(int)  queue;

    n_components = 0;

    not_done = (Smallest_int) 999;

    for_less( poly, 0, polygons->n_items )
        polygon_classes[poly] = not_done;

    for_less( poly, 0, polygons->n_items )
    {
        if( polygon_classes[poly] != not_done )
            continue;

        if( n_components == 999 )
        {
            ++n_components;
            break;
        }

        INITIALIZE_QUEUE( queue );
        INSERT_IN_QUEUE( queue, poly );
        polygon_classes[poly] = (Smallest_int) n_components;
        n_in_class[n_components] = 1;

        while( !IS_QUEUE_EMPTY(queue) )
        {
            REMOVE_FROM_QUEUE( queue, current_poly );
            size = GET_OBJECT_SIZE( *polygons, current_poly );

            for_less( edge, 0, size )
            {
                neigh = polygons->neighbours[
                    POINT_INDEX(polygons->end_indices,current_poly,edge)];
                if( neigh >= 0 &&
                    polygon_classes[neigh] == not_done )
                {
                    polygon_classes[neigh] = (Smallest_int) n_components;
                    ++n_in_class[n_components];
                    INSERT_IN_QUEUE( queue, neigh );
                }
            }
        }

        DELETE_QUEUE( queue );

        ++n_components;
    }

    return( n_components );
}

int   separate_polygons(
    polygons_struct    *polygons,
    int                desired_index,
    object_struct      **out[] )
{
    int                point, ind, p_ind, poly, vertex, size, i, j, tmp;
    int                point_index, *new_point_ids, n_objects, comp, c;
    int                biggest;
    Smallest_int       *poly_classes;
    int                n_components, *n_in_class, *ordered;
    polygons_struct    *new_poly;

    ALLOC( poly_classes, polygons->n_items );
    ALLOC( n_in_class, 1000 );
    ALLOC( ordered, 1000 );

    n_components = make_connected_components( polygons, poly_classes,
                                              n_in_class );

    for_less( i, 0, n_components )
        ordered[i] = i;

    for_less( i, 0, n_components-1 )
    {
        biggest = i;
        for_less( j, i+1, n_components )
        {
            if( n_in_class[ordered[j]] > n_in_class[ordered[biggest]] )
                biggest = j;
        }

        tmp = ordered[i];
        ordered[i] = ordered[biggest];
        ordered[biggest] = tmp;
    }

    ALLOC( new_point_ids, polygons->n_points );

    n_objects = 0;

    for_less( c, 0, n_components )
    {
        if( desired_index >= 0 && c != desired_index )
            continue;

        comp = ordered[c];

        for_less( point, 0, polygons->n_points )
            new_point_ids[point] = -1;

        SET_ARRAY_SIZE( *out, n_objects, n_objects+1,
                        DEFAULT_CHUNK_SIZE);
        (*out)[n_objects] = create_object( POLYGONS );
        new_poly = get_polygons_ptr( (*out)[n_objects] );
        ++n_objects;
        initialize_polygons( new_poly, WHITE, NULL );
        if( desired_index >= 0 )
        {
            new_poly->points = polygons->points;
            new_poly->normals = polygons->normals;
            new_poly->indices = polygons->indices;
            new_poly->n_items = 0;

            for_less( poly, 0, polygons->n_items )
            {
                if( poly_classes[poly] != (Smallest_int) comp )
                    continue;
                size = GET_OBJECT_SIZE( *polygons, poly );
                ++new_poly->n_items;
                for_less( vertex, 0, size )
                {
                    point_index = polygons->indices[
                              POINT_INDEX(polygons->end_indices,poly,vertex)];
                    if( new_point_ids[point_index] < 0 )
                        new_point_ids[point_index] = 0;
                }
            }

            ALLOC( new_poly->end_indices, new_poly->n_items );

            ind = 0;
            for_less( point, 0, polygons->n_points )
            {
                if( new_point_ids[point] >= 0 )
                {
                    new_point_ids[point] = ind;
                    new_poly->points[ind] = new_poly->points[point];
                    new_poly->normals[ind] = new_poly->normals[point];
                    ++ind;
                }
            }

            new_poly->n_points = ind;

            p_ind = 0;
            ind = 0;
            for_less( poly, 0, polygons->n_items )
            {
                if( poly_classes[poly] != (Smallest_int) comp )
                    continue;

                size = GET_OBJECT_SIZE( *polygons, poly );

                for_less( vertex, 0, size )
                {
                    point_index = polygons->indices[
                         POINT_INDEX(polygons->end_indices,poly,vertex)];
                    new_poly->indices[ind] = new_point_ids[point_index];
                    ++ind;
                }

                new_poly->end_indices[p_ind] = ind;
                ++p_ind;
            }
        }
        else
        {
            ind = 0;
            for_less( poly, 0, polygons->n_items )
            {
                if( poly_classes[poly] != (Smallest_int) comp )
                    continue;

                size = GET_OBJECT_SIZE( *polygons, poly );
                for_less( vertex, 0, size )
                {
                    point_index = polygons->indices[
                              POINT_INDEX(polygons->end_indices,poly,vertex)];

                    if( new_point_ids[point_index] < 0 )
                    {
                        new_point_ids[point_index] = new_poly->n_points;
                        ADD_ELEMENT_TO_ARRAY( new_poly->points,
                                              new_poly->n_points,
                                              polygons->points[point_index],
                                              DEFAULT_CHUNK_SIZE );
                        --new_poly->n_points;
                        ADD_ELEMENT_TO_ARRAY( new_poly->normals,
                                              new_poly->n_points,
                                              polygons->normals[point_index],
                                              DEFAULT_CHUNK_SIZE );
                    }

                    ADD_ELEMENT_TO_ARRAY( new_poly->indices, ind,
                                          new_point_ids[point_index],
                                          DEFAULT_CHUNK_SIZE );
                }

                ADD_ELEMENT_TO_ARRAY( new_poly->end_indices, new_poly->n_items,
                                      ind, DEFAULT_CHUNK_SIZE );
            }

            REALLOC( new_poly->points, new_poly->n_points );
            REALLOC( new_poly->normals, new_poly->n_points );
            REALLOC( new_poly->end_indices, new_poly->n_items );
            REALLOC( new_poly->indices, ind );
        }
    }

    FREE( poly_classes );
    FREE( new_point_ids );
    FREE( n_in_class );
    FREE( ordered );

    return( n_objects );
}
