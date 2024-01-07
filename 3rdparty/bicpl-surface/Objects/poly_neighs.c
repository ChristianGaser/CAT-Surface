/* ----------------------------------------------------------------------------
@COPYRIGHT  :
              Copyright 1993,1994,1995 David MacDonald,
              McConnell Brain Imaging Centre,
              Montreal Neurological Institute, McGill University.
              Permission to use, copy, modify, and distribute this
              software and its documentation for any purpose and without
              fee is hereby granted, provided that the above copyright
              notice appear in all copies.  The author and McGill University
              make no representations about the suitability of this
              software for any purpose.  It is provided "as is" without
              express or implied warranty.
---------------------------------------------------------------------------- */

#include "bicpl_internal.h"

#ifndef lint
static char rcsid[] = "$Header: /private-cvsroot/libraries/bicpl/Objects/poly_neighs.c,v 1.23 2005/08/17 22:28:27 bert Exp $";
#endif

#define  SMALL_CHUNK_SIZE          4

static   void   create_polygon_neighbours(
    polygons_struct  *polygons,
    int              neighbours[] );

/* ----------------------------- MNI Header -----------------------------------
@NAME       : check_polygons_neighbours_computed
@INPUT      : polygons
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: Creates the polygons neighbours (neighbouring polygon index for
              each polygon edge), if necessary.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  void  check_polygons_neighbours_computed(
    polygons_struct   *polygons )
{
    if( polygons->neighbours == NULL && polygons->n_items > 0 )
    {
        ALLOC( polygons->neighbours,polygons->end_indices[polygons->n_items-1]);
        create_polygon_neighbours( polygons, polygons->neighbours );
    }
}

BICAPI  void  delete_polygon_point_neighbours(
    polygons_struct  *polygons,
    int              n_point_neighbours[],
    int              *point_neighbours[],
    Smallest_int     interior_flags[],
    int              *point_polygons[] )
{
    int   i;

    FREE( n_point_neighbours );

    for_less( i, 0, polygons->n_points )
        FREE( point_neighbours[i] );
    FREE( point_neighbours );

    if( interior_flags != NULL )
        FREE( interior_flags );

    if( point_polygons != NULL )
    {
        FREE( point_polygons[0] );
        FREE( point_polygons );
    }
}

static  void  insert_neighbours(
    int    n_to_add,
    int    indices[],
    int    n_nodes,
    int    *n_neighbours,
    int    *neighbours[] )
{
    int       p1, p2, first_index, start, i, n_to_insert, n_to_do, last_index;
    BOOLEAN   wrapped;

    p1 = indices[0];
    for_less( first_index, 0, *n_neighbours )
    {
        if( (*neighbours)[first_index] == p1 + n_nodes )
            break;
    }

    wrapped = FALSE;

    if( first_index < *n_neighbours )
    {
        start = first_index-1;
        while( start >= 0 && (*neighbours)[start] < n_nodes )
            --start;
        ++start;

        n_to_insert = first_index - start;
        for_down( i, n_to_add - 1, 0 )
            indices[n_to_insert+i] = indices[i];

        for_less( i, 0, n_to_insert )
        {
            indices[i] = (*neighbours)[start+i];
            (*neighbours)[start+i] = -1;
        }
        (*neighbours)[first_index] = -1;
        n_to_add += n_to_insert;
        if( indices[n_to_add-1] == indices[0] )
        {
            --n_to_add;
            wrapped = TRUE;
        }
    }

    p2 = indices[n_to_add-1];
    for_less( last_index, 0, *n_neighbours )
    {
        if( (*neighbours)[last_index] == p2 )
            break;
    }

    if( last_index == *n_neighbours  )
    {
        if( !wrapped )
            indices[n_to_add-1] += n_nodes;
        last_index = 0;
        n_to_do = *n_neighbours;
    }
    else
    {
        ++last_index;
        n_to_do = *n_neighbours-1;
    }

    for_less( i, 0, n_to_do )
    {
        if( (*neighbours)[(last_index+i)%*n_neighbours] >= 0 )
        {
            indices[n_to_add] = (*neighbours)[(last_index+i)%*n_neighbours];
            ++n_to_add;
        }
    }

    SET_ARRAY_SIZE( *neighbours, *n_neighbours, n_to_add, SMALL_CHUNK_SIZE);

    *n_neighbours = n_to_add;
    for_less( i, 0, n_to_add )
        (*neighbours)[i] = indices[i];
}

BICAPI   void   create_polygon_point_neighbours(
    polygons_struct  *polygons,
    BOOLEAN          across_polygons_flag,
    int              *n_point_neighbours_ptr[],
    int              **point_neighbours_ptr[],
    Smallest_int     *interior_flags_ptr[],
    int              **point_polygons_ptr[] )
{
    int                 edge, i0, i1, size, poly, total_neighbours, p0, p1;
    int                 points_to_add_alloced;
    int                 n_points, n;
    int                 *n_point_neighbours, **point_neighbours;
    int                 **point_polygons, point, index0, index1;
    int                 i, v, indices[MAX_POINTS_PER_POLYGON], ii;
    int                 n_to_add, *points_to_add;
    progress_struct     progress;

    if( across_polygons_flag && point_polygons_ptr != NULL )
    {
        print_error(
                "create_polygon_point_neighbours: conflicting argument.\n" );
        return;
    }

    points_to_add_alloced = 0;
    points_to_add = NULL;

    n_points = polygons->n_points;

    ALLOC( n_point_neighbours, n_points );
    for_less( point, 0, n_points )
        n_point_neighbours[point] = 0;

    ALLOC( point_neighbours, n_points );

    initialize_progress_report( &progress, FALSE, polygons->n_items,
                                "Neighbour-finding" );

    for_less( poly, 0, polygons->n_items )
    {

        size = GET_OBJECT_SIZE( *polygons, poly );
        for_less( v, 0, size )
        {
            indices[v] = polygons->indices[
                                    POINT_INDEX(polygons->end_indices,poly,v)];
        }

        for_less( i0, 0, size )
        {
            p0 = indices[i0];

            if( size + n_point_neighbours[p0] > points_to_add_alloced )
            {
                SET_ARRAY_SIZE( points_to_add, points_to_add_alloced,
                                size + n_point_neighbours[p0],
                                DEFAULT_CHUNK_SIZE);
                points_to_add_alloced = size + n_point_neighbours[p0];
            }

            n_to_add = 0;
            for_less( ii, 0, size-1 )
            {
                i1 = (i0 + 1 + ii) % size;
                if( !across_polygons_flag && ii != 0 && ii != size-2 )
                    continue;

                p1 = indices[i1];
                points_to_add[n_to_add] = p1;
                ++n_to_add;
            }

            insert_neighbours( n_to_add, points_to_add, n_points,
                               &n_point_neighbours[p0], &point_neighbours[p0] );
        }

        update_progress_report( &progress, poly+1 );
    }

    terminate_progress_report( &progress );

    if( points_to_add_alloced > 0 )
        FREE( points_to_add );

    total_neighbours = 0;
    for_less( i, 0, n_points )
    {
        REALLOC( point_neighbours[i], n_point_neighbours[i] );
        total_neighbours += n_point_neighbours[i];
    }

    if( interior_flags_ptr != NULL )
    {
        ALLOC( *interior_flags_ptr, n_points );
        for_less( point, 0, n_points )
        {
            (*interior_flags_ptr)[point] = (Smallest_int)
              (point_neighbours[point][n_point_neighbours[point]-1] < n_points);
        }
    }

    for_less( point, 0, n_points )
    {
        for_less( n, 0, n_point_neighbours[point] )
        {
            if( point_neighbours[point][n] >= n_points )
                point_neighbours[point][n] -= n_points;
        }
    }

    *n_point_neighbours_ptr = n_point_neighbours;   
    *point_neighbours_ptr = point_neighbours;   

    if( point_polygons_ptr == NULL )
        return;

    ALLOC( point_polygons, n_points );
    ALLOC( point_polygons[0], total_neighbours );

    for_less( point, 0, total_neighbours )
        point_polygons[0][point] = -1;

    for_less( point, 1, n_points )
    {
        point_polygons[point] = &point_polygons[point-1]
                                               [n_point_neighbours[point-1]];
    }

    initialize_progress_report( &progress, FALSE, polygons->n_items,
                                "Neighbour-finding step 2" );

    for_less( poly, 0, polygons->n_items )
    {
        size = GET_OBJECT_SIZE( *polygons, poly );

        i1 = polygons->indices[POINT_INDEX(polygons->end_indices,poly,0)];
        for_less( edge, 0, size )
        {
            i0 = i1;
            i1 = polygons->indices[
                  POINT_INDEX(polygons->end_indices,poly,(edge+1)%size)];

            for_less( index0, 0, n_point_neighbours[i0] )
            {
                if( point_neighbours[i0][index0] == i1 )
                    break;
            }

            for_less( index1, 0, n_point_neighbours[i1] )
            {
                if( point_neighbours[i1][index1] == i0 )
                    break;
            }

            if( point_polygons[i0][index0] < 0 )
                point_polygons[i0][index0] = poly;
            else
                point_polygons[i1][index1] = poly;
        }

        update_progress_report( &progress, poly+1 );
    }

    terminate_progress_report( &progress );

    *point_polygons_ptr = point_polygons;   
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : create_polygon_neighbours
@INPUT      : n_polygons
              indices
              end_indices
@OUTPUT     : neighbours
@RETURNS    : 
@DESCRIPTION: Computes the neighbours for each edge of the polygons
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

static   void   create_polygon_neighbours(
    polygons_struct  *polygons,
    int              neighbours[] )
{
    int                 i0, i1, size1, size2, n1, n2;
    int                 poly1, poly2, point1, point2, edge1, edge2;
    int                 *n_point_neighbours, **point_neighbours;
    int                 **point_polygons;
    progress_struct     progress;

    for_less( i0, 0, polygons->end_indices[polygons->n_items-1] )
        neighbours[i0] = -1;

    create_polygon_point_neighbours( polygons, FALSE, &n_point_neighbours,
                                     &point_neighbours, NULL, &point_polygons );

    initialize_progress_report( &progress, FALSE, polygons->n_items,
                                "Neighbour-finding step 2" );

    for_less( point1, 0, polygons->n_points )
    {
        for_less( n1, 0, n_point_neighbours[point1] )
        {
            point2 = point_neighbours[point1][n1];

            if( point2 <= point1 )
                continue;

            poly1 = point_polygons[point1][n1];
            if( poly1 < 0 )
                continue;

            for_less( n2, 0, n_point_neighbours[point2] )
            {
                if( point_neighbours[point2][n2] == point1 )
                    break;
            }

            if( n2 >= n_point_neighbours[point2] )
                handle_internal_error( "create_polygon_neighbours" );

            poly2 = point_polygons[point2][n2];
            if( poly2 < 0 )
                continue;

            size1 = GET_OBJECT_SIZE( *polygons, poly1 );
            for_less( edge1, 0, size1 )
            {
                i0 = polygons->indices[POINT_INDEX(polygons->end_indices,
                                                   poly1,edge1)];
                i1 = polygons->indices[POINT_INDEX(polygons->end_indices,
                                                   poly1,(edge1+1)%size1)];

                if( (i0 == point1 && i1 == point2) ||
                    (i1 == point1 && i0 == point2) )
                    break;
            }

            if( edge1 >= size1 )
                handle_internal_error( "create_polygon_neighbours" );

            size2 = GET_OBJECT_SIZE( *polygons, poly2 );
            for_less( edge2, 0, size2 )
            {
                i0 = polygons->indices[POINT_INDEX(polygons->end_indices,
                                                   poly2,edge2)];
                i1 = polygons->indices[POINT_INDEX(polygons->end_indices,
                                                   poly2,(edge2+1)%size2)];

                if( (i0 == point1 && i1 == point2) ||
                    (i1 == point1 && i0 == point2) )
                    break;
            }

            if( edge2 >= size2 )
                handle_internal_error( "create_polygon_neighbours" );

            neighbours[POINT_INDEX( polygons->end_indices, poly1, edge1 )] =
                                                               poly2;
            neighbours[POINT_INDEX( polygons->end_indices, poly2, edge2 )] =
                                                               poly1;
        }

        update_progress_report( &progress, point1+1 );
    }

    terminate_progress_report( &progress );

    delete_polygon_point_neighbours( polygons,
                                     n_point_neighbours, point_neighbours,
                                     NULL, point_polygons );

#ifdef  DEBUG
    for_less( i0, 0, polygons->end_indices[polygons->n_items-1] )
    {
        if( neighbours[i0] < 0 )
            break;
    }

    if( i0 < polygons->end_indices[polygons->n_items-1] )
        handle_internal_error( "create_polygon_neighbours: topology\n" );
#endif
}
