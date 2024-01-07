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
static char rcsid[] = "$Header: /private-cvsroot/libraries/bicpl/Geometry/path_surface.c,v 1.8 2005/08/17 22:30:25 bert Exp $";
#endif

static  void  follow_path(
    int   polygon1,
    int   polygon2,
    int   end_indices[],
    int   neighbours[],
    int   distances[],
    int   path_length,
    int   path[] );

#define  INVALID_DISTANCE   -1

/* ----------------------------- MNI Header -----------------------------------
@NAME       : find_path_between_polygons
@INPUT      : polygon1
              polygon2
              n_polygons
              end_indices
              visibilities
@OUTPUT     : path_exists
              path_length
              path
@RETURNS    : 
@DESCRIPTION: Finds the shortest path between two polygons, in terms of the
              number of intervening polygons.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI   void   find_path_between_polygons(
    int            polygon1,
    int            polygon2,
    int            n_polygons,
    int            end_indices[],
    Smallest_int   visibilities[],
    int            neighbours[],
    BOOLEAN        *path_exists,
    int            *path_length,
    int            *path[] )
{
    BOOLEAN    found_polygon1;
    int        start_index, end_index, neighbour;
    int        i, n, *distances, distance, curr_index, prev_index;
    int        n_step[2], n_step_alloced[2], *step_list[2];

    ALLOC( distances, n_polygons );

    for_less( i, 0, n_polygons )
        distances[i] = INVALID_DISTANCE;

    n_step[0] = 0;
    n_step[1] = 0;
    n_step_alloced[0] = 0;
    n_step_alloced[1] = 0;
    prev_index = 0;
    curr_index = 1;

    ADD_ELEMENT_TO_ARRAY_WITH_SIZE( step_list[curr_index],
                  n_step_alloced[curr_index], n_step[curr_index],
                  polygon2, DEFAULT_CHUNK_SIZE );
    distances[polygon2] = 0;

    distance = 0;

    found_polygon1 = (polygon1 == polygon2);

    while( !found_polygon1 && n_step[curr_index] > 0 )
    {
        ++distance;
        curr_index = 1 - curr_index;
        prev_index = 1 - prev_index;

        n_step[curr_index] = 0;

        for_less( i, 0, n_step[prev_index] )
        {
            start_index = START_INDEX( end_indices,
                                       step_list[prev_index][i] );
            end_index = end_indices[step_list[prev_index][i]];

            for_less( n, start_index, end_index )
            {
                neighbour = neighbours[n];

                if( neighbour >= 0 &&
                    (visibilities == (Smallest_int *) 0 ||
                     visibilities[neighbour]) &&
                    distances[neighbour] == INVALID_DISTANCE )
                {
                    if( neighbour == polygon1 )
                    {
                        found_polygon1 = TRUE;
                        break;
                    }
                    else
                    {
                        ADD_ELEMENT_TO_ARRAY_WITH_SIZE( step_list[curr_index],
                          n_step_alloced[curr_index], n_step[curr_index],
                          neighbour, DEFAULT_CHUNK_SIZE );
                        distances[neighbour] = distance;
                    }
                }
            }

            if( found_polygon1 )
                break;
        }
    }

    if( n_step_alloced[0] > 0 )
    {
        FREE( step_list[0] );
    }

    if( n_step_alloced[1] > 0 )
    {
        FREE( step_list[1] );
    }

    if( found_polygon1 )
    {
        *path_exists = TRUE;

        *path_length = distance + 1;

        ALLOC( *path, *path_length );

        follow_path( polygon1, polygon2, end_indices,
                     neighbours, distances, *path_length, *path );
    }
    else
    {
        *path_exists = FALSE;
    }

    FREE( distances );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : follow_path
@INPUT      : polygon1
              polygon2
              end_indices
              neighbours
              distances
              path_length
@OUTPUT     : path
@RETURNS    : 
@DESCRIPTION: Fills in the path from polygon1 to polygon2, given the distances
              of every polygon from polygon1.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

static  void  follow_path(
    int   polygon1,
    int   polygon2,
    int   end_indices[],
    int   neighbours[],
    int   distances[],
    int   path_length,
    int   path[] )
{
    int   n, start_index, end_index, path_index, neighbour, current_polygon;

    path[0] = polygon1;
    path[path_length-1] = polygon2;

    current_polygon = polygon1;

    for_less( path_index, 1, path_length-1 )
    {
        start_index = START_INDEX( end_indices, current_polygon );
        end_index = end_indices[current_polygon];

        for_less( n, start_index, end_index )
        {
            neighbour = neighbours[n];

            if( neighbour >= 0 &&
                distances[neighbour] == path_length - path_index - 1 )
            {
                break;
            }
        }

        if( n == end_index )
            handle_internal_error( "follow_path" );

        current_polygon = neighbour;
        path[path_index] = current_polygon;
    }
}
