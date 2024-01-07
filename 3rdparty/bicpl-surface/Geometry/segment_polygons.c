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
static char rcsid[] = "$Header: /private-cvsroot/libraries/bicpl/Geometry/segment_polygons.c,v 1.9 2005/08/17 22:30:25 bert Exp $";
#endif

/* ----------------------------- MNI Header -----------------------------------
@NAME       : remove_invisible_polygons
@INPUT      : polygons
              visibilities
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: Removes any polygon that is not marked as visible from the
              polygons.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  void   remove_invisible_polygons(
    polygons_struct  *polygons,
    Smallest_int     visibilities[] )
{
    int              tmp_n_points, i, n_points, n_items, *indices, *end_indices;
    int              p, poly, n_indices, size, point_index;
    int              *new_ids;
    Point            *points;
    Vector           *normals;
    Colour           col;

    if( visibilities == NULL )
        return;

    n_points = 0;
    n_items = 0;
    n_indices = 0;
    indices = NULL;
    end_indices = NULL;
    points = NULL;
    normals = NULL;

    ALLOC( new_ids, polygons->n_points );

    for_less( i, 0, polygons->n_points )
        new_ids[i] = -1;

    for_less( poly, 0, polygons->n_items )
    {
        if( visibilities[poly] )
        {
            size = GET_OBJECT_SIZE( *polygons, poly );

            for_less( p, 0, size )
            {
                point_index = polygons->indices[
                        POINT_INDEX(polygons->end_indices,poly,p)];

                if( new_ids[point_index] < 0 )
                {
                    new_ids[point_index] = n_points;

                    tmp_n_points = n_points;

                    ADD_ELEMENT_TO_ARRAY( points, tmp_n_points,
                         polygons->points[point_index], DEFAULT_CHUNK_SIZE );

                    ADD_ELEMENT_TO_ARRAY( normals, n_points,
                          polygons->normals[point_index], DEFAULT_CHUNK_SIZE );
                }

                ADD_ELEMENT_TO_ARRAY( indices, n_indices,
                         new_ids[point_index], DEFAULT_CHUNK_SIZE );
            }

            ADD_ELEMENT_TO_ARRAY( end_indices, n_items,
                                  n_indices, DEFAULT_CHUNK_SIZE );
        }
    }

    FREE( new_ids );

    col = polygons->colours[0];

    delete_polygons( polygons );

    ALLOC( polygons->colours, 1 );

    polygons->colour_flag = ONE_COLOUR;
    polygons->colours[0] = col;
    polygons->n_points = n_points;
    polygons->points = points;
    polygons->normals = normals;
    polygons->n_items = n_items;
    polygons->indices = indices;
    polygons->end_indices = end_indices;
}
