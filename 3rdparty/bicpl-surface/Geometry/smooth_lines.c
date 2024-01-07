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
static char rcsid[] = "$Header: /private-cvsroot/libraries/bicpl/Geometry/smooth_lines.c,v 1.12 2005/08/17 22:30:25 bert Exp $";
#endif

/* ----------------------------- MNI Header -----------------------------------
@NAME       : smooth_lines
@INPUT      : lines
              smooth_length
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: Removes segments that are smaller than the smooth length.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  void  smooth_lines(
    lines_struct  *lines,
    Real          smooth_length )
{
    int       *new_ids, n_points, n_items, n_indices, *indices, *end_indices;
    int       l, p, point_index, size;
    BOOLEAN   keep_point;
    Point     prev;
    Real      dist_to_prev;
    Point     *points;
    Colour    colour;

    n_points = 0;
    n_items = 0;
    n_indices = 0;
    indices = NULL;
    end_indices = NULL;
    points = NULL;

    ALLOC( new_ids, lines->n_points );

    for_less( l, 0, lines->n_points )
        new_ids[l] = -1;

    for_less( l, 0, lines->n_items )
    {
        size = GET_OBJECT_SIZE( *lines, l );

        for_less( p, 0, size )
        {
            point_index = lines->indices[
                             POINT_INDEX(lines->end_indices,l,p)];

            if( p == 0 || p == size-1 )
            {
                keep_point = TRUE;
            }
            else
            {
                dist_to_prev = distance_between_points(&prev,
                                         &lines->points[point_index]);

                keep_point = (dist_to_prev >= smooth_length);
            }

            if( keep_point )
            {
                prev = lines->points[point_index];

                if( new_ids[point_index] < 0 )
                {
                    new_ids[point_index] = n_points;

                    ADD_ELEMENT_TO_ARRAY( points, n_points,
                              lines->points[point_index], DEFAULT_CHUNK_SIZE );
                }

                ADD_ELEMENT_TO_ARRAY( indices, n_indices,
                                     new_ids[point_index], DEFAULT_CHUNK_SIZE );
            }
        }

        ADD_ELEMENT_TO_ARRAY( end_indices, n_items,
                              n_indices, DEFAULT_CHUNK_SIZE );
    }

    colour = lines->colours[0];
    delete_lines( lines );

    lines->colour_flag = ONE_COLOUR;
    ALLOC( lines->colours, 1 );
    lines->colours[0] = colour;
    lines->n_points = n_points;
    lines->points = points;
    lines->n_items = n_items;
    lines->indices = indices;
    lines->end_indices = end_indices;

    FREE( new_ids );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : create_line_spline
@INPUT      : lines
              n_curve_segments
@OUTPUT     : new_lines
@RETURNS    : 
@DESCRIPTION: Creates a new line which is a cubic spline interpolation of the
              input lines, with n_curve_segments per input line segment.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  void  create_line_spline(
    lines_struct  *lines,
    int           n_curve_segments,
    lines_struct  *new_lines )
{
    int       c, l, p, point_index1, point_index2, pt_index, line_size, segment;
    BOOLEAN   wrap_around;
    Point     points[4], point;
    Real      u;

    initialize_lines( new_lines, lines->colours[0] );

    new_lines->colours[0] = lines->colours[0];
    new_lines->line_thickness = lines->line_thickness;

    new_lines->n_points = 0;
    new_lines->n_items = 0;

    for_less( l, 0, lines->n_items )
    {
        line_size = GET_OBJECT_SIZE( *lines, l );

        point_index1 = lines->indices[POINT_INDEX(lines->end_indices,l,0)];
        point_index2 = lines->indices[
                         POINT_INDEX(lines->end_indices,l,line_size-1)];

        wrap_around = (point_index1 == point_index2);

        start_new_line( new_lines );

        add_point_to_line( new_lines, &lines->points[point_index1] );

        for_less( p, 0, line_size-1 )
        {
            if( p == 0 )
            {
                if( wrap_around && line_size - 2 >= 0 )
                    pt_index = line_size-2;
                else
                    pt_index = 0;
            }
            else
                pt_index = p - 1;

            points[0] = lines->points[lines->indices[
                          POINT_INDEX(lines->end_indices,l,pt_index)]];

            points[1] = lines->points[lines->indices[
                          POINT_INDEX(lines->end_indices,l,p)]];

            points[2] = lines->points[lines->indices[
                          POINT_INDEX(lines->end_indices,l,p+1)]];

            if( p == line_size-2 )
            {
                if( wrap_around && 1 < line_size )
                    pt_index = 1;
                else
                    pt_index = line_size-1;
            }
            else
                pt_index = p + 2;

            points[3] = lines->points[lines->indices[
                          POINT_INDEX(lines->end_indices,l,pt_index)]];

            for_inclusive( segment, 1, n_curve_segments )
            {
                u = (Real) segment / (Real) n_curve_segments;

                for_less( c, 0, N_DIMENSIONS )
                {
                    Point_coord(point,c) = (Point_coord_type)
                                           cubic_interpolate( u,
                                            (Real) Point_coord(points[0],c),
                                            (Real) Point_coord(points[1],c),
                                            (Real) Point_coord(points[2],c),
                                            (Real) Point_coord(points[3],c) );
                }

                add_point_to_line( new_lines, &point );
            }
        }
    }
}
