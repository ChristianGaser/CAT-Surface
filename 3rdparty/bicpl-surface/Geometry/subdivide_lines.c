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
static char rcsid[] = "$Header: /private-cvsroot/libraries/bicpl/Geometry/subdivide_lines.c,v 1.11 2005/08/17 22:30:25 bert Exp $";
#endif

/* ----------------------------- MNI Header -----------------------------------
@NAME       : subdivide_line
@INPUT      : l
@OUTPUT     : new_n_points
              new_points
              new_n_lines
              new_end_indices
              new_n_indices
              new_indices
@RETURNS    : 
@DESCRIPTION: Subdivides the l'th line segment, adding to the new line.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

static  void  subdivide_line(
    lines_struct      *lines,
    int               l,
    int               *new_n_points,
    Point             *new_points[],
    int               *new_n_lines,
    int               *new_end_indices[],
    int               *new_n_indices,
    int               *new_indices[] )
{
    int     edge, size, p1, p2, midpoint_index;
    Point   midpoint;

    size = GET_OBJECT_SIZE( *lines, l );

    for_less( edge, 0, size-1 )
    {
        p1 = lines->indices[POINT_INDEX(lines->end_indices,l,edge)];
        p2 = lines->indices[POINT_INDEX(lines->end_indices,l,edge+1)];
        midpoint_index = *new_n_points;
        INTERPOLATE_POINTS( midpoint, lines->points[p1],
                            lines->points[p2], 0.5 );
        ADD_ELEMENT_TO_ARRAY( *new_points, *new_n_points,
                              midpoint, DEFAULT_CHUNK_SIZE );

        if( edge == 0 )
        {
            ADD_ELEMENT_TO_ARRAY( *new_indices, *new_n_indices,
                                  p1, DEFAULT_CHUNK_SIZE );
        }

        ADD_ELEMENT_TO_ARRAY( *new_indices, *new_n_indices,
                              midpoint_index, DEFAULT_CHUNK_SIZE );
        ADD_ELEMENT_TO_ARRAY( *new_indices, *new_n_indices,
                              p2, DEFAULT_CHUNK_SIZE );
    }

    ADD_ELEMENT_TO_ARRAY( *new_end_indices, *new_n_lines,
                          *new_n_indices, DEFAULT_CHUNK_SIZE );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : general_subdivide_lines
@INPUT      : lines
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: Subdivides a line structure of any topology.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

static  void  general_subdivide_lines(
    lines_struct  *lines )
{
    int                size, i, new_n_points, new_n_indices, new_n_lines;
    int                *new_indices, *new_end_indices;
    Point              *new_points;

    new_n_points = lines->n_points;
    new_n_lines = 0;
    new_n_indices = 0;
    new_points = NULL;

    size = 0;
    SET_ARRAY_SIZE( new_points, size, new_n_points, DEFAULT_CHUNK_SIZE );

    for_less( i, 0, new_n_points )
        new_points[i] = lines->points[i];

    for_less( i, 0, lines->n_items )
    {
        subdivide_line( lines, i, &new_n_points, &new_points,
                        &new_n_lines, &new_end_indices,
                        &new_n_indices, &new_indices );
    }

    delete_lines( lines );

    ALLOC( lines->colours, 1);
    lines->colours[0] = WHITE;
    lines->colour_flag = ONE_COLOUR;
    lines->n_points = new_n_points;
    lines->points = new_points;
    lines->n_items = new_n_lines;
    lines->end_indices = new_end_indices;
    lines->indices = new_indices;
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : subdivide_closed_curve
@INPUT      : lines
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: Subdivides a lines structure of a single loop topology, into
              similar topology.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

static  void  subdivide_closed_curve(
    lines_struct      *lines )
{
    int     i;
    Point   p1, p2, midpoint;

    REALLOC( lines->points, 2 * lines->n_points );

    for( i = lines->n_points-1;  i >= 0;  --i )
        lines->points[2*i] = lines->points[i];

    for_less( i, 0, lines->n_points )
    {
        p1 = lines->points[2* i];
        p2 = lines->points[2* ((i+1)%lines->n_points)];
        INTERPOLATE_POINTS( midpoint, p1, p2, 0.5 );

        lines->points[2*i+1] = midpoint;
    }

    lines->n_points = 2 * lines->n_points;
    lines->end_indices[0] = lines->n_points + 1;

    REALLOC( lines->indices, lines->n_points + 1 );

    for_less( i, 0, lines->n_points + 1 )
        lines->indices[i] = i % lines->n_points;
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : is_single_closed_curve
@INPUT      : 
@OUTPUT     : 
@RETURNS    : TRUE if single loop
@DESCRIPTION: Checks if lines is a single closed loop, with indices
              0, 1, 2, lines->n_points-1, 0.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  BOOLEAN is_single_closed_curve(
    lines_struct   *lines )
{
    int      i;
    BOOLEAN  single_closed_curve;

    single_closed_curve = FALSE;

    if( lines->n_items == 1 &&
        lines->n_points+1 == lines->end_indices[0] )
    {
        single_closed_curve = TRUE;
        for_less( i, 0, lines->n_points + 1 )
        {
            if( lines->indices[i] != (i % lines->n_points) )
            {
                single_closed_curve = FALSE;
                break;
            }
        }
    }

    return( single_closed_curve );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : subdivide_lines
@INPUT      : lines
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: Subdivides the lines, splitting each segment in half.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  void  subdivide_lines(
    lines_struct  *lines )
{
    if( is_single_closed_curve( lines ) )
        subdivide_closed_curve( lines );
    else
        general_subdivide_lines( lines );
}
