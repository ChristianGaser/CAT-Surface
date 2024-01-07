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
static char rcsid[] = "$Header: /private-cvsroot/libraries/bicpl/Geometry/line_circle.c,v 1.12 2005/08/17 22:30:25 bert Exp $";
#endif

/* ----------------------------- MNI Header -----------------------------------
@NAME       : create_line_circle
@INPUT      : centre
              plane_axis
              x_radius
              y_radius
              n_points
@OUTPUT     : lines
@RETURNS    : 
@DESCRIPTION: Creates an ellipse with the given centre and x and y radius,
              as a closed line structure.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  void  create_line_circle(
    Point            *centre,
    int              plane_axis,
    Real             x_radius,
    Real             y_radius,
    int              n_points,
    lines_struct     *lines )
{
    int      i, a1, a2;
    Real     angle;

    initialize_lines( lines, WHITE );

    lines->n_points = n_points;
    lines->n_items = 1;

    ALLOC( lines->points, n_points );
    ALLOC( lines->end_indices, 1 );
    ALLOC( lines->indices, n_points+1 );

    a1 = (plane_axis+1) % 3;
    a2 = (plane_axis+2) % 3;

    for_less( i, 0, n_points )
    {
        angle = 2.0 * PI * (Real) i / (Real) n_points;
        Point_coord(lines->points[i],a1) = (Point_coord_type)
                   ( (Real) Point_coord(*centre,a1) + x_radius * cos( angle ) );
        Point_coord(lines->points[i],a2) = (Point_coord_type)
                   ( (Real) Point_coord(*centre,a2) + y_radius * sin( angle ) );
        Point_coord(lines->points[i],plane_axis) =
                                           Point_coord(*centre,plane_axis);
    }

    for_less( i, 0, n_points )
        lines->indices[i] = i;

    lines->indices[n_points] = 0;

    lines->end_indices[0] = n_points+1;
}
