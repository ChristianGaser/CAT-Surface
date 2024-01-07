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
static char rcsid[] = "$Header: /private-cvsroot/libraries/bicpl/Geometry/surface_area.c,v 1.8 2005/08/17 22:30:25 bert Exp $";
#endif

/* ----------------------------- MNI Header -----------------------------------
@NAME       : get_polygon_2d_area
@INPUT      : n_points
              points
@OUTPUT     :
@RETURNS    : area
@DESCRIPTION: Computes the area of a 2D polygon (convex or not).
@METHOD     :
@GLOBALS    :
@CALLS      :
@CREATED    :         1993    David MacDonald
@MODIFIED   :
---------------------------------------------------------------------------- */

BICAPI  Real  get_polygon_2d_area(
    int      n_points,
    Point    points[] )
{
    int    i, next_i;
    Real   area;

    area = 0.0;

    for_less( i, 0, n_points )
    {
        next_i = (i + 1) % n_points;
        area += (Real) Point_x(points[i])      * (Real) Point_y(points[next_i])-
                (Real) Point_x(points[next_i]) * (Real) Point_y(points[i]);
    }

    return( FABS( area / 2.0 ) );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : get_polygon_surface_area
@INPUT      : n_points
              points
@OUTPUT     : 
@RETURNS    : surface area
@DESCRIPTION: Returns the surface area of the 3D polygon (convex or not).
              If it is not planar, then the surface area is just an
              approximation.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  Real  get_polygon_surface_area(
    int      n_points,
    Point    points[] )
{
    int    i;
    Vector prev, this, cross_prod, sum;
    Real   surface_area;

    fill_Vector( sum, 0.0, 0.0, 0.0 );
    SUB_VECTORS( this, points[1], points[0] );

    for_less( i, 2, n_points )
    {
        prev = this;
        SUB_VECTORS( this, points[i], points[0] );
        CROSS_VECTORS( cross_prod, prev, this );
        ADD_VECTORS( sum, sum, cross_prod );
    }

    surface_area = MAGNITUDE( sum ) / 2.0;

    return( surface_area );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : get_polygons_surface_area
@INPUT      : polygons
@OUTPUT     : 
@RETURNS    : surface area
@DESCRIPTION: Computes the sum of the surface areas of all the polygons.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  Real  get_polygons_surface_area(
    polygons_struct  *polygons )
{
    int      poly, size;
    Point    points[MAX_POINTS_PER_POLYGON];
    Real     surface_area;

    surface_area = 0.0;

    for_less( poly, 0, polygons->n_items )
    {
        size = get_polygon_points( polygons, poly, points );
        surface_area += get_polygon_surface_area( size, points );
    }

    return( surface_area );
}
