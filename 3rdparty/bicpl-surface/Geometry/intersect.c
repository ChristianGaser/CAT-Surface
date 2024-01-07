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
static char rcsid[] = "$Header: /private-cvsroot/libraries/bicpl/Geometry/intersect.c,v 1.13 2005/08/17 22:30:25 bert Exp $";
#endif

/* ----------------------------- MNI Header -----------------------------------
@NAME       : line_segment_intersects_plane
@INPUT      : p1
              p2
              plane_origin
              plane_normal
@OUTPUT     : intersection_point
@RETURNS    : TRUE if line segment intersects plane
@DESCRIPTION: Tests if a line segment from p1 to p2 intersects the plane,
              passing back the intersection point, if so.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : Jul 12, 1995    David MacDonald - rewrote to be clearer.
---------------------------------------------------------------------------- */

BICAPI  BOOLEAN   line_segment_intersects_plane(
    Point   *p1,
    Point   *p2,
    Point   *plane_origin,
    Vector  *plane_normal,
    Point   *intersection_point )
{
    Real     t, t1, t2;
    Vector   v1, v2;

    SUB_POINTS( v1, *p1, *plane_origin );
    SUB_POINTS( v2, *p2, *plane_origin );

    t1 = DOT_VECTORS( v1, *plane_normal );
    t2 = DOT_VECTORS( v2, *plane_normal );

    if( t1 == t2 )
        return( FALSE );

    if( t1 == 0.0 )
        *intersection_point = *p1;
    else if( t2 == 0.0 )
        *intersection_point = *p2;
    else
    {
        t = t1 / (t1 - t2);
 
        if( t < 0.0 || t > 1.0 )
            return( FALSE );

        INTERPOLATE_POINTS( *intersection_point, *p1, *p2, t );
    }

    return( TRUE );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : get_nearest_point_on_lines
@INPUT      : origin1
              delta1
              origin2
              delta2
@OUTPUT     : nearest_pt
@RETURNS    : TRUE if near point found
@DESCRIPTION: Finds the nearest point on the first line to the second line,
              or returns FALSE if they are parallel.  If it can, passes back
              the intersection point.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  BOOLEAN  get_nearest_point_on_lines(
    Point   *origin1,
    Vector  *delta1,
    Point   *origin2,
    Vector  *delta2,
    Point   *nearest_point )
{
    BOOLEAN   intersects;
    Real      t, bottom;
    Real      d11, d12, d22;
    Real      o11, o12, o21, o22;

    d11 = DOT_VECTORS( *delta1, *delta1 );
    d12 = DOT_VECTORS( *delta1, *delta2 );
    d22 = DOT_VECTORS( *delta2, *delta2 );

    bottom = d11 * d22 - d12 * d12;

    intersects = (bottom != 0.0);

    if( intersects )
    {
        o11 = DOT_POINT_VECTOR( *origin1, *delta1 );
        o12 = DOT_POINT_VECTOR( *origin1, *delta2 );
        o21 = DOT_POINT_VECTOR( *origin2, *delta1 );
        o22 = DOT_POINT_VECTOR( *origin2, *delta2 );

        t = -(d22 * (o11 - o21) + d12 * (o22 - o12)) / bottom;

        GET_POINT_ON_RAY( *nearest_point, *origin1, *delta1, t );
    }

    return( intersects );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : clip_line_to_box
@INPUT      : origin
              direction
              x_min
              x_max
              y_min
              y_max
              z_min
              z_max
@OUTPUT     : t_min
              t_max
@RETURNS    : TRUE if line intersects box
@DESCRIPTION: Determines if the line intersects the box, and if so, passes
              back the min and max multiple of direction vector for the box
              intersection.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  BOOLEAN  clip_line_to_box(
    Point    *origin,
    Vector   *direction,
    Real     x_min,
    Real     x_max,
    Real     y_min,
    Real     y_max,
    Real     z_min,
    Real     z_max,
    Real     *t_min,
    Real     *t_max )
{
    BOOLEAN  first;
    Real     dir, org, t1, t2, min_t, max_t;

    *t_min = 0.0;
    *t_max = -1.0;

    dir = (Real) Vector_x( *direction );
    org = (Real) Point_x( *origin );

    if( dir == 0.0 )
    {
        if( org < x_min || org > x_max )
            return( FALSE );

        dir = (Real) Vector_y( *direction );
        org = (Real) Point_y( *origin );

        if( dir < 0.0 )
        {
            first = FALSE;
            min_t = (y_max - org) / dir;
            max_t = (y_min - org) / dir;
        }
        else if( dir > 0.0 )
        {
            first = FALSE;
            min_t = (y_min - org) / dir;
            max_t = (y_max - org) / dir;
        }
        else
        {
            if( org < y_min || org > y_max )
                return( FALSE );

            first = TRUE;
        }

        dir = (Real) Vector_z( *direction );
        org = (Real) Point_z( *origin );

        if( dir < 0.0 )
        {
            t1 = (z_max - org) / dir;
            t2 = (z_min - org) / dir;

            if( first || t1 > min_t )
                min_t = t1;
            if( first || t2 < max_t )
                max_t = t2;
            first = FALSE;
        }
        else if( dir > 0.0 )
        {
            t1 = (z_min - org) / dir;
            t2 = (z_max - org) / dir;

            if( first || t1 > min_t )
                min_t = t1;
            if( first || t2 < max_t )
                max_t = t2;
            first = FALSE;
        }
        else
        {
            if( org < z_min || org > z_max )
                return( FALSE );
        }

        if( first )
        {
            min_t = 0.0;
            max_t = 0.0;
        }

        *t_min = min_t;
        *t_max = max_t;
        return( min_t <= max_t );
    }

    if( dir > 0.0 )
    {
        min_t = (x_min - org) / dir;
        max_t = (x_max - org) / dir;
    }
    else
    {
        min_t = (x_max - org) / dir;
        max_t = (x_min - org) / dir;
    }

    dir = (Real) Vector_y( *direction );
    org = (Real) Point_y( *origin );

    if( dir == 0.0 )
    {
        if( org < y_min || org > y_max )
            return( FALSE );
    }
    else
    {
        if( dir > 0.0 )
        {
            t1 = (y_min - org) / dir;
            t2 = (y_max - org) / dir;
        }
        else
        {
            t1 = (y_max - org) / dir;
            t2 = (y_min - org) / dir;
        }

        if( t1 > min_t )
        {
            if( t1 > max_t )
                return( FALSE );
            min_t = t1;
        }
        if( t2 < max_t )
        {
            if( t2 < min_t )
                return( FALSE );
            max_t = t2;
        }
    }

    dir = (Real) Vector_z( *direction );
    org = (Real) Point_z( *origin );

    if( dir == 0.0 )
    {
        if( org < z_min || org > z_max )
            return( FALSE );
    }
    else
    {
        if( dir > 0.0 )
        {
            t1 = (z_min - org) / dir;
            t2 = (z_max - org) / dir;
        }
        else
        {
            t1 = (z_max - org) / dir;
            t2 = (z_min - org) / dir;
        }

        if( t1 > min_t )
        {
            if( t1 > max_t )
                return( FALSE );
            min_t = t1;
        }
        if( t2 < max_t )
        {
            if( t2 < min_t )
                return( FALSE );
            max_t = t2;
        }
    }

    *t_min = min_t;
    *t_max = max_t;

    return( min_t <= max_t );
}
