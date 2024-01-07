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
static char rcsid[] = "$Header: /private-cvsroot/libraries/bicpl/Geometry/poly_dist.c,v 1.10 2005/08/17 22:30:25 bert Exp $";
#endif

/* ----------------------------- MNI Header -----------------------------------
@NAME       : sq_distance_between_points
@INPUT      : p1
              p2
@OUTPUT     : 
@RETURNS    : squared distance between the points
@DESCRIPTION: 
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  Real  sq_distance_between_points(
    Point  *p1,
    Point  *p2 )
{
    Real   dx, dy, dz;

    dx = RPoint_x(*p1) - RPoint_x(*p2);
    dy = RPoint_y(*p1) - RPoint_y(*p2);
    dz = RPoint_z(*p1) - RPoint_z(*p2);

    return( dx * dx + dy * dy + dz * dz );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : point_segment_sq_distance
@INPUT      : p
              q1
              q2
@OUTPUT     : closest_point
@RETURNS    : squared distance
@DESCRIPTION: Returns the squared distance from a point to a line segment.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

static  Real  point_segment_sq_distance(
    Point   *p,
    Point   *q1,
    Point   *q2,
    Point   *closest_point )
{
    get_closest_point_on_line_segment( p, q1, q2, closest_point );

    return( sq_distance_between_points( p, closest_point ) );
}

BICAPI  Real  find_point_polygon_distance_sq(
    Point     *point,
    int       n_points,
    Point     poly_points[],
    Point     *closest_point )
{
    int      i, closest;
    Real     n_dot_n, t, closest_dist, dist, d1, d2;
    Vector   offset, o_a, normal;
    Point    seg1_point, seg2_point;

    /*--- first, find closest point on plane of polygon */

    find_polygon_normal( n_points, poly_points, &normal );

    n_dot_n = DOT_VECTORS( normal, normal );
    if( n_dot_n == 0.0 )
    {
        fill_Point( *closest_point, 0.0, 0.0, 0.0 );
        return( 1.0e30 );
    }

    SUB_POINTS( o_a, poly_points[0], *point );

    t = DOT_VECTORS( o_a, normal ) / n_dot_n;

    SCALE_VECTOR( offset, normal, t ); 
    ADD_POINT_VECTOR( *closest_point, *point, offset );

    /*--- if closest point on plane of polygon is within the polygon, then
          it is the closest point of the polygon */

    if( point_within_polygon( closest_point, n_points, poly_points, &normal ) )
        return( DOT_VECTORS( offset, offset ) );

    /*---  find the vertex of the polygon which is closest */

    closest = 0;
    closest_dist = 0.0;

    for_less( i, 0, n_points )
    {
        dist = sq_distance_between_points( point, &poly_points[i] );
        if( i == 0 || dist < closest_dist )
        {
            closest = i;
            closest_dist = dist;
        }
    }

    /*--- the closest point is on one of the two polygon edges touching the
          closest vertex */

    d1 = point_segment_sq_distance( point,
                                    &poly_points[(closest-1+n_points)%n_points],
                                    &poly_points[closest], &seg1_point );
    d2 = point_segment_sq_distance( point,
                                    &poly_points[closest],
                                    &poly_points[(closest+1)%n_points],
                                    &seg2_point );

    if( d1 < d2 )
    {
        *closest_point = seg1_point;
        closest_dist = d1;
    }
     else
    {
        *closest_point = seg2_point;
        closest_dist = d2;
    }

    return( closest_dist );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : find_point_polygon_distance
@INPUT      : point
              n_points
              poly_points
@OUTPUT     : closest_point
@RETURNS    : distance
@DESCRIPTION: Returns the closest distance of a point to a polygon.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  Real  find_point_polygon_distance(
    Point     *point,
    int       n_points,
    Point     poly_points[],
    Point     *closest_point )
{
    return( sqrt( find_point_polygon_distance_sq( point, n_points, poly_points,
                                                  closest_point ) ) );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : find_closest_polygon_point
@INPUT      : point
              polygons
@OUTPUT     : closest_point
@RETURNS    : index of polygon
@DESCRIPTION: Finds the closest point on a polygons struct, returning the
              index of the polygon which has the closest point.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  int  find_closest_polygon_point(
    Point              *point,
    polygons_struct    *polygons,
    Point              *closest_point )
{
    int            poly, size, closest_poly;
    Real           dist, closest_dist;
    Point          poly_points[MAX_POINTS_PER_POLYGON], closest, best;
    object_struct  object;

    closest_dist = 0.0;

    if( polygons->bintree != NULL )
    {
        object.object_type = POLYGONS;
        object.specific.polygons = *polygons;
        (void) find_closest_point_on_object( point, &object, &closest_poly,
                                             &best );
    }
    else
    {
        for_less( poly, 0, polygons->n_items )
        {
            size = get_polygon_points( polygons, poly, poly_points );
            dist = find_point_polygon_distance_sq( point, size, poly_points,
                                                   &closest);
            if( poly == 0 || dist < closest_dist )
            {
                closest_poly = poly;
                closest_dist = dist;
                best = closest;
            }
        }
    }

    *closest_point = best;

    return( closest_poly );
}
