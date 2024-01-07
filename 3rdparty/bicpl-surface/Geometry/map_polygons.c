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
static char rcsid[] = "$Header: /private-cvsroot/libraries/bicpl/Geometry/map_polygons.c,v 1.14 2005/08/17 22:30:25 bert Exp $";
#endif

/* ----------------------------- MNI Header -----------------------------------
@NAME       : get_baricentric
@INPUT      : point
              p1
              p2
              p3
@OUTPUT     : 
@RETURNS    : baricentric coordinate
@DESCRIPTION: Returns the baricentric coordinate of the point within the
              triangle (p1,p2,p3), relative to the p3 point.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

static  Real  get_baricentric(
    Point       *point,
    Point       *p1,
    Point       *p2,
    Point       *p3 )
{
    Real     weight;
    Vector   hor, up, normal, vert, point_offset;

    SUB_POINTS( point_offset, *point, *p1 );
    SUB_POINTS( hor, *p2, *p1 );
    SUB_POINTS( up, *p3, *p1 );

    CROSS_VECTORS( normal, hor, up );
    CROSS_VECTORS( vert, normal, hor );

    weight = DOT_VECTORS( point_offset, vert ) / DOT_VECTORS( up, vert );

    return( weight );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : get_triangle_interpolation_weights
@INPUT      : point
              points
@OUTPUT     : weights
@RETURNS    : 
@DESCRIPTION: Computes the triangle interpolation weights, meaning the weights
              corresponding to each vertex.  If the values at the vertex are
              v[0], v[1], v[2], they can be interpolated at point as
              value at point = v[0] * weights[0] + v[1] * weights[1] +
                               v[2] * weights[2].
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

static  void  get_triangle_interpolation_weights(
    Point       *point,
    Point       points[],
    Real        weights[] )
{
    weights[0] = get_baricentric( point, &points[1], &points[2], &points[0] );
    weights[1] = get_baricentric( point, &points[2], &points[0], &points[1] );
    weights[2] = get_baricentric( point, &points[0], &points[1], &points[2] );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : intersect_lines_2d
@INPUT      : p1
              p2
              q1
              q2
@OUTPUT     : intersect
@RETURNS    : TRUE if lines intersect
@DESCRIPTION: Tests if two 2d lines intersect, i.e., are not parallel, and
              passes back the point of intersection.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

static  BOOLEAN  intersect_lines_2d(
    Real    p1[],
    Real    p2[],
    Real    q1[],
    Real    q2[],
    Real    intersect[] )
{
    BOOLEAN  intersects;
    Real     t, bottom;
    Real     dp[2], dq[2];

    dp[0] = p2[0] - p1[0];
    dp[1] = p2[1] - p1[1];
    dq[0] = q2[0] - q1[0];
    dq[1] = q2[1] - q1[1];

    bottom = dp[0] * dq[1] - dp[1] * dq[0];

    if( bottom == 0.0 )
    {
        intersects = FALSE;
    }
    else
    {
        t = (dq[0] * (p1[1] - q1[1]) + dq[1] * (q1[0] - p1[0])) / bottom;
        intersect[0] = p1[0] + t * (dp[0] - p1[0]);
        intersect[1] = p1[1] + t * (dp[1] - p1[1]);
        intersects = TRUE;
    }

    return( intersects );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : get_two_d_coordinate
@INPUT      : p
              p1
              p2
              q1
              q2
@OUTPUT     : 
@RETURNS    : two coordinate
@DESCRIPTION: Returns a coordinate distance between 0 and 1 for a point in
              a quadrilateral, whose two opposite sides are p1--p2 and q1--q2.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

static  Real  get_two_d_coordinate(
    Point    *p,
    Point    *p1,
    Point    *p2,
    Point    *q1,
    Point    *q2 )
{
    Real     coords[4][2], coord[2], intersect[2];
    Real     intersect_point[2], len2, dx, dy, idx, idy, factor;
    Vector   cross, offset, hor, vert, v1, v2;
    int      i, j;
    Point    points[4];
    BOOLEAN  found;

    points[0] = *p1;
    points[1] = *p2;
    points[2] = *q2;
    points[3] = *q1;

    found = FALSE;
    for_less( i, 0, 3 )
    {
        SUB_POINTS( v1, points[i], points[0] );
        for_less( j, i+1, 4 )
        {
            SUB_POINTS( v2, points[j], points[0] );
            CROSS_VECTORS( cross, v1, v2 );
            if( !null_Vector(&cross) )
            {
                found = TRUE;
                break;
            }
        }

        if( found )
            break;
    }

    if( !found )
        print_error( "get_two_d_coordinate: degenerate polygons\n" );

    create_two_orthogonal_vectors( &cross, &hor, &vert );

    for_less( i, 0, 4 )
    {
        SUB_POINTS( offset, points[i], points[0] );
        coords[i][0] = DOT_VECTORS( hor, offset );
        coords[i][1] = DOT_VECTORS( vert, offset );
    }

    SUB_POINTS( offset, *p, points[0] );
    coord[0] = DOT_VECTORS( hor, offset );
    coord[1] = DOT_VECTORS( vert, offset );

    if( intersect_lines_2d( coords[0], coords[1], coords[3], coords[2],
                            intersect ) )
    {
        if( !intersect_lines_2d( intersect, coord, coords[0], coords[3],
                                 intersect_point ) )
            return( 0.0 );
    }
    else
    {
        intersect_point[0] = coord[0];
        intersect_point[1] = coord[1];
    }

    dx = coords[3][0] - coords[0][0];
    dy = coords[3][1] - coords[0][1];

    len2 = dx * dx + dy * dy;

    idx = intersect_point[0] - coords[0][0];
    idy = intersect_point[1] - coords[0][1];

    factor = (dx * idx + dy * idy) / len2;

    return( factor );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : get_quadrilateral_interpolation_weights
@INPUT      : point
              points  - 4 points
@OUTPUT     : weights
@RETURNS    : 
@DESCRIPTION: Gets the interpolation weights of the 4 quadrilateral vertices,
              by getting the equivalent u and v components on a unit square.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

static  void  get_quadrilateral_interpolation_weights(
    Point       *point,
    Point       points[],
    Real        weights[] )
{
    Real     u, v;

    u = get_two_d_coordinate( point,
                              &points[0], &points[3], &points[1], &points[2] );
    v = get_two_d_coordinate( point,
                              &points[0], &points[1], &points[3], &points[2] );

    weights[0] = (1.0 - u) * (1.0 - v);
    weights[1] =        u  * (1.0 - v);
    weights[2] =        u  *        v;
    weights[3] = (1.0 - u) *        v;
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : get_arbitrary_polygon_interpolation_weights
@INPUT      : point
              n_points
@OUTPUT     : weights
@RETURNS    : 
@DESCRIPTION: Computes the interpolation weights of each vertex of a polygon.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

static  void  get_arbitrary_polygon_interpolation_weights(
    Point       *point,
    int         n_points,
    Point       points[],
    Real        weights[] )
{
    int      i, j;
    Real     sum_weights, alpha, dist;

    for_less( i, 0, n_points )
        weights[i] = 0.0;

    sum_weights = 0.0;

    for_less( i, 0, n_points )
    {
        dist = get_distance_to_line_segment( point,
                                           &points[i], &points[(i+1)%n_points],
                                           &alpha );

        if( dist == 0.0 )
        {
            for_less( j, 0, n_points )
                weights[j] = 0.0;
            weights[i] = 1.0 - alpha;
            weights[(i+1)%n_points] = alpha;
            sum_weights = 1.0;
            break;
        }

        weights[i] += (1.0 - alpha) / dist;
        weights[(i+1)%n_points] += alpha / dist;

        sum_weights += 1.0 / dist;
    }

    for_less( i, 0, n_points )
        weights[i] /= sum_weights;
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : get_polygon_interpolation_weights
@INPUT      : point
              n_points
              points
@OUTPUT     : weights
@RETURNS    : 
@DESCRIPTION: Computes the interpolation weights for the polygon vertices.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  void  get_polygon_interpolation_weights(
    Point       *point,
    int         n_points,
    Point       points[],
    Real        weights[] )
{
    if( n_points == 3 )
    {
        get_triangle_interpolation_weights( point, points, weights );
    }
    else if( n_points == 4 )
    {
        get_quadrilateral_interpolation_weights( point, points, weights );
    }
    else
    {
        get_arbitrary_polygon_interpolation_weights( point, n_points,
                                                     points, weights );
    }
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : map_point_between_polygons
@INPUT      : p1           - polygons struct
              poly_index   - which polygon
              p1_point     - a point in the polygon
              p2           - the destination polygons of same topology
@OUTPUT     : p2_point     - the corresponding point in the second polygon
@RETURNS    : distance of p1_point from the polygons
@DESCRIPTION: Given two polygons of the same topology and a point in the first
              polygon, finds the corresponding point in the second polygon.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  void  map_point_between_polygons(
    polygons_struct  *p1,
    int              poly_index,
    Point            *p1_point,
    polygons_struct  *p2,
    Point            *p2_point )
{
    int     i, size;
    Point   poly1_points[MAX_POINTS_PER_POLYGON];
    Point   poly2_points[MAX_POINTS_PER_POLYGON];
    Point   scaled_point;
    Real    weights[MAX_POINTS_PER_POLYGON];

    size = get_polygon_points( p1, poly_index, poly1_points );

    get_polygon_interpolation_weights( p1_point, size, poly1_points, weights );

    if( get_polygon_points( p2, poly_index, poly2_points ) != size )
        handle_internal_error( "map_point_between_polygons" );

    fill_Point( *p2_point, 0.0, 0.0, 0.0 );

    for_less( i, 0, size )
    {
        SCALE_POINT( scaled_point, poly2_points[i], weights[i] );
        ADD_POINTS( *p2_point, *p2_point, scaled_point );
    }
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : map_point_to_unit_sphere
@INPUT      : p
              point
              unit_sphere
@OUTPUT     : unit_sphere_point
@RETURNS    : distance of point from polygons
@DESCRIPTION: Maps a point on a polygon to the unit sphere of the same topology.
              Returns the distance of the point from the polygons.  Usually,
              we are mapping points on the polygon to points on the sphere,
              so the distance will be 0.  However, we may want map points near
              the polygons, so we will return the distance so that the
              calling program may use the distance to determine if mapping this
              point makes sense.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  Real  map_point_to_unit_sphere(
    polygons_struct  *p,
    Point            *point,
    polygons_struct  *unit_sphere,
    Point            *unit_sphere_point )
{
    Real     mag, dist;
    Point    poly_point;
    Vector   offset;
    int      poly;

    poly = find_closest_polygon_point( point, p, &poly_point );

    dist = distance_between_points( point, &poly_point );

    map_point_between_polygons( p, poly, &poly_point, unit_sphere,
                                unit_sphere_point );

    CONVERT_POINT_TO_VECTOR( offset, *unit_sphere_point );

    mag = MAGNITUDE( offset );

    /*--- project it to the sphere */

    if( mag != 1.0 )
    {
        SCALE_POINT( *unit_sphere_point, *unit_sphere_point, 1.0 / mag );
    }

    return( dist );
}

BICAPI  void  map_unit_sphere_to_point(
    polygons_struct  *unit_sphere,
    Point            *unit_sphere_point,
    polygons_struct  *p,
    Point            *point )
{
    Point    poly_point;
    int      poly;

    poly = find_closest_polygon_point( unit_sphere_point, unit_sphere,
                                       &poly_point );

    map_point_between_polygons( unit_sphere, poly, &poly_point, p,
                                point );
}


/* ----------------------------- MNI Header -----------------------------------
@NAME       : polygon_transform_point
@INPUT      : src_object
              dest_object
              src_point
@OUTPUT     : dest_point
@RETURNS    : 
@DESCRIPTION: Finds the point on the destination polygons corresponding to the
              src_point on the source polygons.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

static  void  polygon_transform_point(
    object_struct   *src_object,
    object_struct   *dest_object,
    Point           *src_point,
    Point           *dest_point )
{
    int     obj_index;
    Point   point;

    (void) find_closest_point_on_object( src_point, src_object, &obj_index,
                                         &point );

    map_point_between_polygons( get_polygons_ptr(src_object), obj_index,
                                &point, get_polygons_ptr(dest_object),
                                dest_point );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : polygon_transform_points
@INPUT      : src_polygons
              dest_polygons
              n_points
              src_points
@OUTPUT     : dest_points
@RETURNS    : 
@DESCRIPTION: Maps each of the src_point to the dest_polygons.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  void  polygon_transform_points(
    polygons_struct   *src_polygons,
    polygons_struct   *dest_polygons,
    int               n_points,
    Point             src_points[],
    Point             dest_points[] )
{
    int             i;
    object_struct   *src_object, *dest_object;

    if( !polygons_are_same_topology( src_polygons, dest_polygons ) )
    {
        print_error(
             "polygon_transform_points: polygons are not same topology.\n" );
        return;
    }

    src_object = create_object( POLYGONS );
    *get_polygons_ptr(src_object) = *src_polygons;
    dest_object = create_object( POLYGONS );
    *get_polygons_ptr(dest_object) = *dest_polygons;

    for_less( i, 0, n_points )
    {
        polygon_transform_point( src_object, dest_object,
                                 &src_points[i], &dest_points[i] );
    }
}

BICAPI  void  map_sphere_to_uv(
    Real    x,
    Real    y,
    Real    z,
    Real    *u,
    Real    *v )
{
    Real    angle_up, angle_around;

    angle_up = acos( z );

    *v = 1.0 - angle_up / PI;

    angle_around = compute_clockwise_rotation( x, -y );

    *u = angle_around / (2.0*PI);
}

BICAPI  void  map_uv_to_sphere(
    Real    u,
    Real    v,
    Real    *x,
    Real    *y,
    Real    *z )
{
    Real   r;

    *z = cos( (1.0 - v) * PI );
    r = sin( (1.0 - v) * PI );

    *x = r * cos( u * 2.0 * PI );
    *y = r * sin( u * 2.0 * PI );
}
