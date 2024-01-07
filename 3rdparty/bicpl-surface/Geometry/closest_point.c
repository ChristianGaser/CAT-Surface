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
static char rcsid[] = "$Header: /private-cvsroot/libraries/bicpl/Geometry/closest_point.c,v 1.10 2005/08/17 22:30:25 bert Exp $";
#endif

#define  MAX_POINTS    300

/* ----------------------------- MNI Header -----------------------------------
@NAME       : get_point_polygon_distance_sq
@INPUT      : point
              polygons
              poly_index
@OUTPUT     : object_point
@RETURNS    : distance
@DESCRIPTION: Returns the closest distance from the point to the polygon,
              and passes back the closest point in the polygon.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

static  Real   get_point_polygon_distance_sq(
    Point            *point,
    polygons_struct  *polygons,
    int              poly_index,
    Point            *object_point )
{
    Point    points[MAX_POINTS];
    int      size;

    if( polygons->visibilities == NULL ||
        polygons->visibilities[poly_index] )
    {
        size = get_polygon_points( polygons, poly_index, points );

        return( find_point_polygon_distance_sq( point, size, points,
                                                object_point ) );
    }
    else
        return( 1.0e60 );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : get_point_quadmesh_distance
@INPUT      : point
              quadmeshs
              obj_index
@OUTPUT     : object_point
@RETURNS    : distance
@DESCRIPTION: Returns the closest distance from the point to the quadmesh,
              and passes back the closest point in the quadmesh.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

static  Real   get_point_quadmesh_distance_sq(
    Point            *point,
    quadmesh_struct  *quadmesh,
    int              obj_index,
    Point            *object_point )
{
    Point    points[4];
    int      i, j, m, n;

    get_quadmesh_n_objects( quadmesh, &m, &n );

    i = obj_index / n;
    j = obj_index % n;

    get_quadmesh_patch( quadmesh, i, j, points );

    return( find_point_polygon_distance_sq( point, 4, points, object_point ) );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : get_line_segment_alpha
@INPUT      : point
              p1     \ endpoints of segment
              p2     /
@OUTPUT     : closest_point
              alpha
@RETURNS    : 
@DESCRIPTION: Finds the closest point on the line segment to the 'point',
              and the position of the point on the segment, ranging from
              0 to 1.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

static  void  get_line_segment_alpha(
    Point  *point,
    Point  *p1,
    Point  *p2,
    Point  *closest_point,
    Real   *alpha )
{
    Vector   p2_minus_p1, p_minus_p1;
    Real     mag;

    SUB_VECTORS( p2_minus_p1, *p2, *p1 );
    SUB_VECTORS( p_minus_p1, *point, *p1 );

    mag = DOT_VECTORS( p2_minus_p1, p2_minus_p1 );

    if( mag == 0.0 )
        *alpha = 0.0;
    else
        *alpha = DOT_VECTORS( p_minus_p1, p2_minus_p1 ) / mag;

    if( *alpha <= 0.0 )
    {
        *closest_point = *p1;
        *alpha = 0.0;
    }
    else if( *alpha >= 1.0 )
    {
        *closest_point = *p2;
        *alpha = 1.0;
    }
    else
    {
        INTERPOLATE_POINTS( *closest_point, *p1, *p2, *alpha );
    }
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : get_closest_point_on_line_segment
@INPUT      : point
              p1     \ endpoints of segment
              p2     /
@OUTPUT     : closest_point
@RETURNS    : 
@DESCRIPTION: Finds the closest point on the line segment to the 'point'.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  void  get_closest_point_on_line_segment(
    Point  *point,
    Point  *p1,
    Point  *p2,
    Point  *closest_point )
{
    Real     t;

    get_line_segment_alpha( point, p1, p2, closest_point, &t );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : get_distance_to_line_segment
@INPUT      : point
              p1     \ endpoints of segment
              p2     /
@OUTPUT     : alpha
@RETURNS    : Distance
@DESCRIPTION: Finds the distance to the line segment, and returns the
              relative position on the segment of the closest point
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  Real  get_distance_to_line_segment(
    Point  *point,
    Point  *p1,
    Point  *p2,
    Real   *alpha )
{
    Point    closest;

    get_line_segment_alpha( point, p1, p2, &closest, alpha );

    return( distance_between_points( point, &closest ) );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : get_point_line_segment_distance_sq
@INPUT      : point
              lines
              obj_index
@OUTPUT     : object_point
@RETURNS    : distance
@DESCRIPTION: Returns the distance to the closest point on the line segment.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

static  Real   get_point_line_segment_distance_sq(
    Point            *point,
    lines_struct     *lines,
    int              obj_index,
    Point            *object_point )
{
    int      line, seg, p1, p2;

    get_line_segment_index( lines, obj_index, &line, &seg );

    p1 = lines->indices[POINT_INDEX(lines->end_indices,line,seg)];
    p2 = lines->indices[POINT_INDEX(lines->end_indices,line,seg+1)];

    get_closest_point_on_line_segment( point, &lines->points[p1],
                                       &lines->points[p2], object_point );

    return( sq_distance_between_points( point, object_point ) );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : get_point_marker_distance
@INPUT      : point
              marker
@OUTPUT     : object_point
@RETURNS    : distance
@DESCRIPTION: Returns the distance from the point to the marker.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

static  Real  get_point_marker_distance(
    Point            *point,
    marker_struct    *marker,
    Point            *object_point )
{
    *object_point = marker->position;

    return( distance_between_points( point, object_point ) );
}

BICAPI  Real  get_point_object_distance_sq(
    Point                 *point,
    object_struct         *object,
    int                   obj_index,
    Point                 *object_point )
{
    Real                  dist;

    if( get_object_type( object ) == POLYGONS )
    {
        dist = get_point_polygon_distance_sq( point, get_polygons_ptr(object),
                                              obj_index, object_point );
    }
    else if( get_object_type( object ) == QUADMESH )
    {
        dist = get_point_quadmesh_distance_sq( point, get_quadmesh_ptr(object),
                                               obj_index, object_point );
    }
    else if( get_object_type( object ) == LINES )
    {
        dist = get_point_line_segment_distance_sq( point, get_lines_ptr(object),
                                                   obj_index, object_point );
    }
    else if( get_object_type( object ) == MARKER )
    {
        dist = get_point_marker_distance( point, get_marker_ptr(object),
                                          object_point );
        dist *= dist;
    }
    else
        dist = 1.0e60;

    return( dist );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : get_point_object_distance
@INPUT      : point
              object
              obj_index
@OUTPUT     : object_point
@RETURNS    : distance
@DESCRIPTION: Finds the closest distance to the given object.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  Real  get_point_object_distance(
    Point                 *point,
    object_struct         *object,
    int                   obj_index,
    Point                 *object_point )
{
    return( sqrt( get_point_object_distance_sq( point, object,
                                                obj_index, object_point ) ) );
}

static  Real   get_point_polygon_vertex_distance(
    Point            *point,
    polygons_struct  *polygons,
    int              poly_index,
    int              *object_vertex )
{
    Point    points[MAX_POINTS];
    Real     dist, best_dist;
    int      i, size, point_index;

    if( polygons->visibilities == (Smallest_int *) 0 ||
        polygons->visibilities[poly_index] )
    {
        size = get_polygon_points( polygons, poly_index, points );

        best_dist = 0.0;
        for_less( i, 0, size )
        {
            point_index = polygons->indices[
                    POINT_INDEX(polygons->end_indices,poly_index,i)];

            dist = distance_between_points( point,
                                            &polygons->points[point_index] );

            if( i == 0 || dist < best_dist )
            {
                best_dist = dist;
                *object_vertex = point_index;
            }
        }
    }
    else
        best_dist = 1.0e30;

    return( best_dist );
}

BICAPI  Real  get_point_object_vertex_distance(
    Point                 *point,
    object_struct         *object,
    int                   obj_index,
    int                   *object_vertex )
{
    Real                  dist;

    if( get_object_type( object ) == POLYGONS )
    {
        dist = get_point_polygon_vertex_distance( point,
                     get_polygons_ptr(object), obj_index, object_vertex );
    }
    else
    {
        print_error(
             "not implemented type in get_point_object_vertex_distance()n" );
        dist = 1.0e30;
    }

    return( dist );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : find_closest_point_on_object
@INPUT      : point
              object
@OUTPUT     : obj_index
              point_on_object
@RETURNS    : distance
@DESCRIPTION: Finds the closest point on the object, and the index of the
              object component where it was found.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  Real  find_closest_point_on_object(
    Point           *point,
    object_struct   *object,
    int             *obj_index,
    Point           *point_on_object )
{
    Real             closest_dist, dist;
    Point            obj_point;
    lines_struct     *lines;
    polygons_struct  *polygons;
    quadmesh_struct  *quadmesh;
    int              i, m, n, n_objects;

    if( obj_index != (int *) NULL )
        *obj_index = -1;

    switch( get_object_type( object ) )
    {
    case LINES:
        lines = get_lines_ptr( object );
        if( lines->n_items == 0 )
        {
            n_objects = 0;
            break;
        }

        if( lines->bintree != (bintree_struct_ptr) NULL )
        {
            return( find_closest_point_in_bintree( point, lines->bintree,
                             object, obj_index, point_on_object ) );
        }
   
        n_objects = lines->end_indices[lines->n_items-1] - lines->n_items;
        break;

    case POLYGONS:
        polygons = get_polygons_ptr( object );
        if( polygons->bintree != (bintree_struct_ptr) NULL )
        {
            return( find_closest_point_in_bintree( point, polygons->bintree,
                             object, obj_index, point_on_object ) );
        }
        n_objects = polygons->n_items;
        break;

    case QUADMESH:
        quadmesh = get_quadmesh_ptr( object );
        if( quadmesh->bintree != (bintree_struct_ptr) NULL )
        {
            return( find_closest_point_in_bintree( point, quadmesh->bintree,
                             object, obj_index, point_on_object ) );
        }
        get_quadmesh_n_objects( quadmesh, &m, &n );
        n_objects = m * n;
        break;

    case MARKER:
        n_objects = 1;
        break;

    default:
        n_objects = 0;
        break;
    }

    closest_dist = 0.0;
    for_less( i, 0, n_objects )
    {
        dist = get_point_object_distance( point, object, i, &obj_point);

        if( i == 0 || dist < closest_dist )
        {
            closest_dist = dist;
            *obj_index = i;
            *point_on_object = obj_point;
        }
    }

    return( closest_dist );
}

BICAPI  Real  find_closest_vertex_on_object(
    Point           *point,
    object_struct   *object,
    int             *vertex_on_object )
{
    Real             closest_dist, dist;
    Point            *points;
    int              i, n_points;

    if( get_object_type( object ) == POLYGONS &&
        get_polygons_ptr( object )->bintree != (bintree_struct_ptr) NULL )
    {
        closest_dist = find_closest_vertex_in_bintree( point,
                get_polygons_ptr( object )->bintree,
                object, vertex_on_object );
    }
    else
    {
        n_points = get_object_points( object, &points );
        closest_dist = 0.0;
        for_less( i, 0, n_points )
        {
            dist = distance_between_points( point, &points[i] );

            if( i == 0 || dist < closest_dist )
            {
                closest_dist = dist;
                *vertex_on_object = i;
            }
        }
    }

    return( closest_dist );
}
