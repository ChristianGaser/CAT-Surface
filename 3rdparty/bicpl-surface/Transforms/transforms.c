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
static char rcsid[] = "$Header: /private-cvsroot/libraries/bicpl/Transforms/transforms.c,v 1.14 2005/08/17 22:26:48 bert Exp $";
#endif

/* ----------------------------- MNI Header -----------------------------------
@NAME       : make_scale_transform
@INPUT      : sx
              sy
              sz
@OUTPUT     : transform
@RETURNS    : 
@DESCRIPTION: Assigns the transform a (possibly anisotropic) scaling transform.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : 1993            David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  void  make_scale_transform( 
    Real        sx,
    Real        sy,
    Real        sz,
    Transform   *transform )
{
    make_identity_transform( transform );

    Transform_elem( *transform, 0, 0 ) = sx;
    Transform_elem( *transform, 1, 1 ) = sy;
    Transform_elem( *transform, 2, 2 ) = sz;
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : set_transform_x_and_z_axes
@INPUT      : x_axis
              z_axis
@OUTPUT     : transform
@RETURNS    : 
@DESCRIPTION: Sets the transform to the orthogonal coordinate system defined
              by the x_axis and z_axis.  To do this, the z axis is used as
              the z axis, and the cross product of the z_axis and the
              x_axis is used as the y_axis.  Then, the x_axis is taken as
              the cross product of the y and z axes, as opposed to the
              one passed in to the function.  In this way, as long as the
              passed in x_axis is not collinear with the z_axis, it need
              not be perpendicular to the z_axis.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : 1993            David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  void  set_transform_x_and_z_axes(
    Transform   *transform,
    Vector      *x_axis,
    Vector      *z_axis )
{
    Vector  n_z, n_y, n_x;

    NORMALIZE_VECTOR( n_z, *z_axis );
    CROSS_VECTORS( n_y, n_z, *x_axis );
    NORMALIZE_VECTOR( n_y, n_y );
    CROSS_VECTORS( n_x, n_z, n_y );
    NORMALIZE_VECTOR( n_x, n_x );

    set_transform_x_axis( transform, &n_x );
    set_transform_y_axis( transform, &n_y );
    set_transform_z_axis( transform, &n_z );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : make_translation_transform
@INPUT      : x_trans
              y_trans
              z_trans
@OUTPUT     : transform
@RETURNS    : 
@DESCRIPTION: Creates a translation transform.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : 1993            David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  void  make_translation_transform(
    Real        x_trans,
    Real        y_trans,
    Real        z_trans,
    Transform   *transform )
{
    make_identity_transform( transform );

    Transform_elem( *transform, 0, 3 ) = x_trans;
    Transform_elem( *transform, 1, 3 ) = y_trans;
    Transform_elem( *transform, 2, 3 ) = z_trans;
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : make_origin_transform
@INPUT      : origin
@OUTPUT     : transform
@RETURNS    : 
@DESCRIPTION: Makes a transformation which has a different origin.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : 1993            David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  void  make_origin_transform(
    Point       *origin,
    Transform   *transform )
{
    make_identity_transform( transform );

    Transform_elem( *transform, 0, 3 ) = (Transform_elem_type) Point_x(*origin);
    Transform_elem( *transform, 1, 3 ) = (Transform_elem_type) Point_y(*origin);
    Transform_elem( *transform, 2, 3 ) = (Transform_elem_type) Point_z(*origin);
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : make_rotation_transform
@INPUT      : radians
              axis
@OUTPUT     : transform
@RETURNS    : 
@DESCRIPTION: Makes a transform which rotates the specified number of radians
              around the given axis (X, Y, or Z).
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : 1993            David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  void  make_rotation_transform(
    Real       radians,
    int        axis,
    Transform  *transform )
{
    int   a1, a2;
    Real  c, s;

    a1 = (axis + 1) % N_DIMENSIONS;
    a2 = (axis + 2) % N_DIMENSIONS;

    make_identity_transform( transform );

    c = cos( (double) radians );
    s = sin( (double) radians );

    Transform_elem( *transform, a1, a1 ) = c;
    Transform_elem( *transform, a1, a2 ) = s;
    Transform_elem( *transform, a2, a1 ) = -s;
    Transform_elem( *transform, a2, a2 ) = c;
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : make_transform_relative_to_point
@INPUT      : point
              transform
@OUTPUT     : rel_transform
@RETURNS    : 
@DESCRIPTION: Creates a transform which is the specified transform performed
              relative to the given point.  For instance to rotate around
              a point, one could pass the point, and a rotation transform
              created by make_rotation_transform(),
              and the resulting transform would perform rotation about that
              point.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : 1993            David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  void  make_transform_relative_to_point(
    Point      *point,
    Transform  *transform,
    Transform  *rel_transform )
{
    Transform  to_origin, to_point;

    make_translation_transform( (Real) Point_x(*point), (Real) Point_y(*point),
                                (Real) Point_z(*point), &to_point );

    make_translation_transform( -(Real) Point_x(*point),
                                -(Real) Point_y(*point),
                                -(Real) Point_z(*point), &to_origin );

    concat_transforms( rel_transform, &to_origin, transform );
    concat_transforms( rel_transform, rel_transform, &to_point );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : make_transform_in_coordinate_system
@INPUT      : origin
              x_axis
              y_axis
              z_axis
              transform
@OUTPUT     : rel_transform
@RETURNS    : 
@DESCRIPTION: Creates a transform which is the input transform relative to
              the specified coordinate system.  For instance if the input
              transform is a scaling transform, then the resulting transform
              will scale points about the 'origin', in the directions of
              'x_axis', 'y_axis', and 'z_axis'.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : 1993            David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  void  make_transform_in_coordinate_system(
    Point      *origin,
    Vector     *x_axis,
    Vector     *y_axis,
    Vector     *z_axis,
    Transform  *transform,
    Transform  *rel_transform )
{
    Transform  to_bases, from_bases;

    make_change_to_bases_transform( origin, x_axis, y_axis, z_axis, &to_bases );
    make_change_from_bases_transform( origin, x_axis, y_axis, z_axis,
                                      &from_bases );

    concat_transforms( rel_transform, &from_bases, transform );
    concat_transforms( rel_transform, rel_transform, &to_bases );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : make_rotation_about_axis
@INPUT      : axis
              angle
@OUTPUT     : transform
@RETURNS    : 
@DESCRIPTION: Creates a transform which corresponds to a rotation about an
              arbitrary vector.  Code comes from Graphics Gems book.
@METHOD     : Graphics Gems, Page ??
@GLOBALS    : 
@CALLS      : 
@CREATED    : 1993            David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  void  make_rotation_about_axis(
    Vector     *axis,
    Real       angle,
    Transform  *transform )
{
    Real    c, s, t;
    Real    txy, txz, tyz, sx, sy, sz;
    Real    x, y, z;
    Vector  unit_axis;

    NORMALIZE_VECTOR( unit_axis, *axis );

    c = cos( (double) -angle );
    s = sin( (double) -angle );
    t = 1.0 - c;

    x = (Real) Vector_x( unit_axis );
    y = (Real) Vector_y( unit_axis );
    z = (Real) Vector_z( unit_axis );

    txy = t * x * y;
    txz = t * x * z;
    tyz = t * y * z;

    sx = s * x;
    sy = s * y;
    sz = s * z;

    Transform_elem( *transform, 0, 0 ) = t * x * x + c;
    Transform_elem( *transform, 0, 1 ) = txy + sz;
    Transform_elem( *transform, 0, 2 ) = txz - sy;
    Transform_elem( *transform, 0, 3 ) = 0.0;

    Transform_elem( *transform, 1, 0 ) = txy - sz;
    Transform_elem( *transform, 1, 1 ) = t * y * y + c;
    Transform_elem( *transform, 1, 2 ) = tyz + sx;
    Transform_elem( *transform, 1, 3 ) = 0.0;

    Transform_elem( *transform, 2, 0 ) = txz + sy;
    Transform_elem( *transform, 2, 1 ) = tyz - sx;
    Transform_elem( *transform, 2, 2 ) = t * z * z + c;
    Transform_elem( *transform, 2, 3 ) = 0.0;

    Transform_elem( *transform, 3, 0 ) = 0.0;
    Transform_elem( *transform, 3, 1 ) = 0.0;
    Transform_elem( *transform, 3, 2 ) = 0.0;
    Transform_elem( *transform, 3, 3 ) = 1.0;
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : convert_2d_transform_to_rotation_translation
@INPUT      : transform
@OUTPUT     : degrees_clockwise
              x_trans
              y_trans
@RETURNS    : 
@DESCRIPTION: Converts a 2d transform to a angle of rotation and origin
              of rotation.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : 1993            David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  void  convert_2d_transform_to_rotation_translation(
    Transform  *transform,
    Real       *degrees_clockwise,
    Real       *x_trans,
    Real       *y_trans )
{
    Real   x, y;

    x = Transform_elem(*transform, X, X );
    y = Transform_elem(*transform, Y, X );

    *degrees_clockwise = RAD_TO_DEG * compute_clockwise_rotation( x, y );
    *x_trans = Transform_elem( *transform, X, 3 );
    *y_trans = Transform_elem( *transform, Y, 3 );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : compute_clockwise_rotation
@INPUT      : x
              y
@OUTPUT     : 
@RETURNS    : Radians
@DESCRIPTION: Computes the clockwise angle of rotation required to rotate
              the x axis to the point (x,y).  Will be value greater than or
              equal to 0 and less than two pi.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : 1993            David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  Real  compute_clockwise_rotation( Real x, Real y )
{
    Real   radians;

    if( x == 0.0 )
    {
        if( y < 0.0 )
            return( PI / 2.0 );
        else if( y > 0.0 )
            return( 3.0 * PI / 2.0 );
        else
            return( 0.0 );
    }
    else if( y == 0.0 )
    {
        if( x > 0.0 )
            return( 0.0 );
        else
            return( PI );
    }
    else
    {
        radians = - (Real) atan2( (double) y, (double) x );

        if( radians < 0.0 )
            radians += 2.0 * PI;

        return( radians );
    }
}

BICAPI  BOOLEAN  is_transform_right_handed(
    Transform   *transform )
{
    Real     volume;
    Vector   x_axis, y_axis, z_axis, offset;

    get_transform_x_axis( transform, &x_axis );
    get_transform_y_axis( transform, &y_axis );
    get_transform_z_axis( transform, &z_axis );

    CROSS_VECTORS( offset, x_axis, y_axis );
    volume = DOT_VECTORS( offset, z_axis );

    return( volume > 0.0 );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : make_identity_transform_2d
@INPUT      : transform
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: Fills in the 2d transform with the identity transform.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : 1993            David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  void  make_identity_transform_2d(
    Transform_2d *transform )
{
    Transform_2d_elem( *transform, 0, 0 ) = 1.0;
    Transform_2d_elem( *transform, 0, 1 ) = 0.0;
    Transform_2d_elem( *transform, 0, 2 ) = 0.0;
    Transform_2d_elem( *transform, 1, 0 ) = 0.0;
    Transform_2d_elem( *transform, 1, 1 ) = 1.0;
    Transform_2d_elem( *transform, 1, 2 ) = 0.0;
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : get_inverse_transform_2d
@INPUT      : transform
@OUTPUT     : inverse
@RETURNS    : 
@DESCRIPTION: Passes back the inverse of the 2d transform.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : 1993            David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  void  get_inverse_transform_2d(
    Transform_2d   *transform,
    Transform_2d   *inverse )
{
    Real  determinant;

    determinant = Transform_2d_elem(*transform,0,0) *
                  Transform_2d_elem(*transform,1,1) -
                  Transform_2d_elem(*transform,1,0) *
                  Transform_2d_elem(*transform,0,1);

    Transform_2d_elem(*inverse,0,0) = Transform_2d_elem(*transform,1,1) /
                                      determinant;
    Transform_2d_elem(*inverse,0,1) = -Transform_2d_elem(*transform,0,1) /
                                      determinant;
    Transform_2d_elem(*inverse,0,2) = (Transform_2d_elem(*transform,1,2) *
                                       Transform_2d_elem(*transform,0,1) -
                                       Transform_2d_elem(*transform,1,1) *
                                       Transform_2d_elem(*transform,0,2)) /
                                      determinant;

    Transform_2d_elem(*inverse,1,0) = -Transform_2d_elem(*transform,1,0) /
                                      determinant;
    Transform_2d_elem(*inverse,1,1) = Transform_2d_elem(*transform,0,0) /
                                      determinant;
    Transform_2d_elem(*inverse,1,2) = (Transform_2d_elem(*transform,1,0) *
                                       Transform_2d_elem(*transform,0,2) -
                                       Transform_2d_elem(*transform,0,0) *
                                       Transform_2d_elem(*transform,1,2)) /
                                      determinant;
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : transform_point_2d
@INPUT      : transform
              x
              y
@OUTPUT     : x_trans
              y_trans
@RETURNS    : 
@DESCRIPTION: Transforms a 2d point by the specified transform.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : 1993            David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  void  transform_point_2d(
    Transform_2d   *transform,
    Real           x,
    Real           y,
    Real           *x_trans,
    Real           *y_trans )
{
    *x_trans = Transform_2d_elem(*transform,0,0) * x +
               Transform_2d_elem(*transform,0,1) * y +
               Transform_2d_elem(*transform,0,2);
    *y_trans = Transform_2d_elem(*transform,1,0) * x +
               Transform_2d_elem(*transform,1,1) * y +
               Transform_2d_elem(*transform,1,2);
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : get_least_squares_transform_2d
@INPUT      : n_points
              x
              y
              x_trans
              y_trans
@OUTPUT     : transform_2d
@RETURNS    : 
@DESCRIPTION: Computes the 2d transform that best transforms the points
              ( x[], y[] ) to ( x_trans[], y_trans[] ), using linear least
              squares.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : 1993            David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  void  get_least_squares_transform_2d(
    int           n_points,
    Real          x[],
    Real          y[],
    Real          x_trans[],
    Real          y_trans[],
    Transform_2d  *transform_2d )
{
    int   p;
    Real  coefs[2+1];
    Real  **coords;

    ALLOC2D( coords, n_points, 2 );

    for_less( p, 0, n_points )
    {
        coords[p][X] = x[p];
        coords[p][Y] = y[p];
    }

    least_squares( n_points, 2, coords, x_trans, coefs );

    Transform_2d_elem( *transform_2d, 0, 0 ) = coefs[1];
    Transform_2d_elem( *transform_2d, 0, 1 ) = coefs[2];
    Transform_2d_elem( *transform_2d, 0, 2 ) = coefs[0];
    
    least_squares( n_points, 2, coords, y_trans, coefs );

    Transform_2d_elem( *transform_2d, 1, 0 ) = coefs[1];
    Transform_2d_elem( *transform_2d, 1, 1 ) = coefs[2];
    Transform_2d_elem( *transform_2d, 1, 2 ) = coefs[0];

    FREE2D( coords );
}
