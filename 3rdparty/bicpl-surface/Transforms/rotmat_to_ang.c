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

/* ----------------------------- MNI Header ----------------------------------- @NAME       : rotmat_to_ang.c
@INPUT      : rot      - rotation matrix (3x3 in num recipes form) calculated
                         by the calling program.
@OUTPUT     : ang      - vector giving rx,ry and rz rotation angles (in 
                         radians). This vector must be defined by the 
                         calling routine.
@RETURNS    : TRUE if ok, FALSE if error.
@DESCRIPTION: this routine extracts the rotation angles from the rotation
              matrix.  The rotation matrix is assumed to be a 3x3 matrix in
              numerical recipes form.  It is locally copied into a 
              4x4 homogeneous matrix for manipulation.

              we assume that the matrix rotation center is (0,0,0).
              we assume that the application of this matrix to a vector
                 is done with rot_mat*vec = premultiplication by matrix

              the resulting angles rx=ang[1],ry=ang[2],rz=ang[3], follow
              the following assumptions:

              -rotations are between -PI to PI
              -rotations are applied 1 - rx, 2 - ry and 3 - rz
              -applying these rotations to an identity matrix will
               result in a matrix equal to `rot' (the input matrix)
              -positive rotations are counter-clockwise when looking
               down the axis, from the positive end towards the origin.
              -I assume a coordinate system:
                          ^ Y
                          |
                          |
                          |
                          |
                          +---------> X
                         /
                        /
                       / Z  (towards the viewer).

              -The problem is posed as:  
                 given a rotation matrix ROT, determine the rotations
                 rx,ry,rz applied in order to give ROT
               solution:
                 assume the rot matrix is equivalent to a normalized
                 orthogonal local coord sys.  i.e.  row 1 of ROT is
                 the local x direction, row 2 is the local y and row 3
                 is the local z.
             
                 (note local is lower case, world is UPPER)

                 1- find RZ that brings local x into the XZ plane, on + size X
                 2- find RY that brings local x*RZ onto X axis
                 3- find RX that brings local z*RZ*RY onto Z axis

                 the required rotations are -RX,-RY and -RZ!

@GLOBALS    : 
@CALLS      : 
@CREATED    : Feb 9, 1992 lc
@MODIFIED   : 
Tue Jun  8 08:44:59 EST 1993 LC
   changes all vec*matrix to matrix*vec.  Std is premultiplication by matrix!
@MODIFIED   : July    4, 1995 D. MacDonald - removed recipes-style code,
                               rewrote completely to handle angles from
                               -PI to PI, instead of just -PI/2 to -PI/2
@MODIFIED   : July   19, 1996 D. MacDonald - now handles left-handed coordinate
                                             systems properly
---------------------------------------------------------------------------- */

#include "bicpl_internal.h"

#ifndef lint
static char rcsid[] = "$Header: /private-cvsroot/libraries/bicpl/Transforms/rotmat_to_ang.c,v 1.21 2005/08/17 22:26:48 bert Exp $";
#endif

#ifdef DEBUG
static  void  are_rotations_equivalent(
    Transform  *rot_trans,
    Real       rx,
    Real       ry,
    Real       rz );
#endif

BICAPI  BOOLEAN  rotmat_to_ang(
    Transform   *rot_trans,
    Real        *ang )
{
    Real       rx, ry, rz, vx, vy, vz, d;
    Vector     x_axis, z_axis, y_axis, cross;
    Transform  z_rot, y_rot;

    /*--- get the x and z axis of the transform, these are all we need
          to compute the three angles, we simply compute the rz, ry, rx
          rotation needed to line up these on (1,0,0) and (0,0,1).  The
          matrix is therefore equivalent to the inverse of this, which
          is the concatenation, in order of the rotations -rx, -ry, -rz */

    get_transform_x_axis( rot_trans, &x_axis );
    get_transform_y_axis( rot_trans, &y_axis );
    get_transform_z_axis( rot_trans, &z_axis );

    /*--- check if it is a left-handed coordinate system, if so
          switch it to right-handed */
 
    CROSS_VECTORS( cross, x_axis, y_axis );
    d = DOT_VECTORS( cross, z_axis );

    if( d < 0.0 )
    {
        print( "rotmat_to_ang: warning, input transform is left-handed.\n" );
        SCALE_VECTOR( x_axis, x_axis, -1.0 );
    }
    else if( d == 0.0 )
    {
        print_error( "rotmat_to_ang: singular system passed in.\n" );
        return( FALSE );
    }
   
    /*--- step one,  find the RZ rotation reqd to bring 
                     the local x into the world XZ plane with a positive X */

    rz = compute_clockwise_rotation( (Real) Vector_x(x_axis),
                                     (Real) Vector_y(x_axis) );

    if( rz >= PI )
        rz -= 2.0 * PI;

    /*--- step two:  find the RY rotation reqd to align 
                     the local x on the world X axis  */

    make_rotation_transform( -rz, Z, &z_rot );

    transform_vector( &z_rot, (Real) Vector_x(x_axis), (Real) Vector_y(x_axis),
                      (Real) Vector_z(x_axis), &vx, &vy, &vz );

    ry = - compute_clockwise_rotation( vx, vz );

    if( ry <= -PI )
        ry += 2.0 * PI;

    /*--- step three, rotate around RX to align the local z with Z */

    make_rotation_transform( -ry, Y, &y_rot );

    transform_vector( &z_rot, (Real) Vector_x(z_axis), (Real) Vector_y(z_axis),
                      (Real) Vector_z(z_axis), &vx, &vy, &vz );
    transform_vector( &y_rot, vx, vy, vz, &vx, &vy, &vz );

    rx = - compute_clockwise_rotation( vz, vy );

    if( rx <= -PI )
        rx += 2.0 * PI;

    /*--- the actual rotations to make up the transform are the negatives
          of these */

    rx = -rx;  /* these are the required rotations */
    ry = -ry;
    rz = -rz;

#ifdef DEBUG
    {
        Transform  test;

        make_identity_transform( &test );
        set_transform_x_axis( &test, &x_axis );
        set_transform_y_axis( &test, &y_axis );
        set_transform_z_axis( &test, &z_axis );

        are_rotations_equivalent( &test, rx, ry, rz );
    }
#endif

    ang[0] = rx;
    ang[1] = ry;
    ang[2] = rz;

    return(TRUE);
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : are_rotations_equivalent
@INPUT      : rot_trans
              rx
              ry
              rz
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: Tests if the rot_transform is equivalent to the three
              successive counter-clockwise rotations.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

#ifdef DEBUG
static  void  are_rotations_equivalent(
    Transform  *rot_trans,
    Real       rx,
    Real       ry,
    Real       rz )
{
    BOOLEAN  error;
    Transform  test, Rx, Ry, Rz;

    /*--- use negatives, since make_rotation_transform assumes clockwise */

    make_rotation_transform( -rx, X, &Rx ) ;
    make_rotation_transform( -ry, Y, &Ry ) ;
    make_rotation_transform( -rz, Z, &Rz ) ;

    concat_transforms( &test, &Rx, &Ry );
    concat_transforms( &test, &test, &Rz );


    error = FALSE;
    for_less( m, 0, N_DIMENSIONS )
        for_less( n, 0, N_DIMENSIONS )
            if( !numerically_close( Transform_elem(test,m,n), 
                                    Transform_elem(*rot_trans,m,n), 1.0e-3 ) )
            {
                error = TRUE;
            }

    if( error )
    {
        for_less( m, 0, N_DIMENSIONS )
        {
            for_less( n, 0, N_DIMENSIONS )
                print( " %g", Transform_elem(test,m,n) );
            print( "\n" );
        }
        print( "\n" );
        for_less( m, 0, N_DIMENSIONS )
        {
            for_less( n, 0, N_DIMENSIONS )
                print( " %g", Transform_elem(*rot_trans,m,n) );
            print( "\n" );
        }

        handle_internal_error( "rotmat_to_ang" );
    }
}
#endif
