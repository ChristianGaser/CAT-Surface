
/* ----------------------------- MNI Header -----------------------------------
@NAME       : 
@INPUT      : 
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: 
@METHOD     : 
@GLOBALS    : 
@CALLS      :  
@CREATED    : Sep. 10, 1996    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */


#ifndef lint
static char rcsid[] = "$Header: /private-cvsroot/libraries/bicpl/Geometry/flatten.c,v 1.6 2005/08/17 22:30:25 bert Exp $";
#endif

#include "bicpl_internal.h"

BICAPI  void  flatten_around_vertex(
    Point     *vertex,
    int       n_neighbours,
    Point     neighbours[],
    BOOLEAN   closed_flag,
    Real      x_flat[],
    Real      y_flat[] )
{
    Real  sum_angles, factor, angle, dist;
    int   i;

    sum_angles = 0.0;

    for_less( i, 0, n_neighbours )
    {
        if( closed_flag || i < n_neighbours-1 )
        {
            angle = get_angle_between_points( &neighbours[i], vertex,
                                              &neighbours[(i+1)%n_neighbours] );
            sum_angles += angle;
        }
    }

    if( sum_angles == 0.0 )
    {
        print_error( "flatten_around_vertex: sum of angles is 0.\n" );
        return;
    }

    if( !closed_flag && sum_angles < 2.0 * PI )
        factor = 1.0;
    else
        factor = 2.0 * PI / sum_angles;

    sum_angles = 0.0;
    for_less( i, 0, n_neighbours )
    {
        angle = sum_angles * factor;

        dist = distance_between_points( vertex, &neighbours[i] );

        x_flat[i] = dist * cos( angle );
        y_flat[i] = dist * sin( angle );

        angle = get_angle_between_points( &neighbours[i], vertex,
                                          &neighbours[(i+1)%n_neighbours] );
        sum_angles += angle;
    }
}

BICAPI  void  flatten_around_vertex_to_sphere(
    Real      radius,
    Point     *vertex,
    int       n_neighbours,
    Point     neighbours[],
    Real      x_sphere[],
    Real      y_sphere[],
    Real      z_sphere[] )
{
    int        n;
    Real       *x_flat, *y_flat, angle, x, y, z, dist;
    Vector     axis;
    Transform  transform;

    ALLOC( x_flat, n_neighbours );
    ALLOC( y_flat, n_neighbours );

    flatten_around_vertex( vertex, n_neighbours, neighbours, TRUE,
                           x_flat, y_flat );

    for_less( n, 0, n_neighbours )
    {
        dist = sqrt( x_flat[n] * x_flat[n] + y_flat[n] * y_flat[n] );
        if( dist / 2.0 > radius )
            dist = 2.0 * radius;
        angle = -2.0 * asin( dist / 2.0 / radius );
        fill_Vector( axis, y_flat[n], -x_flat[n], 0.0 );
        NORMALIZE_VECTOR( axis, axis );
        if( null_Vector( &axis ) )
        {
            x_sphere[n] = 0.0;
            y_sphere[n] = 0.0;
            z_sphere[n] = 0.0;
        }
        else
        {
            make_rotation_about_axis( &axis, angle, &transform );
            transform_point( &transform, 0.0, 0.0, radius, &x, &y, &z );
            z -= radius;
            x_sphere[n] = x;
            y_sphere[n] = y;
            z_sphere[n] = z;
        }
    }

    FREE( x_flat );
    FREE( y_flat );
}
