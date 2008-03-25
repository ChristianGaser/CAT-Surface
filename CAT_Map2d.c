/* Christian Gaser - christian.gaser@uni-jena.de                             */
/* Department of Psychiatry                                                  */
/* University of Jena                                                        */
/*                                                                           */
/* Copyright Christian Gaser, University of Jena.                            */

#include <volume_io/internal_volume_io.h>

Real  compute_clockwise_rotation2( Real x, Real y )
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

public  void  point_to_uv(
    Point            *point,
    Real             *u,
    Real             *v )
{
    Real   x, y, z, theta, phi;

    x = (Real) Point_x(*point);
    y = (Real) Point_y(*point);
    z = (Real) Point_z(*point);

    phi = acos( z );
    theta = compute_clockwise_rotation2( y, x );
    *u = theta / (PI * 2.0);
    *v = phi / PI;
}

public  Real point_to_uv_radius(
    Point            *point,
    Real             *u,
    Real             *v )
{
    Real   x, y, z, theta, phi, scale;

    scale = 100.0;
    
    x = (Real) Point_x(*point)/scale;
    y = (Real) Point_y(*point)/scale;
    z = (Real) Point_z(*point)/scale;

    phi = acos( z );
    theta = compute_clockwise_rotation( y, x );
    *u = theta / (PI * 2.0);
    *v = phi / PI;
    return(sqrt(x*x + y*y + z*z));
}

public void  uv_to_point(
    Real             u,
    Real             v,
    Point            *point )
{
    Real   x, y, z, theta, phi, cos_u;
    Real   sin_u, cos_v, sin_v;

    /* shift theta by 90¡ to obtain correct position of midline */    
    theta = u * PI * 2.0 + PI/2.0;
    phi =   v * PI;

    cos_u = cos( theta );
    sin_u = sin( theta );
    cos_v = cos( phi );
    sin_v = sin( phi );

    z = cos_v;
    x = sin_v * cos_u;
    y = sin_v * sin_u;

    fill_Point( *point, x, y, z );
}
