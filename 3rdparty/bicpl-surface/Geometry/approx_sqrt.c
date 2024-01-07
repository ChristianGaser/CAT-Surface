#include "bicpl_internal.h"

#undef   USE_FAST_APPROX
#define  USE_FAST_APPROX

#ifdef  USE_FAST_APPROX
#define  FAST_SQRT  fast_approx_sqrt
#else
#define  FAST_SQRT  sqrt
#endif

#define  MIN_SQRT_ARG              1.0
#define  MAX_SQRT_ARG              1000.0
#define  N_SQRT_PRECOMPUTE         30000

BICAPI  double  fast_approx_sqrt(
    double  y )
{
    Real     x;
    static   BOOLEAN  first = TRUE;
    static   struct {  double  scale; double trans; }
                                 precomp_sqrt[N_SQRT_PRECOMPUTE], *lookup;

    if( y <= MIN_SQRT_ARG || y >= MAX_SQRT_ARG )
    {
        return( sqrt( y ) );
    }

    if( first )
    {
        int    i;
        Real   next_sqrt, this_sqrt;

        first = FALSE;

        next_sqrt = 0.0;
        for_less( i, 0, N_SQRT_PRECOMPUTE )
        {
            this_sqrt = next_sqrt;
            next_sqrt = sqrt( (Real) (i+1) / (Real) N_SQRT_PRECOMPUTE *
                              MAX_SQRT_ARG );

            precomp_sqrt[i].scale = next_sqrt - this_sqrt;
            precomp_sqrt[i].trans = this_sqrt -
                                    precomp_sqrt[i].scale * (Real) i;
        }
    }

    x = ((Real) N_SQRT_PRECOMPUTE / MAX_SQRT_ARG) * y;

    lookup = &precomp_sqrt[(int) x];

    return( lookup->trans + lookup->scale * x );
}

BICAPI  Real  fast_approx_distance_between_points(
    Point  *p1,
    Point  *p2 )
{
    Real   dx, dy, dz;

    dx = (Real) Point_x(*p2) - (Real) Point_x(*p1);
    dy = (Real) Point_y(*p2) - (Real) Point_y(*p1);
    dz = (Real) Point_z(*p2) - (Real) Point_z(*p1);

    return( FAST_SQRT( dx * dx + dy * dy + dz * dz ) );
}
