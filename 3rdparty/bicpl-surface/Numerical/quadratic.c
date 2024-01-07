#include "bicpl_internal.h"

BICAPI  void  initialize_quadratic(
    int              n_parameters,
    Real             *constant_term,
    float            *linear_terms[],
    float            *square_terms[],
    int              *n_cross_terms[],
    int              **cross_parms[],
    float            **cross_terms[] )
{
    int   parm;

    ALLOC( *linear_terms, n_parameters );
    ALLOC( *square_terms, n_parameters );
    ALLOC( *n_cross_terms, n_parameters );
    ALLOC( *cross_parms, n_parameters );
    ALLOC( *cross_terms, n_parameters );

    *constant_term = 0.0;

    for_less( parm, 0, n_parameters )
    {
        (*linear_terms)[parm] = 0.0f;
        (*square_terms)[parm] = 0.0f;
        (*n_cross_terms)[parm] = 0;
    }
}

BICAPI  void  zero_quadratic(
    int              n_parameters,
    Real             *constant_term,
    float            linear_terms[],
    float            square_terms[],
    int              n_cross_terms[],
    int              *cross_parms[],
    float            *cross_terms[] )
{
    int   parm, n;

    *constant_term = 0.0;

    for_less( parm, 0, n_parameters )
    {
        linear_terms[parm] = 0.0f;
        square_terms[parm] = 0.0f;
        for_less( n, 0, n_cross_terms[parm] )
            cross_terms[parm][n] = 0.0f;
    }
}

BICAPI  void  add_to_quadratic_cross_term(
    int              *n_cross_terms[],
    int              **cross_parms[],
    float            **cross_terms[],
    int              parm1,
    int              parm2,
    Real             value,
    int              alloc_increment )
{
    int   p1, p2, t;

    p1 = MIN( parm1, parm2 );
    p2 = MAX( parm1, parm2 );

    for_less( t, 0, (*n_cross_terms)[p1] )
    {
        if( (*cross_parms)[p1][t] == p2 )
            break;
    }

    if( t >= (*n_cross_terms)[p1] )
    {
        t = (*n_cross_terms)[p1];
        SET_ARRAY_SIZE( (*cross_terms)[p1], t, t+1, alloc_increment );
        SET_ARRAY_SIZE( (*cross_parms)[p1], t, t+1, alloc_increment );
        (*cross_parms)[p1][t] = p2;
        (*cross_terms)[p1][t] = 0.0f;
        ++(*n_cross_terms)[p1];
    }

    (*cross_terms)[p1][t] += (float) value;
}

BICAPI  void  realloc_quadratic_cross_terms(
    int              n_parameters,
    int              n_cross_terms[],
    int              **cross_parms[],
    float            **cross_terms[] )
{
    int    p;

    for_less( p, 0, n_parameters )
    {
        if( n_cross_terms[p] > 0 )
        {
            REALLOC( (*cross_terms)[p], n_cross_terms[p] );
            REALLOC( (*cross_parms)[p], n_cross_terms[p] );
        }
    }
}

BICAPI  void  delete_quadratic(
    int              n_parameters,
    float            linear_terms[],
    float            square_terms[],
    int              n_cross_terms[],
    int              *cross_parms[],
    float            *cross_terms[] )
{
    int  p;

    for_less( p, 0, n_parameters )
    {
        if( n_cross_terms[p] > 0 )
        {
            FREE( cross_parms[p] );
            FREE( cross_terms[p] );
        }
    }

    FREE( linear_terms );
    FREE( square_terms );
    FREE( n_cross_terms );
    FREE( cross_terms );
    FREE( cross_parms );
}

BICAPI  Real   evaluate_quadratic(
    int           n_parameters,
    float         parameters[],
    Real          constant,
    float         linear[],
    float         square[],
    int           n_cross_terms[],
    int           *cross_parms[],
    float         *cross_terms[] )
{
    int    parm, c;
    Real   fit, x, inc;

    fit = constant;

    for_less( parm, 0, n_parameters )
    {
        x = (Real) parameters[parm];

        inc = (Real) linear[parm] + x * (Real) square[parm];

        for_less( c, 0, n_cross_terms[parm] )
            inc += (Real) parameters[cross_parms[parm][c]] *
                   (Real) cross_terms[parm][c];

        fit += x * inc;
    }

    return( fit );
}

BICAPI  void   evaluate_quadratic_deriv(
    int           n_parameters,
    float         parameters[],
    float         linear[],
    float         square[],
    int           n_cross_terms[],
    int           *cross_parms[],
    float         *cross_terms[],
    float         deriv[] )
{
    int    parm, c, n_index;
    float  tx, dx, x;

    for_less( parm, 0, n_parameters )
    {
        x = parameters[parm];

        dx = linear[parm] + 2.0f * x * square[parm];

        for_less( c, 0, n_cross_terms[parm] )
        {
            n_index = cross_parms[parm][c];
            tx = cross_terms[parm][c];

            deriv[n_index] += x * tx;
            dx += parameters[n_index] * tx;
        }

        deriv[parm] += dx;
    }
}
