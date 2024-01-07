#include "bicpl_internal.h"

#define  LSQ_TYPE      Real
#define  MINIMIZE_LSQ  minimize_lsq

#include  "minimize_lsq_include.c"

BICAPI  Real   minimize_lsq(
    int              n_parameters,
    Real             constant_term,
    Real             linear_terms[],
    Real             square_terms[],
    int              n_cross_terms[],
    int              *cross_parms[],
    Real             *cross_terms[],
    Real             max_step_size,
    int              n_iters,
    Real             node_values[] )
{
    return( private_minimize_lsq( n_parameters, constant_term,
                                  linear_terms, square_terms,
                                  n_cross_terms, cross_parms, cross_terms,
                                  max_step_size, n_iters, node_values ) );
}

BICAPI  void  initialize_lsq_terms(
    int              n_parameters,
    Real             *constant_term,
    Real             *linear_terms[],
    Real             *square_terms[],
    int              *n_cross_terms[],
    int              **cross_parms[],
    Real             **cross_terms[] )
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
        (*linear_terms)[parm] = 0.0;
        (*square_terms)[parm] = 0.0;
        (*n_cross_terms)[parm] = 0;
    }
}

BICAPI  void  reset_lsq_terms(
    int              n_parameters,
    Real             *constant_term,
    Real             linear_terms[],
    Real             square_terms[],
    int              n_cross_terms[],
    int              *cross_parms[],
    Real             *cross_terms[] )
{
    int   parm, n;

    *constant_term = 0.0;

    for_less( parm, 0, n_parameters )
    {
        linear_terms[parm] = 0.0;
        square_terms[parm] = 0.0;
        for_less( n, 0, n_cross_terms[parm] )
            cross_terms[parm][n] = 0.0;
    }
}

BICAPI  void  add_to_lsq_terms(
    int              n_parameters,
    Real             *constant_term,
    Real             linear_terms[],
    Real             square_terms[],
    int              n_cross_terms[],
    int              *cross_parms[],
    Real             *cross_terms[],
    int              n_in_list,
    int              list[],
    Real             weights[],
    Real             constant,
    int              alloc_increment )
{
    int   p, q, p1, p2, t;

    *constant_term += constant * constant;

    for_less( p, 0, n_in_list )
    {
        linear_terms[list[p]] += 2.0 * weights[p] * constant;
        square_terms[list[p]] += weights[p] * weights[p];
        for_less( q, p+1, n_in_list )
        {
            p1 = MIN( list[p], list[q] );
            p2 = MAX( list[p], list[q] );
            for_less( t, 0, n_cross_terms[p1] )
            {
                if( cross_parms[p1][t] == p2 )
                    break;
            }

            if( t >= n_cross_terms[p1] )
            {
                t = n_cross_terms[p1];
                SET_ARRAY_SIZE( cross_terms[p1], t, t+1, alloc_increment );
                SET_ARRAY_SIZE( cross_parms[p1], t, t+1, alloc_increment );
                cross_parms[p1][t] = p2;
                cross_terms[p1][t] = 0.0;
                ++n_cross_terms[p1];
            }

            cross_terms[p1][t] += 2.0 * weights[p] * weights[q];
        }
    }
}

BICAPI  void  realloc_lsq_terms(
    int              n_parameters,
    int              n_cross_terms[],
    int              *cross_parms[],
    Real             *cross_terms[] )
{
    int   p;

    for_less( p, 0, n_parameters )
    {
        if( n_cross_terms[p] > 0 )
        {
            REALLOC( cross_terms[p], n_cross_terms[p] );
            REALLOC( cross_parms[p], n_cross_terms[p] );
        }
    }
}

BICAPI  void  delete_lsq_terms(
    int              n_parameters,
    Real             linear_terms[],
    Real             square_terms[],
    int              n_cross_terms[],
    int              *cross_parms[],
    Real             *cross_terms[] )
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
