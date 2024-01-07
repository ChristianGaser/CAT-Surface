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
static char rcsid[] = "$Header: /private-cvsroot/libraries/bicpl/Numerical/gradient_minimize.c,v 1.6 2005/08/17 22:28:59 bert Exp $";
#endif

#define  SMALLEST_STEP_SIZE     1.0e-20

#define  STEP_RATIO   0.5

static  Real  take_step(
    int   n_dims,
    Real  current[],
    Real  current_value,
    Real  (*func) ( Real [], void * ),
    void  *func_data,
    Real  step_sizes[],
    int   n_search_dims,
    Real  parameters[],
    Real  best[],
    int   min_pos[],
    int   max_pos[],
    int   pos[] );

BICAPI  Real  gradient_steps_minimize_function(
    int   n_dims,
    Real  initial_parameters[],
    Real  initial_step_sizes[],
    Real  (*func) ( Real [], void * ),
    void  *func_data,
    int   n_search_dims,
    int   max_iterations,
    Real  tolerance,
    Real  solution[] )
{
    BOOLEAN  done;
    int      iteration, dim, *min_pos, *max_pos, *pos;
    Real     *step_sizes, *parameters, *best, best_value, next_value;

    ALLOC( step_sizes, n_dims );

    for_less( dim, 0, n_dims )
    {
        step_sizes[dim] = initial_step_sizes[dim];
        solution[dim] = initial_parameters[dim];
    }

    best_value = func( solution, func_data );

    ALLOC( parameters, n_dims );
    ALLOC( best, n_dims );
    ALLOC( min_pos, n_dims );
    ALLOC( max_pos, n_dims );
    ALLOC( pos, n_dims );

    done = FALSE;
    iteration = 0;
    while( !done && (max_iterations <= 0 || iteration < max_iterations) )
    {
        next_value = take_step( n_dims, solution, best_value, func, func_data,
                                step_sizes, n_search_dims,
                                parameters, best, min_pos, max_pos, pos );

        if( next_value < best_value &&
            !numerically_close( next_value, best_value, tolerance ) )
        {
            best_value = next_value;
            for_less( dim, 0, n_dims )
            {
                step_sizes[dim] *= 1.0 / STEP_RATIO;
                if( step_sizes[dim] > initial_step_sizes[dim] )
                    step_sizes[dim] = initial_step_sizes[dim];
            }
            done = FALSE;
        }
        else
        {
            done = TRUE;
            for_less( dim, 0, n_dims )
            {
                if( step_sizes[dim] > SMALLEST_STEP_SIZE )
                    done = FALSE;

                step_sizes[dim] *= STEP_RATIO;
                if( step_sizes[dim] < SMALLEST_STEP_SIZE )
                {
                    step_sizes[dim] = SMALLEST_STEP_SIZE;
                }
            }
        }

        ++iteration;
    }

    FREE( parameters );
    FREE( best );
    FREE( min_pos );
    FREE( max_pos );
    FREE( pos );

    return( best_value );
}

static  Real  take_step(
    int   n_dims,
    Real  current[],
    Real  current_value,
    Real  (*func) ( Real [], void * ),
    void  *func_data,
    Real  step_sizes[],
    int   n_search_dims,
    Real  parameters[],
    Real  best[],
    int   min_pos[],
    int   max_pos[],
    int   pos[] )
{
    int   dim, fit_dim, changed_from, n_fit_dims;
    Real  best_value, value;

    for_less( dim, 0, n_dims )
        best[dim] = current[dim];
    best_value = current_value;

    if( n_search_dims < 1 )
        n_search_dims = 1;

    for( fit_dim = 0;  fit_dim < n_dims;  fit_dim += n_search_dims )
    {
        n_fit_dims = MIN( n_search_dims, n_dims - fit_dim );

        for_less( dim, 0, n_dims )
        {
            min_pos[dim] = 0;
            max_pos[dim] = 0;
            parameters[dim] = best[dim];
        }

        for_less( dim, fit_dim, MAX( n_dims, fit_dim + n_fit_dims ) )
        {
            min_pos[dim] = -1;
            max_pos[dim] = 1;
        }

        for_less( dim, 0, n_dims )
            pos[dim] = min_pos[dim];

        changed_from = 0;
        while( changed_from >= 0 )
        {
            for_less( dim, changed_from, n_fit_dims )
            {
                parameters[fit_dim+dim] = current[fit_dim+dim] +
                            step_sizes[fit_dim+dim] * (Real) pos[fit_dim+dim];
            }

            value = func( parameters, func_data );

            if( value < best_value )
            {
                best_value = value;
                for_less( dim, 0, n_dims )
                    best[dim] = parameters[dim];
            }

            changed_from = n_fit_dims - 1;

            do
            {
                ++pos[fit_dim+changed_from];
                if( pos[fit_dim+changed_from] <= max_pos[fit_dim+changed_from] )
                    break;

                pos[fit_dim+changed_from] = min_pos[fit_dim+changed_from];
                --changed_from;
            }
            while( changed_from >= 0 );
        }
    }

    for_less( dim, 0, n_dims )
        current[dim] = best[dim];

    return( best_value );
}
