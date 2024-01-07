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
static char rcsid[] = "$Header: /private-cvsroot/libraries/bicpl/Numerical/amoeba.c,v 1.15 2005/08/17 22:28:58 bert Exp $";
#endif

#define  FLIP_RATIO      1.0
#define  CONTRACT_RATIO  0.5
#define  STRETCH_RATIO   2.0

/* ----------------------------- MNI Header -----------------------------------
@NAME       : get_function_value
@INPUT      : amoeba
              parameters
@OUTPUT     : 
@RETURNS    : function value
@DESCRIPTION: Evaluates the function being minimized by amoeba, by calling
              the user function.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

static  Real  get_function_value(
    amoeba_struct  *amoeba,
    float          parameters[] )
{
    return( (*amoeba->function) ( amoeba->function_data, parameters ) );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : initialize_amoeba
@INPUT      : n_parameters
              initial_parameters
              parameter_deltas
              function
              function_data
              tolerance
@OUTPUT     : amoeba
@RETURNS    : 
@DESCRIPTION: Initializes the amoeba structure to minimize the function.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  void  initialize_amoeba(
    amoeba_struct     *amoeba,
    int               n_parameters,
    Real              initial_parameters[],
    Real              parameter_deltas[],
    amoeba_function   function,
    void              *function_data,
    Real              tolerance )
{
    int    i, j;

    amoeba->n_parameters = n_parameters;
    amoeba->function = function;
    amoeba->function_data = function_data;
    amoeba->tolerance = tolerance;
    amoeba->n_steps_no_improvement = 0;
    ALLOC2D( amoeba->parameters, n_parameters+1, n_parameters );
    ALLOC( amoeba->values, n_parameters+1 );

    ALLOC( amoeba->sum, n_parameters );

    for_less( j, 0, n_parameters )
        amoeba->sum[j] = 0.0;

    for_less( i, 0, n_parameters+1 )
    {
        for_less( j, 0, n_parameters )
        {
            amoeba->parameters[i][j] = (float) initial_parameters[j];
            if( i > 0 && j == i - 1 )
                amoeba->parameters[i][j] = (float) parameter_deltas[j];
            amoeba->sum[j] += (Real) amoeba->parameters[i][j];
        }

        amoeba->values[i] = get_function_value( amoeba, amoeba->parameters[i] );
    }
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : get_amoeba_parameters
@INPUT      : amoeba
@OUTPUT     : parameters
@RETURNS    : function value
@DESCRIPTION: Passes back the current position of the amoeba (best value),
              and returns the function value at that point.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  Real  get_amoeba_parameters(
    amoeba_struct  *amoeba,
    Real           parameters[] )
{
    int   i, j, low;

    low = 0;
    for_less( i, 1, amoeba->n_parameters+1 )
    {
        if( amoeba->values[i] < amoeba->values[low] )
            low = i;
    }

    for_less( j, 0, amoeba->n_parameters )
        parameters[j] = (Real) amoeba->parameters[low][j];

    return( amoeba->values[low] );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : terminate_amoeba
@INPUT      : amoeba
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: Frees the amoeba.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  void  terminate_amoeba(
    amoeba_struct  *amoeba )
{
    FREE2D( amoeba->parameters );
    FREE( amoeba->values );
    FREE( amoeba->sum );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : try_amoeba
@INPUT      : amoeba
              sum
              high
              fac
@OUTPUT     : 
@RETURNS    : value
@DESCRIPTION: Does a modification to the high vertex of the amoeba and
              returns the value of the new point.  If the new point is
              better (smaller value), it replaces the high vertex of the
              amoeba.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

static  Real  try_amoeba(
    amoeba_struct  *amoeba,
    Real           sum[],
    int            high,
    Real           fac )
{
    int    j;
    Real   y_try, fac1, fac2;
    float  *parameters;

    ALLOC( parameters, amoeba->n_parameters );

    fac1 = (1.0 - fac) / (Real) amoeba->n_parameters;
    fac2 = fac - fac1;

    for_less( j, 0, amoeba->n_parameters )
        parameters[j] = (float) (sum[j] * fac1 +
                                 (Real) amoeba->parameters[high][j] * fac2);

    y_try = get_function_value( amoeba, parameters );

    if( y_try < amoeba->values[high] )
    {
        amoeba->values[high] = y_try;
        for_less( j, 0, amoeba->n_parameters )
        {
            sum[j] += (Real) parameters[j] - (Real) amoeba->parameters[high][j];
            amoeba->parameters[high][j] = parameters[j];
        }
    }

    FREE( parameters );

    return( y_try );
}

#define  N_STEPS_NO_IMPROVEMENT  20

/* ----------------------------- MNI Header -----------------------------------
@NAME       : perform_amoeba
@INPUT      : amoeba
@OUTPUT     : 
@RETURNS    : TRUE if numerically significant improvement
@DESCRIPTION: Performs one iteration of an amoeba, returning true if a
              numerically significant improvement has been found recently.
              Even if it returns FALSE, you can keep calling this function,
              since it may be contracting with no improvement, but will 
              eventually shrink small enough to get an improvment.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  BOOLEAN  perform_amoeba(
    amoeba_struct  *amoeba )
{
    int      i, j, low, high, next_high;
    Real     y_try, y_save;
    BOOLEAN  improvement_found;

    improvement_found = TRUE;

    if( amoeba->values[0] > amoeba->values[1] )
    {
        high = 0;
        next_high = 1;
    }
    else
    {
        high = 1;
        next_high = 0;
    }

    low = next_high;

    for_less( i, 2, amoeba->n_parameters+1 )
    {
        if( amoeba->values[i] < amoeba->values[low] )
            low = i;
        else if( amoeba->values[i] > amoeba->values[high] )
        {
            next_high = high;
            high = i;
        }
        else if( amoeba->values[i] > amoeba->values[next_high] )
            next_high = i;
    }

    if( numerically_close( amoeba->values[low], amoeba->values[high],
                           amoeba->tolerance ) )
    {
        ++amoeba->n_steps_no_improvement;
        if( amoeba->n_steps_no_improvement == N_STEPS_NO_IMPROVEMENT )
            improvement_found = FALSE;
    }
    else
        amoeba->n_steps_no_improvement = 0;

    y_try = try_amoeba( amoeba, amoeba->sum, high, -FLIP_RATIO );

    if( y_try <= amoeba->values[low] )
        y_try = try_amoeba( amoeba, amoeba->sum, high, STRETCH_RATIO );
    else if( y_try >= amoeba->values[next_high] )
    {
        y_save = amoeba->values[high];
        y_try = try_amoeba( amoeba, amoeba->sum, high, CONTRACT_RATIO );

        if( y_try >= y_save )
        {
            for_less( i, 0, amoeba->n_parameters+1 )
            {
                if( i != low )
                {
                    for_less( j, 0, amoeba->n_parameters )
                    {
                        amoeba->parameters[i][j] = (amoeba->parameters[i][j] +
                                            amoeba->parameters[low][j]) / 2.0f;
                    }

                    amoeba->values[i] = get_function_value( amoeba,
                                                  amoeba->parameters[i] );
                }
            }

            for_less( j, 0, amoeba->n_parameters )
            {
                amoeba->sum[j] = 0.0;
                for_less( i, 0, amoeba->n_parameters+1 )
                    amoeba->sum[j] += (Real) amoeba->parameters[i][j];
            }
        }
    }

    return( improvement_found );
}
