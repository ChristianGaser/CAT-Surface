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
static char rcsid[] = "$Header: /private-cvsroot/libraries/bicpl/Numerical/least_squares.c,v 1.9 2005/08/17 22:28:59 bert Exp $";
#endif

/* ----------------------------- MNI Header -----------------------------------
@NAME       : least_squares
@INPUT      : n_points
              n_dims
              points
              values
@OUTPUT     : parameters
@RETURNS    : 
@DESCRIPTION: Computes the standard problem of a least squares fit of the
              equation:
              values(points[i]) = parameters[0] +
                                  parameters[1] * points[i][0] +
                                  parameters[2] * points[i][1] +
                                             ...
                                  parameters[n_dims] * points[i][n_dims-1].
              
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  BOOLEAN  least_squares(
    int     n_points,
    int     n_dims,
    Real    **points,
    Real    values[],
    Real    parameters[] )
{
    int                   pt, i;
    linear_least_squares  lsq;
    Real                  *p;
    BOOLEAN               success;

    initialize_linear_least_squares( &lsq, n_dims+1 );

    ALLOC( p, n_dims+1 );

    for_less( pt, 0, n_points )
    {
        p[0] = 1.0;
        for_less( i, 0, n_dims )
            p[i+1] = points[pt][i];

        add_to_linear_least_squares( &lsq, p, values[pt] );
    }

    FREE( p );

    success = get_linear_least_squares_solution( &lsq, parameters );

    delete_linear_least_squares( &lsq );

    return( success );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : initialize_linear_least_squares
@INPUT      : lsq
              n_parameters
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: Initializes the general least squares structure.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  void  initialize_linear_least_squares(
    linear_least_squares   *lsq,
    int                    n_parameters )
{
    int   i, j;

    lsq->n_parameters = n_parameters;

    if( n_parameters > 0 )
    {
        ALLOC2D( lsq->second_derivs, n_parameters, n_parameters );
        ALLOC( lsq->constants, n_parameters );

        for_less( i, 0, n_parameters )
        {
            for_less( j, 0, n_parameters )
                lsq->second_derivs[i][j] = 0.0;
            lsq->constants[i] = 0.0;
        }
    }
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : add_to_linear_least_squares
@INPUT      : lsq
              parameter_coefs[n_parameters]
              constant
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: Adds to the linear least squares.  For instance if you
              are solving for:

               minimize sum ( a * x + b * y + c - z )^2
                a,b,c    n         i       i       i

              Then you would call this function n times, each time with
              parameter_coefs[0] = x[i],
              parameter_coefs[1] = y[i],
              constant = z[i],
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  void  add_to_linear_least_squares(
    linear_least_squares   *lsq,
    Real                   parameter_coefs[],
    Real                   constant )
{
    int   i, j;

    for_less( i, 0, lsq->n_parameters )
    {
        for_less( j, i, lsq->n_parameters )
            lsq->second_derivs[i][j] += parameter_coefs[i] * parameter_coefs[j];

        lsq->constants[i] += parameter_coefs[i] * constant;
    }
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : get_linear_least_squares_solution
@INPUT      : lsq
@OUTPUT     : solution
@RETURNS    : TRUE if successful
@DESCRIPTION: Returns the solution to the linear least squares problem.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  BOOLEAN  get_linear_least_squares_solution(
    linear_least_squares   *lsq,
    Real                   solution[] )
{
    BOOLEAN  solved;
    int      i, j;

    for_less( i, 0, lsq->n_parameters )
    {
        for_less( j, i+1, lsq->n_parameters )
            lsq->second_derivs[j][i] = lsq->second_derivs[i][j];
    }

    solved = solve_linear_system( lsq->n_parameters, lsq->second_derivs,
                                  lsq->constants, solution );

    if( !solved )
    {
        print_error(
             "Could not solve least squares, non-invertible matrix.\n" );

        for_less( i, 0, lsq->n_parameters )
            solution[i] = 0.0;
    }

    return( solved );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : delete_linear_least_squares
@INPUT      : lsq
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: Deletes the linear least squares structure.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  void  delete_linear_least_squares(
    linear_least_squares   *lsq )
{
    FREE2D( lsq->second_derivs );
    FREE( lsq->constants );
}
