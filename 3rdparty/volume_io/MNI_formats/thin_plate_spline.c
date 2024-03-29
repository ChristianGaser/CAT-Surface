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

#include <internal_volume_io.h>

#ifndef lint
static char rcsid[] = "$Header: /private-cvsroot/minc/volume_io/MNI_formats/thin_plate_spline.c,v 1.13.2.1 2004/10/04 20:18:52 bert Exp $";
#endif

#define   INVERSE_FUNCTION_TOLERANCE     0.01
#define   INVERSE_DELTA_TOLERANCE        0.01
#define   MAX_INVERSE_ITERATIONS         20

/* ----------------------------- MNI Header -----------------------------------
@NAME       : thin_plate_spline.c
@INPUT      : 
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: library of routines for warping/mapping transformations.
@METHOD     : The original source code for these routines are taken from
              program VOI, written by Weiqian Dai.
@GLOBALS    : 
@CALLS      : 
@CREATED    : Dec 2, 1991 LC
@MODIFIED   : Mon Apr  5 09:00:54 EST 1993 louis 
                - building new library routines, with prototypes
@MODIFIED   : Wed Jul  14 1993  david 
                - incorporated into libmni.c
              Feb. 28, 1995     D. MacDonald
                - rewrote to get rid of mnewt and floats
---------------------------------------------------------------------------- */


/* ----- structure used by newton root finding ---- */

typedef  struct
{
    Real   **points;
    Real   **weights;
    int    n_points;
    int    n_dims;
} spline_data_struct;

/*------------ private functions -----------------*/

static  void   newton_function(
    void     *function_data,
    Real     parameters[],
    Real     values[],
    Real     **first_derivs );

static  Real  thin_plate_spline_U_deriv(
   Real   pos[],
   Real   landmark[],
   int    n_dims,
   int    deriv_dim );

/* ----------------------------- MNI Header -----------------------------------
@NAME       : evaluate_thin_plate_spline
@INPUT      : n_dims           - dimensionality of the function
              n_values         - number of values of the function
              n_points         - number of defining landmarks
              points[n_points][n_dims]  - landmarks
              weights[n_points+1+n_dims][n_values] - weights for the points
              pos[n_dims]      - position at which to evaluate
@OUTPUT     : values[n_values] - function values at this position
              deriv[n_values][n_dims] - function derivatives at this point
@RETURNS    : 
@DESCRIPTION: Evaluates the thin plate spline at the given point, and, if
              the argument is non-null, the derivatives also.  The thin-plate
              spline takes a point in n_dims dimensional space and returns
              a point in n_values dimensional space.  When used for transforms,
              as in this file, n_values == n_dims, but the code will work for
              the general case where n_values != n_dims.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :
@MODIFIED   : Feb. 27, 1995    David MacDonald 
              modified from some code written by WeiQian Dai later modified by
              Louis Collins
---------------------------------------------------------------------------- */

VIOAPI  void  evaluate_thin_plate_spline(
    int     n_dims,
    int     n_values,
    int     n_points,
    Real    **points,
    Real    **weights,
    Real    pos[],
    Real    values[],
    Real    **derivs )
{
    int       v, d, p;
    Real      dist, dist_deriv;

    /* f(x,y[,z]) =a_{n} + a_{n+1}x + a_{n+1}y + sum_{0}^{n-1}
     *          w_{i}U(|P_{i} - (x,y)|) 
     */

    /* --- initialize derivatives, if desired */

    if( derivs != NULL )
    {
        for_less( v, 0, n_values )
           for_less( d, 0, n_dims )
               derivs[v][d] = 0.0;
    }

    /* --- initialize value of thin plate spline to 0 */

    for_less( v, 0, n_values )
        values[v] = 0.0;

    /* --- for each point, add its contribution to the values and derivs */

    for_less( p, 0, n_points )
    {
        /* --- the thin plate spline weighting function for this point */

        dist = thin_plate_spline_U( pos, points[p], n_dims );

        /* --- add the weighted component to the values */

        for_less( v, 0, n_values )
            values[v] = values[v] + (Real) weights[p][v] * dist;

        /* --- add the weighted component to the derivatives */

        if( derivs != NULL )
        {
            for_less( v, 0, n_values )
            {
                for_less( d, 0, n_dims )
                {
                    dist_deriv = thin_plate_spline_U_deriv( pos, points[p],
                                                            n_dims, d );
                    derivs[v][d] += (Real) weights[p][v] * dist_deriv;
                }
            }
        }
    }

    /* --- add the constant component to the values */

    for_less( v, 0, n_values )
        values[v] += (Real) weights[n_points][v];

    /* --- add the linear components to the values and derivatives */

    for_less( v, 0, n_values )
    {
        for_less( d, 0, n_dims )
        {
            values[v]    += (Real) weights[n_points+1+d][v] * pos[d];
            if( derivs != NULL )
                derivs[v][d] += (Real) weights[n_points+1+d][v];
        }
    }
}


/* ----------------------------- MNI Header -----------------------------------
@NAME       : thin_plate_spline_transform
@INPUT      : 
	      n_dims    - number of dimensions (either 2 or 3).
              n_points  - number of points
              points - array with 'n_points' rows, and 'n_dims' cols,
                       the list of landmarks points in the 'source' volume
	      weights - array with 'n_points+n_dims+1' rows, and 'n_dims' cols,
	               the deformation weights that define the thin plate spline
	      n_points - number of landmark points
              x        - coordinate to transform
              y
              z
@OUTPUT     : x_transformed
              y_transformed
              z_transformed
@RETURNS    : nothing
@DESCRIPTION: Transforms the 1,2, or 3D point given the thin plate spline
              transform
@METHOD     : 
@GLOBALS    : none
@CALLS      : 
@CREATED    : Mon Apr  5 09:00:54 EST 1993
@MODIFIED   : Feb. 27, 1995   D. MacDonald -
                    reorganized to call evaluate_thin_plane_spline()
---------------------------------------------------------------------------- */

VIOAPI  void  thin_plate_spline_transform(
    int     n_dims,
    int     n_points,
    Real    **points,
    Real    **weights,
    Real    x,
    Real    y,
    Real    z,
    Real    *x_transformed,
    Real    *y_transformed,
    Real    *z_transformed )
{
    Real      input_point[N_DIMENSIONS], output_point[N_DIMENSIONS];

    input_point[0] = x;
    input_point[1] = y;
    input_point[2] = z;

    evaluate_thin_plate_spline( n_dims, n_dims, n_points,
                                points, weights, input_point, output_point,
                                NULL );

    *x_transformed = output_point[0];

    if( n_dims >= 2 )
        *y_transformed = output_point[1];

    if( n_dims >= 3 )
        *z_transformed = output_point[2];
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : thin_plate_spline_inverse_transform
@INPUT      : 
	      n_dims    - number of dimensions (either 2 or 3).
              n_points  - number of points
              points - array with 'n_points' rows, and 'n_dims' cols,
                       the list of landmarks points in the 'source' volume
	      weights - array with 'n_points+n_dims+1' rows, and 'n_dims' cols,
	               the deformation weights that define the thin plate spline
	      n_points - number of landmark points
              x        - coordinate to inverse transform
              y
              z
@OUTPUT     : x_transformed
              y_transformed
              z_transformed
@RETURNS    : nothing
@DESCRIPTION: Inverse transforms the 1,2, or 3D point given the thin plate
              spline transform.
@METHOD     : 
@GLOBALS    : none
@CALLS      : 
@CREATED    : Mon Apr  5 09:00:54 EST 1993
@MODIFIED   : Feb. 27, 1995   D. MacDonald -
                    reorganized to call evaluate_thin_plane_spline()
---------------------------------------------------------------------------- */

VIOAPI  void  thin_plate_spline_inverse_transform(
    int     n_dims,
    int     n_points,
    Real    **points,
    Real    **weights,
    Real    x,
    Real    y,
    Real    z,
    Real    *x_transformed,
    Real    *y_transformed,
    Real    *z_transformed )
{
    Real                x_in[N_DIMENSIONS], solution[N_DIMENSIONS];
    spline_data_struct  data;
  
    x_in[X] = x;

    if( n_dims >= 2 )
        x_in[Y] = y;
    else
        x_in[Y] = 0.0;

    if( n_dims >= 3 )
        x_in[Z] = z;
    else
        x_in[Z] = 0.0;

    data.points = points;
    data.weights = weights;
    data.n_points = n_points;
    data.n_dims = n_dims;

    /* --- solve for the root of the function using Newton steps,
           which require a function (newton_function) that evaluates the
           thin plate spline and its derivative at an arbitrary point */

    if( newton_root_find( n_dims, newton_function, (void *) &data,
                          x_in, x_in, solution, INVERSE_FUNCTION_TOLERANCE,
                          INVERSE_DELTA_TOLERANCE, MAX_INVERSE_ITERATIONS ) )
    {
        *x_transformed = solution[0];
        *y_transformed = solution[1];
        *z_transformed = solution[2];
    }
    else
    {
        *x_transformed = x_in[0];
        *y_transformed = x_in[1];
        *z_transformed = x_in[2];
    }
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : newton_function
@INPUT      : function_data
              parameters
@OUTPUT     : values
              first_derivs
@RETURNS    : 
@DESCRIPTION: This function is passed to the newton function root finding
              routine, and evaluates the values and derivatives of the
              thin plate spline.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : Feb. 27, 1995    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

static  void   newton_function(
    void     *function_data,
    Real     parameters[],
    Real     values[],
    Real     **first_derivs )
{
    spline_data_struct *spline_data;

    spline_data = (spline_data_struct *) function_data;

    evaluate_thin_plate_spline( spline_data->n_dims, spline_data->n_dims,
                                spline_data->n_points,
                                spline_data->points, spline_data->weights,
                                parameters, values,
                                first_derivs );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : thin_plate_spline_U
@INPUT      : pos       - position at which to evaluate
              landmark  - landmark
              n_dims      number of dimensions (1,2, or 3)
@OUTPUT     : 
@RETURNS    : U interpolation function of distance between the two args 
@DESCRIPTION: Returns the U interpolation function of the distance between
              points.  In order to correspond to a thin-plate spline, this
              function has a different form in each of the 3 dimensions.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : Feb.   , 1995    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

VIOAPI  Real  thin_plate_spline_U(
    Real   pos[],
    Real   landmark[],
    int    n_dims )
{
    Real r, fu, dx, dy, dz;

    switch( n_dims )
    {
    case  1:
        dx = pos[X] - landmark[X];
        r = FABS( dx );
        fu = r * r * r;
        break;

    case  2:   /* r is actually r^2 */
        dx = pos[X] - landmark[X];
        dy = pos[Y] - landmark[Y];
        r = dx * dx + dy * dy;

        if( r == 0.0 )
            fu = 0.0;
        else
            fu = r * log( r );
        break;

    case  3:
        dx = pos[X] - landmark[X];
        dy = pos[Y] - landmark[Y];
        dz = pos[Z] - landmark[Z];
        r = sqrt( dx * dx + dy * dy + dz * dz );
        fu = r;
        break;

    default:
        handle_internal_error( " impossible error in FU" );
        fu = 0.0;
        break;
    }

    return( fu );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : thin_plate_spline_U_deriv
@INPUT      : pos       - position at which to evaluate
              landmark  - landmark
              n_dims    - number of dimensions (1,2, or 3)
              deriv_dim - dimension to differentiate
@OUTPUT     : 
@RETURNS    : derivative of U interpolation function
@DESCRIPTION: Returns the derivative of the U interpolation function of the
              distance between points (as specified by thin_plate_spline_U()
              above).
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : Feb.   , 1995    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

static  Real  thin_plate_spline_U_deriv(
   Real   pos[],
   Real   landmark[],
   int    n_dims,
   int    deriv_dim )
{
    Real r, r2, deriv, delta[N_DIMENSIONS];

    switch( n_dims )
    {
    case  1:
        delta[X] = pos[X] - landmark[X];
        r = delta[X];
        deriv = 3.0 * r * r;
        break;

    case  2:   /* r2 is r^2 */
        delta[X] = pos[X] - landmark[X];
        delta[Y] = pos[Y] - landmark[Y];
        r2 = delta[X] * delta[X] + delta[Y] * delta[Y];

        if( r2 == 0.0 )
            deriv = 0.0;
        else
            deriv = (1.0 + log( r2 )) * 2.0 * delta[deriv_dim];
        break;

    case  3:
        delta[X] = pos[X] - landmark[X];
        delta[Y] = pos[Y] - landmark[Y];
        delta[Z] = pos[Z] - landmark[Z];
        r = sqrt( delta[X] * delta[X] + delta[Y] * delta[Y] +
                  delta[Z] * delta[Z] );

        if( r == 0.0 )
            deriv = 0.0;
        else
            deriv = delta[deriv_dim] / r;
        break;

    default:
        handle_internal_error( " invalid dimensions error in FU" );
        deriv = 0.0;
        break;
    }

    return( deriv );
}
