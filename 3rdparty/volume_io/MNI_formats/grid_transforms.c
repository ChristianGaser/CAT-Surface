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

#include  <internal_volume_io.h>

#ifndef lint
static char rcsid[] = "$Header: /private-cvsroot/minc/volume_io/MNI_formats/grid_transforms.c,v 1.10.2.1 2004/10/04 20:18:52 bert Exp $";
#endif

#define   DEGREES_CONTINUITY         2    /* Cubic interpolation */
#define   SPLINE_DEGREE         ((DEGREES_CONTINUITY) + 2)

#define   N_COMPONENTS   N_DIMENSIONS /* displacement vector has 3 components */

#define   FOUR_DIMS      4

#define   INVERSE_FUNCTION_TOLERANCE     0.01
#define   INVERSE_DELTA_TOLERANCE        1.0e-5
#define   MAX_INVERSE_ITERATIONS         20

static void   evaluate_grid_volume(
    Volume         volume,
    Real           x,
    Real           y,
    Real           z,
    int            degrees_continuity,
    Real           values[],
    Real           deriv_x[],
    Real           deriv_y[],
    Real           deriv_z[] );

/* ----------------------------- MNI Header -----------------------------------
@NAME       : grid_transform_point
@INPUT      : transform
              x
              y
              z
@OUTPUT     : x_transformed
              y_transformed
              z_transformed
@RETURNS    : 
@DESCRIPTION: Applies a grid transform to the point
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : Feb. 21, 1995    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

VIOAPI  void  grid_transform_point(
    General_transform   *transform,
    Real                x,
    Real                y,
    Real                z,
    Real                *x_transformed,
    Real                *y_transformed,
    Real                *z_transformed )
{
    Real    displacements[N_COMPONENTS];
    Volume  volume;

    /* --- the volume that defines the transform is an offset vector,
           so evaluate the volume at the given position and add the
           resulting offset to the given position */

    volume = (Volume) transform->displacement_volume;

    evaluate_grid_volume( volume, x, y, z, DEGREES_CONTINUITY, displacements,
                          NULL, NULL, NULL );

    *x_transformed = x + displacements[X];
    *y_transformed = y + displacements[Y];
    *z_transformed = z + displacements[Z];
}

#ifdef USE_NEWTONS_METHOD
/* ----------------------------- MNI Header -----------------------------------
@NAME       : forward_function
@INPUT      : function_data  - contains transform info
              parameters     - x,y,z position
@OUTPUT     : values         - where x,y,z, maps to
              derivatives    - the 3 by 3 derivatives of the mapping
@RETURNS    : 
@DESCRIPTION: This function does the same thing as grid_transform_point(),
              but also gets derivatives.  This function is passed to the
              newton function solution routine to perform the inverse mapping
              of the grid transformation.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : Feb.   , 1995    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

private  void  forward_function(
    void   *function_data,
    Real   parameters[],
    Real   values[],
    Real   **derivatives )
{
    int                c;
    General_transform  *transform;
    Real               deriv_x[N_COMPONENTS], deriv_y[N_COMPONENTS];
    Real               deriv_z[N_COMPONENTS];
    Volume             volume;

    transform = (General_transform *) function_data;

    /* --- store the offset vector in values[0-2] */

    volume = transform->displacement_volume;

    evaluate_grid_volume( volume, parameters[X], parameters[Y], parameters[Z],
                          DEGREES_CONTINUITY, values,
                          deriv_x, deriv_y, deriv_z );

    for_less( c, 0, N_COMPONENTS )
    {
        values[c] += parameters[c];   /* to get x',y',z', add offset to x,y,z */

        /*--- given the derivatives of the offset, compute the
              derivatives of (x,y,z) + offset, with respect to x,y,z */

        derivatives[c][X] = deriv_x[c];
        derivatives[c][Y] = deriv_y[c];
        derivatives[c][Z] = deriv_z[c];

        derivatives[c][c] += 1.0;    /* deriv of (x,y,z) w.r.t. x or y or z  */
    }
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : grid_inverse_transform_point
@INPUT      : transform
              x
              y
              z
@OUTPUT     : x_transformed
              y_transformed
              z_transformed
@RETURNS    : 
@DESCRIPTION: Applies the inverse grid transform to the point.  This is done
              by using newton-rhapson steps to find the point which maps to
              the parameters (x,y,z).
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : Feb. 21, 1995    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

/* ---------------------------------------------------------------------------

     There are two different versions of the grid inverse function.  I would
have hoped that my version worked best, since it uses first derivatives and
Newton's method.  However, Louis' version seems to work better, perhaps since
it matches the code he uses in minctracc to generate the grid transforms.

- David MacDonald

---------------------------------------------------------------------------- */

VIOAPI  void  grid_inverse_transform_point(
    General_transform   *transform,
    Real                x,
    Real                y,
    Real                z,
    Real                *x_transformed,
    Real                *y_transformed,
    Real                *z_transformed )
{
    Real   solution[N_DIMENSIONS];
    Real   initial_guess[N_DIMENSIONS];
    Real   desired_values[N_DIMENSIONS];

    /* --- fill in the initial guess */

    initial_guess[X] = x;
    initial_guess[Y] = y;
    initial_guess[Z] = z;

    /* --- define what the desired function values are */

    desired_values[X] = x;
    desired_values[Y] = y;
    desired_values[Z] = z;

    /* --- find the x,y,z that are mapped to the desired values */

    if( newton_root_find( N_DIMENSIONS, forward_function,
                          (void *) transform,
                          initial_guess, desired_values,
                          solution, INVERSE_FUNCTION_TOLERANCE, 
                          INVERSE_DELTA_TOLERANCE, MAX_INVERSE_ITERATIONS ))
    {
        *x_transformed = solution[X];
        *y_transformed = solution[Y];
        *z_transformed = solution[Z];
    }
    else  /* --- if no solution found, not sure what is reasonable to return */
    {
        *x_transformed = x;
        *y_transformed = y;
        *z_transformed = z;
    }
}
#endif

/* ----------------------------- MNI Header -----------------------------------
@NAME       : grid_inverse_transform_point
@INPUT      : transform
              x
              y
              z
@OUTPUT     : x_transformed
              y_transformed
              z_transformed
@RETURNS    : 
@DESCRIPTION: Transforms the point by the inverse of the grid transform.
              Approximates the solution using a simple iterative step
              method.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : 1993?   Louis Collins
@MODIFIED   : 1994    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

VIOAPI  void  grid_inverse_transform_point(
    General_transform   *transform,
    Real                x,
    Real                y,
    Real                z,
    Real                *x_transformed,
    Real                *y_transformed,
    Real                *z_transformed )
{
#define  NUMBER_TRIES  10
    int    tries;
    Real   best_x, best_y, best_z;
    Real   tx, ty, tz;
    Real   gx, gy, gz;
    Real   error_x, error_y, error_z, error, smallest_e;
    Real   ftol;

    ftol = 0.05;

    grid_transform_point( transform, x, y, z, &tx, &ty, &tz );

    tx = x - (tx - x);
    ty = y - (ty - y);
    tz = z - (tz - z);

    grid_transform_point( transform, tx, ty, tz, &gx, &gy, &gz );

    error_x = x - gx;
    error_y = y - gy;
    error_z = z - gz;

    tries = 0;

    error = smallest_e = FABS(error_x) + FABS(error_y) + FABS(error_z);

    best_x = tx;
    best_y = ty;
    best_z = tz;

    while( ++tries < NUMBER_TRIES && smallest_e > ftol )
    {
        tx += 0.67 * error_x;
        ty += 0.67 * error_y;
        tz += 0.67 * error_z;

        grid_transform_point( transform, tx, ty, tz, &gx, &gy, &gz );

        error_x = x - gx;
        error_y = y - gy;
        error_z = z - gz;
    
        error = FABS(error_x) + FABS(error_y) + FABS(error_z);

        if( error < smallest_e )
        {
            smallest_e = error;
            best_x = tx;
            best_y = ty;
            best_z = tz;
        }
    }

    *x_transformed = best_x;
    *y_transformed = best_y;
    *z_transformed = best_z;
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : evaluate_grid_volume
@INPUT      : volume
              voxel
              degrees_continuity
@OUTPUT     : values
              derivs  (if non-NULL)
@RETURNS    : 
@DESCRIPTION: Takes a voxel space position and evaluates the value within
              the volume by nearest_neighbour, linear, quadratic, or
              cubic interpolation.  Rather than use the generic evaluate_volume
              function, this special purpose function is a bit faster.
@CREATED    : Mar. 16, 1995           David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

static  void   evaluate_grid_volume(
    Volume         volume,
    Real           x,
    Real           y,
    Real           z,
    int            degrees_continuity,
    Real           values[],
    Real           deriv_x[],
    Real           deriv_y[],
    Real           deriv_z[] )
{
    Real     voxel[MAX_DIMENSIONS], voxel_vector[MAX_DIMENSIONS];
    int      inc0, inc1, inc2, inc3, inc[MAX_DIMENSIONS], derivs_per_value;
    int      ind0, vector_dim;
    int      start0, start1, start2, start3, inc_so_far;
    int      end0, end1, end2, end3;
    int      v0, v1, v2, v3;
    int      v, d, id, sizes[MAX_DIMENSIONS];
    int      start[MAX_DIMENSIONS];
    int      end[MAX_DIMENSIONS];
    Real     fraction[MAX_DIMENSIONS], bound, pos;
    Real     coefs[SPLINE_DEGREE*SPLINE_DEGREE*SPLINE_DEGREE*N_COMPONENTS];
    Real     values_derivs[N_COMPONENTS + N_COMPONENTS * N_DIMENSIONS];

    convert_world_to_voxel( volume, x, y, z, voxel );

    if( get_volume_n_dimensions(volume) != FOUR_DIMS )
        handle_internal_error( "evaluate_grid_volume" );

    /*--- find which of 4 dimensions is the vector dimension */

    for_less( vector_dim, 0, FOUR_DIMS )
    {
        for_less( d, 0, N_DIMENSIONS )
        {
            if( volume->spatial_axes[d] == vector_dim )
                break;
        }

        if( d == N_DIMENSIONS )
            break;
    }

    get_volume_sizes( volume, sizes );

    bound = (Real) degrees_continuity / 2.0;

    /*--- if near the edges, reduce the degrees of continuity */

    for_less( d, 0, FOUR_DIMS )
    {
        if( d != vector_dim )
        {
            while( degrees_continuity >= -1 &&
                   (voxel[d] < bound  ||
                    voxel[d] > (Real) sizes[d] - 1.0 - bound ||
                    bound == (Real) sizes[d] - 1.0 - bound ) )
            {
                --degrees_continuity;
                if( degrees_continuity == 1 )
                    degrees_continuity = 0;
                bound = (Real) degrees_continuity / 2.0;
            }
        }
    }

    /*--- check to fill in the first derivative */

    if( degrees_continuity < 0 && deriv_x != NULL )
    {
        for_less( v, 0, N_COMPONENTS )
        {
            deriv_x[v] = 0.0;
            deriv_y[v] = 0.0;
            deriv_z[v] = 0.0;
        }
    }

    /*--- if outside */

    if( degrees_continuity == -2 )
    {
        for_less( v, 0, N_COMPONENTS )
            values[v] = 0.0;

        return;
    }

    /*--- determine the starting positions in the volume to grab control
          vertices */

    id = 0;
    for_less( d, 0, FOUR_DIMS )
    {
        if( d != vector_dim )
        {
            pos = voxel[d] - bound;
            start[d] =       FLOOR( pos );
            if( voxel[d] == (Real) sizes[d] - 1.0 - bound )
            {
                --start[d];
                fraction[id] = 1.0;
            }
            else
                fraction[id] = pos - (Real) start[d];

            end[d] = start[d] + degrees_continuity + 2;
            ++id;
        }
    }

    /*--- create the strides */

    start[vector_dim] = 0;
    end[vector_dim] = N_COMPONENTS;

    inc_so_far = N_COMPONENTS;
    for_down( d, FOUR_DIMS-1, 0 )
    {
        if( d != vector_dim )
        {
            inc[d] = inc_so_far;
            inc_so_far *= degrees_continuity + 2;
        }
    }

    /*--- copy stride arrays to variables for speed */

    inc[vector_dim] = 1;

    start0 = start[0];
    start1 = start[1];
    start2 = start[2];
    start3 = start[3];

    end0 = end[0];
    end1 = end[1];
    end2 = end[2];
    end3 = end[3];

    inc0 = inc[0] - inc[1] * (end1 - start1);
    inc1 = inc[1] - inc[2] * (end2 - start2);
    inc2 = inc[2] - inc[3] * (end3 - start3);
    inc3 = inc[3];

    /*--- extract values from volume */

    ind0 = 0;

    for_less( v0, start0, end0 )
    {
        for_less( v1, start1, end1 )
        {
            for_less( v2, start2, end2 )
            {
                for_less( v3, start3, end3 )
                {
                    GET_VALUE_4D_TYPED( coefs[ind0], (Real), volume,
                                        v0, v1, v2, v3 );
                    ind0 += inc3;
                }
                ind0 += inc2;
            }
            ind0 += inc1;
        }
        ind0 += inc0;
    }

    /*--- interpolate values */

    if( degrees_continuity == -1 )
    {
        for_less( v, 0, N_COMPONENTS )
            values[v] = coefs[v];
    }
    else
    {
        evaluate_interpolating_spline( N_DIMENSIONS, fraction,
                                       degrees_continuity + 2, 
                                       N_COMPONENTS, coefs, 0, values_derivs );

        /*--- extract values and derivatives from values_derivs */

        if( deriv_x != NULL )
            derivs_per_value = 8;
        else
            derivs_per_value = 1;

        for_less( v, 0, N_COMPONENTS )
            values[v] = values_derivs[v*derivs_per_value];

        if( deriv_x != NULL )
        {
            for_less( v, 0, N_COMPONENTS )
            {
                id = 0;
                for_less( d, 0, FOUR_DIMS )
                {
                    if( d != vector_dim )
                    {
                        voxel_vector[d] = values_derivs[v*8 + (4>>id)];
                        ++id;
                    }
                    else
                        voxel_vector[d] = 0.0;
                }
             
                convert_voxel_normal_vector_to_world( volume, voxel_vector,
                                    &deriv_x[v], &deriv_y[v], &deriv_z[v] );
            }
        }
    }
}
