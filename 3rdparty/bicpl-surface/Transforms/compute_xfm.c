
/* ----------------------------------------------------------------------------
@COPYRIGHT  :
              Copyright 1993,1994,1995 David MacDonald,
              Peter Neelin, Louis Collins,
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

/* ----------------------------- MNI Header -----------------------------------
@NAME       : compute_xfm.c
@DESCRIPTION: Routine to calculate a General_transform from a pair of tag
              point lists.
@METHOD     : 
@GLOBALS    : 
@CREATED    : August 30, 1993 (Peter Neelin)
@MODIFIED   : $Log: compute_xfm.c,v $
@MODIFIED   : Revision 1.18  2005/08/17 22:26:47  bert
@MODIFIED   : Replace public/private with BICAPI/static
@MODIFIED   :
@MODIFIED   : Revision 1.17  2000/02/06 15:30:50  stever
@MODIFIED   : rearranged header file structure; add gcc -Wall fixes
@MODIFIED   :
@MODIFIED   : Revision 1.16  2000/02/05 21:27:20  stever
@MODIFIED   : change include lines <foo.h> --> <bicpl/foo.h>
@MODIFIED   :
@MODIFIED   : Revision 1.15  1995/07/31 13:46:00  david
@MODIFIED   : check_in_all
@MODIFIED   :
 * Revision 1.14  1995/07/10  18:02:50  david
 * check_in_all
 *
 * Revision 1.13  1995/07/10  14:36:10  david
 * check_in_all
 *
 * Revision 1.12  1995/07/08  03:38:04  david
 * *** empty log message ***
 *
 * Revision 1.11  1995/07/07  18:51:24  david
 * *** empty log message ***
 *
 * Revision 1.10  1995/07/07  18:18:23  david
 * check_in_all
 *
 * Revision 1.9  1995/07/07  18:16:45  david
 * *** empty log message ***
 *
 * Revision 1.8  1995/07/05  14:52:06  david
 * *** empty log message ***
 *
 * Revision 1.7  1995/06/23  14:24:37  david
 * check_in_all
 *
 * Revision 1.6  1995/05/19  21:51:21  david
 * *** empty log message ***
 *
 * Revision 1.5  1995/04/04  03:42:05  david
 * check_in_all
 *
 * Revision 1.4  1995/03/07  18:54:49  david
 * check_in_all
 *
 * Revision 1.3  95/02/27  17:20:20  david
 * *** empty log message ***
 * 
 * Revision 1.2  94/11/25  14:23:13  david
 * check_in_all
 * 
 * Revision 1.1  94/11/04  14:45:47  david
 * Initial revision
 * 
 * Revision 1.4  93/09/10  11:52:58  neelin
 * Fixed 9 and 10 parameter fitting.
 * 
 * Revision 1.3  93/09/08  15:45:00  neelin
 * Added 9 and 10 parameter linear fitting.
 * 
 * Revision 1.2  93/09/01  15:35:17  neelin
 * Changed calls to lfit to lfit_vector.
 * 
 * Revision 1.1  93/09/01  15:25:25  neelin
 * Initial revision
 * 
 * Revision 1.1  93/09/01  13:21:30  neelin
 * Initial revision
 * 
---------------------------------------------------------------------------- */

#include "bicpl_internal.h"

#ifndef lint
static char rcsid[] = "$Header: /private-cvsroot/libraries/bicpl/Transforms/compute_xfm.c,v 1.18 2005/08/17 22:26:47 bert Exp $";
#endif

/* Function declarations */

static void compute_procrustes_transform(int npoints, 
                                          Real **tag_list1, 
                                          Real **tag_list2, 
                                          Trans_type trans_type,
                                          General_transform *transform);
static void compute_arb_param_transform(int npoints, 
                                         Real **tag_list1, 
                                         Real **tag_list2, 
                                         Trans_type trans_type,
                                         General_transform *transform);
static void compute_12param_transform(int npoints, 
                                       Real **tag_list1, 
                                       Real **tag_list2, 
                                       Trans_type trans_type,
                                       General_transform *transform);
static void compute_tps_transform(int npoints, 
                                   Real **tag_list1, 
                                   Real **tag_list2, 
                                   Trans_type trans_type,
                                   General_transform *transform);

static  void  concat_transformation_matrix(
    Transform   *lt, 
    Real        center[],
    Real        translations[],
    Real        scales[],
    Real        shears[],
    Transform   *rotation );

/* ----------------------------- MNI Header -----------------------------------
@NAME       : compute_transform_from_tags
@INPUT      : npoints - number of pairs of tag points
              tag_list1 - first list of tag points
              tag_list2 - second list of tag points
              trans_type - type of transformation to calculate
@OUTPUT     : transform - computed transform
@RETURNS    : (nothing)
@DESCRIPTION: Routine to calculate a general transform from a pair of lists
              of tag points. The transform is from tag_list2 to tag_list1.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : August 30, 1993 (Peter Neelin)
@MODIFIED   : July    4, 1995 D. MacDonald - removed recipes-style code
---------------------------------------------------------------------------- */

BICAPI  void  compute_transform_from_tags(
    int                 npoints, 
    Real                **tag_list1, 
    Real                **tag_list2, 
    Trans_type          trans_type,
    General_transform   *transform)
{
    /* Check number of points for linear transformation */

    if( ((trans_type == TRANS_LSQ6) ||
         (trans_type == TRANS_LSQ7) ||
         (trans_type == TRANS_LSQ9) ||
         (trans_type == TRANS_LSQ10) ||
         (trans_type == TRANS_LSQ12)) &&
        (npoints < MIN_POINTS_LINEAR) )
    {
        create_linear_transform(transform, NULL);
        return;
    }

    /* Check number of points for thin-plate spline transformation */

    if( (trans_type == TRANS_TPS) && (npoints < MIN_POINTS_TPS) )
    {
        create_linear_transform(transform, NULL);
        return;
    }

    switch (trans_type)
    {
    case TRANS_LSQ6:
    case TRANS_LSQ7:
        compute_procrustes_transform( npoints, tag_list1, tag_list2,
                                      trans_type, transform );
        break;

    case TRANS_LSQ9:
    case TRANS_LSQ10:
        compute_arb_param_transform( npoints, tag_list1, tag_list2,
                                     trans_type, transform );
        break;

    case TRANS_LSQ12:
        compute_12param_transform( npoints, tag_list1, tag_list2,
                                   trans_type, transform );
        break;

    case TRANS_TPS:
        compute_tps_transform( npoints, tag_list1, tag_list2,
                               trans_type, transform );
        break;

    default:
        print_error(
            "Invalid transform type in compute_transform_from tags\n" );

        exit(EXIT_FAILURE);
    }
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : compute_procrustes_transform
@INPUT      : npoints - number of pairs of tag points
              tag_list1 - first list of tag points
              tag_list2 - second list of tag points
              trans_type - type of transformation to calculate
@OUTPUT     : transform - computed transform
@RETURNS    : (nothing)
@DESCRIPTION: Routine to calculate a general transform from a pair of lists
              of tag points. The transform is from tag_list2 to tag_list1.
              Uses 6 or 7 parameter procrustes transformation.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : August 30, 1993 (Peter Neelin)
@MODIFIED   : July    4, 1995 D. MacDonald - removed recipes-style code
---------------------------------------------------------------------------- */

static  void  compute_procrustes_transform(
    int                 npoints, 
    Real                **tag_list1, 
    Real                **tag_list2, 
    Trans_type          trans_type,
    General_transform   *transform)
{
    Transform   rotation;
    int         i;
    Real        translation[N_DIMENSIONS];
    Real        centre_of_rotation[N_DIMENSIONS];
    Real        scale, scales[N_DIMENSIONS];
    Real        shears[N_DIMENSIONS];
    Transform   linear_transform;

    /* Do procrustes fit */

    procrustes( npoints, N_DIMENSIONS, tag_list1, tag_list2, translation, 
                centre_of_rotation, &rotation, &scale );

    /* Set scale appropriately */

    if( trans_type == TRANS_LSQ6 )
        scale = 1.0;

    for_less( i, 0, N_DIMENSIONS )
    {
        scales[i] = scale;
        shears[i] = 0.0;
    }

    /* Concatenate translation, scale, shear and rotation */

    concat_transformation_matrix( &linear_transform, centre_of_rotation,
                                  translation, scales, shears, &rotation );

    create_linear_transform( transform, &linear_transform );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : compute_arb_param_transform
@INPUT      : npoints - number of pairs of tag points
              tag_list1 - first list of tag points
              tag_list2 - second list of tag points
              trans_type - type of transformation to calculate
@OUTPUT     : transform - computed transform
@RETURNS    : (nothing)
@DESCRIPTION: Routine to calculate a general transform from a pair of lists
              of tag points. The transform is from tag_list2 to tag_list1.
              Computes a linear transform using arbitrary parameters.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : August 30, 1993 (Louis Collins and Peter Neelin)
@MODIFIED   : July    4, 1995 D. MacDonald - removed recipes-style code
---------------------------------------------------------------------------- */

static  void  compute_arb_param_transform(
    int                 npoints, 
    Real                **tag_list1, 
    Real                **tag_list2, 
    Trans_type          trans_type,
    General_transform   *transform )
{
    Transform  rotation;
    Real       translation[N_DIMENSIONS];
    Real       centre_of_rotation[N_DIMENSIONS], scale;
    Real       scales[N_DIMENSIONS];
    Real       shears[N_DIMENSIONS], angles[N_DIMENSIONS];
    int        i;
    Transform  linear_transform;
  
    if( trans_type != TRANS_LSQ9 && trans_type != TRANS_LSQ10 )
    {
        print_error( "Error in compute_arb_param_transform().\n" );
        print_error( "Trans_type (=%d) requested is not LSQ9 or LSQ10.\n", 
                     trans_type );
        exit(EXIT_FAILURE);
    }

    /* Do procrustes fit */

    procrustes( npoints, N_DIMENSIONS, tag_list1, tag_list2, translation, 
                centre_of_rotation, &rotation, &scale );

    if( !rotmat_to_ang( &rotation, angles ) )
    {
        print_error( "Error in compute_arb_param_transform().\n");
        print_error( "Cannot extract angles from rotation matrix.\n");
        exit(EXIT_FAILURE);
    }
  
    for_less( i, 0, N_DIMENSIONS )
    {
        shears[i] = 0.0;
        scales[i] = scale;
    }

    if( !optimize_simplex( tag_list1, tag_list2,
                           npoints, trans_type,
                           centre_of_rotation, 
                           translation,
                           scales,
                           shears,
                           angles) )
    {
        print_error( "Error in compute_arb_param_transform(),\n");
        print_error( "in call to optimize_simplex().\n");
        exit(EXIT_FAILURE);
    }

    build_transformation_matrix( &linear_transform,
                                 centre_of_rotation, 
                                 translation,
                                 scales,
                                 shears,
                                 angles );


    create_linear_transform( transform, &linear_transform );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : make_rots
@INPUT      : rot_x, rot_y, rot_z - three rotation angles, in radians.
@OUTPUT     : xmat, a linear transform
@RETURNS    : nothing
@DESCRIPTION: to be applied by premultiplication, ie rot*vec = newvec
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : Tue Jun  8 08:44:59 EST 1993 LC
@MODIFIED   : July    4, 1995 D. MacDonald - removed recipes-style code
---------------------------------------------------------------------------- */

static  void   make_rots(
    Transform   *xmat,
    Real        rot_x,
    Real        rot_y,
    Real        rot_z )
{
    Transform  tx, ty, tz;

    make_rotation_transform( -rot_x, X, &tx ) ;
    make_rotation_transform( -rot_y, Y, &ty ) ;
    make_rotation_transform( -rot_z, Z, &tz ) ;

    concat_transforms( xmat, &tx, &ty );
    concat_transforms( xmat, xmat, &tz );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : concat_transformation_matrix
@INPUT      : center, translations, scales, rotations
@OUTPUT     : lt  - a linear transformation matrix
@RETURNS    : nothing
@DESCRIPTION: mat = (c)(sh)(s*r)(-c)(t),
               the matrix is to be  PREmultiplied with a vector (mat*vec)
               when used in the application
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : Thu Jun  3 09:37:56 EST 1993 lc
@MODIFIED   : July    4, 1995 D. MacDonald - removed recipes-style code
---------------------------------------------------------------------------- */

static  void  concat_transformation_matrix(
    Transform   *lt, 
    Real        center[],
    Real        translations[],
    Real        scales[],
    Real        shears[],
    Transform   *rotation )
{
    Transform    T, SH, S, C;

    /* mat = (C)(SH)(S)(R)(-C)(T) */

    /* make (T)(-C) */

    make_translation_transform( translations[X] - center[X],
                                translations[Y] - center[Y],
                                translations[Z] - center[Z], &T );

    /* make shear rotation matrix */

    make_rots( &SH, shears[0], shears[1], shears[2] );

    /* make scaling matrix */

    make_scale_transform( scales[X], scales[Y], scales[Z], &S );

    /* make translation back to centre */

    make_translation_transform( center[X], center[Y], center[Z], &C );

    concat_transforms( lt, &T, rotation );
    concat_transforms( lt, lt, &S );
    concat_transforms( lt, lt, &SH );
    concat_transforms( lt, lt, &C );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : build_transformation_matrix
@INPUT      : center, translations, scales, rotations
@OUTPUT     : lt  - a linear transformation matrix
@RETURNS    : nothing
@DESCRIPTION: mat = (c)(s*r)(-c)(t),
               the matrix is to be  PREmultiplied with a vector (mat*vec)
               when used in the application
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : Thu Jun  3 09:37:56 EST 1993 lc
@MODIFIED   : July    4, 1995 D. MacDonald - removed recipes-style code
---------------------------------------------------------------------------- */

BICAPI  void  build_transformation_matrix(
    Transform   *lt, 
    Real        center[],
    Real        translations[],
    Real        scales[],
    Real        shears[],
    Real        rotations[] )
{
    Transform    R;

    /* make rotation matix */

    make_rots( &R, rotations[0], rotations[1], rotations[2] );

    concat_transformation_matrix( lt, center, translations, scales,
                                  shears, &R );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : compute_12param_transform
@INPUT      : npoints - number of pairs of tag points
              tag_list1 - first list of tag points
              tag_list2 - second list of tag points
              trans_type - type of transformation to calculate
@OUTPUT     : transform - computed transform
@RETURNS    : (nothing)
@DESCRIPTION: Routine to calculate a general transform from a pair of lists
              of tag points. The transform is from tag_list2 to tag_list1.
              Uses 12 parameter linear transformation, and simply solves
              the linear least squares problem for each of the x, y, and
              z components independently.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : August 30, 1993 (Peter Neelin)
@MODIFIED   : 
@MODIFIED   : July    4, 1995 D. MacDonald - removed recipes-style code
---------------------------------------------------------------------------- */

static  void  compute_12param_transform(
    int                 npoints, 
    Real                **tag_list1, 
    Real                **tag_list2, 
    Trans_type          trans_type,
    General_transform   *transform)
{
    Real       *x, solution[N_DIMENSIONS + 1];
    int        d, dim;
    int        point;
    Transform  linear_transform;

    /*--- Check transformation type */

    if( trans_type != TRANS_LSQ12 )
    {
        print_error( "Internal error in compute_12param_transform!\n");
        exit(EXIT_FAILURE);
    }

    /*--- initialize the solution */

    make_identity_transform( &linear_transform );

    ALLOC( x, npoints );

    /*--- for each dimension, find the linear least squares solution */

    for_less( dim, 0, N_DIMENSIONS )
    {
        /*--- Copy the data points into a 1-D array */

        for_less( point, 0, npoints )
            x[point] = tag_list1[point][dim];

        /*--- find the solution */

        least_squares( npoints, N_DIMENSIONS, tag_list2, x, solution );

        /*--- record solution in the linear transform */

        Transform_elem( linear_transform, dim, N_DIMENSIONS ) = solution[0];
        for_less( d, 0, N_DIMENSIONS )
            Transform_elem( linear_transform, dim, d ) = solution[1+d];
    }

    /*--- Create general transform */

    create_linear_transform( transform, &linear_transform );

    /* Free matrices and vectors */

    FREE( x );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : compute_tps_transform
@INPUT      : npoints - number of pairs of tag points
              tag_list1 - first list of tag points
              tag_list2 - second list of tag points
              trans_type - type of transformation to calculate
@OUTPUT     : transform - computed transform
@RETURNS    : (nothing)
@DESCRIPTION: Routine to calculate a general transform from a pair of lists
              of tag points. The transform is from tag_list2 to tag_list1.
              Uses thin-plate spline transform.
@METHOD     : We calculate the tranformation from tag_list1 to tag_list2
              and then invert it so that resampling volumes is faster
              (doing a forward transformation, rather than a iterative
              inverse).
@GLOBALS    : 
@CALLS      : 
@CREATED    : August 30, 1993 (Peter Neelin)
@MODIFIED   : July    4, 1995 D. MacDonald - removed recipes-style code
---------------------------------------------------------------------------- */
 
static  void  compute_tps_transform(
    int                 npoints, 
    Real                **tag_list1, 
    Real                **tag_list2, 
    Trans_type          trans_type,
    General_transform   *transform)
{
    Real               **displacements;
    General_transform  inv_transform;

    /* Check trans_type */

    if (trans_type != TRANS_TPS)
    {
        print_error( "Wrong trans type in compute_tps_transform\n");
        exit(EXIT_FAILURE);
    }

    /* Allocate matrices */

    ALLOC2D( displacements, npoints+1+N_DIMENSIONS, N_DIMENSIONS );

    get_nonlinear_warp( tag_list1, tag_list2, displacements, npoints,
                        N_DIMENSIONS, N_DIMENSIONS );

    /* ---- Create general transform */

    create_thin_plate_transform_real( &inv_transform, N_DIMENSIONS, npoints, 
                                      tag_list1, displacements );

    /* ---- Invert general transform */

    create_inverse_general_transform( &inv_transform, transform );

    /* ---- Free inverse transform */

    delete_general_transform(&inv_transform);

    /* ---- Free displacements */

    FREE2D( displacements );
}
