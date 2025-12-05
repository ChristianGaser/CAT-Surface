/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 *
 */

#ifndef _CAT_MARCHINGCUBES_H_
#define _CAT_MARCHINGCUBES_H_

#include <bicpl/marching.h>
#include "CAT_NiftiLib.h"

#define CHUNK_SIZE    1000000
#define IDX(x, y, z, nx, ny) ((z) * (nx) * (ny) + (y) * (nx) + (x))

void correct_topology(
    float *volume, 
    float thresh, 
    int dims[3], 
    int conn_arr[2], 
    int n_loops);

object_struct *apply_marching_cubes(
    float *input_float, 
    nifti_image *nii_ptr,
    float *label, 
    double min_threshold,
    double pre_fwhm,
    int iter_laplacian,
    double dist_morph,
    int n_median_filter,
    int n_iter,
    double strength_gyri_mask,
    int verbose);

object_struct *apply_marching_cubes_fast(
    float *input_float,
    nifti_image *nii_ptr,
    double min_threshold,
    int verbose);

void extract_isosurface(
    float            *vol,
    int              sizes[3],
    double           min_label,
    double           max_label,
    mat44            nii_mat,
    Marching_cubes_methods method,
    BOOLEAN          binary_flag,
    double           min_threshold,
    double           max_threshold,
    double           valid_low,
    double           valid_high,
    polygons_struct  *polygons,
    int              verbose);

#endif
