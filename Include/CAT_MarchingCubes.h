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

/* argument defaults */
double min_threshold = 0.5;
double post_fwhm = 2.0;
double pre_fwhm = 2.0;
double dist_morph = FLT_MAX;
int n_median_filter = 2;
int verbose = 0;
int n_iter = 10;

void correct_topology(float *volume, float thresh, int dims[3], int conn_arr[2], 
                      int n_loops);
object_struct *apply_marching_cubes(float *input_float, nifti_image *nii_ptr);

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