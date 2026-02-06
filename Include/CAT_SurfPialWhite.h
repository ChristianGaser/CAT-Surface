/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 *
 */

#ifndef _CAT_SURF_PIAL_WHITE_H_
#define _CAT_SURF_PIAL_WHITE_H_

/**
 * @file CAT_SurfPialWhite.h
 * @brief Library for estimating pial and white matter surfaces from central surface.
 *
 * Provides functions for estimating pial and white matter surfaces from a
 * central surface using cortical thickness values and volume label data.
 */

#include <bicpl.h>
#include "CAT_NiftiLib.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Options for pial/white surface estimation.
 */
typedef struct {
    double w1;          /**< Internal smoothness weight (default: 0.05) */
    double w2;          /**< Gradient alignment weight (default: 0.05) */
    double w3;          /**< Balloon force weight (default: 0.05) */
    double sigma;       /**< Displacement smoothing sigma (default: 0.2) */
    int iterations;     /**< Number of deformation iterations (default: 100) */
    int verbose;        /**< Verbose output (default: 0) */
} CAT_PialWhiteOptions;

/**
 * Initialize pial/white estimation options with default values.
 *
 * @param opts Pointer to options structure.
 */
void CAT_PialWhiteOptionsInit(CAT_PialWhiteOptions *opts);

/**
 * Estimate pial and white matter surfaces from central surface.
 *
 * This function:
 * 1. Creates initial pial/white estimates using thickness values
 * 2. Smooths pial surface with curvature-guided blending
 * 3. Performs dual-surface deformation using intensity gradients
 *
 * @param central           Input central surface.
 * @param thickness_values  Per-vertex thickness values.
 * @param labels            NIfTI label volume (tissue classes).
 * @param nii_ptr           NIfTI image header for coordinate transforms.
 * @param pial_out          Output: pial surface (caller must allocate).
 * @param white_out         Output: white surface (caller must allocate).
 * @param opts              Options controlling the algorithm.
 *
 * @return 0 on success, non-zero on error.
 */
int CAT_SurfEstimatePialWhite(
    polygons_struct *central,
    const double *thickness_values,
    float *labels,
    nifti_image *nii_ptr,
    polygons_struct *pial_out,
    polygons_struct *white_out,
    const CAT_PialWhiteOptions *opts
);

#ifdef __cplusplus
}
#endif

#endif /* _CAT_SURF_PIAL_WHITE_H_ */
