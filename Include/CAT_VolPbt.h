/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 *
 */

#ifndef _CAT_VOL_PBT_H_
#define _CAT_VOL_PBT_H_

/**
 * @file CAT_VolPbt.h
 * @brief Projection-based cortical thickness (PBT) estimation library.
 *
 * This module provides the core algorithm for projection-based cortical
 * thickness estimation and percentage position mapping (PPM) as described in:
 *
 *   Dahnke R, Yotter RA, Gaser C.
 *   Cortical thickness and central surface estimation.
 *   Neuroimage. 2013 Jan 15;65:336-48.
 */

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Options structure for controlling PBT estimation behavior.
 */
typedef struct {
    int n_avgs;              /**< Number of averages for distance estimation (default: 2) */
    int n_median_filter;     /**< Iterations of median filter for topology cleanup (default: 2) */
    double range;            /**< Extended range for masking euclidean distance (default: 0.45) */
    double fill_thresh;      /**< Threshold for filling holes in PPM (default: 0.5, 0=disable) */
    double correct_voxelsize;/**< Amount of thickness correction for voxel-size (default: 0.5) */
    int fast;                /**< Fast mode: rougher but quicker estimate (default: 0) */
    int verbose;             /**< Verbose output (default: 0) */
} CAT_PbtOptions;

/**
 * Initialize PBT options with default values.
 *
 * @param opts Pointer to options structure to initialize.
 */
void CAT_PbtOptionsInit(CAT_PbtOptions *opts);

/**
 * Compute projection-based cortical thickness and percentage position map.
 *
 * This function implements the complete PBT pipeline:
 * 1. Distance estimation in WM and CSF by shifting GM borders
 * 2. Thickness estimation via projection_based_thickness()
 * 3. PPM calculation with gyrus/sulcus blending
 * 4. Optional median filtering for topology artifact reduction
 *
 * @param src         Input PVE label image (CSF=1, GM=2, WM=3 with partial volumes).
 * @param GMT_out     Output: Gray matter thickness map (caller must allocate nvox floats).
 * @param PPM_out     Output: Percentage position map (0=WM boundary, 1=CSF boundary).
 * @param dist_CSF_out Optional output: CSF distance map (can be NULL).
 * @param dist_WM_out  Optional output: WM distance map (can be NULL).
 * @param dims        Volume dimensions [nx, ny, nz].
 * @param voxelsize   Voxel sizes in mm [dx, dy, dz].
 * @param opts        Options controlling the algorithm behavior.
 *
 * @return 0 on success, non-zero on error.
 */
int CAT_VolComputePbt(
    const float *src,
    float *GMT_out,
    float *PPM_out,
    float *dist_CSF_out,
    float *dist_WM_out,
    int dims[3],
    double voxelsize[3],
    const CAT_PbtOptions *opts
);

#ifdef __cplusplus
}
#endif

#endif /* _CAT_VOL_PBT_H_ */
