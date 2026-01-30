/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 *
 */

#ifndef _CAT_VOL_LAYER_SMOOTH_H_
#define _CAT_VOL_LAYER_SMOOTH_H_

#include <stdio.h>
#include <math.h>
#include <memory.h>
#include <float.h>

/**
 * smooth_within_cortex - Anisotropic smoothing along cortical layers
 *
 * This function performs Laplace-guided anisotropic smoothing that smooths
 * data along cortical layers (iso-depth contours) while preserving the
 * laminar structure. Smoothing is weighted by both spatial distance and
 * similarity in cortical depth (relative position between WM and CSF).
 *
 * The algorithm:
 * 1. Computes cortical depth field: phi = WMD / (WMD + CSFD)
 * 2. For each voxel, applies Gaussian smoothing weighted by:
 *    - Euclidean distance (standard Gaussian kernel)
 *    - Depth similarity: exp(-(delta_phi)^2 / (2*sigma_depth^2))
 * 3. Optionally extends smoothing slightly outside the cortical band
 *
 * Parameters:
 *  @param data       Input/output data volume to be smoothed (modified in-place)
 *  @param seg        PVE segmentation image (1=CSF, 2=GM, 3=WM, with partial volumes)
 *  @param dims       Volume dimensions [x, y, z]
 *  @param voxelsize  Voxel size in mm [x, y, z]
 *  @param fwhm       Smoothing kernel size in mm (FWHM along cortex)
 *  @param extend     Extension outside cortical band in mm (0 = strict GM only)
 *  @param datatype   NIfTI datatype of input data
 */
void smooth_within_cortex(void *data, float *seg, int dims[3], double voxelsize[3],
                          double fwhm, double extend, int datatype);

/**
 * smooth_within_cortex_float - Float version of layer smoothing
 *
 * Same as smooth_within_cortex but operates directly on float data.
 *
 * Parameters:
 *  @param data       Input/output float data volume (modified in-place)
 *  @param seg        PVE segmentation image
 *  @param dims       Volume dimensions [x, y, z]
 *  @param voxelsize  Voxel size in mm [x, y, z]
 *  @param fwhm       Smoothing kernel size in mm (FWHM along cortex)
 *  @param extend     Extension outside cortical band in mm
 */
void smooth_within_cortex_float(float *data, float *seg, int dims[3], double voxelsize[3],
                                double fwhm, double extend);

/**
 * compute_cortical_depth - Compute cortical depth field (Laplace-like)
 *
 * Computes the relative position within the cortical ribbon:
 *   phi = WMD / (WMD + CSFD)
 * where WMD is distance to WM boundary and CSFD is distance to CSF boundary.
 * Values range from 0 (at WM surface) to 1 (at pial/CSF surface).
 *
 * Parameters:
 *  @param seg        PVE segmentation image
 *  @param depth      Output depth field (must be pre-allocated)
 *  @param dist_WM    Output WM distance (can be NULL if not needed)
 *  @param dist_CSF   Output CSF distance (can be NULL if not needed)
 *  @param dims       Volume dimensions [x, y, z]
 *  @param voxelsize  Voxel size in mm [x, y, z]
 *  @param extend     Extension outside cortical band in mm
 */
void compute_cortical_depth(float *seg, float *depth, float *dist_WM, float *dist_CSF,
                            int dims[3], double voxelsize[3], double extend);

#endif /* _CAT_VOL_LAYER_SMOOTH_H_ */
