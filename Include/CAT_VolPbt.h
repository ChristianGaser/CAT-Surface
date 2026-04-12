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
extern "C"
{
#endif

    /**
     * Options structure for controlling PBT estimation behavior.
     */
    typedef struct
    {
        int n_avgs;               /**< Number of averages for distance estimation (default: 2) */
        int n_median_filter;      /**< Iterations of local median filtering for PPM cleanup (default: 2). The filter is blended in only where a topology-artifact weight map is high: the weight is estimated from the positive residual PPM - smooth(PPM), restricted to sufficiently thick cortex, cleaned morphologically, then smoothed before mixing original and median-filtered PPM. Set to 0 to disable. */
        int median_subsample;     /**< Subsampling for median filter to smooth thickness values */
        double range;             /**< Extended range for masking euclidean distance (default: 0.45) */
        double fill_thresh;       /**< Threshold for filling holes in PPM (default: 0.5, 0=disable) */
        double correct_voxelsize; /**< Amount of thickness correction for voxel-size (default: 0.5) */
        double sulcal_width;      /**< Max distance from CSF boundary (mm) for sulcal PPM correction (default: 2.5, 0=disable) */
        int fast;                 /**< Fast mode: rougher but quicker estimate (default: 0) */
        int verbose;              /**< Verbose output (default: 0) */
    } CAT_PbtOptions;

    /**
     * \brief Blood-vessel correction for float PVE label maps.
     *
     * Implements blood-vessel correction on PVE labels in [0..3] using downcut growing,
     * ring-constrained median inpainting, and class-aware value clamping.
     *
     * \param Yp0              (in/out) float PVE label map
     * \param dims             (in)     dimensions {nx, ny, nz}
     * \param vx_vol           (in)     voxel spacing {sx, sy, sz}; NULL -> {1,1,1}
     */
    void blood_vessel_correction_pve_float(float *Yp0, int dims[3], double vx_vol[3]);

    /**
     * \brief Datatype-generic wrapper for PVE blood-vessel correction.
     *
     * Converts input volume to float, applies `blood_vessel_correction_pve_float()`, then converts back
     * to the specified datatype.
     *
     * \param data             (in/out) PVE label volume
     * \param dims             (in)     dimensions {nx, ny, nz}
     * \param vx_vol           (in)     voxel spacing {sx, sy, sz}; NULL -> {1,1,1}
     * \param datatype         (in)     datatype code (DT_UINT8, DT_UINT16, DT_FLOAT32, etc.)
     */
    void blood_vessel_correction_pve(void *data, int dims[3], double vx_vol[3], int datatype);

    /**
     * Initialize PBT options with default values.
     *
     * @param opts Pointer to options structure to initialize.
     */
    void CAT_PbtOptionsInit(CAT_PbtOptions *opts);

    /**
     * \brief Compute projection-based cortical thickness and percentage position map.
     *
     * This function implements the complete PBT pipeline:
     * 1. Distance estimation in WM and CSF by shifting GM borders
     * 2. Thickness estimation via projection_based_thickness()
     * 3. PPM calculation with gyrus/sulcus blending
     * 4. Optional weighted local median filtering for topology artifact reduction
     *
     * If opts->n_median_filter > 0, the final PPM cleanup is applied only
     * where a topology-artifact likelihood map is high rather than everywhere
     * uniformly. This likelihood is estimated from the positive residual between
     * the PPM and a smoothed PPM, restricted to thicker cortex (GMT > 1.5), then
     * regularized by close/open/dilate operations and smoothing. The resulting
     * soft weight in [0,1] blends the original PPM with a locally median-filtered
     * PPM, so stronger suspected artifacts receive more median-filter influence.
     *
     * @param src         Input PVE label image (CSF=1, GM=2, WM=3 with partial volumes).
     * @param GMT_out     Output: Gray matter thickness map (caller must allocate nvox floats).
     * @param PPM_out     Output: Percentage position map (0=WM boundary, 1=CSF boundary).
     * @param dist_CSF_out Optional output: CSF distance map (can be NULL).
     * @param dist_WM_out  Optional output: WM distance map (can be NULL).
     * @param dims        Volume dimensions [nx, ny, nz].
     * @param voxelsize   Voxel sizes in mm [dx, dy, dz].
     * @param opts        Options controlling the algorithm behavior, including
     *                    n_median_filter for the weighted local PPM cleanup.
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
        const CAT_PbtOptions *opts);

    /**
     * Calculate the thickness of segmented structures using projection-based method.
     *
     * This function estimates the projection-based thickness of segmented structures
     * in a 3D volume, typically used in medical imaging for brain tissue analysis.
     * It requires PVE label images and distance maps for WM and CSF.
     *
     * @param SEG      PVE label image with labels for CSF (1), GM (2), and WM (3).
     * @param WMD      White Matter distance map.
     * @param CSFD     Cerebrospinal Fluid distance map.
     * @param GMT      Output thickness image.
     * @param dims     Array containing the dimensions of the volume [nx, ny, nz].
     * @param voxelsize Array containing the size of each voxel in mm [dx, dy, dz].
     */
    void projection_based_thickness(float *SEG, float *WMD, float *CSFD, float *GMT, int dims[3], double *voxelsize);

    /**
     * Build a smooth "gyri/sulci mask" (gyri ≈ 0, sulci ≈ 1).
     *
     * This mask supports downstream operations that benefit from different behavior in
     * gyral vs. sulcal regions, e.g.:
     *  - Surface extraction: prevent sulcal closure by using a higher isovalue in sulci,
     *    and prevent cutting gyri by using a lower isovalue in gyri.
     *  - Projection-based cortical thickness: use different parameters over gyri vs. sulci.
     *
     * Algorithm (overview):
     *  1. Initial thresholding of the input scalar image at (thresh * max(src)) to form
     *     a coarse tissue mask; then distance-based closing to fill sulcal gaps.
     *  2. Gyri emphasis: a slight dilation followed by a stronger erosion to
     *     preferentially shrink gyri crowns relative to sulcal regions.
     *  3. Smoothing to create a soft (0..1) transition between gyri and sulci.
     *  4. CSF enforcement: voxels clearly below tissue threshold (e.g., CSF) are set to 1
     *     (open sulci), followed by a light final smoothing to avoid hard borders.
     *
     * Conventions:
     *  - Output mask is in [0,1] (float). Values close to 0 indicate gyri crowns;
     *    values close to 1 indicate sulcal fundi.
     *
     * @param src       Input scalar image (e.g., label map float[nx*ny*nz])
     * @param mask      Output mask (0..1) pre-allocated by caller float[nx*ny*nz]
     * @param dims      Volume dimensions {nx, ny, nz}
     * @param voxelsize Voxel spacing in mm {sx, sy, sz}
     * @param thresh    Threshold for src; seeds the initial mask (recommended: 1.5 for label map)
     * @param fwhm      Smoothing FWHM for the main blur step (recommended: 8.0)
     *
     * @note The heuristic constants (closing=5, dilate=2, erode=5, CSF_factor=0.75, final_FWHM=2)
     *       follow the original intent and can be exposed as parameters if needed.
     */
    void smooth_gyri_mask(const float *src, float *mask, int dims[3], double voxelsize[3],
                          double thresh, double fwhm);

#ifdef __cplusplus
}
#endif

#endif /* _CAT_VOL_PBT_H_ */
