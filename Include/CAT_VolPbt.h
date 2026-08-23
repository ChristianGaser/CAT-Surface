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
        double correct_thickness; /**< Additive thickness correction in **mm** (default: 0.25). */
        double sulcal_width;      /**< Max distance from CSF boundary (mm) for sulcal PPM correction (default: 2.5, 0=disable) */
        int pve_distance;         /**< Sub-voxel PVE correction of the WM/CSF distance maps (default: 0=off). The distance transform measures to the nearest source voxel *centre*; with this enabled the partial volume of that source voxel is used to measure to the tissue boundary instead, as in CAT12's cat_vol_pbtsimpleCS4. EXPERIMENTAL: it shifts the thickness calibration, so correct_thickness must be re-derived when enabling it. */
        int sulcal_barrier;       /**< Stop the CSF distance at the sulcal medial surface (default: 0=off). PBT overestimates thickness wherever the classifier lost the CSF in a sulcus: with no CSF boundary to stop at, the distance from one bank propagates through the fused grey matter into the other, and the thickness -- and with it the PPM, and with it the surface -- follows it across. This bounds dist_CSF by the distance to the sulcal midline, located by front collision from geometry alone (CAT_VolSulcalMedialSet), so no intensity image and no per-subject threshold are involved. The bound is applied as a minimum, which makes it self-limiting: where CSF *was* segmented it is nearer than the midline and nothing changes at all, and the correction can only ever reduce an overestimated thickness, never inflate one. */
        double barrier_q;         /**< Shock threshold of the medial set (default: 0, which selects CAT_MEDIAL_SET_Q). Lower is stricter and gives a thinner midline. */
        double barrier_tmin;      /**< Ignore collisions closer than this to the WM boundary, in mm (default: 0.5). ACE's "T > 1" guard; it is what stops the barrier from carving into thin cortex from the inside. */
        double barrier_halfwidth; /**< Half the width the barrier stands for, in mm (default: 0, which selects one voxel). The distance transform measures to the nearest medial *voxel centre* while the grey matter ends at the surface of that set, so the raw distance is short by half the set's width -- and the shock test admits a band about two voxels across, not one, which is why the correction is a whole voxel rather than half. Measured rather than assumed: on a fused-bank phantom one voxel reproduces the thickness of the correctly segmented sulcus exactly (2.19 mm), where half a voxel leaves it 0.25 mm long. It does not move the surface either way -- the 0.5 crossing is identical across the range -- so it trades against thickness accuracy only. */
        int oriented_filter;      /**< Replace the isotropic 3x3x3 medians by sheetness-oriented ones (default: 0=off). An isotropic median penalizes boundary area and therefore removes thin structures whichever side of the label boundary they sit on -- it opens a glued sulcus and closes a cerebellar fissure at the same rate. The oriented variant admits only neighbours lying *in* the plane of the locally detected sheet, so it can no longer close one; where no sheet is detected it is identical to the isotropic median. See CAT_Sheetness.h. */
        double oriented_cutoff;   /**< Admission cutoff of the oriented medians (default: 0, which selects CAT_ORIENTED_MEDIAN_CUTOFF). A neighbour is admitted when sheetness*(dhat.n)^2 < cutoff, so a one-voxel-thick sheet is preserved from a sheetness of 2*cutoff upwards. See CAT_Sheetness.h. */
        double oriented_strength; /**< Overall gain on the sheetness before the oriented filters use it (default: 1.0). Passed through as CAT_SheetnessOpts::gain, so 0 reproduces the isotropic filters exactly and values above 1 amplify a response too weak to cross the hard 0.5 threshold of CAT_VolOrientedMedian(). See CAT_Sheetness.h. */
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
