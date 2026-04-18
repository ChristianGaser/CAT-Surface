/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 */

/**
 * \file CAT_VolMyelinCorrection.h
 * \brief Correct PVE label bias at tissue boundaries using T1w intensities.
 *
 * In heavily myelinated cortical regions (primary motor/somatosensory cortex),
 * deeper cortical layers have high myelin content, making gray matter appear
 * brighter on T1w MRI — closer to white matter intensity.  This causes
 * tissue-classification algorithms to misclassify myelinated GM as WM,
 * pushing the GM/WM boundary outward and producing artificially thin
 * cortical thickness estimates.
 *
 * Analogously, at the pial (GM/CSF) boundary, partial-volume effects or
 * intensity biases can cause GM to be misclassified as CSF, pulling the
 * pial surface inward.
 *
 * The correction operates on the PVE label volume *before* projection-based
 * thickness (PBT) estimation.  It uses the raw T1w intensity image to detect
 * voxels near both tissue boundaries whose intensity and gradient
 * characteristics are inconsistent with their labelled tissue class:
 *
 * **WM/GM boundary** (myelination correction):
 *  1. Erode the WM mask to obtain a "deep WM" core.
 *  2. Compute intensity statistics (mean, std) from the T1w within that core.
 *  3. At each WM-boundary voxel (PVE ~ 2.3-3.0), flag it as suspected
 *     myelination artifact when *all three* conditions hold:
 *       a) T1w intensity is below the deep-WM range,
 *       b) T1w gradient magnitude is weak (no clear tissue boundary),
 *       c) Euclidean distance to the deep-WM core exceeds a threshold.
 *  4. For flagged voxels, shift PVE values toward GM.
 *
 * **GM/CSF boundary** (pial correction):
 *  1. Erode the CSF mask to obtain a "deep CSF" core.
 *  2. Compute intensity statistics (mean, std) from the T1w within that core.
 *  3. At each CSF-boundary voxel (PVE ~ 1.0-1.7), flag it as suspected
 *     GM misclassified as CSF when *all three* conditions hold:
 *       a) T1w intensity is above the deep-CSF range,
 *       b) T1w gradient magnitude is weak,
 *       c) Euclidean distance to the deep-CSF core exceeds a threshold.
 *  4. For flagged voxels, shift PVE values toward GM.
 *
 * Both correction fields are median-filtered before application to suppress
 * noisy, isolated corrections while preserving coherent boundary shifts.
 */

#ifndef _CAT_VOL_MYELIN_CORRECTION_H_
#define _CAT_VOL_MYELIN_CORRECTION_H_

#include "CAT_Vol.h"
#include "CAT_NiftiLib.h"

/**
 * \brief Options controlling boundary correction behavior.
 */
typedef struct
{
    double erosion_mm;       /**< Erosion radius (mm) for deep tissue cores (default: 3.0) */
    double k_intensity;      /**< Std-dev factor for intensity threshold (default: 1.0) */
    double grad_percentile;  /**< Percentile of boundary gradient below which gradient is "weak" (default: 25.0) */
    double dist_mm;          /**< Min distance (mm) from tissue core to be suspicious (default: 2.0) */
    double max_correction;   /**< Maximum PVE shift toward GM (default: 0.8) */
    double min_cluster_mm3;  /**< Minimum flagged-cluster volume in mm^3 (default: 20.0). Clusters smaller than this are discarded to prevent noisy, isolated corrections. */
    double gm_grad_pct;      /**< Percentile of GM gradient magnitude above which gradient is "high" (default: 50.0). Used for double-gradient myelination detection in motor cortex. */
    double max_gm_grad_dist; /**< Max distance (mm) from high-gradient GM for double-gradient detection (default: 3.0). WM-boundary voxels within this distance of elevated-gradient GM tissue are flagged as likely myelinated GM. Set to 0 to disable. */
    int n_median_filter;     /**< Iterations of 3x3x3 median filter on the correction field (default: 2). Removes noisy isolated corrections while preserving coherent boundary shifts. Set to 0 to disable. */
    int correct_wm;          /**< Apply WM/GM (myelination) correction (default: 1) */
    int correct_csf;         /**< Apply GM/CSF (pial) correction (default: 1) */
    int verbose;             /**< Print progress messages (default: 0) */
} CAT_MyelinCorrOptions;

/**
 * \brief Initialize myelination correction options with defaults.
 *
 * \param opts (out) options structure to fill
 */
void CAT_MyelinCorrOptionsInit(CAT_MyelinCorrOptions *opts);

/**
 * \brief Correct PVE labels at both tissue boundaries using T1w intensities.
 *
 * Uses the raw T1w intensity image to identify voxels near the GM/WM
 * and/or GM/CSF boundaries that were likely misclassified, and shifts
 * their PVE values toward GM.
 *
 * At the WM side, this detects myelination-induced WM over-labelling.
 * At the CSF side, this detects GM voxels misclassified as CSF.
 *
 * Which boundaries are corrected is controlled by opts->correct_wm
 * and opts->correct_csf.
 *
 * The correction is applied in-place to the PVE label volume.
 *
 * \param pve        (in/out) PVE label volume (float, values in [0..3])
 * \param t1w        (in)     raw T1w intensity volume (float)
 * \param dims       (in)     volume dimensions {nx, ny, nz}
 * \param voxelsize  (in)     voxel sizes in mm {dx, dy, dz}
 * \param opts       (in)     algorithm options (NULL for defaults)
 * \return 0 on success, non-zero on error
 */
int CAT_VolCorrectMyelination(float *pve, const float *t1w,
                              int dims[3], double voxelsize[3],
                              const CAT_MyelinCorrOptions *opts);

#endif /* _CAT_VOL_MYELIN_CORRECTION_H_ */
