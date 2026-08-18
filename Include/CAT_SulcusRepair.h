/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 *
 */

#ifndef _CAT_SULCUSREPAIR_H_
#define _CAT_SULCUSREPAIR_H_

/*
 * Anatomy-aware repair of a PVE label map, to be run *before* PBT.
 *
 * Three failure modes of a tissue classifier break central-surface extraction,
 * and all three are failures of *evidence*, not of smoothness:
 *
 *  1. Glued sulci.  Two banks of a tight sulcus are labelled as one thick GM
 *     band because no CSF was detected between them.  Typical in the occipital
 *     midline, where cortex is thin and the classifier has least to work with.
 *
 *  2. Broken gyri.  A thin gyral WM blade is interrupted by a small
 *     missegmentation, so the WM distance map -- and with it the thickness and
 *     the central surface -- is wrong along the whole blade.
 *
 *  3. Residual partial-volume error in the band where (1) happens: the label
 *     map says "no CSF anywhere near" while the T1 still shows an intensity
 *     dip across the sulcus.
 *
 * No amount of regularization fixes these, because regularization cannot
 * create evidence -- it can only redistribute what the classifier already
 * committed to.  What all three routines here have in common is that they go
 * back to the *intensity* image and recover evidence the classifier discarded,
 * using a Hessian sheetness filter (see CAT_Sheetness.h) as the shape prior.
 *
 * Every operation is one-sided and gated:
 *  - the CSF recovery can only ever lower a label towards what the intensity
 *    actually supports, never below it, and never within one voxel of the WM
 *    boundary (the "T > 1" guard of ACE);
 *  - the gyral reconnection only fires where WM is present on two opposite
 *    sides, so it cannot fire inside a sulcus;
 *  - the narrow-band refit only touches voxels that show a genuine local
 *    minimum across the sheet.
 *
 * Conventions: the label map is a PVE image in [0..3] with CSF = 1, GM = 2,
 * WM = 3, matching the input of CAT_VolComputePbt().  The intensity image is
 * expected to be bias corrected but may be in arbitrary units; it is
 * internally rescaled onto the same 1..3 axis.
 *
 * References
 * ----------
 * Han et al., Proc SPIE Med Imag 4322:194-203, 2001 (ACE).
 * Han et al., NeuroImage 23(3):997-1012, 2004 (CRUISE).
 * Kim et al., NeuroImage 27(1):210-221, 2005 (CLASP, CSF skeleton).
 * Descoteaux et al., Med Image Anal 10(4):638-651, 2006 (sheetness).
 */

/** \brief Parameters of the pre-PBT repair operations. */
typedef struct
{
    /* shared sheetness scales, in mm */
    double sheet_sigma_min; /**< smallest Gaussian scale */
    double sheet_sigma_max; /**< largest Gaussian scale */
    int sheet_n_scales;     /**< number of log-spaced scales */
    double sheet_strength;  /**< overall gain on the sheetness (default 1.0); passed
                                 through as CAT_SheetnessOpts::gain.  Raise it above 1
                                 when the response is too weak to clear csf_thresh /
                                 wm_thresh; 0 disables every sheetness-gated step. */

    /* (1) sulcal CSF recovery */
    double csf_min_dist;  /**< mm; act only where the label map sees no CSF within this radius */
    double csf_min_wmdist;/**< mm; never carve closer than this to the WM boundary */
    double csf_thresh;    /**< sheetness below this is ignored (default 0.1); the blend
                               ramps from 0 here to csf_strength at twice this value, so
                               it is half the sheetness at which the correction acts at
                               full strength -- the same relation the oriented median has
                               to its cutoff.  Match it to the response the data actually
                               produces; see CAT_Sheetness.h. */
    double csf_strength;  /**< 0..1 blend towards the intensity-implied label */

    /* (2) gyral WM reconnection */
    double wm_thresh;     /**< sheetness below this is ignored (default 0.1); ramps to
                               wm_strength at twice this value, as csf_thresh does */
    double wm_strength;   /**< 0..1 blend towards WM */
    double wm_min_int;    /**< intensity floor on the 1..3 axis (default 2.1); a blade tip
                               is dragged towards GM by partial volume, so this sits just
                               above pure GM rather than half way to WM */
    int wm_max_gap;       /**< how far from existing WM, in voxels, a blade may still be
                               strengthened (default 3) */

    /* (3) narrow-band PVE refit */
    double band_min_dist; /**< mm; refit only outside this distance to detected CSF */
    int band_window;      /**< half-width in voxels of the window for local class means */
    double band_strength; /**< 0..1 blend towards the refitted value */

    int verbose;
} CAT_SulcusRepairOpts;

/**
 * \brief Fill a CAT_SulcusRepairOpts with defaults tuned for 0.5 mm data.
 *
 * \param opts (out) option block to initialize; NULL is ignored
 */
void CAT_SulcusRepairOptionsInit(CAT_SulcusRepairOpts *opts);

/**
 * \brief Rescale an intensity image onto the 1..3 axis of a PVE label map.
 *
 * Estimates the CSF, GM and WM intensity levels as the medians of the voxels
 * that the label map calls pure tissue, then maps the intensity piecewise
 * linearly onto [1,3].  This makes every threshold in this module independent
 * of the units and the scaling of the input image.
 *
 * \param t1        (in)  intensity image, arbitrary units
 * \param label     (in)  PVE label image in [0..3]
 * \param t1n       (out) float[nvox] rescaled image on the 1..3 axis
 * \param dims      (in)  {nx, ny, nz}
 * \return 0 on success, -1 on invalid arguments, -3 if a class is unpopulated
 */
int CAT_VolNormalizeToLabel(const float *t1, const float *label, float *t1n,
                            int dims[3]);

/**
 * \brief Recover sulcal CSF that the classifier missed (opens glued sulci).
 *
 * Runs a dark-sheet (valley) sheetness filter on the intensity image and uses
 * it to pull GM labels down towards the intensity-implied value wherever the
 * label map reports no CSF nearby.  This is the mu_CSF evidence term of ACE,
 * recovered from the image rather than from the classifier: a T1 still shows
 * the intensity dip through a glued sulcus long after the label map has
 * committed to pure GM.
 *
 * The correction is strictly one-sided -- a label is only ever lowered, and
 * only as far as the intensity supports -- and is suppressed within
 * opts->csf_min_wmdist of the WM boundary so that thin cortex is not eaten.
 *
 * \param t1        (in)     bias-corrected intensity image
 * \param label     (in/out) PVE label image in [0..3], corrected in place
 * \param sheetness (out)    float[nvox] dark-sheet response, or NULL to skip
 * \param dims      (in)     {nx, ny, nz}
 * \param voxelsize (in)     voxel spacing in mm
 * \param opts      (in)     parameters; NULL selects the defaults
 * \return 0 on success, negative on error
 */
int CAT_VolRecoverSulcalCSF(const float *t1, float *label, float *sheetness,
                            int dims[3], double voxelsize[3],
                            const CAT_SulcusRepairOpts *opts);

/**
 * \brief Strengthen thin WM blades the classifier under-labelled.
 *
 * Targets the fine white-matter fingers reaching into the gyral crowns, which
 * are one to two voxels across at their far end.  Partial volume pulls their
 * intensity towards GM there, so a classifier that resolves the trunk of a
 * blade correctly still tends to lose its tip -- and losing the last millimetre
 * of a blade corrupts the WM distance map, and with it the thickness and the
 * central surface along the whole gyrus.
 *
 * Runs a bright-sheet (ridge) sheetness filter on the intensity image to find
 * thin bright structures, keeps the voxels that are labelled GM, brighter than
 * opts->wm_min_int, brighter than their own label, and within opts->wm_max_gap
 * voxels of existing WM, then accepts those a geodesic growth (downcut_float)
 * seeded from the existing WM can reach *through the candidate set itself*.
 * Growing through the candidates is what makes the repair follow the blade
 * instead of dilating in all directions, and it is also the connectivity
 * evidence: an isolated bright speck in GM is never reached, the tip of a real
 * blade is reached from the trunk behind it.
 *
 * The operation cannot fire inside a sulcus because a sulcus is a *dark* sheet
 * and the polarity guard of the sheetness filter rejects it outright; the
 * intensity floor and the one-sided blend bound the rest.
 *
 * Must run before PBT: repairing afterwards is too late, because dist_WM and
 * the thickness were already computed on the broken blade.
 *
 * \param t1        (in)     bias-corrected intensity image
 * \param label     (in/out) PVE label image in [0..3], corrected in place
 * \param sheetness (out)    float[nvox] bright-sheet response, or NULL to skip
 * \param dims      (in)     {nx, ny, nz}
 * \param voxelsize (in)     voxel spacing in mm
 * \param opts      (in)     parameters; NULL selects the defaults
 * \return 0 on success, negative on error
 */
int CAT_VolStrengthenWmBlades(const float *t1, float *label, float *sheetness,
                              int dims[3], double voxelsize[3],
                              const CAT_SulcusRepairOpts *opts);

/**
 * \brief Re-estimate partial volume in a narrow band from local class means.
 *
 * Targets the voxels where the label map claims there is no CSF within
 * opts->band_min_dist but the intensity image has a local minimum across the
 * sheet normal.  For those voxels the two bracketing class levels are
 * re-estimated from a local window -- which makes the refit robust against
 * residual bias that the global classifier could not remove -- and the label
 * is recomputed as the linear mixing fraction between them.
 *
 * \param t1        (in)     bias-corrected intensity image
 * \param label     (in/out) PVE label image in [0..3], corrected in place
 * \param dims      (in)     {nx, ny, nz}
 * \param voxelsize (in)     voxel spacing in mm
 * \param opts      (in)     parameters; NULL selects the defaults
 * \return 0 on success, negative on error
 */
int CAT_VolRefinePveNarrowBand(const float *t1, float *label,
                               int dims[3], double voxelsize[3],
                               const CAT_SulcusRepairOpts *opts);

#endif /* _CAT_SULCUSREPAIR_H_ */
