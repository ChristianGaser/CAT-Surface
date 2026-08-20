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

#include "CAT_Sheetness.h"

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
    double sheet_normalize; /**< value the response's p99.9 is scaled to (default
                                 CAT_SHEETNESS_NORMALIZE = 1.0); <= 0 keeps the raw
                                 response.  This is what makes csf_thresh / wm_thresh
                                 mean the same thing on every dataset; disable it only
                                 where the image may contain no sheets at all, since a
                                 percentile anchor would then amplify noise. */
    int sheet_skeleton;     /**< 1 to thin the response to its medial sheet before
                                 any threshold is applied (default 0); passed through as
                                 CAT_SheetnessOpts::skeletonize.  The plate response is
                                 as wide as the Gaussian that produced it, so the scale
                                 that finds a blade or a fundus best also answers several
                                 voxels into the tissue on either side, and the per-voxel
                                 gates below then correct that surrounding tissue too.
                                 Enable it when a large sheet_sigma_max locates the
                                 structures well but the corrections are too broad. */
    double sheet_strength;  /**< overall gain on the sheetness (default 1.0); passed
                                 through as CAT_SheetnessOpts::gain.  Since the response
                                 is now anchored to its own p99.9 the gain is normally
                                 unnecessary -- it used to absorb the per-dataset scale
                                 that the anchor removes, which is why values of 20 or
                                 more were once needed.  It remains as a deliberate
                                 relative adjustment; 0 disables every sheetness-gated
                                 step. */

    /* (1) sulcal CSF recovery */
    double csf_min_dist;  /**< mm; act only where the label map sees no CSF within this radius */
    double csf_min_wmdist;/**< mm; never carve closer than this to the WM boundary */
    double csf_thresh;    /**< sheetness below this is ignored (default 0.3); the blend
                               ramps from 0 here to csf_strength at twice this value, so
                               it is half the sheetness at which the correction acts at
                               full strength -- the same relation the oriented median has
                               to its cutoff.  The response is anchored to its own p99.9
                               by CAT_SheetnessOpts::normalize, so 0.3 means the same
                               thing on every dataset: full strength from 0.6, which is
                               where the p99 of the reference response sits.  See
                               CAT_Sheetness.h. */
    double csf_strength;  /**< 0..1 blend towards the intensity-implied label */

    /* (2) gyral WM reconnection */
    double wm_thresh;     /**< sheetness below this is ignored (default 0.3); ramps to
                               wm_strength at twice this value, as csf_thresh does */
    double wm_strength;   /**< 0..1 blend towards WM */
    double wm_min_int;    /**< intensity floor on the 1..3 axis (default 2.1); a blade tip
                               is dragged towards GM by partial volume, so this sits just
                               above pure GM rather than half way to WM */
    int wm_max_gap;       /**< how far from existing WM, in voxels, a blade may still be
                               strengthened (default 3) */
    double wm_sulcus_guard; /**< 0..1, default 1.0; how strongly a nearby sulcus vetoes
                               the blade strengthening.  A blade tip and the sulcal
                               floor behind it are one voxel apart where the cortex is
                               thin, so raising the tip towards WM closes the sulcus --
                               the failure is common in the occipital lobe, where the
                               banks are already almost touching.  The polarity guard
                               alone does not catch it: the bright ridge is genuinely
                               there, it is the *neighbouring* dark sheet that must not
                               be filled.  A second, dark-sheet pass is therefore run
                               and dilated by one voxel, and the blend weight is damped
                               by (1 - guard * ramp(dark)), reaching zero at twice
                               csf_thresh.  Set 0 to disable the guard entirely. */

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

/*
 * Opening buried sulci directly in a percentage position map (PPM).
 *
 * The pre-PBT repairs above need the intensity image. By the time the central
 * surface is extracted the T1 is gone -- marching cubes sees only the PPM --
 * but the PPM carries the geometry itself, so no intensity is required:
 *
 *   crossing a sulcus:  1 (WM) -> 0.5 -> ~0 (pial) -> 0.5 -> 1 (other bank)
 *   crossing a gyrus:   0 (pial) -> 0.5 -> ~1 (WM blade) -> 0.5 -> 0
 *
 * A sulcus is therefore a *valley* in the PPM and a gyral blade a *ridge*, and
 * the polarity guard of the sheetness filter separates them exactly. A buried
 * sulcus is a valley whose floor never drops below the isovalue: the geometry
 * is still there, only the amplitude is missing, which is why an isovalue-based
 * surface fuses the two banks while the valley remains plainly visible to a
 * shape filter.
 */

/** \brief Parameters for opening buried sulci in a PPM. */
typedef struct
{
    double sigma_min; /**< smallest sheetness scale in mm */
    double sigma_max; /**< largest sheetness scale in mm */
    int n_scales;     /**< number of log-spaced scales */
    double sheet_normalize; /**< value the response's p99.9 is scaled to (default
                                 CAT_SHEETNESS_NORMALIZE); <= 0 keeps it raw */
    int sheet_skeleton;    /**< 1 to thin the response to its medial sheet (default 0);
                                see CAT_SulcusRepairOpts::sheet_skeleton */
    double sheet_strength; /**< overall gain on the response (default 1.0), passed
                                through as CAT_SheetnessOpts::gain.  The automatic
                                noise scale of the sheetness filter is half the
                                largest Hessian norm in the volume, which on real
                                data is set by the cortical ribbon itself -- a thin
                                sulcal valley is far weaker than that, so the raw
                                response usually sits an order of magnitude below
                                the thresholds here and nothing happens at a gain
                                of 1.  This is the same reason CAT_VolSulcusRepair
                                needs its own -sheet-strength, and the value does
                                not carry over: that one is measured on the
                                intensity image, this one on the PPM.  Measure with
                                CAT_VolSheetness -polarity -1 on the PPM and choose
                                the gain that puts p99 near twice thresh. */
    double thresh;    /**< sheetness below this is ignored */
    double band;      /**< only act on values in (isovalue, isovalue + band) */
    double margin;    /**< how far below the isovalue an opened valley is pushed */
    double strength;  /**< 0..1 blend towards that target */
    double cutoff;    /**< admission cutoff of the oriented median the caller may run
                           alongside; <= 0 selects CAT_ORIENTED_MEDIAN_CUTOFF */
    int verbose;
} CAT_PpmSulciOpts;

/**
 * \brief Fill a CAT_PpmSulciOpts with defaults for a 0..1 PPM at 0.5 mm.
 *
 * \param opts (out) option block to initialize; NULL is ignored
 */
void CAT_PpmSulciOptionsInit(CAT_PpmSulciOpts *opts);

/**
 * \brief Push buried sulcal valleys in a PPM below the isovalue.
 *
 * Finds thin low-value sheets -- valleys -- in the map and lowers the ones
 * that sit just above the isovalue until they cross it, so the isosurface
 * separates the two banks instead of bridging them.
 *
 * Three conditions must hold together, and it is the conjunction that makes
 * the operation safe rather than just another filter:
 *
 *   - the value lies in (isovalue, isovalue + band), so only marginal cases
 *     are touched and solid white matter (PPM near 1) never is;
 *   - the dark-sheet response exceeds opts->thresh, so there really is a thin
 *     planar structure and not a blob or noise;
 *   - the value is a genuine local minimum *along the sheet normal*, which a
 *     gyral crown can never satisfy -- crossing a blade the PPM has a maximum,
 *     not a minimum, so the operation cannot cut a gyrus open.
 *
 * The correction is one-sided: it only ever lowers a value.
 *
 * \param ppm       (in/out) percentage position map in [0,1], corrected in place
 * \param sheetness (in)     precomputed dark-sheet response, or NULL to compute
 * \param normal    (in)     precomputed unit sheet normals, or NULL to compute
 * \param dims      (in)     {nx, ny, nz}
 * \param voxelsize (in)     voxel spacing in mm
 * \param isovalue  (in)     threshold the surface will be extracted at
 * \param opts      (in)     parameters; NULL selects the defaults
 * \return 0 on success, negative on error
 */
int CAT_VolOpenPpmSulci(float *ppm, const float *sheetness, const float *normal,
                        int dims[3], double voxelsize[3], double isovalue,
                        const CAT_PpmSulciOpts *opts);

#endif /* _CAT_SULCUSREPAIR_H_ */
