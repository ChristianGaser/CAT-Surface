/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 *
 */

#ifndef _CAT_SHEETNESS_H_
#define _CAT_SHEETNESS_H_

/*
 * Hessian-based detection of thin sheet-like ("plate") structures and the
 * oriented filters built on top of it.
 *
 * Motivation
 * ----------
 * Isotropic regularizers -- a local median, a Potts MRF, total variation --
 * penalize boundary area.  A thin structure has an extreme area-to-volume
 * ratio, so deleting it is always the cheaper labelling ("shrinking bias").
 * That is why one and the same median filter opens glued sulci in one place
 * and closes cerebellar fissures in another: to an isotropic prior the two
 * are the same object and only the sign of the label differs.
 *
 * The escape is to stop using a smoothness prior and use a *shape* prior.
 * From the eigenvalues of the Hessian, |l1| <= |l2| <= |l3|, the local shape
 * is read off directly:
 *
 *   sheet (plate):  |l1| ~ |l2| ~ 0,  |l3| large
 *   tube  (line):   |l1| ~ 0,         |l2| ~ |l3| large
 *   blob:           |l1| ~ |l2| ~ |l3| large
 *
 * A sheetness filter therefore keeps thin sheets and ignores blobs, and it
 * shrinks nothing, because it makes no statement about boundary length.  The
 * eigenvector belonging to l3 is the sheet *normal*, which is what turns an
 * isotropic filter into an oriented one.
 *
 * Sulcal CSF is a dark sheet in a T1 (polarity -1), a gyral white-matter
 * blade is a bright sheet (polarity +1); the same operator finds both.
 *
 * References
 * ----------
 * Sato et al., Med Image Anal 2(2):143-168, 1998.
 * Frangi et al., MICCAI 1998, LNCS 1496:130-137.
 * Descoteaux et al., Med Image Anal 10(4):638-651, 2006 (sheetness form).
 */

/** \brief Parameters of the multi-scale sheetness filter. */
typedef struct
{
    double sigma_min; /**< smallest Gaussian scale in mm */
    double sigma_max; /**< largest Gaussian scale in mm */
    int n_scales;     /**< number of log-spaced scales in [sigma_min, sigma_max] */
    double alpha;     /**< plate-vs-tube sensitivity (R_sheet term) */
    double beta;      /**< blob-vs-plate sensitivity (R_blob term) */
    double c;         /**< structure-vs-noise sensitivity; <= 0 selects automatic */
    double gain;      /**< overall gain on the response (default 1.0); see below.
                           With normalize enabled this is rarely needed: the
                           anchor already puts the response on a fixed scale, so
                           the gain becomes a deliberate relative adjustment
                           rather than a per-dataset calibration. */
    double normalize; /**< scale the response so its p99.9 equals this value
                           (default CAT_SHEETNESS_NORMALIZE = 1.0); <= 0 keeps
                           the raw response.
                           The automatic noise scale c is half the largest
                           Hessian norm in the volume, so the absolute level of
                           the response depends on whatever the strongest
                           structure in that image happens to be.  Anchoring to
                           the map's own p99.9 makes every threshold downstream
                           -- CAT_ORIENTED_MEDIAN_CUTOFF, csf_thresh, wm_thresh
                           -- mean the same thing on every dataset.  Scaling is
                           a single positive factor, so the ranking of voxels,
                           the winning scale and the zero set are unchanged, and
                           the s = 0 invariant still holds exactly.  The hazard
                           is the mirror of the one it fixes: an image with no
                           sheets at all still has a p99.9, and normalizing to
                           it amplifies noise -- disable it there. */
    int skeletonize;  /**< 1 to keep only the medial sheet (default 0).
                           The response of a plate filter at scale sigma is as
                           wide as the Gaussian that produced it, so the scale
                           that detects a structure best also spreads its answer
                           several voxels into the tissue on either side.  Every
                           consumer gates per voxel, so that width is read as
                           "part of a sheet" well beyond the sheet itself and a
                           correction reaches into what surrounds it -- which is
                           what makes a large sigma_max both the most sensitive
                           and the least precise setting.
                           Suppressing everything that is not a maximum along
                           its own normal collapses the band onto its ridge line,
                           a one-voxel medial sheet, leaving the value on the
                           ridge unchanged.  It is Canny's non-maximum
                           suppression lifted to the sheet normal, and it gives
                           the medial surface far more cheaply than a
                           morphological thinning would.
                           Enable it when a large sigma_max finds the structures
                           but the corrections bleed into neighbouring tissue.
                           It runs before the percentile anchor, so p99.9 is
                           then taken over the ridge values. */
    int polarity;     /**< +1 bright sheets, -1 dark sheets, 0 either */
    int verbose;      /**< 1 to report progress */
} CAT_SheetnessOpts;

/*
 * About the gain
 * --------------
 * Every consumer of the sheetness gates on it through a threshold, and the
 * one in CAT_VolOrientedMedian() is hard: a neighbour is admitted unless
 * 1 - s*(dhat.n)^2 > 0.5 fails, and since (dhat.n)^2 <= 1 that cannot happen
 * for s <= 0.5.  A response that peaks at 0.5 therefore leaves the oriented
 * median bit-identical to the isotropic one -- not weakened, but unchanged.
 *
 * Crossing 0.5 is only where it begins to differ.  Preserving a one-voxel-thick
 * sheet takes more: of the 27 neighbours only the 9 in the sheet plane carry the
 * sheet value, and the 8 edge neighbours at (dhat.n)^2 = 1/2 are excluded only
 * once s reaches 1, so below saturation the admitted off-plane neighbours
 * outvote the sheet and the median erases it anyway.  On a 3x3x3 neighbourhood
 * the operator is therefore close to binary: useful at s = 1, inert below it.
 *
 * The automatic noise scale is what usually puts it there.  It is half the
 * largest Hessian norm found anywhere in the volume, so the noise factor at a
 * voxel carrying a fraction r of that maximum is 1 - exp(-2r^2): a voxel needs
 * r > 0.59 before that term alone reaches 0.5.  On a real head the maximum is
 * set by the scalp/air step, which no sulcal dip comes close to, so the
 * cortical response collapses towards zero.  Setting `c` explicitly is the
 * principled fix; the gain is the blunt one, and it is useful because it needs
 * no knowledge of the intensity units.
 *
 * The gain multiplies the response and the result is clamped to [0,1].  It is
 * a linear map fixing zero, so a zero response stays zero at any gain and the
 * invariant that every oriented operator degenerates to its isotropic
 * counterpart where no sheet was found is preserved exactly.  A gain above 1
 * also amplifies whatever noise floor the map has, so raise it while watching
 * the map from CAT_VolSheetness rather than blind.
 */

/**
 * \brief Fill a CAT_SheetnessOpts with defaults suitable for 0.5 mm brain data.
 *
 * \param opts (out) option block to initialize; NULL is ignored
 */
void CAT_SheetnessOptionsInit(CAT_SheetnessOpts *opts);

/**
 * \brief Multi-scale Hessian sheetness (plate) filter.
 *
 * \param src        (in)  scalar volume, e.g. a bias-corrected T1
 * \param sheetness  (out) float[nvox] response in [0,1]; maximum over scales
 * \param normal     (out) float[3*nvox] unit sheet normal, or NULL to skip
 * \param mask       (in)  optional uint8 mask; voxels with 0 are skipped, NULL = all
 * \param dims       (in)  {nx, ny, nz}
 * \param voxelsize  (in)  voxel spacing in mm {dx, dy, dz}
 * \param opts       (in)  filter parameters; NULL selects the defaults
 * \return 0 on success, -1 on invalid arguments, -2 on allocation failure
 */
int CAT_VolSheetness(const float *src, float *sheetness, float *normal,
                     const unsigned char *mask, int dims[3], double voxelsize[3],
                     const CAT_SheetnessOpts *opts);

/** \brief Value the response's p99.9 is scaled to by default; 0 keeps it raw. */
#define CAT_SHEETNESS_NORMALIZE 1.0

/** \brief Default admission cutoff of CAT_VolOrientedMedian(). */
#define CAT_ORIENTED_MEDIAN_CUTOFF 0.10

/**
 * \brief Median filter over a sheetness-oriented neighbourhood.
 *
 * Drop-in replacement for an isotropic 3x3x3 median that cannot close a thin
 * sheet.  A neighbour at offset d is admitted only when
 *
 *     s * (dhat . n)^2 < cutoff
 *
 * with s the local sheetness and n the local sheet normal.  For s = 0 the left
 * side is zero for every offset, so every neighbour is admitted and the
 * operator is bit-identical to the plain isotropic median at any cutoff --
 * behaviour away from thin structures never changes.
 *
 * How to choose the cutoff
 * -----------------------
 * Sort the 26 offsets by (dhat . n)^2: the 6 face neighbours sit at 1, the 12
 * edge neighbours at 1/2 and the 8 corners at 1/3, while the 9 offsets lying in
 * the sheet plane sit at 0 and are therefore admitted always.  A neighbour
 * class at (dhat . n)^2 = q drops out once s >= cutoff / q, so:
 *
 *   s >= cutoff       the 6 face neighbours drop out  (25 of 27 admitted)
 *   s >= 2 * cutoff   the 12 edge neighbours follow   (17 of 27 admitted)
 *   s >= 3 * cutoff   only the sheet plane survives   ( 9 of 27 admitted)
 *
 * The middle row is the one that matters.  A one-voxel-thick sheet puts the
 * sheet value on only the 9 in-plane offsets, so it survives the median exactly
 * when those 9 are at least half of what is admitted -- that is, from
 * s >= 2 * cutoff onwards.  Below it the filter differs from the isotropic
 * median without yet being able to keep a thin sheet.
 *
 * The cutoff is therefore half the sheetness at which a thin sheet starts being
 * preserved, and it should be set from the response the sheetness filter
 * actually produces on the data at hand.  The historical value was 0.5, which
 * put that point at s = 1 -- reachable only where the response saturates, so
 * the filter behaved as an on/off switch confined to a thin rim.
 *
 * With CAT_SheetnessOpts::normalize left at its default the response is scaled
 * so its p99.9 is 1, which is what makes the cutoff data-independent: it is now
 * read directly as a fraction of that anchor, and the same number means the
 * same thing on every dataset.
 *
 * The default of 0.30 is that fraction, derived from the reference response on
 * a 0.5 mm skull-stripped MPRAGE, whose raw dark-sheet map had p99 = 0.20 and
 * p99.9 = 0.33.  Full sheet preservation should begin at the p99 level, i.e. at
 * 0.20 / 0.33 = 0.61 of the anchor, and the cutoff is half of that -- so 0.30.
 * The top percent of the response is preserved outright, the few percent below
 * it get partial anisotropy, and everything else is left isotropic; on the
 * reference scan that changed about 2.6% of brain voxels against the plain
 * median.  Because the anchor now travels with the data, other scans reproduce
 * that behaviour at the same 0.30 rather than needing a new value each time.
 *
 * Reading it the other way round: a cutoff of 0.5 would start preservation
 * exactly at p99.9 (the top 0.1%), and the historical 0.5 *without*
 * normalization put it at a raw s = 1 -- reachable only where the response
 * saturates, so the filter behaved as an on/off switch confined to a thin rim.
 * That failure is what the anchor removes.
 *
 * \param vol       (in/out) volume to filter, in place
 * \param sheetness (in)     float[nvox] in [0,1]; NULL falls back to isotropic
 * \param normal    (in)     float[3*nvox] unit sheet normals; NULL falls back
 * \param mask      (in)     optional uint8 mask; NULL = filter everywhere
 * \param dims      (in)     {nx, ny, nz}
 * \param cutoff    (in)     admission cutoff; <= 0 selects
 *                           CAT_ORIENTED_MEDIAN_CUTOFF
 * \param iters     (in)     number of successive passes
 * \return 0 on success, -1 on invalid arguments, -2 on allocation failure
 */
int CAT_VolOrientedMedian(float *vol, const float *sheetness, const float *normal,
                          const unsigned char *mask, int dims[3], double cutoff,
                          int iters);

/** \brief Default shock threshold of CAT_VolSulcalMedialSet(). */
#define CAT_MEDIAL_SET_Q 0.6

/**
 * \brief Locate sulcal midlines by front collision, from geometry alone.
 *
 * The Hessian filters above need the structure to be *visible*: a sulcus has to
 * leave a dark sheet in the intensity, or a valley in the PPM, before they can
 * answer.  Where a classifier lost the CSF completely there is nothing left to
 * see, and on those subjects every intensity-derived threshold has to be
 * re-tuned -- which is the failure mode this routine exists to avoid.
 *
 * It uses no intensity at all.  Take the distance to the white-matter boundary,
 * propagated through the cortical band, and follow the front outwards: away
 * from any collision it advances one unit per unit, so ||grad T|| = 1.  Where
 * two banks stand back to back the nearest white matter flips from one bank to
 * the other, T turns around, and a centred difference straddling that turn
 * cancels -- ||grad T|| collapses towards 0.  The set where it does is the
 * medial surface between the banks, i.e. the sulcus, and it exists whenever the
 * two banks exist, whatever the intensity looks like between them.
 *
 * Gyral crowns cannot produce a false positive, and that is a property of the
 * construction rather than a guard bolted on: a front leaving a gyral blade has
 * nothing to collide with, because outside the crown is CSF and not more
 * cortex.  Only two facing banks make a shock.  This is the shock detection of
 * ACE (Han et al. 2001), reduced to the Euclidean case that PBT already has the
 * distance map for.
 *
 * A strict q gives a thin set: for a front turning symmetrically about the
 * midline, ||grad T|| is already back to ~1 one voxel away, so a q in this
 * range keeps essentially the midline voxels alone.
 *
 * The result is flat in q over a wide band and then falls off a cliff: on a
 * fused-bank phantom 0.5 and 0.7 place the central surface identically, while
 * 0.3 finds too few collisions to change anything at all.  The default sits in
 * the middle of the flat region rather than at its edge.
 *
 * \param dist_wm   (in)  distance to the WM boundary, in the same units the
 *                        gradient is taken in (voxel steps; PBT builds it that
 *                        way, and the isotropic-step assumption is inherited
 *                        from there)
 * \param mask      (in)  optional uint8 mask limiting the search to the
 *                        cortical band; NULL searches everywhere
 * \param medial    (out) float[nvox], 1 on the medial set and 0 elsewhere
 * \param dims      (in)  {nx, ny, nz}
 * \param q         (in)  shock threshold on ||grad T||; <= 0 selects
 *                        CAT_MEDIAL_SET_Q
 * \param t_min     (in)  ignore voxels closer than this to the WM boundary,
 *                        the "T > 1" guard of ACE, which is what stops the set
 *                        from hugging the inner boundary and carving into thin
 *                        cortex
 * \return number of voxels in the medial set, or negative on error
 */
long CAT_VolSulcalMedialSet(const float *dist_wm, const unsigned char *mask,
                            float *medial, int dims[3], double q, double t_min);

/**
 * \brief Eigen-decomposition of a symmetric 3x3 matrix, sorted by magnitude.
 *
 * Closed-form solution; no iteration, no external dependency.  Eigenvalues are
 * returned sorted so that |eval[0]| <= |eval[1]| <= |eval[2]|, which is the
 * ordering the shape measures above are defined in.
 *
 * \param a     (in)  upper triangle {axx, axy, axz, ayy, ayz, azz}
 * \param eval  (out) eigenvalues sorted by ascending absolute value
 * \param evec3 (out) unit eigenvector of eval[2], or NULL to skip
 */
void CAT_EigenSym3(const double a[6], double eval[3], double evec3[3]);

#endif /* _CAT_SHEETNESS_H_ */
