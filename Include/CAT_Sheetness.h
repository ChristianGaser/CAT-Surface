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
    int polarity;     /**< +1 bright sheets, -1 dark sheets, 0 either */
    int verbose;      /**< 1 to report progress */
} CAT_SheetnessOpts;

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

/**
 * \brief Median filter over a sheetness-oriented neighbourhood.
 *
 * Drop-in replacement for an isotropic 3x3x3 median that cannot close a thin
 * sheet.  A neighbour at offset d is admitted only when
 *
 *     1 - s * (dhat . n)^2 > 0.5
 *
 * with s the local sheetness and n the local sheet normal.  For s = 0 every
 * neighbour is admitted and the operator is bit-identical to the plain
 * isotropic median, so behaviour away from thin structures is unchanged.  For
 * s = 1 only offsets within 45 degrees of the sheet plane survive, i.e. the
 * filter averages *along* the sheet and never across it.
 *
 * \param vol       (in/out) volume to filter, in place
 * \param sheetness (in)     float[nvox] in [0,1]; NULL falls back to isotropic
 * \param normal    (in)     float[3*nvox] unit sheet normals; NULL falls back
 * \param mask      (in)     optional uint8 mask; NULL = filter everywhere
 * \param dims      (in)     {nx, ny, nz}
 * \param iters     (in)     number of successive passes
 * \return 0 on success, -1 on invalid arguments, -2 on allocation failure
 */
int CAT_VolOrientedMedian(float *vol, const float *sheetness, const float *normal,
                          const unsigned char *mask, int dims[3], int iters);

/**
 * \brief Sheetness-oriented (coherence-enhancing) smoothing.
 *
 * Diffuses along the sheet and not across it: the 3x3x3 neighbourhood is
 * weighted by (1 - s) + s * exp(-(dhat . n)^2 / (2 sigma^2)), so tangential
 * neighbours keep full weight while normal neighbours are suppressed in
 * proportion to the local sheetness.  With s = 0 the weights are the plain
 * distance weights and the operator reduces to ordinary local averaging.
 *
 * \param vol       (in/out) volume to filter, in place
 * \param sheetness (in)     float[nvox] in [0,1]; NULL falls back to isotropic
 * \param normal    (in)     float[3*nvox] unit sheet normals; NULL falls back
 * \param mask      (in)     optional uint8 mask; NULL = filter everywhere
 * \param dims      (in)     {nx, ny, nz}
 * \param sigma     (in)     angular width of the anisotropy, e.g. 0.5
 * \param iters     (in)     number of successive passes
 * \return 0 on success, -1 on invalid arguments, -2 on allocation failure
 */
int CAT_VolOrientedSmooth(float *vol, const float *sheetness, const float *normal,
                          const unsigned char *mask, int dims[3], double sigma,
                          int iters);

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
