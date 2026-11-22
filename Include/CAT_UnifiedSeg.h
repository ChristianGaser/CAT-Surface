/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Unified segmentation with geodesic shooting: the generative model of
 * SPM's unified segmentation (Ashburner & Friston 2005; spm_preproc8.m and
 * CAT's cat_spm_preproc8.m / cat_spm_preproc_write8.m) -- a mixture of
 * Gaussians per tissue class, a smooth multiplicative bias field and tissue
 * priors that are deformed onto the image -- with the priors taken from a
 * series of shooting templates and deformed by geodesic shooting on a grid
 * of increasing resolution (CAT_GeodesicShoot.h) instead of SPM's
 * small-deformation model on a tissue probability map.
 *
 * Copyright Christian Gaser, University of Jena.
 */

#ifndef _CAT_UNIFIEDSEG_H_
#define _CAT_UNIFIEDSEG_H_

#include "CAT_GeodesicShoot.h"

/** Tissue classes of the templates: GM, WM, CSF, bone, soft tissue, air. */
#define CAT_USEG_NCLASSES 6

/** Largest number of Gaussians of one class. */
#define CAT_USEG_MAXGAUSS 8

/** Largest number of Gaussians in total. */
#define CAT_USEG_MAXK (CAT_USEG_NCLASSES * CAT_USEG_MAXGAUSS)

/** Regularisation of the affine registration (spm_affine_priors.m). */
enum
{
    CAT_USEG_AFFREG_NONE = 0, /**< no regularisation */
    CAT_USEG_AFFREG_MNI,      /**< European brains to MNI space (default) */
    CAT_USEG_AFFREG_EASTERN,  /**< East Asian brains to MNI space */
    CAT_USEG_AFFREG_SUBJ,     /**< inter-subject registration */
    CAT_USEG_AFFREG_RIGID     /**< almost rigid */
};

/** Basis functions of the bias field. */
enum
{
    CAT_USEG_BIAS_DCT = 0, /**< SPM's discrete cosine transform */
    CAT_USEG_BIAS_SPLINE   /**< cubic B-splines (default) */
};

/** What to do with the affine matrix handed to CAT_UnifiedSeg(). */
enum
{
    CAT_USEG_AFFINE_ESTIMATE = 0, /**< estimate it; a given matrix is one more start */
    CAT_USEG_AFFINE_REFINE,       /**< refine the given matrix, no other start */
    CAT_USEG_AFFINE_FIXED         /**< use the given matrix as it is */
};

/**
 * \brief Settings of CAT_UnifiedSeg().
 *
 * Initialise with CAT_UnifiedSegOptionsInit(), the source of truth for every
 * front-end.
 */
typedef struct
{
    int ngaus[CAT_USEG_NCLASSES]; /**< Gaussians per class (1 1 2 3 4 2, as SPM
                                       and CAT) */
    double samp;        /**< distance of the voxels the parameters are estimated
                             from, in mm (3) */
    int bias_model;     /**< CAT_USEG_BIAS_* (SPLINE) */
    double biasfwhm;    /**< smoothness of the bias field as a FWHM in mm (60);
                             a spline has knots 0.75 times this apart; <= 0:
                             no bias correction */
    double biasreg;     /**< regularisation of the bias field (1e-3) */
    double tol;         /**< stopping criterion of the EM steps, as a
                             log-likelihood change per voxel (1e-4) */
    double wp_reg;      /**< holds the class weights near 1/6: pseudo-counts
                             of this multiple of the prior mass of a class
                             (0.5).  0 estimates them freely as SPM does,
                             a large value fixes them */
    int affreg;         /**< CAT_USEG_AFFREG_* (MNI) */
    int affine_mode;    /**< CAT_USEG_AFFINE_* (ESTIMATE) */
    double noise;       /**< > 0: Gaussian noise is added to the image the
                             parameters are estimated from until its noise SD
                             reaches this fraction of the intensity range (1st
                             to 99th percentile); an image that is already
                             noisier is left alone.  The final segmentation is
                             always computed from the image without it.  0: off
                             (default) */
    unsigned int noise_seed; /**< seed of the added noise (1) */
    int sanlm;          /**< 1: denoise the image with the SANLM filter first;
                             everything that follows uses the filtered image
                             (0) */
    double sanlm_strength; /**< strength of the SANLM filter (1) */
    int sanlm_rician;   /**< 1: Rician noise model of the SANLM filter (0) */
    double mrf;         /**< strength of the Markov random field cleanup of the
                             class maps, as SPM's (0: off) */
    int pve;            /**< 1: estimate partial volumes of CSF, GM and WM for
                             the label map (default); 0: the label is the
                             posterior mean CSF + 2 GM + 3 WM */
    CAT_VolShootOpts shoot; /**< the shooting: regstr, reg_scale, lowres, ...
                                 (CAT_VolShootOptionsInit(), with regstr 0.5) */
    int n_threads;      /**< worker threads; <= 0: the library default */
    int verbose;        /**< 1: stages and parameters, 2: every iteration */
} CAT_UnifiedSegOpts;

/**
 * \brief Estimated parameters and a summary of the run.
 *
 * bias_coef is allocated by CAT_UnifiedSeg(); release it with
 * CAT_UnifiedSegResultFree().
 */
typedef struct
{
    double affine[16];  /**< input-to-template world mapping (SPM's Affine),
                             4x4 row-major */
    double affine_ll;   /**< objective of the affine registration: mutual
                             information of intensity and class in bits */
    int n_gauss;        /**< number of Gaussians */
    int lkp[CAT_USEG_MAXK];   /**< class (0-based) of every Gaussian */
    double mg[CAT_USEG_MAXK]; /**< mixing proportion within its class */
    double mn[CAT_USEG_MAXK]; /**< mean, of the bias-corrected intensity */
    double vr[CAT_USEG_MAXK]; /**< variance */
    double wp[CAT_USEG_NCLASSES]; /**< class weights on the priors */
    int bias_model;     /**< CAT_USEG_BIAS_* of the bias field */
    int bias_dims[3];   /**< number of basis functions per axis */
    double *bias_coef;  /**< their coefficients, x fastest; the bias field is
                             exp() of their sum (NULL without bias correction).
                             A spline has bias_dims - 3 knot spans over the
                             voxels 0 .. dims - 1 of an axis. */
    double ll;          /**< final log-likelihood per voxel used */
    int n_samples;      /**< voxels the parameters were estimated from */
    double intensity_range[2]; /**< 1st and 99th percentile of the image */
    double noise_estimated; /**< noise SD of the image, as a fraction of the
                                 intensity range */
    double noise_added;     /**< SD of the noise added for the estimation,
                                 same unit (0: none) */
    int pve_done;       /**< partial volumes were estimated */
    double pve_mean[3]; /**< intensity of pure CSF, GM and WM */
    double pve_sd;      /**< SD of a pure tissue in the partial volume model */
    double volume[CAT_USEG_NCLASSES]; /**< volume of every class in ml; voxels
                                           without data count with their prior */
    CAT_VolShootInfo shoot;  /**< summary of the shooting */
    double time_total;  /**< seconds, wall clock */
    double time_filter; /**< SANLM filter */
    double time_affine; /**< affine registration */
    double time_fit;    /**< estimation and shooting */
    double time_write;  /**< final segmentation */
} CAT_UnifiedSegResult;

/**
 * \brief Fill a CAT_UnifiedSegOpts with the defaults.
 *
 * \param opts (out) options to initialise
 */
void CAT_UnifiedSegOptionsInit(CAT_UnifiedSegOpts *opts);

/**
 * \brief Release what CAT_UnifiedSeg() allocated inside a result.
 *
 * \param res (in/out) result; NULL is ignored
 */
void CAT_UnifiedSegResultFree(CAT_UnifiedSegResult *res);

/**
 * \brief Unified segmentation of one image with geodesic shooting.
 *
 * Classifies every voxel into the six classes of a shooting template series
 * (GM, WM, CSF, bone, soft tissue, air) with the model of SPM's unified
 * segmentation: every class is a mixture of Gaussians of the bias-corrected
 * intensity, the bias field is the exponential of a few DCT basis functions,
 * and the prior of a class is its template, deformed onto the image and
 * scaled by a class weight.  The steps:
 *
 * 1. Optionally the SANLM filter, and optionally added noise for the image
 *    the parameters are estimated from (opts->noise).
 * 2. Affine registration to the first template by maximising the mutual
 *    information of intensity and class (spm_maff8.m), from several starting
 *    points.
 * 3. On voxels opts->samp apart, the Gaussians, the class weights and the
 *    bias field are estimated by EM and Gauss-Newton as in cat_spm_preproc8.m.
 *    They alternate with the deformation: geodesic shooting on a grid refined
 *    from opts->shoot.lowres to the template resolution, through the template
 *    series, where every Gauss-Newton step is an EM step that matches the
 *    templates to the current class posteriors (CAT_VolShootFit()).
 * 4. The class posteriors of every voxel from the final parameters
 *    (cat_spm_preproc_write8.m), and a label map with partial volumes.
 *
 * Voxels that are zero or not finite carry no data: they take no part in
 * the estimation and get their prior as posterior, so skull-stripped and
 * defaced images can be segmented.
 *
 * Voxel-to-world matrices are 4x4, row-major, NIfTI convention (0-based
 * voxels to mm).  Volumes are stored x fastest.
 *
 * \param src       (in)  image on the input grid
 * \param dims      (in)  input grid dimensions
 * \param mat       (in)  input grid voxel-to-world matrix
 * \param templates (in)  CAT_SHOOT_NTEMPLATES templates, coarse to fine, each
 *                        with CAT_USEG_NCLASSES volumes of tdims that sum to
 *                        one
 * \param tdims     (in)  template dimensions
 * \param tmat      (in)  template voxel-to-world matrix
 * \param mask      (in)  optional (NULL) mask on the input grid: voxels where
 *                        it is zero carry no data (lesions)
 * \param affine    (in)  optional (NULL) input-to-template world mapping, used
 *                        as opts->affine_mode says
 * \param opts      (in)  settings; NULL uses CAT_UnifiedSegOptionsInit()
 * \param prob      (out) optional (NULL), CAT_USEG_NCLASSES volumes of dims,
 *                        one after the other: class posteriors times 255
 * \param label     (out) optional (NULL), dims: label map, 0 background,
 *                        1 CSF, 2 GM, 3 WM, with partial volumes in between
 * \param pve       (out) optional (NULL), 3 volumes of dims: fractions of
 *                        GM, WM and CSF times 255 (the partial volume estimate;
 *                        without opts->pve the posteriors of these classes)
 * \param corrected (out) optional (NULL), dims: bias-corrected image
 * \param bias      (out) optional (NULL), dims: bias field (corrected / image)
 * \param y         (out) optional (NULL), tdims x 3: for every template voxel
 *                        the input voxel (0-based) it pulls from
 * \param iy        (out) optional (NULL), dims x 3: for every input voxel the
 *                        template voxel (0-based) it maps to
 * \param jac       (out) optional (NULL), tdims: Jacobian determinant of y in
 *                        voxel units
 * \param res       (out) optional (NULL) parameters and summary; release with
 *                        CAT_UnifiedSegResultFree()
 * \return 0 on success, 1 on invalid input, 2 on memory failure, 3 if the
 *         estimation failed (no usable voxels, a degenerate model)
 */
int CAT_UnifiedSeg(const float *src, const int dims[3], const double mat[16],
                   const float *const templates[CAT_SHOOT_NTEMPLATES],
                   const int tdims[3], const double tmat[16],
                   const unsigned char *mask, const double affine[16],
                   const CAT_UnifiedSegOpts *opts, unsigned char *prob,
                   float *label, unsigned char *pve, float *corrected,
                   float *bias, float *y, float *iy, float *jac,
                   CAT_UnifiedSegResult *res);

/**
 * \brief Write the parameters of a segmentation as JSON.
 *
 * The file holds what SPM stores in its *_seg8.mat -- the affine matrix, the
 * Gaussians, the class weights and the bias field coefficients -- together
 * with the settings and the summary of the run.
 *
 * \param file  (in) output file
 * \param res   (in) result of CAT_UnifiedSeg()
 * \param opts  (in) the settings it ran with; NULL: the defaults
 * \param image (in) optional (NULL) name of the segmented image
 * \param tpl   (in) optional (NULL) name of the first template
 * \return 0 on success, 1 if the file cannot be written
 */
int CAT_UnifiedSegWriteJson(const char *file, const CAT_UnifiedSegResult *res,
                            const CAT_UnifiedSegOpts *opts, const char *image,
                            const char *tpl);

#endif
