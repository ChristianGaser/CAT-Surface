/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Geodesic shooting of tissue maps onto a series of shooting templates: a C
 * port of CAT12/CAT26's optimised shooting (run_Shooting and regsetup in
 * cat_main_registration.m, Christian Gaser and Robert Dahnke), built on SPM's
 * shooting primitives (3rdparty/spm_shoot, John Ashburner).
 *
 * Copyright Christian Gaser, University of Jena.
 */

#ifndef _CAT_GEODESICSHOOT_H_
#define _CAT_GEODESICSHOOT_H_

/** Number of templates of a shooting series (Template_0 .. Template_4). */
#define CAT_SHOOT_NTEMPLATES 5

/** Largest number of tissue classes that can be registered (GM, WM, CSF). */
#define CAT_SHOOT_MAXCLASSES 3

/**
 * FWHM in mm by which the template of a level is derived from the final one
 * when only that is given: CAT's MNI152NLin2009cAsym series is close to it
 * (Template_k_GS ~ Template_4_GS smoothed by these, R^2 0.96-0.999).
 */
#define CAT_SHOOT_LEVEL_FWHM {6.5, 3.0, 1.5, 1.0, 0.0}

/**
 * \brief Settings of CAT_VolShoot().
 *
 * Initialise with CAT_VolShootOptionsInit().  The defaults are the setting
 * with the best label accuracy on Mindboggle-101 (T1Prep's evaluation):
 * regstr 1 and half of CAT's regularisation.  CAT12/CAT26 itself runs
 * regstr 0.5 and reg_scale 1 (cat_defaults.m: extopts.regstr = 0.5;
 * cat_main_registration.m: nits = 64, ll1th = 0.001).
 */
typedef struct
{
    double regstr;   /**< 1 (default): CAT's optimised shooting, 2.5 mm to the
                          template resolution in five levels, every scheduled
                          iteration run; values in (0, 1) switch to the next
                          level early (CAT's default 0.5 after about three
                          iterations per level); 4: SPM's default schedule,
                          every level at the template resolution */
    int nits;        /**< maximum number of Gauss-Newton iterations (64) */
    double ll1th;    /**< base of the level-switch threshold (0.001) */
    double lowres;   /**< coarsest registration resolution in mm (2.5) */
    double reg_scale; /**< factor on CAT's regularisation (0.5; CAT: 1);
                           smaller values allow more local deformation */
    int hydro;       /**< separate very large ventricles as their own class,
                          as CAT does (1) */
    int hotspot;     /**< filter local Jacobian hot spots of the velocity,
                          as CAT does (1) */
    int n_threads;   /**< worker threads; <= 0 uses the library default */
    int verbose;     /**< 1 prints CAT's iteration log to stdout */
} CAT_VolShootOpts;

/**
 * \brief What CAT_VolShoot() did, for reports and tests.
 */
typedef struct
{
    int iterations;                          /**< Gauss-Newton iterations run */
    int level_iterations[CAT_SHOOT_NTEMPLATES]; /**< iterations per level */
    double level_res[CAT_SHOOT_NTEMPLATES];  /**< registration resolution (mm) */
    double ll_match;                         /**< final matching term per voxel */
    double ll_reg;                           /**< final regularisation term */
    int hydro_levels;                        /**< levels with the ventricle class */
    int hotspot_corrections;                 /**< iterations that filtered u */
    double det_min;                          /**< smallest Jacobian determinant */
    double det_max;                          /**< largest Jacobian determinant */
} CAT_VolShootInfo;

/**
 * \brief An intensity channel matched by sum of squares (experimental).
 *
 * A label or T1 image matched to a template intensity from level \p from on,
 * next to the multinomial matching of the tissue classes or instead of it.
 * The template intensity is mapped to the subject's by a polynomial refitted
 * at every iteration.  With tissue classes, \p weight 1 gives the term the
 * summed curvature of the multinomial term at the first iteration of the
 * level; without them, \p weight multiplies the maximum-likelihood precision
 * 1/sigma^2 of the residuals at that iteration.  Voxels where the subject's
 * image is not positive (outside a skull-stripped brain) do not count, nor
 * do those of a lesion mask passed to CAT_VolShootEx().
 */
typedef struct
{
    const float *img;                       /**< subject intensity, input grid */
    const float *tpl[CAT_SHOOT_NTEMPLATES]; /**< template intensity per level on
                                                 the template grid; tpl[4] is
                                                 required, a NULL level is
                                                 tpl[4] smoothed by
                                                 CAT_SHOOT_LEVEL_FWHM */
    const float *mask;   /**< weight on the template grid (NULL: 1 everywhere) */
    double weight;       /**< see above */
    int degree;          /**< degree of the intensity mapping (1..5) */
    int from;            /**< first level (0..4) with the term */
} CAT_VolShootIntensity;

/**
 * \brief CAT_VolShoot() with an optional intensity channel (experimental).
 *
 * As CAT_VolShoot(), with two extensions for testing: \p n_cls may be 0
 * (\p cls and \p templates NULL) when \p ich is given, which registers by
 * the intensity alone; and a NULL entry of \p templates below the last is
 * the final template smoothed by CAT_SHOOT_LEVEL_FWHM.
 *
 * \param ich (in) intensity channel; NULL runs CAT_VolShoot() unchanged
 * \return as CAT_VolShoot()
 */
int CAT_VolShootEx(const float *const cls[], int n_cls, const int dims[3],
                   const double mat[16], const double affine[16],
                   const float *const templates[CAT_SHOOT_NTEMPLATES],
                   const int tdims[3], const double tmat[16], double reg_vx,
                   const float *lesion, const CAT_VolShootIntensity *ich,
                   const CAT_VolShootOpts *opts, float *y, float *iy,
                   float *jac, CAT_VolShootInfo *info);

/**
 * \brief Label template of a shooting template: 1 + GM + 2 WM.
 *
 * The partial-volume label of CAT and T1Prep (CSF 1, GM 2, WM 3) that the
 * tissue probabilities of a shooting template predict, as the template of a
 * label image registered with CAT_VolShootEx().  Where the template holds
 * neither GM nor WM -- CSF and the space around the brain alike -- it is 1.
 * No brain mask is needed: the subject's background (label 0) takes no part
 * in the matching, so the edge of its skull-strip exerts no force.  A
 * template that is 0 outside a brain mask, with the background matched,
 * was tested on Mindboggle-101 and adds nothing at a matched weight.
 *
 * \param tpl   (in)  template with GM and WM as its first two volumes
 * \param tdims (in)  template dimensions
 * \param label (out) tdims: the label template
 */
void CAT_VolShootLabelTemplate(const float *tpl, const int tdims[3], float *label);

/**
 * \brief Fill a CAT_VolShootOpts with the defaults (best label accuracy).
 *
 * \param opts (out) options to initialise
 */
void CAT_VolShootOptionsInit(CAT_VolShootOpts *opts);

/**
 * \brief Geodesic shooting of tissue maps onto a shooting template series.
 *
 * Reproduces CAT26's run_Shooting: the subject's tissue maps are sampled onto
 * a registration grid that is refined from opts->lowres to the template
 * resolution together with the template series, and at each level the
 * initial velocity of a geodesic is updated by Gauss-Newton iterations of the
 * multinomial matching of SPM's shooting toolbox.  CAT's safeguards -- the
 * ventricle class for very large ventricles, the filtering of Jacobian hot
 * spots and the masking of lesions -- are included.
 *
 * Voxel-to-world matrices are 4x4, row-major, and follow the NIfTI
 * convention: they map 0-based voxel indices to millimetres.  Volumes are
 * stored x fastest.
 *
 * \param cls       (in)  n_cls tissue maps on the input grid (GM, WM[, CSF]),
 *                        probabilities in [0, 1]
 * \param n_cls     (in)  number of tissue classes, 2 or 3; the templates must
 *                        hold at least as many volumes
 * \param dims      (in)  input grid dimensions
 * \param mat       (in)  input grid voxel-to-world matrix
 * \param affine    (in)  input-to-template world mapping (CAT's res.Affine);
 *                        NULL for the identity, i.e. maps already in template
 *                        space such as CAT's rp* or T1Prep's affine maps
 * \param templates (in)  CAT_SHOOT_NTEMPLATES templates, coarse to fine; each
 *                        holds n_cls volumes of tdims, one after the other
 * \param tdims     (in)  template dimensions
 * \param tmat      (in)  template voxel-to-world matrix
 * \param reg_vx    (in)  voxel size of the regulariser in mm; <= 0 uses the
 *                        template's.  CAT takes its TPM's 1.5 mm at every
 *                        template resolution, which equals the template's for
 *                        its 1.5 mm series but weakens the regularisation of a
 *                        1 mm series
 * \param lesion    (in)  optional lesion mask on the input grid (NULL: none)
 * \param opts      (in)  settings; NULL uses CAT_VolShootOptionsInit()
 * \param y         (out) tdims x 3: for every template voxel the input voxel
 *                        (0-based) it pulls from
 * \param iy        (out) optional (NULL), dims x 3: for every input voxel the
 *                        template voxel (0-based) it maps to
 * \param jac       (out) optional (NULL), tdims: Jacobian determinant of y in
 *                        voxel units
 * \param info      (out) optional (NULL) summary of the run
 * \return 0 on success, 1 on invalid input, 2 on memory failure
 */
int CAT_VolShoot(const float *const cls[], int n_cls, const int dims[3],
                 const double mat[16], const double affine[16],
                 const float *const templates[CAT_SHOOT_NTEMPLATES],
                 const int tdims[3], const double tmat[16], double reg_vx,
                 const float *lesion, const CAT_VolShootOpts *opts,
                 float *y, float *iy, float *jac, CAT_VolShootInfo *info);

#endif
