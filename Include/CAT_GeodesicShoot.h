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
 * \brief Settings of CAT_VolShoot().
 *
 * Initialise with CAT_VolShootOptionsInit(); the defaults are CAT12/CAT26's
 * (cat_defaults.m: extopts.regstr = 0.5; cat_main_registration.m: nits = 64,
 * ll1th = 0.001, ll3th = 0.002).
 */
typedef struct
{
    double regstr;   /**< 0.5 (default): CAT's optimised shooting, 2.5 to the
                          template resolution in five levels with early level
                          switching; any value in (0, 1] sets the threshold of
                          that switching (1 = most iterations); 4: SPM's
                          default schedule, every level at the template
                          resolution */
    int nits;        /**< maximum number of Gauss-Newton iterations (64) */
    double ll1th;    /**< base of the level-switch threshold (0.001) */
    double lowres;   /**< coarsest registration resolution in mm (2.5) */
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
 * \brief Fill a CAT_VolShootOpts with CAT26's defaults.
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
 * \param affine    (in)  template-to-input world mapping (CAT's res.Affine);
 *                        NULL for the identity, i.e. maps already in template
 *                        space such as CAT's rp* or T1Prep's affine maps
 * \param templates (in)  CAT_SHOOT_NTEMPLATES templates, coarse to fine; each
 *                        holds n_cls volumes of tdims, one after the other
 * \param tdims     (in)  template dimensions
 * \param tmat      (in)  template voxel-to-world matrix
 * \param tpm_vx    (in)  voxel size of CAT's tissue probability map in mm,
 *                        which scales the regularisation (1.5 for SPM's TPM)
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
                 const int tdims[3], const double tmat[16], double tpm_vx,
                 const float *lesion, const CAT_VolShootOpts *opts,
                 float *y, float *iy, float *jac, CAT_VolShootInfo *info);

#endif
