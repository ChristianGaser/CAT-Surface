/*
 * CAT_WarpDemons.h
 *
 * Demon-based spherical surface registration driver used by CAT_SurfWarpDemon.
 *
 * Implements four demon variants for matching a source cortical surface to a
 * template by deforming the source sphere:
 *   1 - Thirion's passive-force demons (Thirion, Med. Image Anal. 1998)
 *   2 - Accelerated active+passive demons (Wang et al., PMB 2005)
 *   3 - Fast inverse-consistent demons (Yang et al., PMB 2008)
 *   4 - Spherical Demons (Yeo et al., IEEE TMI 2010) - default
 *
 * Method 4 performs a diffeomorphic update: each per-vertex Gauss-Newton
 * velocity step is integrated with a scaling-and-squaring exponential map and
 * composed onto the running warp, mirroring SD_SphericalExpMap.m from the
 * reference Spherical Demons toolbox.
 */

#ifndef CAT_WARP_DEMONS_H
#define CAT_WARP_DEMONS_H

#include <bicpl.h>

#ifdef __cplusplus
extern "C" {
#endif

#define CAT_WARP_DEMONS_MAX_STEPS 3

/**
 * \brief Options controlling demon-based spherical registration.
 *
 * Populate with CAT_WarpDemonsDefaults() and override fields as needed before
 * calling CAT_WarpDemonsRegister().
 */
typedef struct {
    int    n_points;            /* spherical resample resolution (e.g. 20480) */
    int    method;              /* demon variant 1..4; 4 = Spherical Demons */
    int    n_steps;             /* number of feature stages (1..3) */
    int    curvtype[CAT_WARP_DEMONS_MAX_STEPS]; /* curvature type per stage:
                                   0 mean curv (3mm, deg), 1 gaussian,
                                   2 curvedness, 3 shape index,
                                   4 mean curv (rad), 5 sulcal-depth-like */
    int    iters;               /* maximum iterations per stage */
    int    rotate;              /* rigid rotation pre-alignment on stage 0 */
    int    smooth_velocity;     /* low-pass the velocity update (fluid prior) */
    int    smooth_displacement; /* low-pass the displacement field (elastic) */
    int    use_hessian;         /* per-vertex Gauss-Newton 2x2 Hessian update */
    int    use_line_search;     /* adaptive step backtracking on stalled CC */
    int    use_expmap;          /* diffeomorphic scaling-and-squaring (method 4) */
    double fwhm_flow;           /* FWHM for velocity-update smoothing */
    double fwhm_curv;           /* FWHM for the initial curvature smoothing */
    double fwhm_disp;           /* FWHM for displacement-field smoothing */
    double rate;                /* per-iteration multiplier for fwhm_flow */
    double alpha0;              /* demon force normalization weight */
    double max_step_deg;        /* clamp per-iteration step (deg); <=0 disables */
    double sigma_x;             /* Spherical Demons regularization weight */
    double step_factor;         /* global step-size factor */
    int    verbose;             /* print per-iteration progress */
    int    debug;               /* write intermediate debug files */
} CAT_WarpDemonsOptions;

/**
 * \brief Fill an options struct with the default 2-stage Spherical Demons setup.
 *
 * Defaults to method 4 (Spherical Demons) with diffeomorphic integration,
 * n_steps = 2 (stage 0 sulcal depth, stage 1 mean curvature), rotation
 * pre-alignment, and the fluid/elastic priors enabled.
 *
 * \param opt (out) options struct to initialize
 * \return void
 */
void CAT_WarpDemonsDefaults(CAT_WarpDemonsOptions *opt);

/**
 * \brief Register a source surface to a template by warping its sphere.
 *
 * Resamples both surfaces and their spheres to opt->n_points, optionally
 * applies a rigid rotation pre-alignment, then runs the selected demon method
 * for each feature stage, carrying the accumulated warp forward between stages.
 * The deformed source sphere at full input resolution is returned.
 *
 * \param src               (in)  source cortical surface mesh
 * \param src_sphere        (in)  spherical parameterization of src
 * \param trg               (in)  template cortical surface mesh
 * \param trg_sphere        (in)  spherical parameterization of trg
 * \param warped_src_sphere (out) deformed source sphere (caller-allocated,
 *                                overwritten with the registration result)
 * \param opt               (in)  registration options
 * \return OK on success, ERROR on invalid input or failure
 */
Status CAT_WarpDemonsRegister(polygons_struct *src, polygons_struct *src_sphere,
                              polygons_struct *trg, polygons_struct *trg_sphere,
                              polygons_struct *warped_src_sphere,
                              const CAT_WarpDemonsOptions *opt);

#ifdef __cplusplus
}
#endif

#endif
