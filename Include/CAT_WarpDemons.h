/*
 * CAT_WarpDemons.h
 *
 * Spherical Demons surface registration driver used by CAT_SurfSphericalDemon.
 *
 * Matches a source cortical surface to a template by deforming the source
 * sphere (Spherical Demons, Yeo et al., IEEE TMI 2010). Each per-vertex
 * Gauss-Newton velocity step is integrated with a scaling-and-squaring
 * exponential map and composed onto the running warp, mirroring
 * SD_SphericalExpMap.m from the reference Spherical Demons toolbox. The warp is
 * refined coarse-to-fine over a multi-resolution pyramid.
 */

#ifndef CAT_WARP_DEMONS_H
#define CAT_WARP_DEMONS_H

#include <bicpl.h>

#ifdef __cplusplus
extern "C" {
#endif

#define CAT_WARP_DEMONS_MAX_STEPS 4

/**
 * \brief Options controlling demon-based spherical registration.
 *
 * Populate with CAT_WarpDemonsDefaults() and override fields as needed before
 * calling CAT_WarpDemonsRegister().
 */
typedef struct {
    int    n_points;            /* fallback spherical resolution if level_points unset */
    int    level_points[CAT_WARP_DEMONS_MAX_STEPS]; /* resolution per pyramid level
                                   (coarse -> fine); <=0 falls back to n_points */
    int    n_steps;             /* number of multi-resolution levels (1..3) */
    int    curvtype[CAT_WARP_DEMONS_MAX_STEPS]; /* curvature type per level:
                                   0 mean curv (3mm, deg), 1 gaussian,
                                   2 curvedness, 3 shape index,
                                   4 mean curv (rad), 5 sulcal-depth-like */
    int    iters;               /* maximum iterations per level */
    int    rotate;              /* rigid rotation pre-alignment on level 0 */
    int    smooth_velocity;     /* low-pass the velocity update (fluid prior; SD default off) */
    int    smooth_displacement; /* low-pass the displacement field (elastic prior; SD default on) */
    int    use_hessian;         /* per-vertex Gauss-Newton 2x2 Hessian update */
    int    use_line_search;     /* adaptive step backtracking on stalled CC */
    int    use_expmap;          /* diffeomorphic scaling-and-squaring exp map */
    int    use_tangent;         /* per-vertex tangent-plane update (SD) instead of
                                   the global lat-lon chart; requires use_expmap */
    double fwhm_flow;           /* FWHM for velocity-update smoothing (fluid) */
    double fwhm_curv;           /* FWHM for the initial curvature smoothing */
    double fwhm_disp;           /* FWHM for displacement-field smoothing (elastic) */
    double rate;                /* per-iteration multiplier for fwhm_flow */
    double max_step_deg;        /* clamp per-iteration step (deg); <=0 disables */
    double sigma_x;             /* SD regularization weight (= max_step; SD default 2) */
    double step_factor;         /* global step-size factor */
    double *std_map;            /* optional per-vertex std of the (mean-curvature)
                                   feature on the TEMPLATE mesh, length
                                   trg->n_points. When set, the Gauss-Newton data
                                   term is locally weighted by 1/variance
                                   (atlas-style, as in SD template registration);
                                   resampled to each pyramid level. NULL = off. */
    double std_exp;             /* exponent on the precision weight: w = (1/var)^e.
                                   1 = SD's 1/variance; >1 sharpens a low-contrast
                                   std map; 0 = uniform. Only used with std_map. */
    int    verbose;             /* print per-iteration progress */
    int    debug;               /* write intermediate debug files */
} CAT_WarpDemonsOptions;

/**
 * \brief Fill an options struct with the default multi-resolution SD setup.
 *
 * Uses Spherical Demons (Yeo et al.) with diffeomorphic integration and a
 * 2-level coarse-to-fine sulcal-depth pyramid (5120 -> 20480 points), which is
 * the empirical accuracy/runtime sweet spot. Smoothing FWHM is scaled by mesh
 * spacing, anchored to a fixed reference resolution. Rotation pre-alignment and
 * the SD regularization (constant Tikhonov, sigma_x = 2, elastic displacement
 * smoothing) are enabled.
 *
 * \param opt (out) options struct to initialize
 * \return void
 */
void CAT_WarpDemonsDefaults(CAT_WarpDemonsOptions *opt);

/**
 * \brief Register a source surface to a template by warping its sphere.
 *
 * Resamples both surfaces and their spheres to opt->n_points, optionally
 * applies a rigid rotation pre-alignment, then runs Spherical Demons
 * registration for each feature stage, carrying the accumulated warp forward
 * between stages. The deformed source sphere at full input resolution is
 * returned.
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
