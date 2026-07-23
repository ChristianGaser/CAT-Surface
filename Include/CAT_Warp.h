/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 *
 */

#ifndef _CAT_WARP_H_
#define _CAT_WARP_H_

#include <bicpl.h>
#include <float.h>

#define INVERSE_WARPING 0
#define RADIANS(deg) ((PI * (double)(deg)) / 180.0)
#define DEGREES(rad) ((180.0 * (double)(rad)) / PI)

#define MAX_ITER 500
#define ALPHA 1.0    // Reflection coefficient
#define GAMMA 2.0    // Expansion coefficient
#define RHO 0.5      // Contraction coefficient
#define SIGMA 0.5    // Shrinkage coefficient
#define TOL 1e-4     // Convergence tolerance

typedef struct {
    polygons_struct *src;
    polygons_struct *src_sphere;
    polygons_struct *trg_sphere;
    double *orig_trg;  // Precomputed target curvatures
    double *map_trg;   // Preallocated buffer for rotated target curvatures
    double *map_src;   // Precomputed source curvatures
    double *pre_rot;   // Optional 3x3 row-major seed rotation pre-multiplied
                       // into every candidate (NULL = none); lets the
                       // Nelder-Mead refine reside in the residual angles.
} OptimizationParams;

/**
 * \brief Public API for rotate_polygons.
 *
 * This function is part of the CAT-Surface public library interface and is used by command-line tools.
 *
 * \param param (in/out) Parameter of rotate_polygons.
 * \param param (in/out) Parameter of rotate_polygons.
 * \param rotation_matrix (in/out) Parameter of rotate_polygons.
 * \return void (no return value).
 */
void rotate_polygons(polygons_struct *, polygons_struct *, double *rotation_matrix);
/**
 * \brief Public API for rotation_to_matrix.
 *
 * This function is part of the CAT-Surface public library interface and is used by command-line tools.
 *
 * \param param (in/out) Parameter of rotation_to_matrix.
 * \param double (in/out) Parameter of rotation_to_matrix.
 * \param double (in/out) Parameter of rotation_to_matrix.
 * \param gamma (in/out) Parameter of rotation_to_matrix.
 * \return void (no return value).
 */
void rotation_to_matrix(double *, double, double, double gamma);          
void apply_warp(polygons_struct *, polygons_struct *, double *, int *, int);
void apply_uv_warp(polygons_struct *, polygons_struct *, double *, double *, int );
/**
 * \brief Public API for average_xz_surf.
 *
 * This function is part of the CAT-Surface public library interface and is used by command-line tools.
 *
 * \param param (in/out) Parameter of average_xz_surf.
 * \param param (in/out) Parameter of average_xz_surf.
 * \param param (in/out) Parameter of average_xz_surf.
 * \return void (no return value).
 */
void average_xz_surf(polygons_struct *, polygons_struct *, polygons_struct *);
void rotate_polygons_to_atlas(polygons_struct *, polygons_struct *,
             polygons_struct *, polygons_struct *, double, int, double *, int);

/**
 * \brief Exhaustive coarse-to-fine global search for the initial rigid rotation.
 *
 * Mirrors FreeSurfer's MRISrigidBodyAlignGlobal. Instead of trusting a local
 * optimiser (Nelder-Mead) to descend from a seed, it evaluates the alignment cost
 * on a dense grid spanning +/- \p max_degrees over all three rotation angles,
 * re-centres on the best candidate, and halves the span only after a pass finds
 * no improvement - repeating until the span falls below \p min_degrees.
 *
 * Because every candidate within the span is evaluated, the search cannot be
 * trapped by a neighbouring-sulcus local minimum, and its capture range is the
 * full \p max_degrees rather than the width of a single basin. This is the
 * property a simplex search lacks: with a quasi-periodic folding cost, a local
 * method converges to whichever basin its seed lands in, no matter how the seed
 * is chosen.
 *
 * The returned rotation maps SOURCE coordinates onto TARGET coordinates (that is
 * the direction the underlying cost establishes). A caller that deforms the
 * source, such as CAT_WarpDemonsRegister, needs the opposite direction and must
 * apply the inverse (for an orthonormal rotation, the transpose).
 *
 * \param src             (in)  source central surface mesh
 * \param src_sphere      (in)  spherical parameterization of \p src
 * \param trg             (in)  template central surface mesh
 * \param trg_sphere      (in)  spherical parameterization of \p trg
 * \param fwhm            (in)  smoothing FWHM for the curvature cost feature
 * \param curvtype        (in)  curvature type for the cost feature
 * \param max_degrees     (in)  half-width of the initial angular search span
 * \param min_degrees     (in)  stop once the span falls below this
 * \param nangles         (in)  grid samples per axis per pass
 * \param refine          (in)  1 = Nelder-Mead refine of the residual afterwards
 * \param rotation_matrix (out) 9-element row-major rotation for the source sphere
 * \param verbose         (in)  1 for progress output; 0 silent
 * \return void
 */
void rotate_polygons_to_atlas_global(polygons_struct *src,
             polygons_struct *src_sphere, polygons_struct *trg,
             polygons_struct *trg_sphere, double fwhm, int curvtype,
             double max_degrees, double min_degrees, int nangles, int refine,
             double *rotation_matrix, int verbose);

#endif
