/*
 * CAT_CorrectThicknessFolding.h
 *
 * High-level helper to compensate cortical thickness for
 * folding-related variation.
 */

#ifndef CAT_CORRECT_THICKNESS_FOLDING_H
#define CAT_CORRECT_THICKNESS_FOLDING_H

#include <bicpl.h>

#ifdef __cplusplus
extern "C"
{
#endif

    /**
     * \brief Correct thickness values for folding-related variation (weighted).
     *
     * The correction is applied in-place and preserves the original mean. It
     * matches CAT_SurfCorrectThicknessFolding:
     * - uses 4 folding-related curvature measures (curvtype 1..4)
     * - uses linear and squared terms for each measure
     * - applies heat-kernel smoothing with FWHM = 3 mm
     * - removes the projection of thickness onto the resulting design matrix
     *
     * The correction term is scaled by a per-vertex weight using a bounded
     * transfer function, w_i = 1 + |slope| * tanh(f_i), where f_i is the
     * z-scored thickness, so larger thickness values get stronger correction.
     * It is only applied where the mean curvature, measured before it is
     * centred, is positive. With outward normals that is convex (gyral)
     * cortex -- the opposite of FreeSurfer's ?h.curv sign.
     *
     * \param polygons  (in)     surface mesh
     * \param n_vals    (in)     number of thickness values (must equal n_points)
     * \param thickness (in/out) thickness values, corrected in-place
     * \param slope     (in)     weighting strength; its sign is ignored and 0
     *                          reduces this to the unweighted correction
     * \return OK on success, ERROR otherwise
     */
    Status CAT_CorrectThicknessFoldingWeighted(polygons_struct *polygons, int n_vals,
                                               double *thickness, double slope);

    /**
     * \brief Correct thickness values for folding-related variation (unweighted).
     *
     * Same as CAT_CorrectThicknessFoldingWeighted() with slope = 0.
     *
     * \param polygons  (in)     surface mesh
     * \param n_vals    (in)     number of thickness values (must equal n_points)
     * \param thickness (in/out) thickness values, corrected in-place
     * \return OK on success, ERROR otherwise
     */
    Status CAT_CorrectThicknessFolding(polygons_struct *polygons, int n_vals,
                                       double *thickness);

#ifdef __cplusplus
}
#endif

#endif
