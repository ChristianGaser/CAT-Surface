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
     * Correct thickness values for folding-related variation.
     *
     * The correction is applied in-place and preserves the original mean.
     *
     * Parameters are chosen to match the behavior of the command line tool
     * CAT_SurfCorrectThicknessFolding:
     * - Uses 4 folding-related curvature measures (curvtype 1..4)
     * - Uses linear and squared terms for each measure
     * - Applies heatkernel smoothing with FWHM = 3mm
     * - Removes the projection of thickness onto the resulting design matrix
     *
     * @param polygons Surface mesh.
     * @param n_vals Number of thickness values (must equal polygons->n_points).
     * @param thickness Thickness values (in/out).
     *
     * @return OK on success, ERROR otherwise.
     */
    /**
     * Correct thickness values for folding-related variation with optional
     * slope-controlled weighting.
     *
     * The correction term is scaled by a per-vertex weight using a bounded
     * transfer function:
     *
     *   w_i = 1 + |slope| * tanh(f_i)
     *
    * where f_i is z-scored thickness. Larger thickness values therefore get
    * stronger correction and smaller thickness values get weaker correction.
    *
    * With non-zero slope, weighting/correction is only applied for vertices
    * with positive mean curvature (sulcal-focused correction).
     *
     * With slope=0, this reduces to the unweighted correction.
     *
    * @param slope Thickness-based weighting strength; sign is ignored and
    *              magnitude controls strength.
     */
    Status CAT_CorrectThicknessFoldingWeighted(polygons_struct *polygons, int n_vals,
                                               double *thickness, double slope);

    /**
     * Backward-compatible unweighted correction (slope = 0).
     */
    Status CAT_CorrectThicknessFolding(polygons_struct *polygons, int n_vals,
                                       double *thickness);

#ifdef __cplusplus
}
#endif

#endif
