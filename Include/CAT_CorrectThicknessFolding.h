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
extern "C" {
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
 * thickness-dependent weighting.
 *
 * The correction term is scaled by a per-vertex weight:
 *
 *   w_i = max(0, 1 + slope * (t_i - mean(t)))
 *
 * where t_i are the original (uncorrected) thickness values.
 * With slope=0, this reduces to the unweighted correction.
 *
 * @param slope Steepness of weighting (units: 1 / thickness-unit).
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
