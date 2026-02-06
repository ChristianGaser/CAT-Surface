/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 *
 */

#ifndef _CAT_VOL2SURFUTILS_H_
#define _CAT_VOL2SURFUTILS_H_

/**
 * Helper utilities shared by the CAT_Vol2Surf program.
 *
 * These functions implement the mapping function evaluation and the kernel
 * construction used for exponential and gaussian-weighted averaging.
 */

#include "CAT_Math.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Evaluate the mapping function for an array of sampled values.
 *
 * \param val_array Array of sampled values.
 * \param n_val Number of values in val_array.
 * \param map_func Mapping function selector (e.g. F_MEAN, F_MAX, ...).
 * \param kernel Optional weights for weighted functions (F_WAVERAGE, F_EXP).
 * \param index_out If non-NULL, receives the index of the selected value for
 *                  functions that pick an element (e.g. max/min).
 *
 * \return The evaluated value.
 */
double CAT_Vol2SurfEvaluateFunction(const double *val_array, int n_val,
                                   int map_func, const double *kernel,
                                   int *index_out);

/**
 * Build an exponential kernel with half-decay distance exp_half.
 * The kernel is normalized to sum to 1.
 */
void CAT_Vol2SurfBuildExpKernel(const double *length_array, int n,
                               double exp_half, double *kernel_out);

/**
 * Build the gaussian kernel used by the "weighted average" mapping.
 * The kernel is normalized to sum to 1.
 */
void CAT_Vol2SurfBuildGaussianKernel50(int grid_steps, int grid_steps1,
                                      double *kernel_out);

#ifdef __cplusplus
}
#endif

#endif
