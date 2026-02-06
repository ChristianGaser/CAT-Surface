/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 *
 */

#ifndef _CAT_AVERAGE_SURF_H_
#define _CAT_AVERAGE_SURF_H_

/**
 * @file CAT_AverageSurf.h
 * @brief Surface averaging and RMS computation library.
 *
 * Provides functions for averaging multiple surfaces and computing
 * root-mean-square error across surface sets.
 */

#include <bicpl.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Compute the average of multiple surfaces.
 *
 * All surfaces must have the same number of points. The output surface
 * inherits connectivity from the first input surface.
 *
 * @param surfaces      Array of pointers to input surfaces.
 * @param n_surfaces    Number of surfaces in the array.
 * @param avg_out       Output: averaged surface (connectivity copied from surfaces[0]).
 *
 * @return 0 on success, non-zero on error.
 */
int CAT_SurfComputeAverage(
    polygons_struct **surfaces,
    int n_surfaces,
    polygons_struct *avg_out
);

/**
 * Compute per-vertex RMS (root-mean-square) deviation from the average surface.
 *
 * This measures the geometric variability at each vertex across the surface set.
 *
 * @param surfaces      Array of pointers to input surfaces.
 * @param n_surfaces    Number of surfaces (must be >= 2 for meaningful RMS).
 * @param avg_surface   The average surface (computed via CAT_SurfComputeAverage).
 * @param rms_out       Output: per-vertex RMS values (caller must allocate n_points doubles).
 *
 * @return 0 on success, non-zero on error.
 */
int CAT_SurfComputeRMS(
    polygons_struct **surfaces,
    int n_surfaces,
    polygons_struct *avg_surface,
    double *rms_out
);

#ifdef __cplusplus
}
#endif

#endif /* _CAT_SURF_AVERAGE_H_ */
