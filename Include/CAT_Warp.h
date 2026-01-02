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
#define TOL 1e-6     // Convergence tolerance

typedef struct {
    polygons_struct *src;
    polygons_struct *src_sphere;
    polygons_struct *trg_sphere;
    double *orig_trg;  // Precomputed target curvatures
    double *map_trg;   // Preallocated buffer for rotated target curvatures
    double *map_src;   // Precomputed source curvatures
    int distortion_correction;  // Apply sin(theta) weighting for distortion correction (0=off, 1=on)
} OptimizationParams;

void rotate_polygons(polygons_struct *, polygons_struct *, double *rotation_matrix);
void rotation_to_matrix(double *, double, double, double gamma);          
void apply_warp(polygons_struct *, polygons_struct *, double *, int *, int);
void apply_uv_warp(polygons_struct *, polygons_struct *, double *, double *, int );
void average_xz_surf(polygons_struct *, polygons_struct *, polygons_struct *);
void rotate_polygons_to_atlas(polygons_struct *, polygons_struct *,
             polygons_struct *, polygons_struct *, double, int, double *, int, int);

#endif
