/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "CAT_SurfPialWhite.h"
#include "CAT_Surf.h"
#include "CAT_Vol.h"
#include "CAT_Smooth.h"
#include "CAT_Deform.h"
#include "CAT_Curvature.h"

/* Tissue class thresholds */
#ifndef CGM
#define CGM 1.5
#endif
#ifndef GWM
#define GWM 2.5
#endif

void
CAT_PialWhiteOptionsInit(CAT_PialWhiteOptions *opts)
{
    if (!opts) return;
    opts->w1 = 0.05;
    opts->w2 = 0.05;
    opts->w3 = 0.05;
    opts->sigma = 0.2;
    opts->iterations = 100;
    opts->verbose = 0;
}

int
CAT_SurfEstimatePialWhite(
    polygons_struct *central,
    const double *thickness_values,
    float *labels,
    nifti_image *nii_ptr,
    polygons_struct *pial_out,
    polygons_struct *white_out,
    const CAT_PialWhiteOptions *opts
)
{
    int p;
    int n_points;
    double *extents = NULL;
    double *weight = NULL;
    int *n_neighbours = NULL;
    int **neighbours = NULL;
    object_struct **objects_out;
    polygons_struct *polygons_pial;
    polygons_struct *polygons_white;
    polygons_struct *polygons_smoothed = NULL;
    double weights[3];
    double shifting[2] = {-0.25, 0.25};

    if (!central || !thickness_values || !labels || !nii_ptr || 
        !pial_out || !white_out || !opts)
        return -1;

    n_points = central->n_points;

    /* Allocate working arrays */
    extents = (double *)malloc(sizeof(double) * n_points);
    weight = (double *)malloc(sizeof(double) * n_points);
    polygons_smoothed = (polygons_struct *)malloc(sizeof(polygons_struct));

    if (!extents || !weight || !polygons_smoothed) {
        if (extents) free(extents);
        if (weight) free(weight);
        if (polygons_smoothed) free(polygons_smoothed);
        return -2;
    }

    /* Get neighbours for curvature computation */
    get_all_polygon_point_neighbours(central, &n_neighbours, &neighbours);

    /* Compute curvature-based weights (negative mean curvature -> smoothing) */
    get_polygon_vertex_curvatures_cg(central, n_neighbours, neighbours, 3.0, 0.0, weight);
    for (p = 0; p < n_points; p++) {
        weight[p] = fmin(0.0, weight[p]);   /* Only negative curvatures */
        weight[p] = fmax(-90.0, weight[p]); /* Clip at -90 */
        weight[p] /= -90.0;                 /* Normalize to [0..1] */
    }

    /* Initial estimate of pial surface */
    for (p = 0; p < n_points; p++) extents[p] = 0.5;
    objects_out = central_to_new_pial(central, (double *)thickness_values, extents, 
                                       1, 0.5 * opts->sigma, 5, opts->verbose);
    polygons_pial = get_polygons_ptr(objects_out[0]);

    /* Smooth pial surface based on local curvature */
    copy_polygons(polygons_pial, polygons_smoothed);
    smooth_heatkernel(polygons_smoothed, NULL, 5.0);

    /* Blend original and smoothed using curvature weights */
    for (p = 0; p < n_points; p++) {
        Point_x(polygons_pial->points[p]) = weight[p] * Point_x(polygons_smoothed->points[p]) +
                                            (1.0 - weight[p]) * Point_x(polygons_pial->points[p]);
        Point_y(polygons_pial->points[p]) = weight[p] * Point_y(polygons_smoothed->points[p]) +
                                            (1.0 - weight[p]) * Point_y(polygons_pial->points[p]);
        Point_z(polygons_pial->points[p]) = weight[p] * Point_z(polygons_smoothed->points[p]) +
                                            (1.0 - weight[p]) * Point_z(polygons_pial->points[p]);
    }

    /* Initial estimate of white surface */
    for (p = 0; p < n_points; p++) extents[p] = -0.5;
    objects_out = central_to_new_pial(central, (double *)thickness_values, extents,
                                       0, 0.5 * opts->sigma, 5, opts->verbose);
    polygons_white = get_polygons_ptr(objects_out[0]);

    /* Dual-surface deformation */
    weights[0] = opts->w1;
    weights[1] = opts->w2;
    weights[2] = opts->w3;
    surf_deform_dual(polygons_pial, polygons_white, central, labels, nii_ptr,
                     weights, opts->sigma, CGM + shifting[0], GWM + shifting[1],
                     (double *)thickness_values, opts->iterations, opts->verbose);

    /* Copy results to output */
    copy_polygons(polygons_pial, pial_out);
    copy_polygons(polygons_white, white_out);

    /* Cleanup */
    free(extents);
    free(weight);
    delete_polygon_point_neighbours(central, n_neighbours, neighbours, NULL, NULL);
    /* Note: polygons_pial/white are managed by objects_out, don't double-free */
    free(polygons_smoothed->points);
    free(polygons_smoothed->normals);
    free(polygons_smoothed->indices);
    free(polygons_smoothed->end_indices);
    free(polygons_smoothed);

    return 0;
}
