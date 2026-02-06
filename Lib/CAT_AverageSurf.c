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

#include "CAT_AverageSurf.h"

int
CAT_SurfComputeAverage(
    polygons_struct **surfaces,
    int n_surfaces,
    polygons_struct *avg_out
)
{
    int i, p;
    int n_points;

    if (!surfaces || n_surfaces < 1 || !avg_out)
        return -1;

    n_points = surfaces[0]->n_points;

    /* Copy first surface to output (preserves connectivity) */
    copy_polygons(surfaces[0], avg_out);

    /* Add remaining surfaces */
    for (i = 1; i < n_surfaces; i++) {
        if (surfaces[i]->n_points != n_points) {
            fprintf(stderr, "CAT_SurfComputeAverage: point count mismatch\n");
            return -2;
        }
        for (p = 0; p < n_points; p++) {
            Point_x(avg_out->points[p]) += Point_x(surfaces[i]->points[p]);
            Point_y(avg_out->points[p]) += Point_y(surfaces[i]->points[p]);
            Point_z(avg_out->points[p]) += Point_z(surfaces[i]->points[p]);
        }
    }

    /* Compute average */
    for (p = 0; p < n_points; p++) {
        Point_x(avg_out->points[p]) /= (Real)n_surfaces;
        Point_y(avg_out->points[p]) /= (Real)n_surfaces;
        Point_z(avg_out->points[p]) /= (Real)n_surfaces;
    }

    /* Recompute normals */
    compute_polygon_normals(avg_out);

    return 0;
}

int
CAT_SurfComputeRMS(
    polygons_struct **surfaces,
    int n_surfaces,
    polygons_struct *avg_surface,
    double *rms_out
)
{
    int i, p;
    int n_points;
    double dx, dy, dz;
    double sum_sq;

    if (!surfaces || n_surfaces < 2 || !avg_surface || !rms_out)
        return -1;

    n_points = avg_surface->n_points;

    /* Initialize RMS to zero */
    for (p = 0; p < n_points; p++)
        rms_out[p] = 0.0;

    /* Accumulate squared differences */
    for (i = 0; i < n_surfaces; i++) {
        if (surfaces[i]->n_points != n_points) {
            fprintf(stderr, "CAT_SurfComputeRMS: point count mismatch\n");
            return -2;
        }
        for (p = 0; p < n_points; p++) {
            dx = Point_x(surfaces[i]->points[p]) - Point_x(avg_surface->points[p]);
            dy = Point_y(surfaces[i]->points[p]) - Point_y(avg_surface->points[p]);
            dz = Point_z(surfaces[i]->points[p]) - Point_z(avg_surface->points[p]);
            rms_out[p] += dx*dx + dy*dy + dz*dz;
        }
    }

    /* Compute RMS (using n-1 for sample variance) */
    for (p = 0; p < n_points; p++) {
        rms_out[p] = sqrt(rms_out[p] / (double)(n_surfaces - 1));
    }

    return 0;
}
