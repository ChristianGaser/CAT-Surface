/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 *
 */

/*
 * Compensate cortical thickness for folding-related variation.
 *
 * This implements a simple regression-based correction inspired by:
 *
 *   Nagehan Demirci, Timothy S. Coalson, Maria A. Holland,
 *   David C. Van Essen, Matthew F. Glasser
 *   "Compensating Cortical Thickness for Cortical Folding-Related Variation"
 *   https://doi.org/10.1101/2025.05.03.651968
 *
 * For each vertex, we build a design matrix from multiple curvature-based
 * folding measures (curvtype 1..4), including linear and squared terms.
 * The folding measures are smoothed on the surface using a heatkernel
 * smoothing (FWHM = 3mm) and normalized. We then project thickness onto
 * this design matrix (via pseudo-inverse) and subtract the projection.
 * The original mean thickness is preserved.
 */

#include <bicpl.h>
#include <math.h>

#include "CAT_CorrectThicknessFolding.h"
#include "CAT_Curvature.h"
#include "CAT_Math.h"
#include "CAT_Smooth.h"

/**
 * \brief Correct thickness values for folding-related variation (weighted).
 *
 * Builds a regression design matrix from multiple curvature measures and
 * removes the projection from thickness values. Optionally applies a
 * per-vertex weighting based on deviation from the mean thickness.
 *
 * Weighting/correction with non-zero slope is restricted to vertices with
 * positive mean curvature, targeting sulcal regions where over-correction is
 * most likely.
 *
 * \param polygons  (in)  surface mesh
 * \param n_vals    (in)  number of thickness values (must match n_points)
 * \param thickness (in/out) thickness values, corrected in-place
 * \param slope     (in)  thickness-based weighting strength (0.0 disables weighting)
 * \return OK on success, ERROR otherwise
 */
Status
CAT_CorrectThicknessFoldingWeighted(polygons_struct *polygons, int n_vals,
                                    double *thickness, double slope)
{
    double *curvatures, *mean_curvatures;
    double *orig_thickness;
    double **G, **invG, *beta;
    double mean_thickness, std_thickness, std_curvature;
    int i, j;
    int *n_neighbours, **neighbours;

    if (polygons == NULL || thickness == NULL)
        return (ERROR);

    if (n_vals <= 0 || polygons->n_points != n_vals)
        return (ERROR);

    /* Keep original thickness for weighting, then remove mean for regression. */
    orig_thickness = (double *)malloc(sizeof(double) * n_vals);
    for (i = 0; i < n_vals; i++)
        orig_thickness[i] = thickness[i];

    /* remove mean from thickness (regression is done in original units) */
    mean_thickness = get_mean_double(orig_thickness, n_vals, 0);
    std_thickness = get_std_double(orig_thickness, n_vals, 0);
    if (!isfinite(std_thickness) || std_thickness <= 1e-12)
        std_thickness = 1.0;

    for (i = 0; i < n_vals; i++)
        thickness[i] -= mean_thickness;

    get_all_polygon_point_neighbours(polygons, &n_neighbours, &neighbours);

    curvatures = (double *)malloc(sizeof(double) * n_vals);
    mean_curvatures = (double *)malloc(sizeof(double) * n_vals);

    /* add gaussian curvature, curvedness, shape index, mean curvature to G */
    int curvtype[4] = {1, 2, 3, 4};
    int n_curvtypes = sizeof(curvtype) / sizeof(curvtype[0]);
    int n_beta = 1 + 2 * (n_curvtypes);

    ALLOC2D(G, n_vals, n_beta);
    ALLOC2D(invG, n_beta, n_vals);
    beta = (double *)malloc(sizeof(double) * n_beta);

    /* Create design matrix G using linear and squared terms */
    for (j = 0; j < n_curvtypes; j++)
    {
        get_polygon_vertex_curvatures_cg(polygons, n_neighbours,
                                         neighbours, 0.0, curvtype[j],
                                         curvatures);

        /* Smooth foldings with FWHM of 3mm and normalize to mean zero */
        smooth_heatkernel(polygons, curvatures, 3.0);
        normalize_double(curvatures, n_vals);

        /* Add linear term */
        for (i = 0; i < n_vals; i++)
            G[i][2 * j] = curvatures[i];

        /* Rescue mean curvature */
        if (curvtype[j] == 4)
        {
            std_curvature = get_std_double(curvatures, n_vals, 0);
            if (!isfinite(std_curvature) || std_curvature <= 1e-12)
                std_curvature = 1.0;
            for (i = 0; i < n_vals; i++)
                mean_curvatures[i] = curvatures[i];
        }

        /* Add squared term */
        for (i = 0; i < n_vals; i++)
            curvatures[i] *= curvatures[i];
        normalize_double(curvatures, n_vals);
        for (i = 0; i < n_vals; i++)
            G[i][2 * j + 1] = curvatures[i];
    }

    /* also add constant as last column */
    for (i = 0; i < n_vals; i++)
        G[i][n_beta - 1] = 1.0;

    /* Compute pseudo inverse from design matrix */
    (void)pinv(n_vals, n_beta, G, invG);

    /* Get betas */
    for (i = 0; i < n_beta; i++)
    {
        beta[i] = 0.0;
        for (j = 0; j < n_vals; j++)
            beta[i] += invG[i][j] * thickness[j];
    }

    /* Correct thickness by removing effects due to folding */
    for (i = 0; i < n_vals; i++)
    {
        double proj = 0.0;
        double w = 1.0;
        for (j = 0; j < n_beta; j++)
            proj += G[i][j] * beta[j];

        if (slope != 0.0)
        {
            double signed_feature;

            signed_feature = (orig_thickness[i] - mean_thickness) / std_thickness;

            /*
             * Bound the weighting response to avoid extreme correction strength
             * while preserving directionality (larger feature -> larger weight).
             */
            if (mean_curvatures[i] < 0)
                w = 0.0;
            else
                w = 1.0 + fabs(slope) * tanh(signed_feature);

            if (!isfinite(w) || w < 0.0)
                w = 0.0;
        }

        thickness[i] -= w * proj;
    }

    free(beta);
    FREE2D(G);
    FREE2D(invG);

    /* add mean again */
    for (i = 0; i < n_vals; i++)
        thickness[i] += mean_thickness;

    free(curvatures);
    free(mean_curvatures);
    free(orig_thickness);

    /* Cleanup (get_all_polygon_point_neighbours allocates a single flat block for neighbours[0]). */
    if (n_neighbours)
        free(n_neighbours);
    if (neighbours)
    {
        if (neighbours[0])
            free(neighbours[0]);
        free(neighbours);
    }

    return (OK);
}

/**
 * \brief Correct thickness values for folding-related variation (unweighted).
 *
 * Convenience wrapper that applies the weighted correction with slope 0.0.
 *
 * \param polygons  (in)  surface mesh
 * \param n_vals    (in)  number of thickness values (must match n_points)
 * \param thickness (in/out) thickness values, corrected in-place
 * \return OK on success, ERROR otherwise
 */
Status
CAT_CorrectThicknessFolding(polygons_struct *polygons, int n_vals,
                            double *thickness)
{
    return (CAT_CorrectThicknessFoldingWeighted(polygons, n_vals, thickness, 0.0));
}
