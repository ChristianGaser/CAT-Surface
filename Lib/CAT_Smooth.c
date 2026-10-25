/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 *
 */

#include <bicpl.h>
#include "CAT_Smooth.h"

/**
 * \brief Extract all unique vertex neighbors for every point on a polygon mesh.
 *
 * Builds a complete neighbor adjacency structure by iterating through all polygons
 * and collecting the direct neighbors (vertices connected by edges) of each point.
 * Allocates compact neighbor arrays and returns them via output pointers.
 *
 * \param polygons                (in)  polygon mesh structure
 * \param n_point_neighbours_ptr (out)  int[n_points]; count of neighbors for each point
 * \param point_neighbours_ptr   (out)  int*[n_points]; neighbor vertex indices for each point
 */
void get_all_polygon_point_neighbours(polygons_struct *polygons,
                                      int *n_point_neighbours_ptr[],
                                      int **point_neighbours_ptr[])
{
    int p, poly, size, vertex;
    int *n_neighbours, **neighbours;
    int *total_neighbours, total_n_neighbours;
    BOOLEAN interior;

    check_polygons_neighbours_computed(polygons);

    ALLOC(n_neighbours, polygons->n_points);
    ALLOC(neighbours, polygons->n_points);
    for (p = 0; p < polygons->n_points; p++)
        n_neighbours[p] = 0;

    total_n_neighbours = 0;
    for (poly = 0; poly < polygons->n_items; poly++)
    {
        size = GET_OBJECT_SIZE(*polygons, poly);
        for (vertex = 0; vertex < size; vertex++)
        {
            p = polygons->indices[POINT_INDEX(polygons->end_indices, poly, vertex)];
            if (n_neighbours[p] > 0)
                continue;

            n_neighbours[p] = get_neighbours_of_point(polygons,
                                                      poly, vertex, NULL, 0,
                                                      &interior);
            total_n_neighbours += n_neighbours[p];
        }
    }

    ALLOC(total_neighbours, total_n_neighbours);
    total_n_neighbours = 0;
    for (p = 0; p < polygons->n_points; p++)
    {
        neighbours[p] = &total_neighbours[total_n_neighbours];
        total_n_neighbours += n_neighbours[p];
    }

    for (p = 0; p < polygons->n_points; p++)
        n_neighbours[p] = 0;

    for (poly = 0; poly < polygons->n_items; poly++)
    {
        size = GET_OBJECT_SIZE(*polygons, poly);
        for (vertex = 0; vertex < size; vertex++)
        {
            p = polygons->indices[POINT_INDEX(polygons->end_indices, poly, vertex)];

            if (n_neighbours[p] > 0)
                continue;

            n_neighbours[p] = get_neighbours_of_point(polygons,
                                                      poly, vertex, neighbours[p],
                                                      total_n_neighbours,
                                                      &interior);
        }
    }

    *n_point_neighbours_ptr = n_neighbours;
    *point_neighbours_ptr = neighbours;
}

/**
 * \brief Free the neighbour lists returned by get_all_polygon_point_neighbours().
 *
 * The lists of all points are slices of one block that starts at
 * point_neighbours[0]. bicpl's delete_polygon_point_neighbours() frees every
 * list on its own and therefore aborts at the second one.
 *
 * \param n_point_neighbours (in)  int[n_points]; neighbour counts, freed
 * \param point_neighbours   (in)  int*[n_points]; neighbour lists, freed
 */
void delete_all_polygon_point_neighbours(int n_point_neighbours[],
                                         int *point_neighbours[])
{
    FREE(n_point_neighbours);
    if (point_neighbours != NULL)
    {
        FREE(point_neighbours[0]);
        FREE(point_neighbours);
    }
}

/**
 * \brief Evaluate heat diffusion kernel for a given distance.
 *
 * Computes the Gaussian heat kernel: exp(-(x*x) / (2*sigma*sigma))
 * Used in heat diffusion-based smoothing to weight influence based on distance.
 *
 * \param x     (in)  distance value
 * \param sigma (in)  kernel bandwidth parameter
 * \return            kernel weight in range (0, 1]
 */
static double
evaluate_heatkernel(double x, double sigma)
{
    return (exp(-(x * x) / (2.0 * sigma * sigma)));
}

/**
 * \brief Apply weighted averaging for a point using heat kernel distance weighting.
 *
 * Computes a weighted average of values at neighboring points on the mesh,
 * using heat kernel weights based on geodesic distance. Both point positions
 * and scalar values are smoothed.
 *
 * \param n_polygon_pts (in)  total number of points on mesh
 * \param polygon_pts   (in)  Point[n_polygon_pts]; 3D coordinates of all mesh points
 * \param values        (in)  double[n_polygon_pts]; scalar values at each point
 * \param n_neighbours  (in)  number of neighbors to consider
 * \param neighbours    (in)  int[n_neighbours]; neighbor vertex indices
 * \param ptidx         (in)  current point index
 * \param sigma         (in)  heat kernel bandwidth
 * \param smooth_point  (out) Point; smoothed position
 * \param value         (out) double; smoothed scalar value
 */
void heatkernel_blur_points(int n_polygon_pts, Point polygon_pts[],
                            double values[], int n_neighbours, int *neighbours,
                            int ptidx, double sigma, Point *smooth_point,
                            double *value)
{
    double sum[3], weight, sum_weight, point_dist;
    int i, c, n_pts, neigh;

    if (sigma <= 0.0)
        sigma = 1e-20;

    sum[0] = sum[1] = sum[2] = sum_weight = 0.0;

    for (i = 0; i < n_neighbours + 1; i++)
        *value = 1.0;

    for (i = 0; i < n_neighbours + 1; i++)
    {
        if (i > 0)
        { /* neighbouring points */
            neigh = neighbours[i - 1];
            if (values != NULL)
            {
                point_dist = distance_between_points(
                    &polygon_pts[ptidx],
                    &polygon_pts[neigh]);
            }
        }
        else
        { /* center point */
            point_dist = 0.0;
            neigh = ptidx;
        }

        if (values != NULL)
        {
            /* use Gaussian kernel for values */
            weight = evaluate_heatkernel(point_dist, sigma);

            /* only consider values that are not NaN */
            if (!isnan(values[neigh]))
            {
                sum[0] += weight * values[neigh];
                sum_weight += weight;
            }
        }
        else
        {
            /* this is rather based on empirically estimated values to be
               compatible to older versions */
            if (i > 0)
                weight = 1.0 / N_DIMENSIONS;
            else
                weight = 1.0;

            for (c = 0; c < N_DIMENSIONS; c++)
                sum[c] += weight * (double)
                                       Point_coord(polygon_pts[neigh], c);
            sum_weight += weight;
        }
    }

    if (values != NULL)
        *value = sum[0] / sum_weight;
    else
    {
        for (c = 0; c < N_DIMENSIONS; c++)
            Point_coord(*smooth_point, c) = (Point_coord_type)(sum[c] / sum_weight);
    }
}

/**
 * \brief Heat-kernel smoothing of values on a mesh whose points do not move.
 *
 * Smoothing values leaves the geometry, and with it every kernel weight,
 * unchanged, so the weights are computed once instead of in every iteration
 * (about 4x faster).  They and the order of summation are those of
 * heatkernel_blur_points(), so the result is bit-identical to it.
 *
 * \param polygons     (in)     mesh
 * \param values       (in/out) double[n_points]; smoothed in place, NaN kept
 * \param n_neighbours (in)     int[n_points]; neighbour count per point
 * \param neighbours   (in)     int*[n_points]; neighbour indices per point
 * \param sigma        (in)     heat kernel bandwidth
 * \param n_iter       (in)     number of iterations
 */
static void
smooth_values_heatkernel(polygons_struct *polygons, double *values,
                         int *n_neighbours, int **neighbours, double sigma,
                         int n_iter)
{
    double *weights, *smooth_values, centre, sum, sum_weight;
    int *offset, i, j, k, neigh, n_points = polygons->n_points;

    if (sigma <= 0.0)
        sigma = 1e-20;

    offset = (int *)malloc(sizeof(int) * (n_points + 1));
    offset[0] = 0;
    for (i = 0; i < n_points; i++)
        offset[i + 1] = offset[i] + n_neighbours[i];
    weights = (double *)malloc(sizeof(double) * (offset[n_points] + 1));
    for (i = 0; i < n_points; i++)
        for (k = 0; k < n_neighbours[i]; k++)
            weights[offset[i] + k] = evaluate_heatkernel(
                distance_between_points(&polygons->points[i],
                                        &polygons->points[neighbours[i][k]]),
                sigma);
    centre = evaluate_heatkernel(0.0, sigma);

    smooth_values = (double *)malloc(sizeof(double) * n_points);
    for (j = 0; j < n_iter; j++)
    {
        for (i = 0; i < n_points; i++)
        {
            if (isnan(values[i]))
            {
                smooth_values[i] = FNAN;
                continue;
            }
            /* accumulated from 0, as in heatkernel_blur_points(), so that
               even the sign of a zero result is the same */
            sum = sum_weight = 0.0;
            sum += centre * values[i];
            sum_weight += centre;
            for (k = 0; k < n_neighbours[i]; k++)
            {
                /* only consider values that are not NaN */
                neigh = neighbours[i][k];
                if (!isnan(values[neigh]))
                {
                    sum += weights[offset[i] + k] * values[neigh];
                    sum_weight += weights[offset[i] + k];
                }
            }
            smooth_values[i] = sum / sum_weight;
        }
        for (i = 0; i < n_points; i++)
            values[i] = smooth_values[i];
    }

    free(smooth_values);
    free(weights);
    free(offset);
}

/**
 * \brief Apply heat diffusion-based smoothing to polygon mesh scalar data.
 *
 * Iteratively applies Gaussian heat kernel smoothing to scalar values defined
 * on mesh vertices. The smoothing kernel bandwidth (sigma) is automatically
 * calculated from the mean vertex spacing. Multiple iterations can be specified
 * for stronger smoothing effect.
 *
 * \param polygons  (in/out) polygon mesh structure; coordinates smoothed in-place
 * \param values    (in/out) double[n_points]; scalar values; smoothed in-place
 * \param fwhm      (in)     Gaussian smoothing FWHM in mm (converted to sigma internally)
 */
void smooth_heatkernel(polygons_struct *polygons, double *values, double fwhm)
{
    double sigma, value, point_dist, sum_dist;
    Point point, *smooth_pts;
    int n_iter, i, j, size, p1, p2, n;
    int *n_neighbours, **neighbours;
    BOOLEAN values_present;
    progress_struct progress;

    get_all_polygon_point_neighbours(polygons, &n_neighbours, &neighbours);

    /* estimate mean distance between the points to define optimal sigma */
    n = 0;
    sum_dist = 0.0;
    for (i = 0; i < polygons->n_items; i++)
    {
        size = GET_OBJECT_SIZE(*polygons, i);

        for (j = 0; j < size; j++)
        {
            p1 = polygons->indices[POINT_INDEX(polygons->end_indices, i, j)];
            p2 = polygons->indices[POINT_INDEX(polygons->end_indices, i, (j + 1) % size)];

            point_dist = distance_between_points(
                &polygons->points[p1],
                &polygons->points[p2]);
            n++;
            sum_dist += point_dist;
        }
    }

    /* use mean distance as sigma */
    sigma = sum_dist / (double)n;

    /* calculate n_iter with regard to fwhm
       0.541011 equals to 3.0/(8.0*log(2.0))
       see SurfStatSmooth.m in surfstat from Keith Worsley */
    n_iter = ceil(fwhm * fwhm * 0.541011 / (sigma * sigma));
    if (n_iter < 1)
        n_iter = 1;

    /* recalibrate sigma to ensure integer numbers for iterations */
    sigma = (fwhm * 0.7355345) / sqrt((double)n_iter);

    initialize_progress_report(&progress, FALSE, n_iter * polygons->n_points,
                               "Blurring");
    if (values != NULL)
        smooth_values_heatkernel(polygons, values, n_neighbours, neighbours,
                                 sigma, n_iter);
    else
    {
        smooth_pts = (Point *)malloc(sizeof(Point) * polygons->n_points);

        /* diffusion smoothing of the points themselves */
        for (j = 0; j < n_iter; j++)
        {
            for (i = 0; i < polygons->n_points; i++)
            {
                heatkernel_blur_points(polygons->n_points,
                                       polygons->points, NULL,
                                       n_neighbours[i], neighbours[i],
                                       i, sigma, &point, &value);
                smooth_pts[i] = point;
            }
            for (i = 0; i < polygons->n_points; i++)
                polygons->points[i] = smooth_pts[i];
        }
        free(smooth_pts);
    }

    terminate_progress_report(&progress);

    delete_all_polygon_point_neighbours(n_neighbours, neighbours);
}

/**
 * \brief HC Laplacian smoothing of a mesh (Vollmer et al., 1999).
 *
 * Laplacian smoothing with Humphrey's Classes: after each Laplacian step the
 * points are pushed back towards a blend of the original and the previous
 * positions, which keeps the volume from shrinking. Border edges are locked.
 * The mesh is fan-triangulated if it is not already a triangle mesh.
 *
 * \param polygons (in/out) mesh, smoothed in-place
 * \param iter     (in)     number of iterations
 * \param alpha    (in)     weight of the original positions in the correction
 *                          (0..1; 0.1 is typical)
 * \param beta     (in)     weight of a point's own correction against that of
 *                          its neighbours (0..1; 0.5 is typical)
 * \return 0 on success, -1 on an empty mesh or allocation failure
 */
int smooth_laplacian(polygons_struct *polygons,
                     int iter,
                     double alpha,
                     double beta)
{
    if (!polygons || polygons->n_points <= 0 || polygons->n_items <= 0)
        return -1;

    vec3d *V = NULL;
    vec3i *F = NULL;
    int nv = 0, nf = 0, fan = 0;

    /* --- in: BICPL -> (V,F) triangles --- */
    if (polygons_to_tri_arrays(polygons, &V, &F, &nv, &nf, &fan) != 0)
        return -1;

    /* --- simplify --- */
    laplacian_smoothHC(V, F, nv, nf, alpha, beta, iter, true);

    if (nv <= 0 || nf <= 0)
    {
        free(V);
        free(F);
        return -1;
    }

    /* --- out: (V,F) -> BICPL (triangles, exclusive end_indices) --- */
    const int ok = tri_arrays_to_polygons(polygons, V, F, nv, nf);
    free(V);
    free(F);
    if (ok != 0)
        return -1;

    return 0;
}
