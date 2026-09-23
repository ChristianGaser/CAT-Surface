/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 *
 */

#include <bicpl.h>
#include <float.h>

#include <string.h>

#include "CAT_Warp.h"
#include "CAT_Map.h"
#include "CAT_Surf.h"
#include "CAT_SurfUtils.h"
#include "CAT_Curvature.h"
#include "CAT_Smooth.h"
#include "CAT_Resample.h"
#include "CAT_SafeAlloc.h"

/**
 * \brief Apply 3x3 rotation matrix to all vertices of a mesh.
 *
 * Transforms mesh coordinates via matrix multiplication: each vertex is rotated
 * by the 9-element rotation matrix (stored in row-major order). Can modify in-place
 * or write to separate output mesh. Used for surface reparameterization and alignment.
 *
 * \param polygons (in) source mesh
 * \param rotated_polygons (out) output mesh (if NULL, modifies input in-place)
 * \param rotation_matrix (in) 9-element row-major 3x3 rotation matrix
 */
void rotate_polygons(polygons_struct *polygons, polygons_struct *rotated_polygons,
                     double *rotation_matrix)
{
    int i;
    double x, y, z;

    if (rotated_polygons != NULL)
        copy_polygons(polygons, rotated_polygons);

    for (i = 0; i < polygons->n_points; i++)
    {
        x = Point_x(polygons->points[i]) * rotation_matrix[0] + Point_y(polygons->points[i]) * rotation_matrix[1] + Point_z(polygons->points[i]) * rotation_matrix[2];
        y = Point_x(polygons->points[i]) * rotation_matrix[3] + Point_y(polygons->points[i]) * rotation_matrix[4] + Point_z(polygons->points[i]) * rotation_matrix[5];
        z = Point_x(polygons->points[i]) * rotation_matrix[6] + Point_y(polygons->points[i]) * rotation_matrix[7] + Point_z(polygons->points[i]) * rotation_matrix[8];
        if (rotated_polygons != NULL)
        {
            fill_Point(rotated_polygons->points[i], x, y, z);
        }
        else
            fill_Point(polygons->points[i], x, y, z);
    }
}

/**
 * \brief Compute 3x3 rotation matrix from Euler angles (alpha, beta, gamma).
 *
 * Combines three rotation matrices (X, Y, Z axes) via sequential matrix multiplications:
 * R_total = R_z(gamma) * R_y(beta) * R_x(alpha). Outputs 9-element row-major matrix.
 * Standard aerospace/ZYX Euler angle convention.
 *
 * \param rotation_matrix (out) 9-element row-major result matrix
 * \param alpha (in) rotation angle about X-axis (radians)
 * \param beta (in) rotation angle about Y-axis (radians)
 * \param gamma (in) rotation angle about Z-axis (radians)
 */
void rotation_to_matrix(double *rotation_matrix, double alpha, double beta,
                        double gamma)
{
    int i, j, k;
    double sum, rot[9];

    /* rotation matrices */
    double rot_x[9] = {1.0, 0.0, 0.0,
                       0.0, cos(alpha), sin(alpha),
                       0.0, -sin(alpha), cos(alpha)};
    double rot_y[9] = {cos(beta), 0.0, sin(beta),
                       0.0, 1.0, 0.0,
                       -sin(beta), 0.0, cos(beta)};
    double rot_z[9] = {cos(gamma), sin(gamma), 0.0,
                       -sin(gamma), cos(gamma), 0.0,
                       0.0, 0.0, 1.0};

    /* combine x and y rotation */
    for (i = 0; i < 3; i++)
    {
        for (j = 0; j < 3; j++)
        {
            sum = 0.0;
            for (k = 0; k < 3; k++)
                sum += rot_y[i + 3 * k] * rot_x[k + 3 * j];
            rot[i + 3 * j] = sum;
        }
    }

    /* combine with z rotation */
    for (i = 0; i < 3; i++)
    {
        for (j = 0; j < 3; j++)
        {
            sum = 0.0;
            for (k = 0; k < 3; k++)
                sum += rot_z[i + 3 * k] * rot[k + 3 * j];
            rotation_matrix[i + 3 * j] = sum;
        }
    }
}

/**
 * \brief Apply 2D UV-space deformation directly using separated u,v displacement fields.
 *
 * Like apply_warp() but takes separate u and v displacement arrays for flexibility.
 * Better for workflows that compute displacement components independently. Still uses
 * spherical triangle interpolation for mapping to 3D mesh.
 *
 * \param polygons (in/out) mesh modified by warp
 * \param sphere (in) spherical reference surface
 * \param ux (in) u-component displacement field
 * \param vy (in) v-component displacement field (same size as ux)
 * \param inverse (in) 1 for inverse warp; 0 for forward
 */
void apply_uv_warp(polygons_struct *polygons, polygons_struct *sphere, double *ux,
                   double *vy, int inverse)
{
    Point centre, unit_point, *new_points, trans_point;
    polygons_struct unit_sphere;
    double u, v, x, y, z;
    double indx, indy;
    double xo, yo, zo;
    int i, p, ind;

    copy_polygons(sphere, &unit_sphere);

    create_polygons_bintree(polygons, round((double)polygons->n_items *
                                            BINTREE_FACTOR));
    create_polygons_bintree(&unit_sphere,
                            round((double)unit_sphere.n_items *
                                  BINTREE_FACTOR));

    ALLOC(new_points, sphere->n_points);

    for (p = 0; p < polygons->n_points; p++)
    {
        xo = Point_x(sphere->points[p]);
        yo = Point_y(sphere->points[p]);
        zo = Point_z(sphere->points[p]);

        if (inverse)
        {
            u = -ux[p];
            v = -vy[p];
        }
        else
        {
            u = ux[p];
            v = vy[p];
        }

        x = xo * cos(u) + xo * cos(v) + yo * sin(v) + zo * sin(u) - xo;
        y = xo * -sin(u) * sin(u) + yo * cos(u) + zo * cos(u) * sin(u) + xo * -sin(v) + yo * cos(v) - yo;
        z = xo * -sin(u) * cos(u) + yo * -sin(u) + zo * cos(u) * cos(u);

        fill_Point(trans_point, x, y, z);

        map_unit_sphere_to_point(&unit_sphere, &trans_point,
                                 polygons, &new_points[p]);
    }

    for (p = 0; p < polygons->n_points; p++)
        polygons->points[p] = new_points[p];

    /* set radius to 1 */
    for (i = 0; i < unit_sphere.n_points; i++)
        set_vector_length(&unit_sphere.points[i], 1.0);

    compute_polygon_normals(polygons);
    free(new_points);
}

/* Row-major 3x3 matrix product C = A * B (C must not alias A or B). */
static void
mat3_mul(double *C, const double *A, const double *B)
{
    int r, c, k;

    for (r = 0; r < 3; r++)
        for (c = 0; c < 3; c++)
        {
            double s = 0.0;
            for (k = 0; k < 3; k++)
                s += A[r * 3 + k] * B[k * 3 + c];
            C[r * 3 + c] = s;
        }
}

/* Curvature sum-of-squared-differences cost of applying the row-major 3x3
   rotation R to the source sphere: resample the target feature onto the rotated
   source grid and accumulate squared differences against the source feature.

   Note the direction this establishes: the cost is lowest when the SOURCE
   feature at s matches the TARGET feature at R*s, so the minimising R maps
   source coordinates onto target coordinates. Callers that deform the source
   (CAT_WarpDemonsRegister) need the opposite direction and must apply R^-1. */
static double
cost_for_rotation(OptimizationParams *p, const double *R)
{
    polygons_struct rot_src_sphere;
    double sum_sq = 0.0, d;
    int i;

    rotate_polygons(p->src_sphere, &rot_src_sphere, (double *) R);
    resample_values_sphere(p->trg_sphere, &rot_src_sphere,
                           p->orig_trg, p->map_trg, 0, 0);
    for (i = 0; i < p->src->n_points; i++)
    {
        d = p->map_src[i] - p->map_trg[i];
        sum_sq += d * d;
    }
    delete_polygons(&rot_src_sphere);
    return sum_sq;
}

static double compute_cost(double *angles, void *params)
{
    OptimizationParams *p = (OptimizationParams *) params;
    double R[9];

    rotation_to_matrix(R, angles[0], angles[1], angles[2]);

    /* When a seed rotation is supplied, the optimiser searches the residual on
       top of it: net rotation = R_residual * R_seed. */
    if (p->pre_rot != NULL)
    {
        double Rc[9];
        mat3_mul(Rc, R, p->pre_rot);
        return cost_for_rotation(p, Rc);
    }
    return cost_for_rotation(p, R);
}

static void nelder_mead(double **simplex, double *f_values, int n, int max_iter, double tol, OptimizationParams *params, double *optimal_params, int verbose)
{
    int i, j, iter;
    int highest, second_highest, lowest;
    double centroid[n];
    double reflected[n], expanded[n], contracted[n];
    double f_reflected, f_expanded, f_contracted;

    for (iter = 0; iter < max_iter; iter++)
    {
        // Identify the lowest, highest, and second-highest points
        highest = 0;
        lowest = 0;
        second_highest = 1;
        for (i = 0; i <= n; i++)
        {
            if (f_values[i] > f_values[highest])
            {
                second_highest = highest;
                highest = i;
            }
            else if (f_values[i] > f_values[second_highest] && i != highest)
            {
                second_highest = i;
            }
            if (f_values[i] < f_values[lowest])
            {
                lowest = i;
            }
        }

        // Compute the centroid
        for (j = 0; j < n; j++)
        {
            centroid[j] = 0.0;
            for (i = 0; i <= n; i++)
            {
                if (i != highest)
                {
                    centroid[j] += simplex[i][j];
                }
            }
            centroid[j] /= n;
        }

        // Reflection
        for (j = 0; j < n; j++)
        {
            reflected[j] = centroid[j] + ALPHA * (centroid[j] - simplex[highest][j]);
        }
        f_reflected = compute_cost(reflected, params);

        if (f_reflected < f_values[lowest])
        {
            // Expansion
            for (j = 0; j < n; j++)
            {
                expanded[j] = centroid[j] + GAMMA * (reflected[j] - centroid[j]);
            }
            f_expanded = compute_cost(expanded, params);

            if (f_expanded < f_reflected)
            {
                for (j = 0; j < n; j++)
                {
                    simplex[highest][j] = expanded[j];
                }
                f_values[highest] = f_expanded;
            }
            else
            {
                for (j = 0; j < n; j++)
                {
                    simplex[highest][j] = reflected[j];
                }
                f_values[highest] = f_reflected;
            }
        }
        else if (f_reflected < f_values[second_highest])
        {
            for (j = 0; j < n; j++)
            {
                simplex[highest][j] = reflected[j];
            }
            f_values[highest] = f_reflected;
        }
        else
        {
            // Contraction
            for (j = 0; j < n; j++)
            {
                contracted[j] = centroid[j] + RHO * (simplex[highest][j] - centroid[j]);
            }
            f_contracted = compute_cost(contracted, params);

            if (f_contracted < f_values[highest])
            {
                for (j = 0; j < n; j++)
                {
                    simplex[highest][j] = contracted[j];
                }
                f_values[highest] = f_contracted;
            }
            else
            {
                // Shrink the simplex
                for (i = 0; i <= n; i++)
                {
                    if (i != lowest)
                    {
                        for (j = 0; j < n; j++)
                        {
                            simplex[i][j] = simplex[lowest][j] + SIGMA * (simplex[i][j] - simplex[lowest][j]);
                        }
                        f_values[i] = compute_cost(simplex[i], params);
                    }
                }
            }
        }

        // Check for convergence
        double max_diff = 0.0;
        for (i = 0; i <= n; i++)
        {
            double diff = fabs(f_values[i] - f_values[lowest]);
            if (diff > max_diff)
            {
                max_diff = diff;
            }
        }
        if (max_diff < tol)
        {
            break;
        }
    }

    // Copy the optimal parameters from the lowest point
    for (j = 0; j < n; j++)
    {
        optimal_params[j] = simplex[lowest][j];
    }

    if (verbose)
    {
        printf("Optimization completed in %d iterations\n", iter);
        printf("Minimum found at:\n");
        for (j = 0; j < n; j++)
        {
            printf("Angle[%d] = %.6f\n", j, optimal_params[j]);
        }
        printf("Minimum squared difference: %.6f\n", f_values[lowest]);
    }
}

/* Cap on consecutive re-centring passes at one angular scale, bounding runtime
   when the optimum keeps drifting (FreeSurfer relies on a visited-flag cache
   instead; here a hard cap is simpler and equally effective). */
#define ROT_MAX_RECENTER 4

/**
 * \brief Exhaustive coarse-to-fine global search for the initial rigid rotation.
 *
 * See the header declaration for the full description. This mirrors FreeSurfer's
 * MRISrigidBodyAlignGlobal: rather than trusting a local optimiser, it evaluates
 * the cost on a dense grid over all three rotation angles, re-centres on the best
 * candidate, and halves the angular span only once a pass yields no improvement -
 * repeating until the span drops below \p min_degrees. Because every candidate in
 * the span is evaluated, the search cannot be trapped in a neighbouring-sulcus
 * local minimum the way a simplex started from any single seed can, and its
 * capture range is the full \p max_degrees rather than a basin width.
 *
 * \param src             (in)  source central surface mesh
 * \param src_sphere      (in)  spherical parameterization of \p src
 * \param trg             (in)  template central surface mesh
 * \param trg_sphere      (in)  spherical parameterization of \p trg
 * \param fwhm            (in)  smoothing FWHM for the curvature cost feature
 * \param curvtype        (in)  curvature type for the cost feature
 * \param max_degrees     (in)  half-width of the initial angular search span
 * \param min_degrees     (in)  stop once the span falls below this
 * \param nangles         (in)  grid samples per axis per pass (span/nangles step)
 * \param refine          (in)  1 = Nelder-Mead refine of the residual afterwards
 * \param rotation_matrix (out) 9-element row-major rotation for the source sphere
 * \param verbose         (in)  1 for progress output; 0 silent
 */
void
rotate_polygons_to_atlas_global(polygons_struct *src, polygons_struct *src_sphere,
                                polygons_struct *trg, polygons_struct *trg_sphere,
                                double fwhm, int curvtype,
                                double max_degrees, double min_degrees,
                                int nangles, int refine,
                                double *rotation_matrix, int verbose)
{
    double *orig_trg, *map_trg, *map_src;
    double best_ang[3] = {0.0, 0.0, 0.0}, best_cost, cost_id;
    double degrees;
    int    recenter = 0, n_eval = 0, i;

    orig_trg = SAFE_MALLOC(double, trg->n_points);
    map_trg  = SAFE_MALLOC(double, src->n_points);
    map_src  = SAFE_MALLOC(double, src->n_points);

    get_smoothed_curvatures(trg, orig_trg, fwhm, curvtype);
    get_smoothed_curvatures(src, map_src, fwhm, curvtype);

    OptimizationParams params = {src, src_sphere, trg_sphere, orig_trg, map_trg,
                                 map_src, NULL};

    if (nangles < 2)
        nangles = 2;
    if (min_degrees <= 0.0)
        min_degrees = 0.5;
    if (max_degrees < min_degrees)
        max_degrees = min_degrees;

    cost_id = compute_cost(best_ang, &params);
    best_cost = cost_id;
    n_eval++;

    if (verbose)
        fprintf(stdout, "Global rotation search: +/-%.3g deg down to %.3g deg, "
                "%d samples/axis (identity cost %.4g)\n",
                max_degrees, min_degrees, nangles, cost_id);

    for (degrees = max_degrees; degrees >= min_degrees; ) {
        double step = 2.0 * degrees / (double) nangles;
        double ctr[3];
        int improved = 0, ai, bi, gi;

        for (i = 0; i < 3; i++)
            ctr[i] = best_ang[i];

        for (gi = 0; gi <= nangles; gi++)
            for (bi = 0; bi <= nangles; bi++)
                for (ai = 0; ai <= nangles; ai++) {
                    double cand[3], cost;

                    cand[0] = ctr[0] + RADIANS((ai - nangles / 2.0) * step);
                    cand[1] = ctr[1] + RADIANS((bi - nangles / 2.0) * step);
                    cand[2] = ctr[2] + RADIANS((gi - nangles / 2.0) * step);

                    cost = compute_cost(cand, &params);
                    n_eval++;
                    if (cost < best_cost) {
                        best_cost = cost;
                        for (i = 0; i < 3; i++)
                            best_ang[i] = cand[i];
                        improved = 1;
                    }
                }

        if (verbose)
            fprintf(stdout, "  span %+.3g deg (step %.3g): best cost %.4g at "
                    "(%.2f, %.2f, %.2f) deg%s\n", degrees, step, best_cost,
                    DEGREES(best_ang[0]), DEGREES(best_ang[1]),
                    DEGREES(best_ang[2]), improved ? "" : " [no change]");

        /* Follow FreeSurfer: only shrink the span once a full pass finds nothing
           better; otherwise re-centre and re-scan at the same scale (capped). */
        if (improved && ++recenter < ROT_MAX_RECENTER)
            continue;
        recenter = 0;
        degrees *= 0.5;
    }

    if (verbose)
        fprintf(stdout, "Global rotation: (%.2f, %.2f, %.2f) deg, cost %.4g "
                "(identity %.4g), %d evaluations\n",
                DEGREES(best_ang[0]), DEGREES(best_ang[1]), DEGREES(best_ang[2]),
                best_cost, cost_id, n_eval);

    /* Optional local refine of the residual on top of the global optimum. */
    if (refine) {
        double Rseed[9];
        double simplex[4][3] = {
            {0.0, 0.0, 0.0},
            {0.1, 0.0, 0.0},
            {0.0, 0.1, 0.0},
            {0.0, 0.0, 0.1}};
        double *simplex_ptrs[4] = {simplex[0], simplex[1], simplex[2], simplex[3]};
        double f_values[4], resid[3], Rres[9];

        rotation_to_matrix(Rseed, best_ang[0], best_ang[1], best_ang[2]);
        params.pre_rot = Rseed;
        for (i = 0; i < 4; i++)
            f_values[i] = compute_cost(simplex_ptrs[i], &params);
        nelder_mead(simplex_ptrs, f_values, 3, MAX_ITER, TOL, &params, resid,
                    verbose);
        rotation_to_matrix(Rres, resid[0], resid[1], resid[2]);
        mat3_mul(rotation_matrix, Rres, Rseed);
        params.pre_rot = NULL;
    } else
        rotation_to_matrix(rotation_matrix, best_ang[0], best_ang[1],
                           best_ang[2]);

    free(orig_trg);
    free(map_trg);
    free(map_src);
}
