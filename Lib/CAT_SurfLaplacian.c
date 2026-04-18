/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Laplacian streamline-based pial/white surface placement.
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 *
 */

#include <bicpl.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <float.h>

#include "CAT_NiftiLib.h"
#include "CAT_Vol.h"
#include "CAT_Smooth.h"
#include "CAT_SurfLaplacian.h"

/* Voxel classification for the Laplace solver */
#define MASK_BG 0     /* background / outside brain        */
#define MASK_RIBBON 1 /* cortical ribbon  (solve here)     */
#define MASK_PIAL 2   /* pial boundary    (phi = 1 fixed)  */
#define MASK_WHITE 3  /* white boundary   (phi = 0 fixed)  */

/* -------------------------------------------------------------------
 * solve_laplace_ribbon
 * -------------------------------------------------------------------*/

/**
 * \brief Solve the Laplace equation in the cortical ribbon using SOR.
 *
 * Boundary conditions:
 *   phi = 1  where 0.5 < label <= lim_pial  (CSF / pial side)
 *   phi = 0  where label >= lim_white        (WM / white side)
 * The cortical ribbon is the set of voxels with lim_pial < label < lim_white.
 *
 * Uses Successive Over-Relaxation (SOR) with omega = 1.7 for fast
 * convergence.  The initial guess is a linear interpolation of the
 * label value inside the ribbon, which is close to the true solution
 * and accelerates convergence.
 *
 * \param labels      (in)     tissue label volume
 * \param dims        (in)     volume dimensions [nx, ny, nz]
 * \param ribbon_pial (in)     pial-side ribbon boundary (phi = 1 fixed)
 * \param ribbon_white(in)     white-side ribbon boundary (phi = 0 fixed)
 * \param phi         (out)    solution volume (pre-allocated, same dims)
 * \param max_iter    (in)     maximum SOR iterations
 * \param tol         (in)     convergence tolerance (max abs change)
 * \param verbose     (in)     print progress
 */
static void
solve_laplace_ribbon(const float *labels, int dims[3],
                     float ribbon_pial, float ribbon_white,
                     float *phi, int max_iter, float tol,
                     int verbose)
{
    int nx = dims[0], ny = dims[1], nz = dims[2];
    int nxy = nx * ny;
    int nvox = nxy * nz;
    int x, y, z, idx, iter, n;
    unsigned char *mask;
    float denom, sum, gs, new_val, change, max_change;
    const float omega = 1.7f;

    mask = (unsigned char *)calloc(nvox, sizeof(unsigned char));
    if (!mask)
    {
        fprintf(stderr, "solve_laplace_ribbon: memory allocation error\n");
        return;
    }

    /* ---- Classify voxels and initialise phi ---- */
    denom = ribbon_pial - ribbon_white; /* negative, used for linear init */
    for (idx = 0; idx < nvox; idx++)
    {
        float val = labels[idx];
        if (val > 0.5f && val <= ribbon_pial)
        {
            mask[idx] = MASK_PIAL;
            phi[idx] = 1.0f;
        }
        else if (val >= ribbon_white)
        {
            mask[idx] = MASK_WHITE;
            phi[idx] = 0.0f;
        }
        else if (val > ribbon_pial && val < ribbon_white)
        {
            mask[idx] = MASK_RIBBON;
            /* Linear interpolation: phi = 1 at ribbon_pial, phi = 0 at ribbon_white */
            phi[idx] = (val - ribbon_white) / denom;
        }
        else
        {
            mask[idx] = MASK_BG;
            phi[idx] = 0.0f;
        }
    }

    /* ---- SOR iteration ---- */
    for (iter = 0; iter < max_iter; iter++)
    {
        max_change = 0.0f;

        for (z = 1; z < nz - 1; z++)
        {
            for (y = 1; y < ny - 1; y++)
            {
                idx = 1 + y * nx + z * nxy;
                for (x = 1; x < nx - 1; x++, idx++)
                {
                    if (mask[idx] != MASK_RIBBON)
                        continue;

                    /* 6-connectivity average, skipping background */
                    sum = 0.0f;
                    n = 0;

                    if (mask[idx - 1]  > 0) { sum += phi[idx - 1];  n++; }
                    if (mask[idx + 1]  > 0) { sum += phi[idx + 1];  n++; }
                    if (mask[idx - nx] > 0) { sum += phi[idx - nx]; n++; }
                    if (mask[idx + nx] > 0) { sum += phi[idx + nx]; n++; }
                    if (mask[idx - nxy] > 0) { sum += phi[idx - nxy]; n++; }
                    if (mask[idx + nxy] > 0) { sum += phi[idx + nxy]; n++; }

                    if (n == 0)
                        continue;

                    gs = sum / (float)n;
                    new_val = (1.0f - omega) * phi[idx] + omega * gs;
                    change = fabsf(new_val - phi[idx]);
                    if (change > max_change)
                        max_change = change;
                    phi[idx] = new_val;
                }
            }
        }

        if (verbose && ((iter + 1) % 200 == 0 || max_change < tol))
            fprintf(stdout, "\rLaplace solver: iter %4d/%d  max_change %.2e    ",
                    iter + 1, max_iter, max_change);

        if (max_change < tol)
            break;
    }

    if (verbose)
        fprintf(stdout, "\rLaplace solver: converged in %d iterations (max_change %.2e)\n",
                iter < max_iter ? iter + 1 : max_iter,
                max_change);

    free(mask);
}

/* -------------------------------------------------------------------
 * solve_ade_ribbon  — Adaptive Diffusion Equation solver
 * -------------------------------------------------------------------*/

/**
 * \brief Solve the Adaptive Diffusion Equation in the cortical ribbon.
 *
 * Solves  div( f(x) grad(phi) ) = 0  where f(x) is a gray-matter
 * fraction derived from the tissue label volume.  Boundary conditions
 * are identical to the Laplace solver (phi=1 pial, phi=0 white).
 *
 * The diffusivity f is set proportional to the GM fraction:
 *   - pure GM (label == 2.0): f = 1
 *   - partial volume voxels:  f linearly interpolated between the
 *     boundary labels and 1.0, clamped to [epsilon, 1]
 *   - WM/CSF boundaries:      f = large value (perfect conductor)
 *
 * At each face between voxels i and j, the effective conductance is
 * the harmonic mean  2*f_i*f_j / (f_i + f_j).
 *
 * Reference:
 *   Joshi AA et al. (2025). Robust Cortical Thickness Estimation in
 *   the Presence of Partial Volumes using Adaptive Diffusion Equation.
 *   J Neurosci Methods 423:110552.
 *
 * \param labels      (in)     tissue label volume
 * \param dims        (in)     volume dimensions [nx, ny, nz]
 * \param ribbon_pial (in)     pial-side ribbon boundary (phi = 1 fixed)
 * \param ribbon_white(in)     white-side ribbon boundary (phi = 0 fixed)
 * \param phi         (out)    solution volume (pre-allocated)
 * \param max_iter    (in)     maximum SOR iterations
 * \param tol         (in)     convergence tolerance
 * \param verbose     (in)     print progress
 */
static void
solve_ade_ribbon(const float *labels, int dims[3],
                 float ribbon_pial, float ribbon_white,
                 float *phi, int max_iter, float tol,
                 int verbose)
{
    int nx = dims[0], ny = dims[1], nz = dims[2];
    int nxy = nx * ny;
    int nvox = nxy * nz;
    int x, y, z, idx, iter;
    unsigned char *mask;
    float *frac; /* GM fraction / diffusivity per voxel */
    float denom, wsum, wtotal, fi, fj, w, gs, new_val, change, max_change;
    const float omega = 1.7f;
    const float eps = 0.01f;       /* minimum diffusivity (avoids division by 0) */
    const float conductor = 10.0f; /* high conductivity for boundary voxels */

    mask = (unsigned char *)calloc(nvox, sizeof(unsigned char));
    frac = (float *)malloc(sizeof(float) * nvox);
    if (!mask || !frac)
    {
        fprintf(stderr, "solve_ade_ribbon: memory allocation error\n");
        if (mask)
            free(mask);
        if (frac)
            free(frac);
        return;
    }

    /* ---- Classify voxels, compute GM fraction, initialise phi ---- */
    denom = ribbon_pial - ribbon_white;
    for (idx = 0; idx < nvox; idx++)
    {
        float val = labels[idx];
        if (val > 0.5f && val <= ribbon_pial)
        {
            mask[idx] = MASK_PIAL;
            phi[idx] = 1.0f;
            frac[idx] = conductor;
        }
        else if (val >= ribbon_white)
        {
            mask[idx] = MASK_WHITE;
            phi[idx] = 0.0f;
            frac[idx] = conductor;
        }
        else if (val > ribbon_pial && val < ribbon_white)
        {
            mask[idx] = MASK_RIBBON;
            phi[idx] = (val - ribbon_white) / denom;

            /* GM fraction: peak at GM label (2.0), fall off toward boundaries.
             * Map label linearly: 2.0 -> 1.0,  ribbon_pial -> 0.25,
             * ribbon_white -> 0.25.  Clamp to [eps, 1]. */
            float gm_label = 2.0f;
            float half_range = (ribbon_white - ribbon_pial) * 0.5f;
            float dist_from_gm = fabsf(val - gm_label) / half_range;
            frac[idx] = fmaxf(eps, 1.0f - 0.75f * dist_from_gm);
        }
        else
        {
            mask[idx] = MASK_BG;
            phi[idx] = 0.0f;
            frac[idx] = eps;
        }
    }

    /* ---- SOR iteration with weighted Laplacian ---- */
    for (iter = 0; iter < max_iter; iter++)
    {
        max_change = 0.0f;

        for (z = 1; z < nz - 1; z++)
        {
            for (y = 1; y < ny - 1; y++)
            {
                idx = 1 + y * nx + z * nxy;
                for (x = 1; x < nx - 1; x++, idx++)
                {
                    if (mask[idx] != MASK_RIBBON)
                        continue;

                    fi = frac[idx];
                    wsum = 0.0f;
                    wtotal = 0.0f;

                    /* 6-connected neighbors, harmonic mean conductance */
#define ADE_NEIGHBOR(OFF)                            \
    do                                               \
    {                                                \
        int nb = (OFF);                              \
        if (mask[idx + nb] > 0)                      \
        {                                            \
            fj = frac[idx + nb];                     \
            w = 2.0f * fi * fj / (fi + fj + 1e-30f); \
            wsum += w * phi[idx + nb];               \
            wtotal += w;                             \
        }                                            \
    } while (0)

                    ADE_NEIGHBOR(-1);
                    ADE_NEIGHBOR(+1);
                    ADE_NEIGHBOR(-nx);
                    ADE_NEIGHBOR(+nx);
                    ADE_NEIGHBOR(-nxy);
                    ADE_NEIGHBOR(+nxy);
#undef ADE_NEIGHBOR

                    if (wtotal < 1e-20f)
                        continue;

                    gs = wsum / wtotal;
                    new_val = (1.0f - omega) * phi[idx] + omega * gs;
                    change = fabsf(new_val - phi[idx]);
                    if (change > max_change)
                        max_change = change;
                    phi[idx] = new_val;
                }
            }
        }

        if (verbose && ((iter + 1) % 200 == 0 || max_change < tol))
            fprintf(stdout, "\rADE solver: iter %4d/%d  max_change %.2e    ",
                    iter + 1, max_iter, max_change);

        if (max_change < tol)
            break;
    }

    if (verbose)
        fprintf(stdout, "\rADE solver: converged in %d iterations (max_change %.2e)\n",
                iter < max_iter ? iter + 1 : max_iter,
                max_change);

    free(mask);
    free(frac);
}

/* -------------------------------------------------------------------
 * surf_laplacian_pial_white  /  surf_ade_pial_white
 *
 * Both share the same streamline tracing logic; only the PDE solver
 * differs (Laplace vs ADE).
 * -------------------------------------------------------------------*/

/**
 * \brief Internal: solve PDE + trace streamlines.
 *
 * The PDE is always solved on a wide ribbon from CSF (1.0) to WM (3.0)
 * so the phi field has enough room for smooth gradients.  The target
 * isovalues lim_pial and lim_white are mapped to phi stop values:
 *   phi_stop = (target - ribbon_white) / (ribbon_pial - ribbon_white)
 * This means changing lim_pial/lim_white actually moves the surfaces.
 *
 * \param use_ade  0 = Laplace equation,  1 = Adaptive Diffusion Equation
 */
static int
surf_pde_pial_white_internal(polygons_struct *central,
                             float *labels,
                             nifti_image *nii_ptr,
                             float lim_pial,
                             float lim_white,
                             const double *thickness_values,
                             polygons_struct *pial_out,
                             polygons_struct *white_out,
                             int use_ade,
                             int verbose)
{
    int v, step, dims[3], nvox;
    double vx[3];
    int n_points;
    float *phi = NULL;
    float *grad_x = NULL, *grad_y = NULL, *grad_z = NULL;

    /* Wide ribbon boundaries for the PDE solver (CSF to WM) */
    const float ribbon_pial  = 1.5f;   /* CSF label — pial side (phi = 1) */
    const float ribbon_white = 2.5f;   /* WM label  — white side (phi = 0) */

    /* Map target isovalues to phi stop thresholds.
     * phi = 1 at ribbon_pial (1.0), phi = 0 at ribbon_white (3.0).
     * Linear mapping:  phi = (target - ribbon_white) / (ribbon_pial - ribbon_white)
     * Example: lim_pial=1.5 -> phi_stop_pial = (1.5 - 3.0) / (1.0 - 3.0) = 0.75
     *          lim_white=2.5 -> phi_stop_white = (2.5 - 3.0) / (1.0 - 3.0) = 0.25 */
    float phi_stop_pial  = (lim_pial  - ribbon_white) / (ribbon_pial - ribbon_white);
    float phi_stop_white = (lim_white - ribbon_white) / (ribbon_pial - ribbon_white);

    /* Clamp to sensible range */
    if (phi_stop_pial  > 0.99f) phi_stop_pial  = 0.999f;
    if (phi_stop_pial  < 0.01f) phi_stop_pial  = 0.50f;
    if (phi_stop_white > 0.99f) phi_stop_white = 0.50f;
    if (phi_stop_white < 0.01f) phi_stop_white = 0.001f;

    /* Streamline parameters */
    const double step_size = 0.1;  /* mm per integration step          */
    const int    max_steps = 120;  /* max steps = 12 mm travel         */
    const double min_grad  = 1e-6; /* gradient magnitude floor         */

    if (!central || !labels || !nii_ptr || !pial_out || !white_out)
        return -1;

    n_points = central->n_points;
    dims[0] = nii_ptr->nx;
    dims[1] = nii_ptr->ny;
    dims[2] = nii_ptr->nz;
    nvox = dims[0] * dims[1] * dims[2];
    vx[0] = fabs(nii_ptr->dx);
    vx[1] = fabs(nii_ptr->dy);
    vx[2] = fabs(nii_ptr->dz);

    /* ================================================================
     * Step 1 — Solve PDE in the cortical ribbon
     * ================================================================ */
    phi = (float *)calloc(nvox, sizeof(float));
    if (!phi)
    {
        fprintf(stderr, "surf_pde_pial_white: alloc error (phi)\n");
        return -2;
    }

    if (use_ade)
        solve_ade_ribbon(labels, dims, ribbon_pial, ribbon_white,
                         phi, 2000, 1e-5f, verbose);
    else
        solve_laplace_ribbon(labels, dims, ribbon_pial, ribbon_white,
                             phi, 2000, 1e-5f, verbose);

    if (verbose)
        fprintf(stdout, "Phi stop thresholds: pial >= %.3f  white <= %.3f\n",
                phi_stop_pial, phi_stop_white);

    /* ================================================================
     * Step 2 — Compute gradient of φ (in voxel-grid coordinates)
     * ================================================================ */
    grad_x = (float *)malloc(sizeof(float) * nvox);
    grad_y = (float *)malloc(sizeof(float) * nvox);
    grad_z = (float *)malloc(sizeof(float) * nvox);
    if (!grad_x || !grad_y || !grad_z)
    {
        fprintf(stderr, "surf_pde_pial_white: alloc error (gradient)\n");
        free(phi);
        if (grad_x) free(grad_x);
        if (grad_y) free(grad_y);
        if (grad_z) free(grad_z);
        return -3;
    }

    gradient3D(phi, NULL, grad_x, grad_y, grad_z, dims, vx);

    /* ================================================================
     * Step 3 — Prepare output surfaces (copy of central)
     * ================================================================ */
    copy_polygons(central, pial_out);
    copy_polygons(central, white_out);
    compute_polygon_normals(central);

    /* ================================================================
     * Step 4 — Trace streamlines for every vertex
     * ================================================================ */
    int n_pial_ok = 0;
    int n_white_ok = 0;

    for (v = 0; v < n_points; v++)
    {
        double cx, cy, cz;       /* central vertex position   */
        double nx, ny, nz, nlen; /* surface normal (fallback) */
        double pos[3];
        double gx, gy, gz, glen;
        double max_travel, dist;
        int found;
        float phi_val;

        cx = Point_x(central->points[v]);
        cy = Point_y(central->points[v]);
        cz = Point_z(central->points[v]);

        /* Surface normal for fallback when gradient is too weak */
        nx = Point_x(central->normals[v]);
        ny = Point_y(central->normals[v]);
        nz = Point_z(central->normals[v]);
        nlen = sqrt(nx * nx + ny * ny + nz * nz);
        if (nlen > 1e-10)
        {
            nx /= nlen;
            ny /= nlen;
            nz /= nlen;
        }

        /* Limit travel distance using cortical thickness */
        max_travel = thickness_values
                         ? fabs(thickness_values[v]) * 1.0
                         : 4.0;
        if (max_travel < 0.5)  max_travel = 0.5;
        if (max_travel > 6.0)  max_travel = 6.0;

        /* ---- trace toward PIAL (follow +grad phi, toward phi = 1) ---- */
        pos[0] = cx;
        pos[1] = cy;
        pos[2] = cz;
        dist = 0.0;
        found = 0;

        for (step = 0; step < max_steps && dist < max_travel; step++)
        {
            gx = isoval(grad_x, pos[0], pos[1], pos[2], dims, nii_ptr);
            gy = isoval(grad_y, pos[0], pos[1], pos[2], dims, nii_ptr);
            gz = isoval(grad_z, pos[0], pos[1], pos[2], dims, nii_ptr);

            glen = sqrt(gx * gx + gy * gy + gz * gz);
            if (glen < min_grad)
            {
                /* Fallback: follow outward normal */
                gx = nx;
                gy = ny;
                gz = nz;
            }
            else
            {
                gx /= glen;
                gy /= glen;
                gz /= glen;
            }

            pos[0] += step_size * gx;
            pos[1] += step_size * gy;
            pos[2] += step_size * gz;
            dist += step_size;

            /* Stop when phi reaches the target pial isovalue */
            phi_val = isoval(phi, pos[0], pos[1], pos[2], dims, nii_ptr);
            if (phi_val >= phi_stop_pial)
            {
                found = 1;
                break;
            }
        }

        if (found)
        {
            fill_Point(pial_out->points[v], pos[0], pos[1], pos[2]);
            n_pial_ok++;
        }
        else
        {
            /* Fallback: half-thickness along outward normal */
            double t = thickness_values ? 0.5 * thickness_values[v] : 1.5;
            fill_Point(pial_out->points[v],
                       cx + t * nx, cy + t * ny, cz + t * nz);
        }

        /* ---- trace toward WHITE (follow -grad phi, toward phi = 0) ---- */
        pos[0] = cx;
        pos[1] = cy;
        pos[2] = cz;
        dist = 0.0;
        found = 0;

        for (step = 0; step < max_steps && dist < max_travel; step++)
        {
            gx = isoval(grad_x, pos[0], pos[1], pos[2], dims, nii_ptr);
            gy = isoval(grad_y, pos[0], pos[1], pos[2], dims, nii_ptr);
            gz = isoval(grad_z, pos[0], pos[1], pos[2], dims, nii_ptr);

            glen = sqrt(gx * gx + gy * gy + gz * gz);
            if (glen < min_grad)
            {
                /* Fallback: follow inward normal */
                gx = -nx;
                gy = -ny;
                gz = -nz;
            }
            else
            {
                gx /= glen;
                gy /= glen;
                gz /= glen;
            }

            /* Negative gradient direction (toward phi = 0, WM) */
            pos[0] -= step_size * gx;
            pos[1] -= step_size * gy;
            pos[2] -= step_size * gz;
            dist += step_size;

            /* Stop when phi reaches the target white isovalue */
            phi_val = isoval(phi, pos[0], pos[1], pos[2], dims, nii_ptr);
            if (phi_val <= phi_stop_white)
            {
                found = 1;
                break;
            }
        }

        if (found)
        {
            fill_Point(white_out->points[v], pos[0], pos[1], pos[2]);
            n_white_ok++;
        }
        else
        {
            double t = thickness_values ? 0.5 * thickness_values[v] : 1.5;
            fill_Point(white_out->points[v],
                       cx - t * nx, cy - t * ny, cz - t * nz);
        }
    }

    if (verbose)
        fprintf(stdout, "%s streamlines: pial %d/%d  white %d/%d converged\n",
                use_ade ? "ADE" : "Laplacian",
                n_pial_ok, n_points, n_white_ok, n_points);

    /* ================================================================
     * Step 5 — Light Laplacian smoothing for regularisation
     * ================================================================ */
    smooth_laplacian(pial_out, 5, 0.1, 0.5);
    smooth_laplacian(white_out, 5, 0.1, 0.5);

    /* Recompute normals for output */
    compute_polygon_normals(pial_out);
    compute_polygon_normals(white_out);

    /* ================================================================
     * Cleanup
     * ================================================================ */
    free(phi);
    free(grad_x);
    free(grad_y);
    free(grad_z);

    return 0;
}

/* -------------------------------------------------------------------
 * Public wrappers
 * -------------------------------------------------------------------*/

int surf_laplacian_pial_white(polygons_struct *central,
                              float *labels,
                              nifti_image *nii_ptr,
                              float lim_pial,
                              float lim_white,
                              const double *thickness_values,
                              polygons_struct *pial_out,
                              polygons_struct *white_out,
                              int verbose)
{
    return surf_pde_pial_white_internal(central, labels, nii_ptr,
                                        lim_pial, lim_white,
                                        thickness_values,
                                        pial_out, white_out,
                                        0, verbose);
}

int surf_ade_pial_white(polygons_struct *central,
                        float *labels,
                        nifti_image *nii_ptr,
                        float lim_pial,
                        float lim_white,
                        const double *thickness_values,
                        polygons_struct *pial_out,
                        polygons_struct *white_out,
                        int verbose)
{
    return surf_pde_pial_white_internal(central, labels, nii_ptr,
                                        lim_pial, lim_white,
                                        thickness_values,
                                        pial_out, white_out,
                                        1, verbose);
}
