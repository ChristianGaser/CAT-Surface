/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 *
 * Profile-based placement of the pial surface.
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <float.h>
#include <math.h>

#include "CAT_SurfPialProfile.h"

/* Facing-wall search: neighbourhood radius, lateral tolerance and the
 * normal opposition (n_i . n_j < -CONTACT_FACING) that separates the other
 * bank of a sulcus from the same sheet.  A facing wall less than
 * CONTACT_CROSSED behind a vertex has been crossed already: pial walls on the
 * two sides of a gyral blade are several millimetres apart. */
#define CONTACT_RADIUS 2.0
#define CONTACT_LATERAL 1.0
#define CONTACT_FACING 0.3
#define CONTACT_CROSSED 1.0

#define PROFILE_MAX_SAMPLES 256

#define TARGET_NONE 0
#define TARGET_CROSSING 1
#define TARGET_VALLEY 2

/* Trilinear sampler in world coordinates */
typedef struct
{
    const float *vol;
    int dims[3];
    double A[3][4]; /* world -> voxel */
} profile_sampler;

/* Uniform hash grid over the mesh vertices */
typedef struct
{
    int n[3];
    double min[3];
    double cell;
    int *head;
    int *next;
} vertex_grid;

/**
 * \brief Initialize pial profile placement options with defaults.
 *
 * \param opts (out) options structure to initialize
 * \return void
 */
void CAT_PialProfileOptionsInit(CAT_PialProfileOptions *opts)
{
    if (!opts)
        return;
    opts->isovalue = 1.5;
    opts->search_out = 2.0;
    opts->search_in = 1.0;
    opts->sample_step = 0.1;
    opts->valley_depth = 0.05;
    opts->step_fraction = 0.5;
    opts->max_step = 0.25;
    opts->max_offset = 2.0;
    opts->smooth_lambda = 0.5;
    opts->smooth_passes_start = 5;
    opts->smooth_passes_end = 1;
    opts->tangential_weight = 0.3;
    opts->concave_weight = 0.1;
    opts->contact_margin = 0.1;
    opts->fold_angle = 72.5;
    opts->iterations = 60;
    opts->verbose = 0;
}

static void
sampler_init(profile_sampler *S, const float *vol, nifti_image *nii)
{
    int r, c;
    mat44 inv = nifti_mat44_inverse(nii->sto_xyz);

    S->vol = vol;
    S->dims[0] = nii->nx;
    S->dims[1] = nii->ny;
    S->dims[2] = nii->nz;
    for (r = 0; r < 3; r++)
        for (c = 0; c < 4; c++)
            S->A[r][c] = inv.m[r][c];
}

static double
clamp_index(double u, int n)
{
    if (u < 0.0)
        return 0.0;
    if (u > n - 1.001)
        return n - 1.001;
    return u;
}

static double
sample_world(const profile_sampler *S, const double *x)
{
    int nx = S->dims[0], ny = S->dims[1], nz = S->dims[2];
    double u, v, w, fu, fv, fw, c00, c10, c01, c11;
    const float *V = S->vol;
    int i, j, k;
    size_t o, sx, sxy;

    u = S->A[0][0] * x[0] + S->A[0][1] * x[1] + S->A[0][2] * x[2] + S->A[0][3];
    v = S->A[1][0] * x[0] + S->A[1][1] * x[1] + S->A[1][2] * x[2] + S->A[1][3];
    w = S->A[2][0] * x[0] + S->A[2][1] * x[1] + S->A[2][2] * x[2] + S->A[2][3];
    u = clamp_index(u, nx);
    v = clamp_index(v, ny);
    w = clamp_index(w, nz);

    i = (int)u;
    j = (int)v;
    k = (int)w;
    fu = u - i;
    fv = v - j;
    fw = w - k;

    sx = (size_t)nx;
    sxy = (size_t)nx * (size_t)ny;
    o = (size_t)i + sx * (size_t)j + sxy * (size_t)k;

    c00 = V[o] * (1 - fu) + V[o + 1] * fu;
    c10 = V[o + sx] * (1 - fu) + V[o + sx + 1] * fu;
    c01 = V[o + sxy] * (1 - fu) + V[o + sxy + 1] * fu;
    c11 = V[o + sxy + sx] * (1 - fu) + V[o + sxy + sx + 1] * fu;

    return (c00 * (1 - fv) + c10 * fv) * (1 - fw) + (c01 * (1 - fv) + c11 * fv) * fw;
}

/*
 * Target offset along the unit normal n at p.  Inside the boundary the
 * profile is walked outwards and stops at the isovalue crossing or at the
 * bottom of a valley (a rise of valley_depth above the running minimum).
 * If the profile rises right away the vertex is past the valley bottom and
 * the minimum is searched inwards.  Outside the boundary the profile is
 * walked inwards to the crossing.  Without a boundary in range the vertex
 * is held (offset 0) inside and pulled back by search_in outside.
 */
static double
profile_target(const profile_sampler *S, const double *p, const double *n,
               const CAT_PialProfileOptions *o, int *kind)
{
    double I[2 * PROFILE_MAX_SAMPLES + 1], x[3], m;
    int n_in = (int)floor(o->search_in / o->sample_step + 0.5);
    int n_out = (int)floor(o->search_out / o->sample_step + 0.5);
    int q, c, q_min, k;
    const double target = o->isovalue;

    if (n_in > PROFILE_MAX_SAMPLES)
        n_in = PROFILE_MAX_SAMPLES;
    if (n_out > PROFILE_MAX_SAMPLES)
        n_out = PROFILE_MAX_SAMPLES;

    c = n_in; /* index of the vertex itself */
    for (q = -n_in; q <= n_out; q++)
    {
        for (k = 0; k < 3; k++)
            x[k] = p[k] + q * o->sample_step * n[k];
        I[c + q] = sample_world(S, x);
    }

    *kind = TARGET_NONE;

    if (I[c] < target)
    {
        for (q = 1; q <= n_in; q++)
        {
            if (I[c - q] >= target)
            {
                double a = I[c - q + 1], b = I[c - q];
                double fr = (a != b) ? (target - a) / (b - a) : 0.0;
                *kind = TARGET_CROSSING;
                return -(q - 1 + fr) * o->sample_step;
            }
        }
        return -o->search_in;
    }

    m = I[c];
    q_min = 0;
    for (q = 1; q <= n_out; q++)
    {
        if (I[c + q] <= target)
        {
            double a = I[c + q - 1], b = I[c + q];
            double fr = (a != b) ? (a - target) / (a - b) : 0.0;
            *kind = TARGET_CROSSING;
            return (q - 1 + fr) * o->sample_step;
        }
        if (I[c + q] < m)
        {
            m = I[c + q];
            q_min = q;
        }
        else if (I[c + q] > m + o->valley_depth)
        {
            *kind = TARGET_VALLEY;
            break;
        }
    }

    if (*kind != TARGET_VALLEY)
        return 0.0;

    if (q_min > 0)
        return q_min * o->sample_step;

    /* rising right away: find the valley bottom inwards */
    m = I[c];
    for (q = 1; q <= n_in; q++)
    {
        if (I[c - q] < m)
        {
            m = I[c - q];
            q_min = q;
        }
        else if (I[c - q] > m + o->valley_depth)
            break;
    }
    return -q_min * o->sample_step;
}

static int
grid_build(vertex_grid *G, const polygons_struct *P, double cell)
{
    double max[3];
    size_t n_cells, id;
    int i, k, idx[3];

    for (k = 0; k < 3; k++)
    {
        G->min[k] = DBL_MAX;
        max[k] = -DBL_MAX;
    }
    for (i = 0; i < P->n_points; i++)
        for (k = 0; k < 3; k++)
        {
            double x = Point_coord(P->points[i], k);
            if (x < G->min[k])
                G->min[k] = x;
            if (x > max[k])
                max[k] = x;
        }

    G->cell = cell;
    for (k = 0; k < 3; k++)
        G->n[k] = (int)((max[k] - G->min[k]) / cell) + 1;

    n_cells = (size_t)G->n[0] * (size_t)G->n[1] * (size_t)G->n[2];
    G->head = (int *)malloc(sizeof(int) * n_cells);
    G->next = (int *)malloc(sizeof(int) * P->n_points);
    if (!G->head || !G->next)
    {
        free(G->head);
        free(G->next);
        return -1;
    }

    for (id = 0; id < n_cells; id++)
        G->head[id] = -1;
    for (i = 0; i < P->n_points; i++)
    {
        for (k = 0; k < 3; k++)
            idx[k] = (int)((Point_coord(P->points[i], k) - G->min[k]) / cell);
        id = (size_t)idx[0] + (size_t)G->n[0] * ((size_t)idx[1] + (size_t)G->n[1] * (size_t)idx[2]);
        G->next[i] = G->head[id];
        G->head[id] = i;
    }
    return 0;
}

/* Signed distance along n_i to the nearest facing wall (DBL_MAX if none) */
static double
facing_gap(const polygons_struct *P, const vertex_grid *G, int i)
{
    double p[3], n[3], best = DBL_MAX;
    int idx[3], k, da, db, dc;

    for (k = 0; k < 3; k++)
    {
        p[k] = Point_coord(P->points[i], k);
        n[k] = Point_coord(P->normals[i], k);
        idx[k] = (int)((p[k] - G->min[k]) / G->cell);
    }

    for (da = -1; da <= 1; da++)
        for (db = -1; db <= 1; db++)
            for (dc = -1; dc <= 1; dc++)
            {
                int a = idx[0] + da, b = idx[1] + db, c = idx[2] + dc, j;
                if (a < 0 || b < 0 || c < 0 || a >= G->n[0] || b >= G->n[1] || c >= G->n[2])
                    continue;
                j = G->head[(size_t)a + (size_t)G->n[0] * ((size_t)b + (size_t)G->n[1] * (size_t)c)];
                for (; j >= 0; j = G->next[j])
                {
                    double d[3], g, lateral2, n_dot = 0.0;
                    if (j == i)
                        continue;
                    for (k = 0; k < 3; k++)
                        n_dot += n[k] * Point_coord(P->normals[j], k);
                    if (n_dot > -CONTACT_FACING)
                        continue;
                    for (k = 0; k < 3; k++)
                        d[k] = Point_coord(P->points[j], k) - p[k];
                    g = d[0] * n[0] + d[1] * n[1] + d[2] * n[2];
                    lateral2 = d[0] * d[0] + d[1] * d[1] + d[2] * d[2] - g * g;
                    if (lateral2 > CONTACT_LATERAL * CONTACT_LATERAL || g < -CONTACT_CROSSED)
                        continue;
                    if (g < best)
                        best = g;
                }
            }
    return best;
}

static void
face_normals(const polygons_struct *P, double (*fn)[3])
{
    int f, k;

    for (f = 0; f < P->n_items; f++)
    {
        int a = P->indices[3 * f], b = P->indices[3 * f + 1], c = P->indices[3 * f + 2];
        double u[3], v[3], len;
        for (k = 0; k < 3; k++)
        {
            u[k] = Point_coord(P->points[b], k) - Point_coord(P->points[a], k);
            v[k] = Point_coord(P->points[c], k) - Point_coord(P->points[a], k);
        }
        fn[f][0] = u[1] * v[2] - u[2] * v[1];
        fn[f][1] = u[2] * v[0] - u[0] * v[2];
        fn[f][2] = u[0] * v[1] - u[1] * v[0];
        len = sqrt(fn[f][0] * fn[f][0] + fn[f][1] * fn[f][1] + fn[f][2] * fn[f][2]);
        if (len > 0.0)
            for (k = 0; k < 3; k++)
                fn[f][k] /= len;
    }
}

static void
unit_vertex_normals(polygons_struct *P)
{
    int i;

    compute_polygon_normals(P);
    for (i = 0; i < P->n_points; i++)
    {
        double len = sqrt(Point_x(P->normals[i]) * Point_x(P->normals[i]) +
                          Point_y(P->normals[i]) * Point_y(P->normals[i]) +
                          Point_z(P->normals[i]) * Point_z(P->normals[i]));
        if (len > 1e-12)
            SCALE_VECTOR(P->normals[i], P->normals[i], 1.0 / len);
    }
}

/* Smooth the scalar update over the mesh (no damping, so no shrinkage) */
static void
smooth_update(double *d, double *tmp, int n_points, int *n_neighbours,
              int **neighbours, int passes, double lambda)
{
    int pass, i, j;

    for (pass = 0; pass < passes; pass++)
    {
        for (i = 0; i < n_points; i++)
        {
            double acc = 0.0;
            for (j = 0; j < n_neighbours[i]; j++)
                acc += d[neighbours[i][j]];
            tmp[i] = (1.0 - lambda) * d[i] + lambda * acc / n_neighbours[i];
        }
        memcpy(d, tmp, sizeof(double) * n_points);
    }
}

/* Revert the vertices of faces that rotated by more than acos(min_cos) */
static int
revert_folds(polygons_struct *P, double (*fn_prev)[3], double (*fn_cur)[3],
             double (*prev)[3], double min_cos)
{
    int pass, f, q, k, n_reverted = 0;

    for (pass = 0; pass < 3; pass++)
    {
        int changed = 0;
        face_normals(P, fn_cur);
        for (f = 0; f < P->n_items; f++)
        {
            double dt = fn_cur[f][0] * fn_prev[f][0] + fn_cur[f][1] * fn_prev[f][1] +
                        fn_cur[f][2] * fn_prev[f][2];
            if (dt >= min_cos)
                continue;
            for (q = 0; q < 3; q++)
            {
                int v = P->indices[3 * f + q];
                if (Point_x(P->points[v]) == prev[v][0] && Point_y(P->points[v]) == prev[v][1] &&
                    Point_z(P->points[v]) == prev[v][2])
                    continue;
                for (k = 0; k < 3; k++)
                    Point_coord(P->points[v], k) = prev[v][k];
                changed++;
            }
        }
        n_reverted += changed;
        if (!changed)
            break;
    }
    return n_reverted;
}

/**
 * \brief Move a pial surface onto the CSF/GM boundary by profile search.
 *
 * Each iteration: (1) search the profile along every vertex normal for the
 * isovalue crossing or a valley bottom, (2) take step_fraction of that offset,
 * limited to max_step, (3) smooth this update over the mesh, (4) add a normal
 * Laplacian term in concave regions only (a convex crown would be pulled
 * inwards), (5) clamp the outward step to half the gap to a facing sulcal
 * wall, (6) move along the normal with tangential relaxation, (7) limit the
 * offset from the start surface to max_offset and (8) revert vertices of
 * faces that rotated by more than fold_angle.
 *
 * \param pial (in/out) pial surface; start positions in, placed surface out
 * \param vol  (in)     volume data (e.g. label map with CSF=1, GM=2, WM=3)
 * \param nii  (in)     NIfTI header of vol
 * \param opts (in)     options, see CAT_PialProfileOptionsInit()
 * \return 0 on success, -1 on invalid arguments, -2 on allocation failure
 */
int CAT_SurfDeformPialProfile(polygons_struct *pial, const float *vol,
                              nifti_image *nii,
                              const CAT_PialProfileOptions *opts)
{
    int n_points, it, i, j, k, rc = 0;
    int *n_neighbours = NULL, **neighbours = NULL;
    double *d, *tmp, *gap, (*start)[3], (*start_n)[3], (*prev)[3];
    double (*fn_prev)[3], (*fn_cur)[3];
    double min_cos;
    profile_sampler S;

    if (!pial || !vol || !nii || !opts || pial->n_points < 4 ||
        opts->sample_step <= 0.0 || opts->iterations < 0)
        return -1;
    for (i = 0; i < pial->n_items; i++)
        if (pial->end_indices[i] != 3 * (i + 1))
            return -1; /* triangle meshes only */

    n_points = pial->n_points;
    min_cos = cos(opts->fold_angle * PI / 180.0);

    d = (double *)malloc(sizeof(double) * n_points);
    tmp = (double *)malloc(sizeof(double) * n_points);
    gap = (double *)malloc(sizeof(double) * n_points);
    start = malloc(sizeof(double[3]) * n_points);
    start_n = malloc(sizeof(double[3]) * n_points);
    prev = malloc(sizeof(double[3]) * n_points);
    fn_prev = malloc(sizeof(double[3]) * pial->n_items);
    fn_cur = malloc(sizeof(double[3]) * pial->n_items);
    if (!d || !tmp || !gap || !start || !start_n || !prev || !fn_prev || !fn_cur)
    {
        rc = -2;
        goto cleanup;
    }

    sampler_init(&S, vol, nii);
    create_polygon_point_neighbours(pial, TRUE, &n_neighbours, &neighbours, NULL, NULL);

    unit_vertex_normals(pial);
    for (i = 0; i < n_points; i++)
        for (k = 0; k < 3; k++)
        {
            start[i][k] = Point_coord(pial->points[i], k);
            start_n[i][k] = Point_coord(pial->normals[i], k);
        }

    for (it = 0; it < opts->iterations; it++)
    {
        double fr = (opts->iterations > 1) ? (double)it / (opts->iterations - 1) : 1.0;
        int passes = (int)floor(opts->smooth_passes_start +
                                (opts->smooth_passes_end - opts->smooth_passes_start) * fr + 0.5);
        int n_kind[3] = {0, 0, 0}, n_clamped = 0, n_reverted;
        double sum_abs = 0.0;
        vertex_grid G;

        if (it > 0)
            unit_vertex_normals(pial);

        /* (1)-(2) per-vertex target offsets */
        for (i = 0; i < n_points; i++)
        {
            double p[3], n[3], s;
            int kind;
            for (k = 0; k < 3; k++)
            {
                p[k] = Point_coord(pial->points[i], k);
                n[k] = Point_coord(pial->normals[i], k);
            }
            s = opts->step_fraction * profile_target(&S, p, n, opts, &kind);
            n_kind[kind]++;
            d[i] = fmax(-opts->max_step, fmin(opts->max_step, s));
        }

        /* (3) regularize the update, not the accumulated displacement */
        smooth_update(d, tmp, n_points, n_neighbours, neighbours, passes,
                      opts->smooth_lambda);

        /* snapshot for the tangential term and the fold check */
        for (i = 0; i < n_points; i++)
            for (k = 0; k < 3; k++)
                prev[i][k] = Point_coord(pial->points[i], k);
        face_normals(pial, fn_prev);

        /* facing-wall gaps from the positions before this step */
        if (grid_build(&G, pial, CONTACT_RADIUS) != 0)
        {
            rc = -2;
            goto cleanup;
        }
        for (i = 0; i < n_points; i++)
            gap[i] = facing_gap(pial, &G, i);
        free(G.head);
        free(G.next);

        for (i = 0; i < n_points; i++)
        {
            double c[3] = {0.0, 0.0, 0.0}, t[3], n[3], x[3], dot = 0.0, dn, off = 0.0;

            for (k = 0; k < 3; k++)
                n[k] = Point_coord(pial->normals[i], k);
            for (j = 0; j < n_neighbours[i]; j++)
                for (k = 0; k < 3; k++)
                    c[k] += prev[neighbours[i][j]][k];
            for (k = 0; k < 3; k++)
            {
                t[k] = c[k] / n_neighbours[i] - prev[i][k];
                dot += t[k] * n[k];
            }
            for (k = 0; k < 3; k++)
                t[k] -= dot * n[k];

            /* (4) concave regions only: dot > 0 where the neighbours lie outwards */
            dn = d[i] + ((dot > 0.0) ? opts->concave_weight * dot : 0.0);

            /* (5) facing walls meet in the middle */
            if (gap[i] != DBL_MAX && dn > 0.5 * (gap[i] - opts->contact_margin))
            {
                dn = 0.5 * (gap[i] - opts->contact_margin);
                n_clamped++;
            }

            /* (6) move */
            for (k = 0; k < 3; k++)
                x[k] = prev[i][k] + dn * n[k] + opts->tangential_weight * t[k];

            /* (7) stay within max_offset of the start surface */
            for (k = 0; k < 3; k++)
                off += (x[k] - start[i][k]) * start_n[i][k];
            if (off > opts->max_offset)
                for (k = 0; k < 3; k++)
                    x[k] -= (off - opts->max_offset) * start_n[i][k];
            else if (off < -opts->max_offset)
                for (k = 0; k < 3; k++)
                    x[k] -= (off + opts->max_offset) * start_n[i][k];

            fill_Point(pial->points[i], x[0], x[1], x[2]);
            sum_abs += fabs(dn);
        }

        /* (8) folds */
        n_reverted = revert_folds(pial, fn_prev, fn_cur, prev, min_cos);

        if (opts->verbose)
            fprintf(stdout, "\rPial profile: iter %03d | mean |step| %.4f | crossing %d valley %d none %d"
                            " | clamped %d | reverted %d   ",
                    it + 1, sum_abs / n_points, n_kind[TARGET_CROSSING], n_kind[TARGET_VALLEY],
                    n_kind[TARGET_NONE], n_clamped, n_reverted);
    }
    if (opts->verbose && opts->iterations > 0)
        fprintf(stdout, "\n");

    compute_polygon_normals(pial);

cleanup:
    free(d);
    free(tmp);
    free(gap);
    free(start);
    free(start_n);
    free(prev);
    free(fn_prev);
    free(fn_cur);
    if (neighbours)
        delete_polygon_point_neighbours(pial, n_neighbours, neighbours, NULL, NULL);
    return rc;
}
