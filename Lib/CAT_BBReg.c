/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Boundary-Based Registration (BBR) — CAT-native implementation.
 *
 * Reference: Greve & Fischl (2009), NeuroImage 48:63-72.
 *   "Accurate and robust brain image alignment using boundary-based registration"
 *
 * Copyright Christian Gaser, University of Jena.
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>

#include <bicpl.h>

#include "CAT_BBReg.h"
#include "CAT_Vol.h"

/* =========================================================================
 * Internal helpers
 * ========================================================================= */

/* Trilinear sample at continuous voxel coords (vx, vy, vz).
 * Returns NaN if any corner of the interpolation cell is outside [0, dims-1].
 * This is used by BBR where out-of-FOV vertices must be excluded, not clamped. */
static double sample_vox(const float *vol, const int dims[3],
                         double vx, double vy, double vz)
{
    int x0 = (int)floor(vx), y0 = (int)floor(vy), z0 = (int)floor(vz);
    int x1 = x0 + 1,         y1 = y0 + 1,         z1 = z0 + 1;

    /* Strict out-of-bounds → exclude this vertex entirely */
    if (x0 < 0 || x1 >= dims[0] ||
        y0 < 0 || y1 >= dims[1] ||
        z0 < 0 || z1 >= dims[2])
        return (double)NAN;

    double wx = vx - x0, wy = vy - y0, wz = vz - z0;

#define IDX(a,b,c) ((a) + dims[0]*((b) + dims[1]*(c)))
    double v000 = vol[IDX(x0,y0,z0)]; double v100 = vol[IDX(x1,y0,z0)];
    double v010 = vol[IDX(x0,y1,z0)]; double v110 = vol[IDX(x1,y1,z0)];
    double v001 = vol[IDX(x0,y0,z1)]; double v101 = vol[IDX(x1,y0,z1)];
    double v011 = vol[IDX(x0,y1,z1)]; double v111 = vol[IDX(x1,y1,z1)];
#undef IDX

    return v000*(1-wx)*(1-wy)*(1-wz) + v100*wx*(1-wy)*(1-wz) +
           v010*(1-wx)*wy*(1-wz)     + v110*wx*wy*(1-wz)     +
           v001*(1-wx)*(1-wy)*wz     + v101*wx*(1-wy)*wz     +
           v011*(1-wx)*wy*wz         + v111*wx*wy*wz;
}

/* =========================================================================
 * Public API
 * ========================================================================= */

/**
 * \brief Convert 6-DOF rigid parameters to a 4×4 homogeneous matrix (row-major).
 *
 * The resulting matrix maps points in the *fixed* surface (T1) RAS space to
 * the *moving* volume (EPI) RAS space — the direction the BBR cost function
 * applies when transforming surface sample points before EPI interpolation.
 * Use CAT_BBReg_invert_matrix() to get the standard EPI→T1 output matrix.
 *
 *   M = T_trans * Rz * Ry * Rx
 *
 * \param p  (in)  rigid parameters (tx, ty, tz mm; rx, ry, rz radians)
 * \param m  (out) 4×4 matrix, row-major, 16 doubles
 */
void CAT_BBReg_params_to_matrix(const CAT_RigidParams *p, double m[16])
{
    double cx = cos(p->rx), sx = sin(p->rx);
    double cy = cos(p->ry), sy = sin(p->ry);
    double cz = cos(p->rz), sz = sin(p->rz);

    /* Rotation: R = Rz * Ry * Rx (column-major interpretation) */
    double r00 = cz * cy;
    double r01 = cz * sy * sx - sz * cx;
    double r02 = cz * sy * cx + sz * sx;

    double r10 = sz * cy;
    double r11 = sz * sy * sx + cz * cx;
    double r12 = sz * sy * cx - cz * sx;

    double r20 = -sy;
    double r21 = cy * sx;
    double r22 = cy * cx;

    m[0] = r00;
    m[1] = r01;
    m[2] = r02;
    m[3] = p->tx;
    m[4] = r10;
    m[5] = r11;
    m[6] = r12;
    m[7] = p->ty;
    m[8] = r20;
    m[9] = r21;
    m[10] = r22;
    m[11] = p->tz;
    m[12] = 0.0;
    m[13] = 0.0;
    m[14] = 0.0;
    m[15] = 1.0;
}

/**
 * \brief Apply a 4×4 RAS transform to every vertex of a surface (in-place).
 *
 * \param surface (in/out) surface mesh whose vertices are modified
 * \param m       (in)     4×4 transform, row-major
 */
void CAT_BBReg_apply_matrix(polygons_struct *surface, const double m[16])
{
    int i;
    double x, y, z, nx, ny, nz;

    for (i = 0; i < surface->n_points; i++)
    {
        x = Point_x(surface->points[i]);
        y = Point_y(surface->points[i]);
        z = Point_z(surface->points[i]);

        nx = m[0] * x + m[1] * y + m[2] * z + m[3];
        ny = m[4] * x + m[5] * y + m[6] * z + m[7];
        nz = m[8] * x + m[9] * y + m[10] * z + m[11];

        fill_Point(surface->points[i], nx, ny, nz);
    }
}

/**
 * \brief Invert a 4×4 rigid homogeneous matrix analytically.
 *
 * \param m    (in)  4×4 row-major rigid matrix
 * \param inv  (out) 4×4 row-major inverse
 */
void CAT_BBReg_invert_matrix(const double m[16], double inv[16])
{
    /* R^T */
    inv[0] = m[0];
    inv[1] = m[4];
    inv[2] = m[8];
    inv[4] = m[1];
    inv[5] = m[5];
    inv[6] = m[9];
    inv[8] = m[2];
    inv[9] = m[6];
    inv[10] = m[10];

    /* -R^T t */
    inv[3] = -(inv[0] * m[3] + inv[1] * m[7] + inv[2] * m[11]);
    inv[7] = -(inv[4] * m[3] + inv[5] * m[7] + inv[6] * m[11]);
    inv[11] = -(inv[8] * m[3] + inv[9] * m[7] + inv[10] * m[11]);

    inv[12] = 0.0;
    inv[13] = 0.0;
    inv[14] = 0.0;
    inv[15] = 1.0;
}

/**
 * \brief Compute the BBR cost over one or more surfaces.
 *
 * For each non-masked vertex the volume is sampled at WM and GM sides of the
 * surface normal.  When thickness data and a fractional projection factor are
 * provided the per-vertex GM offset is frac * thickness[i]; otherwise the
 * fixed gm_dist is used.
 *
 * \param p               (in)  current 6-DOF rigid parameters
 * \param surfs           (in)  array of CAT_SurfData structures
 * \param n_surfs         (in)  number of elements in surfs[]
 * \param vol             (in)  floating-point volume data (moving image)
 * \param nii_ptr         (in)  NIfTI header of the moving volume
 * \param dims            (in)  volume dimensions [nx, ny, nz]
 * \param wm_dist         (in)  WM sampling offset (mm, positive)
 * \param gm_dist         (in)  absolute GM offset (mm); used when no thickness
 * \param slope           (in)  saturation slope of the cost function
 * \param invert_contrast (in)  non-zero for T2/BOLD
 * \return mean BBR cost (lower = better alignment)
 */
double
CAT_BBReg_cost(const CAT_RigidParams *p,
               const CAT_SurfData *surfs,
               int n_surfs,
               float *vol,
               nifti_image *nii_ptr,
               int dims[3],
               double wm_dist,
               double gm_dist,
               double slope,
               int invert_contrast)
{
    double m[16];
    double cost = 0.0;
    int n_valid = 0;
    int s, i;
    const double eps = 1e-6;

    CAT_BBReg_params_to_matrix(p, m);

    /* Precompute voxel-space inverse once (world → voxel) to avoid
     * recomputing it inside isoval() for every sample point. */
    mat44 sto_ijk = nifti_mat44_inverse(nii_ptr->sto_xyz);
    /* Row-major 3×4 of sto_ijk for fast dot-product */
    double ij[12];
    int r, c;
    for (r = 0; r < 3; r++)
        for (c = 0; c < 4; c++)
            ij[r*4+c] = (double)sto_ijk.m[r][c];

    for (s = 0; s < n_surfs; s++)
    {
        const polygons_struct *surface = surfs[s].surface;
        const float *mask = surfs[s].cortex_mask;
        const float *thick = surfs[s].thickness;
        double frac = surfs[s].gm_proj_frac;

        for (i = 0; i < surface->n_points; i++)
        {
            /* Skip vertices outside cortex label */
            if (mask && mask[i] <= 0.5f)
                continue;

            double vx = Point_x(surface->points[i]);
            double vy = Point_y(surface->points[i]);
            double vz = Point_z(surface->points[i]);

            double nx = Vector_x(surface->normals[i]);
            double ny = Vector_y(surface->normals[i]);
            double nz = Vector_z(surface->normals[i]);

            double nlen = sqrt(nx * nx + ny * ny + nz * nz);
            if (nlen < 1e-10)
                continue;
            nx /= nlen;
            ny /= nlen;
            nz /= nlen;

            /* GM offset: fractional thickness or fixed distance */
            double d_gm = gm_dist;
            if (thick && frac > 0.0)
                d_gm = frac * (double)thick[i];

            /* WM sample point (inward along normal) */
            double wm_x = vx - wm_dist * nx;
            double wm_y = vy - wm_dist * ny;
            double wm_z = vz - wm_dist * nz;

            /* GM sample point (outward along normal) */
            double gm_x = vx + d_gm * nx;
            double gm_y = vy + d_gm * ny;
            double gm_z = vz + d_gm * nz;

            /* Apply rigid transform: surface RAS → EPI RAS */
            double twm_x = m[0]*wm_x + m[1]*wm_y + m[2]*wm_z + m[3];
            double twm_y = m[4]*wm_x + m[5]*wm_y + m[6]*wm_z + m[7];
            double twm_z = m[8]*wm_x + m[9]*wm_y + m[10]*wm_z + m[11];

            double tgm_x = m[0]*gm_x + m[1]*gm_y + m[2]*gm_z + m[3];
            double tgm_y = m[4]*gm_x + m[5]*gm_y + m[6]*gm_z + m[7];
            double tgm_z = m[8]*gm_x + m[9]*gm_y + m[10]*gm_z + m[11];

            /* EPI RAS → voxel (using precomputed inverse) */
            double wvx = ij[0]*twm_x + ij[1]*twm_y + ij[2]*twm_z + ij[3];
            double wvy = ij[4]*twm_x + ij[5]*twm_y + ij[6]*twm_z + ij[7];
            double wvz = ij[8]*twm_x + ij[9]*twm_y + ij[10]*twm_z + ij[11];

            double gvx = ij[0]*tgm_x + ij[1]*tgm_y + ij[2]*tgm_z + ij[3];
            double gvy = ij[4]*tgm_x + ij[5]*tgm_y + ij[6]*tgm_z + ij[7];
            double gvz = ij[8]*tgm_x + ij[9]*tgm_y + ij[10]*tgm_z + ij[11];

            /* sample_vox returns NaN when either point is outside the FOV,
             * ensuring such vertices are excluded rather than clamped to
             * edge values (which would bias the cost). */
            double I_wm = sample_vox(vol, dims, wvx, wvy, wvz);
            double I_gm = sample_vox(vol, dims, gvx, gvy, gvz);

            if (isnan(I_wm) || isnan(I_gm))
                continue;

            double denom = 0.5 * (I_wm + I_gm) + eps;
            /* Percent contrast (×100): slope=0.5 saturates at ~6% contrast.
             * Matching neuroreg/bbregister convention. */
            double c = 100.0 * (I_wm - I_gm) / denom;
            if (invert_contrast)
                c = -c;

            double t = tanh(slope * c);
            cost += 1.0 - t * t;
            n_valid++;
        }
    }

    if (n_valid == 0)
        return 1.0;
    return cost / (double)n_valid;
}

/* Print BBR intensity statistics at a given parameter set (verbose diagnostic). */
static void bbr_print_stats(const CAT_RigidParams *p,
                            const CAT_SurfData *surfs, int n_surfs,
                            float *vol, nifti_image *nii_ptr, int dims[3],
                            double wm_dist, double gm_dist)
{
    double m[16];
    CAT_BBReg_params_to_matrix(p, m);
    mat44 sto_ijk = nifti_mat44_inverse(nii_ptr->sto_xyz);
    double ij[12];
    int r, col;
    for (r = 0; r < 3; r++)
        for (col = 0; col < 4; col++)
            ij[r*4+col] = (double)sto_ijk.m[r][col];

    double sum_wm = 0.0, sum_gm = 0.0;
    int n_valid = 0, n_fov = 0;
    int s, i;

    for (s = 0; s < n_surfs; s++)
    {
        const polygons_struct *surf = surfs[s].surface;
        const float *mask = surfs[s].cortex_mask;
        const float *thick = surfs[s].thickness;
        double frac = surfs[s].gm_proj_frac;

        for (i = 0; i < surf->n_points; i++)
        {
            if (mask && mask[i] <= 0.5f) continue;
            n_fov++;

            double vx = Point_x(surf->points[i]);
            double vy = Point_y(surf->points[i]);
            double vz = Point_z(surf->points[i]);
            double nx = Vector_x(surf->normals[i]);
            double ny = Vector_y(surf->normals[i]);
            double nz = Vector_z(surf->normals[i]);
            double nlen = sqrt(nx*nx + ny*ny + nz*nz);
            if (nlen < 1e-10) continue;
            nx /= nlen; ny /= nlen; nz /= nlen;

            double d_gm = gm_dist;
            if (thick && frac > 0.0) d_gm = frac * (double)thick[i];

            double wm_x = vx - wm_dist*nx, wm_y = vy - wm_dist*ny, wm_z = vz - wm_dist*nz;
            double gm_x = vx + d_gm*nx,   gm_y = vy + d_gm*ny,   gm_z = vz + d_gm*nz;

            double twm_x = m[0]*wm_x + m[1]*wm_y + m[2]*wm_z + m[3];
            double twm_y = m[4]*wm_x + m[5]*wm_y + m[6]*wm_z + m[7];
            double twm_z = m[8]*wm_x + m[9]*wm_y + m[10]*wm_z + m[11];
            double tgm_x = m[0]*gm_x + m[1]*gm_y + m[2]*gm_z + m[3];
            double tgm_y = m[4]*gm_x + m[5]*gm_y + m[6]*gm_z + m[7];
            double tgm_z = m[8]*gm_x + m[9]*gm_y + m[10]*gm_z + m[11];

            double wvx = ij[0]*twm_x+ij[1]*twm_y+ij[2]*twm_z+ij[3];
            double wvy = ij[4]*twm_x+ij[5]*twm_y+ij[6]*twm_z+ij[7];
            double wvz = ij[8]*twm_x+ij[9]*twm_y+ij[10]*twm_z+ij[11];
            double gvx = ij[0]*tgm_x+ij[1]*tgm_y+ij[2]*tgm_z+ij[3];
            double gvy = ij[4]*tgm_x+ij[5]*tgm_y+ij[6]*tgm_z+ij[7];
            double gvz = ij[8]*tgm_x+ij[9]*tgm_y+ij[10]*tgm_z+ij[11];

            double I_wm = sample_vox(vol, dims, wvx, wvy, wvz);
            double I_gm = sample_vox(vol, dims, gvx, gvy, gvz);
            if (isnan(I_wm) || isnan(I_gm)) continue;

            sum_wm += I_wm;
            sum_gm += I_gm;
            n_valid++;
        }
    }

    if (n_valid > 0)
        printf("  mean I_WM=%.1f  mean I_GM=%.1f  "
               "contrast=(WM-GM)/mean=%.1f%%  "
               "valid vertices=%d/%d\n",
               sum_wm / n_valid, sum_gm / n_valid,
               100.0 * (sum_wm - sum_gm) / (0.5*(sum_wm+sum_gm) + 1e-6),
               n_valid, n_fov);
}

/* =========================================================================
 * Contrast auto-detection
 * ========================================================================= */

/**
 * \brief Detect image contrast type by sampling WM and GM intensities.
 */
int CAT_BBReg_detect_contrast(const CAT_RigidParams *p,
                              const CAT_SurfData *surfs, int n_surfs,
                              float *vol, nifti_image *nii_ptr, int dims[3],
                              double wm_dist, double gm_dist,
                              int verbose)
{
    double m[16];
    CAT_BBReg_params_to_matrix(p, m);
    mat44 sto_ijk = nifti_mat44_inverse(nii_ptr->sto_xyz);
    double ij[12];
    int r, col;
    for (r = 0; r < 3; r++)
        for (col = 0; col < 4; col++)
            ij[r*4+col] = (double)sto_ijk.m[r][col];

    double sum_wm = 0.0, sum_gm = 0.0;
    int n_valid = 0;
    int s, i;

    for (s = 0; s < n_surfs; s++)
    {
        const polygons_struct *surf = surfs[s].surface;
        const float *mask  = surfs[s].cortex_mask;
        const float *thick = surfs[s].thickness;
        double frac = surfs[s].gm_proj_frac;

        for (i = 0; i < surf->n_points; i++)
        {
            if (mask && mask[i] <= 0.5f) continue;

            double vx = Point_x(surf->points[i]);
            double vy = Point_y(surf->points[i]);
            double vz = Point_z(surf->points[i]);
            double nx = Vector_x(surf->normals[i]);
            double ny = Vector_y(surf->normals[i]);
            double nz = Vector_z(surf->normals[i]);
            double nlen = sqrt(nx*nx + ny*ny + nz*nz);
            if (nlen < 1e-10) continue;
            nx /= nlen; ny /= nlen; nz /= nlen;

            double d_gm = (thick && frac > 0.0) ? frac * (double)thick[i]
                                                 : gm_dist;

            double wm_x = vx - wm_dist*nx, wm_y = vy - wm_dist*ny;
            double wm_z = vz - wm_dist*nz;
            double gm_x = vx + d_gm*nx,   gm_y = vy + d_gm*ny;
            double gm_z = vz + d_gm*nz;

            /* Apply rigid transform */
            double twm_x = m[0]*wm_x+m[1]*wm_y+m[2]*wm_z+m[3];
            double twm_y = m[4]*wm_x+m[5]*wm_y+m[6]*wm_z+m[7];
            double twm_z = m[8]*wm_x+m[9]*wm_y+m[10]*wm_z+m[11];
            double tgm_x = m[0]*gm_x+m[1]*gm_y+m[2]*gm_z+m[3];
            double tgm_y = m[4]*gm_x+m[5]*gm_y+m[6]*gm_z+m[7];
            double tgm_z = m[8]*gm_x+m[9]*gm_y+m[10]*gm_z+m[11];

            double wvx = ij[0]*twm_x+ij[1]*twm_y+ij[2]*twm_z+ij[3];
            double wvy = ij[4]*twm_x+ij[5]*twm_y+ij[6]*twm_z+ij[7];
            double wvz = ij[8]*twm_x+ij[9]*twm_y+ij[10]*twm_z+ij[11];
            double gvx = ij[0]*tgm_x+ij[1]*tgm_y+ij[2]*tgm_z+ij[3];
            double gvy = ij[4]*tgm_x+ij[5]*tgm_y+ij[6]*tgm_z+ij[7];
            double gvz = ij[8]*tgm_x+ij[9]*tgm_y+ij[10]*tgm_z+ij[11];

            double I_wm = sample_vox(vol, dims, wvx, wvy, wvz);
            double I_gm = sample_vox(vol, dims, gvx, gvy, gvz);
            if (isnan(I_wm) || isnan(I_gm)) continue;

            sum_wm += I_wm;
            sum_gm += I_gm;
            n_valid++;
        }
    }

    if (n_valid < 100)
    {
        if (verbose)
            fprintf(stderr, "Contrast detection: only %d valid vertices"
                    " — defaulting to T1.\n", n_valid);
        return -1;
    }

    double mean_wm  = sum_wm / n_valid;
    double mean_gm  = sum_gm / n_valid;
    double mean_all = 0.5 * (mean_wm + mean_gm);
    double contrast = (mean_all > 1e-6) ?
                      (mean_wm - mean_gm) / mean_all : 0.0;

    /* Require at least 1 % WM/GM contrast to make a reliable call */
    if (fabs(contrast) < 0.01)
    {
        if (verbose)
            fprintf(stderr, "Contrast detection: WM/GM contrast too low"
                    " (%.2f%%) — defaulting to T1.\n", 100.0 * contrast);
        return -1;
    }

    int is_t2 = (contrast < 0.0);  /* WM darker than GM → T2/BOLD */

    if (verbose)
        printf("Contrast detection: mean I_WM=%.1f  mean I_GM=%.1f"
               "  contrast=%.1f%%  → %s\n",
               mean_wm, mean_gm, 100.0 * contrast,
               is_t2 ? "T2/BOLD (inverted)" : "T1/FLAIR");

    return is_t2;
}

/* =========================================================================
 * Powell optimizer (self-contained, no external deps)
 * ========================================================================= */

/* Context struct passed to the 1-D line function. */
typedef struct
{
    CAT_RigidParams base;
    double dir[6]; /* unit direction in parameter space */
    const CAT_SurfData *surfs;
    int n_surfs;
    float *vol;
    nifti_image *nii_ptr;
    int dims[3];
    double wm_dist;
    double gm_dist;
    double slope;
    int invert_contrast;
} LineCtx;

/* Evaluate cost along a 1-D line: params = base + t * dir */
static double line_cost(const LineCtx *ctx, double t)
{
    CAT_RigidParams p;
    p.tx = ctx->base.tx + t * ctx->dir[0];
    p.ty = ctx->base.ty + t * ctx->dir[1];
    p.tz = ctx->base.tz + t * ctx->dir[2];
    p.rx = ctx->base.rx + t * ctx->dir[3];
    p.ry = ctx->base.ry + t * ctx->dir[4];
    p.rz = ctx->base.rz + t * ctx->dir[5];

    return CAT_BBReg_cost(&p, ctx->surfs, ctx->n_surfs,
                          ctx->vol, ctx->nii_ptr,
                          (int *)ctx->dims, ctx->wm_dist, ctx->gm_dist,
                          ctx->slope, ctx->invert_contrast);
}

/* Brent's 1-D minimiser (no derivative).
 * Returns t_min in [ax, cx] and sets *fmin. */
#define BRENT_ITMAX 100
#define BRENT_CGOLD 0.3819660
#define BRENT_ZEPS 1.0e-10

static double brent_minimize(const LineCtx *ctx,
                             double ax, double bx, double cx,
                             double tol, double *fmin)
{
    int iter;
    double a, b, d = 0.0, e = 0.0, etemp, fu, fv, fw, fx;
    double p, q, r, tol1, tol2, u, v, w, x, xm;

    a = (ax < cx) ? ax : cx;
    b = (ax > cx) ? ax : cx;
    x = w = v = bx;
    fw = fv = fx = line_cost(ctx, x);

    for (iter = 0; iter < BRENT_ITMAX; iter++)
    {
        xm = 0.5 * (a + b);
        tol1 = tol * fabs(x) + BRENT_ZEPS;
        tol2 = 2.0 * tol1;

        if (fabs(x - xm) <= (tol2 - 0.5 * (b - a)))
        {
            *fmin = fx;
            return x;
        }

        if (fabs(e) > tol1)
        {
            r = (x - w) * (fx - fv);
            q = (x - v) * (fx - fw);
            p = (x - v) * q - (x - w) * r;
            q = 2.0 * (q - r);
            if (q > 0.0)
                p = -p;
            q = fabs(q);
            etemp = e;
            e = d;

            if (fabs(p) >= fabs(0.5 * q * etemp) ||
                p <= q * (a - x) || p >= q * (b - x))
            {
                e = (x >= xm) ? a - x : b - x;
                d = BRENT_CGOLD * e;
            }
            else
            {
                d = p / q;
                u = x + d;
                if (u - a < tol2 || b - u < tol2)
                    d = (xm - x >= 0.0) ? tol1 : -tol1;
            }
        }
        else
        {
            e = (x >= xm) ? a - x : b - x;
            d = BRENT_CGOLD * e;
        }

        u = (fabs(d) >= tol1) ? x + d : x + ((d >= 0.0) ? tol1 : -tol1);
        fu = line_cost(ctx, u);

        if (fu <= fx)
        {
            if (u < x)
                b = x;
            else
                a = x;
            v = w;
            fv = fw;
            w = x;
            fw = fx;
            x = u;
            fx = fu;
        }
        else
        {
            if (u < x)
                a = u;
            else
                b = u;
            if (fu <= fw || w == x)
            {
                v = w;
                fv = fw;
                w = u;
                fw = fu;
            }
            else if (fu <= fv || v == x || v == w)
            {
                v = u;
                fv = fu;
            }
        }
    }
    *fmin = fx;
    return x;
}

/* Bracket a minimum along the line and call Brent. */
static double line_minimize(const LineCtx *ctx,
                            double step, double *t_out)
{
    double t0 = 0.0, t1 = step;
    double f0 = line_cost(ctx, t0);
    double f1 = line_cost(ctx, t1);
    double fmin;

    /* Make sure we step downhill */
    if (f1 > f0)
    {
        t1 = -step;
        f1 = line_cost(ctx, t1);
    }
    if (f1 > f0)
    {
        *t_out = 0.0;
        return f0;
    } /* no descent direction */

    /* Expand bracket until we find a minimum */
    double t2 = t1 * 1.618;
    double f2 = line_cost(ctx, t2);
    int bracket_iter = 0;
    while (f2 < f1 && bracket_iter++ < 50)
    {
        t0 = t1;
        f0 = f1;
        t1 = t2;
        f1 = f2;
        t2 = t1 + 1.618 * (t1 - t0);
        f2 = line_cost(ctx, t2);
    }

    /* Brent on [t0, t2] with initial guess t1 */
    double tmin = brent_minimize(ctx, t0, t1, t2, 1e-4, &fmin);
    *t_out = tmin;
    return fmin;
}

/* Powell's direction-set method.  Mutates *p toward the minimum. */
static double powell_minimise(CAT_RigidParams *p,
                              const CAT_SurfData *surfs,
                              int n_surfs,
                              float *vol, nifti_image *nii_ptr,
                              int dims[3],
                              double wm_dist, double gm_dist,
                              double slope, int invert_contrast,
                              int max_iter, double tol, int verbose)
{
    const int N = 6;                                        /* number of DOF */
    double steps[6] = {0.5, 0.5, 0.5, 0.005, 0.005, 0.005}; /* mm / rad */
    /* Direction matrix: start as identity */
    double dirs[6][6];
    int i, j, iter;
    double f_curr, f_prev, t_opt, f_opt;
    double delta_max, delta_i;
    int idx_max;

    for (i = 0; i < N; i++)
        for (j = 0; j < N; j++)
            dirs[i][j] = (i == j) ? 1.0 : 0.0;

    LineCtx ctx;
    ctx.surfs = surfs;
    ctx.n_surfs = n_surfs;
    ctx.vol = vol;
    ctx.nii_ptr = nii_ptr;
    ctx.dims[0] = dims[0];
    ctx.dims[1] = dims[1];
    ctx.dims[2] = dims[2];
    ctx.wm_dist = wm_dist;
    ctx.gm_dist = gm_dist;
    ctx.slope = slope;
    ctx.invert_contrast = invert_contrast;

    f_curr = CAT_BBReg_cost(p, surfs, n_surfs, vol, nii_ptr, dims,
                            wm_dist, gm_dist, slope, invert_contrast);

    for (iter = 0; iter < max_iter; iter++)
    {
        f_prev = f_curr;
        delta_max = 0.0;
        idx_max = 0;

        double p_arr[6];
        p_arr[0] = p->tx;
        p_arr[1] = p->ty;
        p_arr[2] = p->tz;
        p_arr[3] = p->rx;
        p_arr[4] = p->ry;
        p_arr[5] = p->rz;

        for (i = 0; i < N; i++)
        {
            /* Set up line context */
            ctx.base.tx = p->tx;
            ctx.base.ty = p->ty;
            ctx.base.tz = p->tz;
            ctx.base.rx = p->rx;
            ctx.base.ry = p->ry;
            ctx.base.rz = p->rz;
            for (j = 0; j < N; j++)
                ctx.dir[j] = dirs[i][j];

            double f_before = f_curr;
            f_opt = line_minimize(&ctx, steps[i], &t_opt);

            /* Update p along direction */
            p->tx += t_opt * dirs[i][0];
            p->ty += t_opt * dirs[i][1];
            p->tz += t_opt * dirs[i][2];
            p->rx += t_opt * dirs[i][3];
            p->ry += t_opt * dirs[i][4];
            p->rz += t_opt * dirs[i][5];
            f_curr = f_opt;

            delta_i = f_before - f_curr;
            if (delta_i > delta_max)
            {
                delta_max = delta_i;
                idx_max = i;
            }
        }

        if (verbose)
            printf("Powell iter %3d  cost=%.6f  delta=%.2e\n",
                   iter + 1, f_curr, f_prev - f_curr);

        /* Convergence check */
        if (2.0 * fabs(f_prev - f_curr) <=
            tol * (fabs(f_prev) + fabs(f_curr) + 1e-20))
            break;

        /* Conjugate direction update (Powell rule):
         * Replace direction idx_max with the average direction taken. */
        double new_dir[6];
        new_dir[0] = p->tx - (p_arr[0] + dirs[0][0] * 0.0); /* placeholder: rebuild below */

        /* Rebuild new_dir = current p - p at start of this iteration */
        double pstart[6];
        pstart[0] = p_arr[0];
        pstart[1] = p_arr[1];
        pstart[2] = p_arr[2];
        pstart[3] = p_arr[3];
        pstart[4] = p_arr[4];
        pstart[5] = p_arr[5];

        double p_now[6];
        p_now[0] = p->tx;
        p_now[1] = p->ty;
        p_now[2] = p->tz;
        p_now[3] = p->rx;
        p_now[4] = p->ry;
        p_now[5] = p->rz;

        double len2 = 0.0;
        for (j = 0; j < N; j++)
        {
            new_dir[j] = p_now[j] - pstart[j];
            len2 += new_dir[j] * new_dir[j];
        }
        if (len2 > 1e-20)
        {
            double inv_len = 1.0 / sqrt(len2);
            for (j = 0; j < N; j++)
                new_dir[j] *= inv_len;

            /* Drop direction idx_max, shift remaining dirs AND their steps */
            for (i = idx_max; i < N - 1; i++)
            {
                for (j = 0; j < N; j++)
                    dirs[i][j] = dirs[i + 1][j];
                steps[i] = steps[i + 1];
            }
            for (j = 0; j < N; j++)
                dirs[N - 1][j] = new_dir[j];

            /* Adaptive step for the new conjugate direction based on its
             * translational vs rotational content, so that the line minimizer
             * explores a physically meaningful range regardless of direction. */
            double trans_sq = 0.0, rot_sq = 0.0;
            for (j = 0; j < 3; j++) trans_sq += new_dir[j] * new_dir[j];
            for (j = 3; j < N; j++) rot_sq   += new_dir[j] * new_dir[j];
            double step_t = (trans_sq > 1e-10) ? 0.5  / sqrt(trans_sq) : 1e30;
            double step_r = (rot_sq   > 1e-10) ? 0.005 / sqrt(rot_sq)  : 1e30;
            steps[N - 1] = (step_t < step_r) ? step_t : step_r;
            if (steps[N - 1] < 1e-8 || steps[N - 1] > 1.0)
                steps[N - 1] = 0.5;
        }
    }

    return f_curr;
}

/**
 * \brief Optimise the BBR cost function using a two-stage strategy.
 *
 * Stage 1 — brute-force grid search over a coarse ±range neighbourhood.
 * Stage 2 — Powell's direction-set method to refine.
 *
 * \param p_init          (in)     initial 6-DOF parameters
 * \param p_best          (out)    optimised 6-DOF parameters
 * \param surface         (in)     WM surface in surface-native RAS space
 * \param vol             (in)     floating-point volume data (moving image)
 * \param nii_ptr         (in)     NIfTI header of the moving volume
 * \param dims            (in)     volume dimensions [nx, ny, nz]
 * \param wm_dist         (in)     WM sampling offset (mm)
 * \param gm_dist         (in)     GM sampling offset (mm)
 * \param slope           (in)     BBR cost saturation slope
 * \param invert_contrast (in)     non-zero for T2/BOLD
 * \param grid_range_mm   (in)     half-width of translation grid (mm)
 * \param grid_range_rad  (in)     half-width of rotation grid (rad)
 * \param grid_steps      (in)     steps per DOF in grid search
 * \param max_iter        (in)     maximum Powell iterations
 * \param tol             (in)     convergence tolerance for Powell
 * \param verbose         (in)     non-zero to print progress
 * \return final BBR cost after optimisation
 */
double
CAT_BBReg_optimise(const CAT_RigidParams *p_init,
                   CAT_RigidParams *p_best,
                   const CAT_SurfData *surfs,
                   int n_surfs,
                   float *vol,
                   nifti_image *nii_ptr,
                   int dims[3],
                   double wm_dist,
                   double gm_dist,
                   double slope,
                   int invert_contrast,
                   double grid_range_mm,
                   double grid_range_rad,
                   int grid_steps,
                   int max_iter,
                   double tol,
                   int verbose)
{
    /* Stage 1: brute-force grid search
     * Only search over translations (tx, ty, tz) and one rotation axis (rz)
     * to keep the search tractable when grid_steps > 1.
     * Full 6-D grid would be (2*grid_steps+1)^6 which is very expensive
     * for grid_steps > 2.  We restrict Stage 1 to translations only. */
    double f_best = DBL_MAX;
    *p_best = *p_init;

    /* Evaluate cost at the starting point so the user can compare */
    double f_init = CAT_BBReg_cost(p_init, surfs, n_surfs, vol, nii_ptr, dims,
                                   wm_dist, gm_dist, slope, invert_contrast);

    if (verbose)
    {
        printf("BBReg initial cost (at identity/p_init): %.6f\n", f_init);
        printf("BBReg Stage 1: brute-force grid search...\n");
    }

    double step_mm = (grid_steps > 0) ? grid_range_mm / grid_steps : 0.0;
    double step_rad = (grid_steps > 0) ? grid_range_rad / grid_steps : 0.0;

    int n_steps = 2 * grid_steps + 1;
    int itx, ity, itz;

    for (itx = 0; itx < n_steps; itx++)
    {
        double dtx = (itx - grid_steps) * step_mm;
        for (ity = 0; ity < n_steps; ity++)
        {
            double dty = (ity - grid_steps) * step_mm;
            for (itz = 0; itz < n_steps; itz++)
            {
                double dtz = (itz - grid_steps) * step_mm;

                CAT_RigidParams p_try = *p_init;
                p_try.tx += dtx;
                p_try.ty += dty;
                p_try.tz += dtz;

                double f = CAT_BBReg_cost(&p_try, surfs, n_surfs, vol, nii_ptr,
                                          dims, wm_dist, gm_dist,
                                          slope, invert_contrast);
                if (f < f_best)
                {
                    f_best = f;
                    *p_best = p_try;
                }
            }
        }
    }

    if (verbose)
        printf("BBReg Stage 1 best cost: %.6f  (tx=%.2f ty=%.2f tz=%.2f)\n",
               f_best, p_best->tx, p_best->ty, p_best->tz);

    /* Stage 2: Powell refinement */
    if (verbose)
        printf("BBReg Stage 2: Powell optimisation...\n");

    double f_final = powell_minimise(p_best, surfs, n_surfs, vol, nii_ptr, dims,
                                     wm_dist, gm_dist, slope, invert_contrast,
                                     max_iter, tol, verbose);

    if (verbose)
    {
        printf("BBReg Stage 2 final cost: %.6f  "
               "(tx=%.3f ty=%.3f tz=%.3f  rx=%.4f ry=%.4f rz=%.4f)\n",
               f_final,
               p_best->tx, p_best->ty, p_best->tz,
               p_best->rx, p_best->ry, p_best->rz);
        bbr_print_stats(p_best, surfs, n_surfs, vol, nii_ptr, dims,
                        wm_dist, gm_dist);
    }

    return f_final;
}

/**
 * \brief Write a 4×4 RAS-to-RAS transform to a plain-text file.
 *
 * Format (compatible with FSL flirt -init / ANTs):
 *
 *   # CAT_BBReg RAS-to-RAS transform
 *   m00 m01 m02 m03
 *   m10 m11 m12 m13
 *   m20 m21 m22 m23
 *   0   0   0   1
 *
 * \param filename (in) output file path
 * \param m        (in) 4×4 row-major matrix
 * \return 0 on success, -1 on I/O error
 */
int CAT_BBReg_write_matrix(const char *filename, const double m[16])
{
    FILE *fp = fopen(filename, "w");
    if (!fp)
    {
        fprintf(stderr, "CAT_BBReg_write_matrix: cannot open '%s'\n", filename);
        return -1;
    }

    fprintf(fp, "# CAT_BBReg RAS-to-RAS transform\n");
    fprintf(fp, "# Generated by CAT-Surface CAT_SurfBBReg\n");
    fprintf(fp, "%.10g %.10g %.10g %.10g\n", m[0], m[1], m[2], m[3]);
    fprintf(fp, "%.10g %.10g %.10g %.10g\n", m[4], m[5], m[6], m[7]);
    fprintf(fp, "%.10g %.10g %.10g %.10g\n", m[8], m[9], m[10], m[11]);
    fprintf(fp, "%.10g %.10g %.10g %.10g\n", m[12], m[13], m[14], m[15]);

    fclose(fp);
    return 0;
}
