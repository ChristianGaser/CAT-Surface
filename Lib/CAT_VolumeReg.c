/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Robust multi-resolution volume-to-volume rigid registration.
 *
 * Algorithm: linearised intensity matching with Tukey biweight M-estimation
 * (Gauss-Newton IRLS), solved via 6×6 normal equations at each resolution
 * level of a Gaussian image pyramid.
 *
 * Reference: Reuter M, Rosas HD, Fischl B (2010). Highly Accurate Inverse
 *   Consistent Registration: A Robust Approach. NeuroImage 53:1181-1196.
 *
 * Copyright Christian Gaser, University of Jena.
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>

#include "CAT_VolumeReg.h"
#include "CAT_BBReg.h"

/* =========================================================================
 * Internal types
 * ========================================================================= */

typedef struct
{
    float *data;
    int nx, ny, nz;
    double sxyz[16]; /* voxel→world, row-major */
    double sijk[16]; /* world→voxel, row-major */
} VLevel;

/* =========================================================================
 * Low-level helpers
 * ========================================================================= */

/* Extract a row-major 16-element copy of nifti mat44 */
static void mat44_to_d16(mat44 m, double out[16])
{
    int r, c;
    for (r = 0; r < 4; r++)
        for (c = 0; c < 4; c++)
            out[r * 4 + c] = (double)m.m[r][c];
}

/* Invert a 4×4 rigid matrix (R^T | -R^T t) */
static void invert_rigid(const double m[16], double inv[16])
{
    inv[0] = m[0]; inv[1] = m[4]; inv[2] = m[8];
    inv[4] = m[1]; inv[5] = m[5]; inv[6] = m[9];
    inv[8] = m[2]; inv[9] = m[6]; inv[10] = m[10];
    inv[3]  = -(inv[0]*m[3] + inv[1]*m[7] + inv[2]*m[11]);
    inv[7]  = -(inv[4]*m[3] + inv[5]*m[7] + inv[6]*m[11]);
    inv[11] = -(inv[8]*m[3] + inv[9]*m[7] + inv[10]*m[11]);
    inv[12] = inv[13] = inv[14] = 0.0; inv[15] = 1.0;
}

/* Trilinear sample; returns NaN if out of bounds */
static double sample_vol(const float *vol, int nx, int ny, int nz,
                         double vx, double vy, double vz)
{
    int x0 = (int)floor(vx), y0 = (int)floor(vy), z0 = (int)floor(vz);
    int x1 = x0 + 1, y1 = y0 + 1, z1 = z0 + 1;
    if (x0 < 0 || x1 >= nx || y0 < 0 || y1 >= ny || z0 < 0 || z1 >= nz)
        return (double)NAN;
    double wx = vx - x0, wy = vy - y0, wz = vz - z0;
#define V(a,b,c) ((double)vol[(a) + nx*((b) + ny*(c))])
    return V(x0,y0,z0)*(1-wx)*(1-wy)*(1-wz) + V(x1,y0,z0)*wx*(1-wy)*(1-wz)
         + V(x0,y1,z0)*(1-wx)*wy*(1-wz)     + V(x1,y1,z0)*wx*wy*(1-wz)
         + V(x0,y0,z1)*(1-wx)*(1-wy)*wz     + V(x1,y0,z1)*wx*(1-wy)*wz
         + V(x0,y1,z1)*(1-wx)*wy*wz         + V(x1,y1,z1)*wx*wy*wz;
#undef V
}

/* 1-D separable Gaussian blur (in-place, kernel built on the fly) */
static void gauss1d(float *buf, float *tmp, int n, int stride, double sigma)
{
    if (sigma <= 0.0) return;
    int r = (int)(ceil(2.5 * sigma));
    if (r < 1) r = 1;
    int klen = 2 * r + 1;
    double *k = (double *)malloc(klen * sizeof(double));
    double s = 0.0;
    int i;
    for (i = 0; i < klen; i++)
    {
        double x = i - r;
        k[i] = exp(-0.5 * x * x / (sigma * sigma));
        s += k[i];
    }
    for (i = 0; i < klen; i++) k[i] /= s;

    for (i = 0; i < n; i++)
    {
        double acc = 0.0, wt = 0.0;
        int j;
        for (j = -r; j <= r; j++)
        {
            int idx = i + j;
            if (idx < 0 || idx >= n) continue;
            acc += k[j + r] * (double)buf[idx * stride];
            wt  += k[j + r];
        }
        tmp[i] = (float)(acc / wt);
    }
    for (i = 0; i < n; i++) buf[i * stride] = tmp[i];
    free(k);
}

/* Smooth a 3-D volume in-place with an isotropic Gaussian (sigma in voxels) */
static void smooth_vol(float *vol, int nx, int ny, int nz, double sigma)
{
    int x, y, z;
    float *tmp = (float *)malloc((nx > ny ? (nx > nz ? nx : nz) : (ny > nz ? ny : nz)) * sizeof(float));
    if (!tmp) return;

    /* x-pass */
    for (z = 0; z < nz; z++)
        for (y = 0; y < ny; y++)
            gauss1d(vol + y * nx + z * nx * ny, tmp, nx, 1, sigma);

    /* y-pass */
    for (z = 0; z < nz; z++)
        for (x = 0; x < nx; x++)
            gauss1d(vol + x + z * nx * ny, tmp, ny, nx, sigma);

    /* z-pass */
    for (y = 0; y < ny; y++)
        for (x = 0; x < nx; x++)
            gauss1d(vol + x + y * nx, tmp, nz, nx * ny, sigma);

    free(tmp);
}

/* Downsample volume by factor 2 in all dimensions (after prior smoothing) */
static float *downsample2(const float *vol, int nx, int ny, int nz,
                          int *out_nx, int *out_ny, int *out_nz)
{
    *out_nx = (nx + 1) / 2;
    *out_ny = (ny + 1) / 2;
    *out_nz = (nz + 1) / 2;
    float *out = (float *)malloc((size_t)(*out_nx) * (*out_ny) * (*out_nz) * sizeof(float));
    if (!out) return NULL;
    int ox, oy, oz;
    for (oz = 0; oz < *out_nz; oz++)
        for (oy = 0; oy < *out_ny; oy++)
            for (ox = 0; ox < *out_nx; ox++)
            {
                int sx = ox * 2, sy = oy * 2, sz = oz * 2;
                if (sx >= nx) sx = nx - 1;
                if (sy >= ny) sy = ny - 1;
                if (sz >= nz) sz = nz - 1;
                out[ox + (*out_nx)*(oy + (*out_ny)*oz)] =
                    vol[sx + nx*(sy + ny*sz)];
            }
    return out;
}

/* Build pyramid level: smooth then downsample 2^(level) times.
 * Level 0 = original; level k = 1/2^k resolution.
 * sxyz_in/sijk_in are the original world<->voxel matrices. */
static VLevel build_level(const float *vol, int nx, int ny, int nz,
                          const double sxyz_in[16],
                          int level)
{
    VLevel lv;
    float *cur = (float *)malloc((size_t)nx * ny * nz * sizeof(float));
    memcpy(cur, vol, (size_t)nx * ny * nz * sizeof(float));
    int cx = nx, cy = ny, cz = nz;

    int l;
    for (l = 0; l < level; l++)
    {
        /* Anti-aliasing: blur by 0.5 voxels before subsampling */
        smooth_vol(cur, cx, cy, cz, 0.5);
        int nx2, ny2, nz2;
        float *next = downsample2(cur, cx, cy, cz, &nx2, &ny2, &nz2);
        free(cur);
        cur = next;
        cx = nx2; cy = ny2; cz = nz2;
    }

    lv.data = cur;
    lv.nx = cx; lv.ny = cy; lv.nz = cz;

    /* Scale the voxel→world matrix: voxel size doubles at each level */
    int r, c;
    double scale = (double)(1 << level);
    for (r = 0; r < 4; r++)
        for (c = 0; c < 4; c++)
        {
            double v = sxyz_in[r * 4 + c];
            /* Scale the rotation/scale columns (0-2) */
            lv.sxyz[r * 4 + c] = (c < 3) ? v * scale : v;
        }

    /* Compute world→voxel by inverting sxyz (rigid, use transpose trick) */
    mat44 m44;
    for (r = 0; r < 4; r++)
        for (c = 0; c < 4; c++)
            m44.m[r][c] = (float)lv.sxyz[r * 4 + c];
    mat44 inv44 = nifti_mat44_inverse(m44);
    for (r = 0; r < 4; r++)
        for (c = 0; c < 4; c++)
            lv.sijk[r * 4 + c] = (double)inv44.m[r][c];

    return lv;
}

/* =========================================================================
 * Robust statistics
 * ========================================================================= */

static int cmp_float(const void *a, const void *b)
{
    float fa = *(const float *)a, fb = *(const float *)b;
    return (fa > fb) - (fa < fb);
}

/* Median of a float array (sorts a copy) */
static float median_f(float *arr, int n)
{
    float *tmp = (float *)malloc((size_t)n * sizeof(float));
    memcpy(tmp, arr, (size_t)n * sizeof(float));
    qsort(tmp, (size_t)n, sizeof(float), cmp_float);
    float med = (n & 1) ? tmp[n / 2] : 0.5f * (tmp[n/2 - 1] + tmp[n/2]);
    free(tmp);
    return med;
}

/* Percentile of a float array (sorts a copy); pct in [0,1] */
static float vol_percentile(const float *data, int n, float pct)
{
    float *tmp = (float *)malloc((size_t)n * sizeof(float));
    if (!tmp) return 0.0f;
    memcpy(tmp, data, (size_t)n * sizeof(float));
    qsort(tmp, (size_t)n, sizeof(float), cmp_float);
    int idx = (int)(pct * (float)(n - 1));
    if (idx < 0)     idx = 0;
    if (idx >= n)    idx = n - 1;
    float v = tmp[idx];
    free(tmp);
    return v;
}

/* Tukey biweight weight: w = (1-(r/sat)^2)^2, 0 if |r| >= sat */
static double tukey_w(double r, double sat)
{
    if (sat <= 0.0 || fabs(r) >= sat) return 0.0;
    double u = r / sat;
    double t = 1.0 - u * u;
    return t * t;
}

/* =========================================================================
 * 6×6 Cholesky solve (positive-definite symmetric ATA)
 * ========================================================================= */

static int cholesky6(double A[6][6], double L[6][6])
{
    int i, j, k;
    for (i = 0; i < 6; i++)
    {
        for (j = 0; j <= i; j++)
        {
            double s = A[i][j];
            for (k = 0; k < j; k++) s -= L[i][k] * L[j][k];
            if (i == j)
            {
                if (s <= 0.0) return -1; /* not positive definite */
                L[i][j] = sqrt(s);
            }
            else
                L[i][j] = s / L[j][j];
        }
        for (j = i + 1; j < 6; j++) L[i][j] = 0.0;
    }
    return 0;
}

/* Solve L L^T x = b */
static void cholesky6_solve(double L[6][6], double b[6], double x[6])
{
    int i, k;
    /* Forward substitution: L y = b */
    double y[6];
    for (i = 0; i < 6; i++)
    {
        double s = b[i];
        for (k = 0; k < i; k++) s -= L[i][k] * y[k];
        y[i] = s / L[i][i];
    }
    /* Back substitution: L^T x = y */
    for (i = 5; i >= 0; i--)
    {
        double s = y[i];
        for (k = i + 1; k < 6; k++) s -= L[k][i] * x[k];
        x[i] = s / L[i][i];
    }
}

/* =========================================================================
 * One registration level: IRLS loop
 * ========================================================================= */

static double register_level(const VLevel *fixed, const VLevel *moving,
                             CAT_RigidParams *p,
                             double sat_k, int max_iter, int verbose,
                             int level)
{
    int iter, x, y, z;
    const int fnx = fixed->nx, fny = fixed->ny, fnz = fixed->nz;
    double residual = 1.0;

    /* Pre-compute 3×3 upper-left of moving's world→voxel matrix (for
     * converting world gradient to voxel-space efficiently later) */

    /* Gradient of moving image in voxel space (precomputed once per level) */
    const int mnx = moving->nx, mny = moving->ny, mnz = moving->nz;
    float *gx_m = (float *)calloc((size_t)mnx * mny * mnz, sizeof(float));
    float *gy_m = (float *)calloc((size_t)mnx * mny * mnz, sizeof(float));
    float *gz_m = (float *)calloc((size_t)mnx * mny * mnz, sizeof(float));
    if (!gx_m || !gy_m || !gz_m) { free(gx_m); free(gy_m); free(gz_m); return 1.0; }

    /* Central-difference gradient in VOXEL space for the moving image */
    for (z = 0; z < mnz; z++)
        for (y = 0; y < mny; y++)
            for (x = 0; x < mnx; x++)
            {
#define MI(xx,yy,zz) ((double)moving->data[  \
    (xx) + mnx*((yy) + mny*(zz))])
                int idx = x + mnx*(y + mny*z);
                gx_m[idx] = (float)(0.5 * (
                    (x < mnx-1 ? MI(x+1,y,z) : MI(x,y,z)) -
                    (x > 0     ? MI(x-1,y,z) : MI(x,y,z))));
                gy_m[idx] = (float)(0.5 * (
                    (y < mny-1 ? MI(x,y+1,z) : MI(x,y,z)) -
                    (y > 0     ? MI(x,y-1,z) : MI(x,y,z))));
                gz_m[idx] = (float)(0.5 * (
                    (z < mnz-1 ? MI(x,y,z+1) : MI(x,y,z)) -
                    (z > 0     ? MI(x,y,z-1) : MI(x,y,z))));
#undef MI
            }

    /* Moving image sto_ijk upper-left 3×3 (world→voxel Jacobian).
     * We need this to convert gradient from voxel to world units. */
    /* Actually we convert gradient to world coords by multiplying the voxel
     * gradient by the world→voxel Jacobian transposed:
     *   g_world[k] = sum_j g_vox[j] * sijk[j][k]   (j,k = 0..2) */
    double Jm[3][3]; /* rows: output vox axis; cols: world axis */
    {
        int r2, c2;
        for (r2 = 0; r2 < 3; r2++)
            for (c2 = 0; c2 < 3; c2++)
                Jm[r2][c2] = moving->sijk[r2 * 4 + c2];
    }

    for (iter = 0; iter < max_iter; iter++)
    {
        /* Build 4×4 rigid matrix from current parameters */
        double M[16];
        CAT_BBReg_params_to_matrix(p, M);

        /* ---- Accumulate residuals for MAD sigma estimation ---- */
        int n_res = 0;
        float *res_arr = (float *)malloc(
            (size_t)fnx * fny * fnz * sizeof(float));
        if (!res_arr) break;

        for (z = 0; z < fnz; z++)
            for (y = 0; y < fny; y++)
                for (x = 0; x < fnx; x++)
                {
                    /* Fixed voxel → world */
                    double fw[3];
                    fw[0] = fixed->sxyz[0]*x + fixed->sxyz[1]*y +
                            fixed->sxyz[2]*z + fixed->sxyz[3];
                    fw[1] = fixed->sxyz[4]*x + fixed->sxyz[5]*y +
                            fixed->sxyz[6]*z + fixed->sxyz[7];
                    fw[2] = fixed->sxyz[8]*x + fixed->sxyz[9]*y +
                            fixed->sxyz[10]*z + fixed->sxyz[11];

                    /* Apply transform: fixed world → moving world */
                    double mw[3];
                    mw[0] = M[0]*fw[0] + M[1]*fw[1] + M[2]*fw[2] + M[3];
                    mw[1] = M[4]*fw[0] + M[5]*fw[1] + M[6]*fw[2] + M[7];
                    mw[2] = M[8]*fw[0] + M[9]*fw[1] + M[10]*fw[2] + M[11];

                    /* Moving world → moving voxel */
                    double mvx = moving->sijk[0]*mw[0] + moving->sijk[1]*mw[1] +
                                 moving->sijk[2]*mw[2] + moving->sijk[3];
                    double mvy = moving->sijk[4]*mw[0] + moving->sijk[5]*mw[1] +
                                 moving->sijk[6]*mw[2] + moving->sijk[7];
                    double mvz = moving->sijk[8]*mw[0] + moving->sijk[9]*mw[1] +
                                 moving->sijk[10]*mw[2] + moving->sijk[11];

                    double I_m = sample_vol(moving->data, mnx, mny, mnz, mvx, mvy, mvz);
                    if (isnan(I_m)) continue;

                    double I_f = (double)fixed->data[x + fnx*(y + fny*z)];
                    res_arr[n_res++] = (float)(I_m - I_f);
                }

        if (n_res < 10) { free(res_arr); break; }

        /* Robust sigma via MAD */
        float med_r = median_f(res_arr, n_res);
        float *abs_dev = (float *)malloc((size_t)n_res * sizeof(float));
        int i;
        for (i = 0; i < n_res; i++)
            abs_dev[i] = fabsf(res_arr[i] - med_r);
        float mad = median_f(abs_dev, n_res);
        free(abs_dev);
        free(res_arr);
        double sigma = 1.4826 * (double)mad;
        if (sigma < 1e-6) sigma = 1e-6;
        double sat = sat_k * sigma;

        /* ---- Build weighted normal equations ---- */
        double ATA[6][6];
        double ATb[6];
        memset(ATA, 0, sizeof(ATA));
        memset(ATb, 0, sizeof(ATb));
        double wsum = 0.0, wres2 = 0.0;
        int n_valid = 0;

        for (z = 0; z < fnz; z++)
            for (y = 0; y < fny; y++)
                for (x = 0; x < fnx; x++)
                {
                    /* Fixed voxel → world */
                    double fw[3];
                    fw[0] = fixed->sxyz[0]*x + fixed->sxyz[1]*y +
                            fixed->sxyz[2]*z + fixed->sxyz[3];
                    fw[1] = fixed->sxyz[4]*x + fixed->sxyz[5]*y +
                            fixed->sxyz[6]*z + fixed->sxyz[7];
                    fw[2] = fixed->sxyz[8]*x + fixed->sxyz[9]*y +
                            fixed->sxyz[10]*z + fixed->sxyz[11];

                    /* Apply transform: fixed world → moving world */
                    double mw[3];
                    mw[0] = M[0]*fw[0] + M[1]*fw[1] + M[2]*fw[2] + M[3];
                    mw[1] = M[4]*fw[0] + M[5]*fw[1] + M[6]*fw[2] + M[7];
                    mw[2] = M[8]*fw[0] + M[9]*fw[1] + M[10]*fw[2] + M[11];

                    /* Moving world → moving voxel */
                    double mvx = moving->sijk[0]*mw[0] + moving->sijk[1]*mw[1] +
                                 moving->sijk[2]*mw[2] + moving->sijk[3];
                    double mvy = moving->sijk[4]*mw[0] + moving->sijk[5]*mw[1] +
                                 moving->sijk[6]*mw[2] + moving->sijk[7];
                    double mvz = moving->sijk[8]*mw[0] + moving->sijk[9]*mw[1] +
                                 moving->sijk[10]*mw[2] + moving->sijk[11];

                    double I_m = sample_vol(moving->data, mnx, mny, mnz, mvx, mvy, mvz);
                    if (isnan(I_m)) continue;

                    double I_f = (double)fixed->data[x + fnx*(y + fny*z)];
                    double ri  = I_m - I_f;

                    double w = tukey_w(ri - (double)med_r, sat);
                    if (w < 1e-12) continue;
                    double w2 = w * w;

                    /* Gradient of moving image at (mvx, mvy, mvz) in WORLD coords.
                     * Interpolate voxel-space gradients then convert via Jacobian. */
                    double gvx_m = sample_vol(gx_m, mnx, mny, mnz, mvx, mvy, mvz);
                    double gvy_m = sample_vol(gy_m, mnx, mny, mnz, mvx, mvy, mvz);
                    double gvz_m = sample_vol(gz_m, mnx, mny, mnz, mvx, mvy, mvz);
                    if (isnan(gvx_m) || isnan(gvy_m) || isnan(gvz_m)) continue;

                    /* Convert voxel gradient → world gradient:
                     * g_world[k] = sum_j g_vox[j] * Jm[j][k]  */
                    double gw[3];
                    int k;
                    for (k = 0; k < 3; k++)
                        gw[k] = gvx_m * Jm[0][k] + gvy_m * Jm[1][k] +
                                gvz_m * Jm[2][k];

                    /* Design matrix row:
                     *   a[0..2] = gw[0..2]                (translation)
                     *   a[3]    = gw[2]*mw[1] - gw[1]*mw[2]  (rot x)
                     *   a[4]    = gw[0]*mw[2] - gw[2]*mw[0]  (rot y)
                     *   a[5]    = gw[1]*mw[0] - gw[0]*mw[1]  (rot z)
                     * Position used is mw[] = position in MOVING world space. */
                    double a[6];
                    a[0] = gw[0];
                    a[1] = gw[1];
                    a[2] = gw[2];
                    a[3] = gw[2] * mw[1] - gw[1] * mw[2];
                    a[4] = gw[0] * mw[2] - gw[2] * mw[0];
                    a[5] = gw[1] * mw[0] - gw[0] * mw[1];

                    /* Accumulate ATA and ATb.
                     * Gauss-Newton: ATA dp = -A^T W r  (descent toward minimum)
                     * so accumulate the negative of A^T W r into ATb. */
                    int r2, c2;
                    for (r2 = 0; r2 < 6; r2++)
                    {
                        ATb[r2] -= w2 * a[r2] * ri;
                        for (c2 = 0; c2 <= r2; c2++)
                            ATA[r2][c2] += w2 * a[r2] * a[c2];
                    }

                    wsum  += w2;
                    wres2 += w2 * ri * ri;
                    n_valid++;
                }

        if (n_valid < 10) break;

        /* Symmetrize ATA */
        {
            int r2, c2;
            for (r2 = 0; r2 < 6; r2++)
                for (c2 = r2 + 1; c2 < 6; c2++)
                    ATA[r2][c2] = ATA[c2][r2];
        }

        /* Regularise diagonal slightly for numerical stability */
        {
            int j;
            double tr = 0.0;
            for (j = 0; j < 6; j++) tr += ATA[j][j];
            double eps_reg = 1e-8 * tr / 6.0;
            for (j = 0; j < 6; j++) ATA[j][j] += eps_reg;
        }

        /* Cholesky solve: ATA * dp = -A^T W r  (stored as positive ATb above) */
        double L[6][6];
        if (cholesky6(ATA, L) != 0) break;
        double dp[6];
        cholesky6_solve(L, ATb, dp);

        /* Update parameters */
        p->tx += dp[0]; p->ty += dp[1]; p->tz += dp[2];
        p->rx += dp[3]; p->ry += dp[4]; p->rz += dp[5];

        residual = (wsum > 0.0) ? sqrt(wres2 / wsum) : 1.0;

        if (verbose > 1)
            printf("  Level %d  iter %2d  residual=%.4f  sigma=%.2f"
                   "  n_valid=%d  dp_t=%.3f dp_r=%.5f\n",
                   level, iter + 1, residual, sigma, n_valid,
                   sqrt(dp[0]*dp[0]+dp[1]*dp[1]+dp[2]*dp[2]),
                   sqrt(dp[3]*dp[3]+dp[4]*dp[4]+dp[5]*dp[5]));

        /* Convergence: small parameter update */
        double dt = sqrt(dp[0]*dp[0] + dp[1]*dp[1] + dp[2]*dp[2]);
        double dr = sqrt(dp[3]*dp[3] + dp[4]*dp[4] + dp[5]*dp[5]);
        if (dt < 0.01 && dr < 1e-4) break;
    }

    free(gx_m); free(gy_m); free(gz_m);
    return residual;
}

/* =========================================================================
 * Public API
 * ========================================================================= */

/**
 * \brief Robust multi-resolution volume-to-volume rigid registration.
 */
double
CAT_VolumeReg_register(float *fixed_vol,
                       nifti_image *fixed_nii,
                       int fixed_dims[3],
                       float *moving_vol,
                       nifti_image *moving_nii,
                       int moving_dims[3],
                       CAT_RigidParams *p_out,
                       int n_levels,
                       double sat_k,
                       int max_iter,
                       int verbose)
{
    double fsxyz[16], msxyz[16];
    mat44_to_d16(fixed_nii->sto_xyz,  fsxyz);
    mat44_to_d16(moving_nii->sto_xyz, msxyz);

    /* Clamp n_levels to a sensible range */
    if (n_levels < 1) n_levels = 1;
    if (n_levels > 6) n_levels = 6;

    /* Build pyramids for fixed and moving */
    VLevel *fpyr = (VLevel *)malloc((size_t)n_levels * sizeof(VLevel));
    VLevel *mpyr = (VLevel *)malloc((size_t)n_levels * sizeof(VLevel));
    if (!fpyr || !mpyr) { free(fpyr); free(mpyr); return 1.0; }

    int l;
    for (l = 0; l < n_levels; l++)
    {
        fpyr[l] = build_level(fixed_vol,  fixed_dims[0],  fixed_dims[1],  fixed_dims[2],  fsxyz, l);
        mpyr[l] = build_level(moving_vol, moving_dims[0], moving_dims[1], moving_dims[2], msxyz, l);
    }

    /* Initialise parameters from caller (or zero if unset) */
    CAT_RigidParams p = *p_out;

    double final_res = 1.0;

    /* Coarse-to-fine: level (n_levels-1) → 0 */
    for (l = n_levels - 1; l >= 0; l--)
    {
        if (verbose)
            printf("VolumeReg level %d  (%dx%dx%d fixed, %dx%dx%d moving)\n",
                   l,
                   fpyr[l].nx, fpyr[l].ny, fpyr[l].nz,
                   mpyr[l].nx, mpyr[l].ny, mpyr[l].nz);

        final_res = register_level(&fpyr[l], &mpyr[l],
                                   &p, sat_k, max_iter, verbose, l);

        if (verbose)
            printf("  → residual=%.4f  tx=%.2f ty=%.2f tz=%.2f"
                   "  rx=%.4f ry=%.4f rz=%.4f\n",
                   final_res, p.tx, p.ty, p.tz, p.rx, p.ry, p.rz);
    }

    *p_out = p;

    /* Free pyramids */
    for (l = 0; l < n_levels; l++) { free(fpyr[l].data); free(mpyr[l].data); }
    free(fpyr); free(mpyr);

    return final_res;
}

/* =========================================================================
 * NMI-based rigid registration  (cross-modal, e.g. EPI ↔ T1w)
 *
 * Algorithm:
 *   1. Build Gaussian pyramid for fixed and moving (shared with IRLS path).
 *   2. At each pyramid level (coarse→fine):
 *      a. Map fixed intensity and sampled moving intensity to joint 2-D
 *         histogram using bilinear partial-volume interpolation.
 *      b. Compute NMI = (H_f + H_m) / H_joint from the histogram.
 *      c. Maximise NMI with Powell's direction-set method + Brent line search.
 * ========================================================================= */

/* ---- Internal types ---- */

/* Context shared between the NMI cost function and the line minimizer */
typedef struct
{
    CAT_RigidParams  base;      /* current base point in parameter space     */
    double           dir[6];    /* unit search direction                      */
    const VLevel    *fixed;
    const VLevel    *moving;
    int              n_bins;    /* histogram resolution per axis              */
    float            flo;       /* fixed intensity → bin: bin = (I-flo)*fscale */
    float            fscale;
    float            mlo;       /* moving intensity → bin                    */
    float            mscale;
} NMILineCtx;

/* ---- Histogram helpers ---- */

/* Map intensity to a fractional bin in [0, n_bins-1] */
static double intensity_to_bin(double intensity, float lo, float scale,
                                int n_bins)
{
    double b = ((double)intensity - (double)lo) * (double)scale;
    if (b < 0.0)                b = 0.0;
    if (b > (double)(n_bins - 1)) b = (double)(n_bins - 1);
    return b;
}

/* Compute NMI from a (nb × nb) joint histogram.
 * Rows = fixed bins, columns = moving bins.
 * Returns NMI ≥ 1; returns 0 on empty or degenerate histogram. */
static double nmi_from_hist(const double *jhist, int nb)
{
    int i, j;
    double total = 0.0;
    for (i = 0; i < nb * nb; i++) total += jhist[i];
    if (total < 1.0) return 0.0;
    double inv = 1.0 / total;

    double hf = 0.0, hm = 0.0, hj = 0.0;

    /* Marginal entropy of fixed */
    for (i = 0; i < nb; i++)
    {
        double pi = 0.0;
        for (j = 0; j < nb; j++) pi += jhist[i * nb + j];
        pi *= inv;
        if (pi > 1e-15) hf -= pi * log(pi);
    }
    /* Marginal entropy of moving */
    for (j = 0; j < nb; j++)
    {
        double pj = 0.0;
        for (i = 0; i < nb; i++) pj += jhist[i * nb + j];
        pj *= inv;
        if (pj > 1e-15) hm -= pj * log(pj);
    }
    /* Joint entropy */
    for (i = 0; i < nb * nb; i++)
    {
        double p = jhist[i] * inv;
        if (p > 1e-15) hj -= p * log(p);
    }

    if (hj < 1e-10) return 0.0;
    return (hf + hm) / hj;   /* NMI ≥ 1; larger = better alignment */
}

/* ---- NMI cost (returns −NMI for minimisation) ---- */

static double nmi_cost_params(const CAT_RigidParams *p,
                              const NMILineCtx *ctx)
{
    double M[16];
    CAT_BBReg_params_to_matrix((CAT_RigidParams *)p, M);

    const VLevel *fixed  = ctx->fixed;
    const VLevel *moving = ctx->moving;
    const int fnx = fixed->nx,  fny = fixed->ny,  fnz = fixed->nz;
    const int mnx = moving->nx, mny = moving->ny, mnz = moving->nz;
    const int nb  = ctx->n_bins;

    double *jhist = (double *)calloc((size_t)nb * nb, sizeof(double));
    if (!jhist) return 0.0;

    int x, y, z;
    for (z = 0; z < fnz; z++)
    for (y = 0; y < fny; y++)
    for (x = 0; x < fnx; x++)
    {
        /* Fixed voxel → world */
        double fw0 = fixed->sxyz[0]*x + fixed->sxyz[1]*y +
                     fixed->sxyz[2]*z + fixed->sxyz[3];
        double fw1 = fixed->sxyz[4]*x + fixed->sxyz[5]*y +
                     fixed->sxyz[6]*z + fixed->sxyz[7];
        double fw2 = fixed->sxyz[8]*x + fixed->sxyz[9]*y +
                     fixed->sxyz[10]*z + fixed->sxyz[11];

        /* Apply rigid transform: fixed world → moving world */
        double mw0 = M[0]*fw0 + M[1]*fw1 + M[2]*fw2 + M[3];
        double mw1 = M[4]*fw0 + M[5]*fw1 + M[6]*fw2 + M[7];
        double mw2 = M[8]*fw0 + M[9]*fw1 + M[10]*fw2 + M[11];

        /* Moving world → moving voxel */
        double mvx = moving->sijk[0]*mw0 + moving->sijk[1]*mw1 +
                     moving->sijk[2]*mw2 + moving->sijk[3];
        double mvy = moving->sijk[4]*mw0 + moving->sijk[5]*mw1 +
                     moving->sijk[6]*mw2 + moving->sijk[7];
        double mvz = moving->sijk[8]*mw0 + moving->sijk[9]*mw1 +
                     moving->sijk[10]*mw2 + moving->sijk[11];

        double I_m = sample_vol(moving->data, mnx, mny, mnz, mvx, mvy, mvz);
        if (isnan(I_m)) continue;

        double I_f = (double)fixed->data[x + fnx * (y + fny * z)];

        /* Fractional bins */
        double bf = intensity_to_bin(I_f, ctx->flo, ctx->fscale, nb);
        double bm = intensity_to_bin(I_m, ctx->mlo, ctx->mscale, nb);

        /* Bilinear partial-volume interpolation into the joint histogram */
        int bf0 = (int)bf, bm0 = (int)bm;
        int bf1 = bf0 + 1, bm1 = bm0 + 1;
        double wf1 = bf - (double)bf0, wm1 = bm - (double)bm0;
        double wf0 = 1.0 - wf1,        wm0 = 1.0 - wm1;

        if (bf0 >= 0 && bf0 < nb)
        {
            if (bm0 >= 0 && bm0 < nb) jhist[bf0 * nb + bm0] += wf0 * wm0;
            if (bm1 >= 0 && bm1 < nb) jhist[bf0 * nb + bm1] += wf0 * wm1;
        }
        if (bf1 >= 0 && bf1 < nb)
        {
            if (bm0 >= 0 && bm0 < nb) jhist[bf1 * nb + bm0] += wf1 * wm0;
            if (bm1 >= 0 && bm1 < nb) jhist[bf1 * nb + bm1] += wf1 * wm1;
        }
    }

    double nmi = nmi_from_hist(jhist, nb);
    free(jhist);
    return -nmi;   /* minimise */
}

/* ---- Brent + Powell line search for NMI ---- */

static double nmi_line_cost(const NMILineCtx *ctx, double t)
{
    CAT_RigidParams p;
    p.tx = ctx->base.tx + t * ctx->dir[0];
    p.ty = ctx->base.ty + t * ctx->dir[1];
    p.tz = ctx->base.tz + t * ctx->dir[2];
    p.rx = ctx->base.rx + t * ctx->dir[3];
    p.ry = ctx->base.ry + t * ctx->dir[4];
    p.rz = ctx->base.rz + t * ctx->dir[5];
    return nmi_cost_params(&p, ctx);
}

/* Brent's 1-D minimiser; same algorithm as in CAT_BBReg.c */
static double nmi_brent(const NMILineCtx *ctx,
                        double ax, double bx, double cx,
                        double tol, double *fmin)
{
    int iter;
    double a, b, d = 0.0, e = 0.0, etemp, fu, fv, fw, fx;
    double p, q, r2, tol1, tol2, u, v, w, x, xm;

    a = (ax < cx) ? ax : cx;
    b = (ax > cx) ? ax : cx;
    x = w = v = bx;
    fw = fv = fx = nmi_line_cost(ctx, x);

    for (iter = 0; iter < 100; iter++)
    {
        xm   = 0.5 * (a + b);
        tol1 = tol * fabs(x) + 1.0e-10;
        tol2 = 2.0 * tol1;
        if (fabs(x - xm) <= (tol2 - 0.5 * (b - a))) { *fmin = fx; return x; }

        if (fabs(e) > tol1)
        {
            r2 = (x - w) * (fx - fv);
            q  = (x - v) * (fx - fw);
            p  = (x - v) * q - (x - w) * r2;
            q  = 2.0 * (q - r2);
            if (q > 0.0) p = -p;
            q    = fabs(q);
            etemp = e; e = d;
            if (fabs(p) >= fabs(0.5 * q * etemp) ||
                p <= q * (a - x) || p >= q * (b - x))
            {
                e = (x >= xm) ? a - x : b - x;
                d = 0.3819660 * e;
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
            d = 0.3819660 * e;
        }

        u  = (fabs(d) >= tol1) ? x + d : x + ((d >= 0.0) ? tol1 : -tol1);
        fu = nmi_line_cost(ctx, u);

        if (fu <= fx)
        {
            if (u < x) b = x; else a = x;
            v = w; fv = fw; w = x; fw = fx; x = u; fx = fu;
        }
        else
        {
            if (u < x) a = u; else b = u;
            if (fu <= fw || w == x) { v = w; fv = fw; w = u; fw = fu; }
            else if (fu <= fv || v == x || v == w) { v = u; fv = fu; }
        }
    }
    *fmin = fx;
    return x;
}

static double nmi_line_minimize(const NMILineCtx *ctx,
                                double step, double *t_out)
{
    double t0 = 0.0, t1 = step;
    double f0 = nmi_line_cost(ctx, t0);
    double f1 = nmi_line_cost(ctx, t1);
    double fmin;

    /* Ensure we step downhill first */
    if (f1 > f0) { t1 = -step; f1 = nmi_line_cost(ctx, t1); }
    if (f1 > f0) { *t_out = 0.0; return f0; }

    /* Expand bracket with golden ratio until we bracket a minimum */
    double t2 = t1 * 1.618, f2 = nmi_line_cost(ctx, t2);
    int it = 0;
    while (f2 < f1 && it++ < 50)
    {
        t0 = t1; f0 = f1;
        t1 = t2; f1 = f2;
        t2 = t1 + 1.618 * (t1 - t0);
        f2 = nmi_line_cost(ctx, t2);
    }

    double tmin = nmi_brent(ctx, t0, t1, t2, 1e-4, &fmin);
    *t_out = tmin;
    return fmin;
}

/* Powell's direction-set minimiser for NMI.
 * Uses the same conjugate-direction update as the BBR Powell in CAT_BBReg.c. */
static double nmi_powell(CAT_RigidParams *p,
                         const NMILineCtx *ctx_tmpl,
                         const double init_steps[6],
                         int max_iter, double tol,
                         int verbose, int level)
{
    const int N = 6;
    double steps[6], dirs[6][6];
    int i, j, iter;

    for (i = 0; i < N; i++) steps[i] = init_steps[i];
    for (i = 0; i < N; i++)
        for (j = 0; j < N; j++)
            dirs[i][j] = (i == j) ? 1.0 : 0.0;

    NMILineCtx ctx = *ctx_tmpl;
    double f_curr = nmi_cost_params(p, &ctx);

    for (iter = 0; iter < max_iter; iter++)
    {
        double f_prev = f_curr;
        double delta_max = 0.0;
        int    idx_max   = 0;
        double p_arr[6];
        p_arr[0] = p->tx; p_arr[1] = p->ty; p_arr[2] = p->tz;
        p_arr[3] = p->rx; p_arr[4] = p->ry; p_arr[5] = p->rz;

        for (i = 0; i < N; i++)
        {
            ctx.base = *p;
            for (j = 0; j < N; j++) ctx.dir[j] = dirs[i][j];

            double f_before = f_curr, t_opt;
            double f_opt = nmi_line_minimize(&ctx, steps[i], &t_opt);

            p->tx += t_opt * dirs[i][0];
            p->ty += t_opt * dirs[i][1];
            p->tz += t_opt * dirs[i][2];
            p->rx += t_opt * dirs[i][3];
            p->ry += t_opt * dirs[i][4];
            p->rz += t_opt * dirs[i][5];
            f_curr = f_opt;

            double delta_i = f_before - f_curr;
            if (delta_i > delta_max) { delta_max = delta_i; idx_max = i; }
        }

        if (verbose > 1)
            printf("  NMI level %d  iter %2d  NMI=%.5f  delta=%.2e\n",
                   level, iter + 1, -f_curr, f_prev - f_curr);

        if (2.0 * fabs(f_prev - f_curr) <=
            tol * (fabs(f_prev) + fabs(f_curr) + 1e-20))
            break;

        /* Conjugate direction update (Powell rule) */
        double new_dir[6], len2 = 0.0;
        for (j = 0; j < N; j++)
        {
            double pnow = (j == 0) ? p->tx : (j == 1) ? p->ty : (j == 2) ? p->tz :
                          (j == 3) ? p->rx : (j == 4) ? p->ry : p->rz;
            new_dir[j] = pnow - p_arr[j];
            len2 += new_dir[j] * new_dir[j];
        }
        if (len2 > 1e-20)
        {
            double inv_len = 1.0 / sqrt(len2);
            for (j = 0; j < N; j++) new_dir[j] *= inv_len;

            for (i = idx_max; i < N - 1; i++)
            {
                for (j = 0; j < N; j++) dirs[i][j] = dirs[i + 1][j];
                steps[i] = steps[i + 1];
            }
            for (j = 0; j < N; j++) dirs[N - 1][j] = new_dir[j];

            /* Adaptive step: scale to physical magnitude of direction */
            double trans_sq = 0.0, rot_sq = 0.0;
            for (j = 0; j < 3; j++) trans_sq += new_dir[j] * new_dir[j];
            for (j = 3; j < N; j++) rot_sq   += new_dir[j] * new_dir[j];
            double step_t = (trans_sq > 1e-10) ?
                            init_steps[0] / sqrt(trans_sq) : 1e30;
            double step_r = (rot_sq   > 1e-10) ?
                            init_steps[3] / sqrt(rot_sq)   : 1e30;
            steps[N - 1] = (step_t < step_r) ? step_t : step_r;
            if (steps[N - 1] < 1e-8 || steps[N - 1] > 10.0)
                steps[N - 1] = init_steps[0];
        }
    }
    return f_curr;
}

/* =========================================================================
 * Public API — NMI registration
 * ========================================================================= */

/**
 * \brief Cross-modal rigid registration using Normalised Mutual Information.
 *
 * Equivalent to mri_coreg (FreeSurfer) which is the initialisation step used
 * by fMRIPrep before bbregister.  A Gaussian image pyramid is used for
 * coarse-to-fine registration; at each level the NMI is computed from a
 * joint intensity histogram (bilinear partial-volume interpolation) and
 * maximised with Powell's direction-set method.
 *
 * Intensity histograms use n_bins bins per axis.  Intensities are normalised
 * from the 1st–99th percentile of each volume so outlier voxels do not
 * distort the histogram.
 *
 * \param fixed_vol   (in)  reference volume data (float, x-fastest)
 * \param fixed_nii   (in)  NIfTI header of the reference volume (T1w)
 * \param fixed_dims  (in)  [nx, ny, nz] of the reference volume
 * \param moving_vol  (in)  moving volume data (float) — e.g. mean EPI/BOLD
 * \param moving_nii  (in)  NIfTI header of the moving volume
 * \param moving_dims (in)  [nx, ny, nz] of the moving volume
 * \param p_out       (out) estimated 6-DOF rigid parameters (fixed→moving)
 * \param n_levels    (in)  Gaussian pyramid levels (3–4 recommended)
 * \param n_bins      (in)  histogram bins per axis (64 recommended)
 * \param max_iter    (in)  Powell iterations per level (20–50)
 * \param verbose     (in)  ≥1 prints level summaries; ≥2 prints per-iteration
 * \return            final NMI value (higher = better alignment)
 */
double
CAT_VolumeReg_register_NMI(float *fixed_vol,
                            nifti_image *fixed_nii,
                            int fixed_dims[3],
                            float *moving_vol,
                            nifti_image *moving_nii,
                            int moving_dims[3],
                            CAT_RigidParams *p_out,
                            int n_levels,
                            int n_bins,
                            int max_iter,
                            int verbose)
{
    int l, j;
    double fsxyz[16], msxyz[16];
    mat44_to_d16(fixed_nii->sto_xyz,  fsxyz);
    mat44_to_d16(moving_nii->sto_xyz, msxyz);

    if (n_levels < 1)  n_levels = 1;
    if (n_levels > 6)  n_levels = 6;
    if (n_bins  < 16)  n_bins   = 16;
    if (n_bins  > 256) n_bins   = 256;

    /* Build Gaussian pyramids */
    VLevel *fpyr = (VLevel *)malloc((size_t)n_levels * sizeof(VLevel));
    VLevel *mpyr = (VLevel *)malloc((size_t)n_levels * sizeof(VLevel));
    if (!fpyr || !mpyr) { free(fpyr); free(mpyr); return 0.0; }

    for (l = 0; l < n_levels; l++)
    {
        fpyr[l] = build_level(fixed_vol,
                              fixed_dims[0],  fixed_dims[1],  fixed_dims[2],
                              fsxyz, l);
        mpyr[l] = build_level(moving_vol,
                              moving_dims[0], moving_dims[1], moving_dims[2],
                              msxyz, l);
    }

    /* Compute intensity normalisation from the coarsest level (fastest).
     * Clip to 1st–99th percentile to suppress outlier voxels. */
    int cl = n_levels - 1;
    int fn = fpyr[cl].nx * fpyr[cl].ny * fpyr[cl].nz;
    int mn = mpyr[cl].nx * mpyr[cl].ny * mpyr[cl].nz;

    float flo = vol_percentile(fpyr[cl].data, fn, 0.01f);
    float fhi = vol_percentile(fpyr[cl].data, fn, 0.99f);
    float mlo = vol_percentile(mpyr[cl].data, mn, 0.01f);
    float mhi = vol_percentile(mpyr[cl].data, mn, 0.99f);

    float fscale = (fhi > flo) ? (float)(n_bins - 1) / (fhi - flo) : 1.0f;
    float mscale = (mhi > mlo) ? (float)(n_bins - 1) / (mhi - mlo) : 1.0f;

    if (verbose)
        printf("NMI intensity range:  fixed [%.1f, %.1f]  moving [%.1f, %.1f]"
               "  bins=%d\n", (double)flo, (double)fhi,
               (double)mlo, (double)mhi, n_bins);

    CAT_RigidParams p = *p_out;
    double final_nmi  = 0.0;

    for (l = n_levels - 1; l >= 0; l--)
    {
        /* Step sizes scale with pyramid level: coarser = larger initial steps */
        double scale = (double)(1 << l);   /* 2^l */
        double init_steps[6];
        init_steps[0] = init_steps[1] = init_steps[2] = 2.0 * scale; /* mm */
        init_steps[3] = init_steps[4] = init_steps[5] = 0.02;        /* rad */

        NMILineCtx ctx;
        ctx.base   = p;
        for (j = 0; j < 6; j++) ctx.dir[j] = 0.0;
        ctx.fixed  = &fpyr[l];
        ctx.moving = &mpyr[l];
        ctx.n_bins = n_bins;
        ctx.flo    = flo;   ctx.fscale = fscale;
        ctx.mlo    = mlo;   ctx.mscale = mscale;

        if (verbose)
            printf("NMI level %d  (%dx%dx%d fixed, %dx%dx%d moving)"
                   "  step_t=%.1f mm\n",
                   l, fpyr[l].nx, fpyr[l].ny, fpyr[l].nz,
                   mpyr[l].nx, mpyr[l].ny, mpyr[l].nz,
                   init_steps[0]);

        double neg_nmi = nmi_powell(&p, &ctx, init_steps,
                                    max_iter, 1e-5, verbose, l);
        final_nmi = -neg_nmi;

        if (verbose)
            printf("  → NMI=%.5f  tx=%.2f ty=%.2f tz=%.2f"
                   "  rx=%.4f ry=%.4f rz=%.4f\n",
                   final_nmi, p.tx, p.ty, p.tz, p.rx, p.ry, p.rz);
    }

    *p_out = p;

    for (l = 0; l < n_levels; l++)
    {
        free(fpyr[l].data);
        free(mpyr[l].data);
    }
    free(fpyr);
    free(mpyr);

    return final_nmi;
}
