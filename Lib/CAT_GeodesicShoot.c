/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Geodesic shooting of tissue maps onto a series of shooting templates: a C
 * port of CAT12/CAT26's optimised shooting (run_Shooting and regsetup in
 * cat_main_registration.m, Christian Gaser and Robert Dahnke), built on SPM's
 * shooting primitives (3rdparty/spm_shoot, John Ashburner): the Gauss-Newton
 * update of spm_shoot_update.m, the geodesic integration of spm_shoot3d.m and
 * the Green's operator of spm_shoot_greens.m are ported here as well.
 *
 * The code follows the MATLAB statement by statement, including its
 * precision (single where MATLAB computes in single), so that results can be
 * compared with CAT26 voxel by voxel.  Where the MATLAB has a quirk that
 * changes results, the quirk is kept and commented.
 *
 * Copyright Christian Gaser, University of Jena.
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <fftw3.h>

#include "CAT_GeodesicShoot.h"
#include "CAT_Vol.h"
#include "CAT_Cat12Utils.h"
#include "spm_shoot.h"

/* MATLAB evaluates without fused multiply-add, and so were SPM's mex files
 * built; contracting here would only move the results away from CAT's. */
#pragma STDC FP_CONTRACT OFF

#define MATLAB_EPS 2.220446049250313e-16 /* eps in MATLAB (double) */
#define SPM_DIFFEO_EPS 1.1921e-7         /* eps as defined in spm_diffeo.c */
#define SLICE_TINY 5e-2                  /* TINY of spm_vol_utils.c */
#define NLEV CAT_SHOOT_NTEMPLATES
#define MAXCLS (CAT_SHOOT_MAXCLASSES + 2) /* + background + ventricle class */
#define MAX_SHOOT_THREADS 16

/* ======================================================================== */
/* Threading                                                                */
/* ======================================================================== */

typedef void (*range_fn)(void *ctx, int ini, int fin);

typedef struct
{
    range_fn fn;
    void *ctx;
    int ini, fin;
} range_args;

static void *
range_worker(void *p)
{
    range_args *a = (range_args *)p;
    a->fn(a->ctx, a->ini, a->fin);
    return NULL;
}

/* Run fn over [0, n) in nthreads contiguous chunks.  Every caller writes only
 * to the elements of its own chunk, so the result does not depend on the
 * number of threads. */
static void
parallel_for(int n, int nthreads, range_fn fn, void *ctx)
{
    range_args args[MAX_SHOOT_THREADS];
    int t;

    if (nthreads > n) nthreads = n;
    if (nthreads > MAX_SHOOT_THREADS) nthreads = MAX_SHOOT_THREADS;
    if (nthreads <= 1)
    {
        if (n > 0) fn(ctx, 0, n);
        return;
    }
    for (t = 0; t < nthreads; t++)
    {
        args[t].fn = fn;
        args[t].ctx = ctx;
        args[t].ini = (int)(((long)t * n) / nthreads);
        args[t].fin = (int)(((long)(t + 1) * n) / nthreads);
    }
    cat_parallel_run(nthreads, range_worker, args, sizeof(range_args));
}

/* ======================================================================== */
/* 4x4 matrices, row-major                                                  */
/* ======================================================================== */

static void
m4_eye(double *A)
{
    int i;
    for (i = 0; i < 16; i++) A[i] = (i % 5 == 0) ? 1.0 : 0.0;
}

static void
m4_mul(const double *A, const double *B, double *C)
{
    double T[16];
    int i, j, k;
    for (i = 0; i < 4; i++)
        for (j = 0; j < 4; j++)
        {
            double s = 0.0;
            for (k = 0; k < 4; k++) s += A[i * 4 + k] * B[k * 4 + j];
            T[i * 4 + j] = s;
        }
    memcpy(C, T, sizeof(T));
}

/* Inverse by Gauss-Jordan elimination with partial pivoting. */
static int
m4_inv(const double *A, double *Ai)
{
    double M[4][8];
    int i, j, k;

    for (i = 0; i < 4; i++)
        for (j = 0; j < 8; j++)
            M[i][j] = (j < 4) ? A[i * 4 + j] : ((j - 4 == i) ? 1.0 : 0.0);
    for (k = 0; k < 4; k++)
    {
        int p = k;
        double piv;
        for (i = k + 1; i < 4; i++)
            if (fabs(M[i][k]) > fabs(M[p][k])) p = i;
        if (fabs(M[p][k]) < 1e-300) return 1;
        if (p != k)
            for (j = 0; j < 8; j++)
            {
                double t = M[k][j];
                M[k][j] = M[p][j];
                M[p][j] = t;
            }
        piv = M[k][k];
        for (j = 0; j < 8; j++) M[k][j] /= piv;
        for (i = 0; i < 4; i++)
            if (i != k)
            {
                double f = M[i][k];
                for (j = 0; j < 8; j++) M[i][j] -= f * M[k][j];
            }
    }
    for (i = 0; i < 4; i++)
        for (j = 0; j < 4; j++) Ai[i * 4 + j] = M[i][j + 4];
    return 0;
}

/* NIfTI matrices map 0-based voxels, SPM's 1-based ones: Mspm = Mnii*T(-1). */
static void
nii2spm(const double *nii, double *spm)
{
    double T[16];
    m4_eye(T);
    T[3] = T[7] = T[11] = -1.0;
    m4_mul(nii, T, spm);
}

/* Least-squares X = A/B for A, B of size 4x8 (MATLAB's mrdivide), solved by
 * the normal equations X = A*B'*inv(B*B').  The corner sets used here are
 * related by an exact affine map, so this is that map. */
static int
mrdivide48(const double A[4][8], const double B[4][8], double *X)
{
    double ABt[16], BBt[16], BBti[16];
    int i, j, k;

    for (i = 0; i < 4; i++)
        for (j = 0; j < 4; j++)
        {
            double s1 = 0.0, s2 = 0.0;
            for (k = 0; k < 8; k++)
            {
                s1 += A[i][k] * B[j][k];
                s2 += B[i][k] * B[j][k];
            }
            ABt[i * 4 + j] = s1;
            BBt[i * 4 + j] = s2;
        }
    if (m4_inv(BBt, BBti)) return 1;
    m4_mul(ABt, BBti, X);
    return 0;
}

/* ======================================================================== */
/* spm_smooth (spm_smoothkern with linear interpolation, spm_conv_vol)       */
/* ======================================================================== */

/* spm_smoothkern(fwhm, -r:r, 1), normalised to unit sum as spm_smooth does. */
static void
smoothkern(double fwhm, int r, double *k)
{
    double s = pow(fwhm / sqrt(8.0 * log(2.0)), 2) + MATLAB_EPS;
    double w1 = 0.5 * sqrt(2.0 / s), w2 = -0.5 / s, w3 = sqrt(s / 2.0 / M_PI);
    double sum = 0.0;
    int i;

    for (i = -r; i <= r; i++)
    {
        double x = (double)i, v;
        v = 0.5 * (erf(w1 * (x + 1)) * (x + 1) + erf(w1 * (x - 1)) * (x - 1) -
                   2 * erf(w1 * x) * x) +
            w3 * (exp(w2 * (x + 1) * (x + 1)) + exp(w2 * (x - 1) * (x - 1)) -
                  2 * exp(w2 * x * x));
        if (v < 0) v = 0;
        k[i + r] = v;
        sum += v;
    }
    for (i = 0; i < 2 * r + 1; i++) k[i] /= sum;
}

typedef struct
{
    double *vol;
    int d[3];
    const double *fx, *fy;
    int nx, ny, xoff, yoff;
} convxy_ctx;

/* convxy of spm_conv_vol.c on whole planes: x then y, zero padded, NaN -> 0
 * in the x pass. */
static void
convxy_range(void *p, int z0, int z1)
{
    convxy_ctx *c = (convxy_ctx *)p;
    int xdim = c->d[0], ydim = c->d[1], x, y, k, z;
    double *buff = (double *)malloc(sizeof(double) * ((xdim > ydim) ? xdim : ydim));

    for (z = z0; z < z1; z++)
    {
        double *out = c->vol + (size_t)z * xdim * ydim;
        for (y = 0; y < ydim; y++)
        {
            for (x = 0; x < xdim; x++)
            {
                buff[x] = out[x + y * xdim];
                if (!isfinite(buff[x])) buff[x] = 0.0;
            }
            for (x = 0; x < xdim; x++)
            {
                double sum1 = 0.0;
                int fstart = ((x - c->xoff >= xdim) ? x - xdim - c->xoff + 1 : 0);
                int fend = ((x - (c->xoff + c->nx) < 0) ? x - c->xoff + 1 : c->nx);
                for (k = fstart; k < fend; k++)
                    sum1 += buff[x - c->xoff - k] * c->fx[k];
                out[x + y * xdim] = sum1;
            }
        }
        for (x = 0; x < xdim; x++)
        {
            for (y = 0; y < ydim; y++) buff[y] = out[x + y * xdim];
            for (y = 0; y < ydim; y++)
            {
                double sum1 = 0.0;
                int fstart = ((y - c->yoff >= ydim) ? y - ydim - c->yoff + 1 : 0);
                int fend = ((y - (c->yoff + c->ny) < 0) ? y - c->yoff + 1 : c->ny);
                for (k = fstart; k < fend; k++)
                    sum1 += buff[y - c->yoff - k] * c->fy[k];
                out[y * xdim + x] = sum1;
            }
        }
    }
    free(buff);
}

typedef struct
{
    const double *in;
    double *outd;
    float *outf;
    int d[3];
    const double *fz;
    int nz, zoff;
} convz_ctx;

/* The z pass of convxyz() in spm_conv_vol.c: over the planes inside the
 * volume, renormalised by the sum of the filter taps used. */
static void
convz_range(void *p, int zo0, int zo1)
{
    convz_ctx *c = (convz_ctx *)p;
    size_t nxy = (size_t)c->d[0] * c->d[1], xy;
    int zdim = c->d[2], zo, k;

    for (zo = zo0; zo < zo1; zo++)
    {
        int z = zo + c->nz + c->zoff - 1; /* the loop index of convxyz */
        int fstart = ((z >= zdim) ? z - zdim + 1 : 0);
        int fend = ((z - c->nz < 0) ? z + 1 : c->nz);
        double sum2 = 0.0;
        for (k = fstart; k < fend; k++) sum2 += c->fz[k];
        for (xy = 0; xy < nxy; xy++)
        {
            double sum1 = 0.0, v;
            for (k = fstart; k < fend; k++)
                sum1 += c->fz[k] * c->in[(size_t)(z - k) * nxy + xy];
            v = (sum2 != 0.0) ? sum1 / sum2 : 0.0;
            if (c->outf) c->outf[(size_t)zo * nxy + xy] = (float)v;
            else c->outd[(size_t)zo * nxy + xy] = v;
        }
    }
}

/* spm_smooth(V, V, fwhm) on a numeric array (voxel units), in place.  The
 * volume is given as double (vol_d) or float (vol_f). */
static void
spm_smooth_vol(double *vol_d, float *vol_f, const int d[3], const double fwhm[3],
               int nthreads)
{
    int r[3], i;
    double *k[3], *buf;
    size_t n = (size_t)d[0] * d[1] * d[2], j;
    convxy_ctx cxy;
    convz_ctx cz;

    for (i = 0; i < 3; i++)
    {
        double s1 = fwhm[i] / sqrt(8.0 * log(2.0));
        r[i] = (int)round(6.0 * s1);
        k[i] = (double *)malloc(sizeof(double) * (2 * r[i] + 1));
        smoothkern(fwhm[i], r[i], k[i]);
    }
    buf = (double *)malloc(sizeof(double) * n);
    for (j = 0; j < n; j++) buf[j] = vol_f ? (double)vol_f[j] : vol_d[j];

    cxy.vol = buf;
    memcpy(cxy.d, d, sizeof(int) * 3);
    cxy.fx = k[0];
    cxy.fy = k[1];
    cxy.nx = 2 * r[0] + 1;
    cxy.ny = 2 * r[1] + 1;
    cxy.xoff = -r[0];
    cxy.yoff = -r[1];
    parallel_for(d[2], nthreads, convxy_range, &cxy);

    cz.in = buf;
    cz.outd = vol_d;
    cz.outf = vol_f;
    memcpy(cz.d, d, sizeof(int) * 3);
    cz.fz = k[2];
    cz.nz = 2 * r[2] + 1;
    cz.zoff = -r[2];
    parallel_for(d[2], nthreads, convz_range, &cz);

    free(buf);
    for (i = 0; i < 3; i++) free(k[i]);
}

/* ======================================================================== */
/* spm_slice_vol(V, A, dim, [1 NaN]): SLICE_1 of spm_vol_utils.c             */
/* ======================================================================== */

/* Trilinear sampling of one plane.  mat is the 4x4 matrix in MATLAB's column
 * order (mat[0..3] is the first column); the plane's points are mat*[x y 0 1]
 * with 1-based x, y, accumulated incrementally as in SLICE_1. */
static void
slice1(const double *vol, const int d2[3], const double mat[16], int xdim1,
       int ydim1, double *image)
{
    double x2, y2, z2, s2, dx3 = mat[0], dy3 = mat[1], dz3 = mat[2], ds3 = mat[3];
    int xdim2 = d2[0], ydim2 = d2[1], zdim2 = d2[2];
    int t = 0;
    double yy;

    x2 = mat[12] + 0 * mat[8];
    y2 = mat[13] + 0 * mat[9];
    z2 = mat[14] + 0 * mat[10];
    s2 = mat[15] + 0 * mat[11];

    for (yy = 1; yy <= ydim1; yy++)
    {
        double xx;
        double x3 = x2 + yy * mat[4];
        double y3 = y2 + yy * mat[5];
        double z3 = z2 + yy * mat[6];
        double s3 = s2 + yy * mat[7];
        for (xx = 1; xx <= xdim1; xx++)
        {
            double x4, y4, z4;
            s3 += ds3;
            x4 = (x3 += dx3) / s3 - 1.0;
            y4 = (y3 += dy3) / s3 - 1.0;
            z4 = (z3 += dz3) / s3 - 1.0;

            if (z4 >= -SLICE_TINY && z4 < zdim2 + SLICE_TINY - 1 &&
                y4 >= -SLICE_TINY && y4 < ydim2 + SLICE_TINY - 1 &&
                x4 >= -SLICE_TINY && x4 < xdim2 + SLICE_TINY - 1)
            {
                double k111, k112, k121, k122, k211, k212, k221, k222;
                double dx1, dx2, dy1, dy2, dz1, dz2;
                int off1, off2, offx, offy, offz, ix4, iy4, iz4;
                const double *p0, *p1;

                ix4 = (int)floor(x4);
                dx1 = x4 - ix4;
                dx2 = 1.0 - dx1;
                iy4 = (int)floor(y4);
                dy1 = y4 - iy4;
                dy2 = 1.0 - dy1;
                iz4 = (int)floor(z4);
                dz1 = z4 - iz4;
                dz2 = 1.0 - dz1;

                ix4 = (ix4 < 0) ? ((offx = 0), 0) : ((ix4 >= xdim2 - 1) ? ((offx = 0), xdim2 - 1) : ((offx = 1), ix4));
                iy4 = (iy4 < 0) ? ((offy = 0), 0) : ((iy4 >= ydim2 - 1) ? ((offy = 0), ydim2 - 1) : ((offy = xdim2), iy4));
                iz4 = (iz4 < 0) ? ((offz = 0), 0) : ((iz4 >= zdim2 - 1) ? ((offz = 0), zdim2 - 1) : ((offz = 1), iz4));

                off1 = ix4 + xdim2 * iy4;
                off2 = off1 + offy;
                p0 = vol + (size_t)iz4 * xdim2 * ydim2;
                p1 = vol + (size_t)(iz4 + offz) * xdim2 * ydim2;
                k222 = p0[off1];
                k122 = p0[off1 + offx];
                k212 = p0[off2];
                k112 = p0[off2 + offx];
                k221 = p1[off1];
                k121 = p1[off1 + offx];
                k211 = p1[off2];
                k111 = p1[off2 + offx];

                image[t] = (((k222 * dx2 + k122 * dx1) * dy2 + (k212 * dx2 + k112 * dx1) * dy1)) * dz2 +
                           (((k221 * dx2 + k121 * dx1) * dy2 + (k211 * dx2 + k111 * dx1) * dy1)) * dz1;
            }
            else
                image[t] = NAN;
            t++;
        }
    }
}

typedef struct
{
    const double *vol;
    int d2[3], d1[3];
    const double *M; /* row-major 4x4: output voxel -> input voxel */
    float *out;
} sample_ctx;

/* out(:,:,i) = single(spm_slice_vol(vol, M*spm_matrix([0 0 i]), d1(1:2),
 * [1,NaN])) for i = 1..d1(3). */
static void
sample_range(void *p, int z0, int z1)
{
    sample_ctx *c = (sample_ctx *)p;
    size_t nxy = (size_t)c->d1[0] * c->d1[1], j;
    double A[16], *plane = (double *)malloc(sizeof(double) * nxy);
    int z, r;

    for (z = z0; z < z1; z++)
    {
        /* A = M*T(0,0,i) in MATLAB's column order */
        for (r = 0; r < 4; r++)
        {
            A[0 + r] = c->M[r * 4 + 0];
            A[4 + r] = c->M[r * 4 + 1];
            A[8 + r] = c->M[r * 4 + 2];
            A[12 + r] = c->M[r * 4 + 0] * 0.0 + c->M[r * 4 + 1] * 0.0 +
                        c->M[r * 4 + 2] * (double)(z + 1) + c->M[r * 4 + 3];
        }
        slice1(c->vol, c->d2, A, c->d1[0], c->d1[1], plane);
        for (j = 0; j < nxy; j++)
            c->out[(size_t)z * nxy + j] = (float)plane[j];
    }
    free(plane);
}

static void
sample_volume(const double *vol, const int d2[3], const double *M,
              const int d1[3], float *out, int nthreads)
{
    sample_ctx c;
    c.vol = vol;
    memcpy(c.d2, d2, sizeof(int) * 3);
    memcpy(c.d1, d1, sizeof(int) * 3);
    c.M = M;
    c.out = out;
    parallel_for(d1[2], nthreads, sample_range, &c);
}

/* ======================================================================== */
/* Quadratic B-spline coefficients with wrapped boundaries (spm_bsplinc)     */
/* ======================================================================== */

/* spm_bsplinc(V, [2 2 2 1 1 1]): quadratic coefficients, periodic. */
static int
bsplinc_wrap2(const float *vol, const int d[3], double *c)
{
    mwSize dm[3];
    int deg[3] = {2, 2, 2}, wrapflags[3] = {1, 1, 1};
    dm[0] = d[0];
    dm[1] = d[1];
    dm[2] = d[2];
    return spm_bsplinc3(vol, dm, deg, wrapflags, c);
}

/* ======================================================================== */
/* Green's operator of the regulariser (spm_shoot_greens) with FFTW          */
/* ======================================================================== */

typedef struct
{
    int d[3];
    double prm[8];
    int valid, scalar;
    size_t n, nc;           /* real voxels, complex half-spectrum entries */
    float *F;               /* scalar: nc; else 9*nc, F[(j*3+i)*nc + w] */
    double *rbuf[3];
    fftw_complex *cbuf[3], *abuf[3];
    fftw_plan fwd, inv;
} Greens;

static void
greens_free(Greens *G)
{
    int i;
    if (G->F) free(G->F);
    for (i = 0; i < 3; i++)
    {
        if (G->rbuf[i]) fftw_free(G->rbuf[i]);
        if (G->cbuf[i]) fftw_free(G->cbuf[i]);
        if (G->abuf[i]) fftw_free(G->abuf[i]);
    }
    if (G->fwd) fftw_destroy_plan(G->fwd);
    if (G->inv) fftw_destroy_plan(G->inv);
    memset(G, 0, sizeof(Greens));
}

/* Allocate buffers and plans for grid d (plans are made in this thread). */
static int
greens_grid(Greens *G, const int d[3])
{
    int i;
    if (G->rbuf[0] && G->d[0] == d[0] && G->d[1] == d[1] && G->d[2] == d[2])
        return 0;
    greens_free(G);
    memcpy(G->d, d, sizeof(int) * 3);
    G->n = (size_t)d[0] * d[1] * d[2];
    G->nc = (size_t)d[2] * d[1] * (d[0] / 2 + 1);
    for (i = 0; i < 3; i++)
    {
        G->rbuf[i] = (double *)fftw_malloc(sizeof(double) * G->n);
        G->cbuf[i] = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * G->nc);
        G->abuf[i] = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * G->nc);
        if (!G->rbuf[i] || !G->cbuf[i] || !G->abuf[i]) return 2;
    }
    G->F = (float *)malloc(sizeof(float) * 9 * G->nc);
    if (!G->F) return 2;
    /* FFTW is row-major: the slowest dimension (z) first */
    G->fwd = fftw_plan_dft_r2c_3d(d[2], d[1], d[0], G->rbuf[0], G->cbuf[0], FFTW_ESTIMATE);
    G->inv = fftw_plan_dft_c2r_3d(d[2], d[1], d[0], G->abuf[0], G->rbuf[0], FFTW_ESTIMATE);
    return (G->fwd && G->inv) ? 0 : 2;
}

typedef struct
{
    Greens *G;
} fft_ctx;

typedef struct
{
    Greens *G;
    const float *K; /* the kernel, 9 blocks of n */
    float *R;       /* its real spectrum, 9 blocks of nc */
    int pair0, pair1, buf;
} kfft_args;

/* The upper triangle of the symmetric 3x3 kernel, pairs in this order. */
static const int kpair[6][2] = {{0, 0}, {1, 0}, {2, 0}, {1, 1}, {2, 1}, {2, 2}};

static void *
kfft_worker(void *p)
{
    kfft_args *a = (kfft_args *)p;
    Greens *G = a->G;
    int q;
    size_t w;
    for (q = a->pair0; q < a->pair1; q++)
    {
        int i = kpair[q][0], j = kpair[q][1];
        const float *Kij = a->K + (size_t)(i + 3 * j) * G->n;
        for (w = 0; w < G->n; w++) G->rbuf[a->buf][w] = (double)Kij[w];
        fftw_execute_dft_r2c(G->fwd, G->rbuf[a->buf], G->cbuf[a->buf]);
        for (w = 0; w < G->nc; w++)
        {
            float v = (float)G->cbuf[a->buf][w][0];
            a->R[(size_t)(i + 3 * j) * G->nc + w] = v;
            a->R[(size_t)(j + 3 * i) * G->nc + w] = v;
        }
    }
    return NULL;
}

typedef struct
{
    Greens *G;
    const float *R;
} kinv_ctx;

/* Per frequency, the inverse of the 3x3 symbol by cofactors, as
 * spm_shoot_greens.m does it (in double, from the single symbol). */
static void
kinv_range(void *p, int w0, int w1)
{
    kinv_ctx *c = (kinv_ctx *)p;
    Greens *G = c->G;
    size_t nc = G->nc, w;
    int i, j;
    for (w = (size_t)w0; w < (size_t)w1; w++)
    {
        double A[3][3], dt;
        for (i = 0; i < 3; i++)
            for (j = 0; j < 3; j++)
                A[i][j] = (double)c->R[(size_t)(i + 3 * j) * nc + w];
        dt = A[0][0] * (A[1][1] * A[2][2] - A[1][2] * A[2][1]) +
             A[0][1] * (A[1][2] * A[2][0] - A[1][0] * A[2][2]) +
             A[0][2] * (A[1][0] * A[2][1] - A[1][1] * A[2][0]);
        dt = 1.0 / dt;
        /* F(:,:,k,r,c) as in spm_shoot_greens.m, stored at (r + 3*c) */
        G->F[(size_t)(0 + 3 * 0) * nc + w] = (float)((A[1][1] * A[2][2] - A[1][2] * A[2][1]) * dt);
        G->F[(size_t)(1 + 3 * 0) * nc + w] = (float)((A[1][2] * A[2][0] - A[1][0] * A[2][2]) * dt);
        G->F[(size_t)(2 + 3 * 0) * nc + w] = (float)((A[1][0] * A[2][1] - A[1][1] * A[2][0]) * dt);
        G->F[(size_t)(0 + 3 * 1) * nc + w] = (float)((A[0][2] * A[2][1] - A[0][1] * A[2][2]) * dt);
        G->F[(size_t)(1 + 3 * 1) * nc + w] = (float)((A[0][0] * A[2][2] - A[0][2] * A[2][0]) * dt);
        G->F[(size_t)(2 + 3 * 1) * nc + w] = (float)((A[0][1] * A[2][0] - A[0][0] * A[2][1]) * dt);
        G->F[(size_t)(0 + 3 * 2) * nc + w] = (float)((A[0][1] * A[1][2] - A[0][2] * A[1][1]) * dt);
        G->F[(size_t)(1 + 3 * 2) * nc + w] = (float)((A[0][2] * A[1][0] - A[0][0] * A[1][2]) * dt);
        G->F[(size_t)(2 + 3 * 2) * nc + w] = (float)((A[0][0] * A[1][1] - A[0][1] * A[1][0]) * dt);
    }
}

/* spm_shoot_greens('kernel', d, prm): the Fourier symbol of the regulariser,
 * inverted per frequency.  As in MATLAB the forward transform is taken in
 * double, rounded to single, and inverted in double from those singles. */
static int
greens_kernel(Greens *G, const int d[3], const double prm[8], int nthreads)
{
    double param[8];
    mwSize dm[5];
    float *K;
    size_t n, nc, w;
    int s, ret;

    ret = greens_grid(G, d);
    if (ret) return ret;
    if (G->valid && memcmp(G->prm, prm, sizeof(double) * 8) == 0) return 0;

    spm_shoot_param(prm, param);
    n = G->n;
    nc = G->nc;
    dm[0] = d[0];
    dm[1] = d[1];
    dm[2] = d[2];
    G->scalar = (param[6] == 0 && param[7] == 0);
    /* kernel() only writes the stencil; the gateway's array starts at 0 */
    K = (float *)calloc(n * (G->scalar ? 1 : 9), sizeof(float));
    if (!K) return 2;
    dm[3] = dm[4] = 3;
    spm_set_bound(0);
    spm_kernel(dm, param, K);

    if (G->scalar)
    {
        for (w = 0; w < n; w++) G->rbuf[0][w] = (double)K[w];
        fftw_execute_dft_r2c(G->fwd, G->rbuf[0], G->cbuf[0]);
        for (w = 0; w < nc; w++)
        {
            float v = (float)G->cbuf[0][w][0];
            if (v < 0) v = 0;
            G->F[w] = 1.0f / v;
        }
        if (prm[3] == 0) G->F[0] = 0;
    }
    else
    {
        float *R = (float *)malloc(sizeof(float) * 9 * nc);
        if (!R)
        {
            free(K);
            return 2;
        }
        /* MATLAB's F(:,:,:,i,j) is K[(i + 3*j)*n]; the symbol is symmetric:
         * the six transforms of its upper triangle, up to three at a time */
        {
            kfft_args args[3];
            kinv_ctx kc;
            int nt = (nthreads < 3) ? ((nthreads < 1) ? 1 : nthreads) : 3, t;
            for (t = 0; t < nt; t++)
            {
                args[t].G = G;
                args[t].K = K;
                args[t].R = R;
                args[t].pair0 = (t * 6) / nt;
                args[t].pair1 = ((t + 1) * 6) / nt;
                args[t].buf = t;
            }
            if (nt > 1) cat_parallel_run(nt, kfft_worker, args, sizeof(kfft_args));
            else kfft_worker(&args[0]);
            kc.G = G;
            kc.R = R;
            parallel_for((int)nc, nthreads, kinv_range, &kc);
        }
        if (prm[3] == 0)
            for (s = 0; s < 9; s++) G->F[(size_t)s * nc] = 0;
        free(R);
    }
    free(K);
    memcpy(G->prm, prm, sizeof(double) * 8);
    G->valid = 1;
    return 0;
}

typedef struct
{
    Greens *G;
    const float *m;
    float *v;
    double prm0[3];
} apply_ctx;

static void
greens_fwd_range(void *p, int i0, int i1)
{
    apply_ctx *c = (apply_ctx *)p;
    Greens *G = c->G;
    int i;
    size_t w;
    for (i = i0; i < i1; i++)
    {
        const float *mi = c->m + (size_t)i * G->n;
        for (w = 0; w < G->n; w++) G->rbuf[i][w] = (double)mi[w];
        fftw_execute_dft_r2c(G->fwd, G->rbuf[i], G->cbuf[i]);
    }
}

static void
greens_inv_range(void *p, int j0, int j1)
{
    apply_ctx *c = (apply_ctx *)p;
    Greens *G = c->G;
    double scale = 1.0 / (double)G->n;
    int i, j;
    size_t w;
    for (j = j0; j < j1; j++)
    {
        fftw_complex *a = G->abuf[j];
        float *vj = c->v + (size_t)j * G->n;
        if (G->scalar)
        {
            double s = 1.0 / (c->prm0[j] * c->prm0[j]);
            for (w = 0; w < G->nc; w++)
            {
                a[w][0] = G->F[w] * G->cbuf[j][w][0] * s;
                a[w][1] = G->F[w] * G->cbuf[j][w][1] * s;
            }
        }
        else
        {
            for (w = 0; w < G->nc; w++)
            {
                double re = 0.0, im = 0.0;
                for (i = 0; i < 3; i++)
                {
                    double f = G->F[(size_t)(j + 3 * i) * G->nc + w];
                    re += f * G->cbuf[i][w][0];
                    im += f * G->cbuf[i][w][1];
                }
                a[w][0] = re;
                a[w][1] = im;
            }
        }
        fftw_execute_dft_c2r(G->inv, a, G->rbuf[j]);
        for (w = 0; w < G->n; w++) vj[w] = (float)(G->rbuf[j][w] * scale);
    }
}

/* v = spm_shoot_greens(m, F, prm): three forward and three inverse
 * transforms, one per component and thread. */
static void
greens_apply(Greens *G, const float *m, float *v, const double prm[8],
             int nthreads)
{
    apply_ctx c;
    c.G = G;
    c.m = m;
    c.v = v;
    c.prm0[0] = prm[0];
    c.prm0[1] = prm[1];
    c.prm0[2] = prm[2];
    parallel_for(3, (nthreads < 3) ? nthreads : 3, greens_fwd_range, &c);
    parallel_for(3, (nthreads < 3) ? nthreads : 3, greens_inv_range, &c);
}

/* ======================================================================== */
/* Shooting state                                                           */
/* ======================================================================== */

typedef struct
{
    int nthreads;
    Greens G;
    /* scratch of spm_shoot3d, sized for the finest grid */
    float *m0, *mt, *vt, *dp, *dJ, *phi2, *Jphi2;
} Shoot;

static void
shoot_free(Shoot *S)
{
    greens_free(&S->G);
    free(S->m0);
    free(S->mt);
    free(S->vt);
    free(S->dp);
    free(S->dJ);
    free(S->phi2);
    free(S->Jphi2);
}

static int
shoot_alloc(Shoot *S, size_t n)
{
    S->m0 = (float *)malloc(sizeof(float) * 3 * n);
    S->mt = (float *)malloc(sizeof(float) * 3 * n);
    S->vt = (float *)malloc(sizeof(float) * 3 * n);
    S->dp = (float *)malloc(sizeof(float) * 3 * n);
    S->dJ = (float *)malloc(sizeof(float) * 9 * n);
    S->phi2 = (float *)malloc(sizeof(float) * 3 * n);
    S->Jphi2 = (float *)malloc(sizeof(float) * 9 * n);
    return (S->m0 && S->mt && S->vt && S->dp && S->dJ && S->phi2 && S->Jphi2) ? 0 : 2;
}

/* spm_diffeo('vel2mom', v, prm) */
static void
gw_vel2mom(const int d[3], float *v, const double prm[8], float *m)
{
    double param[8];
    mwSize dm[3];
    dm[0] = d[0];
    dm[1] = d[1];
    dm[2] = d[2];
    spm_shoot_param(prm, param);
    spm_vel2mom(dm, v, param, m);
}

/* [phi, Jphi] = spm_shoot3d(v0, prm, [N ...]) */
static int
shoot3d(Shoot *S, const int d[3], float *v0, const double prm[8], int N,
        float *phi, float *Jphi)
{
    mwSize dm[4];
    size_t n = (size_t)d[0] * d[1] * d[2];
    float sc;
    int t, ret;

    spm_set_bound(0);
    dm[0] = d[0];
    dm[1] = d[1];
    dm[2] = d[2];
    dm[3] = 3;
    if (N < 1) N = 1;
    ret = greens_kernel(&S->G, d, prm, S->nthreads);
    if (ret) return ret;
    gw_vel2mom(d, v0, prm, S->m0);
    sc = (float)(1.0 / N);
    spm_smalldef_jac1(dm, sc, v0, phi, Jphi);
    for (t = 2; t <= N; t++)
    {
        /* pushc_grads adds into its output; the gateway hands it zeros */
        memset(S->mt, 0, sizeof(float) * 3 * n);
        spm_pushc_grads(dm, dm, phi, Jphi, S->m0, S->mt);
        greens_apply(&S->G, S->mt, S->vt, prm, S->nthreads);
        spm_smalldef_jac1(dm, sc, S->vt, S->dp, S->dJ);
        spm_composition_jacobian(dm, (mwSize)n, S->dp, S->dJ, phi, Jphi, S->phi2, S->Jphi2);
        memcpy(phi, S->phi2, sizeof(float) * 3 * n);
        memcpy(Jphi, S->Jphi2, sizeof(float) * 9 * n);
    }
    return 0;
}

/* dt = spm_diffeo('det', J) */
static void
jac_det(const int d[3], float *J, float *dt)
{
    mwSize dm[3];
    dm[0] = d[0];
    dm[1] = d[1];
    dm[2] = d[2];
    spm_determinant(dm, J, dt);
}

/* ======================================================================== */
/* The template of one level and its derivatives (mnom_derivs)              */
/* ======================================================================== */

typedef struct
{
    int d[3];
    size_t n;
    int ng;           /* number of template classes incl. background */
    float *logmu;     /* ng*n: single(log(mu)) */
    float *mu;        /* ng*n: single(mu) */
    float *D;         /* 3*ng*n: single(d1|d2|d3) of each class */
    float *W;         /* ng*ng*n: single(k==k1 ? max(mu.*(1-mu1),0) : -mu.*mu1) */
} TplDeriv;

static void
tplderiv_free(TplDeriv *T)
{
    free(T->logmu);
    free(T->mu);
    free(T->D);
    free(T->W);
    memset(T, 0, sizeof(TplDeriv));
}

/* The B-spline part of mnom_derivs, which depends only on the template: the
 * MATLAB evaluates it at every iteration, here once per level.  spm_bsplins
 * of the three-slice window at its middle slice equals spm_bsplins of the
 * whole volume at that slice with periodic boundaries.  Slice by slice, so
 * that the B-spline calls (threaded inside) never run concurrently. */
static int
tplderiv_build(TplDeriv *T, const int d[3], int ng, double *const *coef)
{
    size_t nxy = (size_t)d[0] * d[1], i;
    mwSize dm[3];
    int deg[3] = {2, 2, 2}, wrapflags[3] = {1, 1, 1}, x, y, z, k, k1;
    double *x0, *x1, *x2, *g, *d0, *d1, *d2;

    tplderiv_free(T);
    memcpy(T->d, d, sizeof(int) * 3);
    T->n = (size_t)d[0] * d[1] * d[2];
    T->ng = ng;
    T->logmu = (float *)malloc(sizeof(float) * ng * T->n);
    T->mu = (float *)malloc(sizeof(float) * ng * T->n);
    T->D = (float *)malloc(sizeof(float) * 3 * ng * T->n);
    T->W = (float *)malloc(sizeof(float) * ng * ng * T->n);
    x0 = (double *)malloc(sizeof(double) * 3 * nxy);
    g = (double *)malloc(sizeof(double) * 4 * ng * nxy);
    if (!T->logmu || !T->mu || !T->D || !T->W || !x0 || !g)
    {
        free(x0);
        free(g);
        return 2;
    }
    x1 = x0 + nxy;
    x2 = x1 + nxy;
    d0 = g + ng * nxy;
    d1 = d0 + ng * nxy;
    d2 = d1 + ng * nxy;
    dm[0] = d[0];
    dm[1] = d[1];
    dm[2] = d[2];
    for (y = 0; y < d[1]; y++)
        for (x = 0; x < d[0]; x++)
        {
            x0[x + (size_t)d[0] * y] = x + 1;
            x1[x + (size_t)d[0] * y] = y + 1;
        }
    for (z = 0; z < d[2]; z++)
    {
        for (i = 0; i < nxy; i++) x2[i] = z + 1;
        for (k = 0; k < ng; k++)
            spm_bsplins3(coef[k], dm, nxy, x0, x1, x2, deg, wrapflags, g + k * nxy,
                         d0 + k * nxy, d1 + k * nxy, d2 + k * nxy);
        for (i = 0; i < nxy; i++)
        {
            size_t v = (size_t)z * nxy + i;
            double mu[MAXCLS], s = 0.0;
            for (k = 0; k < ng; k++) mu[k] = exp(g[k * nxy + i]);
            for (k = 0; k < ng; k++) s = s + mu[k];
            for (k = 0; k < ng; k++) mu[k] = mu[k] / s;
            for (k = 0; k < ng; k++)
            {
                T->logmu[k * T->n + v] = (float)log(mu[k]);
                T->mu[k * T->n + v] = (float)mu[k];
                T->D[(3 * k + 0) * T->n + v] = (float)d0[k * nxy + i];
                T->D[(3 * k + 1) * T->n + v] = (float)d1[k * nxy + i];
                T->D[(3 * k + 2) * T->n + v] = (float)d2[k * nxy + i];
                for (k1 = 0; k1 < ng; k1++)
                {
                    double w;
                    if (k1 != k) w = -mu[k] * mu[k1];
                    else
                    {
                        w = mu[k] * (1 - mu[k1]);
                        if (w < 0) w = 0;
                    }
                    T->W[(k * ng + k1) * T->n + v] = (float)w;
                }
            }
        }
    }
    free(x0);
    free(g);
    return 0;
}

typedef struct
{
    const TplDeriv *T;
    const float *f1;  /* ng*n: pulled and Jacobian-weighted tissue maps */
    const float *dt;
    float *b, *A;
    double *ll_slice; /* per slice */
} derivs_ctx;

/* The per-voxel part of mnom_derivs, in single as in MATLAB. */
static void
derivs_range(void *p, int z0, int z1)
{
    derivs_ctx *c = (derivs_ctx *)p;
    const TplDeriv *T = c->T;
    size_t n = T->n, nxy = (size_t)T->d[0] * T->d[1];
    int z, k, k1, ng = T->ng;

    for (z = z0; z < z1; z++)
    {
        size_t i0 = (size_t)z * nxy, i;
        double ll = 0.0;
        for (i = i0; i < i0 + nxy; i++)
        {
            float b1 = 0.0f, b2 = 0.0f, b3 = 0.0f;
            float a1 = 0.0f, a2 = 0.0f, a3 = 0.0f, a4 = 0.0f, a5 = 0.0f, a6 = 0.0f;
            float dti = c->dt[i];
            for (k = 0; k < ng; k++)
            {
                float f1 = c->f1[k * n + i];
                float tmp = f1 * T->logmu[k * n + i];
                const float *Dk = T->D + 3 * k * n;
                ll -= tmp;
                tmp = f1 - T->mu[k * n + i] * dti;
                b1 = b1 + tmp * Dk[i];
                b2 = b2 + tmp * Dk[n + i];
                b3 = b3 + tmp * Dk[2 * n + i];
                for (k1 = 0; k1 < ng; k1++)
                {
                    const float *Dk1 = T->D + 3 * k1 * n;
                    float w = T->W[(k * ng + k1) * n + i] * dti;
                    a1 = a1 + w * Dk[i] * Dk1[i];
                    a2 = a2 + w * Dk[n + i] * Dk1[n + i];
                    a3 = a3 + w * Dk[2 * n + i] * Dk1[2 * n + i];
                    a4 = a4 + w * Dk[i] * Dk1[n + i];
                    a5 = a5 + w * Dk[i] * Dk1[2 * n + i];
                    a6 = a6 + w * Dk[n + i] * Dk1[2 * n + i];
                }
            }
            c->b[i] = b1;
            c->b[n + i] = b2;
            c->b[2 * n + i] = b3;
            c->A[i] = a1;
            c->A[n + i] = a2;
            c->A[2 * n + i] = a3;
            c->A[3 * n + i] = a4;
            c->A[4 * n + i] = a5;
            c->A[5 * n + i] = a6;
        }
        c->ll_slice[z] = ll;
    }
}

/* [u, ll1, ll2] = spm_shoot_update(g, f, u, y, dt, prm, bs_args, scale) */
static int
shoot_update(Shoot *S, const TplDeriv *T, const float *f, float *u,
             const float *y, const float *dt, const double prm[8],
             double scale, double *ll1, double *ll2, float *work)
{
    size_t n = T->n, i;
    float *f1 = work;                 /* ng*n */
    float *b = f1 + (size_t)T->ng * n; /* 3n */
    float *A = b + 3 * n;              /* 6n */
    float *m0 = A + 6 * n;             /* 3n */
    float *du = m0 + 3 * n;            /* 3n */
    double *lls, ll = 0.0, l2 = 0.0;
    derivs_ctx dc;
    mwSize dm0[4];
    int z, k;
    float fs = (float)scale;

    /* f1 = spm_diffeo('pullc', f, y) .* dt for every class */
    spm_set_bound(0);
    dm0[0] = T->d[0];
    dm0[1] = T->d[1];
    dm0[2] = T->d[2];
    dm0[3] = T->ng;
    spm_pullc(dm0, (mwSize)n, (mwSize)T->ng, (float *)y, (float *)f, f1);
    for (k = 0; k < T->ng; k++)
        for (i = 0; i < n; i++) f1[k * n + i] = f1[k * n + i] * dt[i];

    lls = (double *)calloc(T->d[2], sizeof(double));
    dc.T = T;
    dc.f1 = f1;
    dc.dt = dt;
    dc.b = b;
    dc.A = A;
    dc.ll_slice = lls;
    parallel_for(T->d[2], S->nthreads, derivs_range, &dc);
    for (z = 0; z < T->d[2]; z++) ll += lls[z];
    free(lls);

    gw_vel2mom(T->d, u, prm, m0);
    for (i = 0; i < 3 * n; i++) l2 += (double)(m0[i] * u[i]);
    for (i = 0; i < 3 * n; i++) b[i] = b[i] + m0[i];
    {
        mwSize dm[3];
        dm[0] = T->d[0];
        dm[1] = T->d[1];
        dm[2] = T->d[2];
        if (spm_shoot_fmg(dm, A, b, prm, 3, 2, NULL, du)) return 2;
    }
    for (i = 0; i < 3 * n; i++) u[i] = u[i] - fs * du[i];
    *ll1 = ll;
    *ll2 = 0.5 * l2;
    return 0;
}

/* ======================================================================== */
/* Levels                                                                   */
/* ======================================================================== */

typedef struct
{
    int rd[3];        /* registration grid */
    double res;       /* tempres2: registration resolution in mm */
    double rred;      /* rred(ti) (first element, as MATLAB indexes it) */
    double Mads[16];  /* registration voxel -> template voxel (1-based) */
    double Mrregs[16]; /* registration voxel -> input voxel (1-based) */
    double Mys[16];   /* this level's voxel -> previous level's voxel */
} Level;

static void
voxel_corners(const int d[3], double C[4][8])
{
    int c;
    for (c = 0; c < 8; c++)
    {
        C[0][c] = (c & 1) ? d[0] : 1;
        C[1][c] = ((c >> 1) & 1) ? d[1] : 1;
        C[2][c] = ((c >> 2) & 1) ? d[2] : 1;
        C[3][c] = 1;
    }
}

/* ======================================================================== */
/* CAT's safeguards                                                         */
/* ======================================================================== */

/* extdt = max([dt(:)', 1/max(eps,dt(isfinite(dt(:))))]) as written in
 * cat_main_registration.m.  The 1/... is MATLAB's matrix division of a
 * scalar by a column vector, the minimum-norm solution v'/(v'*v), not an
 * element-wise reciprocal -- so in effect the test only sees max(dt).  The
 * quirk is kept: it decides when CAT filters the velocity. */
static double
ext_det(const float *dt, size_t n, int *nonfinite)
{
    double mx = -HUGE_VAL, ss = 0.0, mxv = -HUGE_VAL;
    size_t i;
    int any = 0, nf = 0;

    for (i = 0; i < n; i++)
    {
        double v = dt[i];
        if (isnan(v))
        {
            nf = 1;
            continue;
        }
        if (!isfinite(v)) nf = 1;
        if (v > mx) mx = v;
        if (isfinite(v))
        {
            double e = (v > MATLAB_EPS) ? v : MATLAB_EPS;
            ss += e * e;
            any = 1;
        }
    }
    if (any)
        for (i = 0; i < n; i++)
        {
            double v = dt[i];
            if (isfinite(v))
            {
                double e = (v > MATLAB_EPS) ? v : MATLAB_EPS;
                if (e / ss > mxv) mxv = e / ss;
            }
        }
    *nonfinite = nf;
    return (mxv > mx) ? mxv : mx;
}

static void
smooth3X_copy(const float *in, float *out, const int d[3], double s)
{
    memcpy(out, in, sizeof(float) * (size_t)d[0] * d[1] * d[2]);
    cat12_vol_smooth3X(out, d, s);
}

/* The hot-spot filter of run_Shooting (RD202101), first pass: local outliers
 * of the Jacobian are detected relative to its smoothed copy and the
 * velocity is filtered there. */
static void
hotspot_first(float *u, const float *dt, const int d[3], double vr,
              double sched_it, float *t1, float *t2, unsigned char *msk)
{
    size_t n = (size_t)d[0] * d[1] * d[2], i;
    int c;

    smooth3X_copy(dt, t1, d, vr);
    for (i = 0; i < n; i++) t2[i] = (dt[i] / t1[i] > 1.5f) ? 1.0f : 0.0f;
    cat12_vol_smooth3X(t2, d, 2 * vr); /* mxu */
    for (i = 0; i < n; i++) msk[i] = (t2[i] > 0.01f);
    for (c = 0; c < 3; c++)
    {
        float *uc = u + (size_t)c * n;
        for (i = 0; i < n; i++) uc[i] = uc[i] / 10.0f;
        cat12_vol_laplace3R(uc, msk, d, 0.1);
        for (i = 0; i < n; i++) uc[i] = uc[i] * 10.0f;
    }
    {
        double s = sched_it / 8.0;
        if (s > 2) s = 2;
        if (s < 0.2) s = 0.2;
        for (c = 0; c < 3; c++)
        {
            float *uc = u + (size_t)c * n;
            smooth3X_copy(uc, t1, d, s);
            for (i = 0; i < n; i++) uc[i] = t1[i] * t2[i] + uc[i] * (1 - t2[i]);
        }
    }
}

/* The repeated pass of the hot-spot filter (inside CAT's while loop). */
static void
hotspot_repeat(float *u, const float *dt, const int d[3], double vr, int nfilt,
               double dtl, const float *bgexp, float *t1, float *t2,
               float *t3, unsigned char *msk)
{
    size_t n = (size_t)d[0] * d[1] * d[2], i;
    int c;

    smooth3X_copy(dt, t1, d, vr);
    for (i = 0; i < n; i++) t2[i] = (dt[i] / t1[i] > 1.5f) ? 1.0f : 0.0f;
    cat12_vol_smooth3X(t2, d, 2 * vr); /* mxu = smooth(...) > 0.01 below */
    /* mxu2 = dt>(dtl*2) | dt<1/(dtl*2), plus three times |div(u(:,:,:,1))|
     * -- the MATLAB loops over i but always takes the first component */
    cat12_vol_div(u, d, t3);
    for (i = 0; i < n; i++)
    {
        float mxu2 = (dt[i] > dtl * 2 || dt[i] < 1.0 / (dtl * 2)) ? 1.0f : 0.0f;
        int m;
        mxu2 = mxu2 + fabsf(t3[i]);
        mxu2 = mxu2 + fabsf(t3[i]);
        mxu2 = mxu2 + fabsf(t3[i]);
        m = (t2[i] > 0.01f) || (mxu2 > (1.0f - bgexp[i]));
        m = m && (bgexp[i] < 0.5f);
        m = m || !isfinite(dt[i]) || dt[i] > 1000 || dt[i] < 1.0 / 1000;
        msk[i] = (unsigned char)m;
    }
    for (c = 0; c < 3; c++)
    {
        float *uc = u + (size_t)c * n;
        if (nfilt > 0) cat12_vol_smooth3X(uc, d, 1.0);
        for (i = 0; i < n; i++) uc[i] = uc[i] / 10.0f;
        cat12_vol_laplace3R(uc, msk, d, 0.02);
        for (i = 0; i < n; i++) uc[i] = uc[i] * 10.0f;
    }
}

/* ======================================================================== */
/* Main                                                                     */
/* ======================================================================== */

void
CAT_VolShootOptionsInit(CAT_VolShootOpts *opts)
{
    opts->regstr = 1.0;
    opts->nits = 64;
    opts->ll1th = 0.001;
    opts->lowres = 2.5;
    opts->reg_scale = 0.5;
    opts->hydro = 1;
    opts->hotspot = 1;
    opts->n_threads = 0;
    opts->verbose = 0;
}

static double
vec_min_finite(const float *v, size_t n)
{
    double m = HUGE_VAL;
    size_t i;
    for (i = 0; i < n; i++)
        if (isfinite(v[i]) && v[i] < m) m = v[i];
    return m;
}

int
CAT_VolShoot(const float *const cls[], int n_cls, const int dims[3],
             const double mat[16], const double affine[16],
             const float *const templates[CAT_SHOOT_NTEMPLATES],
             const int tdims[3], const double tmat[16], double reg_vx,
             const float *lesion, const CAT_VolShootOpts *opts_in,
             float *y_out, float *iy_out, float *jac_out,
             CAT_VolShootInfo *info)
{
    CAT_VolShootOpts opts;
    Level L[NLEV];
    Shoot S;
    TplDeriv T;
    double M0[16], tmpM[16], Aff[16], iM0[16], iAff[16], P[16];
    double bb[2][3], mmtmp[4][8], vxtpm[4][8], vx3rr[4][8], Q[4][8];
    double tmpres, regres, vx_vol, resfac[NLEV], ll1th;
    double *sched = NULL, rparam[5] = {1e-4, 0.001, 0.2, 0.05, 0.2};
    int *eul_its = NULL, *tmpl_no = NULL, nits, n1 = n_cls, ti = 0, it;
    int prev_ti = -1, nthreads, i, j, k, c, ret = 0;
    size_t n_in = (size_t)dims[0] * dims[1] * dims[2];
    size_t ntpl = (size_t)tdims[0] * tdims[1] * tdims[2], nmax, n = 0;
    float *f = NULL, *u = NULL, *y = NULL, *J = NULL, *dt = NULL, *work = NULL;
    float *t1 = NULL, *t2 = NULL, *t3 = NULL, *bgexp = NULL;
    unsigned char *msk = NULL;
    double *coef[MAXCLS];
    double prm[8], ll[2] = {0, 0}, llo[2];
    int int_N = 1, hydro = 0, ng = 0, old_threads = 0;
    const double scale = 0.8; /* sd.scale */

    memset(&S, 0, sizeof(S));
    memset(&T, 0, sizeof(T));
    for (k = 0; k < MAXCLS; k++) coef[k] = NULL;
    if (opts_in) opts = *opts_in;
    else CAT_VolShootOptionsInit(&opts);
    if (info) memset(info, 0, sizeof(*info));
    if (n_cls < 2 || n_cls > CAT_SHOOT_MAXCLASSES || opts.nits < 3 ||
        !cls || !templates || !y_out || !(opts.reg_scale > 0))
        return 1;

    nthreads = (opts.n_threads > 0) ? opts.n_threads : cat_parallel_nthreads(1 << 20);
    if (nthreads > MAX_SHOOT_THREADS) nthreads = MAX_SHOOT_THREADS;
    S.nthreads = nthreads;
    old_threads = spm_shoot_get_threads();
    spm_shoot_set_threads(nthreads);
    spm_shoot_clear_error();

    nii2spm(mat, M0);
    nii2spm(tmat, tmpM);
    if (affine) memcpy(Aff, affine, sizeof(Aff));
    else m4_eye(Aff);
    if (m4_inv(M0, iM0) || m4_inv(Aff, iAff)) return 1;

    /* regsetup */
    tmpres = fabs(tmpM[0]);
    regres = tmpres; /* reg.opt.rres for regstr 0.5 and 4 */
    /* CAT's vxreg is its TPM's resolution (1.5 mm) whatever the template's;
     * by default the regulariser gets the template's voxel size instead */
    if (!(reg_vx > 0)) reg_vx = tmpres;
    vx_vol = sqrt(M0[0] * M0[0] + M0[4] * M0[4] + M0[8] * M0[8]);
    if (opts.regstr == 4)
    {
        for (i = 0; i < NLEV; i++) resfac[i] = 1.0;
        ll1th = 0;
    }
    else if (opts.regstr > 0 && opts.regstr <= 1)
    {
        double stepsize = (opts.lowres - tmpres) / 4;
        for (i = 0; i < NLEV; i++) resfac[i] = (opts.lowres - i * stepsize) / tmpres;
        ll1th = 0.0010 + 0.10 * (1 - opts.regstr);
    }
    else
        return 1;

    /* tmpbb, getmm and the per-level matrices of run_Shooting */
    {
        double cmin[3] = {HUGE_VAL, HUGE_VAL, HUGE_VAL}, cmax[3] = {-HUGE_VAL, -HUGE_VAL, -HUGE_VAL};
        double C[4][8], t;
        voxel_corners(tdims, C);
        for (c = 0; c < 8; c++)
            for (i = 0; i < 3; i++)
            {
                double v = tmpM[i * 4 + 0] * C[0][c] + tmpM[i * 4 + 1] * C[1][c] +
                           tmpM[i * 4 + 2] * C[2][c] + tmpM[i * 4 + 3];
                if (v < cmin[i]) cmin[i] = v;
                if (v > cmax[i]) cmax[i] = v;
            }
        for (i = 0; i < 3; i++)
        {
            bb[0][i] = cmin[i];
            bb[1][i] = cmax[i];
        }
        if (bb[0][0] < bb[1][0]) /* if tmpbb(1)<tmpbb(2), swap them */
        {
            t = bb[0][0];
            bb[0][0] = bb[1][0];
            bb[1][0] = t;
        }
        for (c = 0; c < 8; c++)
        {
            mmtmp[0][c] = bb[c & 1][0];
            mmtmp[1][c] = bb[(c >> 1) & 1][1];
            mmtmp[2][c] = bb[(c >> 2) & 1][2];
            mmtmp[3][c] = 1.0;
        }
    }
    {
        double itmpM[16];
        if (m4_inv(tmpM, itmpM)) return 1;
        for (c = 0; c < 8; c++)
            for (i = 0; i < 4; i++)
                vxtpm[i][c] = itmpM[i * 4 + 0] * mmtmp[0][c] + itmpM[i * 4 + 1] * mmtmp[1][c] +
                              itmpM[i * 4 + 2] * mmtmp[2][c] + itmpM[i * 4 + 3] * mmtmp[3][c];
        /* M0\inv(res.Affine)*M1rr, M1rr = tmpM since regres == tmpres */
        m4_mul(iM0, iAff, P);
        m4_mul(P, tmpM, P);
        for (c = 0; c < 8; c++)
            for (i = 0; i < 4; i++)
                Q[i][c] = P[i * 4 + 0] * vxtpm[0][c] + P[i * 4 + 1] * vxtpm[1][c] +
                          P[i * 4 + 2] * vxtpm[2][c] + P[i * 4 + 3] * vxtpm[3][c];
    }
    nmax = 0;
    for (i = 0; i < NLEV; i++)
    {
        size_t nl;
        for (j = 0; j < 3; j++)
            L[i].rd[j] = (int)floor(tdims[j] * tmpres / regres / resfac[i]);
        L[i].res = resfac[i] * regres;
        L[i].rred = (regres / resfac[i]) / vx_vol;
        voxel_corners(L[i].rd, vx3rr);
        if (mrdivide48(vxtpm, vx3rr, L[i].Mads) || mrdivide48(Q, vx3rr, L[i].Mrregs))
            return 1;
        if (i == 0) m4_eye(L[i].Mys);
        else
        {
            double iprev[16];
            if (m4_inv(L[i - 1].Mads, iprev)) return 1;
            m4_mul(L[i].Mads, iprev, L[i].Mys);
        }
        nl = (size_t)L[i].rd[0] * L[i].rd[1] * L[i].rd[2];
        if (nl > nmax) nmax = nl;
        if (info) info->level_res[i] = L[i].res;
    }
    if (L[NLEV - 1].rd[0] != tdims[0] || L[NLEV - 1].rd[1] != tdims[1] ||
        L[NLEV - 1].rd[2] != tdims[2])
        return 1;

    /* the schedule of the regularisation and of the Euler steps */
    if (opts.regstr != 4)
    {
        int onits = 24, nn;
        double lam = 0.25, inter = 32, smin, smax, sch[24], step;
        nits = opts.nits;
        for (i = 0; i < onits; i++) sch[i] = (inter - 1) * exp(-lam * i) + 1;
        smin = sch[onits - 1];
        smax = sch[0];
        for (i = 0; i < onits; i++) sch[i] = (sch[i] - smin) / (smax - smin) * (inter - 1) + 1;
        /* interp1(sched, 1:(numel-1)/(nits-1):numel) */
        step = (double)(onits - 1) / (nits - 1);
        nn = (int)floor((onits - 1) / step + 1e-10) + 1;
        sched = (double *)malloc(sizeof(double) * (nn + 1));
        for (i = 0; i < nn; i++)
        {
            /* linear interpolation of sch at the 1-based position xi */
            double xi = 1 + i * step;
            int lo = (int)floor(xi);
            if (lo >= onits)
                sched[i] = sch[onits - 1];
            else
                sched[i] = sch[lo - 1] + (sch[lo] - sch[lo - 1]) * (xi - lo);
        }
        eul_its = (int *)malloc(sizeof(int) * nits);
        for (i = 0; i < nits; i++)
            eul_its[i] = (int)round(i * (8 - 0.5001) / (nits - 1) + 1);
        nits = nn - 1;
        tmpl_no = (int *)malloc(sizeof(int) * (nits + 1));
        for (i = 0; i < nits; i++)
            tmpl_no[i] = (int)floor((double)i / (nits - 1) * (NLEV - 0.51)) + 1;
    }
    else
    {
        double s_end;
        nits = 24;
        sched = (double *)malloc(sizeof(double) * (nits + 1));
        for (i = 0; i <= nits; i++) sched[i] = (32 - 1) * exp(-0.5 * i) + 1;
        s_end = sched[nits];
        for (i = 0; i <= nits; i++) sched[i] /= s_end;
        eul_its = (int *)malloc(sizeof(int) * nits);
        for (i = 0; i < nits; i++)
            eul_its[i] = (int)round(i * (8 - 0.5001) / (nits - 1) + 1);
        tmpl_no = (int *)malloc(sizeof(int) * (nits + 1));
        for (i = 0; i < nits; i++)
            tmpl_no[i] = (int)round((double)i / (nits - 1) * (NLEV - 0.51)) + 1;
    }

    /* buffers on the finest grid */
    f = (float *)malloc(sizeof(float) * MAXCLS * nmax);
    u = (float *)calloc(3 * nmax, sizeof(float));
    y = (float *)malloc(sizeof(float) * 3 * nmax);
    J = (float *)malloc(sizeof(float) * 9 * nmax);
    dt = (float *)malloc(sizeof(float) * nmax);
    work = (float *)malloc(sizeof(float) * (MAXCLS + 15) * nmax);
    t1 = (float *)malloc(sizeof(float) * nmax);
    t2 = (float *)malloc(sizeof(float) * nmax);
    t3 = (float *)malloc(sizeof(float) * nmax);
    bgexp = (float *)malloc(sizeof(float) * nmax);
    msk = (unsigned char *)malloc(nmax);
    for (k = 0; k < MAXCLS; k++) coef[k] = (double *)malloc(sizeof(double) * nmax);
    if (!f || !u || !y || !J || !dt || !work || !t1 || !t2 || !t3 || !bgexp ||
        !msk || shoot_alloc(&S, nmax))
    {
        ret = 2;
        goto done;
    }
    for (k = 0; k < MAXCLS; k++)
        if (!coef[k])
        {
            ret = 2;
            goto done;
        }

    if (opts.verbose)
    {
        if (opts.regstr != 4)
            printf("Optimized Shooting registration with %0.2f:%0.2f:%0.2f mm (regstr=%0.2f)\n",
                   L[0].res, L[1].res - L[0].res, L[NLEV - 1].res, opts.regstr);
        else
            printf("Default Shooting registration with %0.2f mm\n", regres);
    }

    it = 1;
    while (it <= nits)
    {
        int lvl_change = (it == 1 || tmpl_no[it - 1] != tmpl_no[it - 2]);
        int nfilt = 0, nonfinite; /* passes of CAT's hot-spot filter loop */
        double extdt;

        if (lvl_change)
        {
            Level *Lv;
            double *dbuf = (double *)malloc(sizeof(double) * ((n_in > ntpl) ? n_in : ntpl));
            float *fbuf = (float *)malloc(sizeof(float) * ((n_in > ntpl) ? n_in : ntpl));
            float *gbg;
            double fw[3];

            if (!dbuf || !fbuf)
            {
                free(dbuf);
                free(fbuf);
                ret = 2;
                goto done;
            }
            prev_ti = ti;
            ti = tmpl_no[it - 1] - 1;
            Lv = &L[ti];
            n = (size_t)Lv->rd[0] * Lv->rd[1] * Lv->rd[2];

            /* the tissue maps, smoothed a little on the coarse levels */
            fw[0] = fw[1] = fw[2] = Lv->rred * 0.1 * ((NLEV - (ti + 1) > 0) ? NLEV - (ti + 1) : 0);
            for (k = 0; k < n1; k++)
            {
                size_t v;
                memcpy(fbuf, cls[k], sizeof(float) * n_in);
                spm_smooth_vol(NULL, fbuf, dims, fw, nthreads);
                for (v = 0; v < n_in; v++) dbuf[v] = (double)fbuf[v];
                sample_volume(dbuf, dims, Lv->Mrregs, Lv->rd, f + (size_t)k * n, nthreads);
                for (v = 0; v < n; v++)
                    if (!isfinite(f[k * n + v])) f[k * n + v] = 0.0f;
            }
            for (i = 0; i < (int)n; i++) f[(size_t)n1 * n + i] = 1.0f;

            /* the template */
            gbg = t3;
            for (i = 0; i < (int)n; i++) gbg[i] = 1.0f;
            fw[0] = fw[1] = fw[2] = Lv->res / tmpres;
            for (k = 0; k < n1; k++)
            {
                size_t v;
                float *gk = t1;
                double mn;
                for (v = 0; v < ntpl; v++) dbuf[v] = (double)templates[ti][k * ntpl + v];
                spm_smooth_vol(dbuf, NULL, tdims, fw, nthreads);
                sample_volume(dbuf, tdims, Lv->Mads, Lv->rd, gk, nthreads);
                mn = vec_min_finite(gk, n);
                for (v = 0; v < n; v++)
                    if (isnan(gk[v])) gk[v] = (float)mn;
                for (v = 0; v < n; v++) gbg[v] = gbg[v] - gk[v];
                for (v = 0; v < n; v++) gk[v] = logf(gk[v]);
                bsplinc_wrap2(gk, Lv->rd, coef[k]);
            }
            for (i = 0; i < (int)n; i++)
                t1[i] = logf((gbg[i] > (float)MATLAB_EPS) ? gbg[i] : (float)MATLAB_EPS);
            bsplinc_wrap2(t1, Lv->rd, coef[n1]);
            ng = n1 + 1;

            /* RD202101: a very large ventricle becomes a class of its own */
            hydro = 0;
            if (opts.regstr != 4 && opts.hydro && n1 == 2)
            {
                float *Ygwm = (float *)malloc(sizeof(float) * n);
                float *Yven = (float *)malloc(sizeof(float) * n);
                unsigned char *Yb = (unsigned char *)malloc(n);
                unsigned char *Yv = (unsigned char *)malloc(n);
                unsigned char *all = (unsigned char *)malloc(n);
                double sven = 0.0, sb = 0.0, mx = 0.0, one[3] = {1, 1, 1};
                double tr = Lv->res, vxr[3];
                size_t v;

                vxr[0] = vxr[1] = vxr[2] = tr;
                for (v = 0; v < n; v++)
                {
                    t2[v] = f[n + v] + f[v];
                    all[v] = 1;
                }
                /* cat_vol_localstat(single(f{2}+f{1}), ..., 1, 1) */
                cat12_vol_localstat_mean(t2, all, Lv->rd, 1, Ygwm);
                for (v = 0; v < n; v++) Yb[v] = (Ygwm[v] > 0.5f);
                cat12_vol_morph_labdistclose(Yb, Lv->rd, 12.0 / tr, one);
                for (v = 0; v < n; v++) t1[v] = (float)Yb[v];
                cat12_vol_smooth3X(t1, Lv->rd, 6.0 / tr); /* Ybd */
                for (v = 0; v < n; v++)
                    if (t1[v] > mx) mx = t1[v];
                for (v = 0; v < n; v++)
                    Yv[v] = ((float)Yb[v] - Ygwm[v] > 0.9f) && (t1[v] / mx > 0.9);
                cat12_vol_morph_labdistopen(Yv, Lv->rd, 3.0, vxr);
                cat12_vol_morph_dilate(Yv, Lv->rd, (1.0 / Lv->rred > 1) ? 1.0 / Lv->rred : 1.0, one);
                cat12_vol_morph_lab(Yv, Lv->rd);
                for (v = 0; v < n; v++) Yven[v] = (float)Yv[v];
                cat12_vol_smooth3X(Yven, Lv->rd, 4.0 / tr);
                for (v = 0; v < n; v++)
                {
                    Yven[v] = Yven[v] * 1.5f;
                    if (Yven[v] > 1) Yven[v] = 1;
                    sven += Yven[v];
                    sb += Yb[v];
                }
                if (sb > 0 && sven / sb > 0.2)
                {
                    float *Ygven = t2;
                    hydro = 1;
                    if (info) info->hydro_levels++;
                    if (opts.verbose)
                        printf("  Very large ventricles detected (%0.2f%% of TIV). "
                               "Use separate class and adapt schedule.\n", 100 * sven / sb);
                    if (it > 1)
                    {
                        int cnt = 0, mxl = 0;
                        double fac;
                        for (i = 0; i < nits; i++)
                        {
                            if (tmpl_no[i] == tmpl_no[it - 1]) cnt++;
                            if (tmpl_no[i] > mxl) mxl = tmpl_no[i];
                        }
                        fac = (mxl - tmpl_no[it - 1] + 2) * 2;
                        if (cnt < fac) fac = cnt;
                        fac = log(fac) / log(sched[it - 1]);
                        for (i = 0; i <= nits; i++)
                        {
                            double s2 = pow(sched[i], fac);
                            if (s2 > sched[i]) sched[i] = s2;
                        }
                    }
                    for (v = 0; v < n; v++) bgexp[v] = (float)exp(coef[2][v]);
                    memcpy(Ygven, Yven, sizeof(float) * n);
                    cat12_vol_smooth3X(Ygven, Lv->rd, 2.0 / tr);
                    for (v = 0; v < n; v++) Yv[v] = (bgexp[v] * Ygven[v] > 0.5f);
                    cat12_vol_morph_lab(Yv, Lv->rd);
                    cat12_vol_morph_dilate(Yv, Lv->rd, (1.0 / Lv->rred > 1) ? 1.0 / Lv->rred : 1.0, one);
                    cat12_vol_morph_lab(Yv, Lv->rd);
                    for (v = 0; v < n; v++) Ygven[v] = (float)Yv[v];
                    cat12_vol_smooth3X(Ygven, Lv->rd, 4.0 / tr);
                    for (v = 0; v < n; v++) Ygven[v] = bgexp[v] * Ygven[v];
                    /* f{n1+1}: the ventricle, f{n1+2}: background */
                    for (v = 0; v < n; v++)
                    {
                        float vv = Yven[v] - Ygwm[v];
                        f[(size_t)n1 * n + v] = (vv > (float)MATLAB_EPS) ? vv : (float)MATLAB_EPS;
                        f[(size_t)(n1 + 1) * n + v] = 1.0f;
                    }
                    /* g{n1+2} and g{n1+1} are plain logs here, as in CAT */
                    for (v = 0; v < n; v++)
                    {
                        double bg = exp(coef[2][v]) - Ygven[v];
                        double gv = Ygven[v];
                        coef[n1 + 1][v] = log((bg > MATLAB_EPS) ? bg : MATLAB_EPS);
                        coef[n1][v] = log((gv > MATLAB_EPS) ? gv : MATLAB_EPS);
                    }
                    ng = n1 + 2;
                }
                free(Ygwm);
                free(Yven);
                free(Yb);
                free(Yv);
                free(all);
            }

            /* lesions: replace the tissue by the template there */
            if (lesion)
            {
                size_t v;
                double mn;
                float *ls = t2;
                memcpy(fbuf, lesion, sizeof(float) * n_in);
                if (resfac[ti] > 1)
                {
                    fw[0] = fw[1] = fw[2] = (resfac[ti] - 1) * 2;
                    spm_smooth_vol(NULL, fbuf, dims, fw, nthreads);
                }
                for (v = 0; v < n_in; v++) dbuf[v] = (double)fbuf[v];
                sample_volume(dbuf, dims, Lv->Mrregs, Lv->rd, ls, nthreads);
                mn = vec_min_finite(ls, n);
                for (v = 0; v < n; v++)
                    if (isnan(ls[v])) ls[v] = (float)mn;
                for (k = 0; k < n1; k++)
                {
                    float *tk = t1;
                    for (v = 0; v < ntpl; v++) dbuf[v] = (double)templates[ti][k * ntpl + v];
                    sample_volume(dbuf, tdims, Lv->Mads, Lv->rd, tk, nthreads);
                    mn = vec_min_finite(tk, n);
                    for (v = 0; v < n; v++)
                        if (isnan(tk[v])) tk[v] = (float)mn;
                    for (v = 0; v < n; v++)
                        f[k * n + v] = f[k * n + v] * (1 - ls[v]) + tk[v] * ls[v];
                }
            }

            /* the background */
            for (k = 0; k < ng - 1; k++)
            {
                size_t v;
                for (v = 0; v < n; v++)
                {
                    float *fe = &f[(size_t)(ng - 1) * n + v];
                    *fe = *fe - f[k * n + v];
                    if (!isfinite(f[k * n + v])) *fe = 0.00001f;
                }
            }
            free(dbuf);
            free(fbuf);

            /* template derivatives of this level */
            ret = tplderiv_build(&T, Lv->rd, ng, coef);
            if (ret) goto done;
            for (i = 0; i < (int)n; i++) bgexp[i] = (float)exp(coef[2][i]);

            ll[0] = ll[1] = 0;
            if (it == 1)
            {
                size_t v;
                int x, yy, z;
                for (z = 0; z < Lv->rd[2]; z++)
                    for (yy = 0; yy < Lv->rd[1]; yy++)
                        for (x = 0; x < Lv->rd[0]; x++)
                        {
                            v = x + (size_t)Lv->rd[0] * (yy + (size_t)Lv->rd[1] * z);
                            y[v] = (float)(x + 1);
                            y[n + v] = (float)(yy + 1);
                            y[2 * n + v] = (float)(z + 1);
                        }
                for (v = 0; v < 3 * n; v++) u[v] = 0.0f;
                for (v = 0; v < n; v++) dt[v] = 1.0f;
            }
            else if (L[prev_ti].rd[0] != Lv->rd[0] || L[prev_ti].rd[1] != Lv->rd[1] ||
                     L[prev_ti].rd[2] != Lv->rd[2])
            {
                /* new resolution: interpolate u, re-integrate y */
                size_t no = (size_t)L[prev_ti].rd[0] * L[prev_ti].rd[1] * L[prev_ti].rd[2], v;
                double *uo = (double *)malloc(sizeof(double) * no);
                float sc = (float)Lv->Mys[0];
                float *un = work;
                if (!uo)
                {
                    ret = 2;
                    goto done;
                }
                for (c = 0; c < 3; c++)
                {
                    for (v = 0; v < no; v++) uo[v] = (double)u[c * no + v];
                    sample_volume(uo, L[prev_ti].rd, Lv->Mys, Lv->rd, un + (size_t)c * n, nthreads);
                }
                free(uo);
                for (v = 0; v < 3 * n; v++)
                {
                    float uv = un[v] / sc;
                    u[v] = isfinite(uv) ? uv : (float)MATLAB_EPS;
                }
                ret = shoot3d(&S, Lv->rd, u, prm, int_N, y, J);
                if (ret) goto done;
                jac_det(Lv->rd, J, dt);
            }
        }

        /* the regularisation of this iteration */
        {
            double vx3 = reg_vx * reg_vx * reg_vx;
            prm[0] = prm[1] = prm[2] = reg_vx;
            for (i = 0; i < 5; i++) prm[3 + i] = opts.reg_scale * rparam[i] * sched[it - 1] * vx3;
            int_N = eul_its[it - 1];
        }

        llo[0] = ll[0];
        llo[1] = ll[1];
        ret = shoot_update(&S, &T, f, u, y, dt, prm, scale, &ll[0], &ll[1], work);
        if (ret) goto done;
        ret = shoot3d(&S, L[ti].rd, u, prm, int_N, y, J);
        if (ret) goto done;
        jac_det(L[ti].rd, J, dt);
        extdt = ext_det(dt, n, &nonfinite);

        if (opts.verbose)
            printf("%5d | %0.2f | %7.4f%8.4f%8.4f | %8.4f\n", it, L[ti].res,
                   ll[0] / (3.0 * n), ll[1] / (3.0 * n), (ll[0] + ll[1]) / (3.0 * n),
                   sched[it - 1]);
        if (opts.verbose) fflush(stdout);

        if (opts.regstr != 4 && opts.hotspot && (extdt > 20 || nonfinite))
        {
            /* RD202101: correction of local hot spots by smoothing */
            double vr = 3.0 / L[ti].res, dtl = 60;
            if (info) info->hotspot_corrections++;
            hotspot_first(u, dt, L[ti].rd, vr, sched[it - 1], t1, t2, msk);
            ret = shoot3d(&S, L[ti].rd, u, prm, int_N, y, J);
            if (ret) goto done;
            jac_det(L[ti].rd, J, dt);
            extdt = ext_det(dt, n, &nonfinite);
            while ((nonfinite && nfilt < 20) || (extdt > dtl && nfilt < 4))
            {
                double vx3 = reg_vx * reg_vx * reg_vx;
                hotspot_repeat(u, dt, L[ti].rd, vr, nfilt, dtl, bgexp, t1, t2, t3, msk);
                for (i = 0; i <= nits; i++) sched[i] = pow(sched[i], 1.1);
                for (i = 0; i < 5; i++) prm[3 + i] = opts.reg_scale * rparam[i] * sched[it] * vx3;
                ret = shoot3d(&S, L[ti].rd, u, prm, int_N, y, J);
                if (ret) goto done;
                jac_det(L[ti].rd, J, dt);
                nfilt++;
                extdt = ext_det(dt, n, &nonfinite);
            }
        }

        if (info) info->level_iterations[ti]++;

        /* avoid unnecessary iterations: switch to the next level early */
        {
            double nu = 3.0 * n, r1 = ll[0] / nu;
            double rr = ll[0] / ((llo[0] > MATLAB_EPS) ? llo[0] : MATLAB_EPS);
            if (opts.regstr > 0 && opts.regstr != 4 && nfilt == 0 && !hydro && r1 < 0.075 &&
                (ti > 0 || (ti == 0 && r1 < 1 && rr < 1 && rr > (1 - 0.01))) &&
                (r1 < 1 && rr < 1 && rr > (1 - ll1th)) && tmpl_no[it - 1] > 1)
            {
                int nx = nits + 1;
                for (i = 0; i < nits; i++)
                    if (tmpl_no[i] > tmpl_no[it - 1])
                    {
                        nx = i + 1;
                        break;
                    }
                it = (it + 1 > nx) ? it + 1 : nx;
            }
            else
                it++;
        }
        if (info) info->iterations++;
    }

    /* outputs: y in input voxels (0-based), Jacobians, inverse.  The last
     * registration grid has the template's size, but CAT's corner logic
     * (if tmpbb(1)<tmpbb(2), swap) orders it like a LAS image: for any other
     * orientation it is mirrored against the template grid, so the results
     * are reordered through Mads. */
    {
        const double *Mr = L[NLEV - 1].Mrregs, *Ma = L[NLEV - 1].Mads;
        double iMa[16], dlin, dads, dmin = HUGE_VAL, dmax = -HUGE_VAL;
        const int *rd = L[NLEV - 1].rd;
        size_t v;
        int t0, t1, t2;

        if (m4_inv(Ma, iMa))
        {
            ret = 1;
            goto done;
        }
        dlin = Mr[0] * (Mr[5] * Mr[10] - Mr[6] * Mr[9]) - Mr[1] * (Mr[4] * Mr[10] - Mr[6] * Mr[8]) +
               Mr[2] * (Mr[4] * Mr[9] - Mr[5] * Mr[8]);
        dads = Ma[0] * (Ma[5] * Ma[10] - Ma[6] * Ma[9]) - Ma[1] * (Ma[4] * Ma[10] - Ma[6] * Ma[8]) +
               Ma[2] * (Ma[4] * Ma[9] - Ma[5] * Ma[8]);
        for (t2 = 0; t2 < tdims[2]; t2++)
            for (t1 = 0; t1 < tdims[1]; t1++)
                for (t0 = 0; t0 < tdims[0]; t0++)
                {
                    size_t tv = t0 + (size_t)tdims[0] * (t1 + (size_t)tdims[1] * t2), rv;
                    double r[3], a, b, cc;
                    int ri[3];
                    for (j = 0; j < 3; j++)
                    {
                        r[j] = iMa[j * 4 + 0] * (t0 + 1) + iMa[j * 4 + 1] * (t1 + 1) +
                               iMa[j * 4 + 2] * (t2 + 1) + iMa[j * 4 + 3];
                        ri[j] = (int)floor(r[j] + 0.5) - 1;
                        if (ri[j] < 0) ri[j] = 0;
                        if (ri[j] >= rd[j]) ri[j] = rd[j] - 1;
                    }
                    rv = ri[0] + (size_t)rd[0] * (ri[1] + (size_t)rd[1] * ri[2]);
                    a = y[rv];
                    b = y[ntpl + rv];
                    cc = y[2 * ntpl + rv];
                    for (j = 0; j < 3; j++)
                        y_out[j * ntpl + tv] = (float)(Mr[j * 4 + 0] * a + Mr[j * 4 + 1] * b +
                                                       Mr[j * 4 + 2] * cc + Mr[j * 4 + 3] - 1.0);
                    if (jac_out) jac_out[tv] = (float)(dt[rv] * dlin / dads);
                }
        for (v = 0; v < ntpl; v++)
        {
            if (dt[v] < dmin) dmin = dt[v];
            if (dt[v] > dmax) dmax = dt[v];
        }
        if (info)
        {
            info->ll_match = ll[0] / (3.0 * ntpl);
            info->ll_reg = ll[1] / (3.0 * ntpl);
            info->det_min = dmin;
            info->det_max = dmax;
        }
        if (iy_out)
        {
            mwSize dy[3], diy[3];
            float M1[4][3], M2[4][3];
            for (i = 0; i < 3; i++)
            {
                dy[i] = rd[i];
                diy[i] = dims[i];
            }
            /* y's values are registration voxels, which Mrregs maps to input
             * voxels (M1: "mm" -> voxels of the inverse field); the inverse
             * holds registration voxels, which Mads maps to template voxels
             * (M2: voxels -> "mm" of the forward field) */
            for (i = 0; i < 3; i++)
                for (j = 0; j < 4; j++)
                {
                    M1[j][i] = (float)Mr[i * 4 + j];
                    M2[j][i] = (float)Ma[i * 4 + j];
                }
            memset(iy_out, 0, sizeof(float) * 3 * n_in);
            spm_invdef(dy, y, diy, iy_out, M1, M2);
            for (v = 0; v < 3 * n_in; v++) iy_out[v] -= 1.0f;
        }
    }

done:
    free(sched);
    free(eul_its);
    free(tmpl_no);
    free(f);
    free(u);
    free(y);
    free(J);
    free(dt);
    free(work);
    free(t1);
    free(t2);
    free(t3);
    free(bgexp);
    free(msk);
    for (k = 0; k < MAXCLS; k++) free(coef[k]);
    tplderiv_free(&T);
    shoot_free(&S);
    spm_shoot_set_threads(old_threads);
    return ret;
}
