/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Unified segmentation with geodesic shooting.  The generative model and its
 * estimation follow SPM's unified segmentation as CAT runs it
 * (cat_spm_preproc_run.m, cat_spm_preproc8.m, cat_spm_preproc_write8.m and
 * spm_maff8.m; John Ashburner, with the changes of Christian Gaser and Robert
 * Dahnke): a mixture of Gaussians per class, a bias field of DCT basis
 * functions, class weights on the priors.  The priors are the templates of a
 * shooting series, and their deformation is estimated by geodesic shooting
 * on a grid of increasing resolution (CAT_VolShootFit()) instead of SPM's
 * small-deformation model.
 *
 * Copyright Christian Gaser, University of Jena.
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <stdint.h>
#include <stddef.h>
#include <time.h>

#include "CAT_UnifiedSeg.h"
#include "CAT_Nlm.h"
#include "CAT_Vol.h"
#include "spm_shoot.h"

#define KB CAT_USEG_NCLASSES
#define MAXK CAT_USEG_MAXK
#define NLEV CAT_SHOOT_NTEMPLATES
#define USEG_MAX_THREADS 16
#define TPL_TINY 1e-9    /* keeps the log of a prior finite; SPM's tiny is 1e-4 */
#define MOM_TINY 4.9e-32 /* eps*eps of spm_preproc8.m */
#define DBL_EPS 2.220446049250313e-16
#define AFF_SAMP_COARSE 8.0 /* mm between the voxels of the coarse affine search */
#define AFF_SAMP_FINE 5.0   /* and of its refinement */
#define PVE_NTAB 1024       /* intensity bins of the partial volume tables */
#define PVE_NINT 32         /* nodes of the integral over the mixing fraction */
#define HEAD_NTAB_MAX 16384 /* intensity bins of the tables of the head */
#define HEAD_NPAIR 2        /* mixtures of the head: bone-soft, soft-air */
#define CLEAN_TH_FINAL 0.05 /* threshold of the masks of clean_gwc */

/* template classes */
#define C_GM 0
#define C_WM 1
#define C_CSF 2
#define C_BONE 3
#define C_SOFT 4
#define C_AIR 5

/* ======================================================================== */
/* Small helpers                                                            */
/* ======================================================================== */

static double
wall_time(void)
{
#if defined(_WIN32) || defined(_WIN64)
    return (double)clock() / CLOCKS_PER_SEC;
#else
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + 1e-9 * (double)ts.tv_nsec;
#endif
}

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

/* Run fn over [0, n) in contiguous chunks.  Every caller writes only to the
 * elements of its own units and keeps sums per unit, which are added in the
 * order of the units afterwards: results do not depend on the thread count. */
static void
par_for(int n, int nthreads, range_fn fn, void *ctx)
{
    range_args args[USEG_MAX_THREADS];
    int t;

    if (nthreads > n) nthreads = n;
    if (nthreads > USEG_MAX_THREADS) nthreads = USEG_MAX_THREADS;
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

/* Solve the n x n system A x = b in place (x in b) by Gaussian elimination
 * with partial pivoting; A is destroyed. */
static int
solve_dense(double *A, double *b, int n)
{
    int i, j, k;
    for (k = 0; k < n; k++)
    {
        int p = k;
        for (i = k + 1; i < n; i++)
            if (fabs(A[i * n + k]) > fabs(A[p * n + k])) p = i;
        if (fabs(A[p * n + k]) < 1e-300) return 1;
        if (p != k)
        {
            double t;
            for (j = 0; j < n; j++)
            {
                t = A[k * n + j];
                A[k * n + j] = A[p * n + j];
                A[p * n + j] = t;
            }
            t = b[k];
            b[k] = b[p];
            b[p] = t;
        }
        for (i = k + 1; i < n; i++)
        {
            double f = A[i * n + k] / A[k * n + k];
            for (j = k; j < n; j++) A[i * n + j] -= f * A[k * n + j];
            b[i] -= f * b[k];
        }
    }
    for (k = n - 1; k >= 0; k--)
    {
        double s = b[k];
        for (j = k + 1; j < n; j++) s -= A[k * n + j] * b[j];
        b[k] = s / A[k * n + k];
    }
    return 0;
}

static int
m4_inv(const double *A, double *Ai)
{
    double M[16], col[4];
    int i, j;
    for (j = 0; j < 4; j++)
    {
        memcpy(M, A, sizeof(M));
        for (i = 0; i < 4; i++) col[i] = (i == j) ? 1.0 : 0.0;
        if (solve_dense(M, col, 4)) return 1;
        for (i = 0; i < 4; i++) Ai[i * 4 + j] = col[i];
    }
    return 0;
}

static int
is_valid(float v)
{
    return isfinite(v) && v != 0.0f;
}

/* The k-th smallest of a (0-based) by Hoare selection; a is reordered. */
static float
select_float(float *a, ptrdiff_t n, ptrdiff_t k)
{
    ptrdiff_t lo = 0, hi = n - 1;
    while (lo < hi)
    {
        float p = a[lo + (hi - lo) / 2];
        ptrdiff_t i = lo, j = hi;
        while (i <= j)
        {
            while (a[i] < p) i++;
            while (a[j] > p) j--;
            if (i <= j)
            {
                float t = a[i];
                a[i] = a[j];
                a[j] = t;
                i++;
                j--;
            }
        }
        if (k <= j) hi = j;
        else if (k >= i) lo = i;
        else break;
    }
    return a[k];
}

/* splitmix64 and a normal deviate from it (Box-Muller) */
static uint64_t
rng_next(uint64_t *s)
{
    uint64_t z = (*s += 0x9E3779B97F4A7C15ULL);
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    return z ^ (z >> 31);
}

static double
rng_normal(uint64_t *s)
{
    double u1 = ((double)(rng_next(s) >> 11) + 0.5) * (1.0 / 9007199254740992.0);
    double u2 = ((double)(rng_next(s) >> 11) + 0.5) * (1.0 / 9007199254740992.0);
    return sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2);
}

/* Inverse of the standard normal distribution function (Acklam). */
static double
norm_inv(double p)
{
    static const double a[6] = {-3.969683028665376e+01, 2.209460984245205e+02,
                                -2.759285104469687e+02, 1.383577518672690e+02,
                                -3.066479806614716e+01, 2.506628277459239e+00};
    static const double b[5] = {-5.447609879822406e+01, 1.615858368580409e+02,
                                -1.556989798598866e+02, 6.680131188771972e+01,
                                -1.328068155288572e+01};
    static const double c[6] = {-7.784894002430293e-03, -3.223964580411365e-01,
                                -2.400758277161838e+00, -2.549732539343734e+00,
                                4.374664141464968e+00, 2.938163982698783e+00};
    static const double d[4] = {7.784695709041462e-03, 3.224671290700398e-01,
                                2.445134137142996e+00, 3.754408661907416e+00};
    double q, r;

    if (p < 0.02425)
    {
        q = sqrt(-2.0 * log(p));
        return (((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) /
               ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1.0);
    }
    if (p > 1.0 - 0.02425)
    {
        q = sqrt(-2.0 * log(1.0 - p));
        return -(((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) /
               ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1.0);
    }
    q = p - 0.5;
    r = q * q;
    return (((((a[0] * r + a[1]) * r + a[2]) * r + a[3]) * r + a[4]) * r + a[5]) * q /
           (((((b[0] * r + b[1]) * r + b[2]) * r + b[3]) * r + b[4]) * r + 1.0);
}

/* ======================================================================== */
/* Intensity range, noise                                                   */
/* ======================================================================== */

/* 1st and 99th percentile of the voxels that carry data. */
static int
robust_range(const float *v, size_t n, const unsigned char *mask, double range[2])
{
    size_t step = n / 2000000 + 1, i;
    ptrdiff_t m = 0;
    float *buf = (float *)malloc(sizeof(float) * (n / step + 1));

    if (!buf) return 2;
    for (i = 0; i < n; i += step)
        if (is_valid(v[i]) && (!mask || mask[i])) buf[m++] = v[i];
    if (m < 100)
    {
        free(buf);
        return 3;
    }
    range[0] = select_float(buf, m, (ptrdiff_t)(0.01 * (double)(m - 1)));
    range[1] = select_float(buf, m, (ptrdiff_t)(0.99 * (double)(m - 1)));
    free(buf);
    return (range[1] > range[0]) ? 0 : 3;
}

/* Noise SD from the finest diagonal detail of a Haar transform: over a
 * 2x2x2 block the detail (+ - - + - + + -)/sqrt(8) has the variance of the
 * noise and cancels every signal up to its mixed third derivative.  The SD
 * is the median absolute detail / 0.6745 over the blocks brighter than
 * thresh (tissue, not air, whose noise is Rayleigh distributed). */
static double
estimate_noise(const float *v, const int d[3], const unsigned char *mask, double thresh)
{
    size_t nxy = (size_t)d[0] * d[1];
    size_t nmax = (size_t)(d[0] / 2) * (size_t)(d[1] / 2) * (size_t)(d[2] / 2);
    ptrdiff_t m = 0;
    float *buf, med;
    int x, y, z, c;

    if (nmax < 100) return 0.0;
    buf = (float *)malloc(sizeof(float) * nmax);
    if (!buf) return 0.0;
    for (z = 0; z + 1 < d[2]; z += 2)
        for (y = 0; y + 1 < d[1]; y += 2)
            for (x = 0; x + 1 < d[0]; x += 2)
            {
                size_t i = x + (size_t)d[0] * y + nxy * z;
                double sum = 0.0, det = 0.0;
                int ok = 1;
                for (c = 0; c < 8; c++)
                {
                    size_t j = i + (c & 1) + (size_t)d[0] * ((c >> 1) & 1) + nxy * ((c >> 2) & 1);
                    int par = (c & 1) + ((c >> 1) & 1) + ((c >> 2) & 1);
                    if (!is_valid(v[j]) || (mask && !mask[j]))
                    {
                        ok = 0;
                        break;
                    }
                    sum += v[j];
                    det += (par & 1) ? -(double)v[j] : (double)v[j];
                }
                if (ok && sum / 8.0 > thresh) buf[m++] = (float)fabs(det / sqrt(8.0));
            }
    if (m < 100)
    {
        free(buf);
        return 0.0;
    }
    med = select_float(buf, m, m / 2);
    free(buf);
    return (double)med / 0.6745;
}

typedef struct
{
    const float *in;
    float *out;
    const int *d;
    double sd;
    uint64_t seed;
} noise_ctx;

/* Every slice has its own stream, so the noise does not depend on threads. */
static void
noise_range(void *p, int z0, int z1)
{
    noise_ctx *c = (noise_ctx *)p;
    size_t nxy = (size_t)c->d[0] * c->d[1], i;
    int z;

    for (z = z0; z < z1; z++)
    {
        uint64_t s = c->seed * 0x2545F4914F6CDD1DULL + (uint64_t)(z + 1) * 0xD6E8FEB86659FD93ULL;
        const float *in = c->in + nxy * z;
        float *out = c->out + nxy * z;
        (void)rng_next(&s);
        for (i = 0; i < nxy; i++)
            out[i] = is_valid(in[i]) ? (float)(in[i] + c->sd * rng_normal(&s)) : in[i];
    }
}

/* ======================================================================== */
/* Templates                                                                */
/* ======================================================================== */

typedef struct
{
    int d[3];
    size_t n;
    const float *const *tpl; /* the series */
    float *lt;               /* KB*n: log(template + tiny), class fastest */
    int level;               /* level held by lt, -1: none */
    double bg1[KB], bg2[KB]; /* log of the mean of its bottom and top slice */
} Tpl;

typedef struct
{
    Tpl *T;
    const float *src;
    double pmin[KB]; /* the floor of every class */
} tplload_ctx;

/* The priors of template voxel i (n voxels per class) without the floor
 * pmin of each class, summing to one. */
static void
tpl_prior(const float *src, size_t n, size_t i, const double pmin[KB], double v[KB])
{
    double sum = 0.0;
    int k;

    for (k = 0; k < KB; k++)
    {
        v[k] = (double)src[k * n + i] - pmin[k];
        if (!(v[k] > 0)) v[k] = 0.0;
        sum += v[k];
    }
    if (!(sum > 0)) sum = 1.0;
    for (k = 0; k < KB; k++) v[k] /= sum;
}

static void
tplload_range(void *p, int z0, int z1)
{
    tplload_ctx *c = (tplload_ctx *)p;
    size_t nxy = (size_t)c->T->d[0] * c->T->d[1], i;
    double v[KB];
    int k;

    for (i = nxy * z0; i < nxy * z1; i++)
    {
        tpl_prior(c->src, c->T->n, i, c->pmin, v);
        for (k = 0; k < KB; k++) c->T->lt[KB * i + k] = (float)log(v[k] + TPL_TINY);
    }
}

/* Load the template of a level.  Beyond its volume the prior is the mean of
 * its bottom slice (below) or of its top slice (anywhere else), as bg1 and
 * bg2 of spm_load_priors8.m -- but without the brain classes, which SPM
 * keeps there and which then label neck and shoulders as CSF.
 *
 * A class that the template rules out gets no floor.  SPM adds 1e-4 to every
 * prior, and shooting templates come with a floor of their own (the smallest
 * value of a class, 1e-4 in CAT's series).  With such a floor the classes
 * outside the brain act as outlier classes inside it: white matter that is
 * four to five SD brighter than its Gaussian -- a gyrus next to a coil
 * element, where the bias field is steeper than its basis allows -- has a
 * likelihood ratio of more than 1e4 in favour of the broad Gaussians of soft
 * tissue, becomes soft tissue, and leaves a hole in the brain.  So the floor
 * of every class is taken off, in the volume and in the slices that stand for
 * what lies beyond it, and TPL_TINY only keeps the logarithm finite. */
static void
tpl_load(Tpl *T, int level, int nthreads)
{
    tplload_ctx c;
    size_t nxy = (size_t)T->d[0] * T->d[1], i;
    double v[KB], s1[KB], s2[KB];
    int k;

    if (T->level == level) return;
    c.T = T;
    c.src = T->tpl[level];
    for (k = 0; k < KB; k++)
    {
        const float *t = c.src + (size_t)k * T->n;
        float mn = t[0];
        for (i = 1; i < T->n; i++)
            if (t[i] < mn) mn = t[i];
        c.pmin[k] = (mn > 0) ? (double)mn : 0.0;
    }
    par_for(T->d[2], nthreads, tplload_range, &c);
    for (k = 0; k < KB; k++) s1[k] = s2[k] = 0.0;
    for (i = 0; i < nxy; i++)
    {
        tpl_prior(c.src, T->n, i, c.pmin, v);
        for (k = 0; k < KB; k++) s1[k] += v[k];
        tpl_prior(c.src, T->n, T->n - nxy + i, c.pmin, v);
        for (k = 0; k < KB; k++) s2[k] += v[k];
    }
    for (k = 0; k < KB; k++)
    {
        T->bg1[k] = log(s1[k] / (double)nxy + TPL_TINY);
        T->bg2[k] = log(s2[k] / (double)nxy + TPL_TINY);
        /* the templates hold the whole brain: there is none outside them */
        if (k == C_GM || k == C_WM || k == C_CSF) T->bg1[k] = T->bg2[k] = log(TPL_TINY);
    }
    T->level = level;
}

/* The log template at template voxel (x, y, z), 0-based, with the quadratic
 * B-spline kernel on its samples, as spm_sample_priors8.m does (degree 2 on
 * coefficients of degree 1): smooth, with continuous derivatives g (may be
 * NULL).  Outside the volume it is the constant bg1 or bg2. */
static void
tpl_sample(const Tpl *T, double x, double y, double z, double a[KB], double g[KB][3])
{
    const double c[3] = {x, y, z};
    double w[3][3], dw[3][3];
    size_t ix[3][3];
    int i, j, k, kx, ky, kz;

    if (!(x >= 0.0 && x <= T->d[0] - 1 && y >= 0.0 && y <= T->d[1] - 1 && z >= 0.0 &&
          z <= T->d[2] - 1))
    {
        const double *bg = (z < 0.0) ? T->bg1 : T->bg2;
        for (k = 0; k < KB; k++)
        {
            a[k] = bg[k];
            if (g) g[k][0] = g[k][1] = g[k][2] = 0.0;
        }
        return;
    }
    for (i = 0; i < 3; i++)
    {
        double ci = c[i], t;
        int ic;
        ic = (int)floor(ci + 0.5);
        t = ci - ic;
        w[i][0] = 0.5 * (0.5 - t) * (0.5 - t);
        w[i][1] = 0.75 - t * t;
        w[i][2] = 0.5 * (0.5 + t) * (0.5 + t);
        dw[i][0] = -(0.5 - t);
        dw[i][1] = -2.0 * t;
        dw[i][2] = (0.5 + t);
        for (j = 0; j < 3; j++)
        {
            int id = ic - 1 + j;
            if (id < 0) id = 0;
            if (id > T->d[i] - 1) id = T->d[i] - 1;
            ix[i][j] = (size_t)id;
        }
    }
    for (k = 0; k < KB; k++)
    {
        a[k] = 0.0;
        if (g) g[k][0] = g[k][1] = g[k][2] = 0.0;
    }
    for (kz = 0; kz < 3; kz++)
        for (ky = 0; ky < 3; ky++)
        {
            size_t base = (size_t)T->d[0] * (ix[1][ky] + (size_t)T->d[1] * ix[2][kz]);
            double wyz = w[1][ky] * w[2][kz];
            for (kx = 0; kx < 3; kx++)
            {
                const float *p = T->lt + KB * (base + ix[0][kx]);
                double ww = w[0][kx] * wyz;
                if (g)
                {
                    double gx = dw[0][kx] * wyz, gy = w[0][kx] * dw[1][ky] * w[2][kz];
                    double gz = w[0][kx] * w[1][ky] * dw[2][kz];
                    for (k = 0; k < KB; k++)
                    {
                        a[k] += ww * p[k];
                        g[k][0] += gx * p[k];
                        g[k][1] += gy * p[k];
                        g[k][2] += gz * p[k];
                    }
                }
                else
                    for (k = 0; k < KB; k++) a[k] += ww * p[k];
            }
        }
}

/* b = softmax(a + lw), lb = log(b); lw may be NULL */
static void
softmax_k(const double a[KB], const double *lw, double b[KB], double lb[KB])
{
    double mx = -HUGE_VAL, s = 0.0, t[KB];
    int k;
    for (k = 0; k < KB; k++)
    {
        t[k] = a[k] + (lw ? lw[k] : 0.0);
        if (t[k] > mx) mx = t[k];
    }
    for (k = 0; k < KB; k++)
    {
        b[k] = exp(t[k] - mx);
        s += b[k];
    }
    for (k = 0; k < KB; k++)
    {
        b[k] /= s;
        if (lb) lb[k] = t[k] - mx - log(s);
    }
}

/* ======================================================================== */
/* Affine registration (spm_maff8.m)                                        */
/* ======================================================================== */

/* mean and inverse covariance of the zooms and shears: spm_affine_priors.m */
static void
affine_prior(int typ, double mu[6], double isig[36])
{
    static const double mu_mni[6] = {0.0667, 0.0039, 0.0008, 0.0333, 0.0071, 0.1071};
    static const double is_mni[36] = {
        0.0902, -0.0345, -0.0106, -0.0025, -0.0005, -0.0163,
        -0.0345, 0.7901, 0.3883, 0.0041, -0.0103, -0.0116,
        -0.0106, 0.3883, 2.2599, 0.0113, 0.0396, -0.0060,
        -0.0025, 0.0041, 0.0113, 0.0925, 0.0471, -0.0440,
        -0.0005, -0.0103, 0.0396, 0.0471, 0.2964, -0.0062,
        -0.0163, -0.0116, -0.0060, -0.0440, -0.0062, 0.1144};
    static const double mu_east[6] = {0.0719, -0.0040, -0.0032, 0.1416, 0.0601, 0.2578};
    static const double is_east[36] = {
        0.0757, 0.0220, -0.0224, -0.0049, 0.0304, -0.0327,
        0.0220, 0.3125, -0.1555, 0.0280, -0.0012, -0.0284,
        -0.0224, -0.1555, 1.9727, 0.0196, -0.0019, 0.0122,
        -0.0049, 0.0280, 0.0196, 0.0576, -0.0282, -0.0200,
        0.0304, -0.0012, -0.0019, -0.0282, 0.2128, -0.0275,
        -0.0327, -0.0284, 0.0122, -0.0200, -0.0275, 0.0511};
    static const double is_subj[36] = {
        0.8876, 0.0784, 0.0784, -0.1749, 0.0784, -0.1749,
        0.0784, 5.3894, 0.2655, 0.0784, 0.2655, 0.0784,
        0.0784, 0.2655, 5.3894, 0.0784, 0.2655, 0.0784,
        -0.1749, 0.0784, 0.0784, 0.8876, 0.0784, -0.1749,
        0.0784, 0.2655, 0.2655, 0.0784, 5.3894, 0.0784,
        -0.1749, 0.0784, 0.0784, -0.1749, 0.0784, 0.8876};
    int i;

    for (i = 0; i < 6; i++) mu[i] = 0.0;
    for (i = 0; i < 36; i++) isig[i] = 0.0;
    switch (typ)
    {
    case CAT_USEG_AFFREG_MNI:
        for (i = 0; i < 6; i++) mu[i] = mu_mni[i];
        for (i = 0; i < 36; i++) isig[i] = 1e4 * is_mni[i];
        break;
    case CAT_USEG_AFFREG_EASTERN:
        for (i = 0; i < 6; i++) mu[i] = mu_east[i];
        for (i = 0; i < 36; i++) isig[i] = 1e4 * is_east[i];
        break;
    case CAT_USEG_AFFREG_SUBJ:
        for (i = 0; i < 36; i++) isig[i] = 1e3 * is_subj[i];
        break;
    case CAT_USEG_AFFREG_RIGID:
        for (i = 0; i < 6; i++) isig[i * 7] = 1e8;
        break;
    default:
        break;
    }
}

/* Eigenvalues w and eigenvectors (columns of V) of a symmetric 3x3 matrix
 * by Jacobi rotations. */
static void
sym3_eig(const double A[9], double w[3], double V[9])
{
    double a[3][3], v[3][3];
    int sweep, p, q, k;

    for (p = 0; p < 3; p++)
        for (q = 0; q < 3; q++)
        {
            a[p][q] = A[p * 3 + q];
            v[p][q] = (p == q) ? 1.0 : 0.0;
        }
    for (sweep = 0; sweep < 60; sweep++)
    {
        double off = a[0][1] * a[0][1] + a[0][2] * a[0][2] + a[1][2] * a[1][2];
        if (off < 1e-300) break;
        for (p = 0; p < 2; p++)
            for (q = p + 1; q < 3; q++)
            {
                double theta, t, c, s;
                if (fabs(a[p][q]) < 1e-300) continue;
                theta = (a[q][q] - a[p][p]) / (2.0 * a[p][q]);
                t = ((theta >= 0) ? 1.0 : -1.0) / (fabs(theta) + sqrt(theta * theta + 1.0));
                c = 1.0 / sqrt(t * t + 1.0);
                s = t * c;
                for (k = 0; k < 3; k++)
                {
                    double akp = a[k][p], akq = a[k][q];
                    a[k][p] = c * akp - s * akq;
                    a[k][q] = s * akp + c * akq;
                }
                for (k = 0; k < 3; k++)
                {
                    double apk = a[p][k], aqk = a[q][k];
                    a[p][k] = c * apk - s * aqk;
                    a[q][k] = s * apk + c * aqk;
                }
                for (k = 0; k < 3; k++)
                {
                    double vkp = v[k][p], vkq = v[k][q];
                    v[k][p] = c * vkp - s * vkq;
                    v[k][q] = s * vkp + c * vkq;
                }
            }
    }
    for (p = 0; p < 3; p++)
    {
        w[p] = a[p][p];
        for (q = 0; q < 3; q++) V[p * 3 + q] = v[p][q];
    }
}

/* B = V diag(f(w)) V' with fun 0: exp, 1: sqrt, 2: 1/sqrt, 3: log/2 */
static void
sym3_fun(const double A[9], int fun, double B[9])
{
    double w[3], V[9];
    int i, j, k;

    sym3_eig(A, w, V);
    for (k = 0; k < 3; k++)
    {
        double x = w[k];
        if (fun == 0) w[k] = exp(x);
        else
        {
            if (x < 1e-300) x = 1e-300;
            w[k] = (fun == 1) ? sqrt(x) : ((fun == 2) ? 1.0 / sqrt(x) : 0.5 * log(x));
        }
    }
    for (i = 0; i < 3; i++)
        for (j = 0; j < 3; j++)
        {
            double s = 0.0;
            for (k = 0; k < 3; k++) s += V[i * 3 + k] * w[k] * V[j * 3 + k];
            B[i * 3 + j] = s;
        }
}

/* P2M of spm_maff8.m: translations, the logarithm of the rotation and of
 * the symmetric part (zooms and shears) to a 4x4 matrix. */
static void
aff_P2M(const double P[12], double M[16])
{
    double S[9] = {0, P[3], P[4], -P[3], 0, P[5], -P[4], -P[5], 0};
    double Vs[9] = {P[6], P[7], P[8], P[7], P[9], P[10], P[8], P[10], P[11]};
    double R[9], V[9], S2[9], th2 = P[3] * P[3] + P[4] * P[4] + P[5] * P[5];
    double th = sqrt(th2), f1, f2;
    int i, j, k;

    /* Rodrigues: expm(S) = I + sin(th)/th S + (1 - cos(th))/th^2 S^2 */
    if (th < 1e-8)
    {
        f1 = 1.0 - th2 / 6.0;
        f2 = 0.5 - th2 / 24.0;
    }
    else
    {
        f1 = sin(th) / th;
        f2 = (1.0 - cos(th)) / th2;
    }
    for (i = 0; i < 3; i++)
        for (j = 0; j < 3; j++)
        {
            double s = 0.0;
            for (k = 0; k < 3; k++) s += S[i * 3 + k] * S[k * 3 + j];
            S2[i * 3 + j] = s;
        }
    for (i = 0; i < 9; i++) R[i] = ((i % 4 == 0) ? 1.0 : 0.0) + f1 * S[i] + f2 * S2[i];
    sym3_fun(Vs, 0, V);
    m4_eye(M);
    for (i = 0; i < 3; i++)
    {
        for (j = 0; j < 3; j++)
        {
            double s = 0.0;
            for (k = 0; k < 3; k++) s += V[i * 3 + k] * R[k * 3 + j];
            M[i * 4 + j] = s;
        }
        M[i * 4 + 3] = P[i];
    }
}

/* M2P of spm_maff8.m, by the polar decomposition J = V R. */
static int
aff_M2P(const double M[16], double P[12])
{
    double J[9], JJt[9], Vi[9], R[9], lV[9], tr, c, th, f;
    int i, j, k;

    for (i = 0; i < 3; i++)
        for (j = 0; j < 3; j++) J[i * 3 + j] = M[i * 4 + j];
    if (J[0] * (J[4] * J[8] - J[5] * J[7]) - J[1] * (J[3] * J[8] - J[5] * J[6]) +
        J[2] * (J[3] * J[7] - J[4] * J[6]) <= 0)
        return 1; /* a reflection has no logarithm */
    for (i = 0; i < 3; i++)
        for (j = 0; j < 3; j++)
        {
            double s = 0.0;
            for (k = 0; k < 3; k++) s += J[i * 3 + k] * J[j * 3 + k];
            JJt[i * 3 + j] = s;
        }
    sym3_fun(JJt, 2, Vi);
    sym3_fun(JJt, 3, lV);
    for (i = 0; i < 3; i++)
        for (j = 0; j < 3; j++)
        {
            double s = 0.0;
            for (k = 0; k < 3; k++) s += Vi[i * 3 + k] * J[k * 3 + j];
            R[i * 3 + j] = s;
        }
    tr = R[0] + R[4] + R[8];
    c = 0.5 * (tr - 1.0);
    if (c > 1.0) c = 1.0;
    if (c < -1.0) c = -1.0;
    th = acos(c);
    f = (th < 1e-8) ? 0.5 : th / (2.0 * sin(th));
    P[0] = M[3];
    P[1] = M[7];
    P[2] = M[11];
    /* lR = -logm(R); P(4:6) = lR(2,1), lR(3,1), lR(3,2) */
    P[3] = -f * (R[3] - R[1]);
    P[4] = -f * (R[6] - R[2]);
    P[5] = -f * (R[7] - R[5]);
    P[6] = lV[0];
    P[7] = lV[3];
    P[8] = lV[6];
    P[9] = lV[4];
    P[10] = lV[7];
    P[11] = lV[8];
    return 0;
}

/* Voxels of the image, a few mm apart, with intensities in 256 bins. */
typedef struct
{
    int d[3], sk[3];
    size_t n;
    unsigned char *g, *msk;
    double ff; /* fudge factor for the dependence of neighbouring voxels */
} AffLat;

static void
afflat_free(AffLat *A)
{
    free(A->g);
    free(A->msk);
    memset(A, 0, sizeof(AffLat));
}

static double
fudge_factor(double fwhm, const double vx[3], const int sk[3])
{
    double s = (fwhm + (vx[0] + vx[1] + vx[2]) / 3.0) / sqrt(8.0 * log(2.0)), p = 1.0;
    int i;
    for (i = 0; i < 3; i++)
    {
        double r = s / vx[i] / sk[i];
        p *= 4.0 * M_PI * r * r + 1.0;
    }
    return sqrt(p);
}

/* loadbuf of spm_maff8.m */
static int
afflat_build(AffLat *A, const float *src, const unsigned char *mask, const int dims[3],
             const double vx[3], double samp)
{
    size_t nxy = (size_t)dims[0] * dims[1], i, nv = 0;
    double mn = HUGE_VAL, mx = -HUGE_VAL, h[4000], cum = 0.0, lo, hi;
    int x, y, z, k, ilo = -1, ihi = -1;

    memset(A, 0, sizeof(AffLat));
    for (k = 0; k < 3; k++)
    {
        A->sk[k] = (int)floor(samp / vx[k] + 0.5);
        if (A->sk[k] < 1) A->sk[k] = 1;
        A->d[k] = (dims[k] - 1) / A->sk[k] + 1;
    }
    A->n = (size_t)A->d[0] * A->d[1] * A->d[2];
    A->g = (unsigned char *)calloc(A->n, 1);
    A->msk = (unsigned char *)calloc(A->n, 1);
    if (!A->g || !A->msk) return 2;

    for (z = 0; z < A->d[2]; z++)
        for (y = 0; y < A->d[1]; y++)
            for (x = 0; x < A->d[0]; x++)
            {
                size_t v = (size_t)x * A->sk[0] + (size_t)dims[0] * ((size_t)y * A->sk[1]) +
                           nxy * ((size_t)z * A->sk[2]);
                i = x + (size_t)A->d[0] * (y + (size_t)A->d[1] * z);
                if (!is_valid(src[v]) || (mask && !mask[v])) continue;
                A->msk[i] = 1;
                nv++;
                if (src[v] < mn) mn = src[v];
                if (src[v] > mx) mx = src[v];
            }
    if (nv < 1000 || !(mx > mn)) return 3;

    /* the range between the 0.05% and 99.95% quantiles becomes 0..255 */
    memset(h, 0, sizeof(h));
    for (z = 0; z < A->d[2]; z++)
        for (y = 0; y < A->d[1]; y++)
            for (x = 0; x < A->d[0]; x++)
            {
                size_t v = (size_t)x * A->sk[0] + (size_t)dims[0] * ((size_t)y * A->sk[1]) +
                           nxy * ((size_t)z * A->sk[2]);
                i = x + (size_t)A->d[0] * (y + (size_t)A->d[1] * z);
                if (!A->msk[i]) continue;
                k = (int)floor((src[v] - mn) / (mx - mn) * 3999.0 + 0.5);
                h[(k < 0) ? 0 : ((k > 3999) ? 3999 : k)] += 1.0;
            }
    for (k = 0; k < 4000; k++)
    {
        cum += h[k] / (double)nv;
        if (ilo < 0 && cum > 0.0005) ilo = k;
        if (ihi < 0 && cum > 0.9995) ihi = k;
    }
    if (ihi < 0) ihi = 3999;
    if (ihi <= ilo) ihi = ilo + 1;
    lo = mn + (mx - mn) * ilo / 3999.0;
    hi = mn + (mx - mn) * ihi / 3999.0;
    for (z = 0; z < A->d[2]; z++)
        for (y = 0; y < A->d[1]; y++)
            for (x = 0; x < A->d[0]; x++)
            {
                size_t v = (size_t)x * A->sk[0] + (size_t)dims[0] * ((size_t)y * A->sk[1]) +
                           nxy * ((size_t)z * A->sk[2]);
                double q;
                i = x + (size_t)A->d[0] * (y + (size_t)A->d[1] * z);
                if (!A->msk[i]) continue;
                q = floor((src[v] - lo) / (hi - lo) * 255.0 + 0.5);
                A->g[i] = (unsigned char)((q < 0) ? 0 : ((q > 255) ? 255 : q));
            }
    return 0;
}

#define AFF_NH (256 * KB)
#define AFF_ACC (AFF_NH + 1 + 144 + 12)

typedef struct
{
    const AffLat *A;
    const Tpl *T;
    double Tm[16];       /* input voxel -> template voxel, 0-based */
    float *b;            /* KB*n priors */
    unsigned char *msk1; /* voxel is used with this transform */
    const double *h1;    /* 256 x KB */
    double *acc;         /* per slice: histogram, ll, Alpha, Beta */
    int mode;            /* 0: priors, 1: histogram, 2: derivatives */
} maff_ctx;

static void
maff_range(void *p, int z0, int z1)
{
    maff_ctx *c = (maff_ctx *)p;
    const AffLat *A = c->A;
    const double *Tm = c->Tm;
    int x, y, z, k, r, s;

    for (z = z0; z < z1; z++)
    {
        double *acc = c->acc + (size_t)z * AFF_ACC;
        if (c->mode) memset(acc, 0, sizeof(double) * AFF_ACC);
        for (y = 0; y < A->d[1]; y++)
            for (x = 0; x < A->d[0]; x++)
            {
                size_t i = x + (size_t)A->d[0] * (y + (size_t)A->d[1] * z);
                double xv[3], t[3], a[KB], b[KB], g[KB][3];
                if (!A->msk[i]) continue;
                if (c->mode == 1)
                {
                    const double *h = c->h1 + (size_t)A->g[i] * KB;
                    const float *bi = c->b + KB * i;
                    double q[KB], sq = DBL_EPS;
                    if (!c->msk1[i]) continue;
                    for (k = 0; k < KB; k++)
                    {
                        q[k] = h[k] * bi[k];
                        sq += q[k];
                    }
                    acc[AFF_NH] += log(sq);
                    for (k = 0; k < KB; k++) acc[(size_t)A->g[i] * KB + k] += q[k] / sq;
                    continue;
                }
                xv[0] = (double)x * A->sk[0];
                xv[1] = (double)y * A->sk[1];
                xv[2] = (double)z * A->sk[2];
                for (r = 0; r < 3; r++)
                    t[r] = Tm[r * 4] * xv[0] + Tm[r * 4 + 1] * xv[1] + Tm[r * 4 + 2] * xv[2] + Tm[r * 4 + 3];
                /* nothing below the template is used (the neck) */
                if (t[2] < 0.0)
                {
                    if (c->mode == 0) c->msk1[i] = 0;
                    continue;
                }
                if (c->mode == 0)
                {
                    tpl_sample(c->T, t[0], t[1], t[2], a, NULL);
                    softmax_k(a, NULL, b, NULL);
                    for (k = 0; k < KB; k++) c->b[KB * i + k] = (float)b[k];
                    c->msk1[i] = 1;
                }
                else
                {
                    const double *h = c->h1 + (size_t)A->g[i] * KB;
                    double mi = DBL_EPS, dm[3] = {0, 0, 0}, da[3] = {0, 0, 0}, Av[12];
                    double *Al = acc + AFF_NH + 1, *Be = Al + 144;
                    tpl_sample(c->T, t[0], t[1], t[2], a, g);
                    softmax_k(a, NULL, b, NULL);
                    for (k = 0; k < KB; k++)
                        for (r = 0; r < 3; r++) da[r] += b[k] * g[k][r];
                    for (k = 0; k < KB; k++)
                    {
                        mi += h[k] * b[k];
                        for (r = 0; r < 3; r++) dm[r] += h[k] * b[k] * (g[k][r] - da[r]);
                    }
                    for (r = 0; r < 3; r++) dm[r] /= mi;
                    for (s = 0; s < 3; s++)
                        for (r = 0; r < 3; r++) Av[s * 3 + r] = dm[r] * xv[s];
                    for (r = 0; r < 3; r++) Av[9 + r] = dm[r];
                    for (r = 0; r < 12; r++)
                    {
                        Be[r] -= Av[r];
                        for (s = 0; s < 12; s++) Al[r * 12 + s] += Av[r] * Av[s];
                    }
                }
            }
    }
}

/* Smooth the joint histogram along the intensity and divide it by its
 * marginals: h1 = P(g,k) / (P(g) P(k)). */
static void
maff_smooth_hist(const double *h0, const double *krn, int kr, double *h1)
{
    double tot = 0.0, rs[256], cs[KB];
    int g, k, j;

    for (g = 0; g < AFF_NH; g++) tot += h0[g] + DBL_EPS;
    for (k = 0; k < KB; k++) cs[k] = 0.0;
    for (g = 0; g < 256; g++)
    {
        rs[g] = 0.0;
        for (k = 0; k < KB; k++)
        {
            double s = 0.0;
            for (j = -kr; j <= kr; j++)
                if (g + j >= 0 && g + j < 256)
                    s += (h0[(g + j) * KB + k] + DBL_EPS) / tot * krn[j + kr];
            h1[g * KB + k] = s;
            rs[g] += s;
            cs[k] += s;
        }
    }
    for (g = 0; g < 256; g++)
        for (k = 0; k < KB; k++) h1[g * KB + k] /= rs[g] * cs[k];
}

/* The 12 x 12 derivatives of the voxel-to-voxel matrix (3 x 4, column by
 * column) with respect to the parameters, numerically as in spm_maff8.m. */
static void
maff_derivs(const double itmat[16], const double mat[16], const double P[12], double R[144])
{
    double M[16], M0[16], M1[16], P1[12];
    const double dp = 1e-7;
    int i, r, c;

    aff_P2M(P, M);
    m4_mul(itmat, M, M0);
    m4_mul(M0, mat, M0);
    for (i = 0; i < 12; i++)
    {
        memcpy(P1, P, sizeof(P1));
        P1[i] += dp;
        aff_P2M(P1, M);
        m4_mul(itmat, M, M1);
        m4_mul(M1, mat, M1);
        for (c = 0; c < 4; c++)
            for (r = 0; r < 3; r++)
                R[(c * 3 + r) * 12 + i] = (M1[r * 4 + c] - M0[r * 4 + c]) / dp;
    }
}

/* Affine registration of the image to the template by maximising the mutual
 * information of intensity and class: affreg of spm_maff8.m.  M holds the
 * starting estimate and the result; the objective is returned. */
static double
maff8(const AffLat *A, const Tpl *T, const double mat[16], const double tmat[16], double M[16],
      int regtyp, int maxit, int nthreads, int *err)
{
    double mu[12], isig[36], Alpha0[144], sol[12], sol1[12], dsol[12], itmat[16];
    double *h0, *h1, *acc, krn[25], ll = -HUGE_VAL;
    const int kr = 12, nsearch = 12;
    maff_ctx c;
    int iter, search, subit, i, j, z, done = 0;

    *err = 0;
    affine_prior(regtyp, mu + 6, isig);
    for (i = 0; i < 6; i++) mu[i] = 0.0;
    memset(Alpha0, 0, sizeof(Alpha0));
    for (i = 0; i < 6; i++) Alpha0[i * 13] = 0.00001 * A->ff;
    for (i = 0; i < 6; i++)
        for (j = 0; j < 6; j++) Alpha0[(6 + i) * 12 + 6 + j] = isig[i * 6 + j] * A->ff;
    if (m4_inv(tmat, itmat) || aff_M2P(M, sol))
    {
        *err = 1;
        return ll;
    }
    /* spm_smoothkern(4, x, 0): a Gaussian of 4 bins FWHM, integrated per bin */
    {
        double s = pow(4.0 / sqrt(8.0 * log(2.0)), 2) + DBL_EPS, w1 = 1.0 / sqrt(2.0 * s);
        for (i = -kr; i <= kr; i++)
            krn[i + kr] = 0.5 * (erf(w1 * (i + 0.5)) - erf(w1 * (i - 0.5)));
    }
    h0 = (double *)malloc(sizeof(double) * AFF_NH);
    h1 = (double *)malloc(sizeof(double) * AFF_NH);
    acc = (double *)malloc(sizeof(double) * AFF_ACC * A->d[2]);
    c.b = (float *)malloc(sizeof(float) * KB * A->n);
    c.msk1 = (unsigned char *)malloc(A->n);
    if (!h0 || !h1 || !acc || !c.b || !c.msk1)
    {
        *err = 2;
        goto done;
    }
    c.A = A;
    c.T = T;
    c.h1 = h1;
    c.acc = acc;
    for (i = 0; i < AFF_NH; i++) h1[i] = 1.0;
    memset(dsol, 0, sizeof(dsol));

    for (iter = 1; iter <= maxit && !done; iter++)
    {
        double stepsize = 1.0, ll1 = 0.0, Alpha[144], Beta[12], R[144], RtA[144];
        for (search = 1; search <= nsearch; search++)
        {
            double penalty = 0.0, ll0 = 0.0, ssh = 0.0, Mi[16], sh = 0.0;
            for (i = 0; i < 12; i++) sol1[i] = (iter > 1) ? sol[i] - stepsize * dsol[i] : sol[i];
            for (i = 0; i < 12; i++)
                for (j = 0; j < 12; j++)
                    penalty += 0.5 * (sol1[i] - mu[i]) * Alpha0[i * 12 + j] * (sol1[j] - mu[j]);
            aff_P2M(sol1, Mi);
            m4_mul(itmat, Mi, c.Tm);
            m4_mul(c.Tm, mat, c.Tm);
            c.mode = 0;
            par_for(A->d[2], nthreads, maff_range, &c);

            ll1 = 0.0;
            for (subit = 1; subit <= 32; subit++)
            {
                double lls = 0.0;
                c.mode = 1;
                par_for(A->d[2], nthreads, maff_range, &c);
                for (i = 0; i < AFF_NH; i++) h0[i] = DBL_EPS;
                for (z = 0; z < A->d[2]; z++)
                {
                    const double *az = acc + (size_t)z * AFF_ACC;
                    for (i = 0; i < AFF_NH; i++) h0[i] += az[i];
                    lls += az[AFF_NH];
                }
                maff_smooth_hist(h0, krn, kr, h1);
                if (subit % 4 == 0)
                {
                    ll0 = ll1;
                    ll1 = lls;
                    ssh = 0.0;
                    for (i = 0; i < AFF_NH; i++) ssh += h0[i];
                    if ((ll1 - ll0) / ssh < 1e-5) break;
                }
            }
            ssh = 0.0;
            for (i = 0; i < AFF_NH; i++)
            {
                ssh += h0[i];
                sh += h0[i] * log(h1[i]);
            }
            ll1 = (sh - penalty) / ssh / log(2.0);

            if (iter == 1) break;
            if (fabs(ll1 - ll) < 1e-4)
            {
                done = 1;
                break;
            }
            if (ll1 < ll)
            {
                stepsize *= 0.5;
                if (search == nsearch) done = 2; /* give up, keep sol */
            }
            else
                break;
        }
        if (done == 2) break;
        if (!done || ll1 >= ll)
        {
            ll = ll1;
            memcpy(sol, sol1, sizeof(sol));
        }
        if (done) break;

        /* first and approximate second derivatives for the next step */
        c.mode = 2;
        par_for(A->d[2], nthreads, maff_range, &c);
        memset(Alpha, 0, sizeof(Alpha));
        memset(Beta, 0, sizeof(Beta));
        for (z = 0; z < A->d[2]; z++)
        {
            const double *az = acc + (size_t)z * AFF_ACC + AFF_NH + 1;
            for (i = 0; i < 144; i++) Alpha[i] += az[i];
            for (i = 0; i < 12; i++) Beta[i] += az[144 + i];
        }
        maff_derivs(itmat, mat, sol, R);
        for (i = 0; i < 12; i++)
            for (j = 0; j < 12; j++)
            {
                double s = 0.0;
                for (z = 0; z < 12; z++) s += R[z * 12 + i] * Alpha[z * 12 + j];
                RtA[i * 12 + j] = s;
            }
        for (i = 0; i < 12; i++)
        {
            double s = 0.0;
            for (z = 0; z < 12; z++) s += R[z * 12 + i] * Beta[z];
            dsol[i] = s;
            for (j = 0; j < 12; j++)
            {
                s = 0.0;
                for (z = 0; z < 12; z++) s += RtA[i * 12 + z] * R[z * 12 + j];
                Alpha[i * 12 + j] = s + Alpha0[i * 12 + j];
                dsol[i] += Alpha0[i * 12 + j] * (sol[j] - mu[j]);
            }
        }
        if (solve_dense(Alpha, dsol, 12)) break;
    }
    aff_P2M(sol, M);

done:
    free(h0);
    free(h1);
    free(acc);
    free(c.b);
    free(c.msk1);
    return ll;
}

/* ======================================================================== */
/* The voxels the parameters are estimated from, and the bias field         */
/* ======================================================================== */

typedef struct
{
    int d[3], sk[3];
    size_t n, nm;       /* voxels, and those that are used */
    float *f;           /* intensity of the image the parameters come from */
    float *fseg;        /* intensity of the image that is segmented */
    unsigned char *msk;
    float *b, *lb;      /* KB*n: priors and their logs, class fastest */
    float *lbf, *bf;    /* log bias field and bias field */
    float *wt1, *wt2;   /* weights of the bias update */
    float *ty;          /* 3*n: template voxel (0-based) of every voxel */
    int model;          /* CAT_USEG_BIAS_* */
    int nb[3], nbt;     /* basis functions per axis and in total */
    double *B[3];       /* d[i] x nb[i] */
    double *T, *C;      /* coefficients and the diagonal of their prior */
    double llrb;        /* log prior of the bias + sum of its logs */
    int dims[3];        /* the image the lattice is taken from */
    double vx[3], biasfwhm, biasreg;
} Lattice;

static void
lattice_free(Lattice *L)
{
    int i;
    free(L->f);
    free(L->fseg);
    free(L->msk);
    free(L->b);
    free(L->lb);
    free(L->lbf);
    free(L->bf);
    free(L->wt1);
    free(L->wt2);
    free(L->ty);
    for (i = 0; i < 3; i++) free(L->B[i]);
    free(L->T);
    free(L->C);
    memset(L, 0, sizeof(Lattice));
}

/* The K basis functions of the bias field along an axis of N voxels, at the
 * voxels 0, step, 2 step, ...: B[p*K + k].
 * DCT: spm_dctmtx(N, K, pos).
 * Spline: uniform cubic B-splines with K - 3 knot spans over the axis; four
 * of them are not zero at any position, and they sum to one. */
static void
bias_basis(int model, int N, int K, int npos, int step, double *B)
{
    int p, k;

    if (model == CAT_USEG_BIAS_DCT)
    {
        for (p = 0; p < npos; p++)
        {
            int n = p * step;
            B[(size_t)p * K] = 1.0 / sqrt((double)N);
            for (k = 1; k < K; k++)
                B[(size_t)p * K + k] = sqrt(2.0 / N) * cos(M_PI * (2.0 * n + 1.0) * k / (2.0 * N));
        }
        return;
    }
    memset(B, 0, sizeof(double) * (size_t)npos * K);
    for (p = 0; p < npos; p++)
    {
        int ns = K - 3, j;
        double u = (N > 1) ? (double)(p * step) / (N - 1) * ns : 0.0, t;
        j = (int)floor(u);
        if (j > ns - 1) j = ns - 1;
        t = u - j;
        B[(size_t)p * K + j] = (1 - t) * (1 - t) * (1 - t) / 6.0;
        B[(size_t)p * K + j + 1] = (3 * t * t * t - 6 * t * t + 4) / 6.0;
        B[(size_t)p * K + j + 2] = (-3 * t * t * t + 3 * t * t + 3 * t + 1) / 6.0;
        B[(size_t)p * K + j + 3] = t * t * t / 6.0;
    }
}

/* Number of basis functions per axis, the basis on the lattice and the
 * prior of the coefficients, from L->model, biasfwhm and biasreg.
 *
 * DCT, as spm_preproc8.m: 2 * field of view / FWHM functions, with a prior
 * precision that grows with the frequency like the inverse of a Gaussian.
 *
 * Spline: the same model with basis functions of compact support.  A cubic
 * B-spline with knot spacing h is close to a Gaussian of FWHM 1.36 h, so
 * h = 0.75 FWHM describes fields as smooth as the DCT does, and independent
 * coefficients with the precision biasreg * h^3 (the volume one of them
 * stands for) have the prior that SPM gives a constant field. */
static int
lattice_basis(Lattice *L)
{
    double ff = fudge_factor(0.0, L->vx, L->sk), sd[3], h[3], s;
    int x, y, z, k;

    for (k = 0; k < 3; k++)
    {
        if (L->model == CAT_USEG_BIAS_DCT)
        {
            sd[k] = L->vx[k] * L->dims[k] / L->biasfwhm;
            L->nb[k] = (int)ceil(sd[k] * 2.0);
            if (L->nb[k] < 1) L->nb[k] = 1;
            if (L->nb[k] > L->d[k]) L->nb[k] = L->d[k];
        }
        else
        {
            double len = L->vx[k] * (L->dims[k] - 1);
            int ns = (int)ceil(len / (0.75 * L->biasfwhm));
            if (ns < 1) ns = 1;
            L->nb[k] = ns + 3;
            h[k] = len / ns;
        }
    }
    L->nbt = L->nb[0] * L->nb[1] * L->nb[2];
    L->C = (double *)malloc(sizeof(double) * L->nbt);
    if (!L->C) return 2;
    for (k = 0; k < 3; k++)
    {
        L->B[k] = (double *)malloc(sizeof(double) * (size_t)L->d[k] * L->nb[k]);
        if (!L->B[k]) return 2;
        bias_basis(L->model, L->dims[k], L->nb[k], L->d[k], L->sk[k], L->B[k]);
    }
    s = L->biasreg * ff;
    for (z = 0; z < L->nb[2]; z++)
        for (y = 0; y < L->nb[1]; y++)
            for (x = 0; x < L->nb[0]; x++)
            {
                double c;
                if (L->model == CAT_USEG_BIAS_DCT)
                {
                    double kx = exp(-(double)(x * x) / (sd[0] * sd[0])) / sqrt(L->vx[0]);
                    double ky = exp(-(double)(y * y) / (sd[1] * sd[1])) / sqrt(L->vx[1]);
                    double kz = exp(-(double)(z * z) / (sd[2] * sd[2])) / sqrt(L->vx[2]);
                    c = s / (kx * ky * kz * kx * ky * kz);
                }
                else
                    c = s * h[0] * h[1] * h[2];
                L->C[x + L->nb[0] * (y + L->nb[1] * z)] = c;
            }
    return 0;
}

typedef struct
{
    const Lattice *L;
    const double *T;
    float *lbf, *bf;
    double *sum;      /* per slice */
    double *fd;       /* instead of lbf and bf: the field times wt */
    const float *wt;
} bfield_ctx;

/* The field of the coefficients T on the lattice, slice by slice: the log
 * bias field and the bias field, or (fd) the field times a weight. */
static void
bfield_range(void *p, int z0, int z1)
{
    bfield_ctx *c = (bfield_ctx *)p;
    const Lattice *L = c->L;
    int nb0 = L->nb[0], nb1 = L->nb[1], nb2 = L->nb[2], x, y, z, i, j, k;
    double *t1 = (double *)malloc(sizeof(double) * ((size_t)nb0 * nb1 + (size_t)L->d[0] * nb1));
    double *tmp = t1 + (size_t)nb0 * nb1;

    if (!t1) return;
    for (z = z0; z < z1; z++)
    {
        const double *b2 = L->B[2] + (size_t)z * nb2;
        double sum = 0.0;
        for (i = 0; i < nb0 * nb1; i++)
        {
            double s = 0.0;
            for (k = 0; k < nb2; k++) s += c->T[i + (size_t)nb0 * nb1 * k] * b2[k];
            t1[i] = s;
        }
        for (x = 0; x < L->d[0]; x++)
        {
            const double *b0 = L->B[0] + (size_t)x * nb0;
            for (j = 0; j < nb1; j++)
            {
                double s = 0.0;
                for (i = 0; i < nb0; i++) s += b0[i] * t1[i + nb0 * j];
                tmp[(size_t)x * nb1 + j] = s;
            }
        }
        for (y = 0; y < L->d[1]; y++)
        {
            const double *b1 = L->B[1] + (size_t)y * nb1;
            for (x = 0; x < L->d[0]; x++)
            {
                size_t v = x + (size_t)L->d[0] * (y + (size_t)L->d[1] * z);
                double s = 0.0;
                for (j = 0; j < nb1; j++) s += tmp[(size_t)x * nb1 + j] * b1[j];
                if (c->fd)
                {
                    c->fd[v] = s * c->wt[v];
                    continue;
                }
                c->lbf[v] = (float)s;
                c->bf[v] = (float)exp(s);
                if (L->msk[v]) sum += s;
            }
        }
        if (c->sum) c->sum[z] = sum;
    }
    free(t1);
}

typedef struct
{
    const Lattice *L;
    const double *fd; /* a field on the lattice, or */
    const float *ff;  /* the same as float */
    int square;       /* with the squares of the basis functions */
    double *s2;       /* per slice: nb0*nb1 */
} badj_ctx;

/* The transpose of bfield_range() within the slices: s2 = (B1 x B0)' field */
static void
badj_range(void *p, int z0, int z1)
{
    badj_ctx *c = (badj_ctx *)p;
    const Lattice *L = c->L;
    int nb0 = L->nb[0], nb1 = L->nb[1], m = nb0 * nb1, x, y, z, i, j;
    double *r = (double *)malloc(sizeof(double) * nb0);

    if (!r) return;
    for (z = z0; z < z1; z++)
    {
        double *s2 = c->s2 + (size_t)z * m;
        memset(s2, 0, sizeof(double) * m);
        for (y = 0; y < L->d[1]; y++)
        {
            const double *b1 = L->B[1] + (size_t)y * nb1;
            size_t row = (size_t)L->d[0] * (y + (size_t)L->d[1] * z);
            for (i = 0; i < nb0; i++) r[i] = 0.0;
            for (x = 0; x < L->d[0]; x++)
            {
                const double *b0 = L->B[0] + (size_t)x * nb0;
                double v = c->fd ? c->fd[row + x] : (double)c->ff[row + x];
                if (v == 0.0) continue;
                if (c->square)
                    for (i = 0; i < nb0; i++) r[i] += v * b0[i] * b0[i];
                else
                    for (i = 0; i < nb0; i++) r[i] += v * b0[i];
            }
            for (j = 0; j < nb1; j++)
            {
                double bj = c->square ? b1[j] * b1[j] : b1[j];
                if (bj == 0.0) continue;
                for (i = 0; i < nb0; i++) s2[i + nb0 * j] += bj * r[i];
            }
        }
    }
    free(r);
}

/* g = B' field for a field on the lattice (given as double or as float),
 * B the basis of all three axes; with square, of the squared basis
 * functions (the diagonal of B' diag(field) B).  s2 is scratch of
 * d[2]*nb0*nb1.  The slices are added in their order. */
static void
basis_adjoint(const Lattice *L, const double *fd, const float *ff, int square, double *s2,
              double *g, int nthreads)
{
    badj_ctx c;
    int m = L->nb[0] * L->nb[1], nb2 = L->nb[2], z, k, a;

    c.L = L;
    c.fd = fd;
    c.ff = ff;
    c.square = square;
    c.s2 = s2;
    par_for(L->d[2], nthreads, badj_range, &c);
    for (a = 0; a < L->nbt; a++) g[a] = 0.0;
    for (z = 0; z < L->d[2]; z++)
    {
        const double *sz = s2 + (size_t)z * m;
        for (k = 0; k < nb2; k++)
        {
            double b = L->B[2][(size_t)z * nb2 + k];
            if (square) b *= b;
            if (b == 0.0) continue;
            for (a = 0; a < m; a++) g[a + m * k] += b * sz[a];
        }
    }
}

/* Evaluate the bias field of L->T and its part of the log-likelihood. */
static int
lattice_bias(Lattice *L, int nthreads)
{
    bfield_ctx c;
    double *sum, ll = 0.0;
    size_t i;
    int z;

    if (L->nbt == 0)
    {
        for (i = 0; i < L->n; i++)
        {
            L->lbf[i] = 0.0f;
            L->bf[i] = 1.0f;
        }
        L->llrb = 0.0;
        return 0;
    }
    sum = (double *)calloc(L->d[2], sizeof(double));
    if (!sum) return 2;
    c.L = L;
    c.T = L->T;
    c.lbf = L->lbf;
    c.bf = L->bf;
    c.sum = sum;
    c.fd = NULL;
    c.wt = NULL;
    par_for(L->d[2], nthreads, bfield_range, &c);
    for (z = 0; z < L->d[2]; z++) ll += sum[z];
    for (i = 0; i < (size_t)L->nbt; i++) ll -= 0.5 * L->C[i] * L->T[i] * L->T[i];
    L->llrb = ll;
    free(sum);
    return 0;
}

/* The lattice of voxels samp mm apart with the image est, and clean, the one
 * that is segmented; voxels without data in clean, masked out, or below the
 * template take no part. */
static int
lattice_build(Lattice *L, const float *est, const float *clean, const unsigned char *mask,
              const int dims[3], const double vx[3], const double Tvox[16], double zmin,
              double samp, int model, double biasfwhm, double biasreg, int nthreads)
{
    size_t nxy = (size_t)dims[0] * dims[1], i;
    int x, y, z, k, ret;

    memset(L, 0, sizeof(Lattice));
    L->model = model;
    L->biasfwhm = biasfwhm;
    L->biasreg = biasreg;
    for (k = 0; k < 3; k++)
    {
        L->dims[k] = dims[k];
        L->vx[k] = vx[k];
        L->sk[k] = (int)floor(samp / vx[k] + 0.5);
        if (L->sk[k] < 1) L->sk[k] = 1;
        L->d[k] = (dims[k] - 1) / L->sk[k] + 1;
    }
    L->n = (size_t)L->d[0] * L->d[1] * L->d[2];
    L->f = (float *)calloc(L->n, sizeof(float));
    L->fseg = (float *)calloc(L->n, sizeof(float));
    L->msk = (unsigned char *)calloc(L->n, 1);
    L->b = (float *)malloc(sizeof(float) * KB * L->n);
    L->lb = (float *)malloc(sizeof(float) * KB * L->n);
    L->lbf = (float *)malloc(sizeof(float) * L->n);
    L->bf = (float *)malloc(sizeof(float) * L->n);
    L->wt1 = (float *)malloc(sizeof(float) * L->n);
    L->wt2 = (float *)malloc(sizeof(float) * L->n);
    L->ty = (float *)malloc(sizeof(float) * 3 * L->n);
    if (!L->f || !L->fseg || !L->msk || !L->b || !L->lb || !L->lbf || !L->bf || !L->wt1 || !L->wt2 ||
        !L->ty)
        return 2;

    /* only voxels above zmin, in template voxels along z, are used: SPM
     * leaves out everything less than 5 mm above the bottom of its priors,
     * assuming that the affine registration is close */
    for (z = 0; z < L->d[2]; z++)
        for (y = 0; y < L->d[1]; y++)
            for (x = 0; x < L->d[0]; x++)
            {
                double xv = (double)x * L->sk[0], yv = (double)y * L->sk[1], zv = (double)z * L->sk[2];
                size_t v = (size_t)xv + (size_t)dims[0] * (size_t)yv + nxy * (size_t)zv;
                double tz = Tvox[8] * xv + Tvox[9] * yv + Tvox[10] * zv + Tvox[11];
                i = x + (size_t)L->d[0] * (y + (size_t)L->d[1] * z);
                if (!is_valid(clean[v]) || (mask && !mask[v]) || !(tz > zmin)) continue;
                L->msk[i] = 1;
                L->f[i] = est[v];
                L->fseg[i] = clean[v];
                L->nm++;
            }
    if (L->nm < 1000) return 3;

    /* the basis of the bias field */
    L->nbt = 0;
    if (biasfwhm > 0)
    {
        if ((ret = lattice_basis(L))) return ret;
        L->T = (double *)calloc(L->nbt, sizeof(double));
        if (!L->T) return 2;
    }
    return lattice_bias(L, nthreads);
}


/* ======================================================================== */
/* Mixture of Gaussians                                                     */
/* ======================================================================== */

typedef struct
{
    int K;
    int lkp[MAXK];
    double mg[MAXK], mn[MAXK], vr[MAXK];
    double wp[KB];
    double vr0;   /* Wishart-style prior of the variances */
    double vrmin; /* numerical floor of a variance */
    double wp_reg; /* keeps the class weights near 1/KB, see mstep() */
} Gmm;

#define ACC_N (1 + 3 * MAXK + 2 * KB)

typedef struct
{
    Lattice *L;
    const Gmm *G;
    int mode;    /* 0: log-likelihood, 1: + moments, 2: + weights of the bias update */
    double *acc; /* per slice: ll, mom0, mom1, mom2, mgm, prior mass */
    double cg[MAXK], hg[MAXK], lwp[KB];
} estep_ctx;

/* The E-step over the lattice: responsibilities of the Gaussians from the
 * priors, the class weights and the bias-corrected intensity (latent() of
 * spm_preproc8.m), with what the M-steps need of them. */
static void
estep_range(void *p, int z0, int z1)
{
    estep_ctx *c = (estep_ctx *)p;
    Lattice *L = c->L;
    const Gmm *G = c->G;
    size_t nxy = (size_t)L->d[0] * L->d[1], i;
    int z, k, g, K = G->K;

    for (z = z0; z < z1; z++)
    {
        double *acc = c->acc + (size_t)z * ACC_N;
        double *m0 = acc + 1, *m1 = m0 + MAXK, *m2 = m1 + MAXK, *mgm = m2 + MAXK;
        double *pm = mgm + KB;
        memset(acc, 0, sizeof(double) * ACC_N);
        for (i = nxy * z; i < nxy * (z + 1); i++)
        {
            const float *b = L->b + KB * i, *lb = L->lb + KB * i;
            double cr, s = 0.0, ls, lq[MAXK], q[MAXK], mx = -HUGE_VAL, sp = 0.0;
            if (!L->msk[i])
            {
                if (c->mode == 2) L->wt1[i] = L->wt2[i] = 0.0f;
                continue;
            }
            cr = (double)L->f[i] * L->bf[i];
            for (k = 0; k < KB; k++) s += G->wp[k] * b[k];
            ls = log(s);
            for (g = 0; g < K; g++)
            {
                double d = cr - G->mn[g];
                lq[g] = c->cg[g] - c->hg[g] * d * d + c->lwp[G->lkp[g]] + lb[G->lkp[g]] - ls;
                if (lq[g] > mx) mx = lq[g];
            }
            for (g = 0; g < K; g++)
            {
                q[g] = exp(lq[g] - mx);
                sp += q[g];
            }
            acc[0] += mx + log(sp);
            if (c->mode == 0) continue;
            if (c->mode == 1)
            {
                for (k = 0; k < KB; k++)
                {
                    mgm[k] += b[k] / s;
                    pm[k] += b[k];
                }
                for (g = 0; g < K; g++)
                {
                    double qg = q[g] / sp;
                    m0[g] += qg;
                    m1[g] += qg * cr;
                    m2[g] += qg * cr * cr;
                }
            }
            else
            {
                double w1 = 0.0, w2 = 0.0;
                for (g = 0; g < K; g++)
                {
                    double qg = q[g] / sp;
                    w1 += qg * (G->mn[g] - cr) / G->vr[g];
                    w2 += qg / G->vr[g];
                }
                L->wt1[i] = (float)(-(1.0 + cr * w1));
                L->wt2[i] = (float)(cr * cr * w2 + 1.0);
            }
        }
    }
}

/* Run the E-step; returns the log-likelihood including the bias terms.
 * mode 1 leaves the moments summed over the slices in acc[0..ACC_N). */
static double
estep(Lattice *L, const Gmm *G, int mode, double *acc, int nthreads)
{
    estep_ctx c;
    double ll = 0.0;
    int z, g, i;

    c.L = L;
    c.G = G;
    c.mode = mode;
    c.acc = acc;
    for (g = 0; g < G->K; g++)
    {
        c.cg[g] = log(G->mg[g]) - 0.5 * log(2.0 * M_PI * G->vr[g]);
        c.hg[g] = 0.5 / G->vr[g];
    }
    for (g = 0; g < KB; g++) c.lwp[g] = log(G->wp[g]);
    par_for(L->d[2], nthreads, estep_range, &c);
    for (z = 1; z < L->d[2]; z++)
    {
        const double *az = acc + (size_t)z * ACC_N;
        if (mode == 1)
            for (i = 0; i < ACC_N; i++) acc[i] += az[i];
        else
            acc[0] += az[0];
    }
    ll = acc[0] + L->llrb;
    return ll;
}

/* The M-step of the Gaussians and of the class weights (US eq. 23-27).
 *
 * A class weight w_k scales the prior of its class, b_k -> w_k b_k / sum.
 * SPM estimates them freely: w_k = N_k / sum_i (b_ik / sum_j w_j b_ij), with
 * N_k the voxels the class got.  That is right for priors made with the same
 * model, such as SPM's tissue probability map.  The shooting templates stem
 * from CAT's segmentations, which give white matter what a single Gaussian
 * leaves to grey matter (partial volumes); free weights then turn white
 * matter down by a factor of three to four against grey matter, and the
 * deformation follows.  So the weights are held near 1/KB by wp_reg times
 * the prior mass M_k of every class as pseudo-counts, which leaves room for
 * a brain with more CSF than the template without giving that freedom to the
 * mismatch of the two models: w_k = (N_k + r M_k) / (mgm_k + r M_k KB). */
static void
mstep(Gmm *G, const double *acc)
{
    const double *m0 = acc + 1, *m1 = m0 + MAXK, *m2 = m1 + MAXK, *mgm = m2 + MAXK;
    const double *pm = mgm + KB;
    double cls[KB], s = 0.0;
    int g, k;

    for (k = 0; k < KB; k++) cls[k] = 0.0;
    for (g = 0; g < G->K; g++) cls[G->lkp[g]] += m0[g] + 2.0 * MOM_TINY;
    for (g = 0; g < G->K; g++)
    {
        double n0 = m0[g] + MOM_TINY;
        G->mg[g] = (n0 + MOM_TINY) / cls[G->lkp[g]];
        G->mn[g] = m1[g] / (n0 + MOM_TINY);
        G->vr[g] = (m2[g] - m1[g] * m1[g] / n0 + G->vr0) / (n0 + 1.0);
        if (!(G->vr[g] > G->vrmin)) G->vr[g] = G->vrmin;
    }
    for (k = 0; k < KB; k++)
    {
        double r = G->wp_reg * pm[k] + 1.0; /* + SPM's wp_reg of one voxel */
        G->wp[k] = (cls[k] + r) / (mgm[k] + r * KB);
        s += G->wp[k];
    }
    for (k = 0; k < KB; k++) G->wp[k] /= s;
}

/* One Gaussian per class from the moments under the priors alone. */
static void
gmm_init(Gmm *G, const Lattice *L)
{
    double mm0[KB], mm1[KB], mm2[KB], s0 = 0.0, s1 = 0.0, s2 = 0.0, vr1 = 0.0, sm = 0.0;
    size_t i;
    int k;

    for (k = 0; k < KB; k++) mm0[k] = mm1[k] = mm2[k] = 0.0;
    for (i = 0; i < L->n; i++)
    {
        const float *b = L->b + KB * i;
        double cr, f;
        if (!L->msk[i]) continue;
        f = L->f[i];
        cr = f * L->bf[i];
        s0 += 1.0;
        s1 += f;
        s2 += f * f;
        for (k = 0; k < KB; k++)
        {
            mm0[k] += b[k];
            mm1[k] += b[k] * cr;
            mm2[k] += b[k] * cr * cr;
        }
    }
    G->vr0 = (s2 / s0 - (s1 / s0) * (s1 / s0)) / (KB * KB);
    G->vrmin = 1e-8 * (s2 / s0 - (s1 / s0) * (s1 / s0));
    G->K = KB;
    for (k = 0; k < KB; k++)
    {
        G->lkp[k] = k;
        G->mg[k] = 1.0;
        G->mn[k] = mm1[k] / (mm0[k] + MOM_TINY);
        G->wp[k] = 1.0 / KB;
        vr1 += mm2[k] - mm1[k] * mm1[k] / (mm0[k] + MOM_TINY);
        sm += mm0[k];
    }
    vr1 = (vr1 + G->vr0) / (sm + 1.0);
    for (k = 0; k < KB; k++) G->vr[k] = vr1;
}

/* Replace the single Gaussian of every class by ngaus of them: SPM's crude
 * heuristic, with the means spread over the quantiles of the Gaussian they
 * replace instead of random deviates, so that the result is reproducible. */
static void
gmm_split(Gmm *G, const int ngaus[KB])
{
    Gmm N = *G;
    int k, j, K = 0;

    for (k = 0; k < KB; k++)
    {
        int kk = ngaus[k];
        double w = 1.0 / (1.0 + exp(-(kk - 1) * 0.25)) - 0.5;
        for (j = 0; j < kk; j++)
        {
            double zq = (kk > 1) ? norm_inv((j + 0.5) / kk) : 0.0;
            N.lkp[K] = k;
            N.mn[K] = sqrt(G->vr[k]) * zq * w + G->mn[k];
            N.vr[K] = G->vr[k] * (1.0 - w);
            N.mg[K] = 1.0 / kk;
            K++;
        }
    }
    N.K = K;
    *G = N;
}

/* ======================================================================== */
/* Gauss-Newton update of the bias field                                    */
/* ======================================================================== */

/* ======================================================================== */
/* State of a segmentation                                                  */
/* ======================================================================== */

typedef struct
{
    const float *orig;  /* the image as it came: affine registration */
    const float *clean; /* the image that is segmented (filtered if asked for) */
    const float *est;   /* the image the parameters are estimated from */
    const unsigned char *mask;
    int dims[3];
    double mat[16], vx[3];
    int tdims[3];
    double tmat[16];
    double Tvox[16]; /* input voxel -> template voxel under the affine, 0-based */
    CAT_UnifiedSegOpts o;
    int nthreads;
    Tpl tpl;
    Lattice L;
    Gmm G;
    double ll;
    int initialised;
    double *acc;     /* E-step accumulators */
    float *ireg;     /* image on the registration grid */
    size_t nreg;
    double range[2];
    double t_em, t_bias, t_prior, t_fill; /* seconds spent in the model */
    int n_bias, n_cg;                     /* bias updates and their CG iterations */
} State;

/* H p = B' diag(wt2) B p + C p: the approximate Hessian of the bias update
 * applied to a vector, without forming it. */
static void
bias_hessian(const Lattice *L, const double *pv, double *field, double *s2, double *Hp,
             int nthreads)
{
    bfield_ctx c;
    int i;

    c.L = L;
    c.T = pv;
    c.lbf = c.bf = NULL;
    c.sum = NULL;
    c.fd = field;
    c.wt = L->wt2;
    par_for(L->d[2], nthreads, bfield_range, &c);
    basis_adjoint(L, field, NULL, 0, s2, Hp, nthreads);
    for (i = 0; i < L->nbt; i++) Hp[i] += L->C[i] * pv[i];
}

/* One Gauss-Newton step of the bias field with a line search (US eq. 34):
 * T -= (B' W2 B + C) \ (B' w1 + C T), with the weights w1 and W2 of the
 * E-step.  SPM forms the matrix and factorises it, which is cubic in the
 * number of coefficients; here the system is solved by conjugate gradients
 * with the diagonal as preconditioner, and the matrix is only ever applied,
 * through the separable basis. */
static int
bias_step(State *S)
{
    Lattice *L = &S->L;
    int n = L->nbt, i, it, search, ret = 0, ok = 0;
    double *w, *x, *r, *z, *pv, *Hp, *dg, *oldT, *field, *s2;
    double armijo = 1.0, ll_prev, rz, rr0, rr;

    if (n == 0) return 0;
    w = (double *)malloc(sizeof(double) * 7 * (size_t)n);
    field = (double *)malloc(sizeof(double) * L->n);
    s2 = (double *)malloc(sizeof(double) * (size_t)L->d[2] * L->nb[0] * L->nb[1]);
    if (!w || !field || !s2)
    {
        ret = 2;
        goto done;
    }
    x = w;
    r = x + n;
    z = r + n;
    pv = z + n;
    Hp = pv + n;
    dg = Hp + n;
    oldT = dg + n;

    ll_prev = estep(L, &S->G, 2, S->acc, S->nthreads);
    basis_adjoint(L, NULL, L->wt1, 0, s2, r, S->nthreads);
    basis_adjoint(L, NULL, L->wt2, 1, s2, dg, S->nthreads);
    rr0 = 0.0;
    for (i = 0; i < n; i++)
    {
        r[i] += L->C[i] * L->T[i];
        dg[i] += L->C[i];
        x[i] = 0.0;
        z[i] = r[i] / dg[i];
        pv[i] = z[i];
        rr0 += r[i] * r[i];
    }
    rz = 0.0;
    for (i = 0; i < n; i++) rz += r[i] * z[i];
    /* a Gauss-Newton step with a line search does not need the system solved
     * to more than three digits */
    for (it = 0; it < 100 && rr0 > 0; it++)
    {
        double pHp = 0.0, alpha, rz1 = 0.0;
        bias_hessian(L, pv, field, s2, Hp, S->nthreads);
        for (i = 0; i < n; i++) pHp += pv[i] * Hp[i];
        if (!(pHp > 0)) break;
        alpha = rz / pHp;
        rr = 0.0;
        for (i = 0; i < n; i++)
        {
            x[i] += alpha * pv[i];
            r[i] -= alpha * Hp[i];
            rr += r[i] * r[i];
        }
        if (rr < 1e-6 * rr0) break;
        for (i = 0; i < n; i++)
        {
            z[i] = r[i] / dg[i];
            rz1 += r[i] * z[i];
        }
        for (i = 0; i < n; i++) pv[i] = z[i] + (rz1 / rz) * pv[i];
        rz = rz1;
    }

    memcpy(oldT, L->T, sizeof(double) * n);
    for (search = 0; search < 12; search++)
    {
        double ll;
        for (i = 0; i < n; i++) L->T[i] = oldT[i] - armijo * x[i];
        if ((ret = lattice_bias(L, S->nthreads))) goto done;
        ll = estep(L, &S->G, 0, S->acc, S->nthreads);
        if (ll >= ll_prev)
        {
            S->ll = ll;
            ok = 1;
            break;
        }
        armijo *= 0.5;
    }
    if (!ok)
    {
        memcpy(L->T, oldT, sizeof(double) * n);
        ret = lattice_bias(L, S->nthreads);
    }
    S->n_bias++;
    S->n_cg += it;
    if (S->o.verbose > 1)
        printf("    bias: ll %.6f, step %g%s, %d CG iterations\n", S->ll / (double)L->nm, armijo,
               ok ? "" : " (rejected)", it);

done:
    free(w);
    free(field);
    free(s2);
    return ret;
}

/* Alternate the EM of the Gaussians with updates of the bias field, at most
 * n_inner times: the inner loop of spm_preproc8.m.  With split the Gaussians
 * are those of gmm_init() and get their final number after the first pass. */
static int
estimate(State *S, int n_inner, int split)
{
    Lattice *L = &S->L;
    double tol = S->o.tol * (double)L->nm, ll = -HUGE_VAL, ooll = -HUGE_VAL;
    int iter1, subit, ret;

    for (iter1 = 1; iter1 <= n_inner; iter1++)
    {
        double t0 = wall_time();
        for (subit = 1; subit <= 20; subit++)
        {
            double ll_prev = ll;
            ll = estep(L, &S->G, 1, S->acc, S->nthreads);
            mstep(&S->G, S->acc);
            if (!isfinite(ll)) return 3;
            if (subit > 2 && ll - ll_prev < tol) break;
        }
        S->t_em += wall_time() - t0;
        S->ll = ll;
        if (S->o.verbose > 1)
            printf("    EM %d: ll %.6f (%d steps)\n", iter1, ll / (double)L->nm, subit > 20 ? 20 : subit);
        /* the last pass ends with the Gaussians, which then belong to the
         * bias field that is kept */
        if (iter1 == n_inner) break;
        if (iter1 > 1 && !(ll - ooll > 2.0 * tol) && !(split && iter1 == 2)) break;
        ooll = ll;
        t0 = wall_time();
        if ((ret = bias_step(S))) return ret;
        S->t_bias += wall_time() - t0;
        ll = S->ll;
        if (split && iter1 == 1) gmm_split(&S->G, S->o.ngaus);
    }
    return 0;
}

/* ======================================================================== */
/* Priors on the lattice                                                    */
/* ======================================================================== */

typedef struct
{
    State *S;
} prior_ctx;

/* The priors at the lattice voxels: the template at L->ty, or where that is
 * not defined at the position under the affine transform. */
static void
prior_range(void *p, int z0, int z1)
{
    State *S = ((prior_ctx *)p)->S;
    Lattice *L = &S->L;
    const double *T = S->Tvox;
    int x, y, z, k, r;

    for (z = z0; z < z1; z++)
        for (y = 0; y < L->d[1]; y++)
            for (x = 0; x < L->d[0]; x++)
            {
                size_t i = x + (size_t)L->d[0] * (y + (size_t)L->d[1] * z);
                double t[3], a[KB], b[KB], lb[KB];
                if (!L->msk[i]) continue;
                for (r = 0; r < 3; r++) t[r] = L->ty[r * L->n + i];
                if (!isfinite(t[0]) || !isfinite(t[1]) || !isfinite(t[2]))
                {
                    double xv = (double)x * L->sk[0], yv = (double)y * L->sk[1], zv = (double)z * L->sk[2];
                    for (r = 0; r < 3; r++)
                        t[r] = T[r * 4] * xv + T[r * 4 + 1] * yv + T[r * 4 + 2] * zv + T[r * 4 + 3];
                }
                tpl_sample(&S->tpl, t[0], t[1], t[2], a, NULL);
                softmax_k(a, NULL, b, lb);
                for (k = 0; k < KB; k++)
                {
                    L->b[KB * i + k] = (float)b[k];
                    L->lb[KB * i + k] = (float)lb[k];
                }
            }
}

static void
lattice_priors(State *S)
{
    prior_ctx c;
    c.S = S;
    par_for(S->L.d[2], S->nthreads, prior_range, &c);
}

/* L->ty from the deformation of a level: the inverse of y on the lattice.
 * y maps registration voxels on the template side to those on the input
 * side, to_input the latter to input voxels, to_template the former to
 * template voxels (all 1-based). */
static void
lattice_invdef(State *S, const CAT_VolShootLevel *lv)
{
    Lattice *L = &S->L;
    double Lm[16], M1d[16];
    float M1[4][3], M2[4][3];
    mwSize dy[3], dl[3];
    size_t i;
    int r, c;

    /* input voxel (1-based) -> lattice index (1-based): (v - 1)/sk + 1 */
    m4_eye(Lm);
    for (r = 0; r < 3; r++)
    {
        Lm[r * 5] = 1.0 / L->sk[r];
        Lm[r * 4 + 3] = 1.0 - 1.0 / L->sk[r];
        dy[r] = (mwSize)lv->dims[r];
        dl[r] = (mwSize)L->d[r];
    }
    m4_mul(Lm, lv->to_input, M1d);
    for (r = 0; r < 3; r++)
        for (c = 0; c < 4; c++)
        {
            M1[c][r] = (float)M1d[r * 4 + c];
            M2[c][r] = (float)lv->to_template[r * 4 + c];
        }
    memset(L->ty, 0, sizeof(float) * 3 * L->n);
    spm_invdef(dy, (float *)lv->y, dl, L->ty, M1, M2);
    for (i = 0; i < 3 * L->n; i++) L->ty[i] -= 1.0f;
}

/* ======================================================================== */
/* The data model of the shooting                                           */
/* ======================================================================== */

typedef struct
{
    State *S;
    const CAT_VolShootLevel *lv;
    float *f;
    double cg[MAXK], hg[MAXK];
} reg_ctx;

/* The image on the registration grid of a level, trilinearly from the
 * voxels that carry data; NaN where fewer than half of the weights do. */
static void
regimage_range(void *p, int z0, int z1)
{
    reg_ctx *c = (reg_ctx *)p;
    State *S = c->S;
    const int *rd = c->lv->dims, *d = S->dims;
    const double *M = c->lv->to_input;
    size_t nxy = (size_t)d[0] * d[1];
    int x, y, z, r, q;

    for (z = z0; z < z1; z++)
        for (y = 0; y < rd[1]; y++)
            for (x = 0; x < rd[0]; x++)
            {
                size_t v = x + (size_t)rd[0] * (y + (size_t)rd[1] * z);
                double s[3], fr[3], sum = 0.0, wsum = 0.0;
                int i0[3], in = 1;
                for (r = 0; r < 3; r++)
                {
                    s[r] = M[r * 4] * (x + 1) + M[r * 4 + 1] * (y + 1) + M[r * 4 + 2] * (z + 1) +
                           M[r * 4 + 3] - 1.0;
                    if (!(s[r] >= -0.5 && s[r] <= d[r] - 0.5)) in = 0;
                }
                S->ireg[v] = NAN;
                if (!in) continue;
                for (r = 0; r < 3; r++)
                {
                    if (s[r] < 0) s[r] = 0;
                    if (s[r] > d[r] - 1) s[r] = d[r] - 1;
                    i0[r] = (int)floor(s[r]);
                    if (i0[r] > d[r] - 2) i0[r] = (d[r] > 1) ? d[r] - 2 : 0;
                    fr[r] = s[r] - i0[r];
                }
                for (q = 0; q < 8; q++)
                {
                    int ox = q & 1, oy = (q >> 1) & 1, oz = (q >> 2) & 1;
                    double w = (ox ? fr[0] : 1 - fr[0]) * (oy ? fr[1] : 1 - fr[1]) * (oz ? fr[2] : 1 - fr[2]);
                    size_t j;
                    if (w <= 0) continue;
                    if ((ox && d[0] < 2) || (oy && d[1] < 2) || (oz && d[2] < 2)) continue;
                    j = (size_t)(i0[0] + ox) + (size_t)d[0] * (i0[1] + oy) + nxy * (i0[2] + oz);
                    if (!is_valid(S->clean[j]) || (S->mask && !S->mask[j])) continue;
                    sum += w * S->est[j];
                    wsum += w;
                }
                if (wsum > 0.5) S->ireg[v] = (float)(sum / wsum);
            }
}

/* The likelihoods of the classes on the registration grid, summing to one. */
static void
regfill_range(void *p, int z0, int z1)
{
    reg_ctx *c = (reg_ctx *)p;
    State *S = c->S;
    const Lattice *L = &S->L;
    const Gmm *G = &S->G;
    const int *rd = c->lv->dims;
    const double *M = c->lv->to_input;
    size_t n = (size_t)rd[0] * rd[1] * rd[2];
    int x, y, z, r, q, k, g;

    for (z = z0; z < z1; z++)
        for (y = 0; y < rd[1]; y++)
            for (x = 0; x < rd[0]; x++)
            {
                size_t v = x + (size_t)rd[0] * (y + (size_t)rd[1] * z);
                double l[3], fr[3], lbf = 0.0, cr, lq[MAXK], mx = -HUGE_VAL, pk[KB], sp = 0.0;
                int i0[3];
                if (!isfinite(S->ireg[v]))
                {
                    for (k = 0; k < KB; k++) c->f[k * n + v] = 1.0f / KB;
                    continue;
                }
                /* the bias field, trilinearly from the lattice */
                for (r = 0; r < 3; r++)
                {
                    l[r] = (M[r * 4] * (x + 1) + M[r * 4 + 1] * (y + 1) + M[r * 4 + 2] * (z + 1) +
                            M[r * 4 + 3] - 1.0) / L->sk[r];
                    if (l[r] < 0) l[r] = 0;
                    if (l[r] > L->d[r] - 1) l[r] = L->d[r] - 1;
                    i0[r] = (int)floor(l[r]);
                    if (i0[r] > L->d[r] - 2) i0[r] = (L->d[r] > 1) ? L->d[r] - 2 : 0;
                    fr[r] = l[r] - i0[r];
                }
                for (q = 0; q < 8; q++)
                {
                    int ox = q & 1, oy = (q >> 1) & 1, oz = (q >> 2) & 1;
                    double w = (ox ? fr[0] : 1 - fr[0]) * (oy ? fr[1] : 1 - fr[1]) * (oz ? fr[2] : 1 - fr[2]);
                    if (w <= 0) continue;
                    if ((ox && L->d[0] < 2) || (oy && L->d[1] < 2) || (oz && L->d[2] < 2)) continue;
                    lbf += w * L->lbf[(size_t)(i0[0] + ox) + (size_t)L->d[0] *
                                      ((i0[1] + oy) + (size_t)L->d[1] * (i0[2] + oz))];
                }
                cr = (double)S->ireg[v] * exp(lbf);
                for (g = 0; g < G->K; g++)
                {
                    double dd = cr - G->mn[g];
                    lq[g] = c->cg[g] - c->hg[g] * dd * dd;
                    if (lq[g] > mx) mx = lq[g];
                }
                for (k = 0; k < KB; k++) pk[k] = 0.0;
                for (g = 0; g < G->K; g++) pk[G->lkp[g]] += exp(lq[g] - mx);
                for (k = 0; k < KB; k++)
                {
                    pk[k] += 1e-6; /* keeps the posterior of a class off zero */
                    sp += pk[k];
                }
                for (k = 0; k < KB; k++) c->f[k * n + v] = (float)(pk[k] / sp);
            }
}

static void print_model(const State *S);

/* CAT_VolShootModel.update: the priors under the current deformation, the
 * parameters of the image given these priors, and from them the class
 * likelihoods the next Gauss-Newton step registers the templates to. */
static int
model_update(void *ctx, const CAT_VolShootLevel *lv, float *f, double *logw)
{
    State *S = (State *)ctx;
    size_t n = (size_t)lv->dims[0] * lv->dims[1] * lv->dims[2];
    reg_ctx c;
    double t0;
    int ret, g, k;

    c.S = S;
    c.lv = lv;
    c.f = f;
    if (lv->start)
    {
        tpl_load(&S->tpl, lv->level, S->nthreads);
        if (n > S->nreg)
        {
            free(S->ireg);
            S->ireg = (float *)malloc(sizeof(float) * n);
            if (!S->ireg) return 2;
            S->nreg = n;
        }
        par_for(lv->dims[2], S->nthreads, regimage_range, &c);
    }
    t0 = wall_time();
    lattice_invdef(S, lv);
    lattice_priors(S);
    S->t_prior += wall_time() - t0;
    if (!S->initialised)
    {
        gmm_init(&S->G, &S->L);
        S->G.wp_reg = S->o.wp_reg;
        ret = estimate(S, 16, 1);
        S->initialised = 1;
    }
    else
        ret = estimate(S, 2, 0);
    if (ret) return ret;
    if (S->o.verbose > 1)
        printf("  level %d, iteration %d: ll %.6f\n", lv->level, lv->iteration, S->ll / (double)S->L.nm);
    if (S->o.verbose > 1 && lv->start) print_model(S);

    for (g = 0; g < S->G.K; g++)
    {
        c.cg[g] = log(S->G.mg[g]) - 0.5 * log(2.0 * M_PI * S->G.vr[g]);
        c.hg[g] = 0.5 / S->G.vr[g];
    }
    t0 = wall_time();
    par_for(lv->dims[2], S->nthreads, regfill_range, &c);
    S->t_fill += wall_time() - t0;
    for (k = 0; k < KB; k++) logw[k] = log(S->G.wp[k]);
    return 0;
}

/* ======================================================================== */
/* Partial volume model                                                     */
/* ======================================================================== */

/* Tables over the bias-corrected intensity of a model in which a voxel is
 * pure CSF, GM or WM, or a mixture of two of them with any fraction t: the
 * likelihood of each of these six cases and the first two moments of the
 * fraction of the second tissue of every mixture.  A mixture with fraction t
 * has the mean (1-t) m_a + t m_b; the variance is that of a pure tissue. */
typedef struct
{
    int ok;
    double m[3], v; /* CSF, GM, WM */
    double y0, dy;
    float *dens;     /* PVE_NTAB x 6: C, G, W, CG, GW, CW */
    float *frac;    /* PVE_NTAB x 6: E[t] and E[t^2] of CG, GW, CW */
} Pve;

static const int pve_pair[3][2] = {{0, 1}, {1, 2}, {0, 2}};
static const int pve_cls[3] = {C_CSF, C_GM, C_WM};

static int
pve_tables(Pve *P)
{
    double lo = P->m[0], hi = P->m[0], sd = sqrt(P->v);
    int i, k, t;

    for (k = 1; k < 3; k++)
    {
        if (P->m[k] < lo) lo = P->m[k];
        if (P->m[k] > hi) hi = P->m[k];
    }
    P->y0 = lo - 5.0 * sd;
    P->dy = (hi - lo + 10.0 * sd) / (PVE_NTAB - 1);
    if (!P->dens) P->dens = (float *)malloc(sizeof(float) * PVE_NTAB * 6);
    if (!P->frac) P->frac = (float *)malloc(sizeof(float) * PVE_NTAB * 6);
    if (!P->dens || !P->frac) return 2;
    for (i = 0; i < PVE_NTAB; i++)
    {
        double y = P->y0 + i * P->dy;
        for (k = 0; k < 3; k++)
        {
            double d = y - P->m[k], s0 = 0.0, s1 = 0.0, s2 = 0.0;
            const int *pr = pve_pair[k];
            P->dens[i * 6 + k] = (float)(exp(-0.5 * d * d / P->v) / sqrt(2.0 * M_PI * P->v));
            for (t = 0; t < PVE_NINT; t++)
            {
                double tt = (t + 0.5) / PVE_NINT;
                double dm = y - ((1.0 - tt) * P->m[pr[0]] + tt * P->m[pr[1]]);
                double e = exp(-0.5 * dm * dm / P->v);
                s0 += e;
                s1 += tt * e;
                s2 += tt * tt * e;
            }
            P->dens[i * 6 + 3 + k] = (float)(s0 / PVE_NINT / sqrt(2.0 * M_PI * P->v));
            P->frac[i * 6 + k] = (float)((s0 > 0) ? s1 / s0 : 0.5);
            P->frac[i * 6 + 3 + k] = (float)((s0 > 0) ? s2 / s0 : 1.0 / 3.0);
        }
    }
    return 0;
}

/* Posterior of the six cases for a brain voxel with intensity y and the
 * priors pr of the three tissues (summing to one), with the moments t1 and
 * t2 of the fraction of every mixture.  Two halves of a voxel drawn
 * independently from the priors give a pure voxel of tissue a with
 * probability pr_a^2 and a mixture of a and b with 2 pr_a pr_b.  Returns 0 if
 * no case has a likelihood (y is beyond every tissue). */
static int
pve_posterior(const Pve *P, double y, const double pr[3], double post[6], double t1[3],
              double t2[3])
{
    double pos = (y - P->y0) / P->dy, w, s = 0.0;
    int i, k;

    if (pos < 0) pos = 0;
    if (pos > PVE_NTAB - 1) pos = PVE_NTAB - 1;
    i = (int)floor(pos);
    if (i > PVE_NTAB - 2) i = PVE_NTAB - 2;
    w = pos - i;
    for (k = 0; k < 3; k++)
    {
        double lp = (1 - w) * P->dens[i * 6 + k] + w * P->dens[(i + 1) * 6 + k];
        double lm = (1 - w) * P->dens[i * 6 + 3 + k] + w * P->dens[(i + 1) * 6 + 3 + k];
        t1[k] = (1 - w) * P->frac[i * 6 + k] + w * P->frac[(i + 1) * 6 + k];
        t2[k] = (1 - w) * P->frac[i * 6 + 3 + k] + w * P->frac[(i + 1) * 6 + 3 + k];
        post[k] = pr[k] * pr[k] * lp;
        post[3 + k] = 2.0 * pr[pve_pair[k][0]] * pr[pve_pair[k][1]] * lm;
        s += post[k] + post[3 + k];
    }
    if (!(s > 0)) return 0;
    for (k = 0; k < 6; k++) post[k] /= s;
    return 1;
}

/* Expected fractions of CSF, GM and WM of a brain voxel. */
static void
pve_fractions(const Pve *P, double y, const double pr[3], double fr[3])
{
    double post[6], t1[3], t2[3];
    int k;

    if (!pve_posterior(P, y, pr, post, t1, t2))
    {
        /* beyond every tissue: the nearest pure one that has a prior, so
         * that a tissue the cleanup took out of the priors stays out */
        int kb = -1, pass;
        for (pass = 0; pass < 2 && kb < 0; pass++)
            for (k = 0; k < 3; k++)
                if ((pass || pr[k] > 0) && (kb < 0 || fabs(y - P->m[k]) < fabs(y - P->m[kb])))
                    kb = k;
        for (k = 0; k < 3; k++) fr[k] = (k == kb) ? 1.0 : 0.0;
        return;
    }
    for (k = 0; k < 3; k++) fr[k] = post[k];
    for (k = 0; k < 3; k++)
    {
        fr[pve_pair[k][0]] += post[3 + k] * (1.0 - t1[k]);
        fr[pve_pair[k][1]] += post[3 + k] * t1[k];
    }
}

/* The log weights of the Gaussians and of the classes, for posteriors. */
typedef struct
{
    double cg[MAXK], hg[MAXK], lwp[KB];
} GmmLogs;

static void
gmm_logs(const Gmm *G, GmmLogs *gl)
{
    int k, g;

    for (g = 0; g < G->K; g++)
    {
        gl->cg[g] = log(G->mg[g]) - 0.5 * log(2.0 * M_PI * G->vr[g]);
        gl->hg[g] = 0.5 / G->vr[g];
    }
    for (k = 0; k < KB; k++) gl->lwp[k] = log(G->wp[k]);
}

/* The posterior of the segmentation that lattice voxel i is brain; cr gets
 * its bias-corrected intensity. */
static double
lattice_brain(const State *S, const GmmLogs *gl, size_t i, double *cr)
{
    const Lattice *L = &S->L;
    const Gmm *G = &S->G;
    const float *lb = L->lb + KB * i;
    double lq[MAXK], mx = -HUGE_VAL, q[KB], sq = 0.0;
    int k, g;

    *cr = (double)L->f[i] * L->bf[i];
    for (g = 0; g < G->K; g++)
    {
        double d = *cr - G->mn[g];
        lq[g] = gl->cg[g] - gl->hg[g] * d * d + gl->lwp[G->lkp[g]] + lb[G->lkp[g]];
        if (lq[g] > mx) mx = lq[g];
    }
    for (k = 0; k < KB; k++) q[k] = 0.0;
    for (g = 0; g < G->K; g++) q[G->lkp[g]] += exp(lq[g] - mx);
    for (k = 0; k < KB; k++) sq += q[k];
    return (q[C_CSF] + q[C_GM] + q[C_WM]) / sq;
}

/* The variance that the image the parameters were estimated from has more
 * than the one that is segmented -- added noise, or what the filter took
 * out -- inside the brain and in bias-corrected intensities: the mean squared
 * difference of the two on the lattice, weighted by the posterior that a
 * voxel is brain.  Over the whole image the filter takes out about twice as
 * much (as SD: 18.7 against 10.5 in a 0.75 mm image, 10.1 against 5.7 and
 * 9.4 against 4.8 in two of 1 mm), because it also smooths the fine structure
 * of the scalp. */
static double
pve_excess_variance(const State *S)
{
    const Lattice *L = &S->L;
    double s = 0.0, w = 0.0;
    GmmLogs gl;
    size_t i;

    if (S->est == S->clean) return 0.0;
    gmm_logs(&S->G, &gl);
    for (i = 0; i < L->n; i++)
    {
        double cr, pb, d;
        if (!L->msk[i]) continue;
        pb = lattice_brain(S, &gl, i, &cr);
        d = ((double)L->f[i] - L->fseg[i]) * L->bf[i];
        s += pb * d * d;
        w += pb;
    }
    return (w > 0) ? s / w : 0.0;
}

/* One EM step for the means of the pure tissues on the lattice: with the
 * posterior of the six cases, a pure voxel pulls the mean of its tissue to
 * its intensity y and a mixture pulls (1-t) m_a + t m_b there, a
 * least-squares problem in the three means.  Voxels count by the posterior
 * of the segmentation that they are brain. */
static int
pve_em_step(const State *S, const Pve *P, double m[3])
{
    const Lattice *L = &S->L;
    double A[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0}, rhs[3] = {0, 0, 0};
    GmmLogs gl;
    size_t i;
    int k;

    gmm_logs(&S->G, &gl);
    for (i = 0; i < L->n; i++)
    {
        const float *lb = L->lb + KB * i;
        double cr, pb, pr[3], sb = 0.0;
        double post[6], t1[3], t2[3];
        if (!L->msk[i]) continue;
        pb = lattice_brain(S, &gl, i, &cr);
        if (!(pb > 0.01)) continue;
        for (k = 0; k < 3; k++)
        {
            pr[k] = exp(lb[pve_cls[k]] + gl.lwp[pve_cls[k]]);
            sb += pr[k];
        }
        for (k = 0; k < 3; k++) pr[k] /= sb;
        if (!pve_posterior(P, cr, pr, post, t1, t2)) continue;
        for (k = 0; k < 3; k++)
        {
            int a = pve_pair[k][0], b = pve_pair[k][1];
            double w = pb * post[3 + k];
            A[k * 4] += pb * post[k];
            rhs[k] += pb * post[k] * cr;
            A[a * 4] += w * (1.0 - 2.0 * t1[k] + t2[k]);
            A[b * 4] += w * t2[k];
            A[a * 3 + b] += w * (t1[k] - t2[k]);
            A[b * 3 + a] += w * (t1[k] - t2[k]);
            rhs[a] += w * (1.0 - t1[k]) * cr;
            rhs[b] += w * t1[k] * cr;
        }
    }
    if (solve_dense(A, rhs, 3)) return 1;
    for (k = 0; k < 3; k++) m[k] = rhs[k];
    return 0;
}

/* The intensities of pure CSF, GM and WM and the variance of a pure tissue.
 * The Gaussians of the segmentation are widened by partial volumes: the
 * narrower of those of grey and white matter is taken as the variance of a
 * pure tissue, and the means start at a Gaussian of each tissue and are
 * refined by EM steps of the partial volume model.  If these carry a mean
 * more than half way to another tissue the model does not fit, and the means
 * of the segmentation are kept.
 *
 * The Gaussian of pure CSF is left out of the variance: it is dark and
 * therefore narrow, and as the width of all three tissues it turned the
 * spread of white matter into partial volumes (label RMSE against CAT's p0
 * higher in 5 of 9 images, by up to 0.019, lower in none).  The image that is
 * segmented has less noise than the one the parameters come from (added
 * noise, or what the filter took out), so its tables get the variance less
 * pve_excess_variance(). */
static int
pve_setup(State *S, Pve *P)
{
    const Gmm *G = &S->G;
    double vmin = HUGE_VAL, dmin, m0[3], m[3];
    int it, k, g, ret, ok = 1;

    memset(P, 0, sizeof(Pve));
    /* grey and white matter: their largest Gaussian.  CSF: of its Gaussians
     * that hold at least a tenth of it and lie beyond grey matter as seen
     * from white matter, the one farthest from grey matter -- the others
     * hold its partial volumes with grey matter and are often the larger
     * ones.  Without such a Gaussian (another contrast), its largest. */
    for (k = 2; k >= 0; k--)
    {
        double best = -1.0;
        int gb = -1, pass;
        for (pass = (k == 0) ? 0 : 1; pass < 2 && gb < 0; pass++)
            for (g = 0; g < G->K; g++)
            {
                double score = G->mg[g];
                if (G->lkp[g] != pve_cls[k]) continue;
                if (pass == 0)
                {
                    double away = (G->mn[g] - P->m[1]) * (P->m[1] - P->m[2]);
                    if (!(G->mg[g] >= 0.1) || !(away > 0)) continue;
                    score = fabs(G->mn[g] - P->m[1]);
                }
                if (score > best)
                {
                    best = score;
                    gb = g;
                }
            }
        if (gb < 0) return 0;
        P->m[k] = G->mn[gb];
        if (k > 0 && G->vr[gb] < vmin) vmin = G->vr[gb];
        m0[k] = P->m[k];
    }
    dmin = fabs(m0[1] - m0[0]);
    if (fabs(m0[2] - m0[1]) < dmin) dmin = fabs(m0[2] - m0[1]);
    if (!(dmin > 0) || !isfinite(vmin)) return 0; /* tissues not separated */

    P->v = vmin;
    if (P->v < 1e-4 * dmin * dmin) P->v = 1e-4 * dmin * dmin;
    for (it = 0; it < 20 && ok; it++)
    {
        double change = 0.0;
        if ((ret = pve_tables(P))) return ret;
        if (pve_em_step(S, P, m)) ok = 0;
        for (k = 0; k < 3 && ok; k++)
        {
            if (!(fabs(m[k] - m0[k]) < 0.5 * dmin)) ok = 0;
            if (fabs(m[k] - P->m[k]) > change) change = fabs(m[k] - P->m[k]);
        }
        if (!ok) break;
        for (k = 0; k < 3; k++) P->m[k] = m[k];
        if (change < 1e-3 * dmin) break;
    }
    if (!ok)
        for (k = 0; k < 3; k++) P->m[k] = m0[k];
    P->v = vmin - pve_excess_variance(S);
    if (!(P->v > 1e-4 * dmin * dmin)) P->v = 1e-4 * dmin * dmin;
    if ((ret = pve_tables(P))) return ret;
    P->ok = 1;
    return 0;
}


/* ======================================================================== */
/* Partial volumes of the head                                              */
/* ======================================================================== */

/* Bone, soft tissue and air are mixtures of Gaussians themselves -- cortical
 * bone and marrow, muscle and fat, background noise -- so one intensity per
 * class, as for the brain, would read muscle as half fat and half air.  Here
 * a pure voxel of a head class follows the Gaussians of its class, and a
 * voxel that holds the fraction t of class b and 1-t of class a follows one
 * Gaussian of each: mean (1-t) m_g + t m_h, variance (1-t) v_g + t v_h, the
 * pairs weighted by their mixing proportions and t uniform.  The mixtures
 * are bone with soft tissue and soft tissue with air.  Tables over the
 * bias-corrected intensity hold the likelihood of every mixture and the
 * expected fraction t; their bins are a quarter of the narrowest Gaussian
 * wide, and the integral over t has nodes close enough that its mean moves
 * by a quarter of the narrower SD between them. */
typedef struct
{
    int ok, ntab;
    double y0, dy;
    float *dens; /* ntab x HEAD_NPAIR: likelihood of the mixture */
    float *frac; /* ntab x HEAD_NPAIR: E[t] */
} HeadPve;

static const int head_pair[HEAD_NPAIR][2] = {{C_BONE, C_SOFT}, {C_SOFT, C_AIR}};

static int
head_pve_setup(const State *S, HeadPve *H)
{
    const Gmm *G = &S->G;
    double lo = HUGE_VAL, hi = -HUGE_VAL, sdmin = HUGE_VAL, *den = NULL, *num = NULL;
    int g, h, p, i, it;

    memset(H, 0, sizeof(*H));
    for (g = 0; g < G->K; g++)
    {
        double sd = sqrt(G->vr[g]);
        if (G->lkp[g] < C_BONE) continue;
        if (!(sd > 0) || !isfinite(G->mn[g])) return 0;
        if (G->mn[g] - 6.0 * sd < lo) lo = G->mn[g] - 6.0 * sd;
        if (G->mn[g] + 6.0 * sd > hi) hi = G->mn[g] + 6.0 * sd;
        if (sd < sdmin) sdmin = sd;
    }
    if (!(hi > lo)) return 0;
    H->y0 = lo;
    H->dy = 0.25 * sdmin;
    H->ntab = (int)ceil((hi - lo) / H->dy) + 1;
    if (H->ntab > HEAD_NTAB_MAX)
    {
        H->ntab = HEAD_NTAB_MAX;
        H->dy = (hi - lo) / (H->ntab - 1);
    }
    H->dens = (float *)malloc(sizeof(float) * H->ntab * HEAD_NPAIR);
    H->frac = (float *)malloc(sizeof(float) * H->ntab * HEAD_NPAIR);
    den = (double *)malloc(sizeof(double) * H->ntab);
    num = (double *)malloc(sizeof(double) * H->ntab);
    if (!H->dens || !H->frac || !den || !num)
    {
        free(den);
        free(num);
        return 2;
    }
    for (p = 0; p < HEAD_NPAIR; p++)
    {
        for (i = 0; i < H->ntab; i++) den[i] = num[i] = 0.0;
        for (g = 0; g < G->K; g++)
            for (h = 0; h < G->K; h++)
            {
                double dm, smin, w;
                int nt;
                if (G->lkp[g] != head_pair[p][0] || G->lkp[h] != head_pair[p][1]) continue;
                dm = G->mn[h] - G->mn[g];
                smin = sqrt((G->vr[g] < G->vr[h]) ? G->vr[g] : G->vr[h]);
                nt = (int)ceil(4.0 * fabs(dm) / smin);
                if (nt < PVE_NINT) nt = PVE_NINT;
                if (nt > 8192) nt = 8192;
                w = G->mg[g] * G->mg[h] / nt;
                for (it = 0; it < nt; it++)
                {
                    double t = (it + 0.5) / nt, mu = G->mn[g] + t * dm;
                    double var = (1.0 - t) * G->vr[g] + t * G->vr[h], sd = sqrt(var);
                    double c = w / sqrt(2.0 * M_PI * var);
                    int i0 = (int)floor((mu - 6.0 * sd - H->y0) / H->dy);
                    int i1 = (int)ceil((mu + 6.0 * sd - H->y0) / H->dy);
                    if (i0 < 0) i0 = 0;
                    if (i1 > H->ntab - 1) i1 = H->ntab - 1;
                    for (i = i0; i <= i1; i++)
                    {
                        double d = H->y0 + i * H->dy - mu, e = c * exp(-0.5 * d * d / var);
                        den[i] += e;
                        num[i] += t * e;
                    }
                }
            }
        for (i = 0; i < H->ntab; i++)
        {
            H->dens[i * HEAD_NPAIR + p] = (float)den[i];
            H->frac[i * HEAD_NPAIR + p] = (float)((den[i] > 0) ? num[i] / den[i] : 0.5);
        }
    }
    free(den);
    free(num);
    H->ok = 1;
    return 0;
}

/* Expected fractions of bone, soft tissue and air of a head voxel with
 * intensity y, the likelihoods L of the six classes (their Gaussians) and
 * the priors pr of the three head classes (summing to one).  As for the
 * brain, a pure voxel has the prior pr_a^2 and a mixture 2 pr_a pr_b; the
 * mixture p is only a case where bit p of allowed is set.  Returns 0 if no
 * case has a likelihood. */
static int
head_fractions(const HeadPve *H, double y, const double L[KB], const double pr[3],
               unsigned allowed, double fr[3])
{
    double pos = (y - H->y0) / H->dy, w, post[3 + HEAD_NPAIR], t[HEAD_NPAIR], s = 0.0;
    int i, k, p;

    if (pos < 0) pos = 0;
    if (pos > H->ntab - 1) pos = H->ntab - 1;
    i = (int)floor(pos);
    if (i > H->ntab - 2) i = H->ntab - 2;
    w = pos - i;
    for (k = 0; k < 3; k++)
    {
        post[k] = pr[k] * pr[k] * L[C_BONE + k];
        s += post[k];
    }
    for (p = 0; p < HEAD_NPAIR; p++)
    {
        int a = head_pair[p][0] - C_BONE, b = head_pair[p][1] - C_BONE;
        double lm = (1 - w) * H->dens[i * HEAD_NPAIR + p] + w * H->dens[(i + 1) * HEAD_NPAIR + p];
        t[p] = (1 - w) * H->frac[i * HEAD_NPAIR + p] + w * H->frac[(i + 1) * HEAD_NPAIR + p];
        post[3 + p] = (allowed & (1u << p)) ? 2.0 * pr[a] * pr[b] * lm : 0.0;
        s += post[3 + p];
    }
    if (!(s > 0)) return 0;
    for (k = 0; k < 3; k++) fr[k] = post[k] / s;
    for (p = 0; p < HEAD_NPAIR; p++)
    {
        fr[head_pair[p][0] - C_BONE] += post[3 + p] / s * (1.0 - t[p]);
        fr[head_pair[p][1] - C_BONE] += post[3 + p] / s * t[p];
    }
    return 1;
}

/* ======================================================================== */
/* Cleanup of the brain classes (clean_gwc of spm_preproc_write8.m)         */
/* ======================================================================== */

/* The cleanup grows a brain mask from white matter: two erosions, then
 * conditional dilations through grey and white matter, then through CSF,
 * where every step thresholds the mask, multiplies it with the posteriors of
 * the classes it may grow into and smooths it with a 3x3x3 kernel.  Grey and
 * white matter outside the first mask and CSF outside the second are removed.
 * The steps run on the bounding box of the brain classes -- outside it the
 * products are zero, so the result is the same -- and stop once a step no
 * longer changes the thresholded product, after which every further step
 * would repeat it.  Within a step, a slice is only thresholded again if its
 * smoothed mask changed, and only smoothed again if the product changed in
 * it or a neighbouring slice; what is skipped would come out unchanged. */
typedef struct
{
    const unsigned char *P; /* KB volumes of the image, posteriors times 255 */
    size_t n;               /* voxels of the image */
    int d[3];               /* image dimensions */
    int o[3], bd[3];        /* bounding box: first voxel and dimensions */
    unsigned char *a;       /* bounding box: thresholded product */
    unsigned char *b;       /* bounding box: smoothed */
    int ncls;               /* GM + WM (2) or GM + WM + CSF (3) */
    double th;              /* threshold of b / 255 */
    unsigned char *chg;     /* per slice of the box: the product changed */
    unsigned char *chgb;    /* per slice of the box: the smoothed mask changed */
    double k[3];            /* the separable kernel */
} clean_ctx;

/* a = (b / 255 > th) (P_GM + P_WM [+ P_CSF]), saturated at 255 */
static void
clean_threshold_range(void *p, int z0, int z1)
{
    clean_ctx *c = (clean_ctx *)p;
    size_t nxy = (size_t)c->d[0] * c->d[1], bxy = (size_t)c->bd[0] * c->bd[1];
    int x, y, z, k;

    for (z = z0; z < z1; z++)
    {
        unsigned char chg = 0;
        if (!c->chgb[z])
        {
            c->chg[z] = 0;
            continue;
        }
        for (y = 0; y < c->bd[1]; y++)
            for (x = 0; x < c->bd[0]; x++)
            {
                size_t i = x + (size_t)c->bd[0] * y + bxy * z;
                size_t v = (x + c->o[0]) + (size_t)c->d[0] * (y + c->o[1]) + nxy * (z + c->o[2]);
                int val = 0;
                unsigned char nv;
                if ((double)c->b[i] / 255.0 > c->th)
                    for (k = 0; k < c->ncls; k++) val += c->P[k * c->n + v];
                nv = (unsigned char)((val > 255) ? 255 : val);
                chg |= (nv != c->a[i]);
                c->a[i] = nv;
            }
        c->chg[z] = chg;
    }
}

/* b = a smoothed as spm_conv_vol does it, in its order of summation: x and
 * y truncated at the edges of the image, z normalised by the weights used
 * there, rounded to uint8; inside the image but outside the box a is zero. */
static void
clean_smooth_range(void *p, int z0, int z1)
{
    clean_ctx *c = (clean_ctx *)p;
    int nx = c->bd[0], ny = c->bd[1], nz = c->bd[2], x, y, z, j;
    size_t bxy = (size_t)nx * ny;
    double *sl = (double *)malloc(sizeof(double) * (3 * bxy + (size_t)((nx > ny) ? nx : ny)));
    double *buf;
    int have[3] = {-1000, -1000, -1000};

    if (!sl) return;
    buf = sl + 3 * bxy;
    for (z = z0; z < z1; z++)
    {
        double out_w;
        int zz;
        unsigned char chgb = 0;
        if (!(c->chg[z] || (z > 0 && c->chg[z - 1]) || (z < nz - 1 && c->chg[z + 1])))
        {
            c->chgb[z] = 0;
            continue;
        }
        /* x and y smoothed slices z-1, z, z+1 of the box, in a ring */
        for (zz = z - 1; zz <= z + 1; zz++)
        {
            double *s;
            int slot = ((zz % 3) + 3) % 3;
            if (zz < 0 || zz >= nz || have[slot] == zz) continue;
            s = sl + (size_t)slot * bxy;
            for (y = 0; y < ny; y++)
            {
                const unsigned char *row = c->a + (size_t)nx * y + bxy * zz;
                for (x = 0; x < nx; x++)
                {
                    double v = 0.0;
                    if (x < nx - 1) v += row[x + 1] * c->k[0];
                    v += row[x] * c->k[1];
                    if (x > 0) v += row[x - 1] * c->k[2];
                    s[x + (size_t)nx * y] = v;
                }
            }
            for (x = 0; x < nx; x++)
            {
                for (y = 0; y < ny; y++) buf[y] = s[x + (size_t)nx * y];
                for (y = 0; y < ny; y++)
                {
                    double v = 0.0;
                    if (y < ny - 1) v += buf[y + 1] * c->k[0];
                    v += buf[y] * c->k[1];
                    if (y > 0) v += buf[y - 1] * c->k[2];
                    s[x + (size_t)nx * y] = v;
                }
            }
            have[slot] = zz;
        }
        /* z: beyond the box but inside the image the product is zero; at
         * the edge of the image spm_conv_vol divides by the weights used */
        out_w = 0.0;
        if (z < nz - 1 || c->o[2] + nz < c->d[2]) out_w += c->k[0];
        out_w += c->k[1];
        if (z > 0 || c->o[2] > 0) out_w += c->k[2];
        for (j = 0; j < (int)bxy; j++)
        {
            double v = 0.0;
            if (z < nz - 1) v += c->k[0] * sl[(size_t)((z + 1) % 3) * bxy + j];
            v += c->k[1] * sl[(size_t)(z % 3) * bxy + j];
            if (z > 0) v += c->k[2] * sl[(size_t)((z - 1) % 3) * bxy + j];
            v /= out_w;
            if (v < 0) v = 0.0;
            else if (v > 255.0) v = 255.0;
            {
                unsigned char nb = (unsigned char)floor(v + 0.5);
                chgb |= (nb != c->b[(size_t)j + bxy * z]);
                c->b[(size_t)j + bxy * z] = nb;
            }
        }
        c->chgb[z] = chgb;
    }
    free(sl);
}

/* Steps of clean_gwc from the current b; returns whether it stopped early. */
static void
clean_steps(clean_ctx *c, int nsteps, int nfirst, double th_first, double th, int nthreads)
{
    int j, z;

    for (j = 0; j < nsteps; j++)
    {
        int chg = 0;
        c->th = (j < nfirst) ? th_first : th;
        /* a new threshold or a new loop: every slice */
        if (j == 0 || j == nfirst)
            for (z = 0; z < c->bd[2]; z++) c->chgb[z] = 1;
        par_for(c->bd[2], nthreads, clean_threshold_range, c);
        par_for(c->bd[2], nthreads, clean_smooth_range, c);
        for (z = 0; z < c->bd[2]; z++) chg |= c->chg[z];
        /* the same product under the same threshold: a fixed point */
        if (!chg && j > nfirst) break;
    }
}

/* SPM's cleanup of the brain classes for posteriors P (KB volumes times
 * 255): keep gets bit 0 where grey and white matter stay and bit 1 where CSF
 * stays.  Bits 2 and 3 keep grey and white matter, and CSF, in the priors of
 * the partial volume model: wherever the masks reach, also where SPM's floor
 * of 5% takes them out of the posteriors (that would take the partial volumes
 * out of CSF next to cortex and ventricles), and nowhere else.  Beyond the
 * reach of the mask of grey and white matter a voxel cannot hold either of
 * them, whatever its posteriors: dura and vessels lie between CSF and grey
 * matter in intensity, the posteriors give them to a brighter Gaussian of
 * CSF, and the partial volume model, which has no such Gaussian, read them as
 * mixtures of CSF and grey matter -- a rim of grey matter around the brain,
 * separated from the cortex by CSF (on a 0.75 mm image 18 ml of the grey
 * matter fraction lay more than 1.1 mm from the cortex, 3 ml now).
 *
 * level 1 dilates through voxels whose mask exceeds 0.15, level 2 0.2.  SPM
 * runs 2 erosions and 32 steps of each dilation whatever the voxel size; here
 * they are 2 and 32 mm, which is the same at 1 mm.  At 0.5 mm SPM's counts
 * reach half as far and cost 1.1 ml more brain on a 0.5 mm image (CAT's
 * cat_main_clean_gwc.m scales them for the same reason).  With SPM's counts
 * the masks are those of clean_gwc in MATLAB, voxel for voxel. */
static int
clean_gwc(const State *S, const unsigned char *P, int level, unsigned char *keep)
{
    clean_ctx c;
    const int *d = S->dims;
    size_t n = (size_t)d[0] * d[1] * d[2], nxy = (size_t)d[0] * d[1], nb, i;
    int lo[3] = {d[0], d[1], d[2]}, hi[3] = {-1, -1, -1}, x, y, z, k, nero, nit;
    double sm, vx = cbrt(S->vx[0] * S->vx[1] * S->vx[2]), th1 = (level >= 2) ? 0.2 : 0.15;

    memset(keep, 4 | 8, n); /* nothing to remove outside the box */
    for (z = 0; z < d[2]; z++)
        for (y = 0; y < d[1]; y++)
            for (x = 0; x < d[0]; x++)
            {
                size_t v = x + (size_t)d[0] * y + nxy * z;
                if (!(P[C_GM * n + v] || P[C_WM * n + v] || P[C_CSF * n + v])) continue;
                if (x < lo[0]) lo[0] = x;
                if (x > hi[0]) hi[0] = x;
                if (y < lo[1]) lo[1] = y;
                if (y > hi[1]) hi[1] = y;
                if (z < lo[2]) lo[2] = z;
                if (z > hi[2]) hi[2] = z;
            }
    if (hi[0] < 0) return 0;

    memset(&c, 0, sizeof(c));
    c.P = P;
    c.n = n;
    for (k = 0; k < 3; k++)
    {
        c.d[k] = d[k];
        c.o[k] = lo[k];
        c.bd[k] = hi[k] - lo[k] + 1;
    }
    nb = (size_t)c.bd[0] * c.bd[1] * c.bd[2];
    c.a = (unsigned char *)calloc(nb, 1);
    c.b = (unsigned char *)malloc(nb);
    c.chg = (unsigned char *)calloc((size_t)c.bd[2], 1);
    c.chgb = (unsigned char *)calloc((size_t)c.bd[2], 1);
    if (!c.a || !c.b || !c.chg || !c.chgb)
    {
        free(c.a);
        free(c.b);
        free(c.chg);
        free(c.chgb);
        return 2;
    }
    /* [0.75 1 0.75] scaled so that the 3x3x3 kernel sums to one */
    sm = pow(15.625, 1.0 / 3.0);
    c.k[0] = c.k[2] = 0.75 / sm;
    c.k[1] = 1.0 / sm;
    nero = (int)floor(2.0 / vx + 0.5);
    nit = (int)floor(32.0 / vx + 0.5);
    if (nero < 1) nero = 1;
    if (nit < nero + 1) nit = nero + 1;

    /* grey and white matter, from white matter */
    for (z = 0; z < c.bd[2]; z++)
        for (y = 0; y < c.bd[1]; y++)
            for (x = 0; x < c.bd[0]; x++)
                c.b[x + (size_t)c.bd[0] * (y + (size_t)c.bd[1] * z)] =
                    P[C_WM * n + (x + lo[0]) + (size_t)d[0] * (y + lo[1]) + nxy * (z + lo[2])];
    c.ncls = 2;
    clean_steps(&c, nit, nero, 0.6, th1, S->nthreads);
    for (z = 0; z < c.bd[2]; z++)
        for (y = 0; y < c.bd[1]; y++)
            for (x = 0; x < c.bd[0]; x++)
            {
                size_t b = x + (size_t)c.bd[0] * (y + (size_t)c.bd[1] * z);
                size_t v = (x + lo[0]) + (size_t)d[0] * (y + lo[1]) + nxy * (z + lo[2]);
                double gw = (P[C_GM * n + v] + P[C_WM * n + v]) / 255.0;
                int reach = (double)c.b[b] / 255.0 > CLEAN_TH_FINAL;
                unsigned char kv = 0;
                if ((reach ? gw : 0.0) > CLEAN_TH_FINAL) kv |= 1;
                if (reach) kv |= 4;
                keep[v] = kv;
            }

    /* and CSF, from that mask */
    c.ncls = 3;
    for (i = 0; i < nb; i++) c.a[i] = 0;
    clean_steps(&c, nit, 0, th1, th1, S->nthreads);
    for (z = 0; z < c.bd[2]; z++)
        for (y = 0; y < c.bd[1]; y++)
            for (x = 0; x < c.bd[0]; x++)
            {
                size_t b = x + (size_t)c.bd[0] * (y + (size_t)c.bd[1] * z);
                size_t v = (x + lo[0]) + (size_t)d[0] * (y + lo[1]) + nxy * (z + lo[2]);
                double gwc = ((keep[v] & 1) ? (P[C_GM * n + v] + P[C_WM * n + v]) / 255.0 : 0.0) +
                             P[C_CSF * n + v] / 255.0;
                int reach = (double)c.b[b] / 255.0 > CLEAN_TH_FINAL;
                if ((reach ? gwc : 0.0) > CLEAN_TH_FINAL) keep[v] |= 2;
                if (reach) keep[v] |= 8;
            }
    free(c.a);
    free(c.b);
    free(c.chg);
    free(c.chgb);
    return 0;
}

/* ======================================================================== */
/* Final segmentation (cat_spm_preproc_write8.m)                            */
/* ======================================================================== */

typedef struct
{
    State *S;
    const float *iy;      /* 3 volumes: template voxel of every voxel */
    const Pve *P;
    const HeadPve *H;
    double *Bf[3];        /* bias basis at every voxel of an axis */
    int final;            /* 0: posteriors into Q or prob only; 1: everything */
    int head_pass;        /* after the final pass: the fractions of the head */
    unsigned char *hard;  /* final pass: the class of largest posterior */
    int from_prob;        /* final pass: the posteriors are those in prob (MRF) */
    const unsigned char *keep; /* cleanup: bit 0 keeps GM and WM, bit 1 CSF */
    unsigned char *prob;  /* KB volumes */
    float *Q;             /* KB volumes for the MRF, or NULL */
    float *label;
    unsigned char *pve;   /* KB volumes */
    float *corrected, *bias;
    double *vol;          /* per slice: KB sums of the posteriors */
    double *pvol;         /* per slice: KB sums of the fractions */
    double cg[MAXK], hg[MAXK], lwp[KB];
} write_ctx;

/* The posteriors of every voxel, and in the final pass the label map, the
 * fractions, the corrected image and the bias field.  Without the MRF and
 * the cleanup one pass does everything.  Otherwise a first pass writes the
 * posteriors (Q for the MRF, prob for the cleanup), and the final pass takes
 * them from prob (after the MRF) or computes them again, removes what the
 * cleanup removed and renormalises, with 1e-4 added to soft tissue as
 * cat_main_clean_gwc.m does, so that no voxel is left without a class.
 *
 * The fractions of the head need the classes of the neighbours, so a pass of
 * their own follows the final one: a voxel can only be a mixture of two head
 * classes where both are the class of largest posterior of some voxel of its
 * 3x3x3 neighbourhood.  Partial volumes arise at boundaries; without that
 * condition the priors, which are smooth and flat below the templates, let
 * mixtures of fat and air explain muscle, a pure tissue of intermediate
 * intensity (in the neck of a 0.5 mm image, 40% of the volume that the
 * posteriors gave to bone went to air). */
static void
write_range(void *p, int z0, int z1)
{
    write_ctx *c = (write_ctx *)p;
    State *S = c->S;
    const Lattice *L = &S->L;
    const Gmm *G = &S->G;
    const int *d = S->dims;
    size_t n = (size_t)d[0] * d[1] * d[2], nxy = (size_t)d[0] * d[1];
    int nb0 = L->nb[0], nb1 = L->nb[1], nb2 = L->nb[2], x, y, z, i, j, k, g, r;
    double *t1 = NULL, *tmp = NULL;

    if (L->nbt)
    {
        t1 = (double *)malloc(sizeof(double) * ((size_t)nb0 * nb1 + (size_t)d[0] * nb1));
        if (!t1) return;
        tmp = t1 + (size_t)nb0 * nb1;
    }
    for (z = z0; z < z1; z++)
    {
        double *vol = c->vol ? c->vol + (size_t)z * KB : NULL;
        double *pvol = c->pvol ? c->pvol + (size_t)z * KB : NULL;
        for (k = 0; k < KB; k++)
        {
            if (c->head_pass)
            {
                if (k >= C_BONE) pvol[k] = 0.0;
                continue;
            }
            if (vol) vol[k] = 0.0;
            if (pvol) pvol[k] = 0.0;
        }
        if (L->nbt)
        {
            const double *b2 = c->Bf[2] + (size_t)z * nb2;
            for (i = 0; i < nb0 * nb1; i++)
            {
                double s = 0.0;
                for (k = 0; k < nb2; k++) s += L->T[i + (size_t)nb0 * nb1 * k] * b2[k];
                t1[i] = s;
            }
            for (x = 0; x < d[0]; x++)
            {
                const double *b0 = c->Bf[0] + (size_t)x * nb0;
                for (j = 0; j < nb1; j++)
                {
                    double s = 0.0;
                    for (i = 0; i < nb0; i++) s += b0[i] * t1[i + nb0 * j];
                    tmp[(size_t)x * nb1 + j] = s;
                }
            }
        }
        for (y = 0; y < d[1]; y++)
            for (x = 0; x < d[0]; x++)
            {
                size_t v = x + (size_t)d[0] * y + nxy * z;
                double bf = 1.0, t[3], a[KB], pr[KB], q[KB], lk[KB], sq = 0.0, cr = 0.0, pb, ph;
                double fr[3], hf[3], lab;
                float fv = S->clean[v];
                int ok = is_valid(fv) && (!S->mask || S->mask[v]);

                if (L->nbt)
                {
                    const double *b1 = c->Bf[1] + (size_t)y * nb1;
                    double s = 0.0;
                    for (j = 0; j < nb1; j++) s += tmp[(size_t)x * nb1 + j] * b1[j];
                    bf = exp(s);
                }
                if (c->final && !c->head_pass && c->bias) c->bias[v] = (float)bf;
                if (c->final && !c->head_pass && c->corrected)
                    c->corrected[v] = isfinite(fv) ? (float)(fv * bf) : 0.0f;

                /* the prior */
                for (r = 0; r < 3; r++) t[r] = c->iy[r * n + v];
                if (!isfinite(t[0]) || !isfinite(t[1]) || !isfinite(t[2]))
                    for (r = 0; r < 3; r++)
                        t[r] = S->Tvox[r * 4] * x + S->Tvox[r * 4 + 1] * y + S->Tvox[r * 4 + 2] * z +
                               S->Tvox[r * 4 + 3];
                tpl_sample(&S->tpl, t[0], t[1], t[2], a, NULL);
                softmax_k(a, c->lwp, pr, NULL);

                /* the likelihoods of the classes and the posterior; without
                 * data the posterior is the prior */
                for (k = 0; k < KB; k++) lk[k] = 0.0;
                if (ok)
                {
                    cr = (double)fv * bf;
                    for (g = 0; g < G->K; g++)
                    {
                        double dd = cr - G->mn[g];
                        lk[G->lkp[g]] += exp(c->cg[g] - c->hg[g] * dd * dd);
                    }
                }
                if (c->final && c->from_prob)
                    for (k = 0; k < KB; k++) q[k] = c->prob[k * n + v] / 255.0;
                else if (ok)
                {
                    for (k = 0; k < KB; k++)
                    {
                        q[k] = (lk[k] + DBL_EPS) * pr[k];
                        sq += q[k];
                    }
                    for (k = 0; k < KB; k++) q[k] /= sq;
                }
                else
                    for (k = 0; k < KB; k++) q[k] = pr[k];
                if (!c->final)
                {
                    for (k = 0; k < KB; k++)
                    {
                        if (c->Q) c->Q[k * n + v] = (float)q[k];
                        if (c->prob) c->prob[k * n + v] = (unsigned char)floor(255.0 * q[k] + 0.5);
                    }
                    continue;
                }
                if (c->keep)
                {
                    /* the posteriors lose what the cleanup removed; the priors
                     * of the partial volume model only what its masks did not
                     * reach, not what fell below its 5% */
                    if (!(c->keep[v] & 1)) q[C_GM] = q[C_WM] = 0.0;
                    if (!(c->keep[v] & 2)) q[C_CSF] = 0.0;
                    if (!(c->keep[v] & 4)) pr[C_GM] = pr[C_WM] = 0.0;
                    if (!(c->keep[v] & 8)) pr[C_CSF] = 0.0;
                    q[C_SOFT] += 1e-4;
                    sq = 0.0;
                    for (k = 0; k < KB; k++) sq += q[k];
                    for (k = 0; k < KB; k++) q[k] /= sq;
                }
                if (c->head_pass)
                {
                    double ht[3], sh = pr[C_BONE] + pr[C_SOFT] + pr[C_AIR];
                    unsigned allowed = 0;
                    int seen[KB] = {0, 0, 0, 0, 0, 0}, done = 0, dx, dy, dz;
                    ph = q[C_BONE] + q[C_SOFT] + q[C_AIR];
                    for (dz = -1; dz <= 1; dz++)
                        for (dy = -1; dy <= 1; dy++)
                            for (dx = -1; dx <= 1; dx++)
                            {
                                int xx = x + dx, yy = y + dy, zz = z + dz;
                                if (xx < 0 || yy < 0 || zz < 0 || xx >= d[0] || yy >= d[1] || zz >= d[2])
                                    continue;
                                seen[c->hard[xx + (size_t)d[0] * yy + nxy * zz]] = 1;
                            }
                    for (k = 0; k < HEAD_NPAIR; k++)
                        if (seen[head_pair[k][0]] && seen[head_pair[k][1]]) allowed |= 1u << k;
                    if (ok && ph > 0 && sh > 0)
                    {
                        for (k = 0; k < 3; k++) ht[k] = pr[C_BONE + k] / sh;
                        done = head_fractions(c->H, cr, lk, ht, allowed, hf);
                    }
                    if (!done)
                        for (k = 0; k < 3; k++) hf[k] = (ph > 0) ? q[C_BONE + k] / ph : 0.0;
                    for (k = 0; k < 3; k++)
                    {
                        pvol[C_BONE + k] += ph * hf[k];
                        if (c->pve)
                            c->pve[(C_BONE + k) * n + v] = (unsigned char)floor(255.0 * ph * hf[k] + 0.5);
                    }
                    continue;
                }
                for (k = 0; k < KB; k++)
                {
                    if (c->prob) c->prob[k * n + v] = (unsigned char)floor(255.0 * q[k] + 0.5);
                    vol[k] += q[k];
                }
                if (c->hard)
                {
                    int kb = 0;
                    for (k = 1; k < KB; k++)
                        if (q[k] > q[kb]) kb = k;
                    c->hard[v] = (unsigned char)kb;
                }

                /* fractions of CSF, GM and WM given that the voxel is brain,
                 * and of bone, soft tissue and air given that it is head */
                if (!c->label && !c->pve) continue;
                pb = q[C_CSF] + q[C_GM] + q[C_WM];
                ph = q[C_BONE] + q[C_SOFT] + q[C_AIR];
                if (c->P && c->P->ok && ok && pb > 0)
                {
                    double pt[3], sb = pr[C_CSF] + pr[C_GM] + pr[C_WM];
                    for (k = 0; k < 3; k++) pt[k] = (sb > 0) ? pr[pve_cls[k]] / sb : 0.0;
                    pve_fractions(c->P, cr, pt, fr);
                }
                else
                    for (k = 0; k < 3; k++) fr[k] = (pb > 0) ? q[pve_cls[k]] / pb : 0.0;
                /* the head: its posteriors, unless the head pass follows */
                for (k = 0; k < 3; k++) hf[k] = (ph > 0) ? q[C_BONE + k] / ph : 0.0;
                lab = fr[0] + 2.0 * fr[1] + 3.0 * fr[2];
                if (c->label) c->label[v] = (float)(pb * lab);
                pvol[C_CSF] += pb * fr[0];
                pvol[C_GM] += pb * fr[1];
                pvol[C_WM] += pb * fr[2];
                if (!c->hard)
                    for (k = 0; k < 3; k++) pvol[C_BONE + k] += ph * hf[k];
                if (c->pve)
                {
                    c->pve[C_GM * n + v] = (unsigned char)floor(255.0 * pb * fr[1] + 0.5);
                    c->pve[C_WM * n + v] = (unsigned char)floor(255.0 * pb * fr[2] + 0.5);
                    c->pve[C_CSF * n + v] = (unsigned char)floor(255.0 * pb * fr[0] + 0.5);
                    for (k = 0; k < 3 && !c->hard; k++)
                        c->pve[(C_BONE + k) * n + v] = (unsigned char)floor(255.0 * ph * hf[k] + 0.5);
                }
            }
    }
    free(t1);
}

typedef struct
{
    const int *d;
    const float *Q;
    unsigned char *P;
    float G, w[3];
    int parity;
} mrf_ctx;

/* One colour of a red-black sweep of spm_mrf.c with a diagonal connectivity
 * matrix: q = p exp(G a) / sum, a the mean of the six neighbours' q. */
static void
mrf_range(void *p, int z0, int z1)
{
    mrf_ctx *c = (mrf_ctx *)p;
    const int *d = c->d;
    size_t nxy = (size_t)d[0] * d[1], n = nxy * d[2];
    int x, y, z, k;

    for (z = z0; z < z1; z++)
        for (y = 0; y < d[1]; y++)
            for (x = (y + z + c->parity) & 1; x < d[0]; x += 2)
            {
                size_t v = x + (size_t)d[0] * y + nxy * z;
                float a[KB], e[KB], mx = -3.4028e+38f, se = 0.0f;
                for (k = 0; k < KB; k++)
                {
                    const unsigned char *q = c->P + k * n + v;
                    a[k] = 0.0f;
                    if (z > 0) a[k] += q[-(ptrdiff_t)nxy] * c->w[2];
                    if (z < d[2] - 1) a[k] += q[nxy] * c->w[2];
                    if (y > 0) a[k] += q[-(ptrdiff_t)d[0]] * c->w[1];
                    if (y < d[1] - 1) a[k] += q[d[0]] * c->w[1];
                    if (x > 0) a[k] += q[-1] * c->w[0];
                    if (x < d[0] - 1) a[k] += q[1] * c->w[0];
                    e[k] = c->G * a[k] / (255.0f * 6.0f);
                    if (e[k] > mx) mx = e[k];
                }
                for (k = 0; k < KB; k++)
                {
                    e[k] = (float)exp((double)(e[k] - mx)) * c->Q[k * n + v];
                    se += e[k];
                }
                se = 255.0f / se;
                for (k = 0; k < KB; k++) c->P[k * n + v] = (unsigned char)(e[k] * se + 0.5f);
            }
}

/* The MRF cleanup of cat_spm_preproc_write8.m: ten sweeps from P = 0.  The
 * two colours of a sweep are independent of each other within themselves,
 * so every colour runs in parallel over the slices. */
static void
mrf_cleanup(const State *S, const float *Q, unsigned char *P)
{
    mrf_ctx c;
    size_t n = (size_t)S->dims[0] * S->dims[1] * S->dims[2];
    int it, k;

    c.d = S->dims;
    c.Q = Q;
    c.P = P;
    c.G = (float)S->o.mrf;
    for (k = 0; k < 3; k++) c.w[k] = (float)(1.0 / (S->vx[k] * S->vx[k]));
    memset(P, 0, KB * n);
    for (it = 0; it < 10; it++)
        for (c.parity = 0; c.parity < 2; c.parity++)
            par_for(S->dims[2], S->nthreads, mrf_range, &c);
}

/* The posteriors of every voxel, the label map, the fractions and the
 * corrected image: the MRF of cat_spm_preproc_write8.m, then its cleanup. */
static int
write_segmentation(State *S, const float *iy, const Pve *P, const HeadPve *H, unsigned char *prob,
                   float *label, unsigned char *pve, float *corrected, float *bias,
                   CAT_UnifiedSegResult *res)
{
    write_ctx c;
    const Lattice *L = &S->L;
    size_t n = (size_t)S->dims[0] * S->dims[1] * S->dims[2];
    unsigned char *own = NULL, *keep = NULL;
    double voxvol;
    int g, k, z, ret = 0;

    memset(&c, 0, sizeof(c));
    c.S = S;
    c.iy = iy;
    c.P = P;
    c.H = H;
    c.vol = (double *)calloc((size_t)S->dims[2] * KB, sizeof(double));
    c.pvol = (double *)calloc((size_t)S->dims[2] * KB, sizeof(double));
    if (!c.vol || !c.pvol)
    {
        ret = 2;
        goto done;
    }
    for (k = 0; k < 3 && L->nbt; k++)
    {
        c.Bf[k] = (double *)malloc(sizeof(double) * (size_t)S->dims[k] * L->nb[k]);
        if (!c.Bf[k])
        {
            ret = 2;
            goto done;
        }
        bias_basis(L->model, S->dims[k], L->nb[k], S->dims[k], 1, c.Bf[k]);
    }
    for (g = 0; g < S->G.K; g++)
    {
        c.cg[g] = log(S->G.mg[g]) - 0.5 * log(2.0 * M_PI * S->G.vr[g]);
        c.hg[g] = 0.5 / S->G.vr[g];
    }
    for (k = 0; k < KB; k++) c.lwp[k] = log(S->G.wp[k]);
    tpl_load(&S->tpl, NLEV - 1, S->nthreads);

    if (S->o.mrf > 0 || S->o.cleanup > 0)
    {
        if (!prob) prob = own = (unsigned char *)malloc(KB * n);
        if (S->o.mrf > 0) c.Q = (float *)malloc(sizeof(float) * KB * n);
        if (!prob || (S->o.mrf > 0 && !c.Q))
        {
            ret = 2;
            goto done;
        }
        /* the MRF starts from zero; the cleanup needs the posteriors */
        c.prob = (S->o.mrf > 0) ? NULL : prob;
        par_for(S->dims[2], S->nthreads, write_range, &c);
        if (c.Q)
        {
            mrf_cleanup(S, c.Q, prob);
            free(c.Q);
            c.Q = NULL;
            c.from_prob = 1;
        }
        if (S->o.cleanup > 0)
        {
            double t0 = wall_time();
            keep = (unsigned char *)malloc(n);
            if (!keep || (ret = clean_gwc(S, prob, S->o.cleanup, keep)))
            {
                if (!ret) ret = 2;
                goto done;
            }
            if (S->o.verbose) printf("Cleanup: %.1f s\n", wall_time() - t0);
        }
    }
    c.final = 1;
    c.keep = keep;
    c.prob = prob;
    c.label = label;
    c.pve = pve;
    c.corrected = corrected;
    c.bias = bias;
    if (H && H->ok && (label || pve))
    {
        c.hard = (unsigned char *)malloc(n);
        if (!c.hard)
        {
            ret = 2;
            goto done;
        }
    }
    par_for(S->dims[2], S->nthreads, write_range, &c);
    if (c.hard)
    {
        c.head_pass = 1;
        par_for(S->dims[2], S->nthreads, write_range, &c);
    }

    if (res)
    {
        voxvol = fabs(S->mat[0] * (S->mat[5] * S->mat[10] - S->mat[6] * S->mat[9]) -
                      S->mat[1] * (S->mat[4] * S->mat[10] - S->mat[6] * S->mat[8]) +
                      S->mat[2] * (S->mat[4] * S->mat[9] - S->mat[5] * S->mat[8]));
        for (k = 0; k < KB; k++)
        {
            double s = 0.0, sp = 0.0;
            for (z = 0; z < S->dims[2]; z++)
            {
                s += c.vol[(size_t)z * KB + k];
                sp += c.pvol[(size_t)z * KB + k];
            }
            res->volume[k] = s * voxvol / 1000.0;
            res->pve_volume[k] = (label || pve) ? sp * voxvol / 1000.0 : 0.0;
        }
    }

done:
    free(c.vol);
    free(c.pvol);
    free(c.Q);
    free(c.hard);
    free(own);
    free(keep);
    for (k = 0; k < 3; k++) free(c.Bf[k]);
    return ret;
}

/* ======================================================================== */
/* Affine registration from several starting points                         */
/* ======================================================================== */

#define AFF_SLAB 130.0   /* mm of the head below its top that anchor a start */
#define AFF_NSTART 16    /* 3 zooms x 5 starts, and a given matrix */
#define AFF_NREFINE 3    /* candidates that are refined */

/* Centroid in mm of the upper AFF_SLAB * zoom mm of the points with weight
 * w > 0: the top is the 99.5% quantile of their height, so that a few bright
 * voxels above the head do not count.  The amount of neck and shoulders in
 * the field of view, which moves the centre of mass of the whole head by
 * centimetres, does not enter. */
static int
slab_centroid(const float *xyz, const float *w, ptrdiff_t n, double slab, double c[3])
{
    float *zs = (float *)malloc(sizeof(float) * (size_t)(n > 0 ? n : 1));
    double top, sw = 0.0;
    ptrdiff_t i;

    c[0] = c[1] = c[2] = 0.0;
    if (!zs || n < 10)
    {
        free(zs);
        return 1;
    }
    for (i = 0; i < n; i++) zs[i] = xyz[3 * i + 2];
    top = select_float(zs, n, (ptrdiff_t)(0.995 * (double)(n - 1)));
    free(zs);
    for (i = 0; i < n; i++)
        if (xyz[3 * i + 2] > top - slab)
        {
            c[0] += (double)w[i] * xyz[3 * i];
            c[1] += (double)w[i] * xyz[3 * i + 1];
            c[2] += (double)w[i] * xyz[3 * i + 2];
            sw += w[i];
        }
    if (!(sw > 0)) return 1;
    c[0] /= sw;
    c[1] /= sw;
    c[2] /= sw;
    return 0;
}

/* The anchors of the starting estimates: the centroid of the top of the head
 * in the image (voxels brighter than the mean) and in the template
 * (everything but air), which is larger by about the zoom of MNI space. */
static int
head_anchors(const State *S, const AffLat *A, double cs[3], double ct[3])
{
    const float *air = S->tpl.tpl[0] + (size_t)C_AIR * S->tpl.n;
    size_t nxy = (size_t)S->dims[0] * S->dims[1], cnt = 0;
    double mean = 0.0;
    float *xyz, *w;
    ptrdiff_t m = 0;
    int x, y, z, r, ret;

    for (z = 0; z < A->d[2]; z++)
        for (y = 0; y < A->d[1]; y++)
            for (x = 0; x < A->d[0]; x++)
                if (A->msk[x + (size_t)A->d[0] * (y + (size_t)A->d[1] * z)])
                {
                    mean += S->orig[(size_t)x * A->sk[0] + (size_t)S->dims[0] * ((size_t)y * A->sk[1]) +
                                     nxy * ((size_t)z * A->sk[2])];
                    cnt++;
                }
    mean /= (double)cnt;
    xyz = (float *)malloc(sizeof(float) * 3 * ((A->n > S->tpl.n) ? A->n : S->tpl.n));
    w = (float *)malloc(sizeof(float) * ((A->n > S->tpl.n) ? A->n : S->tpl.n));
    if (!xyz || !w)
    {
        free(xyz);
        free(w);
        return 2;
    }
    for (z = 0; z < A->d[2]; z++)
        for (y = 0; y < A->d[1]; y++)
            for (x = 0; x < A->d[0]; x++)
            {
                double xv[3] = {(double)x * A->sk[0], (double)y * A->sk[1], (double)z * A->sk[2]};
                if (!A->msk[x + (size_t)A->d[0] * (y + (size_t)A->d[1] * z)]) continue;
                if (!(S->orig[(size_t)xv[0] + (size_t)S->dims[0] * (size_t)xv[1] + nxy * (size_t)xv[2]] > mean))
                    continue;
                for (r = 0; r < 3; r++)
                    xyz[3 * m + r] = (float)(S->mat[r * 4] * xv[0] + S->mat[r * 4 + 1] * xv[1] +
                                             S->mat[r * 4 + 2] * xv[2] + S->mat[r * 4 + 3]);
                w[m++] = 1.0f;
            }
    ret = slab_centroid(xyz, w, m, AFF_SLAB, cs);
    m = 0;
    for (z = 0; z < S->tdims[2] && !ret; z += 2)
        for (y = 0; y < S->tdims[1]; y += 2)
            for (x = 0; x < S->tdims[0]; x += 2)
            {
                float wa = 1.0f - air[x + (size_t)S->tdims[0] * (y + (size_t)S->tdims[1] * z)];
                if (!(wa > 0.5f)) continue;
                for (r = 0; r < 3; r++)
                    xyz[3 * m + r] = (float)(S->tmat[r * 4] * x + S->tmat[r * 4 + 1] * y +
                                             S->tmat[r * 4 + 2] * z + S->tmat[r * 4 + 3]);
                w[m++] = wa;
            }
    if (!ret) ret = slab_centroid(xyz, w, m, AFF_SLAB * 1.12, ct);
    free(xyz);
    free(w);
    return ret ? 3 : 0;
}

/* CAT26's start (cat_vol_set_com.m): the centre of mass of the voxels
 * brighter than the median of those at least as bright as the median of the
 * image, which cat_run_job.m puts onto [0 -20 -15] mm, the centre of mass of
 * MNI space, before its first affine registration.  Taken from the voxels of
 * the coarse lattice, zeros included as there. */
static const double com_mni[3] = {0.0, -20.0, -15.0};

static int
com_anchor(const State *S, const AffLat *A, double cs[3])
{
    size_t nxy = (size_t)S->dims[0] * S->dims[1];
    float *val = (float *)malloc(sizeof(float) * (A->n ? A->n : 1));
    double sw = 0.0, m1, m2;
    ptrdiff_t m = 0, i, j;
    int x, y, z, r;

    cs[0] = cs[1] = cs[2] = 0.0;
    if (!val) return 2;
    for (z = 0; z < A->d[2]; z++)
        for (y = 0; y < A->d[1]; y++)
            for (x = 0; x < A->d[0]; x++)
            {
                size_t v = (size_t)x * A->sk[0] + (size_t)S->dims[0] * ((size_t)y * A->sk[1]) +
                           nxy * ((size_t)z * A->sk[2]);
                if (isfinite(S->orig[v]) && (!S->mask || S->mask[v])) val[m++] = S->orig[v];
            }
    if (m < 100)
    {
        free(val);
        return 3;
    }
    m1 = select_float(val, m, m / 2);
    for (i = 0, j = 0; i < m; i++)
        if (val[i] >= m1) val[j++] = val[i];
    m2 = select_float(val, j, j / 2);
    free(val);
    for (z = 0; z < A->d[2]; z++)
        for (y = 0; y < A->d[1]; y++)
            for (x = 0; x < A->d[0]; x++)
            {
                double xv[3] = {(double)x * A->sk[0], (double)y * A->sk[1], (double)z * A->sk[2]};
                size_t v = (size_t)xv[0] + (size_t)S->dims[0] * (size_t)xv[1] + nxy * (size_t)xv[2];
                if (!(isfinite(S->orig[v]) && S->orig[v] > m2) || (S->mask && !S->mask[v])) continue;
                for (r = 0; r < 3; r++)
                    cs[r] += S->mat[r * 4] * xv[0] + S->mat[r * 4 + 1] * xv[1] + S->mat[r * 4 + 2] * xv[2] +
                             S->mat[r * 4 + 3];
                sw += 1.0;
            }
    if (!(sw > 0)) return 3;
    for (r = 0; r < 3; r++) cs[r] /= sw;
    return 0;
}

/* M maps cs to ct with a zoom and a rotation by pitch radians about x. */
static void
start_matrix(const double cs[3], const double ct[3], double zoom, double pitch, double M[16])
{
    double R[9] = {1, 0, 0, 0, cos(pitch), -sin(pitch), 0, sin(pitch), cos(pitch)};
    int r, c;

    m4_eye(M);
    for (r = 0; r < 3; r++)
    {
        double t = ct[r];
        for (c = 0; c < 3; c++)
        {
            M[r * 4 + c] = zoom * R[r * 3 + c];
            t -= zoom * R[r * 3 + c] * cs[c];
        }
        M[r * 4 + 3] = t;
    }
}

/* Largest distance in mm between the images of the points 60 mm around p
 * under two matrices: are these two solutions the same? */
static double
affine_distance(const double A[16], const double B[16], const double p[3])
{
    double dmax = 0.0;
    int c, r;

    for (c = 0; c < 8; c++)
    {
        double q[3] = {p[0] + ((c & 1) ? 60.0 : -60.0), p[1] + ((c & 2) ? 60.0 : -60.0),
                       p[2] + ((c & 4) ? 60.0 : -60.0)};
        double d2 = 0.0;
        for (r = 0; r < 3; r++)
        {
            double d = (A[r * 4] - B[r * 4]) * q[0] + (A[r * 4 + 1] - B[r * 4 + 1]) * q[1] +
                       (A[r * 4 + 2] - B[r * 4 + 2]) * q[2] + A[r * 4 + 3] - B[r * 4 + 3];
            d2 += d * d;
        }
        if (d2 > dmax) dmax = d2;
    }
    return sqrt(dmax);
}

/* The affine registration of cat_spm_preproc_run.m, with a wider search.
 * SPM runs a coarse, almost rigid registration from the origin of the header
 * and from the centre of the field of view and refines the better one.  Both
 * starts fail when the head is far from where the header says and tilted, a
 * head can be caught at too small a zoom (MNI space is larger than most
 * heads), and the coarse objective does not always rank the solutions as the
 * fine one does.  So the coarse registration starts from the header, from
 * the top of the head put onto that of the template (straight and tilted by
 * 20 degrees either way) and from CAT26's centre of mass, each with three
 * zooms; the best AFF_NREFINE different solutions are refined, and the best
 * of those is kept.  The centre of mass changed the result in one of 54
 * shifted, tilted and cropped images (6 subjects), from 5.7 to 2.7 mm off. */
static int
affine_register(State *S, const double *init, double Aff[16], double *ll_out)
{
    static const double zooms[3] = {1.0, 1.2, 1.4};
    AffLat A;
    double M[AFF_NSTART][16], ll[AFF_NSTART], cand[AFF_NREFINE][16], cs[3], ct[3], cm[3];
    double bestll = -HUGE_VAL;
    int ns = 0, nc = 0, i, k, zi, err = 0, ret, used[AFF_NSTART], use_com;

    tpl_load(&S->tpl, 0, S->nthreads);
    if ((ret = afflat_build(&A, S->orig, S->mask, S->dims, S->vx, AFF_SAMP_COARSE))) goto done;
    A.ff = fudge_factor(16.0, S->vx, A.sk);
    if ((ret = head_anchors(S, &A, cs, ct))) goto done;
    use_com = !com_anchor(S, &A, cm);
    if (init) memcpy(M[ns++], init, sizeof(double) * 16);
    if (!init || S->o.affine_mode == CAT_USEG_AFFINE_ESTIMATE)
        for (zi = 0; zi < 3; zi++)
        {
            /* the header: a zoom about the anchor of the template */
            start_matrix(ct, ct, zooms[zi], 0.0, M[ns++]);
            start_matrix(cs, ct, zooms[zi], 0.0, M[ns++]);
            start_matrix(cs, ct, zooms[zi], 20.0 * M_PI / 180.0, M[ns++]);
            start_matrix(cs, ct, zooms[zi], -20.0 * M_PI / 180.0, M[ns++]);
            if (use_com) start_matrix(cm, com_mni, zooms[zi], 0.0, M[ns++]);
        }
    for (i = 0; i < ns; i++)
    {
        ll[i] = maff8(&A, &S->tpl, S->mat, S->tmat, M[i], S->o.affreg, 80, S->nthreads, &err);
        if (err) ll[i] = -HUGE_VAL;
        used[i] = 0;
        if (S->o.verbose > 1) printf("  affine start %2d: %.4f bits\n", i + 1, ll[i]);
    }
    afflat_free(&A);

    /* the best solutions that differ by more than 5 mm around the head */
    for (k = 0; k < AFF_NREFINE; k++)
    {
        int b = -1, same;
        do
        {
            b = -1;
            for (i = 0; i < ns; i++)
                if (!used[i] && isfinite(ll[i]) && (b < 0 || ll[i] > ll[b])) b = i;
            if (b < 0) break;
            used[b] = 1;
            same = 0;
            for (i = 0; i < nc; i++)
                if (affine_distance(M[b], cand[i], cs) < 5.0) same = 1;
        } while (same);
        if (b < 0) break;
        memcpy(cand[nc++], M[b], sizeof(double) * 16);
    }
    if (nc == 0)
    {
        ret = 3;
        goto done;
    }
    /* closer to rigid first, then with the regularisation of the image */
    if ((ret = afflat_build(&A, S->orig, S->mask, S->dims, S->vx, AFF_SAMP_FINE))) goto done;
    for (k = 0; k < nc; k++)
    {
        double l;
        A.ff = fudge_factor(16.0, S->vx, A.sk);
        l = maff8(&A, &S->tpl, S->mat, S->tmat, cand[k], S->o.affreg, 80, S->nthreads, &err);
        A.ff = fudge_factor(0.0, S->vx, A.sk);
        if (!err) l = maff8(&A, &S->tpl, S->mat, S->tmat, cand[k], S->o.affreg, 80, S->nthreads, &err);
        if (err == 2)
        {
            ret = 2;
            goto done;
        }
        if (S->o.verbose) printf("  affine candidate %d: %.4f bits%s\n", k + 1, l, err ? " (failed)" : "");
        if (!err && l > bestll)
        {
            bestll = l;
            memcpy(Aff, cand[k], sizeof(double) * 16);
        }
    }
    if (!isfinite(bestll))
    {
        ret = 3;
        goto done;
    }
    *ll_out = bestll;

done:
    afflat_free(&A);
    return ret;
}

/* ======================================================================== */
/* Public interface                                                         */
/* ======================================================================== */

/**
 * \brief Fill a CAT_UnifiedSegOpts with the defaults.
 *
 * \param opts (out) options to initialise
 */
void
CAT_UnifiedSegOptionsInit(CAT_UnifiedSegOpts *opts)
{
    static const int ngaus[KB] = {1, 1, 2, 3, 4, 2};
    int k;

    for (k = 0; k < KB; k++) opts->ngaus[k] = ngaus[k];
    opts->samp = 3.0;
    opts->wp_reg = 0.5;
    opts->bias_model = CAT_USEG_BIAS_SPLINE;
    opts->biasfwhm = 60.0;
    opts->biasreg = 1e-3;
    opts->tol = 1e-4;
    opts->affreg = CAT_USEG_AFFREG_MNI;
    opts->affine_mode = CAT_USEG_AFFINE_ESTIMATE;
    opts->noise = 0.0;
    opts->noise_seed = 1;
    opts->sanlm = 0;
    opts->sanlm_strength = 1.0;
    opts->sanlm_rician = 0;
    opts->mrf = 0.0;
    opts->pve = 1;
    opts->pve_head = 0;
    opts->cleanup = 1;
    CAT_VolShootOptionsInit(&opts->shoot);
    opts->shoot.regstr = 0.5;
    opts->n_threads = 0;
    opts->verbose = 0;
}

/**
 * \brief Release what CAT_UnifiedSeg() allocated inside a result.
 *
 * \param res (in/out) result; NULL is ignored
 */
void
CAT_UnifiedSegResultFree(CAT_UnifiedSegResult *res)
{
    if (!res) return;
    free(res->bias_coef);
    res->bias_coef = NULL;
}

/* The image that is segmented and the one the parameters are estimated from.
 *
 * With the SANLM filter the two differ: the segmentation is computed from
 * the filtered image, the parameters are estimated from the image as it
 * came.  Estimated from the filtered image, the Gaussian of white matter
 * gets too narrow against that of grey matter, which partial volumes keep
 * wide, and the result is worse than without the filter (mean Dice of nine
 * images 0.860 against 0.866, and 0.870 this way).
 *
 * opts->noise adds Gaussian noise to the estimation image until its noise
 * reaches that level; an image that is noisier already is left alone. */
static int
prepare_images(State *S, const float *src, float **clean_own, float **est_own,
               CAT_UnifiedSegResult *res)
{
    size_t n = (size_t)S->dims[0] * S->dims[1] * S->dims[2], i;
    double t0 = wall_time(), sd, target, add = 0.0, range;
    int ret;

    S->orig = S->clean = S->est = src;
    if ((ret = robust_range(src, n, S->mask, S->range))) return ret;
    range = S->range[1] - S->range[0];
    sd = estimate_noise(src, S->dims, S->mask, S->range[0] + 0.1 * range);
    if (S->o.sanlm)
    {
        *clean_own = (float *)malloc(sizeof(float) * n);
        if (!*clean_own) return 2;
        memcpy(*clean_own, src, sizeof(float) * n);
        sanlm(*clean_own, 3, 1, S->o.sanlm_rician, S->o.sanlm_strength, S->dims);
        /* the filter must not create or remove voxels without data */
        for (i = 0; i < n; i++)
            if (!is_valid(src[i]) || !is_valid((*clean_own)[i])) (*clean_own)[i] = src[i];
        S->clean = *clean_own;
        if (res) res->time_filter = wall_time() - t0;
        if (S->o.verbose) printf("SANLM filter: %.1f s\n", wall_time() - t0);
    }
    target = S->o.noise * range;
    if (target > sd) add = sqrt(target * target - sd * sd);
    if (add > 0)
    {
        noise_ctx c;
        *est_own = (float *)malloc(sizeof(float) * n);
        if (!*est_own) return 2;
        c.in = src;
        c.out = *est_own;
        c.d = S->dims;
        c.sd = add;
        c.seed = (uint64_t)S->o.noise_seed;
        par_for(S->dims[2], S->nthreads, noise_range, &c);
        S->est = *est_own;
    }
    if (res)
    {
        res->intensity_range[0] = S->range[0];
        res->intensity_range[1] = S->range[1];
        res->noise_estimated = sd / range;
        res->noise_added = add / range;
    }
    if (S->o.verbose)
        printf("Intensity range %g .. %g, noise %.2f%% of it%s\n", S->range[0], S->range[1],
               100.0 * sd / range, (add > 0) ? "" : ", none added");
    if (S->o.verbose && add > 0)
        printf("  noise of %.2f%% added for the estimation\n", 100.0 * add / range);
    return 0;
}

static void
print_model(const State *S)
{
    static const char *name[KB] = {"GM", "WM", "CSF", "bone", "soft", "air"};
    int g;

    printf("  class  weight    mix       mean         sd\n");
    for (g = 0; g < S->G.K; g++)
        printf("  %-5s  %.4f  %.4f  %10.4g  %9.4g\n", name[S->G.lkp[g]], S->G.wp[S->G.lkp[g]],
               S->G.mg[g], S->G.mn[g], sqrt(S->G.vr[g]));
}

/**
 * \brief Unified segmentation of one image with geodesic shooting.
 *
 * SPM's unified segmentation as CAT runs it, with the priors taken from a
 * shooting template series and deformed by geodesic shooting: the
 * Gaussians, the class weights and the bias field are estimated on voxels
 * opts->samp apart (cat_spm_preproc8.m) and alternate with Gauss-Newton
 * steps of the deformation, each of which is an EM step that matches the
 * templates to the class posteriors (CAT_VolShootFit()).  The class
 * posteriors of every voxel follow from the final parameters
 * (cat_spm_preproc_write8.m), the label map from a partial volume model.
 *
 * \param src       (in)  image on the input grid
 * \param dims      (in)  input grid dimensions
 * \param mat       (in)  input grid voxel-to-world matrix
 * \param templates (in)  CAT_SHOOT_NTEMPLATES templates of CAT_USEG_NCLASSES
 *                        volumes
 * \param tdims     (in)  template dimensions
 * \param tmat      (in)  template voxel-to-world matrix
 * \param mask      (in)  optional (NULL): zero marks voxels without data
 * \param affine    (in)  optional (NULL) input-to-template world mapping
 * \param opts_in   (in)  settings; NULL uses CAT_UnifiedSegOptionsInit()
 * \param prob      (out) optional (NULL), CAT_USEG_NCLASSES volumes:
 *                        posteriors times 255
 * \param label     (out) optional (NULL): label map with partial volumes
 * \param pve       (out) optional (NULL), 3 volumes: fractions of GM, WM and
 *                        CSF times 255
 * \param corrected (out) optional (NULL): bias-corrected image
 * \param bias      (out) optional (NULL): bias field
 * \param y         (out) optional (NULL), tdims x 3: deformation
 * \param iy        (out) optional (NULL), dims x 3: its inverse
 * \param jac       (out) optional (NULL), tdims: Jacobian determinants
 * \param res       (out) optional (NULL) parameters and summary
 * \return 0 on success, 1 on invalid input, 2 on memory failure, 3 if the
 *         estimation failed
 */
int
CAT_UnifiedSeg(const float *src, const int dims[3], const double mat[16],
               const float *const templates[CAT_SHOOT_NTEMPLATES],
               const int tdims[3], const double tmat[16],
               const unsigned char *mask, const double affine[16],
               const CAT_UnifiedSegOpts *opts_in, unsigned char *prob,
               float *label, unsigned char *pve, float *corrected,
               float *bias, float *y, float *iy, float *jac,
               CAT_UnifiedSegResult *res)
{
    State S;
    Pve P;
    HeadPve H;
    CAT_VolShootModel model;
    CAT_VolShootInfo info;
    float *clean_own = NULL, *est_own = NULL, *y_own = NULL, *iy_own = NULL;
    double Aff[16], itmat[16], t0 = wall_time(), t1;
    size_t n, ntpl;
    int k, r, ret = 0;

    memset(&S, 0, sizeof(S));
    memset(&P, 0, sizeof(P));
    memset(&H, 0, sizeof(H));
    if (res) memset(res, 0, sizeof(*res));
    if (!src || !dims || !mat || !templates || !tdims || !tmat) return 1;
    if (opts_in) S.o = *opts_in;
    else CAT_UnifiedSegOptionsInit(&S.o);
    for (k = 0; k < NLEV; k++)
        if (!templates[k]) return 1;
    for (k = 0; k < KB; k++)
        if (S.o.ngaus[k] < 1 || S.o.ngaus[k] > CAT_USEG_MAXGAUSS) return 1;
    if (!(S.o.samp > 0) || !(S.o.tol > 0) || !(S.o.biasreg >= 0) || S.o.noise < 0 || S.o.mrf < 0)
        return 1;
    if (!(S.o.wp_reg >= 0) || S.o.cleanup < 0 || S.o.cleanup > 2) return 1;
    if (S.o.bias_model != CAT_USEG_BIAS_DCT && S.o.bias_model != CAT_USEG_BIAS_SPLINE) return 1;
    if (S.o.affine_mode != CAT_USEG_AFFINE_ESTIMATE && !affine) return 1;
    for (k = 0; k < 3; k++)
        if (dims[k] < 2 || tdims[k] < 2) return 1;

    S.nthreads = (S.o.n_threads > 0) ? S.o.n_threads : cat_parallel_nthreads(1 << 20);
    if (S.nthreads > USEG_MAX_THREADS) S.nthreads = USEG_MAX_THREADS;
    S.o.shoot.n_threads = S.nthreads;
    S.o.shoot.verbose = (S.o.verbose > 1);
    S.mask = mask;
    memcpy(S.dims, dims, sizeof(S.dims));
    memcpy(S.tdims, tdims, sizeof(S.tdims));
    memcpy(S.mat, mat, sizeof(S.mat));
    memcpy(S.tmat, tmat, sizeof(S.tmat));
    for (k = 0; k < 3; k++)
        S.vx[k] = sqrt(mat[k] * mat[k] + mat[4 + k] * mat[4 + k] + mat[8 + k] * mat[8 + k]);
    n = (size_t)dims[0] * dims[1] * dims[2];
    ntpl = (size_t)tdims[0] * tdims[1] * tdims[2];
    if (m4_inv(tmat, itmat)) return 1;

    memcpy(S.tpl.d, tdims, sizeof(S.tpl.d));
    S.tpl.n = ntpl;
    S.tpl.tpl = templates;
    S.tpl.level = -1;
    S.tpl.lt = (float *)malloc(sizeof(float) * KB * ntpl);
    if (!S.tpl.lt)
    {
        ret = 2;
        goto done;
    }

    /* 1. the images */
    if ((ret = prepare_images(&S, src, &clean_own, &est_own, res))) goto done;

    /* 2. affine registration */
    t1 = wall_time();
    if (S.o.affine_mode == CAT_USEG_AFFINE_FIXED)
        memcpy(Aff, affine, sizeof(Aff));
    else
    {
        double ll = 0.0;
        if (S.o.verbose) printf("Affine registration\n");
        if ((ret = affine_register(&S, affine, Aff, &ll))) goto done;
        if (res) res->affine_ll = ll;
    }
    m4_mul(itmat, Aff, S.Tvox);
    m4_mul(S.Tvox, mat, S.Tvox);
    if (res)
    {
        memcpy(res->affine, Aff, sizeof(Aff));
        res->time_affine = wall_time() - t1;
    }

    /* 3. the parameters of the image and the deformation */
    t1 = wall_time();
    {
        /* 5 mm in template voxels, for 0-based coordinates */
        double tvz = sqrt(tmat[2] * tmat[2] + tmat[6] * tmat[6] + tmat[10] * tmat[10]);
        if ((ret = lattice_build(&S.L, S.est, S.clean, mask, dims, S.vx, S.Tvox, 5.0 / tvz - 1.0,
                                 S.o.samp, S.o.bias_model, S.o.biasfwhm, S.o.biasreg,
                                 S.nthreads)))
            goto done;
    }
    S.acc = (double *)malloc(sizeof(double) * ACC_N * S.L.d[2]);
    if (!y) y = y_own = (float *)malloc(sizeof(float) * 3 * ntpl);
    if (!iy) iy = iy_own = (float *)malloc(sizeof(float) * 3 * n);
    if (!S.acc || !y || !iy)
    {
        ret = 2;
        goto done;
    }
    if (S.o.verbose)
        printf("Estimation from %zu voxels (every %d x %d x %d), %d x %d x %d %s bias functions\n",
               S.L.nm, S.L.sk[0], S.L.sk[1], S.L.sk[2], S.L.nb[0], S.L.nb[1], S.L.nb[2],
               (S.L.model == CAT_USEG_BIAS_DCT) ? "DCT" : "spline");
    model.n_cls = KB;
    model.ctx = &S;
    model.update = model_update;
    ret = CAT_VolShootFit(&model, dims, mat, Aff, templates, tdims, tmat, 0.0, &S.o.shoot,
                          y, iy, jac, &info);
    if (ret) goto done;

    /* the parameters once more, with the final template and deformation */
    tpl_load(&S.tpl, NLEV - 1, S.nthreads);
    {
        Lattice *L = &S.L;
        size_t nxy = (size_t)dims[0] * dims[1];
        int x, yy, z;
        for (z = 0; z < L->d[2]; z++)
            for (yy = 0; yy < L->d[1]; yy++)
                for (x = 0; x < L->d[0]; x++)
                {
                    size_t l = x + (size_t)L->d[0] * (yy + (size_t)L->d[1] * z);
                    size_t v = (size_t)x * L->sk[0] + (size_t)dims[0] * ((size_t)yy * L->sk[1]) +
                               nxy * ((size_t)z * L->sk[2]);
                    for (r = 0; r < 3; r++) L->ty[r * L->n + l] = iy[r * n + v];
                }
    }
    lattice_priors(&S);
    if ((ret = estimate(&S, 8, 0))) goto done;
    if (S.o.verbose)
    {
        printf("Shooting: %d iterations, Jacobian determinants %.3f .. %.3f\n", info.iterations,
               info.det_min, info.det_max);
        printf("Log-likelihood per voxel %.6f\n", S.ll / (double)S.L.nm);
        printf("Model updates: priors %.1f, EM %.1f, bias %.1f (%d steps, %d CG iterations), "
               "likelihoods %.1f s\n", S.t_prior, S.t_em, S.t_bias, S.n_bias, S.n_cg, S.t_fill);
        print_model(&S);
    }
    if (res)
    {
        res->shoot = info;
        res->n_gauss = S.G.K;
        for (k = 0; k < S.G.K; k++)
        {
            res->lkp[k] = S.G.lkp[k];
            res->mg[k] = S.G.mg[k];
            res->mn[k] = S.G.mn[k];
            res->vr[k] = S.G.vr[k];
        }
        for (k = 0; k < KB; k++) res->wp[k] = S.G.wp[k];
        res->ll = S.ll / (double)S.L.nm;
        res->n_samples = (int)S.L.nm;
        res->bias_model = S.L.model;
        for (k = 0; k < 3; k++) res->bias_dims[k] = S.L.nb[k];
        if (S.L.nbt)
        {
            res->bias_coef = (double *)malloc(sizeof(double) * S.L.nbt);
            if (!res->bias_coef)
            {
                ret = 2;
                goto done;
            }
            memcpy(res->bias_coef, S.L.T, sizeof(double) * S.L.nbt);
        }
        res->time_fit = wall_time() - t1;
    }

    /* 4. the segmentation of every voxel */
    t1 = wall_time();
    if (S.o.pve && (label || pve))
    {
        if ((ret = pve_setup(&S, &P))) goto done;
        if (res && P.ok)
        {
            res->pve_done = 1;
            for (k = 0; k < 3; k++) res->pve_mean[k] = P.m[k];
            res->pve_sd = sqrt(P.v);
        }
        if (S.o.verbose && P.ok)
            printf("Partial volumes: CSF %.4g, GM %.4g, WM %.4g, SD %.4g\n", P.m[0], P.m[1], P.m[2],
                   sqrt(P.v));
        if (S.o.verbose && !P.ok) printf("Partial volumes: tissues not separated, skipped\n");
    }
    if (S.o.pve_head && (label || pve))
    {
        if ((ret = head_pve_setup(&S, &H))) goto done;
        if (res) res->pve_head_done = H.ok;
        if (S.o.verbose && !H.ok) printf("Partial volumes of the head: skipped\n");
    }
    ret = write_segmentation(&S, iy, &P, &H, prob, label, pve, corrected, bias, res);
    if (res)
    {
        res->time_write = wall_time() - t1;
        res->time_total = wall_time() - t0;
    }
    if (S.o.verbose && res)
        printf("Time: filter %.1f, affine %.1f, estimation %.1f, segmentation %.1f, total %.1f s\n",
               res->time_filter, res->time_affine, res->time_fit, res->time_write, res->time_total);

done:
    if (ret && res) CAT_UnifiedSegResultFree(res);
    free(P.dens);
    free(P.frac);
    free(H.dens);
    free(H.frac);
    lattice_free(&S.L);
    free(S.acc);
    free(S.ireg);
    free(S.tpl.lt);
    free(clean_own);
    free(est_own);
    free(y_own);
    free(iy_own);
    return ret;
}

/* ======================================================================== */
/* JSON                                                                     */
/* ======================================================================== */

static void
json_string(FILE *fp, const char *s)
{
    fputc('"', fp);
    for (; s && *s; s++)
    {
        if (*s == '"' || *s == '\\') fputc('\\', fp);
        if ((unsigned char)*s < 0x20) fprintf(fp, "\\u%04x", (unsigned char)*s);
        else fputc(*s, fp);
    }
    fputc('"', fp);
}

static void
json_array(FILE *fp, const double *v, int n)
{
    int i;
    fputc('[', fp);
    for (i = 0; i < n; i++)
    {
        if (isfinite(v[i])) fprintf(fp, "%s%.9g", i ? ", " : "", v[i]);
        else fprintf(fp, "%snull", i ? ", " : "");
    }
    fputc(']', fp);
}

/**
 * \brief Write the parameters of a segmentation as JSON.
 *
 * What SPM stores in its *_seg8.mat -- the affine matrix, the Gaussians,
 * the class weights and the bias field coefficients -- with the settings
 * and the summary of the run.
 *
 * \param file  (in) output file
 * \param res   (in) result of CAT_UnifiedSeg()
 * \param opts  (in) the settings it ran with; NULL: the defaults
 * \param image (in) optional (NULL) name of the segmented image
 * \param tpl   (in) optional (NULL) name of the first template
 * \return 0 on success, 1 if the file cannot be written
 */
int
CAT_UnifiedSegWriteJson(const char *file, const CAT_UnifiedSegResult *res,
                        const CAT_UnifiedSegOpts *opts, const char *image,
                        const char *tpl)
{
    static const char *affreg[] = {"none", "mni", "eastern", "subj", "rigid"};
    static const char *cls[KB] = {"GM", "WM", "CSF", "bone", "soft tissue", "air"};
    CAT_UnifiedSegOpts d;
    double v[MAXK];
    FILE *fp;
    int i, k, nb;

    if (!file || !res) return 1;
    if (!opts)
    {
        CAT_UnifiedSegOptionsInit(&d);
        opts = &d;
    }
    fp = fopen(file, "w");
    if (!fp) return 1;
    fprintf(fp, "{\n  \"tool\": \"CAT_VolUnifiedSeg\",\n  \"image\": ");
    json_string(fp, image ? image : "");
    fprintf(fp, ",\n  \"template\": ");
    json_string(fp, tpl ? tpl : "");
    fprintf(fp, ",\n  \"classes\": [");
    for (k = 0; k < KB; k++)
    {
        if (k) fprintf(fp, ", ");
        json_string(fp, cls[k]);
    }
    fprintf(fp, "],\n  \"affine\": [\n");
    for (i = 0; i < 4; i++)
    {
        fprintf(fp, "    ");
        json_array(fp, res->affine + 4 * i, 4);
        fprintf(fp, "%s\n", (i < 3) ? "," : "");
    }
    fprintf(fp, "  ],\n  \"affine_mutual_information\": %.6g,\n", res->affine_ll);
    for (k = 0; k < res->n_gauss; k++) v[k] = res->lkp[k] + 1;
    fprintf(fp, "  \"lkp\": ");
    json_array(fp, v, res->n_gauss);
    fprintf(fp, ",\n  \"mg\": ");
    json_array(fp, res->mg, res->n_gauss);
    fprintf(fp, ",\n  \"mn\": ");
    json_array(fp, res->mn, res->n_gauss);
    fprintf(fp, ",\n  \"vr\": ");
    json_array(fp, res->vr, res->n_gauss);
    fprintf(fp, ",\n  \"wp\": ");
    json_array(fp, res->wp, KB);
    fprintf(fp, ",\n  \"ll\": %.9g,\n  \"n_samples\": %d,\n", res->ll, res->n_samples);
    nb = res->bias_coef ? res->bias_dims[0] * res->bias_dims[1] * res->bias_dims[2] : 0;
    fprintf(fp, "  \"bias\": {\"model\": \"%s\", \"fwhm\": %g, \"reg\": %g, "
                "\"dims\": [%d, %d, %d], \"coef\": ",
            (res->bias_model == CAT_USEG_BIAS_DCT) ? "dct" : "spline", opts->biasfwhm,
            opts->biasreg, nb ? res->bias_dims[0] : 0, nb ? res->bias_dims[1] : 0,
            nb ? res->bias_dims[2] : 0);
    json_array(fp, res->bias_coef, nb);
    fprintf(fp, "},\n  \"intensity_range\": ");
    json_array(fp, res->intensity_range, 2);
    fprintf(fp, ",\n  \"noise\": {\"estimated\": %.6g, \"added\": %.6g, \"target\": %g, \"seed\": %u},\n",
            res->noise_estimated, res->noise_added, opts->noise, opts->noise_seed);
    fprintf(fp, "  \"pve\": {\"done\": %s, \"mean\": ", res->pve_done ? "true" : "false");
    json_array(fp, res->pve_mean, 3);
    fprintf(fp, ", \"sd\": %.9g, \"head\": %s},\n  \"volume_ml\": ", res->pve_sd,
            res->pve_head_done ? "true" : "false");
    json_array(fp, res->volume, KB);
    fprintf(fp, ",\n  \"volume_pve_ml\": ");
    json_array(fp, res->pve_volume, KB);
    fprintf(fp, ",\n  \"tiv_ml\": %.9g,\n", res->volume[C_GM] + res->volume[C_WM] + res->volume[C_CSF]);
    fprintf(fp, "  \"shooting\": {\"iterations\": %d, \"level_iterations\": [%d, %d, %d, %d, %d], "
                "\"det_min\": %.6g, \"det_max\": %.6g, \"hotspot_corrections\": %d},\n",
            res->shoot.iterations, res->shoot.level_iterations[0], res->shoot.level_iterations[1],
            res->shoot.level_iterations[2], res->shoot.level_iterations[3],
            res->shoot.level_iterations[4], res->shoot.det_min, res->shoot.det_max,
            res->shoot.hotspot_corrections);
    for (k = 0; k < KB; k++) v[k] = opts->ngaus[k];
    fprintf(fp, "  \"options\": {\"ngaus\": ");
    json_array(fp, v, KB);
    fprintf(fp, ", \"samp\": %g, \"tol\": %g, \"affreg\": \"%s\", \"sanlm\": %s, \"mrf\": %g, "
                "\"cleanup\": %d, \"pve\": %s, \"pve_head\": %s, \"regstr\": %g, \"reg_scale\": %g, "
                "\"lowres\": %g},\n",
            opts->samp, opts->tol,
            affreg[(opts->affreg >= 0 && opts->affreg <= 4) ? opts->affreg : 0],
            opts->sanlm ? "true" : "false", opts->mrf, opts->cleanup, opts->pve ? "true" : "false",
            opts->pve_head ? "true" : "false", opts->shoot.regstr, opts->shoot.reg_scale,
            opts->shoot.lowres);
    fprintf(fp, "  \"time_s\": {\"filter\": %.2f, \"affine\": %.2f, \"estimation\": %.2f, "
                "\"segmentation\": %.2f, \"total\": %.2f}\n}\n",
            res->time_filter, res->time_affine, res->time_fit, res->time_write, res->time_total);
    return fclose(fp) ? 1 : 0;
}
