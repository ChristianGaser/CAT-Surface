/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Head-motion correction of a 4D series, reported in mcflirt's conventions.
 *
 * Copyright Christian Gaser, University of Jena.
 */

#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#if defined(_WIN32) || defined(_WIN64)
#include <direct.h>
#endif

#include "CAT_Realign.h"
#include "CAT_BBReg.h"
#include "CAT_Vol.h"

/* radius of the sphere of mcflirt's RMS deviation (mm) */
#define RMS_RADIUS 80.0

/* =========================================================================
 * 4x4 row-major matrix helpers
 * ========================================================================= */

static void mat4_identity(double M[16])
{
    int i;

    for (i = 0; i < 16; i++)
        M[i] = (i % 5 == 0) ? 1.0 : 0.0;
}

/* C = A * B; C may not alias A or B */
static void mat4_mul(const double A[16], const double B[16], double C[16])
{
    int r, c, k;

    for (r = 0; r < 4; r++)
        for (c = 0; c < 4; c++)
        {
            double s = 0.0;
            for (k = 0; k < 4; k++)
                s += A[4 * r + k] * B[4 * k + c];
            C[4 * r + c] = s;
        }
}

/* inverse of an affine matrix (last row 0 0 0 1) */
static void affine_inverse(const double A[16], double inv[16])
{
    double a = A[0], b = A[1], c = A[2], d = A[4], e = A[5], f = A[6];
    double g = A[8], h = A[9], k = A[10];
    double det = a * (e * k - f * h) - b * (d * k - f * g) + c * (d * h - e * g);
    double R[9];
    int r;

    R[0] = (e * k - f * h) / det;
    R[1] = (c * h - b * k) / det;
    R[2] = (b * f - c * e) / det;
    R[3] = (f * g - d * k) / det;
    R[4] = (a * k - c * g) / det;
    R[5] = (c * d - a * f) / det;
    R[6] = (d * h - e * g) / det;
    R[7] = (b * g - a * h) / det;
    R[8] = (a * e - b * d) / det;
    for (r = 0; r < 3; r++)
    {
        inv[4 * r] = R[3 * r];
        inv[4 * r + 1] = R[3 * r + 1];
        inv[4 * r + 2] = R[3 * r + 2];
        inv[4 * r + 3] = -(R[3 * r] * A[3] + R[3 * r + 1] * A[7] +
                           R[3 * r + 2] * A[11]);
    }
    inv[12] = inv[13] = inv[14] = 0.0;
    inv[15] = 1.0;
}

static void vox2world(const nifti_image *nii, double V[16])
{
    int r, c;

    for (r = 0; r < 4; r++)
        for (c = 0; c < 4; c++)
            V[4 * r + c] = (double)nii->sto_xyz.m[r][c];
}

/* =========================================================================
 * mcflirt's conventions
 * ========================================================================= */

/**
 * \brief Centre of mass of a volume in FSL coordinates, as FSL computes it.
 *
 * FSL's cog() weights every voxel with its intensity minus the minimum of the
 * volume; without the subtraction a reference with an offset background
 * moves the centre, and with it every translation of the .par file (by
 * 0.05 mm for an offset of 1000 in a test).
 *
 * \param vol  (in)  volume (float, x fastest)
 * \param nii  (in)  its header
 * \param dims (in)  [nx, ny, nz]
 * \param com  (out) centre of mass, FSL mm
 */
void CAT_VolRealignCentreOfMass(const float *vol, const nifti_image *nii,
                                const int dims[3], double com[3])
{
    size_t n = (size_t)dims[0] * dims[1] * dims[2], i;
    double sw = 0.0, sx = 0.0, sy = 0.0, sz = 0.0, F[16], c[3];
    float vmin = vol[0];
    int x, y, z, r;

    for (i = 1; i < n; i++)
        if (vol[i] < vmin)
            vmin = vol[i];
    i = 0;
    for (z = 0; z < dims[2]; z++)
        for (y = 0; y < dims[1]; y++)
            for (x = 0; x < dims[0]; x++, i++)
            {
                double w = (double)vol[i] - (double)vmin;
                sw += w;
                sx += w * x;
                sy += w * y;
                sz += w * z;
            }
    if (sw > 0.0)
    {
        c[0] = sx / sw;
        c[1] = sy / sw;
        c[2] = sz / sw;
    }
    else
    {
        c[0] = 0.5 * (dims[0] - 1);
        c[1] = 0.5 * (dims[1] - 1);
        c[2] = 0.5 * (dims[2] - 1);
    }

    vox2fsl_matrix(nii, F);
    for (r = 0; r < 3; r++)
        com[r] = F[4 * r] * c[0] + F[4 * r + 1] * c[1] + F[4 * r + 2] * c[2] +
                 F[4 * r + 3];
}

/**
 * \brief mcflirt's six parameters of an FSL transform.
 *
 * The angles follow FSL's rotmat2euler() for R = Rx Ry Rz, including its
 * branch for cos(ry) ~ 0; the translation is R * centre + t - centre, which
 * is what decompose_aff() reports for a rigid matrix.
 *
 * \param fsl    (in)  4x4 row-major rigid transform in FSL coordinates
 * \param centre (in)  rotation centre, FSL mm
 * \param par    (out) rx ry rz (rad) tx ty tz (mm)
 */
void CAT_VolRealignFslParams(const double fsl[16], const double centre[3],
                             double par[6])
{
    double cy = sqrt(fsl[0] * fsl[0] + fsl[1] * fsl[1]);
    int i;

    if (cy < 1e-4)
    {
        par[0] = atan2(-fsl[9], fsl[5]);
        par[1] = atan2(-fsl[2], 0.0);
        par[2] = 0.0;
    }
    else
    {
        par[0] = atan2(fsl[6] / cy, fsl[10] / cy);
        par[1] = atan2(-fsl[2], cy);
        par[2] = atan2(fsl[1] / cy, fsl[0] / cy);
    }
    for (i = 0; i < 3; i++)
        par[3 + i] = fsl[4 * i] * centre[0] + fsl[4 * i + 1] * centre[1] +
                     fsl[4 * i + 2] * centre[2] + fsl[4 * i + 3] - centre[i];
}

/**
 * \brief Jenkinson's RMS deviation between two transforms.
 *
 * \param a      (in) 4x4 row-major transform
 * \param b      (in) 4x4 row-major transform
 * \param centre (in) centre of the 80 mm sphere, in the coordinates of a and b
 * \return RMS deviation (mm)
 */
double CAT_VolRealignRms(const double a[16], const double b[16],
                         const double centre[3])
{
    double b_inv[16], D[16], t[3], tr = 0.0;
    int r, c;

    affine_inverse(b, b_inv);
    mat4_mul(a, b_inv, D);
    for (r = 0; r < 3; r++)
        D[5 * r] -= 1.0;
    for (r = 0; r < 3; r++)
        t[r] = D[4 * r + 3] + D[4 * r] * centre[0] + D[4 * r + 1] * centre[1] +
               D[4 * r + 2] * centre[2];
    for (r = 0; r < 3; r++)
        for (c = 0; c < 3; c++)
            tr += D[4 * r + c] * D[4 * r + c];
    return sqrt(t[0] * t[0] + t[1] * t[1] + t[2] * t[2] +
                RMS_RADIUS * RMS_RADIUS / 5.0 * tr);
}

/* =========================================================================
 * Estimation
 * ========================================================================= */

/**
 * \brief Set the default options of the motion correction.
 *
 * \param opts (out) options to initialise
 */
void CAT_VolRealignOptionsInit(CAT_VolRealignOpts *opts)
{
    opts->refvol = -1;
    CAT_VolumeRegOptionsInit(&opts->reg);
    opts->reg.n_levels = 3;
    opts->verbose = 0;
}

/* One direction of the chain from the middle volume: volumes first,
 * first + step, ..., up to but excluding last, each starting from the
 * estimate of the previous one. */
typedef struct
{
    const float *series;
    const nifti_image *nii;
    const int *dims;
    size_t nvox;
    const CAT_VolPyramid *ref_pyr;
    const CAT_VolumeRegOpts *reg;
    CAT_RigidParams *p;
    double *scale;
    int first, last, step;
    int verbose;
    int failed;
} chain_args;

static void *chain_worker(void *arg)
{
    chain_args *a = (chain_args *)arg;
    int prev = a->first - a->step, t;

    for (t = a->first; t != a->last; t += a->step)
    {
        CAT_VolPyramid *mp = CAT_VolumeReg_pyramid(a->series + (size_t)t * a->nvox,
                                                   a->nii, a->dims, a->reg->n_levels);
        CAT_RigidParams p = a->p[prev];
        double s = a->scale[prev];

        if (!mp)
        {
            a->failed = 1;
            return NULL;
        }
        CAT_VolumeReg_register_pyramids(a->ref_pyr, mp, &p, &s, a->reg);
        CAT_VolumeReg_pyramid_free(mp);
        a->p[t] = p;
        a->scale[t] = s;
        if (a->verbose)
            printf("volume %4d  t %8.4f %8.4f %8.4f  r %9.5f %9.5f %9.5f"
                   "  scale %.4f\n", t, p.tx, p.ty, p.tz, p.rx, p.ry, p.rz, s);
        prev = t;
    }
    return NULL;
}

/**
 * \brief Estimate the motion of every volume of a series.
 *
 * The registration of a volume maps the reference into it (CAT_VolumeReg's
 * fixed -> moving direction); the results are inverted to mcflirt's volume ->
 * reference direction, converted to FSL coordinates, and decomposed into the
 * parameters and RMS deviations.
 *
 * \param series     (in)  nt volumes of dims[0]*dims[1]*dims[2] floats
 * \param series_nii (in)  header of the series
 * \param dims       (in)  [nx, ny, nz] of one volume
 * \param nt         (in)  number of volumes, at least 1
 * \param ref        (in)  separate reference volume, or NULL for volume
 *                         opts->refvol of the series
 * \param ref_nii    (in)  header of the separate reference
 * \param ref_dims   (in)  [nx, ny, nz] of the separate reference
 * \param opts       (in)  options; NULL for CAT_VolRealignOptionsInit()
 * \param res        (out) estimates; release with CAT_VolRealignFree()
 * \return 0 on success, -1 on error
 */
int CAT_VolRealignEstimate(const float *series, const nifti_image *series_nii,
                           const int dims[3], int nt,
                           const float *ref, const nifti_image *ref_nii,
                           const int ref_dims[3],
                           const CAT_VolRealignOpts *opts,
                           CAT_VolRealignResult *res)
{
    CAT_VolRealignOpts defaults;
    size_t nvox = (size_t)dims[0] * dims[1] * dims[2];
    CAT_VolPyramid *ref_pyr;
    CAT_RigidParams *p;
    chain_args args[2];
    double com[3], centre[3], I[16];
    int c, t, ri;

    memset(res, 0, sizeof(*res));
    if (!opts)
    {
        CAT_VolRealignOptionsInit(&defaults);
        opts = &defaults;
    }
    if (nt < 1)
    {
        fprintf(stderr, "CAT_VolRealignEstimate: no volumes\n");
        return -1;
    }

    if (ref)
        ri = -1;
    else
    {
        ri = (opts->refvol >= 0 && opts->refvol < nt) ? opts->refvol : nt / 2;
        if (opts->refvol >= nt)
            fprintf(stderr, "CAT_VolRealignEstimate: reference volume %d does not"
                            " exist, using %d\n", opts->refvol, ri);
        ref = series + (size_t)ri * nvox;
        ref_nii = series_nii;
        ref_dims = dims;
    }

    res->nt = nt;
    res->ref_index = ri;
    res->ras = (double *)malloc(sizeof(double) * 16 * nt);
    res->fsl = (double *)malloc(sizeof(double) * 16 * nt);
    res->par = (double *)malloc(sizeof(double) * 6 * nt);
    res->rms_abs = (double *)malloc(sizeof(double) * nt);
    res->rms_rel = (double *)malloc(sizeof(double) * (nt > 1 ? nt - 1 : 1));
    res->scale = (double *)malloc(sizeof(double) * nt);
    p = (CAT_RigidParams *)calloc((size_t)nt, sizeof(CAT_RigidParams));
    ref_pyr = CAT_VolumeReg_pyramid(ref, ref_nii, ref_dims, opts->reg.n_levels);
    if (!res->ras || !res->fsl || !res->par || !res->rms_abs || !res->rms_rel ||
        !res->scale || !p || !ref_pyr)
    {
        fprintf(stderr, "CAT_VolRealignEstimate: out of memory\n");
        CAT_VolumeReg_pyramid_free(ref_pyr);
        free(p);
        CAT_VolRealignFree(res);
        return -1;
    }

    /* the middle volume first: the reference itself is the identity, any
       other starts from it and fits its intensity scale */
    c = (ri >= 0) ? ri : nt / 2;
    res->scale[c] = 1.0;
    if (ri < 0)
    {
        CAT_VolPyramid *mp = CAT_VolumeReg_pyramid(series + (size_t)c * nvox,
                                                   series_nii, dims,
                                                   opts->reg.n_levels);
        if (!mp)
        {
            fprintf(stderr, "CAT_VolRealignEstimate: out of memory\n");
            CAT_VolumeReg_pyramid_free(ref_pyr);
            free(p);
            CAT_VolRealignFree(res);
            return -1;
        }
        res->scale[c] = 0.0;
        CAT_VolumeReg_register_pyramids(ref_pyr, mp, &p[c], &res->scale[c],
                                        &opts->reg);
        CAT_VolumeReg_pyramid_free(mp);
    }
    if (opts->verbose)
        printf("volume %4d  t %8.4f %8.4f %8.4f  r %9.5f %9.5f %9.5f"
               "  scale %.4f%s\n", c, p[c].tx, p[c].ty, p[c].tz, p[c].rx,
               p[c].ry, p[c].rz, res->scale[c], ri >= 0 ? "  (reference)" : "");

    /* then outward in both directions, one thread each */
    for (t = 0; t < 2; t++)
    {
        args[t].series = series;
        args[t].nii = series_nii;
        args[t].dims = dims;
        args[t].nvox = nvox;
        args[t].ref_pyr = ref_pyr;
        args[t].reg = &opts->reg;
        args[t].p = p;
        args[t].scale = res->scale;
        args[t].verbose = opts->verbose;
        args[t].failed = 0;
    }
    args[0].first = c + 1;
    args[0].last = nt;
    args[0].step = 1;
    args[1].first = c - 1;
    args[1].last = -1;
    args[1].step = -1;
    cat_parallel_run(2, chain_worker, args, sizeof(chain_args));
    CAT_VolumeReg_pyramid_free(ref_pyr);
    if (args[0].failed || args[1].failed)
    {
        fprintf(stderr, "CAT_VolRealignEstimate: out of memory\n");
        free(p);
        CAT_VolRealignFree(res);
        return -1;
    }

    /* volume -> reference, in world space and in FSL's convention */
    mat4_identity(I);
    for (t = 0; t < nt; t++)
    {
        double M[16];

        if (t == ri)
        {
            memcpy(res->ras + 16 * t, I, sizeof(I));
            memcpy(res->fsl + 16 * t, I, sizeof(I));
            continue;
        }
        CAT_BBReg_params_to_matrix(&p[t], M);
        CAT_BBReg_invert_matrix(M, res->ras + 16 * t);
        ras_to_fsl_matrix(res->ras + 16 * t, series_nii, ref_nii, res->fsl + 16 * t);
    }
    free(p);

    /* mcflirt's parameters and RMS deviations */
    CAT_VolRealignCentreOfMass(ref, ref_nii, ref_dims, com);
    centre[0] = 0.5 * (ref_dims[0] - 1) * fabs(ref_nii->dx);
    centre[1] = 0.5 * (ref_dims[1] - 1) * fabs(ref_nii->dy);
    centre[2] = 0.5 * (ref_dims[2] - 1) * fabs(ref_nii->dz);
    res->rms_abs_mean = res->rms_rel_mean = 0.0;
    for (t = 0; t < nt; t++)
    {
        CAT_VolRealignFslParams(res->fsl + 16 * t, com, res->par + 6 * t);
        res->rms_abs[t] = CAT_VolRealignRms(res->fsl + 16 * t, I, centre);
        res->rms_abs_mean += res->rms_abs[t] / nt;
        if (t > 0)
        {
            res->rms_rel[t - 1] = CAT_VolRealignRms(res->fsl + 16 * t,
                                                    res->fsl + 16 * (t - 1), centre);
            res->rms_rel_mean += res->rms_rel[t - 1] / (nt - 1);
        }
    }
    return 0;
}

/* =========================================================================
 * Resampling
 * ========================================================================= */

/* Trilinear sample at voxel position (x, y, z); 0 outside [0, n-1].  The
 * last voxel along an axis is interpolated from itself, so an identity
 * transform reproduces the volume. */
static float sample_edge(const float *vol, int nx, int ny, int nz,
                         double x, double y, double z)
{
    const double eps = 1e-6;
    int x0, y0, z0, x1, y1, z1;
    double wx, wy, wz;

    if (x < -eps || y < -eps || z < -eps ||
        x > nx - 1 + eps || y > ny - 1 + eps || z > nz - 1 + eps)
        return 0.0f;
    x = (x < 0.0) ? 0.0 : (x > nx - 1 ? nx - 1 : x);
    y = (y < 0.0) ? 0.0 : (y > ny - 1 ? ny - 1 : y);
    z = (z < 0.0) ? 0.0 : (z > nz - 1 ? nz - 1 : z);
    x0 = (int)x; y0 = (int)y; z0 = (int)z;
    x1 = (x0 < nx - 1) ? x0 + 1 : x0;
    y1 = (y0 < ny - 1) ? y0 + 1 : y0;
    z1 = (z0 < nz - 1) ? z0 + 1 : z0;
    wx = x - x0; wy = y - y0; wz = z - z0;
#define V(a, b, c) ((double)vol[(a) + (size_t)nx * ((b) + (size_t)ny * (c))])
    return (float)((V(x0, y0, z0) * (1 - wx) + V(x1, y0, z0) * wx) * (1 - wy) * (1 - wz) +
                   (V(x0, y1, z0) * (1 - wx) + V(x1, y1, z0) * wx) * wy * (1 - wz) +
                   (V(x0, y0, z1) * (1 - wx) + V(x1, y0, z1) * wx) * (1 - wy) * wz +
                   (V(x0, y1, z1) * (1 - wx) + V(x1, y1, z1) * wx) * wy * wz);
#undef V
}

typedef struct
{
    const float *series;
    const int *dims;
    const int *ref_dims;
    const double *vox_map; /* nt voxel-to-voxel maps, reference -> volume */
    float *out;
    int t0, t1;
} resample_args;

static void *resample_worker(void *arg)
{
    resample_args *a = (resample_args *)arg;
    size_t nvox = (size_t)a->dims[0] * a->dims[1] * a->dims[2];
    size_t nref = (size_t)a->ref_dims[0] * a->ref_dims[1] * a->ref_dims[2];
    int t, x, y, z;

    for (t = a->t0; t < a->t1; t++)
    {
        const double *A = a->vox_map + 16 * t;
        const float *vol = a->series + (size_t)t * nvox;
        float *o = a->out + (size_t)t * nref;
        size_t i = 0;

        for (z = 0; z < a->ref_dims[2]; z++)
            for (y = 0; y < a->ref_dims[1]; y++)
                for (x = 0; x < a->ref_dims[0]; x++, i++)
                    o[i] = sample_edge(vol, a->dims[0], a->dims[1], a->dims[2],
                                       A[0] * x + A[1] * y + A[2] * z + A[3],
                                       A[4] * x + A[5] * y + A[6] * z + A[7],
                                       A[8] * x + A[9] * y + A[10] * z + A[11]);
    }
    return NULL;
}

/**
 * \brief Resample every volume into the reference grid.
 *
 * \param series     (in)  the series given to CAT_VolRealignEstimate()
 * \param series_nii (in)  its header
 * \param dims       (in)  [nx, ny, nz] of one volume
 * \param nt         (in)  number of volumes
 * \param ref_nii    (in)  header of the reference grid
 * \param ref_dims   (in)  [nx, ny, nz] of the reference grid
 * \param res        (in)  estimates from CAT_VolRealignEstimate()
 * \param out        (out) nt volumes on the reference grid
 * \return 0 on success, -1 if out of memory
 */
int CAT_VolRealignResample(const float *series, const nifti_image *series_nii,
                           const int dims[3], int nt,
                           const nifti_image *ref_nii, const int ref_dims[3],
                           const CAT_VolRealignResult *res, float *out)
{
    double Vs[16], Vs_inv[16], Vr[16], *vox_map;
    resample_args args[64];
    int nthreads, t;

    vox_map = (double *)malloc(sizeof(double) * 16 * nt);
    if (!vox_map)
    {
        fprintf(stderr, "CAT_VolRealignResample: out of memory\n");
        return -1;
    }

    /* reference voxel -> reference world -> volume world -> volume voxel */
    vox2world(series_nii, Vs);
    affine_inverse(Vs, Vs_inv);
    vox2world(ref_nii, Vr);
    for (t = 0; t < nt; t++)
    {
        double ref2vol[16], tmp[16];

        affine_inverse(res->ras + 16 * t, ref2vol);
        mat4_mul(ref2vol, Vr, tmp);
        mat4_mul(Vs_inv, tmp, vox_map + 16 * t);
    }

    nthreads = cat_parallel_nthreads(nt);
    if (nthreads > (int)(sizeof(args) / sizeof(args[0])))
        nthreads = (int)(sizeof(args) / sizeof(args[0]));
    for (t = 0; t < nthreads; t++)
    {
        args[t].series = series;
        args[t].dims = dims;
        args[t].ref_dims = ref_dims;
        args[t].vox_map = vox_map;
        args[t].out = out;
        args[t].t0 = (int)((long)nt * t / nthreads);
        args[t].t1 = (int)((long)nt * (t + 1) / nthreads);
    }
    cat_parallel_run(nthreads, resample_worker, args, sizeof(resample_args));
    free(vox_map);
    return 0;
}

/* =========================================================================
 * Output
 * ========================================================================= */

static int make_dir(const char *path)
{
#if defined(_WIN32) || defined(_WIN64)
    if (_mkdir(path) == 0 || errno == EEXIST)
        return 0;
#else
    if (mkdir(path, 0755) == 0 || errno == EEXIST)
        return 0;
#endif
    fprintf(stderr, "CAT_VolRealignWrite: cannot create directory '%s'\n", path);
    return -1;
}

/* write n values, one "%g" per line */
static int write_values_file(const char *file, const double *v, int n)
{
    FILE *fp = fopen(file, "w");
    int i;

    if (!fp)
    {
        fprintf(stderr, "CAT_VolRealignWrite: cannot write '%s'\n", file);
        return -1;
    }
    for (i = 0; i < n; i++)
        fprintf(fp, "%g\n", v[i]);
    return fclose(fp) == 0 ? 0 : -1;
}

/**
 * \brief Write mcflirt's text outputs next to an output base name.
 *
 * \param base   (in) output base name (mcflirt's -out)
 * \param res    (in) estimates from CAT_VolRealignEstimate()
 * \param mats   (in) non-zero: write <base>.mat/MAT_%04d
 * \param plots  (in) non-zero: write <base>.par
 * \param rmsabs (in) non-zero: write <base>_abs.rms and <base>_abs_mean.rms
 * \param rmsrel (in) non-zero: write <base>_rel.rms and <base>_rel_mean.rms
 * \return 0 on success, -1 if a file could not be written
 */
int CAT_VolRealignWrite(const char *base, const CAT_VolRealignResult *res,
                        int mats, int plots, int rmsabs, int rmsrel)
{
    size_t len = strlen(base) + 32;
    char *name = (char *)malloc(len);
    int t, ret = 0;

    if (!name)
        return -1;

    if (mats)
    {
        snprintf(name, len, "%s.mat", base);
        ret = make_dir(name);
        for (t = 0; ret == 0 && t < res->nt; t++)
        {
            snprintf(name, len, "%s.mat/MAT_%04d", base, t);
            ret = write_fsl_matrix(name, res->fsl + 16 * t);
        }
    }
    if (ret == 0 && plots)
    {
        FILE *fp;

        snprintf(name, len, "%s.par", base);
        fp = fopen(name, "w");
        if (!fp)
        {
            fprintf(stderr, "CAT_VolRealignWrite: cannot write '%s'\n", name);
            ret = -1;
        }
        else
        {
            for (t = 0; t < res->nt; t++)
            {
                const double *q = res->par + 6 * t;
                fprintf(fp, "%g  %g  %g  %g  %g  %g  \n",
                        q[0], q[1], q[2], q[3], q[4], q[5]);
            }
            if (fclose(fp) != 0)
                ret = -1;
        }
    }
    if (ret == 0 && rmsabs)
    {
        snprintf(name, len, "%s_abs.rms", base);
        ret = write_values_file(name, res->rms_abs, res->nt);
        snprintf(name, len, "%s_abs_mean.rms", base);
        if (ret == 0)
            ret = write_values_file(name, &res->rms_abs_mean, 1);
    }
    if (ret == 0 && rmsrel)
    {
        snprintf(name, len, "%s_rel.rms", base);
        ret = write_values_file(name, res->rms_rel, res->nt - 1);
        snprintf(name, len, "%s_rel_mean.rms", base);
        if (ret == 0)
            ret = write_values_file(name, &res->rms_rel_mean, 1);
    }
    free(name);
    return ret;
}

/**
 * \brief File name of the realigned series for an output base name.
 *
 * \param base  (in)  output base name
 * \param image (out) file name
 * \param size  (in)  size of image in bytes
 */
void CAT_VolRealignImageName(const char *base, char *image, size_t size)
{
    size_t n = strlen(base);
    const char *type = getenv("FSLOUTPUTTYPE");

    if ((n > 4 && strcmp(base + n - 4, ".nii") == 0) ||
        (n > 7 && strcmp(base + n - 7, ".nii.gz") == 0))
        snprintf(image, size, "%s", base);
    else
        snprintf(image, size, "%s%s", base,
                 (type && strcmp(type, "NIFTI") == 0) ? ".nii" : ".nii.gz");
}

/**
 * \brief Release the arrays of a result and zero it.
 *
 * \param res (in/out) result of CAT_VolRealignEstimate()
 */
void CAT_VolRealignFree(CAT_VolRealignResult *res)
{
    free(res->ras);
    free(res->fsl);
    free(res->par);
    free(res->rms_abs);
    free(res->rms_rel);
    free(res->scale);
    memset(res, 0, sizeof(*res));
}
