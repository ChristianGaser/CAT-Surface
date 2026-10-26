/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * C ports of the CAT12/CAT26 MATLAB helpers that the safeguards of CAT's
 * geodesic-shooting registration (run_Shooting in cat_main_registration.m)
 * call: cat_vol_smooth3X, cat_vol_div, cat_vol_morph ('l', 'd', 'o', 'ldc',
 * 'ldo'), cat_vol_localstat (mean), cat_vol_laplace3R and cat_vbdist
 * (Christian Gaser and Robert Dahnke).
 *
 * The originals crop to a bounding box, pad, recurse and threshold in a
 * particular order, and what happens at the volume border depends on it, so
 * the ports follow them step by step rather than computing the textbook
 * operation.  Volumes are stored in MATLAB order (x fastest).
 *
 * Copyright Christian Gaser, University of Jena.
 */

#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "CAT_Cat12Utils.h"
#include "CAT_SafeAlloc.h"
#include "CAT_Vol.h"

/* The CAT12 mex files are built without fused multiply-add (checked against
   cat_vol_interp3f); the exact ports below must not contract either. */
#if defined(__clang__)
#pragma STDC FP_CONTRACT OFF
#elif defined(__GNUC__)
#pragma GCC optimize("fp-contract=off")
#endif

#define NVOX(d) ((size_t)(d)[0] * (size_t)(d)[1] * (size_t)(d)[2])

/* ------------------------------------------------------------------------ */
/* Small helpers                                                             */
/* ------------------------------------------------------------------------ */

/** \brief Voxel size, {1,1,1} for NULL (the default of the originals). */
static void
get_vx(const double *vx, double out[3])
{
    int k;

    for (k = 0; k < 3; k++)
        out[k] = vx ? vx[k] : 1.0;
}

/** \brief Non-zero becomes 1. */
static void
binarize(unsigned char *vol, size_t n)
{
    size_t i;

    for (i = 0; i < n; i++)
        vol[i] = vol[i] != 0;
}

/** \brief Logical not of a 0/1 volume. */
static void
invert(unsigned char *vol, size_t n)
{
    size_t i;

    for (i = 0; i < n; i++)
        vol[i] = !vol[i];
}

/**
 * \brief Neighbour test of the CAT12 mex files (cat_vbdist.c, cat_vol_gradient3.c).
 *
 * The originals add a linear offset dx + dy*nx + dz*nx*ny to the index of the
 * voxel (u, v, w) and use the result if it lies in the volume and its
 * subscripts (ind2sub) differ from (u, v, w) by at most one in every
 * direction.  For nx >= 3 and ny >= 3 that is a bounds check; in thinner
 * volumes the offset can wrap into the next row or slice and still pass.  The
 * subscripts of the offset index are found here by carrying, which gives the
 * same answer as ind2sub without a division.
 *
 * \param u  (in) x subscript of the voxel
 * \param v  (in) y subscript of the voxel
 * \param w  (in) z subscript of the voxel
 * \param dx (in) x part of the offset (-1, 0, 1)
 * \param dy (in) y part of the offset
 * \param dz (in) z part of the offset
 * \param d  (in) {nx, ny, nz}
 * \return 1 if the offset index is a neighbour, otherwise 0
 */
static int
nb_ok(int u, int v, int w, int dx, int dy, int dz, const int d[3])
{
    int nu = u + dx, nv = v + dy, nw = w + dz;

    while (nu < 0)
    {
        nu += d[0];
        nv--;
    }
    while (nu >= d[0])
    {
        nu -= d[0];
        nv++;
    }
    while (nv < 0)
    {
        nv += d[1];
        nw--;
    }
    while (nv >= d[1])
    {
        nv -= d[1];
        nw++;
    }
    /* nw outside [0, nz) is an index outside [0, nvox) */
    return nw >= 0 && nw < d[2] && abs(nu - u) <= 1 && abs(nv - v) <= 1 &&
           abs(nw - w) <= 1;
}

/* ------------------------------------------------------------------------ */
/* Crop and pad (cat_vol_resize 'reduceBrain' / 'dereduceBrain')             */
/* ------------------------------------------------------------------------ */

/** Box of a crop, 0-based and inclusive. */
typedef struct
{
    int lo[3];    /**< first voxel of the box */
    int hi[3];    /**< last voxel of the box */
    int cdims[3]; /**< size of the box */
} crop_box;

/**
 * \brief Box of cat_vol_resize(Y, 'reduceBrain', vx, d, M).
 *
 * The bounding box of the non-zero voxels of M, enlarged by
 * max(1, round(d/vx)) voxels on each side and clipped to the volume.  An empty
 * M gives the whole volume.
 *
 * \param M    (in)  mask
 * \param dims (in)  {nx, ny, nz}
 * \param vx   (in)  voxel size
 * \param d    (in)  border in mm
 * \param bb   (out) box
 */
static void
reduce_brain(const unsigned char *M, const int dims[3], const double vx[3], double d,
             crop_box *bb)
{
    int first[3], last[3], dd[3], k, x, y, z, found = 0;
    size_t i = 0;

    for (k = 0; k < 3; k++)
    {
        double r = round(d / vx[k]);
        dd[k] = (r > 1.0) ? (int)r : 1;
        first[k] = INT_MAX;
        last[k] = -1;
    }
    for (z = 0; z < dims[2]; z++)
        for (y = 0; y < dims[1]; y++)
            for (x = 0; x < dims[0]; x++, i++)
                if (M[i])
                {
                    found = 1;
                    if (x < first[0]) first[0] = x;
                    if (x > last[0]) last[0] = x;
                    if (y < first[1]) first[1] = y;
                    if (y > last[1]) last[1] = y;
                    if (z < first[2]) first[2] = z;
                    if (z > last[2]) last[2] = z;
                }
    for (k = 0; k < 3; k++)
    {
        if (found)
        {
            bb->lo[k] = (first[k] - dd[k] > 0) ? first[k] - dd[k] : 0;
            bb->hi[k] = (last[k] + dd[k] < dims[k] - 1) ? last[k] + dd[k] : dims[k] - 1;
        }
        else
        {
            bb->lo[k] = 0;
            bb->hi[k] = dims[k] - 1;
        }
        bb->cdims[k] = bb->hi[k] - bb->lo[k] + 1;
    }
}

/** \brief Copy of the box of a binary volume. */
static unsigned char *
crop_u8(const unsigned char *vol, const int dims[3], const crop_box *bb)
{
    unsigned char *out = SAFE_MALLOC(unsigned char, NVOX(bb->cdims));
    size_t o = 0;
    int y, z;

    for (z = bb->lo[2]; z <= bb->hi[2]; z++)
        for (y = bb->lo[1]; y <= bb->hi[1]; y++, o += bb->cdims[0])
            memcpy(out + o, vol + ((size_t)z * dims[1] + y) * dims[0] + bb->lo[0],
                   bb->cdims[0]);
    return out;
}

/** \brief Copy of the box of a float volume. */
static float *
crop_f32(const float *vol, const int dims[3], const crop_box *bb)
{
    float *out = SAFE_MALLOC(float, NVOX(bb->cdims));
    size_t o = 0;
    int y, z;

    for (z = bb->lo[2]; z <= bb->hi[2]; z++)
        for (y = bb->lo[1]; y <= bb->hi[1]; y++, o += bb->cdims[0])
            memcpy(out + o, vol + ((size_t)z * dims[1] + y) * dims[0] + bb->lo[0],
                   bb->cdims[0] * sizeof(float));
    return out;
}

/** \brief 'dereduceBrain' of a binary volume: the box into a volume of zeros. */
static void
uncrop_u8(const unsigned char *cvol, const int dims[3], const crop_box *bb,
          unsigned char *vol)
{
    size_t o = 0;
    int y, z;

    memset(vol, 0, NVOX(dims));
    for (z = bb->lo[2]; z <= bb->hi[2]; z++)
        for (y = bb->lo[1]; y <= bb->hi[1]; y++, o += bb->cdims[0])
            memcpy(vol + ((size_t)z * dims[1] + y) * dims[0] + bb->lo[0], cvol + o,
                   bb->cdims[0]);
}

/** \brief 'dereduceBrain' of a float volume: the box into a volume of zeros. */
static void
uncrop_f32(const float *cvol, const int dims[3], const crop_box *bb, float *vol)
{
    size_t o = 0, i, n = NVOX(dims);
    int y, z;

    for (i = 0; i < n; i++)
        vol[i] = 0.0f;
    for (z = bb->lo[2]; z <= bb->hi[2]; z++)
        for (y = bb->lo[1]; y <= bb->hi[1]; y++, o += bb->cdims[0])
            memcpy(vol + ((size_t)z * dims[1] + y) * dims[0] + bb->lo[0], cvol + o,
                   bb->cdims[0] * sizeof(float));
}

/** \brief A binary volume in the middle of a volume of zeros, nv voxels wider per side. */
static unsigned char *
pad_u8(const unsigned char *vol, const int dims[3], const int nv[3], int pdims[3])
{
    unsigned char *out;
    int k, y, z;

    for (k = 0; k < 3; k++)
        pdims[k] = dims[k] + 2 * nv[k];
    out = SAFE_CALLOC(unsigned char, NVOX(pdims));
    for (z = 0; z < dims[2]; z++)
        for (y = 0; y < dims[1]; y++)
            memcpy(out + ((size_t)(z + nv[2]) * pdims[1] + y + nv[1]) * pdims[0] + nv[0],
                   vol + ((size_t)z * dims[1] + y) * dims[0], dims[0]);
    return out;
}

/** \brief vol | the middle of a padded volume (inverse of pad_u8()). */
static void
unpad_or_u8(const unsigned char *pvol, const int pdims[3], const int nv[3],
            const int dims[3], unsigned char *vol)
{
    int x, y, z;

    for (z = 0; z < dims[2]; z++)
        for (y = 0; y < dims[1]; y++)
        {
            const unsigned char *p =
                pvol + ((size_t)(z + nv[2]) * pdims[1] + y + nv[1]) * pdims[0] + nv[0];
            unsigned char *q = vol + ((size_t)z * dims[1] + y) * dims[0];
            for (x = 0; x < dims[0]; x++)
                q[x] = q[x] || p[x];
        }
}

/** \brief The middle of a padded volume (inverse of pad_u8()). */
static void
unpad_u8(const unsigned char *pvol, const int pdims[3], const int nv[3],
         const int dims[3], unsigned char *vol)
{
    int y, z;

    for (z = 0; z < dims[2]; z++)
        for (y = 0; y < dims[1]; y++)
            memcpy(vol + ((size_t)z * dims[1] + y) * dims[0],
                   pvol + ((size_t)(z + nv[2]) * pdims[1] + y + nv[1]) * pdims[0] + nv[0],
                   dims[0]);
}

/* ------------------------------------------------------------------------ */
/* Distance transform (cat_vbdist)                                           */
/* ------------------------------------------------------------------------ */

/** \brief Distance from voxel (u, v, w) to voxel j, in double as cat_vbdist.c computes it. */
static float
vb_exact(int u, int v, int w, int j, int nx, int nxy, float s1, float s2, float s3)
{
    int nw = j / nxy, nr = j - nw * nxy, nv = nr / nx, nu = nr - nv * nx;
    double a = (double)(u - nu) * s1, b = (double)(v - nv) * s2, c = (double)(w - nw) * s3;

    return (float)sqrt(a * a + b * b + c * c);
}

/**
 * \brief Euclidean distance transform of CAT12's [D, I] = cat_vbdist(V, M, vx).
 *
 * Exact port of cat_vbdist.c: a forward and a backward sweep over a 13-voxel
 * half neighbourhood (plus the voxel itself), each voxel of M taking the
 * nearest object index of the neighbour that minimises the neighbour's
 * distance plus the step length, then the exact distance to that object
 * voxel.  Same step lengths, same mixed single/double arithmetic, same
 * first-minimum tie break.
 *
 * \param vol  (in)  float volume, object where >= 0.5
 * \param mask (in)  voxels to compute (non-zero); NULL is all voxels
 * \param dims (in)  {nx, ny, nz}
 * \param vx   (in)  voxel size; NULL is {1,1,1}
 * \param D    (out) distances; FLT_MAX where not computed
 * \param idx  (out) 0-based index of the nearest object voxel; NULL to skip
 */
void
cat12_vbdist(const float *vol, const unsigned char *mask, const int dims[3],
             const double vx[3], float *D, unsigned int *idx)
{
    const int nx = dims[0], nxy = dims[0] * dims[1], nL = nxy * dims[2];
    /* NI of cat_vbdist.c: 0, -1, -x+1, -x, -x-1, -xy+1, -xy, -xy-1, -xy+x+1,
       -xy+x, -xy+x-1, -xy-x+1, -xy-x, -xy-x-1 */
    static const int OX[14] = {0, -1, 1, 0, -1, 1, 0, -1, 1, 0, -1, 1, 0, -1};
    static const int OY[14] = {0, 0, -1, -1, -1, 0, 0, 0, 1, 1, 1, -1, -1, -1};
    static const int OZ[14] = {0, 0, 0, 0, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1};
    double S[3];
    float s1, s2, s3, s12, s13, s23, s123;
    int NI[14];
    float NDST[14], DN[14];
    unsigned int *I;
    int i, n, pass;

    get_vx(vx, S);
    s1 = (float)fabs(S[0]);
    s2 = (float)fabs(S[1]);
    s3 = (float)fabs(S[2]);
    /* the mixed single/double expressions of cat_vbdist.c */
    s12 = (float)sqrt((double)s1 * s1 + s2 * s2);
    s13 = (float)sqrt((double)s1 * s1 + s3 * s3);
    s23 = (float)sqrt((double)s2 * s2 + s3 * s3);
    s123 = (float)sqrt((double)s12 * s12 + s3 * s3);

    for (n = 0; n < 14; n++)
        NI[n] = OX[n] + OY[n] * nx + OZ[n] * nxy;
    NDST[0] = 0.0f;
    NDST[1] = s1;
    NDST[2] = s12;
    NDST[3] = s2;
    NDST[4] = s12;
    NDST[5] = s13;
    NDST[6] = s3;
    NDST[7] = s13;
    NDST[8] = s123;
    NDST[9] = s23;
    NDST[10] = s123;
    NDST[11] = s123;
    NDST[12] = s23;
    NDST[13] = s123;

    I = idx ? idx : SAFE_MALLOC(unsigned int, nL);
    for (i = 0; i < nL; i++)
    {
        /* round(V) >= 1 in the original (unsigned char) */
        D[i] = (vol[i] >= 0.5f) ? 0.0f : FLT_MAX;
        I[i] = (unsigned int)i;
    }

    /* forward sweep (neighbours before i), then backward (after i) */
    for (pass = 0; pass < 2; pass++)
    {
        const int sgn = pass ? -1 : 1;
        int ii, u = 0, v = 0, w = 0;

        if (pass)
        {
            u = dims[0] - 1;
            v = dims[1] - 1;
            w = dims[2] - 1;
        }
        for (ii = 0; ii < nL; ii++)
        {
            float mn = FLT_MAX;
            int mi = 0;

            /* subscripts (u, v, w) of i, stepped along with it */
            if (ii > 0 && !pass)
            {
                if (++u == dims[0])
                {
                    u = 0;
                    if (++v == dims[1])
                    {
                        v = 0;
                        w++;
                    }
                }
            }
            else if (ii > 0)
            {
                if (--u < 0)
                {
                    u = dims[0] - 1;
                    if (--v < 0)
                    {
                        v = dims[1] - 1;
                        w--;
                    }
                }
            }
            i = pass ? nL - 1 - ii : ii;
            if (!(D[i] > 0.0f) || (mask && !mask[i]))
                continue;
            if (u > 0 && v > 0 && w > 0 && u < dims[0] - 1 && v < dims[1] - 1 &&
                w < dims[2] - 1)
            {
                /* away from the border every neighbour exists */
                for (n = 0; n < 14; n++)
                    DN[n] = D[i + sgn * NI[n]] + NDST[n];
            }
            else
                for (n = 0; n < 14; n++)
                {
                    int ni = i + sgn * NI[n];
                    if (!nb_ok(u, v, w, sgn * OX[n], sgn * OY[n], sgn * OZ[n], dims))
                        ni = i;
                    DN[n] = D[ni] + NDST[n];
                }
            for (n = 0; n < 14; n++)
                if (DN[n] > 0.0f && mn > DN[n])
                {
                    mn = DN[n];
                    mi = n;
                }
            if (mi > 0)
            {
                I[i] = I[i + sgn * NI[mi]];
                D[i] = vb_exact(u, v, w, (int)I[i], nx, nxy, s1, s2, s3);
            }
        }
    }

    if (!idx)
        free(I);
}

/* ------------------------------------------------------------------------ */
/* Morphology (cat_vol_morph)                                                */
/* ------------------------------------------------------------------------ */

/**
 * \brief Binary dilation with a box of 2*r+1 voxels per axis, zeros outside.
 *
 * convn(vol, ones(2*r+1), 'same') > 0.  The isotropic case is libCAT's
 * morph_box(); otherwise one running-window pass per axis.
 */
static void
box_dilate(unsigned char *vol, const int dims[3], const int r[3])
{
    int k;

    if (r[0] == r[1] && r[1] == r[2])
    {
        int d[3] = {dims[0], dims[1], dims[2]};
        morph_box(vol, vol, d, r[0], MORPH_BOX_DILATE);
        return;
    }
    for (k = 0; k < 3; k++)
    {
        const int len = dims[k], rk = r[k];
        const size_t stride = (k == 0) ? 1 : (k == 1) ? (size_t)dims[0]
                                                      : (size_t)dims[0] * dims[1];
        const int o1 = (k == 0) ? 1 : 0, o2 = (k == 2) ? 1 : 2;
        unsigned char *line;
        int a, b, t;

        if (rk <= 0)
            continue;
        line = SAFE_MALLOC(unsigned char, len);
        for (b = 0; b < dims[o2]; b++)
            for (a = 0; a < dims[o1]; a++)
            {
                size_t base;
                int cnt = 0;

                if (k == 0)
                    base = ((size_t)b * dims[1] + a) * dims[0];
                else if (k == 1)
                    base = (size_t)b * dims[0] * dims[1] + a;
                else
                    base = (size_t)b * dims[0] + a;
                for (t = 0; t < len; t++)
                    line[t] = vol[base + t * stride];
                for (t = 0; t < rk && t < len; t++)
                    cnt += line[t];
                for (t = 0; t < len; t++)
                {
                    if (t + rk < len)
                        cnt += line[t + rk];
                    if (t - rk - 1 >= 0)
                        cnt -= line[t - rk - 1];
                    vol[base + t * stride] = cnt > 0;
                }
            }
        free(line);
    }
}

/**
 * \brief 'cdilate' of cat_vol_morph: box dilation inside the bounding box.
 *
 * \param vol  (in/out) 0/1 volume
 * \param dims (in)     {nx, ny, nz}
 * \param nn   (in)     size in mm, rounded for the early exit and the crop only
 * \param vx   (in)     voxel size
 */
static void
m_cdilate(unsigned char *vol, const int dims[3], double nn, const double vx[3])
{
    double n = round(nn);
    crop_box bb;
    unsigned char *cv;
    int r[3], k;

    if (!(n > 0.0))
        return;
    reduce_brain(vol, dims, vx, n + 1.0, &bb);
    cv = crop_u8(vol, dims, &bb);
    for (k = 0; k < 3; k++)
        r[k] = (int)round(nn / vx[k]);
    box_dilate(cv, bb.cdims, r);
    uncrop_u8(cv, dims, &bb, vol);
    free(cv);
}

/** \brief 'cerode' of cat_vol_morph: ~cdilate(~vol, round(nn)). */
static void
m_cerode(unsigned char *vol, const int dims[3], double nn, const double vx[3])
{
    double n = round(nn);

    if (!(n > 0.0))
        return;
    invert(vol, NVOX(dims));
    m_cdilate(vol, dims, n, vx);
    invert(vol, NVOX(dims));
}

/**
 * \brief 'cclose' of cat_vol_morph (its test == 2 branch).
 *
 * Crop to the bounding box plus one voxel, pad by ceil(nn/vx) voxels of zeros,
 * cdilate and cerode with round(nn), cut out the box again and paste it into a
 * volume of zeros.
 */
static void
m_cclose(unsigned char *vol, const int dims[3], double nn, const double vx[3])
{
    double n = round(nn);
    crop_box bb;
    unsigned char *cv, *pv;
    int nv[3], pd[3], k;

    if (!(n > 0.0))
        return;
    for (k = 0; k < 3; k++)
        nv[k] = (int)ceil(nn / vx[k]);
    reduce_brain(vol, dims, vx, 1.0, &bb);
    cv = crop_u8(vol, dims, &bb);
    pv = pad_u8(cv, bb.cdims, nv, pd);
    m_cdilate(pv, pd, n, vx);
    m_cerode(pv, pd, n, vx);
    unpad_u8(pv, pd, nv, bb.cdims, cv);
    uncrop_u8(cv, dims, &bb, vol);
    free(pv);
    free(cv);
}

/** \brief 'copen' of cat_vol_morph: ~cclose(~vol, round(nn)). */
static void
m_copen(unsigned char *vol, const int dims[3], double nn, const double vx[3])
{
    double n = round(nn);

    if (!(n > 0.0))
        return;
    invert(vol, NVOX(dims));
    m_cclose(vol, dims, n, vx);
    invert(vol, NVOX(dims));
}

/** Tap of the distance kernel: offset into the volume and weight. */
typedef struct
{
    int dx, dy, dz;
    float w;
    int order;
} dist_tap;

/** \brief Sort taps by decreasing weight, ties in kernel order. */
static int
cmp_tap(const void *a, const void *b)
{
    const dist_tap *p = (const dist_tap *)a, *q = (const dist_tap *)b;

    if (p->w > q->w)
        return -1;
    if (p->w < q->w)
        return 1;
    return p->order - q->order;
}

/**
 * \brief 'distdilate' with round(nn) <= 5: thresholded convolution with a distance kernel.
 *
 * The kernel of cat_vol_morph: d = max(0, cat_vbdist(centre, vx) - 0.5) with
 * d(1) = d(end), then max(0, nn - d), each step rounded to single as MATLAB
 * does for single - double.  The result is
 * min(1, convn(single(vol), kernel, 'same')) >= 0.5; the weights are added in
 * double from the largest down and the sum stops once it reaches 0.5, which
 * gives the same answer as MATLAB's single sum unless that sum lies within a
 * float rounding of 0.5 (not the case for the kernels run_Shooting uses).
 *
 * \param vol  (in/out) 0/1 volume (the crop)
 * \param dims (in)     {nx, ny, nz}
 * \param n    (in)     round(nn), half size of the kernel
 * \param nn   (in)     distance in mm
 * \param vx   (in)     voxel size
 */
static void
conv_threshold(unsigned char *vol, const int dims[3], int n, double nn, const double vx[3])
{
    const int ks = 2 * n + 1, kd[3] = {ks, ks, ks};
    const size_t nk = (size_t)ks * ks * ks, nv = NVOX(dims);
    float *K = SAFE_CALLOC(float, nk), *D = SAFE_MALLOC(float, nk);
    dist_tap *taps = SAFE_MALLOC(dist_tap, nk);
    unsigned char *cand, *out;
    int ntap = 0, a, b, c, x, y, z, t, r[3] = {n, n, n};
    size_t i;

    K[((size_t)n * ks + n) * ks + n] = 1.0f;
    cat12_vbdist(K, NULL, kd, vx, D, NULL);
    for (i = 0; i < nk; i++)
    {
        float d = (float)((double)D[i] - 0.5);
        K[i] = (d > 0.0f) ? d : 0.0f;
    }
    K[0] = K[nk - 1];
    for (i = 0; i < nk; i++)
    {
        float w = (float)(nn - (double)K[i]);
        K[i] = (w > 0.0f) ? w : 0.0f;
    }

    /* convn 'same' with an odd kernel: out(p) = sum_q vol(p + n - q) K(q) */
    for (c = 0, i = 0; c < ks; c++)
        for (b = 0; b < ks; b++)
            for (a = 0; a < ks; a++, i++)
                if (K[i] > 0.0f)
                {
                    taps[ntap].dx = n - a;
                    taps[ntap].dy = n - b;
                    taps[ntap].dz = n - c;
                    taps[ntap].w = K[i];
                    taps[ntap].order = (int)i;
                    ntap++;
                }
    qsort(taps, ntap, sizeof(dist_tap), cmp_tap);

    /* only voxels within the box of the kernel can reach 0.5 */
    cand = SAFE_MALLOC(unsigned char, nv);
    memcpy(cand, vol, nv);
    box_dilate(cand, dims, r);
    out = SAFE_CALLOC(unsigned char, nv);
    for (z = 0, i = 0; z < dims[2]; z++)
        for (y = 0; y < dims[1]; y++)
            for (x = 0; x < dims[0]; x++, i++)
            {
                double s = 0.0;

                if (!cand[i])
                    continue;
                for (t = 0; t < ntap; t++)
                {
                    const int X = x + taps[t].dx, Y = y + taps[t].dy, Z = z + taps[t].dz;
                    if (X < 0 || Y < 0 || Z < 0 || X >= dims[0] || Y >= dims[1] ||
                        Z >= dims[2])
                        continue;
                    if (vol[((size_t)Z * dims[1] + Y) * dims[0] + X])
                    {
                        s += taps[t].w;
                        if (s >= 0.5)
                            break;
                    }
                }
                out[i] = (float)s >= 0.5f;
            }
    memcpy(vol, out, nv);

    free(out);
    free(cand);
    free(taps);
    free(D);
    free(K);
}

/**
 * \brief 'distdilate' of cat_vol_morph.
 *
 * Crop to the bounding box plus round((round(nn)+1)/vx) voxels; for
 * round(nn) > 5 keep cat_vbdist(vol) <= nn, otherwise conv_threshold().
 */
static void
m_distdilate(unsigned char *vol, const int dims[3], double nn, const double vx[3])
{
    double n = round(nn);
    crop_box bb;
    unsigned char *cv;
    size_t nc, i;

    if (!(n > 0.0))
        return;
    reduce_brain(vol, dims, vx, n + 1.0, &bb);
    cv = crop_u8(vol, dims, &bb);
    nc = NVOX(bb.cdims);
    if (n > 5.0)
    {
        float *f = SAFE_MALLOC(float, nc), *D = SAFE_MALLOC(float, nc);
        for (i = 0; i < nc; i++)
            f[i] = (float)cv[i];
        cat12_vbdist(f, NULL, bb.cdims, vx, D, NULL);
        for (i = 0; i < nc; i++)
            cv[i] = (double)D[i] <= nn;
        free(D);
        free(f);
    }
    else
        conv_threshold(cv, bb.cdims, (int)n, nn, vx);
    uncrop_u8(cv, dims, &bb, vol);
    free(cv);
}

/** \brief 'disterode' of cat_vol_morph: ~distdilate(~vol, nn). */
static void
m_disterode(unsigned char *vol, const int dims[3], double nn, const double vx[3])
{
    if (!(round(nn) > 0.0))
        return;
    invert(vol, NVOX(dims));
    m_distdilate(vol, dims, nn, vx);
    invert(vol, NVOX(dims));
}

/**
 * \brief 'distclose' of cat_vol_morph.
 *
 * Crop to the bounding box plus round((round(nn)+1)/vx) voxels, pad by
 * ceil(nn/vx) voxels of zeros, close (two cat_vbdist calls for round(nn) > 5,
 * otherwise distdilate and disterode) and OR the middle with the input.
 */
static void
m_distclose(unsigned char *vol, const int dims[3], double nn, const double vx[3])
{
    double n = round(nn);
    crop_box bb;
    unsigned char *cv, *pv;
    int nv[3], pd[3], k;
    size_t np, i;

    if (!(n > 0.0))
        return;
    for (k = 0; k < 3; k++)
        nv[k] = (int)ceil(nn / vx[k]);
    reduce_brain(vol, dims, vx, n + 1.0, &bb);
    cv = crop_u8(vol, dims, &bb);
    pv = pad_u8(cv, bb.cdims, nv, pd);
    np = NVOX(pd);
    if (n > 5.0)
    {
        float *f = SAFE_MALLOC(float, np), *D = SAFE_MALLOC(float, np);
        unsigned char *msk = SAFE_MALLOC(unsigned char, np);

        /* background beyond the dilation ... */
        for (i = 0; i < np; i++)
            f[i] = (float)pv[i];
        cat12_vbdist(f, NULL, pd, vx, D, NULL);
        for (i = 0; i < np; i++)
        {
            const int far = (double)D[i] > nn;
            f[i] = (float)far;
            msk[i] = !far;
        }
        /* ... and the dilated region that lies at least nn from it */
        cat12_vbdist(f, msk, pd, vx, D, NULL);
        for (i = 0; i < np; i++)
            pv[i] = (double)D[i] >= nn;
        free(msk);
        free(D);
        free(f);
    }
    else
    {
        m_distdilate(pv, pd, nn, vx);
        m_disterode(pv, pd, nn, vx);
    }
    unpad_or_u8(pv, pd, nv, bb.cdims, cv);
    uncrop_u8(cv, dims, &bb, vol);
    free(pv);
    free(cv);
}

/** \brief 'distopen' of cat_vol_morph: ~distclose(~vol, nn). */
static void
m_distopen(unsigned char *vol, const int dims[3], double nn, const double vx[3])
{
    if (!(round(nn) > 0.0))
        return;
    invert(vol, NVOX(dims));
    m_distclose(vol, dims, nn, vx);
    invert(vol, NVOX(dims));
}

/** Component of a labelling: size and rank of its first voxel. */
typedef struct
{
    size_t size;
    int label;
} lab_comp;

/** \brief Sort components by decreasing size, ties by label (stable sort). */
static int
cmp_comp(const void *a, const void *b)
{
    const lab_comp *p = (const lab_comp *)a, *q = (const lab_comp *)b;

    if (p->size > q->size)
        return -1;
    if (p->size < q->size)
        return 1;
    return p->label - q->label;
}

/**
 * \brief 'lab' of cat_vol_morph with n = [nmax th].
 *
 * spm_bwlabel(vol, 6) numbers the 6-connected components in the order of
 * their first voxel; the components are sorted by size (stably) and the first
 * is kept, plus, for nmax > 1, components 2 .. min(lim, ncomp, nmax) that are
 * larger than th (relative for 0 < th < 1, absolute for th < 0 or th > 1),
 * where lim is the last component that passes that size test.
 */
static void
m_lab(unsigned char *vol, const int dims[3], double nmax, double th)
{
    const int nx = dims[0], ny = dims[1], nz = dims[2];
    const size_t nv = NVOX(dims);
    int *lab = SAFE_CALLOC(int, nv);
    size_t *queue = SAFE_MALLOC(size_t, nv);
    unsigned char *keep;
    lab_comp *comp = NULL;
    int ncomp = 0, cap = 0, r;
    size_t i;

    /* breadth-first labelling in memory order */
    for (i = 0; i < nv; i++)
    {
        size_t head = 0, tail = 0;

        if (!vol[i] || lab[i])
            continue;
        if (ncomp == cap)
        {
            cap = cap ? 2 * cap : 64;
            comp = (lab_comp *)realloc(comp, cap * sizeof(lab_comp));
            if (!comp)
            {
                fprintf(stderr, "Memory allocation error in m_lab\n");
                exit(EXIT_FAILURE);
            }
        }
        ncomp++;
        lab[i] = ncomp;
        queue[tail++] = i;
        while (head < tail)
        {
            const size_t p = queue[head++];
            const int x = (int)(p % nx), y = (int)((p / nx) % ny), z = (int)(p / ((size_t)nx * ny));
            size_t q[6];
            int nq = 0, j;

            if (x > 0) q[nq++] = p - 1;
            if (x < nx - 1) q[nq++] = p + 1;
            if (y > 0) q[nq++] = p - nx;
            if (y < ny - 1) q[nq++] = p + nx;
            if (z > 0) q[nq++] = p - (size_t)nx * ny;
            if (z < nz - 1) q[nq++] = p + (size_t)nx * ny;
            for (j = 0; j < nq; j++)
                if (vol[q[j]] && !lab[q[j]])
                {
                    lab[q[j]] = ncomp;
                    queue[tail++] = q[j];
                }
        }
        comp[ncomp - 1].size = tail;
        comp[ncomp - 1].label = ncomp;
    }
    free(queue);
    if (ncomp == 0)
    {
        free(lab);
        return;
    }

    qsort(comp, ncomp, sizeof(lab_comp), cmp_comp);
    keep = SAFE_CALLOC(unsigned char, ncomp + 1);
    keep[comp[0].label] = 1;
    if (nmax > 1.0)
    {
        const int rel = th > 0.0 && th < 1.0, abso = th < 0.0 || th > 1.0;
        double snum = 0.0, upper;
        int lim = 0;

        for (r = 0; r < ncomp; r++)
            snum += (double)comp[r].size;
        for (r = 0; r < ncomp; r++)
            if (rel ? ((double)comp[r].size / snum > th) : (fabs((double)comp[r].size) > th))
                lim = r + 1;
        if (lim > 0)
        {
            upper = (double)lim;
            if ((double)ncomp < upper)
                upper = (double)ncomp;
            if (nmax < upper)
                upper = nmax;
            for (r = 2; (double)r <= upper; r++)
            {
                const double sz = (double)comp[r - 1].size;
                if ((rel && sz / snum > th) || (abso && sz > fabs(th)))
                    keep[comp[r - 1].label] = 1;
            }
        }
    }
    for (i = 0; i < nv; i++)
        vol[i] = vol[i] && keep[lab[i]];
    free(keep);
    free(comp);
    free(lab);
}

/**
 * \brief Largest component of CAT12's cat_vol_morph(vol, 'l').
 *
 * 6-connected labelling in memory order, ties go to the first component.
 *
 * \param vol  (in/out) binary volume; non-zero is foreground, written as 0/1
 * \param dims (in)     {nx, ny, nz}
 */
void
cat12_vol_morph_lab(unsigned char *vol, const int dims[3])
{
    binarize(vol, NVOX(dims));
    m_lab(vol, dims, 1.0, 0.0);
}

/**
 * \brief Several large components of CAT12's cat_vol_morph(vol, 'l', [nmax th]).
 *
 * See m_lab(); returns the union of the kept components instead of their ranks.
 *
 * \param vol  (in/out) binary volume; non-zero is foreground, written as 0/1
 * \param dims (in)     {nx, ny, nz}
 * \param nmax (in)     largest number of components to keep
 * \param th   (in)     relative (0 < th < 1) or absolute (th < 0, th > 1) size limit
 */
void
cat12_vol_morph_lab_multi(unsigned char *vol, const int dims[3], double nmax, double th)
{
    binarize(vol, NVOX(dims));
    m_lab(vol, dims, nmax, th);
}

/**
 * \brief Dilation of CAT12's cat_vol_morph(vol, 'd', n, vx): m_cdilate().
 *
 * \param vol  (in/out) binary volume; non-zero is foreground, written as 0/1
 * \param dims (in)     {nx, ny, nz}
 * \param n    (in)     size in mm
 * \param vx   (in)     voxel size in mm; NULL is {1,1,1}
 */
void
cat12_vol_morph_dilate(unsigned char *vol, const int dims[3], double n, const double vx[3])
{
    double v[3];

    get_vx(vx, v);
    binarize(vol, NVOX(dims));
    m_cdilate(vol, dims, n, v);
}

/**
 * \brief Opening of CAT12's cat_vol_morph(vol, 'o', n, vx): m_copen().
 *
 * \param vol  (in/out) binary volume; non-zero is foreground, written as 0/1
 * \param dims (in)     {nx, ny, nz}
 * \param n    (in)     size in mm
 * \param vx   (in)     voxel size in mm; NULL is {1,1,1}
 */
void
cat12_vol_morph_open(unsigned char *vol, const int dims[3], double n, const double vx[3])
{
    double v[3];

    get_vx(vx, v);
    binarize(vol, NVOX(dims));
    m_copen(vol, dims, n, v);
}

/**
 * \brief Closing and hole filling of CAT12's cat_vol_morph(vol, 'ldc', n, vx).
 *
 * m_distclose(), then ~lab(~vol).  The 'lab' inside gets n as its n(1), which
 * for n > 1 makes MATLAB return single, but with n(2) = 0 it still keeps only
 * the largest component.
 *
 * \param vol  (in/out) binary volume; non-zero is foreground, written as 0/1
 * \param dims (in)     {nx, ny, nz}
 * \param n    (in)     distance in mm
 * \param vx   (in)     voxel size in mm; NULL is {1,1,1}
 */
void
cat12_vol_morph_labdistclose(unsigned char *vol, const int dims[3], double n,
                             const double vx[3])
{
    double v[3];

    get_vx(vx, v);
    binarize(vol, NVOX(dims));
    m_distclose(vol, dims, n, v);
    /* ~lab(~vol): fill all background components but the largest */
    invert(vol, NVOX(dims));
    m_lab(vol, dims, n, 0.0);
    invert(vol, NVOX(dims));
}

/**
 * \brief Opening and largest component of CAT12's cat_vol_morph(vol, 'ldo', n, vx).
 *
 * m_distopen(), then the largest 6-connected component.
 *
 * \param vol  (in/out) binary volume; non-zero is foreground, written as 0/1
 * \param dims (in)     {nx, ny, nz}
 * \param n    (in)     distance in mm
 * \param vx   (in)     voxel size in mm; NULL is {1,1,1}
 */
void
cat12_vol_morph_labdistopen(unsigned char *vol, const int dims[3], double n,
                            const double vx[3])
{
    double v[3];

    get_vx(vx, v);
    binarize(vol, NVOX(dims));
    m_distopen(vol, dims, n, v);
    m_lab(vol, dims, n, 0.0);
}

/* ------------------------------------------------------------------------ */
/* Smoothing (cat_vol_smooth3X)                                              */
/* ------------------------------------------------------------------------ */

/** \brief NaN and +-Inf become 0. */
static void
zero_nonfinite(float *vol, size_t n)
{
    size_t i;

    for (i = 0; i < n; i++)
        if (isnan(vol[i]) || isinf(vol[i]))
            vol[i] = 0.0f;
}

/**
 * \brief MATLAB's smooth3(S, 'gaussian', sz, sd).
 *
 * Kernel exp(-(x^2+y^2+z^2)/(2 sd^2)) normalised to sum 1 and rounded to
 * single (convn of single data with a double kernel), replicate padding,
 * 'valid' convolution.  Accumulated in double.
 *
 * \param S    (in)  float volume
 * \param dims (in)  {nx, ny, nz}
 * \param sz   (in)  odd kernel size (1, 3 or 5)
 * \param sd   (in)  standard deviation in voxels
 * \param out  (out) float volume, must not alias S
 */
static void
smooth3_gauss(const float *S, const int dims[3], int sz, double sd, float *out)
{
    const int r = (sz - 1) / 2, nk = sz * sz * sz;
    double *h = SAFE_MALLOC(double, nk), sum = 0.0;
    float *hk = SAFE_MALLOC(float, nk);
    int *cx = SAFE_MALLOC(int, (size_t)dims[0] * sz);
    int *cy = SAFE_MALLOC(int, (size_t)dims[1] * sz);
    int *cz = SAFE_MALLOC(int, (size_t)dims[2] * sz);
    int a, b, c, x, y, z, k, i;
    size_t o = 0;

    /* gaussian3 of smooth3.m: x runs along the columns (second dimension) */
    for (c = 0, i = 0; c < sz; c++)
        for (b = 0; b < sz; b++)
            for (a = 0; a < sz; a++, i++)
            {
                const double X = b - r, Y = a - r, Z = c - r;
                h[i] = exp(-(X * X + Y * Y + Z * Z) / (2.0 * sd * sd));
            }
    for (i = 0; i < nk; i++)
        sum += h[i];
    for (i = 0; i < nk; i++)
        hk[i] = (float)(h[i] / sum);

    /* replicate padding: clamped source coordinate of every output and tap */
    for (k = 0; k < sz; k++)
    {
        for (x = 0; x < dims[0]; x++)
        {
            int s = x + r - k;
            cx[x * sz + k] = s < 0 ? 0 : (s >= dims[0] ? dims[0] - 1 : s);
        }
        for (y = 0; y < dims[1]; y++)
        {
            int s = y + r - k;
            cy[y * sz + k] = s < 0 ? 0 : (s >= dims[1] ? dims[1] - 1 : s);
        }
        for (z = 0; z < dims[2]; z++)
        {
            int s = z + r - k;
            cz[z * sz + k] = s < 0 ? 0 : (s >= dims[2] ? dims[2] - 1 : s);
        }
    }

    for (z = 0; z < dims[2]; z++)
        for (y = 0; y < dims[1]; y++)
            for (x = 0; x < dims[0]; x++, o++)
            {
                double acc = 0.0;
                for (c = 0, i = 0; c < sz; c++)
                {
                    const size_t zo = (size_t)cz[z * sz + c] * dims[1];
                    for (b = 0; b < sz; b++)
                    {
                        const float *row = S + (zo + cy[y * sz + b]) * dims[0];
                        for (a = 0; a < sz; a++, i++)
                            acc += (double)hk[i] * (double)row[cx[x * sz + a]];
                    }
                }
                out[o] = (float)acc;
            }

    free(cz);
    free(cy);
    free(cx);
    free(hk);
    free(h);
}

/**
 * \brief cat_vol_interp3f(F, X, Y, Z, 'linear') at one point.
 *
 * Exact copy of the arithmetic of the mex (Brian Amberg's ba_interp3): X runs
 * along the second dimension and Y along the first (meshgrid convention),
 * coordinates are 1-based, the eight indices are clamped to the volume.
 *
 * \param F  (in) float volume
 * \param fd (in) its dimensions
 * \param X  (in) coordinate along the second dimension
 * \param Y  (in) coordinate along the first dimension
 * \param Z  (in) coordinate along the third dimension
 * \return interpolated value
 */
static float
interp3f_linear(const float *F, const int fd[3], float X, float Y, float Z)
{
    const int M = fd[0], N = fd[1], O = fd[2];
    const float xf = floorf(X), yf = floorf(Y), zf = floorf(Z);
    const float dx = X - xf, dy = Y - yf, dz = Z - zf;
    const float wx0 = (float)(1.0 - dx), wx1 = dx;
    const float wy0 = (float)(1.0 - dy), wy1 = dy;
    const float wz0 = (float)(1.0 - dz), wz1 = dz;
    const int xi = (int)(xf - 1), yi = (int)(yf - 1), zi = (int)(zf - 1);
    const int x0 = xi < 0 ? 0 : (xi >= N ? N - 1 : xi);
    const int x1 = xi + 1 < 0 ? 0 : (xi + 1 >= N ? N - 1 : xi + 1);
    const int y0 = yi < 0 ? 0 : (yi >= M ? M - 1 : yi);
    const int y1 = yi + 1 < 0 ? 0 : (yi + 1 >= M ? M - 1 : yi + 1);
    const int z0 = zi < 0 ? 0 : (zi >= O ? O - 1 : zi);
    const int z1 = zi + 1 < 0 ? 0 : (zi + 1 >= O ? O - 1 : zi + 1);
#define F3(yy, xx, zz) F[(yy) + (size_t)M * ((xx) + (size_t)N * (zz))]
    return wz0 * (wy0 * (wx0 * F3(y0, x0, z0) + wx1 * F3(y0, x1, z0)) +
                  wy1 * (wx0 * F3(y1, x0, z0) + wx1 * F3(y1, x1, z0))) +
           wz1 * (wy0 * (wx0 * F3(y0, x0, z1) + wx1 * F3(y0, x1, z1)) +
                  wy1 * (wx0 * F3(y1, x0, z1) + wx1 * F3(y1, x1, z1)));
#undef F3
}

/**
 * \brief reduceRes of cat_vol_smooth3X: half resolution.
 *
 * Odd sizes are padded by repeating the last plane, then cat_vol_interp3f
 * samples the volume at 1.5:2:end (the mean of each 2x2x2 block).
 *
 * \param S    (in)  float volume
 * \param dims (in)  its dimensions
 * \param rd   (out) dimensions of the result, ceil(dims/2)
 * \return newly allocated volume of size rd
 */
static float *
reduce_res(const float *S, const int dims[3], int rd[3])
{
    int pd[3], k, x, y, z;
    float *P, *R;
    size_t o = 0;

    for (k = 0; k < 3; k++)
    {
        pd[k] = dims[k] + (dims[k] % 2);
        rd[k] = pd[k] / 2;
    }
    P = SAFE_MALLOC(float, NVOX(pd));
    for (z = 0; z < pd[2]; z++)
        for (y = 0; y < pd[1]; y++)
            for (x = 0; x < pd[0]; x++, o++)
            {
                const int sx = x < dims[0] ? x : dims[0] - 1;
                const int sy = y < dims[1] ? y : dims[1] - 1;
                const int sz = z < dims[2] ? z : dims[2] - 1;
                P[o] = S[((size_t)sz * dims[1] + sy) * dims[0] + sx];
            }
    R = SAFE_MALLOC(float, NVOX(rd));
    o = 0;
    for (z = 0; z < rd[2]; z++)
        for (y = 0; y < rd[1]; y++)
            for (x = 0; x < rd[0]; x++, o++)
                R[o] = interp3f_linear(P, pd, 1.5f + 2.0f * y, 1.5f + 2.0f * x,
                                       1.5f + 2.0f * z);
    free(P);
    return R;
}

/**
 * \brief dereduceRes of cat_vol_smooth3X: back to the original size.
 *
 * cat_vol_interp3f at 0.75:0.5:dims/2+0.25 of the reduced volume, indices
 * clamped at the border.
 *
 * \param R    (in)  reduced volume
 * \param rd   (in)  its dimensions
 * \param dims (in)  original dimensions
 * \param out  (out) volume of size dims
 */
static void
dereduce_res(const float *R, const int rd[3], const int dims[3], float *out)
{
    int x, y, z;
    size_t o = 0;

    for (z = 0; z < dims[2]; z++)
        for (y = 0; y < dims[1]; y++)
            for (x = 0; x < dims[0]; x++, o++)
                out[o] = interp3f_linear(R, rd, 0.75f + 0.5f * y, 0.75f + 0.5f * x,
                                         0.75f + 0.5f * z);
}

/**
 * \brief Gaussian smoothing of CAT12's cat_vol_smooth3X(S, s).
 *
 * The branches of the original: blend with smooth3(S,'gaussian',3,0.5) for
 * s < 0.5, smooth3(S,'gaussian',3,s) up to s = 1, half resolution
 * (reduce_res()), s/2 and back (dereduce_res()) above, the 5x5x5 kernel once
 * a dimension is <= 6.  NaN and Inf are zeroed before and after.
 *
 * \param vol  (in/out) float volume, smoothed in place
 * \param dims (in)     {nx, ny, nz}
 * \param s    (in)     filter size in voxels
 */
void
cat12_vol_smooth3X(float *vol, const int dims[3], double s)
{
    const size_t n = NVOX(dims);
    size_t i;

    zero_nonfinite(vol, n);
    if (s > 0.0 && s < 0.5)
    {
        float *t = SAFE_MALLOC(float, n);
        smooth3_gauss(vol, dims, 3, 0.5, t);
        /* single * double is computed in double and rounded to single */
        for (i = 0; i < n; i++)
            vol[i] = (float)((double)t[i] * s) + (float)((double)vol[i] * (1.0 - s));
        free(t);
    }
    else if (s >= 0.5 && s <= 1.0)
    {
        float *t = SAFE_MALLOC(float, n);
        smooth3_gauss(vol, dims, 3, s, t);
        memcpy(vol, t, n * sizeof(float));
        free(t);
    }
    else if (s > 1.0 && dims[0] > 6 && dims[1] > 6 && dims[2] > 6)
    {
        int rd[3];
        float *R = reduce_res(vol, dims, rd);
        cat12_vol_smooth3X(R, rd, s / 2.0);
        dereduce_res(R, rd, dims, vol);
        free(R);
    }
    else if (s >= 1.0)
    {
        float *t = SAFE_MALLOC(float, n);
        smooth3_gauss(vol, dims, 5, s, t);
        memcpy(vol, t, n * sizeof(float));
        free(t);
    }
    else if (isnan(s))
    {
        /* smooth3(S,'gaussian',1,NaN) has the kernel NaN */
        for (i = 0; i < n; i++)
            vol[i] = 0.0f;
    }
    /* s == 0 and s < 0 (kernel 1) leave the volume unchanged */
    zero_nonfinite(vol, n);
}

/* ------------------------------------------------------------------------ */
/* Divergence (cat_vol_div)                                                  */
/* ------------------------------------------------------------------------ */

/**
 * \brief One component of cat_vol_gradient3: (I(i+1) - I(i-1)) / 2 along an axis.
 *
 * A missing neighbour (see nb_ok()) is replaced by the voxel itself, which
 * gives half the one-sided difference at the border.
 *
 * \param I    (in)  float volume
 * \param dims (in)  {nx, ny, nz}
 * \param axis (in)  0, 1 or 2
 * \param G    (out) gradient component, must not alias I
 */
static void
grad_axis(const float *I, const int dims[3], int axis, float *G)
{
    const int step = (axis == 0) ? 1 : (axis == 1) ? dims[0] : dims[0] * dims[1];
    const int dx = axis == 0, dy = axis == 1, dz = axis == 2;
    int i = 0, u, v, w;

    for (w = 0; w < dims[2]; w++)
        for (v = 0; v < dims[1]; v++)
            for (u = 0; u < dims[0]; u++, i++)
            {
                int n1 = i - step, n2 = i + step;

                if (!nb_ok(u, v, w, -dx, -dy, -dz, dims))
                    n1 = i;
                if (!nb_ok(u, v, w, dx, dy, dz, dims))
                    n2 = i;
                G[i] = (float)((double)(I[n2] - I[n1]) / 2.0);
            }
}

/**
 * \brief Divergence of CAT12's cat_vol_div(Ym) with the defaults (vx_vol = vx_volr = 1.5).
 *
 * NaN filling by cat12_vbdist(), crop by reduce_brain(), max(1/3, Ym), two
 * passes of grad_axis() with the division by 1.5 in between, the cubic
 * 'dereduceV' at the grid points and the paste back into zeros, all in the
 * single-precision steps of the original.
 *
 * \param vol  (in)  float volume
 * \param dims (in)  {nx, ny, nz}
 * \param div  (out) float volume; may be the same array as vol
 */
void
cat12_vol_div(const float *vol, const int dims[3], float *div)
{
    static const double vx15[3] = {1.5, 1.5, 1.5};
    const size_t nL = NVOX(dims);
    const float third = (float)(1.0 / 3.0);
    float *Ym = SAFE_MALLOC(float, nL), *T, *g, *p, *Yd;
    unsigned char *Ynan = SAFE_MALLOC(unsigned char, nL), *M;
    crop_box bb;
    size_t i, nc;
    int any_nan = 0, any_inf = 0, ax;

    for (i = 0; i < nL; i++)
    {
        Ym[i] = vol[i];
        Ynan[i] = isnan(vol[i]) != 0;
        any_nan |= Ynan[i];
    }

    /* NaN voxels near finite ones take the nearest finite value:
       [D,I] = cat_vbdist(single(~Ynan), cat_vol_morph(~Ynan,'d',2)); Ym(D<2) = Ym(I(D<2)) */
    if (any_nan)
    {
        static const double one[3] = {1.0, 1.0, 1.0};
        unsigned char *msk = SAFE_MALLOC(unsigned char, nL);
        float *src = SAFE_MALLOC(float, nL), *D = SAFE_MALLOC(float, nL);
        unsigned int *I = SAFE_MALLOC(unsigned int, nL);

        for (i = 0; i < nL; i++)
        {
            msk[i] = !Ynan[i];
            src[i] = (float)msk[i];
        }
        m_cdilate(msk, dims, 2.0, one);
        cat12_vbdist(src, msk, dims, NULL, D, I);
        for (i = 0; i < nL; i++)
            if (D[i] < 2.0f)
                Ym[i] = vol[I[i]];
        free(I);
        free(D);
        free(src);
        free(msk);
    }

    /* crop to the bounding box of Ym ~= 0 & ~Ynan plus round(2/1.5) voxel;
       'reduceV' from 1.5 to 1.5 mm is a no-op */
    M = SAFE_MALLOC(unsigned char, nL);
    for (i = 0; i < nL; i++)
        M[i] = Ym[i] != 0.0f && !Ynan[i];
    reduce_brain(M, dims, vx15, 2.0, &bb);
    free(M);
    T = crop_f32(Ym, dims, &bb);
    free(Ym);
    nc = NVOX(bb.cdims);

    /* gradient of max(1/3, Ym) (NaN -> 1/3), divided by vx_volr = 1.5, and
       the derivative of each component along its own axis */
    for (i = 0; i < nc; i++)
        if (!((double)T[i] > 1.0 / 3.0))
            T[i] = third;
    g = SAFE_MALLOC(float, nc);
    p = SAFE_MALLOC(float, nc);
    Yd = SAFE_CALLOC(float, nc);
    for (ax = 0; ax < 3; ax++)
    {
        grad_axis(T, bb.cdims, ax, g);
        for (i = 0; i < nc; i++)
            g[i] = (float)((double)g[i] / 1.5);
        grad_axis(g, bb.cdims, ax, p);
        /* single(px) + single(py) + single(pz), left to right */
        for (i = 0; i < nc; i++)
            Yd[i] = (ax == 0) ? p[i] : Yd[i] + p[i];
    }
    free(p);
    free(g);
    free(T);

    /* 'dereduceV' with 'cubic' at the grid points: NaN -> 0, then the identity
       unless a +-Inf sits in the 4x4x4 window of cat_vol_interp3f */
    for (i = 0; i < nc; i++)
    {
        if (isnan(Yd[i]))
            Yd[i] = 0.0f;
        else if (isinf(Yd[i]))
            any_inf = 1;
    }
    if (any_inf)
    {
        const int *cd = bb.cdims;
        float *t = SAFE_MALLOC(float, nc);
        int x, y, z, a, b, c;

        memcpy(t, Yd, nc * sizeof(float));
        for (z = 0, i = 0; z < cd[2]; z++)
            for (y = 0; y < cd[1]; y++)
                for (x = 0; x < cd[0]; x++, i++)
                    for (c = -1; c <= 2; c++)
                        for (b = -1; b <= 2; b++)
                            for (a = -1; a <= 2; a++)
                            {
                                int X = x + a, Y = y + b, Z = z + c;
                                if (a == 0 && b == 0 && c == 0)
                                    continue;
                                X = X < 0 ? 0 : (X >= cd[0] ? cd[0] - 1 : X);
                                Y = Y < 0 ? 0 : (Y >= cd[1] ? cd[1] - 1 : Y);
                                Z = Z < 0 ? 0 : (Z >= cd[2] ? cd[2] - 1 : Z);
                                if (isinf(Yd[((size_t)Z * cd[1] + Y) * cd[0] + X]))
                                    t[i] = NAN;
                            }
        free(Yd);
        Yd = t;
    }

    /* 'dereduceBrain' (zeros outside the box) and Ydiv(Ynan) = 0 */
    uncrop_f32(Yd, dims, &bb, div);
    for (i = 0; i < nL; i++)
        if (Ynan[i])
            div[i] = 0.0f;
    free(Yd);
    free(Ynan);
}

/* ------------------------------------------------------------------------ */
/* Local mean (cat_vol_localstat) and Laplace filter (cat_vol_laplace3R)     */
/* ------------------------------------------------------------------------ */

/**
 * \brief Local mean of CAT12's cat_vol_localstat(V, M, nb, 1) (iter 1, mskval 0, filter0 0).
 *
 * The loops of cat_vol_localstat.c: the sphere of radius nb voxels inside the
 * volume, finite non-zero voxels of the mask only, a single-precision sum in
 * the order x offset, y offset, z offset, divided by the count.
 *
 * \param vol  (in)  float volume
 * \param mask (in)  mask; NULL is all voxels
 * \param dims (in)  {nx, ny, nz}
 * \param nb   (in)  radius in voxels, 1..10
 * \param out  (out) float volume, must not alias vol
 */
void
cat12_vol_localstat_mean(const float *vol, const unsigned char *mask,
                         const int dims[3], int nb, float *out)
{
    const size_t nL = NVOX(dims);
    const int nk = (2 * nb + 1) * (2 * nb + 1) * (2 * nb + 1);
    unsigned char *MB = SAFE_MALLOC(unsigned char, nL);
    float *NV = SAFE_MALLOC(float, nk);
    int x, y, z, i, j, k;
    size_t ind;

    /* mask without NaN, Inf and (filter0 = 0) zeros */
    for (ind = 0; ind < nL; ind++)
        MB[ind] = (!mask || mask[ind]) && !isnan(vol[ind]) && !isinf(vol[ind]) &&
                  vol[ind] != 0.0f;

    for (z = 0, ind = 0; z < dims[2]; z++)
        for (y = 0; y < dims[1]; y++)
            for (x = 0; x < dims[0]; x++, ind++)
            {
                float m = 0.0f, cnt = 0.0f;
                int n = 0, t;

                out[ind] = 0.0f;
                if (!MB[ind])
                    continue;
                /* neighbours within the sphere, in the loop order of the mex */
                for (i = -nb; i <= nb; i++)
                    for (j = -nb; j <= nb; j++)
                        for (k = -nb; k <= nb; k++)
                        {
                            size_t ni;
                            if (x + i < 0 || x + i >= dims[0] || y + j < 0 ||
                                y + j >= dims[1] || z + k < 0 || z + k >= dims[2])
                                continue;
                            ni = ((size_t)(z + k) * dims[1] + (y + j)) * dims[0] + (x + i);
                            if (MB[ni] && sqrtf((float)(i * i + j * j + k * k)) <= (float)nb)
                                NV[n++] = vol[ni];
                        }
                for (t = 0; t < n; t++)
                {
                    m += NV[t];
                    cnt++;
                }
                out[ind] = m / cnt;
            }
    free(NV);
    free(MB);
}

/**
 * \brief Laplace filter of CAT12's cat_vol_laplace3R(V, M, TH): libCAT's laplace3R().
 *
 * laplace3R() is the same algorithm; TH is rounded to single first, as the mex
 * reads it.
 *
 * \param vol  (in/out) float volume, filtered in place
 * \param mask (in)     filter region
 * \param dims (in)     {nx, ny, nz}
 * \param th   (in)     largest change of a sweep that ends the iteration
 */
void
cat12_vol_laplace3R(float *vol, const unsigned char *mask, const int dims[3], double th)
{
    int d[3] = {dims[0], dims[1], dims[2]};

    /* the mex reads TH as single */
    laplace3R(vol, (unsigned char *)mask, d, (double)(float)th);
}
