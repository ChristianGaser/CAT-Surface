/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h>

#include "CAT_Sheetness.h"
#include "CAT_Vol.h"
#include "CAT_Math.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* FWHM = SQRT8LN2 * sigma; smooth3() is parameterized by FWHM in mm. */
#define SQRT8LN2 2.3548200450309493

/**
 * \brief Fill a CAT_SheetnessOpts with defaults suitable for 0.5 mm brain data.
 *
 * The scale range brackets the structures this filter is meant to find: a
 * sulcal CSF sheet or a gyral WM blade that survives at all is between one and
 * three voxels thick at 0.5 mm.  Scales far above that would respond to the
 * cortical ribbon itself, which is not what we are looking for.
 *
 * \param opts (out) option block to initialize; NULL is ignored
 */
void CAT_SheetnessOptionsInit(CAT_SheetnessOpts *opts)
{
    if (!opts)
        return;
    opts->sigma_min = 0.3;
    opts->sigma_max = 1.0;
    opts->n_scales = 3;
    opts->alpha = 0.5;
    opts->beta = 0.5;
    opts->c = -1.0; /* automatic: half the maximum Hessian norm */
    opts->gain = 1.0;
    opts->normalize = CAT_SHEETNESS_NORMALIZE;
    opts->polarity = 0;
    opts->verbose = 0;
}

/**
 * \brief Eigen-decomposition of a symmetric 3x3 matrix, sorted by magnitude.
 *
 * Uses the closed-form trigonometric solution for the characteristic
 * polynomial of a symmetric 3x3 matrix.  The eigenvector of the dominant
 * eigenvalue is recovered from the null space of (A - lambda I) by taking the
 * cross product of the two most linearly independent rows; the pair with the
 * largest cross-product norm is the numerically safest choice and also handles
 * the rank-1 case, where two rows are parallel.
 *
 * \param a     (in)  upper triangle {axx, axy, axz, ayy, ayz, azz}
 * \param eval  (out) eigenvalues sorted by ascending absolute value
 * \param evec3 (out) unit eigenvector of eval[2], or NULL to skip
 */
void CAT_EigenSym3(const double a[6], double eval[3], double evec3[3])
{
    const double axx = a[0], axy = a[1], axz = a[2];
    const double ayy = a[3], ayz = a[4], azz = a[5];
    double p1, q, p2, p, r, phi, e[3], t;
    int i, j;

    p1 = axy * axy + axz * axz + ayz * ayz;
    q = (axx + ayy + azz) / 3.0;

    if (p1 <= 1e-30)
    {
        /* already diagonal */
        e[0] = axx;
        e[1] = ayy;
        e[2] = azz;
    }
    else
    {
        double b[6];
        p2 = (axx - q) * (axx - q) + (ayy - q) * (ayy - q) + (azz - q) * (azz - q) + 2.0 * p1;
        p = sqrt(p2 / 6.0);
        if (p < 1e-30)
            p = 1e-30;

        /* B = (A - q I) / p */
        b[0] = (axx - q) / p;
        b[1] = axy / p;
        b[2] = axz / p;
        b[3] = (ayy - q) / p;
        b[4] = ayz / p;
        b[5] = (azz - q) / p;

        /* r = det(B) / 2 */
        r = 0.5 * (b[0] * (b[3] * b[5] - b[4] * b[4]) - b[1] * (b[1] * b[5] - b[4] * b[2]) + b[2] * (b[1] * b[4] - b[3] * b[2]));

        if (r <= -1.0)
            phi = M_PI / 3.0;
        else if (r >= 1.0)
            phi = 0.0;
        else
            phi = acos(r) / 3.0;

        e[0] = q + 2.0 * p * cos(phi);
        e[2] = q + 2.0 * p * cos(phi + 2.0 * M_PI / 3.0);
        e[1] = 3.0 * q - e[0] - e[2];
    }

    /* sort by ascending absolute value (insertion sort over three elements) */
    for (i = 1; i < 3; i++)
    {
        t = e[i];
        j = i - 1;
        while (j >= 0 && fabs(e[j]) > fabs(t))
        {
            e[j + 1] = e[j];
            j--;
        }
        e[j + 1] = t;
    }
    eval[0] = e[0];
    eval[1] = e[1];
    eval[2] = e[2];

    if (!evec3)
        return;

    /* null space of (A - eval[2] I) */
    {
        const double l = eval[2];
        double row[3][3];
        double best[3] = {0.0, 0.0, 0.0};
        double best_norm = -1.0;
        int pairs[3][2] = {{0, 1}, {0, 2}, {1, 2}};

        row[0][0] = axx - l; row[0][1] = axy;     row[0][2] = axz;
        row[1][0] = axy;     row[1][1] = ayy - l; row[1][2] = ayz;
        row[2][0] = axz;     row[2][1] = ayz;     row[2][2] = azz - l;

        for (i = 0; i < 3; i++)
        {
            const double *u = row[pairs[i][0]];
            const double *v = row[pairs[i][1]];
            double c[3], n;

            c[0] = u[1] * v[2] - u[2] * v[1];
            c[1] = u[2] * v[0] - u[0] * v[2];
            c[2] = u[0] * v[1] - u[1] * v[0];
            n = c[0] * c[0] + c[1] * c[1] + c[2] * c[2];
            if (n > best_norm)
            {
                best_norm = n;
                best[0] = c[0];
                best[1] = c[1];
                best[2] = c[2];
            }
        }

        if (best_norm > 1e-30)
        {
            double n = sqrt(best_norm);
            evec3[0] = best[0] / n;
            evec3[1] = best[1] / n;
            evec3[2] = best[2] / n;
        }
        else
        {
            /* isotropic neighbourhood: no meaningful normal */
            evec3[0] = 0.0;
            evec3[1] = 0.0;
            evec3[2] = 0.0;
        }
    }
}

/**
 * \brief Multi-scale Hessian sheetness (plate) filter.
 *
 * For every scale the volume is Gaussian-smoothed, the Hessian is formed from
 * central differences in mm units and gamma-normalized by sigma^2 (the correct
 * normalization for second derivatives, so that responses are comparable
 * across scales), and the shape measures
 *
 *   R_sheet = |l2| / |l3|                          (0 for a plate, 1 for a tube)
 *   R_blob  = | 2|l3| - |l2| - |l1| | / |l3|       (2 for a plate, 0 for a blob)
 *   R_noise = sqrt(l1^2 + l2^2 + l3^2)
 *
 * are combined into
 *
 *   S = exp(-R_sheet^2 / 2a^2) * (1 - exp(-R_blob^2 / 2b^2))
 *                              * (1 - exp(-R_noise^2 / 2c^2))
 *
 * The polarity guard rejects the wrong sign of l3: a bright sheet on a darker
 * background is a ridge (l3 < 0), a dark sheet such as sulcal CSF inside GM is
 * a valley (l3 > 0).  This is what keeps a WM-blade detector from responding
 * to sulci and vice versa.
 *
 * The maximum over scales is kept, together with the sheet normal at the
 * winning scale, and opts->gain is applied to the result.  See the header for
 * why a gain is worth having: every consumer gates on the response through a
 * threshold, the one in CAT_VolOrientedMedian() is hard at 0.5, and the
 * automatic noise scale routinely leaves the cortical response well below it.
 *
 * \param src        (in)  scalar volume, e.g. a bias-corrected T1
 * \param sheetness  (out) float[nvox] response in [0,1]; maximum over scales
 * \param normal     (out) float[3*nvox] unit sheet normal, or NULL to skip
 * \param mask       (in)  optional uint8 mask; voxels with 0 are skipped, NULL = all
 * \param dims       (in)  {nx, ny, nz}
 * \param voxelsize  (in)  voxel spacing in mm {dx, dy, dz}
 * \param opts       (in)  filter parameters; NULL selects the defaults
 * \return 0 on success, -1 on invalid arguments, -2 on allocation failure
 */
int CAT_VolSheetness(const float *src, float *sheetness, float *normal,
                     const unsigned char *mask, int dims[3], double voxelsize[3],
                     const CAT_SheetnessOpts *opts)
{
    CAT_SheetnessOpts defaults;
    const int nx = dims ? dims[0] : 0;
    const int ny = dims ? dims[1] : 0;
    const int nz = dims ? dims[2] : 0;
    const int xy = nx * ny;
    const int nvox = xy * nz;
    float *work = NULL;
    double vx, vy, vz;
    double a2, b2;
    int i, s;

    if (!src || !sheetness || !dims || !voxelsize || nvox <= 0)
        return -1;
    if (nx < 3 || ny < 3 || nz < 3)
        return -1;

    if (!opts)
    {
        CAT_SheetnessOptionsInit(&defaults);
        opts = &defaults;
    }

    vx = fabs(voxelsize[0]);
    vy = fabs(voxelsize[1]);
    vz = fabs(voxelsize[2]);
    if (vx <= 0.0 || vy <= 0.0 || vz <= 0.0)
        return -1;

    work = (float *)malloc(sizeof(float) * (size_t)nvox);
    if (!work)
        return -2;

    for (i = 0; i < nvox; i++)
        sheetness[i] = 0.0f;
    if (normal)
    {
        for (i = 0; i < 3 * nvox; i++)
            normal[i] = 0.0f;
    }

    a2 = 2.0 * opts->alpha * opts->alpha;
    b2 = 2.0 * opts->beta * opts->beta;
    if (a2 <= 0.0)
        a2 = 2.0 * 0.25;
    if (b2 <= 0.0)
        b2 = 2.0 * 0.25;

    for (s = 0; s < (opts->n_scales < 1 ? 1 : opts->n_scales); s++)
    {
        double sigma, fwhm[3], gamma, c2;
        double hmax = 0.0;
        int x, y, z;

        if (opts->n_scales <= 1)
            sigma = opts->sigma_min;
        else
            sigma = opts->sigma_min *
                    pow(opts->sigma_max / opts->sigma_min,
                        (double)s / (double)(opts->n_scales - 1));
        if (sigma <= 0.0)
            sigma = 0.3;

        memcpy(work, src, sizeof(float) * (size_t)nvox);
        fwhm[0] = fwhm[1] = fwhm[2] = SQRT8LN2 * sigma;
        smooth3(work, dims, voxelsize, fwhm, 0, DT_FLOAT32);

        /* gamma-normalization for second derivatives */
        gamma = sigma * sigma;

        /* Automatic noise scale: half the largest norm of the diagonal Hessian
           entries found at this scale.  The diagonal alone is enough for a
           scale factor and avoids a second pass over the mixed derivatives.
           This is the usual Frangi heuristic and is what makes the filter
           independent of the intensity units of the input. */
        if (opts->c > 0.0)
        {
            c2 = 2.0 * opts->c * opts->c;
        }
        else
        {
            for (z = 1; z < nz - 1; z++)
                for (y = 1; y < ny - 1; y++)
                    for (x = 1; x < nx - 1; x++)
                    {
                        const int idx = x + y * nx + z * xy;
                        double h[6], n2;
                        if (mask && !mask[idx])
                            continue;
                        h[0] = gamma * (work[idx + 1] - 2.0 * work[idx] + work[idx - 1]) / (vx * vx);
                        h[3] = gamma * (work[idx + nx] - 2.0 * work[idx] + work[idx - nx]) / (vy * vy);
                        h[5] = gamma * (work[idx + xy] - 2.0 * work[idx] + work[idx - xy]) / (vz * vz);
                        n2 = h[0] * h[0] + h[3] * h[3] + h[5] * h[5];
                        if (n2 > hmax)
                            hmax = n2;
                    }
            hmax = sqrt(hmax);
            if (hmax <= 0.0)
                hmax = 1.0;
            c2 = 2.0 * (0.5 * hmax) * (0.5 * hmax);
        }
        if (c2 <= 0.0)
            c2 = 1.0;

        for (z = 1; z < nz - 1; z++)
            for (y = 1; y < ny - 1; y++)
                for (x = 1; x < nx - 1; x++)
                {
                    const int idx = x + y * nx + z * xy;
                    double h[6], eval[3], evec[3];
                    double l1, l2, l3, rs, rb, rn, S;

                    if (mask && !mask[idx])
                        continue;

                    /* second derivatives in mm^-2, gamma-normalized */
                    h[0] = gamma * (work[idx + 1] - 2.0 * work[idx] + work[idx - 1]) / (vx * vx);
                    h[3] = gamma * (work[idx + nx] - 2.0 * work[idx] + work[idx - nx]) / (vy * vy);
                    h[5] = gamma * (work[idx + xy] - 2.0 * work[idx] + work[idx - xy]) / (vz * vz);
                    h[1] = gamma * (work[idx + 1 + nx] - work[idx + 1 - nx] - work[idx - 1 + nx] + work[idx - 1 - nx]) / (4.0 * vx * vy);
                    h[2] = gamma * (work[idx + 1 + xy] - work[idx + 1 - xy] - work[idx - 1 + xy] + work[idx - 1 - xy]) / (4.0 * vx * vz);
                    h[4] = gamma * (work[idx + nx + xy] - work[idx + nx - xy] - work[idx - nx + xy] + work[idx - nx - xy]) / (4.0 * vy * vz);

                    CAT_EigenSym3(h, eval, normal ? evec : NULL);

                    l1 = eval[0];
                    l2 = eval[1];
                    l3 = eval[2];

                    /* a plate needs a dominant third eigenvalue */
                    if (fabs(l3) < 1e-12)
                        continue;

                    /* polarity: +1 wants a ridge (l3 < 0), -1 a valley (l3 > 0) */
                    if (opts->polarity > 0 && l3 >= 0.0)
                        continue;
                    if (opts->polarity < 0 && l3 <= 0.0)
                        continue;

                    rs = fabs(l2) / fabs(l3);
                    rb = fabs(2.0 * fabs(l3) - fabs(l2) - fabs(l1)) / fabs(l3);
                    rn = l1 * l1 + l2 * l2 + l3 * l3;

                    S = exp(-(rs * rs) / a2) *
                        (1.0 - exp(-(rb * rb) / b2)) *
                        (1.0 - exp(-rn / c2));

                    if (S > (double)sheetness[idx])
                    {
                        sheetness[idx] = (float)S;
                        if (normal)
                        {
                            normal[3 * idx + 0] = (float)evec[0];
                            normal[3 * idx + 1] = (float)evec[1];
                            normal[3 * idx + 2] = (float)evec[2];
                        }
                    }
                }

        if (opts->verbose)
            fprintf(stderr, "  sheetness scale %d/%d: sigma = %.2f mm\n",
                    s + 1, opts->n_scales, sigma);
    }

    /* Percentile normalization, applied before the gain.
     *
     * The automatic noise scale c is half the largest Hessian norm in the
     * volume, so the absolute level of the response is set by whatever the
     * strongest structure in that particular image happens to be -- the skull,
     * a vessel, the ventricle wall.  A thin sulcal sheet is far weaker than
     * any of those, which is why the raw response on one dataset can sit an
     * order of magnitude below the same fixed threshold that works on another,
     * and why every consumer used to need a hand-tuned gain before it did
     * anything at all.
     *
     * Anchoring the map to its own p99.9 removes that dependence: the response
     * still ranks voxels exactly as before -- the scaling is a single positive
     * factor, so the ordering, the winning scale and the zero set are all
     * untouched -- but the *value* at a given rank is now the same number on
     * every dataset, and a threshold expressed against it means the same thing
     * everywhere.  p99.9 rather than the maximum because a maximum is a single
     * voxel and behaves like one; p99.9 sits at the top of the genuine sulcal
     * response and is stable.
     *
     * Zeros are excluded, so the anchor is a percentile of the voxels that
     * actually responded rather than of the volume, which would otherwise move
     * with the size of the background.
     *
     * The hazard is the mirror image of the one it fixes: an image containing
     * no sheets at all still has a p99.9, and normalizing to it would amplify
     * noise into a full-scale response.  Set normalize <= 0 to keep the raw
     * response when that matters. */
    if (opts->normalize > 0.0)
    {
        double pct[2] = {50.0, 99.9}, val[2] = {0.0, 0.0};

        get_prctile(sheetness, nvox, val, pct, 1, DT_FLOAT32);

        if (val[1] > 1e-12)
        {
            const double scale = opts->normalize / val[1];

            for (i = 0; i < nvox; i++)
                sheetness[i] = (float)((double)sheetness[i] * scale);

            if (opts->verbose)
                fprintf(stderr, "  normalized: raw p99.9 = %.4f -> %.2f "
                                "(scale %.1fx)\n",
                        val[1], opts->normalize, scale);
        }
        else if (opts->verbose)
        {
            fprintf(stderr, "  normalization skipped: the response is empty.\n");
        }
    }

    /* Overall gain, applied once the maximum over scales is complete.  A
       positive gain is monotonic and therefore cannot change which scale won,
       so one final pass is equivalent to scaling inside the loop and cheaper.
       The clamp keeps the map in [0,1]; because the map is linear in the
       response and fixes zero, a voxel with no response stays at zero for any
       gain and every oriented operator still degenerates exactly to its
       isotropic counterpart there.  A negative gain zeroes the map, which is
       the fully isotropic case. */
    if (opts->gain != 1.0 || opts->normalize > 0.0)
    {
        for (i = 0; i < nvox; i++)
        {
            double g = (double)sheetness[i] * opts->gain;
            if (g < 0.0)
                g = 0.0;
            if (g > 1.0)
                g = 1.0;
            sheetness[i] = (float)g;
        }
    }

    free(work);
    return 0;
}

/* comparison for qsort over floats */
static int cmp_float(const void *a, const void *b)
{
    const float fa = *(const float *)a;
    const float fb = *(const float *)b;
    if (fa < fb)
        return -1;
    if (fa > fb)
        return 1;
    return 0;
}

/**
 * \brief Median filter over a sheetness-oriented neighbourhood.
 *
 * See the header for the admission rule and for how the cutoff maps onto the
 * three neighbour classes.  The important property is the degenerate case:
 * with sheetness 0 the test reduces to 0 < cutoff for every offset, so all 27
 * neighbours are admitted whatever the cutoff, the operator is exactly the
 * isotropic median it replaces, and behaviour away from thin structures does
 * not change at all.
 *
 * \param vol       (in/out) volume to filter, in place
 * \param sheetness (in)     float[nvox] in [0,1]; NULL falls back to isotropic
 * \param normal    (in)     float[3*nvox] unit sheet normals; NULL falls back
 * \param mask      (in)     optional uint8 mask; NULL = filter everywhere
 * \param dims      (in)     {nx, ny, nz}
 * \param cutoff    (in)     admission cutoff; <= 0 selects the default
 * \param iters     (in)     number of successive passes
 * \return 0 on success, -1 on invalid arguments, -2 on allocation failure
 */
int CAT_VolOrientedMedian(float *vol, const float *sheetness, const float *normal,
                          const unsigned char *mask, int dims[3], double cutoff,
                          int iters)
{
    const int nx = dims ? dims[0] : 0;
    const int ny = dims ? dims[1] : 0;
    const int nz = dims ? dims[2] : 0;
    const int xy = nx * ny;
    const int nvox = xy * nz;
    float *in = NULL;
    int it, x, y, z, i, j, k;

    if (!vol || !dims || nvox <= 0)
        return -1;
    if (nx < 3 || ny < 3 || nz < 3 || iters < 1)
        return 0;

    if (cutoff <= 0.0)
        cutoff = CAT_ORIENTED_MEDIAN_CUTOFF;

    in = (float *)malloc(sizeof(float) * (size_t)nvox);
    if (!in)
        return -2;

    for (it = 0; it < iters; it++)
    {
        memcpy(in, vol, sizeof(float) * (size_t)nvox);

        for (z = 1; z < nz - 1; z++)
            for (y = 1; y < ny - 1; y++)
                for (x = 1; x < nx - 1; x++)
                {
                    const int idx = x + y * nx + z * xy;
                    float buf[27];
                    int n = 0;
                    double s = 0.0, nvec[3] = {0.0, 0.0, 0.0};

                    if (mask && !mask[idx])
                        continue;

                    if (sheetness && normal)
                    {
                        s = (double)sheetness[idx];
                        nvec[0] = (double)normal[3 * idx + 0];
                        nvec[1] = (double)normal[3 * idx + 1];
                        nvec[2] = (double)normal[3 * idx + 2];
                        if (s < 0.0)
                            s = 0.0;
                        if (s > 1.0)
                            s = 1.0;
                        /* no usable normal -> stay isotropic */
                        if (nvec[0] == 0.0 && nvec[1] == 0.0 && nvec[2] == 0.0)
                            s = 0.0;
                    }

                    for (k = -1; k <= 1; k++)
                        for (j = -1; j <= 1; j++)
                            for (i = -1; i <= 1; i++)
                            {
                                double d2, cosine;

                                if (i == 0 && j == 0 && k == 0)
                                {
                                    buf[n++] = in[idx];
                                    continue;
                                }

                                d2 = (double)(i * i + j * j + k * k);
                                cosine = ((double)i * nvec[0] + (double)j * nvec[1] +
                                          (double)k * nvec[2]);
                                cosine = cosine * cosine / d2;

                                /* s = 0 admits every offset at any cutoff, which is
                                   what keeps this identical to the isotropic median
                                   away from thin structures */
                                if (s * cosine < cutoff)
                                    buf[n++] = in[idx + i + j * nx + k * xy];
                            }

                    qsort(buf, (size_t)n, sizeof(float), cmp_float);
                    vol[idx] = (n & 1) ? buf[n / 2]
                                       : 0.5f * (buf[n / 2 - 1] + buf[n / 2]);
                }
    }

    free(in);
    return 0;
}
