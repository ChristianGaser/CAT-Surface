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

#include "CAT_PpmSulci.h"
#include "CAT_Sheetness.h"
#include "CAT_Vol.h"
#include "CAT_Math.h"

/**
 * \brief Fill a CAT_PpmSulciOpts with defaults for a 0..1 PPM at 0.5 mm.
 *
 * The band of 0.25 keeps the operation to values in (0.5, 0.75] for the usual
 * isovalue: a buried sulcus sits just above the threshold, whereas solid white
 * matter is close to 1 and must never be touched.  The margin of 0.05 puts an
 * opened valley just far enough below the isovalue for marching cubes to
 * separate the banks without carving a wide gap.
 *
 * \param opts (out) option block to initialize; NULL is ignored
 */
void CAT_PpmSulciOptionsInit(CAT_PpmSulciOpts *opts)
{
    if (!opts)
        return;

    opts->sigma_factor = 0.9;
    opts->sigma_min = 0.3;
    opts->sigma_max = 3.0;
    opts->n_scales = 3;
    opts->sheet_normalize = CAT_SHEETNESS_NORMALIZE;
    opts->sheet_skeleton = 0;
    opts->sheet_strength = 30.0;
    opts->thresh = 0.3;
    opts->band = 0.25;
    opts->margin = 0.05;
    opts->strength = 0.0;
    opts->offset = 0.6;
    /* Negative: the raising half uses the same offset, i.e. the signed map
       is applied whole.  See the header for why the two halves are worth
       separating and what the asymmetry measures. */
    opts->offset_gyri = -1.0;
    opts->cutoff = 0.4;
    opts->verbose = 0;
}

/**
 * \brief Push buried sulcal valleys in a PPM below the isovalue.
 *
 * See the header for the three gating conditions.  The local-minimum test is
 * the important one: it is evaluated along the sheet normal rather than over
 * the whole neighbourhood, because that is what distinguishes a sulcal valley
 * floor from an ordinary dark voxel, and it is what makes the operation
 * structurally unable to fire on a gyral crown -- crossing a blade the PPM has
 * a maximum along the normal, never a minimum.
 *
 * \param ppm       (in/out) percentage position map in [0,1], corrected in place
 * \param sheetness (in)     precomputed dark-sheet response, or NULL to compute
 * \param normal    (in)     precomputed unit sheet normals, or NULL to compute
 * \param dims      (in)     {nx, ny, nz}
 * \param voxelsize (in)     voxel spacing in mm
 * \param isovalue  (in)     threshold the surface will be extracted at
 * \param opts      (in)     parameters; NULL selects the defaults
 * \return 0 on success, negative on error
 */
int CAT_VolOpenPpmSulci(float *ppm, const float *sheetness, const float *normal,
                        int dims[3], double voxelsize[3], double isovalue,
                        const CAT_PpmSulciOpts *opts)
{
    CAT_PpmSulciOpts defaults;
    CAT_SheetnessOpts sopts;
    const int nx = dims ? dims[0] : 0;
    const int ny = dims ? dims[1] : 0;
    const int nz = dims ? dims[2] : 0;
    const int xy = nx * ny;
    const int nvox = xy * nz;
    const float *S = sheetness;
    const float *nrm = normal;
    float *own_S = NULL, *own_nrm = NULL;
    int x, y, z, rc = 0;
    long n_changed = 0;

    if (!ppm || !dims || !voxelsize || nvox <= 0)
        return -1;
    if (nx < 3 || ny < 3 || nz < 3)
        return -1;

    if (!opts)
    {
        CAT_PpmSulciOptionsInit(&defaults);
        opts = &defaults;
    }

    if (!S || !nrm)
    {
        own_S = (float *)malloc(sizeof(float) * (size_t)nvox);
        own_nrm = (float *)malloc(sizeof(float) * 3 * (size_t)nvox);
        if (!own_S || !own_nrm)
        {
            free(own_S);
            free(own_nrm);
            return -2;
        }

        CAT_SheetnessOptionsInit(&sopts);
        sopts.sigma_min = opts->sigma_min;
        sopts.sigma_max = opts->sigma_max;
        sopts.n_scales = opts->n_scales;
        sopts.gain = opts->sheet_strength;
        sopts.normalize = opts->sheet_normalize;
        sopts.skeletonize = opts->sheet_skeleton;
        sopts.polarity = -1; /* a sulcus is a valley in the PPM */
        sopts.verbose = opts->verbose;

        if (opts->verbose)
            fprintf(stderr, "Open buried sulci: dark-sheet filter on the PPM.\n");

        rc = CAT_VolSheetness(ppm, own_S, own_nrm, NULL, dims, voxelsize, &sopts);
        if (rc != 0)
        {
            free(own_S);
            free(own_nrm);
            return rc;
        }
        S = own_S;
        nrm = own_nrm;
    }

    for (z = 1; z < nz - 1; z++)
        for (y = 1; y < ny - 1; y++)
            for (x = 1; x < nx - 1; x++)
            {
                const int idx = x + y * nx + z * xy;
                double s, w, target, n0, n1, n2;
                int px, py, pz, mx, my, mz;

                /* only marginal values: above the isovalue but not solid WM */
                if (!((double)ppm[idx] > isovalue &&
                      (double)ppm[idx] < isovalue + opts->band))
                    continue;

                s = (double)S[idx];
                if (s <= opts->thresh)
                    continue;

                n0 = (double)nrm[3 * idx + 0];
                n1 = (double)nrm[3 * idx + 1];
                n2 = (double)nrm[3 * idx + 2];
                if (n0 == 0.0 && n1 == 0.0 && n2 == 0.0)
                    continue;

                px = x + (int)floor(n0 + 0.5);
                py = y + (int)floor(n1 + 0.5);
                pz = z + (int)floor(n2 + 0.5);
                mx = x - (int)floor(n0 + 0.5);
                my = y - (int)floor(n1 + 0.5);
                mz = z - (int)floor(n2 + 0.5);

                if (px == x && py == y && pz == z)
                    continue;
                if (px < 0 || px >= nx || py < 0 || py >= ny || pz < 0 || pz >= nz ||
                    mx < 0 || mx >= nx || my < 0 || my >= ny || mz < 0 || mz >= nz)
                    continue;

                /* a valley floor, not a ridge: this is what a gyral crown can
                   never satisfy, so the test doubles as the gyrus guard */
                if (!(ppm[idx] < ppm[px + py * nx + pz * xy] &&
                      ppm[idx] < ppm[mx + my * nx + mz * xy]))
                    continue;

                w = opts->strength * (s - opts->thresh) / (1.0 - opts->thresh);
                if (w < 0.0)
                    w = 0.0;
                if (w > 1.0)
                    w = 1.0;

                target = isovalue - opts->margin;
                if (target < 0.0)
                    target = 0.0;

                /* one-sided: only ever lower a value */
                if (target >= (double)ppm[idx])
                    continue;

                ppm[idx] = (float)((1.0 - w) * (double)ppm[idx] + w * target);
                n_changed++;
            }

    if (opts->verbose)
        fprintf(stderr, "  opened %ld voxels below the isovalue.\n", n_changed);

    free(own_S);
    free(own_nrm);
    return 0;
}

/**
 * \brief Median cortical thickness read out of a PPM, in mm.
 *
 * See the header for why this works: the map climbs by 1 across the ribbon, so
 * its gradient is 1/thickness.  The band is restricted to the interior because
 * the map flattens at both ends, and the statistic is a median so that the
 * glued sulci -- where the map is stretched across two banks and the gradient
 * reads too small -- cannot drag the estimate down.
 *
 * \param ppm       (in) percentage position map in [0,1]
 * \param dims      (in) {nx, ny, nz}
 * \param voxelsize (in) voxel spacing in mm
 * \return median thickness in mm, or 0 when the band is too small to measure
 */
double CAT_PpmMedianThickness(const float *ppm, int dims[3], double voxelsize[3])
{
    const int nx = dims ? dims[0] : 0;
    const int ny = dims ? dims[1] : 0;
    const int nz = dims ? dims[2] : 0;
    const int xy = nx * ny;
    double vx, vy, vz, med = 0.0;
    double *buf = NULL;
    int x, y, z, n = 0, cap;

    if (!ppm || !dims || !voxelsize || nx < 3 || ny < 3 || nz < 3)
        return 0.0;

    vx = fabs(voxelsize[0]);
    vy = fabs(voxelsize[1]);
    vz = fabs(voxelsize[2]);
    if (vx <= 0.0 || vy <= 0.0 || vz <= 0.0)
        return 0.0;

    /* subsampled by two per axis: ample for a median, and it keeps the
       temporary well under the size of the volume itself */
    cap = (nx / 2 + 1) * (ny / 2 + 1) * (nz / 2 + 1);
    buf = (double *)malloc(sizeof(double) * (size_t)cap);
    if (!buf)
        return 0.0;

    for (z = 2; z < nz - 2; z += 2)
        for (y = 2; y < ny - 2; y += 2)
            for (x = 2; x < nx - 2; x += 2)
            {
                const int i = x + y * nx + z * xy;
                double gx, gy, gz, g;

                /* the interior of the ribbon only: outside it the map is flat
                   and its gradient carries no thickness information */
                if (!(ppm[i] > 0.15f && ppm[i] < 0.85f))
                    continue;

                gx = 0.5 * ((double)ppm[i + 1] - (double)ppm[i - 1]) / vx;
                gy = 0.5 * ((double)ppm[i + nx] - (double)ppm[i - nx]) / vy;
                gz = 0.5 * ((double)ppm[i + xy] - (double)ppm[i - xy]) / vz;

                g = sqrt(gx * gx + gy * gy + gz * gz);
                if (g < 1e-6)
                    continue;

                if (n < cap)
                    buf[n++] = 1.0 / g;
            }

    if (n >= 100)
        med = get_median_double(buf, n, 0);

    free(buf);
    return med;
}
