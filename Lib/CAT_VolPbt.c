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
#include <math.h>
#include <float.h>

#include "CAT_VolPbt.h"
#include "CAT_Vol.h"
#include "CAT_Math.h"

/* Tissue class thresholds (from CAT_Vol.h) */
#ifndef CGM
#define CGM 1.5
#endif
#ifndef GWM
#define GWM 2.5
#endif

void CAT_PbtOptionsInit(CAT_PbtOptions *opts)
{
    if (!opts)
        return;
    opts->n_avgs = 2;
    opts->n_median_filter = 2;
    opts->range = 0.45;
    opts->fill_thresh = 0.5;
    opts->correct_voxelsize = 0.5;
    opts->fast = 0;
    opts->verbose = 0;
}

int CAT_VolComputePbt(
    const float *src,
    float *GMT_out,
    float *PPM_out,
    float *dist_CSF_out,
    float *dist_WM_out,
    int dims[3],
    double voxelsize[3],
    const CAT_PbtOptions *opts)
{
    int i, j;
    int nvox;
    int n_avgs, n_median_filter;
    double range, fill_thresh, correct_voxelsize;
    int verbose;
    float mean_vx_size;
    double add_value, sum_dist;
    double s[3], threshold[2], prctile[2];
    int replace = 0;

    unsigned char *mask = NULL;
    float *input = NULL;
    float *dist_CSF = NULL;
    float *dist_WM = NULL;
    float *GMT = NULL;
    float *GMT1 = NULL;
    float *GMT2 = NULL;
    float *PPM = NULL;
    float *gyrus_mask = NULL;
    float *src_copy = NULL;

    if (!src || !GMT_out || !PPM_out || !dims || !voxelsize || !opts)
        return -1;

    nvox = dims[0] * dims[1] * dims[2];
    mean_vx_size = (voxelsize[0] + voxelsize[1] + voxelsize[2]) / 3.0f;

    /* Copy options (handle fast mode) */
    n_avgs = opts->n_avgs;
    n_median_filter = opts->n_median_filter;
    range = opts->range;
    fill_thresh = opts->fill_thresh;
    correct_voxelsize = opts->correct_voxelsize;
    verbose = opts->verbose;

    if (opts->fast)
    {
        n_avgs /= 2;
        n_median_filter = 0;
        fill_thresh = 0.0;
    }
    if (n_avgs < 1)
        n_avgs = 1;

    /* Allocate working arrays */
    mask = (unsigned char *)malloc(sizeof(unsigned char) * nvox);
    input = (float *)malloc(sizeof(float) * nvox);
    dist_CSF = (float *)malloc(sizeof(float) * nvox);
    dist_WM = (float *)malloc(sizeof(float) * nvox);
    GMT = (float *)malloc(sizeof(float) * nvox);
    GMT1 = (float *)malloc(sizeof(float) * nvox);
    GMT2 = (float *)malloc(sizeof(float) * nvox);
    PPM = (float *)malloc(sizeof(float) * nvox);
    gyrus_mask = (float *)malloc(sizeof(float) * nvox);
    src_copy = (float *)malloc(sizeof(float) * nvox);

    if (!mask || !input || !dist_CSF || !dist_WM || !GMT || !GMT1 || !GMT2 || !PPM || !gyrus_mask || !src_copy)
    {
        if (mask)
            free(mask);
        if (input)
            free(input);
        if (dist_CSF)
            free(dist_CSF);
        if (dist_WM)
            free(dist_WM);
        if (GMT)
            free(GMT);
        if (GMT1)
            free(GMT1);
        if (GMT2)
            free(GMT2);
        if (PPM)
            free(PPM);
        if (gyrus_mask)
            free(gyrus_mask);
        if (src_copy)
            free(src_copy);
        return -2;
    }

    /* Copy source and initialize distances */
    for (i = 0; i < nvox; i++)
    {
        src_copy[i] = src[i];
        dist_CSF[i] = 0.0f;
        dist_WM[i] = 0.0f;
    }

    /* Median-filtering of input */
    localstat3(src_copy, NULL, dims, 1, F_MEDIAN, 1, 1, DT_FLOAT32);

    /* Distance estimation loop */
    for (j = 0; j < n_avgs; j++)
    {
        add_value = ((double)j + 1.0) / ((double)n_avgs + 1.0) - 0.5;

        /* CSF distance map */
        for (i = 0; i < nvox; i++)
        {
            input[i] = (src_copy[i] < (CGM + add_value)) ? 1.0f : 0.0f;
            mask[i] = (src_copy[i] < (GWM + range)) ? 1 : 0;
        }
        if (verbose && (j == 0))
            fprintf(stderr, "Estimate CSF distance map.\n");
        euclidean_distance(input, mask, dims, NULL, replace);
        for (i = 0; i < nvox; i++)
            dist_CSF[i] += input[i];

        /* WM distance map */
        for (i = 0; i < nvox; i++)
        {
            input[i] = (src_copy[i] > (GWM + add_value)) ? 1.0f : 0.0f;
            mask[i] = (src_copy[i] > (CGM - range)) ? 1 : 0;
        }
        if (verbose && (j == 0))
            fprintf(stderr, "Estimate WM distance map.\n");
        euclidean_distance(input, mask, dims, NULL, replace);
        for (i = 0; i < nvox; i++)
            dist_WM[i] += input[i];
    }

    /* Average distances */
    if (n_avgs > 1)
    {
        for (i = 0; i < nvox; i++)
        {
            dist_CSF[i] /= (float)n_avgs;
            dist_WM[i] /= (float)n_avgs;
        }
    }

    /* Thickness estimation using sulci reconstruction */
    if (verbose)
        fprintf(stderr, "Estimate thickness map.\n");
    for (i = 0; i < nvox; i++)
        input[i] = roundf(src_copy[i]);
    projection_based_thickness(input, dist_WM, dist_CSF, GMT1, dims, voxelsize);

    /* Gyri reconstruction (inverse of src) */
    for (i = 0; i < nvox; i++)
        input[i] = roundf(4.0f - src_copy[i]);
    projection_based_thickness(input, dist_CSF, dist_WM, GMT2, dims, voxelsize);

    /* Combine sulci/gyri estimates */
    for (i = 0; i < nvox; i++)
    {
        sum_dist = dist_WM[i] + dist_CSF[i];
        GMT1[i] = fminf(sum_dist, fmaxf(0.0f, GMT1[i] - 0.125f * (GMT1[i] < sum_dist)));

        /* Limit GMT2 to thick regions */
        GMT2[i] = (GMT2[i] > 0.0f) * fmaxf(GMT2[i], 1.75f / mean_vx_size);
        GMT2[i] = fminf(sum_dist, fmaxf(0.0f, GMT2[i] - 0.125f * (GMT2[i] < sum_dist)));
    }

    /* Fill values using Euclidean distance */
    if (verbose)
        fprintf(stderr, "Fill values using Euclidean distance approach\n");
    euclidean_distance(GMT1, NULL, dims, NULL, 1);
    euclidean_distance(GMT2, NULL, dims, NULL, 1);

    /* Iterative median filter */
    median_subsample3(GMT1, dims, voxelsize, 4, 2.0, DT_FLOAT32);
    median_subsample3(GMT2, dims, voxelsize, 4, 2.0, DT_FLOAT32);

    /* Use minimum of thickness measures */
    for (i = 0; i < nvox; i++)
        GMT[i] = fminf(GMT1[i], GMT2[i]);

    for (i = 0; i < nvox; i++)
        mask[i] = (GMT[i] > 1.0f) ? 1 : 0;
    median3(GMT, mask, dims, 3, DT_FLOAT32);

    /* Re-estimate CSF distance */
    for (i = 0; i < nvox; i++)
    {
        if ((src_copy[i] > CGM) && (src_copy[i] < GWM) && (GMT[i] > 1e-15f))
            dist_CSF[i] = fminf(dist_CSF[i], GMT[i] - dist_WM[i]);
    }
    clip_data(dist_CSF, nvox, 0.0, 1e15, DT_FLOAT32);

    /* PPM estimation */
    for (i = 0; i < nvox; i++)
        PPM[i] = (src_copy[i] >= GWM) ? 1.0f : 0.0f;

    if (verbose)
        fprintf(stderr, "Estimate percentage position map.\n");
    smooth_gyri_mask(src_copy, gyrus_mask, dims, voxelsize, CGM, 8.0);

    for (i = 0; i < nvox; i++)
    {
        if ((src_copy[i] > CGM) && (src_copy[i] < GWM) && (GMT[i] > 1e-15f))
        {
            float PPM_sulci = fmaxf(0.0f, GMT[i] - dist_WM[i]) / GMT[i];
            float PPM_gyri = fmaxf(0.0f, GMT1[i] - dist_WM[i]) / GMT1[i];
            PPM[i] = gyrus_mask[i] * PPM_sulci + (1.0f - gyrus_mask[i]) * PPM_gyri;
        }
    }
    clip_data(PPM, nvox, 0.0, 1.0, DT_FLOAT32);

    /* Fill holes */
    if (fill_thresh > 0.0)
        fill_holes(PPM, dims, fill_thresh, 1.0, DT_FLOAT32);

    /* Median-filtering of PPM */
    if (!opts->fast)
        localstat3(PPM, NULL, dims, 1, F_MEDIAN, 2, 1, DT_FLOAT32);

    /* Voxel size correction */
    for (i = 0; i < nvox; i++)
    {
        GMT[i] += correct_voxelsize;
        GMT[i] *= mean_vx_size;
    }

    /* Median filter preprocessing for topology artifact reduction */
    if (n_median_filter)
    {
        float *vol_smoothed = (float *)malloc(sizeof(float) * nvox);
        if (vol_smoothed)
        {
            for (i = 0; i < nvox; i++)
                vol_smoothed[i] = PPM[i];

            s[0] = s[1] = s[2] = 4.0;
            smooth3(vol_smoothed, dims, voxelsize, s, 0, DT_FLOAT32);

            for (i = 0; i < nvox; i++)
                vol_smoothed[i] = PPM[i] - vol_smoothed[i];

            prctile[0] = 0.1;
            prctile[1] = 99.0;
            get_prctile(vol_smoothed, nvox, threshold, prctile, 1, DT_FLOAT32);

            for (i = 0; i < nvox; i++)
                vol_smoothed[i] = ((vol_smoothed[i] > threshold[1]) && (GMT[i] > 1.5f)) ? 1.0f : 0.0f;

            morph_close(vol_smoothed, dims, 1, 0.5, DT_FLOAT32);
            morph_open(vol_smoothed, dims, 1, 0.0, 0, DT_FLOAT32);
            morph_dilate(vol_smoothed, dims, 3, 0.0, DT_FLOAT32);

            s[0] = s[1] = s[2] = 3.0;
            smooth3(vol_smoothed, dims, voxelsize, s, 0, DT_FLOAT32);

            for (i = 0; i < nvox; i++)
                input[i] = PPM[i];
            localstat3(input, NULL, dims, 1, F_MEDIAN, n_median_filter, 1, DT_FLOAT32);

            for (i = 0; i < nvox; i++)
                PPM[i] = (1.0f - vol_smoothed[i]) * PPM[i] + vol_smoothed[i] * input[i];

            free(vol_smoothed);
        }
    }

    clip_data(PPM, nvox, 0.0, 1.0, DT_FLOAT32);

    /* Copy results to output */
    for (i = 0; i < nvox; i++)
    {
        GMT_out[i] = GMT[i];
        PPM_out[i] = PPM[i];
    }
    if (dist_CSF_out)
    {
        for (i = 0; i < nvox; i++)
            dist_CSF_out[i] = dist_CSF[i];
    }
    if (dist_WM_out)
    {
        for (i = 0; i < nvox; i++)
            dist_WM_out[i] = dist_WM[i];
    }

    /* Cleanup */
    free(mask);
    free(input);
    free(dist_CSF);
    free(dist_WM);
    free(GMT);
    free(GMT1);
    free(GMT2);
    free(PPM);
    free(gyrus_mask);
    free(src_copy);

    return 0;
}

/**
 * pmax - Calculate a conditional maximum value from a set of voxels.
 *
 * This function is used in the projection_based_thickness process to find the maximum value among
 * the voxels that are in the range of White Matter Distance (WMD), considering certain constraints.
 * It first finds the pure maximum based on several criteria and then calculates the mean of the
 * highest values under the same constraints.
 *
 * Parameters:
 *  - GMT: Array of thickness/WMD values of neighbours.
 *  - PPM: Array of projection values.
 *  - SEG: Array of segmentation values.
 *  - ND: Array of Euclidean distances.
 *  - WMD: White Matter Distance for the current voxel.
 *  - SEGI: Segmentation value of the current voxel.
 *  - sA: Size of the arrays (number of elements to consider).
 *
 * Returns:
 *  The calculated maximum value under the specified conditions.
 *
 * Notes:
 *  The function applies several constraints based on segmentation and distance measures to determine
 *  the relevant maximum value. This includes checking the range of projection, upper and lower distance
 *  boundaries, and segmentation-based conditions.
 */
static float
pmax(const float *GMT, const float *PPM, const float *SEG, const float *ND, const float WMD, const float SEGI, const int sA)
{
    float maximum = WMD;
    int i;

    // Calculate the pure maximum under specified conditions
    for (i = 0; i <= sA; i++)
    {
        if ((GMT[i] < FLT_MAX) && (maximum < GMT[i]) &&           /* thickness/WMD of neighbours should be larger */
            (SEG[i] >= 1.0) && (SEGI > 1.2 && SEGI <= 2.75) &&    /* projection range */
            (((PPM[i] - ND[i] * 1.2) <= WMD)) &&                  /* upper boundary - maximum distance */
            (((PPM[i] - ND[i] * 0.5) > WMD) || (SEG[i] < 1.5)) && /* lower boundary - minimum distance - corrected values outside */
            ((((SEGI * MAX(1.0, MIN(1.2, SEGI - 1.5))) >= SEG[i])) || (SEG[i] < 1.5)))
        { /* for high values will project data over sulcal gaps */
            maximum = GMT[i];
        }
    }

    // Calculate the mean of the highest values under the same conditions
    float maximum2 = maximum, m2n = 0.0;
    for (i = 0; i <= sA; i++)
    {
        if ((GMT[i] < FLT_MAX) && ((maximum - 1) < GMT[i]) &&
            (SEG[i] >= 1.0) && (SEGI > 1.2 && SEGI <= 2.75) &&
            (((PPM[i] - ND[i] * 1.2) <= WMD)) &&
            (((PPM[i] - ND[i] * 0.5) > WMD) || (SEG[i] < 1.5)) &&
            ((((SEGI * MAX(1.0, MIN(1.2, SEGI - 1.5))) >= SEG[i])) || (SEG[i] < 1.5)))
        {
            maximum2 += GMT[i];
            m2n++;
        }
    }
    if (m2n > 0.0)
        maximum = (maximum2 - maximum) / m2n;

    return maximum;
}

/**
 * \brief Public API for projection_based_thickness.
 *
 * This function is part of the CAT-Surface public library interface and is used by command-line tools.
 *
 * \param SEG (in/out) Parameter of projection_based_thickness.
 * \param WMD (in/out) Parameter of projection_based_thickness.
 * \param CSFD (in/out) Parameter of projection_based_thickness.
 * \param GMT (in/out) Parameter of projection_based_thickness.
 * \param dims (in/out) Parameter of projection_based_thickness.
 * \param voxelsize (in/out) Parameter of projection_based_thickness.
 * \return void (no return value).
 */
void projection_based_thickness(float *SEG, float *WMD, float *CSFD, float *GMT, int dims[3], double *voxelsize)
{
    // Initialization and pre-processing
    const int nvox = dims[0] * dims[1] * dims[2];
    const int x = dims[0], y = dims[1], xy = x * y;
    const float s2 = sqrt(2.0), s3 = sqrt(3.0);
    const int NI[14] = {0, -1, -x + 1, -x, -x - 1, -xy + 1, -xy, -xy - 1, -xy + x + 1, -xy + x, -xy + x - 1, -xy - x + 1, -xy - x, -xy - x - 1}; // Neighbour index offsets
    const float ND[14] = {0.0, 1.0, s2, 1.0, s2, s2, 1.0, s2, s3, s2, s3, s3, s2, s3};                                                           // Neighbour distances
    enum
    {
        sN = (int)(sizeof(NI) / sizeof(NI[0]))
    }; // Number of neighbours

    // Variables for processing
    float DN[sN], DI[sN], GMTN[sN], WMDN[sN], SEGN[sN], DNm;
    float GMTi, CSFDi;
    int i, n, ni, u, v, w, nu, nv, nw, count_WM = 0, count_CSF = 0;

    /* GMT should be allocated */
    if (!GMT)
    {
        GMT = (float *)malloc(sizeof(float) * nvox);
        if (!GMT)
        {
            fprintf(stderr, "Memory allocation error\n");
            exit(EXIT_FAILURE);
        }
    }

    // Initial distance checks and assignment
    for (i = 0; i < nvox; i++)
    {
        // Initial GMT value and WM/CSF counts
        GMT[i] = WMD[i] + 0.0;

        // proof distance input
        if (SEG[i] >= GWM)
            count_WM++;
        if (SEG[i] <= CGM)
            count_CSF++;
    }

    // Error checks for WM and CSF voxels
    if (count_WM == 0)
    {
        fprintf(stderr, "ERROR: no WM voxels\n");
        exit(EXIT_FAILURE);
    }
    if (count_CSF == 0)
    {
        fprintf(stderr, "ERROR: no CSF voxels\n");
        exit(EXIT_FAILURE);
    }

    // Forward thickness calculation
    for (i = 0; i < nvox; i++)
    {
        // Process only GM voxels
        if (SEG[i] > CSF && SEG[i] < WM)
        {
            // Neighbourhood processing
            ind2sub(i, &u, &v, &w, xy, x);

            /* read neighbour values */
            for (n = 0; n < sN; n++)
            {
                ni = i + NI[n];
                ind2sub(ni, &nu, &nv, &nw, xy, x);

                if ((ni < 0) || (ni >= nvox) || (abs(nu - u) > 1) || (abs(nv - v) > 1) || (abs(nw - w) > 1))
                    ni = i;

                GMTN[n] = GMT[ni];
                WMDN[n] = WMD[ni];
                SEGN[n] = SEG[ni];
            }

            /* find minimum distance within the neighbourhood */
            DNm = pmax(GMTN, WMDN, SEGN, ND, WMD[i], SEG[i], sN);
            GMT[i] = DNm;
        }
    }

    // Backward search for thickness correction
    for (i = nvox - 1; i >= 0; i--)
    {
        // Process only GM voxels
        if (SEG[i] > CSF && SEG[i] < WM)
        {
            ind2sub(i, &u, &v, &w, xy, x);

            /* read neighbour values */
            for (n = 0; n < sN; n++)
            {
                ni = i - NI[n];
                ind2sub(ni, &nu, &nv, &nw, xy, x);

                if ((ni < 0) || (ni >= nvox) || (abs(nu - u) > 1) || (abs(nv - v) > 1) || (abs(nw - w) > 1))
                    ni = i;

                GMTN[n] = GMT[ni];
                WMDN[n] = WMD[ni];
                SEGN[n] = SEG[ni];
            }

            /* find minimum distance within the neighbourhood */
            DNm = pmax(GMTN, WMDN, SEGN, ND, WMD[i], SEG[i], sN);
            if ((GMT[i] < DNm) && (DNm > 0))
                GMT[i] = DNm;
        }
    }

    // Post-processing to refine GMT values
    for (i = 0; i < nvox; i++)
    {
        if (SEG[i] < CGM || SEG[i] > GWM)
            GMT[i] = 0.0;
    }

    // Final GMT adjustment based on CSFD and WMD
    for (i = 0; i < nvox; i++)
    {
        if (SEG[i] >= CGM && SEG[i] <= GWM)
        {
            GMTi = CSFD[i] + WMD[i];
            CSFDi = GMT[i] - WMD[i];

            /*if ( CSFD[i]>CSFDi ) CSFD[i] = CSFDi;
            else GMT[i]  = GMTi;*/
        }
    }
}

/**
 * \brief Build a smooth "gyri/sulci mask" (gyri ≈ 0, sulci ≈ 1).
 *
 * This mask supports downstream operations that benefit from different behavior in
 * gyral vs. sulcal regions, e.g.:
 *  - **Surface extraction:** prevent sulcal closure by using a higher isovalue in sulci,
 *    and prevent cutting gyri by using a lower isovalue in gyri.
 *  - **Projection-based cortical thickness:** use different parameters over gyri vs. sulci.
 *
 * **Algorithm (overview)**
 *  1) **Initial thresholding** of the input scalar image `src` at `thresh * max(src)` to form
 *     a coarse tissue mask; then **distance-based closing** to fill sulcal gaps.
 *  2) **Gyri emphasis:** a slight **dilation** followed by a stronger **erosion** to
 *     preferentially shrink gyri crowns relative to sulcal regions.
 *  3) **Smoothing** to create a soft (0..1) transition between gyri and sulci.
 *  4) **CSF enforcement:** voxels clearly below tissue threshold (e.g., CSF) are set to 1
 *     (open sulci), followed by a light final smoothing to avoid hard borders.
 *
 * **Conventions**
 *  - Output `mask` is in **[0,1]** (float). Values close to **0** indicate gyri crowns;
 *    values close to **1** indicate sulcal fundi.
 *
 * \param src       (in)  float[nx*ny*nz] input scalar image (e.g., label map)
 * \param mask      (out) float[nx*ny*nz] output mask (0..1); pre-allocated by caller
 * \param dims      (in)  {nx, ny, nz}
 * \param voxelsize (in)  voxel spacing (e.g., mm) as {sx, sy, sz}
 * \param thresh    (in)  threshold for src; seeds the initial mask (recommended value for label map 1.5)
 * \param fwhm      (in)  smoothing FWHM for the main blur step (recommended value FWHM=8.0)
 *
 * \note The heuristic constants (closing 5, dilate 2, erode 5, CSF factor 0.75, final FWHM 2)
 *       follow the original intent and can be exposed as parameters if needed.
 */
void smooth_gyri_mask(const float *src, float *mask,
                      int dims[3], double voxelsize[3],
                      double thresh, double fwhm)
{
    const int nx = dims[0], ny = dims[1], nz = dims[2];
    const int nvox = nx * ny * nz;
    int i;

    /* Check inputs */
    if (!src || !mask || nvox <= 0)
        return;

    /* mask = (src > thresh) ? 1 : 0 */
    for (i = 0; i < nvox; ++i)
        mask[i] = (src[i] > thresh) ? 1.0f : 0.0f;

    /* Close sulcal gaps with a modest distance closing radius */
    dist_close_float(mask, dims, voxelsize, /*dist=*/5.0, /*th=*/0.5, NULL);

    /* Gyri emphasis: slight dilation then stronger erosion */
    /* This sequence tends to suppress gyri crowns relative to sulci. */
    dist_dilate_float(mask, dims, voxelsize, /*dist=*/2.0, /*th=*/0.5, NULL);
    dist_erode_float(mask, dims, voxelsize, /*dist=*/5.0, /*th=*/0.5, NULL);

    /* Smooth transition between gyri (≈0) and sulci (≈1) */
    double fwhm3[3];
    fwhm3[0] = fwhm3[1] = fwhm3[2] = fwhm;
    smooth3(mask, dims, voxelsize, fwhm3, /*mask=*/0, DT_FLOAT32);

    /* Force obvious CSF to 1 (open sulci), then very light smoothing */
    const float t_csf = 0.75f * thresh; /* slightly below main threshold */
    for (i = 0; i < nvox; ++i)
        if (src[i] < t_csf)
            mask[i] = 1.0f;

    fwhm3[0] = fwhm3[1] = fwhm3[2] = fwhm / 4.0;
    smooth3(mask, dims, voxelsize, fwhm3, /*mask=*/0, DT_FLOAT32);
}
