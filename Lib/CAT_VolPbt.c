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

void
CAT_PbtOptionsInit(CAT_PbtOptions *opts)
{
    if (!opts) return;
    opts->n_avgs = 2;
    opts->n_median_filter = 2;
    opts->range = 0.45;
    opts->fill_thresh = 0.5;
    opts->correct_voxelsize = 0.5;
    opts->fast = 0;
    opts->verbose = 0;
}

int
CAT_VolComputePbt(
    const float *src,
    float *GMT_out,
    float *PPM_out,
    float *dist_CSF_out,
    float *dist_WM_out,
    int dims[3],
    double voxelsize[3],
    const CAT_PbtOptions *opts
)
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

    if (opts->fast) {
        n_avgs /= 2;
        n_median_filter = 0;
        fill_thresh = 0.0;
    }
    if (n_avgs < 1) n_avgs = 1;

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

    if (!mask || !input || !dist_CSF || !dist_WM || !GMT || !GMT1 || !GMT2 || !PPM || !gyrus_mask || !src_copy) {
        if (mask) free(mask);
        if (input) free(input);
        if (dist_CSF) free(dist_CSF);
        if (dist_WM) free(dist_WM);
        if (GMT) free(GMT);
        if (GMT1) free(GMT1);
        if (GMT2) free(GMT2);
        if (PPM) free(PPM);
        if (gyrus_mask) free(gyrus_mask);
        if (src_copy) free(src_copy);
        return -2;
    }

    /* Copy source and initialize distances */
    for (i = 0; i < nvox; i++) {
        src_copy[i] = src[i];
        dist_CSF[i] = 0.0f;
        dist_WM[i] = 0.0f;
    }

    /* Median-filtering of input */
    localstat3(src_copy, NULL, dims, 1, F_MEDIAN, 1, 1, DT_FLOAT32);

    /* Distance estimation loop */
    for (j = 0; j < n_avgs; j++) {
        add_value = ((double)j + 1.0) / ((double)n_avgs + 1.0) - 0.5;

        /* CSF distance map */
        for (i = 0; i < nvox; i++) {
            input[i] = (src_copy[i] < (CGM + add_value)) ? 1.0f : 0.0f;
            mask[i] = (src_copy[i] < (GWM + range)) ? 1 : 0;
        }
        if (verbose && (j == 0)) fprintf(stderr, "Estimate CSF distance map.\n");
        euclidean_distance(input, mask, dims, NULL, replace);
        for (i = 0; i < nvox; i++)
            dist_CSF[i] += input[i];

        /* WM distance map */
        for (i = 0; i < nvox; i++) {
            input[i] = (src_copy[i] > (GWM + add_value)) ? 1.0f : 0.0f;
            mask[i] = (src_copy[i] > (CGM - range)) ? 1 : 0;
        }
        if (verbose && (j == 0)) fprintf(stderr, "Estimate WM distance map.\n");
        euclidean_distance(input, mask, dims, NULL, replace);
        for (i = 0; i < nvox; i++)
            dist_WM[i] += input[i];
    }

    /* Average distances */
    if (n_avgs > 1) {
        for (i = 0; i < nvox; i++) {
            dist_CSF[i] /= (float)n_avgs;
            dist_WM[i] /= (float)n_avgs;
        }
    }

    /* Thickness estimation using sulci reconstruction */
    if (verbose) fprintf(stderr, "Estimate thickness map.\n");
    for (i = 0; i < nvox; i++) input[i] = roundf(src_copy[i]);
    projection_based_thickness(input, dist_WM, dist_CSF, GMT1, dims, voxelsize);

    /* Gyri reconstruction (inverse of src) */
    for (i = 0; i < nvox; i++) input[i] = roundf(4.0f - src_copy[i]);
    projection_based_thickness(input, dist_CSF, dist_WM, GMT2, dims, voxelsize);

    /* Combine sulci/gyri estimates */
    for (i = 0; i < nvox; i++) {
        sum_dist = dist_WM[i] + dist_CSF[i];
        GMT1[i] = fminf(sum_dist, fmaxf(0.0f, GMT1[i] - 0.125f * (GMT1[i] < sum_dist)));
        GMT2[i] = (GMT2[i] > 0.0f) * fmaxf(GMT2[i], 1.75f / mean_vx_size);
        GMT2[i] = fminf(sum_dist, fmaxf(0.0f, GMT2[i] - 0.125f * (GMT2[i] < sum_dist)));
    }

    /* Fill values using Euclidean distance */
    if (verbose) fprintf(stderr, "Fill values using Euclidean distance approach\n");
    euclidean_distance(GMT1, NULL, dims, NULL, 1);
    euclidean_distance(GMT2, NULL, dims, NULL, 1);

    /* Iterative median filter */
    median_subsample3(GMT1, dims, voxelsize, 4, 2.0, DT_FLOAT32);
    median_subsample3(GMT2, dims, voxelsize, 4, 2.0, DT_FLOAT32);

    /* Use minimum of thickness measures */
    for (i = 0; i < nvox; i++) GMT[i] = fminf(GMT1[i], GMT2[i]);

    for (i = 0; i < nvox; i++)
        mask[i] = (GMT[i] > 1.0f) ? 1 : 0;
    median3(GMT, mask, dims, 3, DT_FLOAT32);

    /* Re-estimate CSF distance */
    for (i = 0; i < nvox; i++) {
        if ((src_copy[i] > CGM) && (src_copy[i] < GWM) && (GMT[i] > 1e-15f))
            dist_CSF[i] = fminf(dist_CSF[i], GMT[i] - dist_WM[i]);
    }
    clip_data(dist_CSF, nvox, 0.0, 1e15, DT_FLOAT32);

    /* PPM estimation */
    for (i = 0; i < nvox; i++)
        PPM[i] = (src_copy[i] >= GWM) ? 1.0f : 0.0f;

    if (verbose) fprintf(stderr, "Estimate percentage position map.\n");
    smooth_gyri_mask(src_copy, gyrus_mask, dims, voxelsize, CGM, 8.0);

    for (i = 0; i < nvox; i++) {
        if ((src_copy[i] > CGM) && (src_copy[i] < GWM) && (GMT[i] > 1e-15f)) {
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
    for (i = 0; i < nvox; i++) {
        GMT[i] += correct_voxelsize;
        GMT[i] *= mean_vx_size;
    }

    /* Median filter preprocessing for topology artifact reduction */
    if (n_median_filter) {
        float *vol_smoothed = (float *)malloc(sizeof(float) * nvox);
        if (vol_smoothed) {
            for (i = 0; i < nvox; i++) vol_smoothed[i] = PPM[i];

            s[0] = s[1] = s[2] = 4.0;
            smooth3(vol_smoothed, dims, voxelsize, s, 0, DT_FLOAT32);

            for (i = 0; i < nvox; i++) vol_smoothed[i] = PPM[i] - vol_smoothed[i];

            prctile[0] = 0.1; prctile[1] = 99.0;
            get_prctile(vol_smoothed, nvox, threshold, prctile, 1, DT_FLOAT32);

            for (i = 0; i < nvox; i++)
                vol_smoothed[i] = ((vol_smoothed[i] > threshold[1]) && (GMT[i] > 1.5f)) ? 1.0f : 0.0f;

            morph_close(vol_smoothed, dims, 1, 0.5, DT_FLOAT32);
            morph_open(vol_smoothed, dims, 1, 0.0, 0, DT_FLOAT32);
            morph_dilate(vol_smoothed, dims, 3, 0.0, DT_FLOAT32);

            s[0] = s[1] = s[2] = 3.0;
            smooth3(vol_smoothed, dims, voxelsize, s, 0, DT_FLOAT32);

            for (i = 0; i < nvox; i++) input[i] = PPM[i];
            localstat3(input, NULL, dims, 1, F_MEDIAN, n_median_filter, 1, DT_FLOAT32);

            for (i = 0; i < nvox; i++)
                PPM[i] = (1.0f - vol_smoothed[i]) * PPM[i] + vol_smoothed[i] * input[i];

            free(vol_smoothed);
        }
    }

    clip_data(PPM, nvox, 0.0, 1.0, DT_FLOAT32);

    /* Copy results to output */
    for (i = 0; i < nvox; i++) {
        GMT_out[i] = GMT[i];
        PPM_out[i] = PPM[i];
    }
    if (dist_CSF_out) {
        for (i = 0; i < nvox; i++)
            dist_CSF_out[i] = dist_CSF[i];
    }
    if (dist_WM_out) {
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
