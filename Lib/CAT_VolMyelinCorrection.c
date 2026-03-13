/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 */

/**
 * \file CAT_VolMyelinCorrection.c
 * \brief Correct PVE label bias at tissue boundaries using T1w intensities.
 *
 * See CAT_VolMyelinCorrection.h for the rationale and algorithm overview.
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <string.h>

#include "CAT_VolMyelinCorrection.h"
#include "CAT_Vol.h"
#include "CAT_Math.h"

/* ------------------------------------------------------------------ */

/**
 * \brief Initialize myelination correction options with defaults.
 *
 * \param opts (out) options structure to fill
 */
void CAT_MyelinCorrOptionsInit(CAT_MyelinCorrOptions *opts)
{
    if (!opts)
        return;
    opts->erosion_mm = 3.0;
    opts->k_intensity = 1.0;
    opts->grad_percentile = 25.0;
    opts->dist_mm = 2.0;
    opts->max_correction = 0.5;
    opts->min_cluster_mm3 = 5.0;
    opts->n_median_filter = 1;
    opts->correct_wm = 1;
    opts->correct_csf = 1;
    opts->verbose = 0;
}

/* ------------------------------------------------------------------ */

/**
 * \brief Correct one tissue boundary (WM/GM or GM/CSF).
 *
 * Builds a deep tissue core, computes intensity statistics, identifies
 * misclassified boundary voxels via three-criterion conjunction, and
 * writes correction magnitudes into the output array.
 *
 * For the WM side (side == +1):
 *   - Core: PVE > 2.9, eroded
 *   - Band: PVE in (2.3, 3.0] outside the core
 *   - Flag: T1w < core_mean - k*core_std  (too dark for WM)
 *   - Shift: PVE toward GM (decrease)
 *
 * For the CSF side (side == -1):
 *   - Core: PVE < 1.1 and PVE > 0, eroded
 *   - Band: PVE in [1.0, 1.7) outside the core
 *   - Flag: T1w > core_mean + k*core_std  (too bright for CSF)
 *   - Shift: PVE toward GM (increase)
 *
 * \param pve         (in)     PVE label volume
 * \param t1w         (in)     raw T1w intensity volume
 * \param grad_mag    (in)     T1w gradient magnitude (precomputed)
 * \param dims        (in)     volume dimensions
 * \param voxelsize   (in)     voxel sizes in mm
 * \param correction  (out)    per-voxel correction magnitude (added to)
 * \param side        (in)     +1 for WM/GM, -1 for GM/CSF
 * \param opts        (in)     algorithm options
 * \return number of flagged voxels, or negative on error
 */
static int
correct_boundary(const float *pve, const float *t1w,
                 const float *grad_mag, int dims[3], double voxelsize[3],
                 float *correction, int side,
                 const CAT_MyelinCorrOptions *opts)
{
    int nvox = dims[0] * dims[1] * dims[2];
    int i, n_core, n_band, n_flagged;
    float *tissue_core = NULL;
    unsigned char *core_mask = NULL;
    unsigned char *band_mask = NULL;
    unsigned char *brain_mask = NULL;
    float *dist_to_core = NULL;
    double core_mean, core_std, intensity_thresh, grad_thresh;

    const char *label = (side > 0) ? "WM/GM" : "GM/CSF";

    /* --- Allocate --- */
    tissue_core = (float *)calloc(nvox, sizeof(float));
    core_mask = (unsigned char *)calloc(nvox, sizeof(unsigned char));
    band_mask = (unsigned char *)calloc(nvox, sizeof(unsigned char));
    brain_mask = (unsigned char *)calloc(nvox, sizeof(unsigned char));
    dist_to_core = (float *)calloc(nvox, sizeof(float));

    if (!tissue_core || !core_mask || !band_mask || !brain_mask || !dist_to_core)
    {
        free(tissue_core);
        free(core_mask);
        free(band_mask);
        free(brain_mask);
        free(dist_to_core);
        fprintf(stderr, "Memory allocation error in correct_boundary (%s)\n", label);
        return -2;
    }

    /* ===== Step 1: Build deep tissue core ===== */
    if (side > 0)
    {
        /* WM core: PVE > 2.9 */
        for (i = 0; i < nvox; i++)
            tissue_core[i] = (pve[i] > 2.9f) ? 1.0f : 0.0f;
    }
    else
    {
        /* CSF core: PVE > 0 and PVE < 1.1 (exclude background = 0) */
        for (i = 0; i < nvox; i++)
            tissue_core[i] = (pve[i] > 0.0f && pve[i] < 1.1f) ? 1.0f : 0.0f;
    }

    dist_erode_float(tissue_core, dims, voxelsize, opts->erosion_mm, 0.5, NULL);
    keep_largest_cluster(tissue_core, 0.5, dims, DT_FLOAT32, 0, 1, 18);

    n_core = 0;
    for (i = 0; i < nvox; i++)
    {
        core_mask[i] = (tissue_core[i] > 0.5f) ? 1 : 0;
        if (core_mask[i])
            n_core++;
    }

    if (n_core < 100)
    {
        if (opts->verbose)
            fprintf(stderr, "  %s: tissue core too small (%d voxels), skipping.\n",
                    label, n_core);
        free(tissue_core);
        free(core_mask);
        free(band_mask);
        free(brain_mask);
        free(dist_to_core);
        return 0;
    }

    if (opts->verbose)
        fprintf(stderr, "  %s: tissue core = %d voxels\n", label, n_core);

    /* ===== Step 2: Core intensity statistics ===== */
    core_mean = get_masked_mean_array((void *)t1w, nvox, core_mask, DT_FLOAT32);
    core_std = get_masked_std_array((void *)t1w, nvox, core_mask, DT_FLOAT32);

    if (opts->verbose)
        fprintf(stderr, "  %s: core T1w mean = %.2f, std = %.2f\n",
                label, core_mean, core_std);

    /* Intensity threshold:
     *   WM side: flag if T1w < core_mean - k*std  (too dark for WM)
     *   CSF side: flag if T1w > core_mean + k*std  (too bright for CSF)
     */
    if (side > 0)
        intensity_thresh = core_mean - opts->k_intensity * core_std;
    else
        intensity_thresh = core_mean + opts->k_intensity * core_std;

    /* ===== Step 3: Distance from tissue core ===== */
    for (i = 0; i < nvox; i++)
    {
        dist_to_core[i] = tissue_core[i];
        brain_mask[i] = (pve[i] > 0.0f) ? 1 : 0;
    }
    euclidean_distance(dist_to_core, brain_mask, dims, voxelsize, 0);

    /* ===== Step 4: Boundary band ===== */
    n_band = 0;
    if (side > 0)
    {
        /* WM boundary band: PVE in (2.3, 3.0] outside core */
        for (i = 0; i < nvox; i++)
        {
            band_mask[i] = (pve[i] > 2.3f && pve[i] <= 3.0f &&
                            !core_mask[i])
                               ? 1
                               : 0;
            if (band_mask[i])
                n_band++;
        }
    }
    else
    {
        /* CSF boundary band: PVE in [1.0, 1.7) outside core */
        for (i = 0; i < nvox; i++)
        {
            band_mask[i] = (pve[i] >= 1.0f && pve[i] < 1.7f &&
                            !core_mask[i])
                               ? 1
                               : 0;
            if (band_mask[i])
                n_band++;
        }
    }

    if (n_band < 10)
    {
        if (opts->verbose)
            fprintf(stderr, "  %s: boundary band too small (%d voxels), skipping.\n",
                    label, n_band);
        free(tissue_core);
        free(core_mask);
        free(band_mask);
        free(brain_mask);
        free(dist_to_core);
        return 0;
    }

    /* ===== Step 5: Gradient threshold from band statistics ===== */
    {
        double *band_grads = (double *)malloc(sizeof(double) * n_band);
        if (!band_grads)
        {
            free(tissue_core);
            free(core_mask);
            free(band_mask);
            free(brain_mask);
            free(dist_to_core);
            return -2;
        }

        int cnt = 0;
        for (i = 0; i < nvox; i++)
        {
            if (band_mask[i])
                band_grads[cnt++] = (double)grad_mag[i];
        }

        double prctile_thresholds[2], prctile_vals[2];
        prctile_vals[0] = opts->grad_percentile;
        prctile_vals[1] = 100.0;
        get_prctile_double(band_grads, cnt, prctile_thresholds, prctile_vals, 0);
        grad_thresh = prctile_thresholds[0];

        if (opts->verbose)
            fprintf(stderr, "  %s: gradient threshold (%.0f%% of band) = %.4f\n",
                    label, opts->grad_percentile, grad_thresh);

        free(band_grads);
    }

    /* ===== Step 6: Flag voxels and compute correction ===== */
    n_flagged = 0;
    for (i = 0; i < nvox; i++)
    {
        if (!band_mask[i])
            continue;

        float t1_val = t1w[i];
        float gm_val = grad_mag[i];
        float dist_val = dist_to_core[i];

        /* Conjunction of three criteria */
        int crit_intensity, crit_gradient, crit_distance;

        if (side > 0)
            crit_intensity = (t1_val < intensity_thresh); /* too dark for WM */
        else
            crit_intensity = (t1_val > intensity_thresh); /* too bright for CSF */

        crit_gradient = (gm_val < grad_thresh);
        crit_distance = (dist_val > opts->dist_mm);

        if (crit_intensity && crit_gradient && crit_distance)
        {
            /* Weight by how strongly each criterion is satisfied [0, 1] */
            double w_int, w_grad, w_dist;

            if (side > 0)
                w_int = fmin(1.0, (intensity_thresh - t1_val) /
                                      (opts->k_intensity * core_std + 1e-10));
            else
                w_int = fmin(1.0, (t1_val - intensity_thresh) /
                                      (opts->k_intensity * core_std + 1e-10));

            w_grad = fmin(1.0, (grad_thresh - gm_val) /
                                   (grad_thresh + 1e-10));
            w_dist = fmin(1.0, (dist_val - opts->dist_mm) /
                                   (opts->dist_mm + 1e-10));

            double w_combined = cbrt(w_int * w_grad * w_dist);

            double pve_excess, shift;
            if (side > 0)
            {
                /* WM side: how far above GM (2.0) */
                pve_excess = pve[i] - GM;
                shift = opts->max_correction * w_combined *
                        fmin(1.0, pve_excess / (WM - GM));
            }
            else
            {
                /* CSF side: how far below GM (2.0) */
                pve_excess = GM - pve[i];
                shift = opts->max_correction * w_combined *
                        fmin(1.0, pve_excess / (GM - CSF));
            }

            correction[i] += (float)shift;
            n_flagged++;
        }
    }

    if (opts->verbose)
        fprintf(stderr, "  %s: flagged %d voxels (%.1f%% of band)\n",
                label, n_flagged, 100.0 * n_flagged / (double)n_band);

    /* ===== Step 7: Remove small isolated clusters of flagged voxels ===== */
    if (n_flagged > 0 && opts->min_cluster_mm3 > 0.0)
    {
        /* Convert to voxel count threshold */
        double vox_vol = voxelsize[0] * voxelsize[1] * voxelsize[2];
        int min_vox = (int)ceil(opts->min_cluster_mm3 / vox_vol);
        if (min_vox < 2)
            min_vox = 2;

        /* Build binary float volume of flagged voxels */
        float *flag_vol = (float *)calloc(nvox, sizeof(float));
        if (flag_vol)
        {
            for (i = 0; i < nvox; i++)
                flag_vol[i] = (correction[i] > 0.0f) ? 1.0f : 0.0f;

            /* Remove clusters smaller than min_vox (18-connected) */
            keep_largest_cluster(flag_vol, 0.5, dims, DT_FLOAT32,
                                 min_vox, 1, 18);

            /* Zero out corrections for voxels removed by cluster filter */
            int n_removed = 0;
            for (i = 0; i < nvox; i++)
            {
                if (correction[i] > 0.0f && flag_vol[i] < 0.5f)
                {
                    correction[i] = 0.0f;
                    n_removed++;
                }
            }
            n_flagged -= n_removed;

            if (opts->verbose)
                fprintf(stderr, "  %s: removed %d isolated voxels "
                                "(min cluster = %d vox / %.1f mm^3)\n",
                        label, n_removed, min_vox, opts->min_cluster_mm3);

            free(flag_vol);
        }
    }

    /* Cleanup */
    free(tissue_core);
    free(core_mask);
    free(band_mask);
    free(brain_mask);
    free(dist_to_core);

    return n_flagged;
}

/* ------------------------------------------------------------------ */

/**
 * \brief Correct PVE labels at both tissue boundaries using T1w intensities.
 *
 * See the header file for the full algorithm description.
 *
 * \param pve        (in/out) PVE label volume (float, values in [0..3])
 * \param t1w        (in)     raw T1w intensity volume (float)
 * \param dims       (in)     volume dimensions {nx, ny, nz}
 * \param voxelsize  (in)     voxel sizes in mm {dx, dy, dz}
 * \param opts       (in)     algorithm options (NULL for defaults)
 * \return 0 on success, non-zero on error
 */
int CAT_VolCorrectMyelination(float *pve, const float *t1w,
                              int dims[3], double voxelsize[3],
                              const CAT_MyelinCorrOptions *opts)
{
    CAT_MyelinCorrOptions defaults;
    int nvox, i, rc;
    float *grad_mag = NULL;
    float *correction_wm = NULL;
    float *correction_csf = NULL;

    if (!pve || !t1w || !dims || !voxelsize)
        return -1;

    if (!opts)
    {
        CAT_MyelinCorrOptionsInit(&defaults);
        opts = &defaults;
    }

    if (!opts->correct_wm && !opts->correct_csf)
        return 0;

    nvox = dims[0] * dims[1] * dims[2];

    /* Compute T1w gradient magnitude (shared by both boundaries) */
    grad_mag = (float *)calloc(nvox, sizeof(float));
    if (!grad_mag)
    {
        fprintf(stderr, "Memory allocation error in CAT_VolCorrectMyelination\n");
        return -2;
    }
    gradient3D((float *)t1w, grad_mag, NULL, NULL, NULL, dims, voxelsize);

    /* ===== WM/GM boundary (myelination correction) ===== */
    if (opts->correct_wm)
    {
        correction_wm = (float *)calloc(nvox, sizeof(float));
        if (!correction_wm)
        {
            free(grad_mag);
            return -2;
        }

        if (opts->verbose)
            fprintf(stderr, "Myelin correction: WM/GM boundary\n");

        rc = correct_boundary(pve, t1w, grad_mag, dims, voxelsize,
                              correction_wm, +1, opts);
        if (rc < 0)
        {
            free(grad_mag);
            free(correction_wm);
            return rc;
        }

        if (rc > 0)
        {
            /* Median-filter and apply: shift PVE toward GM (decrease) */
            if (opts->n_median_filter > 0)
                median3(correction_wm, NULL, dims,
                        opts->n_median_filter, DT_FLOAT32);
            for (i = 0; i < nvox; i++)
            {
                if (correction_wm[i] > 0.0f)
                {
                    float new_pve = pve[i] - correction_wm[i];
                    if (new_pve < GM)
                        new_pve = GM;
                    pve[i] = new_pve;
                }
            }
        }
        free(correction_wm);
    }

    /* ===== GM/CSF boundary (pial correction) ===== */
    if (opts->correct_csf)
    {
        correction_csf = (float *)calloc(nvox, sizeof(float));
        if (!correction_csf)
        {
            free(grad_mag);
            return -2;
        }

        if (opts->verbose)
            fprintf(stderr, "Myelin correction: GM/CSF boundary\n");

        rc = correct_boundary(pve, t1w, grad_mag, dims, voxelsize,
                              correction_csf, -1, opts);
        if (rc < 0)
        {
            free(grad_mag);
            free(correction_csf);
            return rc;
        }

        if (rc > 0)
        {
            /* Median-filter and apply: shift PVE toward GM (increase) */
            if (opts->n_median_filter > 0)
                median3(correction_csf, NULL, dims,
                        opts->n_median_filter, DT_FLOAT32);
            for (i = 0; i < nvox; i++)
            {
                if (correction_csf[i] > 0.0f)
                {
                    float new_pve = pve[i] + correction_csf[i];
                    if (new_pve > GM)
                        new_pve = GM;
                    pve[i] = new_pve;
                }
            }
        }
        free(correction_csf);
    }

    free(grad_mag);
    return 0;
}
