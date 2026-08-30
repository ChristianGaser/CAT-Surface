/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 */

/**
 * \file CAT_VolBoundaryOffset.c
 * \brief Measure how far the GM/WM boundary has been displaced by myelination.
 *
 * See CAT_VolBoundaryOffset.h for the rationale and algorithm overview.
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <string.h>

#include "CAT_VolBoundaryOffset.h"
#include "CAT_Vol.h"
#include "CAT_Math.h"

/* The GM/WM isovalue the label map is read at, and how far along the ray the
 * label is still allowed to cross it before the profile is rejected as not
 * running through this boundary. */
#define BO_ISO 2.5f
#define BO_ANCHOR_MM 1.0

/* Displacements below this are rounding residue, not measurements. */
#define BO_EPS_MM 1e-4f

/* ------------------------------------------------------------------ */

/**
 * \brief Initialize boundary-offset options with defaults.
 *
 * \param opts (out) options structure to fill
 */
void CAT_BoundaryOffsetOptionsInit(CAT_BoundaryOffsetOpts *opts)
{
    if (!opts)
        return;
    opts->ref_fwhm = 10.0;
    opts->erosion_mm = 3.0;
    opts->t_lo = 0.25;
    opts->t_hi = 0.75;
    opts->search_mm = 4.0;
    opts->step_mm = 0.25;
    opts->width_pct = 88.0;
    opts->gain = 0.5;
    opts->max_offset_mm = 1.5;
    opts->smooth_fwhm = 8.0;
    opts->fill_ribbon = 1;
    opts->verbose = 0;
}

/* ------------------------------------------------------------------ */

/**
 * \brief Local tissue-intensity reference by normalized convolution.
 *
 * See the header for what this buys and why the alternative does not work.
 *
 * \param values    (in)  volume to average, typically the T1w
 * \param control   (in)  1 where the voxel is a control point for this tissue
 * \param dims      (in)  volume dimensions
 * \param voxelsize (in)  voxel sizes in mm
 * \param fwhm_mm   (in)  smoothing FWHM in mm
 * \param fallback  (in)  value used where no control point carries any weight
 * \param ref       (out) local reference, one value per voxel
 * \return 0 on success, non-zero on allocation failure
 */
int CAT_VolLocalTissueReference(const float *values,
                                const unsigned char *control,
                                int dims[3], double voxelsize[3],
                                double fwhm_mm, double fallback, float *ref)
{
    int nvox = dims[0] * dims[1] * dims[2];
    int i;
    float *num, *den;
    double fwhm[3];

    if (!values || !control || !ref)
        return -1;

    num = (float *)calloc(nvox, sizeof(float));
    den = (float *)calloc(nvox, sizeof(float));
    if (!num || !den)
    {
        free(num);
        free(den);
        return -2;
    }

    for (i = 0; i < nvox; i++)
    {
        num[i] = control[i] ? values[i] : 0.0f;
        den[i] = control[i] ? 1.0f : 0.0f;
    }

    fwhm[0] = fwhm[1] = fwhm[2] = fwhm_mm;

    /* Smoothed on a coarse grid: at this FWHM a full-resolution Gaussian
     * would be hundreds of taps wide per axis, and the quantity is smooth by
     * construction. */
    smooth_subsample3(num, dims, voxelsize, fwhm, 0, 3.0, DT_FLOAT32);
    smooth_subsample3(den, dims, voxelsize, fwhm, 0, 3.0, DT_FLOAT32);

    /* Below this weight the quotient is carried by the far tail of the kernel
     * and amplifies rounding rather than measuring anything. */
    for (i = 0; i < nvox; i++)
        ref[i] = (den[i] > 1e-3f) ? (num[i] / den[i]) : (float)fallback;

    free(num);
    free(den);
    return 0;
}

/* ------------------------------------------------------------------ */

/**
 * \brief Trilinear sample of a volume at a fractional grid position.
 *
 * \param vol  (in) volume
 * \param dims (in) volume dimensions
 * \param x    (in) fractional grid coordinate
 * \param y    (in) fractional grid coordinate
 * \param z    (in) fractional grid coordinate
 * \param ok   (out) set to 0 when the position falls outside the volume
 * \return interpolated value, or 0 when outside
 */
static double
sample_trilinear(const float *vol, int dims[3],
                 double x, double y, double z, int *ok)
{
    int x0 = (int)floor(x), y0 = (int)floor(y), z0 = (int)floor(z);
    double fx = x - x0, fy = y - y0, fz = z - z0;
    int sx = dims[0], sy = dims[1], sz = dims[2];
    int i, j, k;
    double acc = 0.0;

    if (x0 < 0 || y0 < 0 || z0 < 0 ||
        x0 + 1 >= sx || y0 + 1 >= sy || z0 + 1 >= sz)
    {
        *ok = 0;
        return 0.0;
    }
    *ok = 1;

    for (k = 0; k < 2; k++)
        for (j = 0; j < 2; j++)
            for (i = 0; i < 2; i++)
            {
                double w = (i ? fx : 1.0 - fx) *
                           (j ? fy : 1.0 - fy) *
                           (k ? fz : 1.0 - fz);
                acc += w * (double)vol[(x0 + i) + (y0 + j) * sx +
                                       (z0 + k) * sx * sy];
            }
    return acc;
}

/**
 * \brief First upward crossing of a level in a sampled profile.
 *
 * Linear interpolation between the bracketing samples gives the crossing to
 * better than the sampling step.
 *
 * \param prof  (in) profile samples
 * \param valid (in) 1 where the sample is inside the volume
 * \param n     (in) number of samples
 * \param s0    (in) arc length of the first sample, in mm
 * \param ds    (in) sample spacing in mm
 * \param level (in) level to cross
 * \param found (out) set to 1 when a crossing was located
 * \return arc length of the crossing in mm
 */
static double
first_crossing(const double *prof, const unsigned char *valid, int n,
               double s0, double ds, double level, int *found)
{
    int k;

    *found = 0;
    for (k = 1; k < n; k++)
    {
        if (!valid[k] || !valid[k - 1])
            continue;
        if (prof[k - 1] < level && prof[k] >= level)
        {
            double denom = prof[k] - prof[k - 1];
            double frac = (denom > 1e-12) ? (level - prof[k - 1]) / denom : 0.0;
            *found = 1;
            return s0 + ((double)(k - 1) + frac) * ds;
        }
    }
    return 0.0;
}

/* ------------------------------------------------------------------ */

/**
 * \brief Estimate the myelination-induced GM/WM boundary displacement.
 *
 * See the header file for the full algorithm description.
 *
 * \param label     (in)  PVE label volume (float, values in [0..3])
 * \param t1w       (in)  T1w intensity volume on the same grid
 * \param offset    (out) displacement in mm, 0 outside the transition band
 * \param width     (out) optional raw transition width in mm; NULL to skip
 * \param dims      (in)  volume dimensions
 * \param voxelsize (in)  voxel sizes in mm
 * \param opts      (in)  algorithm options (NULL for defaults)
 * \return number of band voxels with a valid profile, or negative on error
 */
long CAT_VolBoundaryOffset(const float *label, const float *t1w,
                           float *offset, float *width,
                           int dims[3], double voxelsize[3],
                           const CAT_BoundaryOffsetOpts *opts)
{
    CAT_BoundaryOffsetOpts defaults;
    int nvox, i, n_samp;
    long n_band = 0, n_valid = 0;
    float *wm_core = NULL, *ref_wm = NULL, *ref_gm = NULL, *tnorm = NULL;
    float *gx = NULL, *gy = NULL, *gz = NULL, *w_map = NULL;
    unsigned char *ctrl = NULL, *band = NULL;
    double *prof = NULL, *lprof = NULL, *widths = NULL;
    unsigned char *pvalid = NULL;
    double wm_mean, gm_mean, w_ref = 0.0, w_thresh = 0.0;
    int rc = 0;

    if (!label || !t1w || !offset || !dims || !voxelsize)
        return -1;

    if (!opts)
    {
        CAT_BoundaryOffsetOptionsInit(&defaults);
        opts = &defaults;
    }

    nvox = dims[0] * dims[1] * dims[2];
    n_samp = (int)(2.0 * opts->search_mm / opts->step_mm) + 1;

    wm_core = (float *)calloc(nvox, sizeof(float));
    ref_wm = (float *)calloc(nvox, sizeof(float));
    ref_gm = (float *)calloc(nvox, sizeof(float));
    tnorm = (float *)calloc(nvox, sizeof(float));
    gx = (float *)calloc(nvox, sizeof(float));
    gy = (float *)calloc(nvox, sizeof(float));
    gz = (float *)calloc(nvox, sizeof(float));
    w_map = (float *)calloc(nvox, sizeof(float));
    ctrl = (unsigned char *)calloc(nvox, sizeof(unsigned char));
    band = (unsigned char *)calloc(nvox, sizeof(unsigned char));
    prof = (double *)calloc(n_samp, sizeof(double));
    lprof = (double *)calloc(n_samp, sizeof(double));
    pvalid = (unsigned char *)calloc(n_samp, sizeof(unsigned char));

    if (!wm_core || !ref_wm || !ref_gm || !tnorm || !gx || !gy || !gz ||
        !w_map || !ctrl || !band || !prof || !lprof || !pvalid)
    {
        rc = -2;
        goto done;
    }

    memset(offset, 0, sizeof(float) * nvox);
    if (width)
        memset(width, 0, sizeof(float) * nvox);

    /* ===== 1. Local GM and WM intensity levels ===== */

    /* WM control: the deep core, so the myelinated band does not calibrate
     * the reference it is going to be measured against. */
    for (i = 0; i < nvox; i++)
        wm_core[i] = (label[i] > 2.9f) ? 1.0f : 0.0f;
    dist_erode_float(wm_core, dims, voxelsize, opts->erosion_mm, 0.5, NULL);

    for (i = 0; i < nvox; i++)
        ctrl[i] = (wm_core[i] > 0.5f) ? 1 : 0;
    wm_mean = get_masked_mean_array((void *)t1w, nvox, ctrl, DT_FLOAT32);
    if (CAT_VolLocalTissueReference(t1w, ctrl, dims, voxelsize,
                                    opts->ref_fwhm, wm_mean, ref_wm) != 0)
    {
        rc = -2;
        goto done;
    }

    /* GM control: the interior of the PVE ramp needs no erosion. */
    for (i = 0; i < nvox; i++)
        ctrl[i] = (label[i] > 1.75f && label[i] < 2.25f) ? 1 : 0;
    gm_mean = get_masked_mean_array((void *)t1w, nvox, ctrl, DT_FLOAT32);
    if (CAT_VolLocalTissueReference(t1w, ctrl, dims, voxelsize,
                                    opts->ref_fwhm, gm_mean, ref_gm) != 0)
    {
        rc = -2;
        goto done;
    }

    if (opts->verbose)
        fprintf(stderr, "Boundary offset: global GM %.2f, WM %.2f\n",
                gm_mean, wm_mean);

    /* Contrast-normalized intensity: 0 at the local GM level, 1 at the local
     * WM level.  A multiplicative bias moves both ends together and cancels. */
    for (i = 0; i < nvox; i++)
    {
        double c = (double)ref_wm[i] - (double)ref_gm[i];
        tnorm[i] = (c > 1e-6) ? (float)(((double)t1w[i] - ref_gm[i]) / c) : 0.0f;
    }

    /* ===== 2. Cortical normal from the label map ===== */
    {
        float *lab_copy = (float *)malloc(sizeof(float) * nvox);
        if (!lab_copy)
        {
            rc = -2;
            goto done;
        }
        memcpy(lab_copy, label, sizeof(float) * nvox);
        gradient3D(lab_copy, NULL, gx, gy, gz, dims, voxelsize);
        free(lab_copy);
    }

    /* The band is the voxels adjacent to the isosurface, found by looking for
     * a sign change across a face -- not the voxels whose label happens to lie
     * near the isovalue.  The distinction matters: a sharp transition produces
     * few intermediate labels and a wide one produces many, so selecting on
     * the label value samples myelinated cortex far more densely than healthy
     * cortex and biases the reference width the whole method is measured
     * against.  A geometric band samples proportionally to surface area. */
    {
        static const int off6[6][3] = {{-1, 0, 0}, {1, 0, 0}, {0, -1, 0},
                                       {0, 1, 0}, {0, 0, -1}, {0, 0, 1}};
        int x, y, z, k;
        int sx = dims[0], sy = dims[1], sz = dims[2];

        for (z = 1; z < sz - 1; z++)
            for (y = 1; y < sy - 1; y++)
                for (x = 1; x < sx - 1; x++)
                {
                    int side_here;

                    i = x + y * sx + z * sx * sy;
                    if (label[i] <= 0.0f)
                        continue;
                    side_here = (label[i] >= BO_ISO);

                    for (k = 0; k < 6; k++)
                    {
                        int j = (x + off6[k][0]) +
                                (y + off6[k][1]) * sx +
                                (z + off6[k][2]) * sx * sy;

                        if (label[j] <= 0.0f)
                            continue;
                        if (((label[j] >= BO_ISO) ? 1 : 0) != side_here)
                        {
                            band[i] = 1;
                            n_band++;
                            break;
                        }
                    }
                }
    }

    if (n_band < 100)
    {
        if (opts->verbose)
            fprintf(stderr, "Boundary offset: transition band too small "
                            "(%ld voxels).\n", n_band);
        rc = 0;
        goto done;
    }

    widths = (double *)calloc(n_band, sizeof(double));
    if (!widths)
    {
        rc = -2;
        goto done;
    }

    /* ===== 3. Read the profile along the normal at each band voxel ===== */
    {
        int x, y, z, k;
        int sx = dims[0], sy = dims[1];

        for (z = 0; z < dims[2]; z++)
            for (y = 0; y < dims[1]; y++)
                for (x = 0; x < dims[0]; x++)
                {
                    double nx, ny, nz, gmag, s_lo, s_hi, s_lab, w;
                    int f_lo, f_hi, f_lab;

                    i = x + y * sx + z * sx * sy;
                    if (!band[i])
                        continue;

                    gmag = sqrt((double)gx[i] * gx[i] +
                                (double)gy[i] * gy[i] +
                                (double)gz[i] * gz[i]);
                    if (gmag < 1e-6)
                        continue;

                    /* Unit normal in physical space, pointing toward WM. */
                    nx = gx[i] / gmag;
                    ny = gy[i] / gmag;
                    nz = gz[i] / gmag;

                    for (k = 0; k < n_samp; k++)
                    {
                        double s = -opts->search_mm + k * opts->step_mm;
                        double px = x + s * nx / voxelsize[0];
                        double py = y + s * ny / voxelsize[1];
                        double pz = z + s * nz / voxelsize[2];
                        int ok, ok2;

                        prof[k] = sample_trilinear(tnorm, dims, px, py, pz,
                                                   &ok);
                        lprof[k] = sample_trilinear(label, dims, px, py, pz,
                                                    &ok2);
                        pvalid[k] = (unsigned char)(ok && ok2);
                    }

                    /* The profile has to actually cross the boundary, and
                     * cross it here.  Where the normal comes out tangential --
                     * at the edge of the mask, or wherever the label gradient
                     * is dominated by something other than this boundary --
                     * the ray runs along the surface instead of through it and
                     * returns a meaningless width. */
                    s_lab = first_crossing(lprof, pvalid, n_samp,
                                           -opts->search_mm, opts->step_mm,
                                           BO_ISO, &f_lab);
                    if (!f_lab || fabs(s_lab) > BO_ANCHOR_MM)
                        continue;

                    s_lo = first_crossing(prof, pvalid, n_samp,
                                          -opts->search_mm, opts->step_mm,
                                          opts->t_lo, &f_lo);
                    s_hi = first_crossing(prof, pvalid, n_samp,
                                          -opts->search_mm, opts->step_mm,
                                          opts->t_hi, &f_hi);

                    if (!f_lo || !f_hi)
                        continue;

                    w = s_hi - s_lo;
                    if (w <= 0.0)
                        continue;

                    w_map[i] = (float)w;
                    widths[n_valid++] = w;
                }
    }

    if (n_valid < 100)
    {
        if (opts->verbose)
            fprintf(stderr, "Boundary offset: too few usable profiles "
                            "(%ld of %ld band voxels).\n", n_valid, n_band);
        rc = 0;
        goto done;
    }

    /* ===== 4. Smooth the width within the boundary sheet =====
     * Before anything is thresholded, not after.  Myelination is a
     * centimetre-scale, stereotyped phenomenon while the per-voxel profile is
     * noisy, and the band is a thin surface, so a normalized convolution
     * restricted to it averages along the boundary and not across it.
     *
     * Doing this first is what makes the threshold below meaningful: on a
     * 0.75 mm subject it takes the spread of the healthy population from 0.121
     * to 0.080 mm and the p99 of the width from 4.06 to 2.01 mm -- the extreme
     * per-voxel values were noise, not structure.  It also removes the reason
     * to worry about the non-negativity clamp rectifying noise into a
     * systematic inflation, because the noise is gone before the clamp. */
    if (opts->smooth_fwhm > 0.0)
    {
        float *sm = (float *)calloc(nvox, sizeof(float));
        unsigned char *valid_mask = (unsigned char *)calloc(nvox, 1);

        if (!sm || !valid_mask)
        {
            free(sm);
            free(valid_mask);
            rc = -2;
            goto done;
        }
        for (i = 0; i < nvox; i++)
            valid_mask[i] = (w_map[i] > 0.0f) ? 1 : 0;

        if (CAT_VolLocalTissueReference(w_map, valid_mask, dims, voxelsize,
                                        opts->smooth_fwhm, 0.0, sm) == 0)
        {
            for (i = 0; i < nvox; i++)
                w_map[i] = valid_mask[i] ? sm[i] : 0.0f;
        }
        free(sm);
        free(valid_mask);

        n_valid = 0;
        for (i = 0; i < nvox; i++)
            if (w_map[i] > 0.0f)
                widths[n_valid++] = (double)w_map[i];
    }

    {
        double pth[2], pv[2];
        pv[0] = 50.0;
        pv[1] = 100.0;
        get_prctile_double(widths, (int)n_valid, pth, pv, 0);
        w_ref = pth[0];
    }

    /* ===== 5. Where the healthy population ends =====
     * A location alone displaces the entire cortex, and that is the single
     * thing most worth knowing about this method: whatever central percentile
     * is chosen, the healthy population straddles it, so every voxel above
     * acquires a positive offset.  With the reference at p25 that was 99.8% of
     * the boundary on a real subject, at 0.05-0.1 mm, because p25 sits *below*
     * the mode of a distribution whose healthy peak is one voxel wide.
     *
     * The threshold is therefore an upper percentile, which fixes the share of
     * the boundary treated as myelinated.  That is the one parameterization
     * that transfers: a location plus a multiple of the spread does not,
     * because the spread depends on resolution and noise rather than on
     * anatomy -- the same setting selected 12.0% of the boundary on a 0.75 mm
     * scan and 1.4% on a 1x1x1.25 mm one.  A fixed share is also what the
     * anatomy suggests, heavily myelinated cortex (M1, S1, V1, A1, MT) being
     * a roughly fixed tenth of the sheet in everyone.
     *
     * Fixing the share does not fix the correction: the displacement is the
     * excess *beyond* the threshold, so a brain whose widths are tightly
     * grouped is barely corrected at all however large the share.  Measured on
     * two subjects at p88, the mean displacement came out at 0.19 mm and
     * 0.13 mm respectively, and it moves by under 10% across p85-p92. */
    {
        double pth[2], pv[2];

        pv[0] = opts->width_pct;
        pv[1] = 100.0;
        get_prctile_double(widths, (int)n_valid, pth, pv, 0);
        w_thresh = pth[0];
    }

    if (opts->verbose)
        fprintf(stderr, "Boundary offset: %ld profiles, width p50 = %.3f mm, "
                        "threshold (p%.0f) = %.3f mm\n",
                n_valid, w_ref, opts->width_pct, w_thresh);

    /* ===== 6. Excess beyond the healthy population -> displacement ===== */
    for (i = 0; i < nvox; i++)
    {
        double d;

        if (w_map[i] <= 0.0f)
            continue;
        d = opts->gain * ((double)w_map[i] - w_thresh);
        if (d < BO_EPS_MM)
            continue;
        if (d > opts->max_offset_mm)
            d = opts->max_offset_mm;
        offset[i] = (float)d;
    }

    if (opts->verbose)
    {
        double sum = 0.0, mx = 0.0;
        long n = 0;
        for (i = 0; i < nvox; i++)
            if (offset[i] > 0.0f)
            {
                sum += offset[i];
                if (offset[i] > mx)
                    mx = offset[i];
                n++;
            }
        fprintf(stderr, "Boundary offset: %ld of %ld sheet voxels displaced "
                        "(%.1f%%), mean %.3f mm, max %.3f mm\n",
                n, n_valid, 100.0 * (double)n / (double)n_valid,
                n ? sum / n : 0.0, mx);
    }

    if (width)
        memcpy(width, w_map, sizeof(float) * nvox);

    /* ===== 7. Carry the displacement through the tissue =====
     * The profile can only be read where the boundary is, but the quantity
     * being corrected is the thickness of the whole cortical column, and it is
     * read at the central surface -- half a ribbon away from any of these
     * voxels.  Each voxel therefore takes the displacement measured at the
     * nearest point of the sheet.
     *
     * The seeds are every sheet voxel that produced a profile, not just the
     * displaced ones: seeding from the nonzero offsets alone would hand
     * healthy cortex the correction belonging to the myelinated patch next to
     * it. */
    if (opts->fill_ribbon)
    {
        float *seed = (float *)calloc(nvox, sizeof(float));
        float *filled = (float *)calloc(nvox, sizeof(float));
        unsigned char *inside = (unsigned char *)calloc(nvox, 1);

        if (!seed || !filled || !inside)
        {
            free(seed);
            free(filled);
            free(inside);
            rc = -2;
            goto done;
        }

        for (i = 0; i < nvox; i++)
        {
            seed[i] = (w_map[i] > 0.0f) ? 1.0f : 0.0f;
            inside[i] = (label[i] > 0.0f) ? 1 : 0;
        }

        euclidean_distance_src(seed, inside, dims, voxelsize, 0,
                               offset, filled);

        for (i = 0; i < nvox; i++)
            offset[i] = inside[i] ? filled[i] : 0.0f;

        free(seed);
        free(filled);
        free(inside);
    }

    rc = 0;

done:
    free(wm_core);
    free(ref_wm);
    free(ref_gm);
    free(tnorm);
    free(gx);
    free(gy);
    free(gz);
    free(w_map);
    free(ctrl);
    free(band);
    free(prof);
    free(lprof);
    free(pvalid);
    free(widths);

    return (rc != 0) ? rc : n_valid;
}
