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
#include "CAT_Nlm.h"
#include "CAT_Sheetness.h"

/* Tissue class thresholds (from CAT_Vol.h) */
#ifndef CGM
#define CGM 1.5
#endif
#ifndef GWM
#define GWM 2.5
#endif

/*
 * Thickness corrections are defined in mm, not in voxels.
 *
 * All distances inside the algorithm are computed in voxel units, so a
 * correction expressed in voxels silently scales with the grid PBT is run on:
 * the same constant removes twice as much cortex at 1 mm as it does at 0.5 mm.
 * Since the bias these constants compensate is a physical property of the
 * segmentation (where it places the GM/WM and GM/CSF borders), it must be a
 * fixed number of millimetres.  Constants defined here are converted to voxels
 * where they are applied to voxel-space quantities.
 *
 * PBT_SHRINK_MM and the default of opts->correct_thickness reproduce the former
 * voxel-based values (0.125 and 0.5 voxel) at the 0.5 mm grid that PBT is
 * normally run on, so results at 0.5 mm are unchanged.
 */
#define PBT_SHRINK_MM 0.0625 /* was 0.125 voxel */
#define PBT_CORRECT_MM 0.25  /* was 0.5 voxel   */
#define PBT_GMT2_MIN_MM 1.75 /* lower bound of the gyral thickness estimate */

/* Forward declarations of static helpers */
static void correct_ppm_sulci(const float *src, float *PPM, float *GMT,
                              const float *dist_CSF, const float *dist_WM,
                              int dims[3], double voxelsize[3],
                              double sulcal_width);
static double estimate_pve_width(const float *src, int dims[3]);

void CAT_PbtOptionsInit(CAT_PbtOptions *opts)
{
    if (!opts)
        return;
    opts->n_avgs = 2;
    opts->n_median_filter = 2;
    opts->range = 0.45;
    opts->median_subsample = 4;
    opts->fill_thresh = 0.5;
    opts->correct_thickness = PBT_CORRECT_MM;
    opts->sulcal_width = 2.5;
    opts->pve_distance = 0;
    opts->sulcal_barrier = 0;
    opts->gyral_barrier = 0;
    opts->barrier_q = 0.0;
    opts->barrier_tmin = 0.5;
    opts->barrier_halfwidth = 0.0;
    opts->oriented_filter = 0;
    opts->oriented_strength = 1.0;
    opts->oriented_cutoff = 0.0; /* 0 selects CAT_ORIENTED_MEDIAN_CUTOFF */
    opts->fast = 0;
    opts->verbose = 0;
}

/**
 * \brief Compute projection-based cortical thickness and percentage position map.
 *
 * Implements the full PBT pipeline consisting of averaged WM/CSF distance
 * estimation, sulcal and gyral thickness estimation, PPM construction, and an
 * optional weighted local median-filter cleanup of the final PPM.
 *
 * If opts->n_median_filter > 0, the median filter is not applied globally.
 * Instead, a topology-artifact likelihood map is estimated from the positive
 * residual PPM - smooth(PPM), restricted to sufficiently thick cortex
 * (GMT > 1.5), regularized morphologically, then smoothed. This soft weight
 * map blends the original PPM with a locally median-filtered PPM so that only
 * likely topology-artifact regions receive strong filtering.
 *
 * \param src            (in)  input PVE label image (CSF=1, GM=2, WM=3)
 * \param GMT_out        (out) output gray matter thickness map
 * \param PPM_out        (out) output percentage position map
 * \param dist_CSF_out   (out) optional output CSF distance map, or NULL
 * \param dist_WM_out    (out) optional output WM distance map, or NULL
 * \param dims           (in)  volume dimensions [nx, ny, nz]
 * \param voxelsize      (in)  voxel sizes in mm [dx, dy, dz]
 * \param opts           (in)  algorithm options, including n_median_filter for
 *                             weighted local PPM cleanup
 * \return 0 on success, non-zero on error
 */
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
    int n_avgs, n_median_filter, subsample;
    double range, fill_thresh, correct_thickness;
    int verbose;
    float mean_vx_size, shrink, gmt2_min;
    double add_value, sum_dist;
    double s[3], threshold[2], prctile[2];
    int replace = 0;

    double pve_width = 1.0;
    float *src_val = NULL;

    /* sheetness field driving the oriented replacements of the isotropic
       medians; NULL keeps every filter isotropic and the result unchanged */
    float *sheet = NULL;
    float *sheet_nrm = NULL;

    /* sulcal medial surface used as a barrier for the CSF distance, and its
       dual: the gyral medial surface bounding the WM distance */
    float *medial = NULL;
    float *gyral = NULL;

    unsigned char *mask = NULL;
    float *input = NULL;
    float *dist_CSF = NULL;
    float *dist_WM = NULL;
    float *GMT = NULL;
    float *GMT1 = NULL;
    float *GMT2 = NULL;
    float *PPM = NULL;
    float *src_copy = NULL;

    if (!src || !GMT_out || !PPM_out || !dims || !voxelsize || !opts)
        return -1;

    nvox = dims[0] * dims[1] * dims[2];
    mean_vx_size = (voxelsize[0] + voxelsize[1] + voxelsize[2]) / 3.0f;
    if (mean_vx_size <= 0.0f)
        return -1;

    /* mm-defined constants expressed in the voxel units used internally */
    shrink = (float)(PBT_SHRINK_MM / mean_vx_size);
    gmt2_min = (float)(PBT_GMT2_MIN_MM / mean_vx_size);

    /* Copy options (handle fast mode) */
    n_avgs = opts->n_avgs;
    n_median_filter = opts->n_median_filter;
    range = opts->range;
    fill_thresh = opts->fill_thresh;
    correct_thickness = opts->correct_thickness;
    verbose = opts->verbose;
    subsample = opts->median_subsample;

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
    src_copy = (float *)malloc(sizeof(float) * nvox);

    if (!mask || !input || !dist_CSF || !dist_WM || !GMT || !GMT1 || !GMT2 || !PPM || !src_copy)
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

    /* Optional orientation field, computed once on the untouched label map and
     * reused by every oriented filter below.  Polarity 0 accepts both signs of
     * the dominant curvature, because both structures we must not destroy are
     * thin sheets: the dark CSF sheet of a tight sulcus and the bright WM
     * blade of a thin gyrus.
     *
     * On a label map the filter also responds along ordinary tissue
     * boundaries, which are planar too.  That is harmless and in fact wanted
     * here: orienting the median along a boundary instead of across it is
     * exactly what makes it edge preserving, so the field does double duty as
     * a thin-structure detector and as an edge orientation field. */
    if (opts->oriented_filter)
    {
        CAT_SheetnessOpts sopts;

        sheet = (float *)malloc(sizeof(float) * nvox);
        sheet_nrm = (float *)malloc(sizeof(float) * 3 * nvox);
        if (!sheet || !sheet_nrm)
        {
            free(sheet);
            free(sheet_nrm);
            sheet = NULL;
            sheet_nrm = NULL;
            fprintf(stderr, "Warning: no memory for the oriented filters, "
                            "falling back to isotropic ones.\n");
        }
        else
        {
            CAT_SheetnessOptionsInit(&sopts);
            sopts.polarity = 0;
            sopts.gain = opts->oriented_strength;
            sopts.verbose = verbose;
            if (verbose)
                fprintf(stderr, "Estimate sheetness field for oriented filtering.\n");
            if (CAT_VolSheetness(src_copy, sheet, sheet_nrm, NULL, dims, voxelsize,
                                 &sopts) != 0)
            {
                free(sheet);
                free(sheet_nrm);
                sheet = NULL;
                sheet_nrm = NULL;
            }
        }
    }

    /* Median-filtering of input */
    if (sheet)
        CAT_VolOrientedMedian(src_copy, sheet, sheet_nrm, NULL, dims,
                              opts->oriented_cutoff, 1);
    else
        localstat3(src_copy, NULL, dims, 1, F_MEDIAN, 1, 1, DT_FLOAT32);

    /* Optional sub-voxel correction of the distance maps.
     *
     * euclidean_distance() measures to the nearest source voxel *centre*, but
     * the boundary the threshold stands for does not run through that centre:
     * a source voxel overshoots the threshold, and how far it lies beyond the
     * boundary is encoded in its partial volume. Converting the overshoot into
     * a distance needs the width of the PVE ramp, because a label difference
     * of 1 spans one voxel only if the ramp is one voxel wide -- which it is
     * not when the label map was interpolated to a finer grid than it was
     * acquired on. This mirrors the correction in CAT12's
     * cat_vol_pbtsimpleCS4.m, where the ramp width is opt.pvet.
     *
     * Off by default: it changes the calibration of the thickness, so
     * opts->correct_thickness has to be re-derived when it is switched on.
     */
    if (opts->pve_distance)
    {
        src_val = (float *)malloc(sizeof(float) * nvox);
        if (!src_val)
        {
            free(mask); free(input); free(dist_CSF); free(dist_WM);
            free(GMT); free(GMT1); free(GMT2); free(PPM); free(src_copy);
            return -2;
        }
        pve_width = estimate_pve_width(src_copy, dims);
        if (verbose)
            fprintf(stderr, "PVE ramp width: %.2f voxel.\n", pve_width);
    }

    /* Distance estimation loop */
    for (j = 0; j < n_avgs; j++)
    {
        double thr_CSF, thr_WM;

        add_value = ((double)j + 1.0) / ((double)n_avgs + 1.0) - 0.5;
        thr_CSF = CGM + add_value;
        thr_WM = GWM + add_value;

        /* CSF distance map */
        for (i = 0; i < nvox; i++)
        {
            input[i] = (src_copy[i] < thr_CSF) ? 1.0f : 0.0f;
            mask[i] = (src_copy[i] < (GWM + range)) ? 1 : 0;
        }
        if (verbose && (j == 0))
            fprintf(stderr, "Estimate CSF distance map.\n");
        euclidean_distance_src(input, mask, dims, NULL, replace,
                               src_val ? src_copy : NULL, src_val);
        if (src_val)
        {
            /* the CSF source sits below the threshold, so it lies this far
               beyond the boundary, capped at half a voxel */
            for (i = 0; i < nvox; i++)
            {
                float d = (float)((thr_CSF - src_val[i]) * pve_width);
                d = fminf(0.5f, fmaxf(0.0f, d));
                dist_CSF[i] += fmaxf(0.0f, input[i] - d);
            }
        }
        else
        {
            for (i = 0; i < nvox; i++)
                dist_CSF[i] += input[i];
        }

        /* WM distance map */
        for (i = 0; i < nvox; i++)
        {
            input[i] = (src_copy[i] > thr_WM) ? 1.0f : 0.0f;
            mask[i] = (src_copy[i] > (CGM - range)) ? 1 : 0;
        }
        if (verbose && (j == 0))
            fprintf(stderr, "Estimate WM distance map.\n");
        euclidean_distance_src(input, mask, dims, NULL, replace,
                               src_val ? src_copy : NULL, src_val);
        if (src_val)
        {
            /* mirrored: the WM source sits above the threshold */
            for (i = 0; i < nvox; i++)
            {
                float d = (float)((src_val[i] - thr_WM) * pve_width);
                d = fminf(0.5f, fmaxf(0.0f, d));
                dist_WM[i] += fmaxf(0.0f, input[i] - d);
            }
        }
        else
        {
            for (i = 0; i < nvox; i++)
                dist_WM[i] += input[i];
        }
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

    /* Sulcal barrier.
     *
     * Where the classifier lost the CSF in a sulcus there is no boundary for
     * the CSF distance to stop at, so the front from one bank runs through the
     * fused grey matter and into the other.  dist_CSF then measures most of the
     * way across the sulcus, GMT follows it, and the PPM never drops below the
     * isovalue -- the buried sulcus is created here, in the distance map, long
     * before marching cubes is asked to draw it.
     *
     * The midline the front should have stopped at is still recoverable, and
     * from geometry rather than intensity: it is where the fronts from the two
     * banks collide.  Bounding dist_CSF by the distance to that surface puts
     * the missing boundary back.
     *
     * Applying it as a minimum is what makes this safe to leave on.  Where CSF
     * was segmented properly it is nearer than the midline, the minimum keeps
     * the real value, and the result is bit-identical to not running this at
     * all.  Where it was lost the bound takes over.  Either way the distance
     * can only shrink, so an overestimated thickness can be corrected but a
     * correct one can never be inflated. */
    if (opts->sulcal_barrier)
    {
        float *dmed = NULL;
        double tmin_vox = opts->barrier_tmin / (double)mean_vx_size;
        double half_vox = (opts->barrier_halfwidth > 0.0)
                              ? opts->barrier_halfwidth / (double)mean_vx_size
                              : 1.0;
        long n_medial;

        medial = (float *)malloc(sizeof(float) * nvox);
        dmed = (float *)malloc(sizeof(float) * nvox);

        if (!medial || !dmed)
        {
            free(medial);
            free(dmed);
            medial = NULL;
            fprintf(stderr, "Warning: no memory for the sulcal barrier.\n");
        }
        else
        {
            /* the cortical band the fronts travel through */
            for (i = 0; i < nvox; i++)
                mask[i] = (src_copy[i] > CGM && src_copy[i] < GWM) ? 1 : 0;

            n_medial = CAT_VolSulcalMedialSet(dist_WM, mask, medial, dims,
                                              opts->barrier_q, tmin_vox);

            if (n_medial <= 0)
            {
                free(medial);
                free(dmed);
                medial = NULL;
                if (verbose)
                    fprintf(stderr, "Sulcal barrier: no collisions found.\n");
            }
            else
            {
                long n_capped = 0;

                /* distance to the medial surface, in voxel units to match the
                   distances it is compared against */
                for (i = 0; i < nvox; i++)
                    dmed[i] = medial[i];
                euclidean_distance(dmed, NULL, dims, NULL, 0);

                for (i = 0; i < nvox; i++)
                {
                    double bound;

                    if (!(src_copy[i] > CGM && src_copy[i] < GWM))
                        continue;

                    bound = (double)dmed[i] - half_vox;
                    if (bound < 0.0)
                        bound = 0.0;

                    if (bound < (double)dist_CSF[i])
                    {
                        dist_CSF[i] = (float)bound;
                        n_capped++;
                    }
                }

                if (verbose)
                    fprintf(stderr, "Sulcal barrier: %ld medial voxels, "
                                    "capped dist_CSF at %ld of them.\n",
                            n_medial, n_capped);
                free(dmed);
            }
        }
    }

    /* Gyral barrier: the exact dual of the block above.
     *
     * A thin white-matter blade lies between two sulci.  Where the classifier
     * lost it the gyrus is grey matter all the way through, so dist_WM measures
     * out to whatever white matter is left further along the gyrus rather than
     * to the blade that should have been here -- the thickness inflates and the
     * PPM at the gyral core never rises, which is the blade disappearing from
     * the surface.
     *
     * Running the fronts inward from the pial boundary instead of outward from
     * the white matter turns the same collision test into a blade detector:
     * where the blade survives, the white matter splits the band and no fronts
     * meet; where it was lost, they collide exactly where it should have been.
     * Nothing else changes -- it is the same routine with the other distance
     * map -- and it is one-sided in the opposite direction, raising the PPM at
     * the core and so strengthening a finger, never thinning one. */
    if (opts->gyral_barrier)
    {
        float *dgyr = NULL;
        double tmin_vox = opts->barrier_tmin / (double)mean_vx_size;
        double half_vox = (opts->barrier_halfwidth > 0.0)
                              ? opts->barrier_halfwidth / (double)mean_vx_size
                              : 1.0;
        long n_gyral;

        gyral = (float *)malloc(sizeof(float) * nvox);
        dgyr = (float *)malloc(sizeof(float) * nvox);

        if (!gyral || !dgyr)
        {
            free(gyral);
            free(dgyr);
            gyral = NULL;
            fprintf(stderr, "Warning: no memory for the gyral barrier.\n");
        }
        else
        {
            for (i = 0; i < nvox; i++)
                mask[i] = (src_copy[i] > CGM && src_copy[i] < GWM) ? 1 : 0;

            n_gyral = CAT_VolSulcalMedialSet(dist_CSF, mask, gyral, dims,
                                             opts->barrier_q, tmin_vox);

            if (n_gyral <= 0)
            {
                free(gyral);
                free(dgyr);
                gyral = NULL;
                if (verbose)
                    fprintf(stderr, "Gyral barrier: no collisions found.\n");
            }
            else
            {
                long n_capped = 0;

                for (i = 0; i < nvox; i++)
                    dgyr[i] = gyral[i];
                euclidean_distance(dgyr, NULL, dims, NULL, 0);

                for (i = 0; i < nvox; i++)
                {
                    double bound;

                    if (!(src_copy[i] > CGM && src_copy[i] < GWM))
                        continue;

                    bound = (double)dgyr[i] - half_vox;
                    if (bound < 0.0)
                        bound = 0.0;

                    if (bound < (double)dist_WM[i])
                    {
                        dist_WM[i] = (float)bound;
                        n_capped++;
                    }
                }

                if (verbose)
                    fprintf(stderr, "Gyral barrier: %ld medial voxels, "
                                    "capped dist_WM at %ld of them.\n",
                            n_gyral, n_capped);
                free(dgyr);
            }
        }
    }

    /* Thickness estimation using sulci reconstruction */
    if (verbose)
        fprintf(stderr, "Estimate thickness map.\n");
    for (i = 0; i < nvox; i++)
        input[i] = roundf(src_copy[i]);

    /* The projection has to stop at the barrier too.  Capping the distance is
     * not enough on its own: projection_based_thickness() carries a thickness
     * along connected grey matter, and with the banks fused it would hand the
     * value straight across the sulcus again.  Marking the medial voxels as CSF
     * in the copy the projection reads stops it there, using the mechanism the
     * routine already has -- it only projects through voxels between CSF and
     * WM.  This touches a scratch buffer, never the segmentation. */
    if (medial)
    {
        for (i = 0; i < nvox; i++)
            if (medial[i] > 0.0f)
                input[i] = (float)CSF;
    }
    if (gyral)
    {
        for (i = 0; i < nvox; i++)
            if (gyral[i] > 0.0f)
                input[i] = (float)WM;
    }

    projection_based_thickness(input, dist_WM, dist_CSF, GMT1, dims, voxelsize);

    /* Gyri reconstruction (inverse of src) */
    for (i = 0; i < nvox; i++)
        input[i] = roundf(4.0f - src_copy[i]);
    projection_based_thickness(input, dist_CSF, dist_WM, GMT2, dims, voxelsize);

    /* Combine sulci/gyri estimates */
    for (i = 0; i < nvox; i++)
    {
        sum_dist = dist_WM[i] + dist_CSF[i];
        GMT1[i] = fminf(sum_dist, fmaxf(0.0f, GMT1[i] - shrink * (GMT1[i] < sum_dist)));

        /* Limit GMT2 to thick regions */
        GMT2[i] = (GMT2[i] > 0.0f) * fmaxf(GMT2[i], gmt2_min);
        GMT2[i] = fminf(sum_dist, fmaxf(0.0f, GMT2[i] - shrink * (GMT2[i] < sum_dist)));
    }

    /* Fill values using Euclidean distance */
    if (verbose)
        fprintf(stderr, "Fill values using Euclidean distance approach\n");
    euclidean_distance(GMT1, NULL, dims, NULL, 1);
    euclidean_distance(GMT2, NULL, dims, NULL, 1);

    /* Iterative median filter */
    median_subsample3(GMT1, dims, voxelsize, subsample, 2.0, DT_FLOAT32);
    median_subsample3(GMT2, dims, voxelsize, subsample, 2.0, DT_FLOAT32);

    /* Use the smaller thickness measure */
    for (i = 0; i < nvox; i++)
    {
        GMT[i] = fminf(GMT1[i], GMT2[i]);    
        mask[i] = (GMT[i] > 1.0f) ? 1 : 0;
    }
    median3(GMT, mask, dims, 3, DT_FLOAT32);

    /* Re-estimate CSF distance */
    for (i = 0; i < nvox; i++)
    {
        if ((src_copy[i] > CGM) && (src_copy[i] < GWM) && (GMT[i] > 1e-15f))
            dist_CSF[i] = fminf(dist_CSF[i], GMT[i] - dist_WM[i]);
    }
    clip_data(dist_CSF, nvox, 0.0, 1e15, DT_FLOAT32);

    if (verbose)
        fprintf(stderr, "Estimate percentage position map.\n");

    /* PPM estimation */
    for (i = 0; i < nvox; i++)
        PPM[i] = (src_copy[i] >= GWM) ? 1.0f : 0.0f;

    for (i = 0; i < nvox; i++)
    {
        if ((src_copy[i] > CGM) && (src_copy[i] < GWM) && (GMT[i] > 1e-15f))
            PPM[i] = fmaxf(0.0f, GMT[i] - dist_WM[i]) / GMT[i];
    }
    clip_data(PPM, nvox, 0.0, 1.0, DT_FLOAT32);

    /* Fill holes */
    if (fill_thresh > 0.0)
        fill_holes(PPM, dims, fill_thresh, 1.0, DT_FLOAT32);

    /* Use the maximum between the median and the PPM for values above the
       isovalue of 0.5 (which are rather gyral), and the minimum otherwise,
       to strengthen gyri and weaken sulci. */
    if (!opts->fast)
    {
        /* Rescue unfiltered PPM */
        memcpy(src_copy, PPM, sizeof(float) * nvox);

        /* Median-filtering of PPM with use of euclidean distance */
        if (sheet)
            CAT_VolOrientedMedian(PPM, sheet, sheet_nrm, NULL, dims,
                                  opts->oriented_cutoff, 2);
        else
            localstat3(PPM, NULL, dims, 1, F_MEDIAN, 2, 1, DT_FLOAT32);

        /* Use either minimum or maximum of median and PPM w.r.t. threshold of 0.5 */
        for (i = 0; i < nvox; i++)
            PPM[i] = (src_copy[i] > 0.5) ? fmaxf(src_copy[i], PPM[i]) : fminf(src_copy[i], PPM[i]);
    }

    /* Correct PPM and GMT in sulci using CSF proximity to prevent sulcal gluing */
    if (opts->sulcal_width > 0.0)
    {
        if (verbose)
            fprintf(stderr, "Correct PPM and GMT in sulci using CSF proximity (width=%.1f mm).\n",
                    opts->sulcal_width);
        correct_ppm_sulci(src, PPM, GMT, dist_CSF, dist_WM, dims, voxelsize,
                          opts->sulcal_width);
    }

    /* Convert thickness from voxels to mm, then apply the mm-defined
       correction.  Adding the correction after the conversion is what keeps it
       independent of the resolution PBT was run on. */
    for (i = 0; i < nvox; i++)
        GMT[i] = GMT[i] * mean_vx_size + (float)correct_thickness;

    /* Median filter preprocessing for topology artifact reduction
     * A topology-artifact likelihood map is estimated from the positive
     * residual PPM - smooth(PPM), restricted to sufficiently thick cortex
     * (GMT > 1.5), regularized morphologically, then smoothed. This soft weight
     * map blends the original PPM with a locally median-filtered PPM so that only
     * likely topology-artifact regions receive strong filtering.
     */
    if (n_median_filter)
    {
        if (verbose)
            fprintf(stderr, "Local median filtering.\n");
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
            if (sheet)
                CAT_VolOrientedMedian(input, sheet, sheet_nrm, NULL, dims,
                                      opts->oriented_cutoff, n_median_filter);
            else
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
    free(src_copy);
    if (src_val)
        free(src_val);
    if (sheet)
        free(sheet);
    if (sheet_nrm)
        free(sheet_nrm);
    if (medial)
        free(medial);
    if (gyral)
        free(gyral);

    return 0;
}

/**
 * \brief Estimate the width of the partial-volume ramp, in voxels.
 *
 * Port of estimatePVEsize() from CAT12's cat_vol_pbtsimpleCS4.m (the accurate
 * branch). For both tissue transitions it measures how far apart the two pure
 * tissue cores are: the distance from a PVE voxel out to the lower core plus
 * the distance in to the upper core spans the ramp plus the voxel itself, so
 * one is subtracted. Only PVE voxels contribute -- inside either core one of
 * the two distances is zero because the voxel is either in the source set or
 * outside the mask, so the sum is zero and the percentile call skips it.
 *
 * The result is the median over the volume and is floored at one voxel: a ramp
 * cannot be sharper than the sampling.
 *
 * \param src  (in) PVE label image (CSF=1, GM=2, WM=3)
 * \param dims (in) volume dimensions {nx, ny, nz}
 * \return ramp width in voxels, >= 1.0
 */
static double estimate_pve_width(const float *src, int dims[3])
{
    const int nvox = dims[0] * dims[1] * dims[2];
    float *a, *b, *pvet;
    unsigned char *m;
    double threshold[2], prctile[2], width;
    int i;

    a = (float *)malloc(sizeof(float) * nvox);
    b = (float *)malloc(sizeof(float) * nvox);
    pvet = (float *)malloc(sizeof(float) * nvox);
    m = (unsigned char *)malloc(sizeof(unsigned char) * nvox);

    if (!a || !b || !pvet || !m)
    {
        if (a) free(a);
        if (b) free(b);
        if (pvet) free(pvet);
        if (m) free(m);
        return 1.0;
    }

    /* CSF/GM transition */
    for (i = 0; i < nvox; i++)
    {
        a[i] = (src[i] < 1.1f) ? 1.0f : 0.0f;
        m[i] = (src[i] < 1.9f) ? 1 : 0;
    }
    euclidean_distance(a, m, dims, NULL, 0);
    for (i = 0; i < nvox; i++)
    {
        b[i] = (src[i] > 1.9f) ? 1.0f : 0.0f;
        m[i] = (src[i] > 1.1f) ? 1 : 0;
    }
    euclidean_distance(b, m, dims, NULL, 0);
    for (i = 0; i < nvox; i++)
        pvet[i] = a[i] + b[i];

    /* GM/WM transition */
    for (i = 0; i < nvox; i++)
    {
        a[i] = (src[i] < 2.1f) ? 1.0f : 0.0f;
        m[i] = (src[i] < 2.9f) ? 1 : 0;
    }
    euclidean_distance(a, m, dims, NULL, 0);
    for (i = 0; i < nvox; i++)
    {
        b[i] = (src[i] > 2.9f) ? 1.0f : 0.0f;
        m[i] = (src[i] > 2.1f) ? 1 : 0;
    }
    euclidean_distance(b, m, dims, NULL, 0);
    for (i = 0; i < nvox; i++)
        pvet[i] = fmaxf(pvet[i], a[i] + b[i]);

    prctile[0] = 50.0;
    prctile[1] = 50.0;
    get_prctile(pvet, nvox, threshold, prctile, 1, DT_FLOAT32);

    width = threshold[0] - 1.0;
    if (!(width >= 1.0))
        width = 1.0;

    free(a);
    free(b);
    free(pvet);
    free(m);

    return width;
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
            (SEG[i] >= 1.0) && (SEGI > 1.2 && SEGI <= 2.75) &&   /* projection range */
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
    const int NI[14] = {0, -1, -x + 1, -x, -x - 1, -xy + 1, -xy, -xy - 1,
                        -xy + x + 1, -xy + x, -xy + x - 1, -xy - x + 1, -xy - x, -xy - x - 1}; // Neighbour index offsets
    const float ND[14] = {0.0, 1.0, s2, 1.0, s2, s2, 1.0, s2, s3, s2, s3, s3, s2, s3};         // Neighbour distances
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

    /* We use GMT as temporary variable to get ORNLM-filtered SEG */
    for (i = 0; i < nvox; i++)
        GMT[i] = SEG[i];
    ornlm(GMT, 3, 1, 0.005, 0.005, dims);
    /* Finally use minimum of original and filtered output */
    for (i = 0; i < nvox; i++)
        SEG[i] = MIN(GMT[i], SEG[i]);

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
}

/**
 * \brief Correct PPM in sulci using CSF tissue values to prevent sulcal gluing.
 *
 * Segmentation errors or blood vessels can produce elevated PPM values inside
 * sulci. When the central surface is extracted at isovalue 0.5, these elevated
 * values cause opposing sulcal banks to merge ("glued sulci"). This function
 * attenuates PPM where the original PVE label indicates CSF proximity,
 * effectively carving a thin separator through sulcal fundi.
 *
 * Additionally, GMT (cortical thickness) is corrected in the same region.
 * Where CSF reconstruction failed, dist_CSF is overestimated, leading to
 * inflated thickness values. For near-CSF voxels the thickness is capped at
 * 2 * dist_WM (the cortex cannot be thicker than twice the WM distance when
 * the CSF boundary is nearby), blended with the same quadratic weight.
 *
 * For cortical voxels in the GM/CSF partial-volume zone (CGM <= src < GM)
 * that lie close to the CSF boundary (dist_CSF < sulcal_width), PPM is
 * clamped to a tissue-composition-based maximum: max_ppm = (src - 1) / 2,
 * which maps CSF (1.0) -> 0.0 and GM (2.0) -> 0.5. A smooth quadratic
 * falloff prevents hard edges.
 *
 * \param src          (in)     PVE label image (CSF=1, GM=2, WM=3)
 * \param PPM          (in/out) percentage position map to correct
 * \param GMT          (in/out) gray matter thickness to correct (voxel units)
 * \param dist_CSF     (in)     distance-to-CSF boundary (in voxel units)
 * \param dist_WM      (in)     distance-to-WM boundary (in voxel units)
 * \param dims         (in)     volume dimensions {nx, ny, nz}
 * \param voxelsize    (in)     voxel spacing in mm {dx, dy, dz}
 * \param sulcal_width (in)     max distance from CSF to correct (mm)
 */
static void correct_ppm_sulci(const float *src, float *PPM, float *GMT,
                              const float *dist_CSF, const float *dist_WM,
                              int dims[3], double voxelsize[3],
                              double sulcal_width)
{
    const int nvox = dims[0] * dims[1] * dims[2];
    const float mean_vx = (float)(voxelsize[0] + voxelsize[1] + voxelsize[2]) / 3.0f;
    /* Convert sulcal_width from mm to voxel units */
    const float sw_vox = (float)(sulcal_width / mean_vx);
    int i;

    if (sw_vox <= 0.0f)
        return;

    for (i = 0; i < nvox; i++)
    {
        /* Only act on GM/CSF partial-volume zone with non-zero PPM */
        if (src[i] >= CGM && src[i] < GM && PPM[i] > 0.0f)
        {
            float csf_prox = dist_CSF[i]; /* distance to CSF in voxel units */

            if (csf_prox < sw_vox)
            {
                /* Smooth quadratic falloff: strongest near CSF */
                float weight = 1.0f - csf_prox / sw_vox;
                weight *= weight;

                /* Expected max PPM from tissue composition:
                 * CSF (1.0) -> 0.0, GM (2.0) -> 0.5 */
                float max_ppm = (src[i] - 1.0f) * 0.5f;
                float corrected_ppm = fminf(PPM[i], max_ppm);
                PPM[i] = (1.0f - weight) * PPM[i] + weight * corrected_ppm;

                /* Cap thickness: near CSF the cortex cannot be thicker
                 * than approximately 2 * dist_WM */
                if (GMT[i] > 0.0f && dist_WM[i] > 0.0f)
                {
                    float max_gmt = 2.0f * dist_WM[i];
                    float corrected_gmt = fminf(GMT[i], max_gmt);
                    GMT[i] = (1.0f - weight) * GMT[i] + weight * corrected_gmt;
                }
            }
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

/**
 * \brief Blood-vessel correction for PVE label maps.
 *
 * Implements blood-vessel correction for a PVE map in the range [0..3].
 * A safe WM seed region is estimated by thresholding and distance-based opening,
 * then expanded with downcut region growing constrained by transformed intensities.
 * Identified vessel-like WM outliers are inpainted from surrounding non-vessel
 * tissue and finally clamped to class-aware limits estimated from a one-voxel
 * ring around the vessel mask.
 *
 * This implementation uses CAT-Surface library tools only:
 * - `downcut_float()` for constrained region growing
 * - `dist_open_float()` / `dist_dilate_float()` for distance morphology
 * - `get_median_double()` for local inpainting statistics
 *
 * \param Yp0              (in/out) float PVE label image in [0..3]
 * \param dims             (in)     dimensions {nx, ny, nz}
 * \param vx_vol           (in)     voxel spacing {sx, sy, sz}; NULL -> {1,1,1}
 */
void blood_vessel_correction_pve_float(float *Yp0, int dims[3],
                                       double vx_vol[3])
{
    const int nvox = dims[0] * dims[1] * dims[2];
    double vx_local[3] = {1.0, 1.0, 1.0};
    float *F, *YwmA, *YwmB, *Ywm, *Yd, *Yp0s, *Ynn, local_replace;
    unsigned char *Ymsk, *Yring, *Yfill, *nanMsk, *fillMsk, *brainMsk;
    float *Icon;
    int i;

    if (!Yp0 || !dims)
    {
        fprintf(stderr, "Invalid NULL input in blood_vessel_correction_pve_float\n");
        exit(EXIT_FAILURE);
    }

    if (vx_vol)
    {
        vx_local[0] = vx_vol[0];
        vx_local[1] = vx_vol[1];
        vx_local[2] = vx_vol[2];
    }

    F = (float *)malloc((size_t)nvox * sizeof(float));
    YwmA = (float *)malloc((size_t)nvox * sizeof(float));
    YwmB = (float *)malloc((size_t)nvox * sizeof(float));
    Ywm = (float *)malloc((size_t)nvox * sizeof(float));
    Yd = (float *)malloc((size_t)nvox * sizeof(float));
    Yp0s = (float *)malloc((size_t)nvox * sizeof(float));
    Ynn = (float *)malloc((size_t)nvox * sizeof(float));
    Ymsk = (unsigned char *)malloc((size_t)nvox * sizeof(unsigned char));
    Yring = (unsigned char *)malloc((size_t)nvox * sizeof(unsigned char));
    Yfill = (unsigned char *)malloc((size_t)nvox * sizeof(unsigned char));
    nanMsk = (unsigned char *)malloc((size_t)nvox * sizeof(unsigned char));
    fillMsk = (unsigned char *)malloc((size_t)nvox * sizeof(unsigned char));
    brainMsk = (unsigned char *)malloc((size_t)nvox * sizeof(unsigned char));
    Icon = (float *)malloc((size_t)nvox * sizeof(float));

    if (!F || !YwmA || !YwmB || !Ywm || !Yd || !Yp0s || !Ynn ||
        !Ymsk || !Yring || !Yfill || !nanMsk || !fillMsk || !brainMsk || !Icon)
    {
        free(F);
        free(YwmA);
        free(YwmB);
        free(Ywm);
        free(Yd);
        free(Yp0s);
        free(Ynn);
        free(Ymsk);
        free(Yring);
        free(Yfill);
        free(nanMsk);
        free(fillMsk);
        free(brainMsk);
        free(Icon);
        fprintf(stderr, "Memory allocation error in blood_vessel_correction_pve_float\n");
        exit(EXIT_FAILURE);
    }

    /*
     * F = max(0, Yp0 - 1); F(Yp0 <= 1.1) = inf;
     */
    for (i = 0; i < nvox; ++i)
    {
        float fval = Yp0[i] - 1.0f;
        if (fval < 0.0f)
            fval = 0.0f;
        if (Yp0[i] <= 1.1f)
            fval = FLT_MAX;
        F[i] = fval;
    }

    /*
     * Build a brain ROI mask covering GM+WM (Yp0 > 1.0) dilated by 5 mm.
     * This is passed to the distance-morphology functions so they can skip
     * voxels deep in the background/CSF, giving a significant speedup.
     */
    for (i = 0; i < nvox; ++i)
        Ywm[i] = (Yp0[i] > 1.0f) ? 1.0f : 0.0f;
    dist_dilate_float(Ywm, dims, vx_local, 5.0, 0.5, NULL);
    for (i = 0; i < nvox; ++i)
        brainMsk[i] = (Ywm[i] > 0.0f) ? 1 : 0;

    /*
     * Ywm = morph(Yp0>2.5,'ldo',2,vx) | morph(Yp0>2.75,'ldo',1,vx)
     * Approximated with distance-based opening for logical masks.
     */
    for (i = 0; i < nvox; ++i)
    {
        YwmA[i] = (Yp0[i] > 2.5f) ? 1.0f : 0.0f;
        YwmB[i] = (Yp0[i] > 2.75f) ? 1.0f : 0.0f;
    }
    double vx1[3] = {1.0, 1.0, 1.0};
    dist_open_float(YwmA, dims, vx_local, 2.0, 0.5, brainMsk);
    dist_open_float(YwmB, dims, vx1, 1.0, 0.5, brainMsk);
    keep_largest_cluster(YwmA, 0.5, dims, DT_FLOAT32, 0, 1, 18);
    keep_largest_cluster(YwmB, 0.5, dims, DT_FLOAT32, 0, 1, 18);
    for (i = 0; i < nvox; ++i)
        Ywm[i] = (YwmA[i] > 0.0f || YwmB[i] > 0.0f) ? 1.0f : 0.0f;

    /*
     * [~,Yd] = cat_vol_downcut(single(Ywm), F, -0.001)
     */
    downcut_float(Ywm, F, Yd, dims, -0.001, vx_local, NULL);

    /*
     * Ymsk = Yd > 1000000 & Yp0 > 2.1
     */
    float th_wm = 2.1f;
    for (i = 0; i < nvox; ++i)
    {
        Ywm[i] = Yp0[i] > th_wm;
        Ymsk[i] = (Yd[i] > 1000000.0f && Ywm[i]) ? 1 : 0;
        if (Ymsk[i] > 0)
            Ywm[i] = 0;
    }
    morph_close(Ywm, dims, 2, 0.5, DT_FLOAT32);
    for (i = 0; i < nvox; ++i)
    {
        if (Ymsk[i])
            Ymsk[i] = (Ywm[i]) ? 0 : Ymsk[i];
        if (Ymsk[i])
            Yp0[i] = NAN;
    }

    for (i = 0; i < nvox; ++i)
        Ywm[i] = (float)Ymsk[i];
    dist_dilate_float(Ywm, dims, vx_local, 1.0, 0.5, brainMsk);

    /*
     * Ring-constrained inpainting.
     * 1) set vessel mask to NaN
     * 2) iteratively inpaint inward with masked median from touching neighbours
     * 3) fallback nearest-value propagation for remaining NaNs
     * 4) class-aware clamp inside vessel mask using 75th percentile of ring values
     */
    for (i = 0; i < nvox; ++i)
    {
        const int is_nan = isnan(Yp0[i]) ? 1 : 0;
        Yring[i] = (Ywm[i] > 0.0f && !Ymsk[i] && !is_nan) ? 1 : 0;
        Yfill[i] = (!is_nan && !Ymsk[i]) ? 1 : 0;
        nanMsk[i] = (is_nan && Ymsk[i]) ? 1 : 0;
    }

    {
        const int nx = dims[0];
        const int ny = dims[1];
        const int nz = dims[2];
        const int xy = nx * ny;
        const int maxIter = 200;
        int iter = 0;

        while (iter < maxIter)
        {
            int changed = 0;
            int has_nan = 0;
            int z, y, x;

            iter++;
            memset(fillMsk, 0, (size_t)nvox * sizeof(unsigned char));

            for (i = 0; i < nvox; ++i)
            {
                if (nanMsk[i])
                    has_nan = 1;
            }
            if (!has_nan)
                break;

            for (z = 0; z < nz; ++z)
            {
                const int zmin = (z > 0) ? (z - 1) : z;
                const int zmax = (z + 1 < nz) ? (z + 1) : z;

                for (y = 0; y < ny; ++y)
                {
                    const int ymin = (y > 0) ? (y - 1) : y;
                    const int ymax = (y + 1 < ny) ? (y + 1) : y;

                    for (x = 0; x < nx; ++x)
                    {
                        const int idx = x + y * nx + z * xy;
                        int zz, yy, xx;

                        if (!nanMsk[idx])
                            continue;

                        for (zz = zmin; zz <= zmax && !fillMsk[idx]; ++zz)
                        {
                            for (yy = ymin; yy <= ymax && !fillMsk[idx]; ++yy)
                            {
                                const int xmin = (x > 0) ? (x - 1) : x;
                                const int xmax = (x + 1 < nx) ? (x + 1) : x;
                                for (xx = xmin; xx <= xmax; ++xx)
                                {
                                    const int nidx = xx + yy * nx + zz * xy;
                                    if (Yfill[nidx])
                                    {
                                        fillMsk[idx] = 1;
                                        break;
                                    }
                                }
                            }
                        }
                    }
                }
            }

            /*
             * Compute median of non-NaN
             * 26-connected neighbours for each fillMsk voxel.
             */

            for (i = 0; i < nvox; ++i)
            {
                if (!fillMsk[i])
                    continue;
                {
                    const int ix = i % nx;
                    const int iy = (i / nx) % ny;
                    const int iz = i / xy;
                    double vals[27];
                    int cnt = 0;
                    int zz, yy, xx;
                    const int zmin_l = (iz > 0) ? (iz - 1) : iz;
                    const int zmax_l = (iz + 1 < nz) ? (iz + 1) : iz;
                    const int ymin_l = (iy > 0) ? (iy - 1) : iy;
                    const int ymax_l = (iy + 1 < ny) ? (iy + 1) : iy;
                    const int xmin_l = (ix > 0) ? (ix - 1) : ix;
                    const int xmax_l = (ix + 1 < nx) ? (ix + 1) : ix;

                    for (zz = zmin_l; zz <= zmax_l; ++zz)
                        for (yy = ymin_l; yy <= ymax_l; ++yy)
                            for (xx = xmin_l; xx <= xmax_l; ++xx)
                            {
                                const int nidx = xx + yy * nx + zz * xy;
                                if (!isnan(Yp0[nidx]))
                                {
                                    vals[cnt++] = (double)Yp0[nidx];
                                }
                            }

                    if (cnt > 0)
                        Yp0s[i] = (float)get_median_double(vals, cnt, 0);
                    else
                        Yp0s[i] = NAN;
                }
            }

            for (i = 0; i < nvox; ++i)
            {
                if (fillMsk[i] && !isnan(Yp0s[i]))
                {
                    Yp0[i] = Yp0s[i];
                    Yfill[i] = 1;
                    nanMsk[i] = 0;
                    changed = 1;
                }
            }

            if (!changed)
                break;
        }

        /* Fallback for remaining NaNs: nearest value from non-vessel support. */

        for (i = 0; i < nvox; ++i)
        {
            Ynn[i] = (!Ymsk[i] && !isnan(Yp0[i])) ? (Yp0[i] + 1.0f) : 0.0f;
            Icon[i] = 0.0f;
        }

        double dd_nn[2] = {1.0, 0.0};
        downcut_float(Ynn, Icon, NULL, dims, 0.0, vx_local, dd_nn);

        for (i = 0; i < nvox; ++i)
        {
            if (nanMsk[i])
            {
                if (Ynn[i] > 0.0f)
                    Yp0[i] = Ynn[i] - 1.0f;
                else
                    Yp0[i] = 2.0f;
                nanMsk[i] = 0;
            }
        }
    }

    free(F);
    free(YwmA);
    free(YwmB);
    free(Ywm);
    free(Yd);
    free(Yp0s);
    free(Ynn);
    free(Ymsk);
    free(Yring);
    free(Yfill);
    free(nanMsk);
    free(fillMsk);
    free(brainMsk);
    free(Icon);
}

/**
 * \brief Datatype-generic blood-vessel correction wrapper for PVE labels.
 *
 * Converts input PVE labels to float, runs `blood_vessel_correction_pve_float()`, and converts back
 * to the requested datatype.
 *
 * \param data             (in/out) PVE label volume data
 * \param Ygmt             (in)     optional local gray-matter thickness map (NULL -> scalar replacement)
 * \param dims             (in)     dimensions {nx, ny, nz}
 * \param vx_vol           (in)     voxel spacing {sx, sy, sz}; NULL -> {1,1,1}
 * \param datatype         (in)     datatype code (DT_UINT8, DT_UINT16, DT_FLOAT32, etc.)
 */
void blood_vessel_correction_pve(void *data, int dims[3], double vx_vol[3], int datatype)
{
    const int nvox = dims[0] * dims[1] * dims[2];
    float *buffer;

    if (!data || !dims)
    {
        fprintf(stderr, "Invalid NULL input in blood_vessel_correction_pve\n");
        exit(EXIT_FAILURE);
    }

    buffer = (float *)malloc((size_t)nvox * sizeof(float));
    if (!buffer)
    {
        fprintf(stderr, "Memory allocation error in blood_vessel_correction_pve\n");
        exit(EXIT_FAILURE);
    }

    convert_input_type_float(data, buffer, nvox, datatype);
    blood_vessel_correction_pve_float(buffer, dims, vx_vol);
    convert_output_type_float(data, buffer, nvox, datatype);

    free(buffer);
}
