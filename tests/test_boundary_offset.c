#include "minunit.h"
#include <stdio.h>
#include "CAT_VolBoundaryOffset.h"
#include "CAT_Vol.h"

#include <stdlib.h>
#include <string.h>
#include <math.h>

#define N 64
#define NVOX (N * N * N)
#define IDX(x, y, z) ((x) + (y) * N + (z) * N * N)

static int dims3[3] = {N, N, N};
static double vx1[3] = {1.0, 1.0, 1.0};

/* ------------------------------------------------------------------ */
/* the local tissue reference                                          */
/* ------------------------------------------------------------------ */

static void test_reference_recovers_a_constant(void)
{
    float *val = (float *)malloc(sizeof(float) * NVOX);
    float *ref = (float *)malloc(sizeof(float) * NVOX);
    unsigned char *ctrl = (unsigned char *)calloc(NVOX, 1);
    int x, y, z, ok = 1;

    for (x = 0; x < NVOX; x++)
        val[x] = 100.0f;
    for (z = 8; z < 56; z++)
        for (y = 8; y < 56; y++)
            for (x = 8; x < 56; x++)
                ctrl[IDX(x, y, z)] = 1;

    MU_ASSERT("reference runs",
              CAT_VolLocalTissueReference(val, ctrl, dims3, vx1, 10.0, 0.0,
                                          ref) == 0);

    for (z = 16; z < 48; z++)
        for (y = 16; y < 48; y++)
            for (x = 16; x < 48; x++)
                if (fabsf(ref[IDX(x, y, z)] - 100.0f) > 0.5f)
                    ok = 0;
    MU_ASSERT("a constant field comes back as itself", ok);

    free(val);
    free(ref);
    free(ctrl);
}

/* The point of the construction: it follows a bias field rather than
 * averaging it away, which is what lets a fixed threshold be used on the
 * ratio of two such references. */
static void test_reference_tracks_a_ramp(void)
{
    float *val = (float *)malloc(sizeof(float) * NVOX);
    float *ref = (float *)malloc(sizeof(float) * NVOX);
    unsigned char *ctrl = (unsigned char *)calloc(NVOX, 1);
    int x, y, z;
    double err_local, err_global, mean = 100.0;

    for (z = 0; z < N; z++)
        for (y = 0; y < N; y++)
            for (x = 0; x < N; x++)
                val[IDX(x, y, z)] = 100.0f + 0.5f * (float)x;

    for (z = 8; z < 56; z++)
        for (y = 8; y < 56; y++)
            for (x = 8; x < 56; x++)
                ctrl[IDX(x, y, z)] = 1;

    mean = get_masked_mean_array((void *)val, NVOX, ctrl, DT_FLOAT32);
    CAT_VolLocalTissueReference(val, ctrl, dims3, vx1, 10.0, mean, ref);

    err_local = err_global = 0.0;
    for (z = 20; z < 44; z++)
        for (y = 20; y < 44; y++)
            for (x = 20; x < 44; x++)
            {
                double truth = val[IDX(x, y, z)];
                err_local += fabs(ref[IDX(x, y, z)] - truth);
                err_global += fabs(mean - truth);
            }

    MU_ASSERT("the local reference follows the ramp far better than the "
              "global mean does", err_local < 0.2 * err_global);

    free(val);
    free(ref);
    free(ctrl);
}

/* ------------------------------------------------------------------ */
/* the boundary offset                                                 */
/* ------------------------------------------------------------------ */

/* A layered slab whose GM/WM transition is a sharp step everywhere except
 * inside a disc, where it is spread into a ramp -- the signature of a
 * myelinated deep cortical layer.  The label map is derived from the
 * intensity, as a classifier would derive it, so the label boundary sits at
 * the intensity midpoint in both regions and only the *width* differs. */
static void
build_slab(float *label, float *t1w, int myelinated, double ramp_mm)
{
    int x, y, z, i;

    memset(label, 0, sizeof(float) * NVOX);
    memset(t1w, 0, sizeof(float) * NVOX);

    for (z = 8; z < 56; z++)
        for (y = 8; y < 56; y++)
            for (x = 0; x < N; x++)
            {
                double t1, half, d;
                int dy = y - 32, dz = z - 32;
                int inside = myelinated && (dy * dy + dz * dz <= 100);

                i = IDX(x, y, z);

                /* Intensity: WM 110 inside, GM 60 in the ribbon, CSF 20.
                 * The GM/WM transition is centred on x = 39.5 and is either
                 * one voxel wide or `ramp_mm` wide. */
                half = inside ? (ramp_mm / 2.0) : 0.5;
                d = (39.5 - x) / half; /* +1 deep in WM, -1 out in GM */
                if (d > 1.0)
                    d = 1.0;
                if (d < -1.0)
                    d = -1.0;

                if (x < 46)
                    t1 = 85.0 + 25.0 * d; /* 110 .. 60 */
                else if (x < 50)
                    t1 = 20.0;
                else
                    t1 = 0.0;

                if (x >= 50)
                    continue;

                t1w[i] = (float)t1;

                /* PVE label follows the intensity, as a classifier does. */
                if (x < 46)
                {
                    double lab = 2.0 + (t1 - 60.0) / 50.0;
                    if (lab < 2.0)
                        lab = 2.0;
                    if (lab > 3.0)
                        lab = 3.0;
                    label[i] = (float)lab;
                }
                else
                    label[i] = 1.0f;
            }
}

/* The patch here is ~14% of this phantom's boundary, where heavily myelinated
 * cortex is nearer a tenth of a real hemisphere's, so the default percentile
 * would clip it.  These tests are about the measurement, not about what share
 * of a brain is myelinated, so they open the threshold to match the phantom. */
static void
phantom_opts(CAT_BoundaryOffsetOpts *opts)
{
    CAT_BoundaryOffsetOptionsInit(opts);
    opts->width_pct = 80.0;
}

static void
offset_stats(const float *offset, double *inside_mean, double *outside_mean)
{
    int x, y, z, i;
    double si = 0.0, so = 0.0;
    long ni = 0, no = 0;

    for (z = 8; z < 56; z++)
        for (y = 8; y < 56; y++)
            for (x = 0; x < N; x++)
            {
                int dy = y - 32, dz = z - 32;
                i = IDX(x, y, z);
                if (offset[i] == 0.0f)
                    continue;
                if (dy * dy + dz * dz <= 100)
                {
                    si += offset[i];
                    ni++;
                }
                else
                {
                    so += offset[i];
                    no++;
                }
            }
    *inside_mean = ni ? si / ni : 0.0;
    *outside_mean = no ? so / no : 0.0;
}

/* A slab with no myelination anywhere has one transition width, so the
 * excess over the derived reference is zero and nothing is displaced. */
static void test_uniform_boundary_gives_no_offset(void)
{
    float *label = (float *)malloc(sizeof(float) * NVOX);
    float *t1w = (float *)malloc(sizeof(float) * NVOX);
    float *offset = (float *)malloc(sizeof(float) * NVOX);
    int i, n_nonzero = 0;

    build_slab(label, t1w, 0, 0.0);
    MU_ASSERT("offset runs",
              CAT_VolBoundaryOffset(label, t1w, offset, NULL, dims3, vx1,
                                    NULL) > 0);

    for (i = 0; i < NVOX; i++)
        if (offset[i] > 0.05f)
            n_nonzero++;

    MU_ASSERT("a uniformly sharp boundary is not displaced", n_nonzero == 0);

    free(label);
    free(t1w);
    free(offset);
}

/* A widened transition inside the disc is picked out, and the cortex around
 * it is left alone. */
static void test_a_widened_transition_is_found(void)
{
    float *label = (float *)malloc(sizeof(float) * NVOX);
    float *t1w = (float *)malloc(sizeof(float) * NVOX);
    float *offset = (float *)malloc(sizeof(float) * NVOX);
    float *width = (float *)malloc(sizeof(float) * NVOX);
    CAT_BoundaryOffsetOpts opts;
    double in_mean, out_mean;

    phantom_opts(&opts);
    build_slab(label, t1w, 1, 3.0);
    MU_ASSERT("offset runs",
              CAT_VolBoundaryOffset(label, t1w, offset, width, dims3, vx1,
                                    &opts) > 0);

    offset_stats(offset, &in_mean, &out_mean);

    MU_ASSERT("the myelinated patch is displaced", in_mean > 0.2);
    MU_ASSERT("the healthy cortex around it is not",
              out_mean < 0.25 * in_mean);

    free(label);
    free(t1w);
    free(offset);
    free(width);
}

/* Myelination can only place the boundary too far out, so the correction is
 * one-sided by construction. */
static void test_offset_is_one_sided_and_bounded(void)
{
    float *label = (float *)malloc(sizeof(float) * NVOX);
    float *t1w = (float *)malloc(sizeof(float) * NVOX);
    float *offset = (float *)malloc(sizeof(float) * NVOX);
    CAT_BoundaryOffsetOpts opts;
    int i, ok = 1;

    phantom_opts(&opts);
    build_slab(label, t1w, 1, 3.0);
    CAT_VolBoundaryOffset(label, t1w, offset, NULL, dims3, vx1, &opts);

    for (i = 0; i < NVOX; i++)
        if (offset[i] < 0.0f || offset[i] > (float)opts.max_offset_mm)
            ok = 0;

    MU_ASSERT("the displacement never goes negative or past the clamp", ok);

    free(label);
    free(t1w);
    free(offset);
}

/* A wider ramp means a larger displacement -- the estimate is a measurement
 * in millimetres, not a score. */
static void test_offset_grows_with_the_ramp(void)
{
    float *label = (float *)malloc(sizeof(float) * NVOX);
    float *t1w = (float *)malloc(sizeof(float) * NVOX);
    float *offset = (float *)malloc(sizeof(float) * NVOX);
    CAT_BoundaryOffsetOpts opts;
    double narrow_in, wide_in, dummy;

    phantom_opts(&opts);

    build_slab(label, t1w, 1, 2.0);
    CAT_VolBoundaryOffset(label, t1w, offset, NULL, dims3, vx1, &opts);
    offset_stats(offset, &narrow_in, &dummy);

    build_slab(label, t1w, 1, 4.0);
    CAT_VolBoundaryOffset(label, t1w, offset, NULL, dims3, vx1, &opts);
    offset_stats(offset, &wide_in, &dummy);

    MU_ASSERT("a wider transition yields a larger displacement",
              wide_in > narrow_in + 0.1);

    free(label);
    free(t1w);
    free(offset);
}

/* With noise the measured widths spread out, and a *central* percentile then
 * leaves the healthy population straddling the reference: every voxel above it
 * acquires a positive offset and the whole cortex ends up displaced.  On a real
 * 0.75 mm subject a reference at p25 displaced 99.8% of the boundary.  Reading
 * the threshold as an upper percentile instead confines the correction to the
 * tail, and fixes the share of the boundary it can touch. */
static void test_an_upper_percentile_confines_the_correction(void)
{
    float *label = (float *)malloc(sizeof(float) * NVOX);
    float *t1w = (float *)malloc(sizeof(float) * NVOX);
    float *offset = (float *)malloc(sizeof(float) * NVOX);
    CAT_BoundaryOffsetOpts opts;
    unsigned int rng = 7u;
    long n_k0 = 0, n_dk0 = 0, n_def = 0, n_ddef = 0;
    int x, y, z, i, k;
    double spec_k0, spec_def;

    build_slab(label, t1w, 1, 3.0);
    for (i = 0; i < NVOX; i++)
        if (label[i] > 0.0f)
        {
            double u = 0.0;
            for (k = 0; k < 12; k++)
            {
                rng = rng * 1103515245u + 12345u;
                u += (double)((rng >> 16) & 0x7fff) / 32767.0;
            }
            t1w[i] += (float)(2.0 * (u - 6.0));
        }

    CAT_BoundaryOffsetOptionsInit(&opts);
    opts.width_pct = 25.0; /* a central percentile, as first shipped */
    CAT_VolBoundaryOffset(label, t1w, offset, NULL, dims3, vx1, &opts);
    for (z = 8; z < 56; z++)
        for (y = 8; y < 56; y++)
            for (x = 0; x < N; x++)
            {
                int dy = y - 32, dz = z - 32;
                i = IDX(x, y, z);
                if (offset[i] <= 0.0f)
                    continue;
                n_k0++;
                if (dy * dy + dz * dz <= 100)
                    n_dk0++;
            }

    phantom_opts(&opts);
    CAT_VolBoundaryOffset(label, t1w, offset, NULL, dims3, vx1, &opts);
    for (z = 8; z < 56; z++)
        for (y = 8; y < 56; y++)
            for (x = 0; x < N; x++)
            {
                int dy = y - 32, dz = z - 32;
                i = IDX(x, y, z);
                if (offset[i] <= 0.0f)
                    continue;
                n_def++;
                if (dy * dy + dz * dz <= 100)
                    n_ddef++;
            }

    MU_ASSERT("both settings correct something", n_k0 > 0 && n_def > 0);
    MU_ASSERT("an upper percentile corrects strictly less", n_def < n_k0);

    /* What it buys is specificity: a larger share of what is corrected is
     * inside the patch that actually has a widened transition. */
    spec_k0 = (double)n_dk0 / (double)n_k0;
    spec_def = (double)n_ddef / (double)n_def;
    MU_ASSERT("and a larger share of it is where the widening is",
              spec_def > spec_k0 * 1.2);

    free(label);
    free(t1w);
    free(offset);
}

int main(void)
{
    MU_RUN_TEST(test_reference_recovers_a_constant);
    MU_RUN_TEST(test_reference_tracks_a_ramp);
    MU_RUN_TEST(test_uniform_boundary_gives_no_offset);
    MU_RUN_TEST(test_a_widened_transition_is_found);
    MU_RUN_TEST(test_offset_is_one_sided_and_bounded);
    MU_RUN_TEST(test_offset_grows_with_the_ramp);
    MU_RUN_TEST(test_an_upper_percentile_confines_the_correction);
    printf("%d tests run, %d failed\n", tests_run, tests_failed);
    return tests_failed ? 1 : 0;
}
