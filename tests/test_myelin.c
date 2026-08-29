#include "minunit.h"
#include <stdio.h>
#include "CAT_VolMyelinCorrection.h"
#include "CAT_Vol.h"

#include <stdlib.h>
#include <string.h>
#include <math.h>

#define N 48
#define NVOX (N * N * N)
#define IDX(x, y, z) ((x) + (y) * N + (z) * N * N)

static int dims3[3] = {N, N, N};
static double vx1[3] = {1.0, 1.0, 1.0};

/* ------------------------------------------------------------------ */
/* the cluster filter                                                  */
/* ------------------------------------------------------------------ */

/* Two cubes of different size, far enough apart to stay disconnected at
 * 18-connectivity. */
static void
two_cubes(float *vol, int small_side)
{
    int x, y, z;

    memset(vol, 0, sizeof(float) * NVOX);

    for (z = 8; z < 14; z++) /* 6^3 = 216 voxels */
        for (y = 8; y < 14; y++)
            for (x = 8; x < 14; x++)
                vol[IDX(x, y, z)] = 1.0f;

    for (z = 30; z < 30 + small_side; z++)
        for (y = 30; y < 30 + small_side; y++)
            for (x = 30; x < 30 + small_side; x++)
                vol[IDX(x, y, z)] = 1.0f;
}

static int
count_nonzero(const float *vol)
{
    int i, n = 0;
    for (i = 0; i < NVOX; i++)
        if (vol[i] != 0.0f)
            n++;
    return n;
}

static void test_cluster_keep_largest(void)
{
    float *vol = (float *)malloc(sizeof(float) * NVOX);

    /* min_size >= 0 keeps the largest cluster only -- the behaviour every
     * other caller in the tree relies on. */
    two_cubes(vol, 3);
    keep_largest_cluster(vol, 0.5, dims3, DT_FLOAT32, 0, 1, 18);
    MU_ASSERT("min_size 0 keeps only the largest cluster",
              count_nonzero(vol) == 216);

    /* and it stays keep-largest for any non-negative magnitude */
    two_cubes(vol, 3);
    keep_largest_cluster(vol, 0.5, dims3, DT_FLOAT32, 5, 1, 18);
    MU_ASSERT("a positive min_size does not become a size floor",
              count_nonzero(vol) == 216);

    free(vol);
}

static void test_cluster_size_floor(void)
{
    float *vol = (float *)malloc(sizeof(float) * NVOX);

    /* A negative min_size keeps every cluster of at least |min_size| voxels,
     * so both cubes survive a floor below the smaller one. */
    two_cubes(vol, 3); /* 216 and 27 voxels */
    keep_largest_cluster(vol, 0.5, dims3, DT_FLOAT32, -10, 1, 18);
    MU_ASSERT("a floor below both sizes keeps both clusters",
              count_nonzero(vol) == 216 + 27);

    /* a floor between the two sizes drops the smaller one only */
    two_cubes(vol, 3);
    keep_largest_cluster(vol, 0.5, dims3, DT_FLOAT32, -100, 1, 18);
    MU_ASSERT("a floor between the two sizes drops the smaller cluster",
              count_nonzero(vol) == 216);

    /* the floor is inclusive: a cluster of exactly |min_size| survives */
    two_cubes(vol, 3);
    keep_largest_cluster(vol, 0.5, dims3, DT_FLOAT32, -27, 1, 18);
    MU_ASSERT("a cluster of exactly |min_size| voxels survives",
              count_nonzero(vol) == 216 + 27);

    /* one above and it does not */
    two_cubes(vol, 3);
    keep_largest_cluster(vol, 0.5, dims3, DT_FLOAT32, -28, 1, 18);
    MU_ASSERT("a cluster one voxel below the floor is dropped",
              count_nonzero(vol) == 216);

    free(vol);
}

/* ------------------------------------------------------------------ */
/* myelination phantom                                                 */
/* ------------------------------------------------------------------ */

/* A layered slab along x inside a margin, with two disjoint patches of
 * myelinated cortex: deep GM that carries WM-like labels but an intensity
 * between GM and WM.  This is the M1 / occipital situation -- two separate
 * regions, neither of them the largest connected structure in the brain. */

/* Deliberately different radii, so that "keep the largest cluster" and "keep
 * every cluster above a size floor" cannot agree by accident. */
#define PATCH_A_R 5
#define PATCH_B_R 3

static int
in_patch(int x, int y, int z, int cy, int cz, int r)
{
    int dy = y - cy, dz = z - cz;
    return (x >= 22 && x < 26 && dy * dy + dz * dz <= r * r);
}

static int
in_patch_a(int x, int y, int z)
{
    return in_patch(x, y, z, 15, 15, PATCH_A_R);
}

static int
in_patch_b(int x, int y, int z)
{
    return in_patch(x, y, z, 33, 33, PATCH_B_R);
}

static void
build_phantom(float *pve, float *t1w)
{
    int x, y, z;

    memset(pve, 0, sizeof(float) * NVOX);
    memset(t1w, 0, sizeof(float) * NVOX);

    for (z = 8; z < 40; z++)
        for (y = 8; y < 40; y++)
            for (x = 0; x < N; x++)
            {
                int i = IDX(x, y, z);

                if (x >= 10 && x < 26)
                {
                    pve[i] = 3.0f; /* WM */
                    t1w[i] = 110.0f;
                }
                else if (x >= 26 && x < 32)
                {
                    pve[i] = 2.0f; /* GM */
                    t1w[i] = 60.0f;
                }
                else if (x >= 32 && x < 36)
                {
                    pve[i] = 1.0f; /* CSF */
                    t1w[i] = 20.0f;
                }

                /* Myelinated cortex mislabelled as WM: still labelled 3.0,
                 * but too dark to be white matter. */
                if (in_patch_a(x, y, z) || in_patch_b(x, y, z))
                    t1w[i] = 85.0f;
            }
}

static int
count_corrected_in_patch(const float *corr, int (*inside)(int, int, int))
{
    int x, y, z, n = 0;

    for (z = 8; z < 40; z++)
        for (y = 8; y < 40; y++)
            for (x = 0; x < N; x++)
                if (inside(x, y, z) && corr[IDX(x, y, z)] != 0.0f)
                    n++;
    return n;
}

static void
phantom_opts(CAT_MyelinCorrOptions *opts)
{
    CAT_MyelinCorrOptionsInit(opts);
    opts->correct_csf = 0;
    /* The phantom is noise-free, so its band gradients are 0 almost
     * everywhere and any percentile of them is 0 too.  Take the gradient
     * criterion out of the way; this test is about which flagged voxels
     * survive, not about how they are flagged. */
    opts->grad_percentile = 100.0;
    opts->max_gm_grad_dist = 0.0;
    opts->n_median_filter = 0;
}

/* Both patches have to survive the cluster filter.  A keep-largest filter
 * passes only one of them, which is what used to happen here. */
static void test_myelin_keeps_every_patch(void)
{
    float *pve = (float *)malloc(sizeof(float) * NVOX);
    float *t1w = (float *)malloc(sizeof(float) * NVOX);
    float *corr = (float *)malloc(sizeof(float) * NVOX);
    CAT_MyelinCorrOptions opts;
    int n_a, n_b;

    build_phantom(pve, t1w);
    phantom_opts(&opts);

    MU_ASSERT("correction runs",
              CAT_VolCorrectMyelination(pve, t1w, dims3, vx1, corr, &opts) == 0);

    n_a = count_corrected_in_patch(corr, in_patch_a);
    n_b = count_corrected_in_patch(corr, in_patch_b);

    MU_ASSERT("the larger myelinated patch is corrected", n_a > 0);
    MU_ASSERT("the smaller myelinated patch is corrected too", n_b > 0);

    free(pve);
    free(t1w);
    free(corr);
}

/* The exported field is the shift that was really applied, and it is
 * one-sided on the WM boundary: labels may only move toward GM. */
static void test_myelin_correction_field(void)
{
    float *pve = (float *)malloc(sizeof(float) * NVOX);
    float *pve0 = (float *)malloc(sizeof(float) * NVOX);
    float *t1w = (float *)malloc(sizeof(float) * NVOX);
    float *corr = (float *)malloc(sizeof(float) * NVOX);
    CAT_MyelinCorrOptions opts;
    int i, n_nonzero = 0, exact = 1, one_sided = 1;

    build_phantom(pve, t1w);
    memcpy(pve0, pve, sizeof(float) * NVOX);
    phantom_opts(&opts);

    CAT_VolCorrectMyelination(pve, t1w, dims3, vx1, corr, &opts);

    for (i = 0; i < NVOX; i++)
    {
        if (corr[i] != pve[i] - pve0[i])
            exact = 0;
        if (corr[i] > 0.0f)
            one_sided = 0;
        if (corr[i] != 0.0f)
            n_nonzero++;
    }

    MU_ASSERT("something was corrected", n_nonzero > 0);
    MU_ASSERT("the field equals the applied shift", exact);
    MU_ASSERT("the WM correction only moves labels toward GM", one_sided);

    free(pve);
    free(pve0);
    free(t1w);
    free(corr);
}

/* Healthy cortex must come back untouched: where the intensities agree with
 * the labels there is nothing to correct, and the field is all zero. */
static void test_myelin_null_case(void)
{
    float *pve = (float *)malloc(sizeof(float) * NVOX);
    float *pve0 = (float *)malloc(sizeof(float) * NVOX);
    float *t1w = (float *)malloc(sizeof(float) * NVOX);
    float *corr = (float *)malloc(sizeof(float) * NVOX);
    CAT_MyelinCorrOptions opts;
    int x, y, z, i, identical = 1, all_zero = 1;

    build_phantom(pve, t1w);

    /* undo the two myelinated patches */
    for (z = 8; z < 40; z++)
        for (y = 8; y < 40; y++)
            for (x = 10; x < 26; x++)
                t1w[IDX(x, y, z)] = 110.0f;

    memcpy(pve0, pve, sizeof(float) * NVOX);
    phantom_opts(&opts);

    CAT_VolCorrectMyelination(pve, t1w, dims3, vx1, corr, &opts);

    for (i = 0; i < NVOX; i++)
    {
        if (pve[i] != pve0[i])
            identical = 0;
        if (corr[i] != 0.0f)
            all_zero = 0;
    }

    MU_ASSERT("healthy labels come back bit-identical", identical);
    MU_ASSERT("the correction field is empty", all_zero);

    free(pve);
    free(pve0);
    free(t1w);
    free(corr);
}

int main(void)
{
    MU_RUN_TEST(test_cluster_keep_largest);
    MU_RUN_TEST(test_cluster_size_floor);
    MU_RUN_TEST(test_myelin_keeps_every_patch);
    MU_RUN_TEST(test_myelin_correction_field);
    MU_RUN_TEST(test_myelin_null_case);
    printf("%d tests run, %d failed\n", tests_run, tests_failed);
    return tests_failed ? 1 : 0;
}
