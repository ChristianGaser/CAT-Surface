#include "minunit.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "CAT_Vol.h"

/*
 * morph_dilate_geodesic() grows a mask step by step, confined to a region.
 * It visits only the front, so the reference here is the obvious dense
 * version: one full pass per step.  Both have to agree exactly, including at
 * the volume border and with the alternating 6/26 schedule, which is what the
 * callers (a ventricle fill, a brainstem cut) rely on.
 */

static int DIMS[3] = {19, 17, 23};

/** \brief Dense reference: one full pass per step. */
static void
reference(unsigned char *mask, const unsigned char *region, int dims[3],
          int niter, int alternate)
{
    int nvox = dims[0] * dims[1] * dims[2];
    unsigned char *tmp = (unsigned char *)malloc(nvox);
    int step, x, y, z, dx, dy, dz;

    for (step = 0; step < niter; step++)
    {
        memcpy(tmp, mask, nvox);
        for (x = 0; x < dims[0]; x++)
            for (y = 0; y < dims[1]; y++)
                for (z = 0; z < dims[2]; z++)
                {
                    int idx = z * dims[0] * dims[1] + y * dims[0] + x;
                    if (tmp[idx])
                        continue;
                    if (region && !region[idx])
                        continue;
                    for (dx = -1; dx <= 1; dx++)
                        for (dy = -1; dy <= 1; dy++)
                            for (dz = -1; dz <= 1; dz++)
                            {
                                int man = abs(dx) + abs(dy) + abs(dz);
                                int X = x + dx, Y = y + dy, Z = z + dz;
                                if (man == 0)
                                    continue;
                                if (alternate && (step % 2 == 0) && man != 1)
                                    continue;
                                if (X < 0 || Y < 0 || Z < 0 || X >= dims[0] ||
                                    Y >= dims[1] || Z >= dims[2])
                                    continue;
                                if (tmp[Z * dims[0] * dims[1] + Y * dims[0] + X])
                                    mask[idx] = 1;
                            }
                }
    }
    free(tmp);
}

/** \brief Random mask/region pair; every voxel of the border is foreground. */
static void
fill_random(unsigned char *mask, unsigned char *region, int nvox, int seed,
            int touch_border, int dims[3])
{
    int i, x, y;
    srand(seed);
    for (i = 0; i < nvox; i++)
    {
        mask[i] = (rand() % 1000) < 4;
        region[i] = (rand() % 100) < 85;
    }
    if (touch_border)
        for (x = 0; x < dims[0]; x++)
            for (y = 0; y < dims[1]; y++)
                mask[y * dims[0] + x] = 1;      /* the whole z = 0 plane */
}

static void
test_matches_the_dense_version(void)
{
    int nvox = DIMS[0] * DIMS[1] * DIMS[2];
    unsigned char *mask = (unsigned char *)malloc(nvox);
    unsigned char *region = (unsigned char *)malloc(nvox);
    unsigned char *want = (unsigned char *)malloc(nvox);
    int iters[5] = {1, 2, 3, 5, 10};
    int seed, k, alt, border, mismatch = 0;

    for (seed = 1; seed <= 3; seed++)
        for (border = 0; border <= 1; border++)
            for (alt = 0; alt <= 1; alt++)
                for (k = 0; k < 5; k++)
                {
                    fill_random(mask, region, nvox, seed, border, DIMS);
                    memcpy(want, mask, nvox);
                    reference(want, region, DIMS, iters[k], alt);
                    morph_dilate_geodesic(mask, region, DIMS, iters[k], alt);
                    if (memcmp(mask, want, nvox) != 0)
                        mismatch++;
                }
    MU_ASSERT("geodesic dilation matches the dense version", mismatch == 0);
    free(mask);
    free(region);
    free(want);
}

static void
test_region_null_grows_everywhere(void)
{
    int nvox = DIMS[0] * DIMS[1] * DIMS[2];
    unsigned char *mask = (unsigned char *)malloc(nvox);
    unsigned char *want = (unsigned char *)malloc(nvox);
    unsigned char *region = (unsigned char *)malloc(nvox);

    fill_random(mask, region, nvox, 7, 0, DIMS);
    memcpy(want, mask, nvox);
    reference(want, NULL, DIMS, 3, 1);
    morph_dilate_geodesic(mask, NULL, DIMS, 3, 1);
    MU_ASSERT("without a region it is a plain dilation",
              memcmp(mask, want, nvox) == 0);
    free(mask);
    free(want);
    free(region);
}

static void
test_edge_cases(void)
{
    int nvox = DIMS[0] * DIMS[1] * DIMS[2];
    unsigned char *mask = (unsigned char *)calloc(nvox, 1);
    unsigned char *region = (unsigned char *)malloc(nvox);
    int i, sum;

    memset(region, 1, nvox);
    morph_dilate_geodesic(mask, region, DIMS, 5, 1);
    for (i = 0, sum = 0; i < nvox; i++)
        sum += mask[i];
    MU_ASSERT("an empty mask stays empty", sum == 0);

    mask[nvox / 2] = 3;                       /* any non-zero is foreground */
    morph_dilate_geodesic(mask, region, DIMS, 0, 1);
    MU_ASSERT("niter 0 is a no-op", mask[nvox / 2] == 3);

    morph_dilate_geodesic(mask, region, DIMS, 1, 0);
    MU_ASSERT("a non-zero seed is normalised to 1", mask[nvox / 2] == 1);

    memset(mask, 0, nvox);
    mask[nvox / 2] = 1;
    memset(region, 0, nvox);
    morph_dilate_geodesic(mask, region, DIMS, 4, 1);
    for (i = 0, sum = 0; i < nvox; i++)
        sum += mask[i];
    MU_ASSERT("an empty region blocks all growth", sum == 1);
    free(mask);
    free(region);
}

int main(void)
{
    MU_RUN_TEST(test_matches_the_dense_version);
    MU_RUN_TEST(test_region_null_grows_everywhere);
    MU_RUN_TEST(test_edge_cases);
    printf("%d tests run, %d failed\n", tests_run, tests_failed);
    return tests_failed ? 1 : 0;
}
