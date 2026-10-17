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

/*
 * morph_box() is a dilation / erosion / closing / opening with a cube of side
 * 2r+1, done by sliding-window counts.  The reference is the definition: look
 * at every voxel of the cube; voxels outside the volume are background.  That
 * is what scipy.ndimage's maximum_filter / minimum_filter do with cval=0, and
 * what the T1Prep callers rely on, so the two have to agree exactly.
 */

/** \brief The definition: any (dilate) or all (erode) of the (2r+1)^3 cube. */
static void
box_reference(const unsigned char *src, unsigned char *dst, int dims[3], int r,
              int erode)
{
    int x, y, z, dx, dy, dz;

    for (z = 0; z < dims[2]; z++)
        for (y = 0; y < dims[1]; y++)
            for (x = 0; x < dims[0]; x++)
            {
                int any = 0, all = 1;
                for (dz = -r; dz <= r; dz++)
                    for (dy = -r; dy <= r; dy++)
                        for (dx = -r; dx <= r; dx++)
                        {
                            int X = x + dx, Y = y + dy, Z = z + dz;
                            int v = X >= 0 && Y >= 0 && Z >= 0 && X < dims[0] &&
                                    Y < dims[1] && Z < dims[2] &&
                                    src[(size_t)Z * dims[0] * dims[1] +
                                        (size_t)Y * dims[0] + X] != 0;
                            any |= v;
                            all &= v;
                        }
                dst[(size_t)z * dims[0] * dims[1] + (size_t)y * dims[0] + x] =
                    (unsigned char)(erode ? all : any);
            }
}

/** \brief The reference of each of the four operations. */
static void
box_reference_op(const unsigned char *src, unsigned char *dst, int dims[3], int r,
                 int op)
{
    size_t nvox = (size_t)dims[0] * dims[1] * dims[2];
    unsigned char *tmp = (unsigned char *)malloc(nvox);

    switch (op)
    {
    case MORPH_BOX_DILATE: box_reference(src, dst, dims, r, 0); break;
    case MORPH_BOX_ERODE:  box_reference(src, dst, dims, r, 1); break;
    case MORPH_BOX_CLOSE:
        box_reference(src, tmp, dims, r, 0);
        box_reference(tmp, dst, dims, r, 1);
        break;
    default:
        box_reference(src, tmp, dims, r, 1);
        box_reference(tmp, dst, dims, r, 0);
        break;
    }
    free(tmp);
}

/** \brief Random mask with the given foreground density in per mille. */
static void
fill_density(unsigned char *mask, size_t nvox, int seed, int per_mille)
{
    size_t i;
    srand(seed);
    for (i = 0; i < nvox; i++)
        mask[i] = (rand() % 1000) < per_mille;
}

static void
test_box_matches_the_definition(void)
{
    /* thin and smaller-than-the-cube shapes included: the border logic is
       where a sliding window goes wrong */
    static const int shapes[6][3] = {{19, 17, 23}, {1, 9, 6}, {3, 3, 3},
                                     {40, 2, 11}, {7, 31, 1}, {2, 1, 1}};
    static const int radii[6] = {1, 2, 3, 5, 9, 40};
    static const int density[3] = {5, 400, 970}; /* sparse, half, almost full */
    int s, k, d, op, mismatch = 0, runs = 0;

    for (s = 0; s < 6; s++)
    {
        int dims[3] = {shapes[s][0], shapes[s][1], shapes[s][2]};
        size_t nvox = (size_t)dims[0] * dims[1] * dims[2];
        unsigned char *mask = (unsigned char *)malloc(nvox);
        unsigned char *want = (unsigned char *)malloc(nvox);
        unsigned char *got = (unsigned char *)malloc(nvox);

        for (k = 0; k < 6; k++)
        {
            double cube = (2.0 * radii[k] + 1) * (2.0 * radii[k] + 1) * (2.0 * radii[k] + 1);
            if ((double)nvox * cube > 3e7)
                continue; /* keep the dense reference quick */
            for (d = 0; d < 3; d++)
                for (op = MORPH_BOX_DILATE; op <= MORPH_BOX_OPEN; op++)
                {
                    fill_density(mask, nvox, 11 * s + 3 * k + d + 1, density[d]);
                    box_reference_op(mask, want, dims, radii[k], op);
                    memset(got, 7, nvox);
                    morph_box(mask, got, dims, radii[k], op);
                    runs++;
                    if (memcmp(got, want, nvox) != 0)
                        mismatch++;
                }
        }
        free(mask);
        free(want);
        free(got);
    }
    MU_ASSERT("box morphology matches the definition", mismatch == 0 && runs > 100);
}

static void
test_box_in_place(void)
{
    int dims[3] = {19, 17, 23};
    size_t nvox = (size_t)dims[0] * dims[1] * dims[2];
    unsigned char *mask = (unsigned char *)malloc(nvox);
    unsigned char *want = (unsigned char *)malloc(nvox);
    unsigned char *work = (unsigned char *)malloc(nvox);
    int op, r, mismatch = 0;

    for (op = MORPH_BOX_DILATE; op <= MORPH_BOX_OPEN; op++)
        for (r = 1; r <= 4; r++)
        {
            fill_density(mask, nvox, 5 * op + r, 450);
            morph_box(mask, want, dims, r, op);
            memcpy(work, mask, nvox);
            morph_box(work, work, dims, r, op); /* src is dst */
            if (memcmp(work, want, nvox) != 0)
                mismatch++;
        }
    MU_ASSERT("dst may be src", mismatch == 0);
    free(mask);
    free(want);
    free(work);
}

static void
test_box_edge_cases(void)
{
    int dims[3] = {19, 17, 23};
    size_t nvox = (size_t)dims[0] * dims[1] * dims[2], i, sum;
    unsigned char *src = (unsigned char *)calloc(nvox, 1);
    unsigned char *dst = (unsigned char *)malloc(nvox);

    src[nvox / 2] = 200; /* any non-zero is foreground */
    memset(dst, 9, nvox);
    morph_box(src, dst, dims, 0, MORPH_BOX_DILATE);
    MU_ASSERT("radius 0 is a 0/1 copy",
              dst[nvox / 2] == 1 && dst[0] == 0 && dst[nvox - 1] == 0);

    morph_box(src, dst, dims, 1, MORPH_BOX_DILATE);
    for (i = 0, sum = 0; i < nvox; i++)
        sum += dst[i];
    MU_ASSERT("one voxel dilated by 1 is a 3x3x3 cube", sum == 27);

    morph_box(src, dst, dims, 100, MORPH_BOX_DILATE);
    for (i = 0, sum = 0; i < nvox; i++)
        sum += dst[i];
    MU_ASSERT("a cube larger than the volume fills it", sum == nvox);

    morph_box(src, dst, dims, 100, MORPH_BOX_ERODE);
    for (i = 0, sum = 0; i < nvox; i++)
        sum += dst[i];
    MU_ASSERT("a cube larger than the volume erodes everything", sum == 0);

    memset(src, 1, nvox);
    morph_box(src, dst, dims, 2, MORPH_BOX_DILATE);
    for (i = 0, sum = 0; i < nvox; i++)
        sum += dst[i];
    MU_ASSERT("dilating a full volume keeps it full (no wrap)", sum == nvox);

    morph_box(src, dst, dims, 2, MORPH_BOX_ERODE);
    for (i = 0, sum = 0; i < nvox; i++)
        sum += dst[i];
    MU_ASSERT("eroding a full volume peels r voxels off every face",
              sum == (size_t)(dims[0] - 4) * (dims[1] - 4) * (dims[2] - 4));

    morph_box(src, dst, dims, 2, MORPH_BOX_CLOSE);
    for (i = 0, sum = 0; i < nvox; i++)
        sum += dst[i];
    /* the erosion half of the closing sees the outside as background, so a
       volume that touches its border loses r voxels there (as scipy does) */
    MU_ASSERT("closing a full volume peels the border like an erosion",
              sum == (size_t)(dims[0] - 4) * (dims[1] - 4) * (dims[2] - 4));

    memset(dst, 9, nvox);
    morph_box(src, dst, dims, 2, 99);
    MU_ASSERT("an unknown operation leaves dst alone", dst[0] == 9 && dst[nvox - 1] == 9);
    free(src);
    free(dst);
}

/*
 * morph_dilate(), morph_close() and morph_open() repeat a 3x3x3 step niter
 * times, which is a (2*niter+1)^3 cube.  An earlier version took the 27-voxel
 * sum with an 8-bit convolution that rounded sum / 3, so a voxel whose only
 * foreground neighbour was one voxel (1/3 -> 0) was never dilated into -- an
 * isolated voxel disappeared -- and an erosion accepted 26 of 27.  These tests
 * hold the steps to the definition.
 */

/** \brief Closing as libCAT defines it: the volume sits in an infinite background. */
static void
close_reference(const unsigned char *src, unsigned char *dst, int dims[3], int r)
{
    int pd[3] = {dims[0] + 2 * r, dims[1] + 2 * r, dims[2] + 2 * r};
    size_t pn = (size_t)pd[0] * pd[1] * pd[2];
    unsigned char *a = (unsigned char *)calloc(pn, 1);
    unsigned char *b = (unsigned char *)malloc(pn);
    unsigned char *c = (unsigned char *)malloc(pn);
    int x, y, z;

    for (z = 0; z < dims[2]; z++)
        for (y = 0; y < dims[1]; y++)
            for (x = 0; x < dims[0]; x++)
                a[(size_t)(z + r) * pd[0] * pd[1] + (size_t)(y + r) * pd[0] + (x + r)] =
                    src[(size_t)z * dims[0] * dims[1] + (size_t)y * dims[0] + x] != 0;
    box_reference(a, b, pd, r, 0);
    box_reference(b, c, pd, r, 1);
    for (z = 0; z < dims[2]; z++)
        for (y = 0; y < dims[1]; y++)
            for (x = 0; x < dims[0]; x++)
                dst[(size_t)z * dims[0] * dims[1] + (size_t)y * dims[0] + x] =
                    c[(size_t)(z + r) * pd[0] * pd[1] + (size_t)(y + r) * pd[0] + (x + r)];
    free(a);
    free(b);
    free(c);
}

/** \brief Zero every voxel within `margin` of the volume border. */
static void
clear_margin(unsigned char *mask, int dims[3], int margin)
{
    int x, y, z;
    for (z = 0; z < dims[2]; z++)
        for (y = 0; y < dims[1]; y++)
            for (x = 0; x < dims[0]; x++)
                if (x < margin || y < margin || z < margin || x >= dims[0] - margin ||
                    y >= dims[1] - margin || z >= dims[2] - margin)
                    mask[(size_t)z * dims[0] * dims[1] + (size_t)y * dims[0] + x] = 0;
}

static void
test_morph_dilate_is_exact(void)
{
    int dims[3] = {21, 19, 17};
    size_t nvox = (size_t)dims[0] * dims[1] * dims[2], i;
    static const int density[3] = {5, 300, 930};
    static const int niters[4] = {1, 2, 3, 5};
    unsigned char *mask = (unsigned char *)malloc(nvox);
    unsigned char *want = (unsigned char *)malloc(nvox);
    unsigned char *got = (unsigned char *)malloc(nvox);
    float *fl = (float *)malloc(nvox * sizeof(float));
    int d, k, mismatch = 0, runs = 0, sum;

    for (d = 0; d < 3; d++)
        for (k = 0; k < 4; k++)
        {
            fill_density(mask, nvox, 100 + 4 * d + k, density[d]);
            box_reference(mask, want, dims, niters[k], 0);

            memcpy(got, mask, nvox);
            morph_dilate(got, dims, niters[k], 0.5, DT_UINT8);
            runs++;
            if (memcmp(got, want, nvox) != 0)
                mismatch++;

            for (i = 0; i < nvox; i++)
                fl[i] = mask[i] ? 1.0f : 0.0f;
            morph_dilate(fl, dims, niters[k], 0.5, DT_FLOAT32);
            for (i = 0; i < nvox; i++)
                got[i] = fl[i] != 0.0f;
            if (memcmp(got, want, nvox) != 0)
                mismatch++;
        }
    MU_ASSERT("morph_dilate is the exact (2r+1)^3 dilation", mismatch == 0 && runs == 12);

    memset(mask, 0, nvox);
    mask[nvox / 2] = 1;
    morph_dilate(mask, dims, 1, 0.5, DT_UINT8);
    for (i = 0, sum = 0; i < nvox; i++)
        sum += mask[i];
    MU_ASSERT("an isolated voxel dilates to a 3x3x3 cube", sum == 27);
    free(mask);
    free(want);
    free(got);
    free(fl);
}

static void
test_morph_close_open_are_exact(void)
{
    int dims[3] = {21, 19, 17};
    size_t nvox = (size_t)dims[0] * dims[1] * dims[2];
    static const int density[3] = {5, 300, 930};
    unsigned char *mask = (unsigned char *)malloc(nvox);
    unsigned char *want = (unsigned char *)malloc(nvox);
    unsigned char *got = (unsigned char *)malloc(nvox);
    int d, r, mismatch_close = 0, mismatch_open = 0;

    for (d = 0; d < 3; d++)
        for (r = 1; r <= 3; r++)
        {
            /* closing: any mask, the foreground may touch the border */
            fill_density(mask, nvox, 200 + 3 * d + r, density[d]);
            close_reference(mask, want, dims, r);
            memcpy(got, mask, nvox);
            morph_close(got, dims, r, 0.5, DT_UINT8);
            if (memcmp(got, want, nvox) != 0)
                mismatch_close++;

            /* opening: foreground kept r voxels from the border, where the
               border convention of the erosion does not come into it */
            fill_density(mask, nvox, 300 + 3 * d + r, density[d]);
            clear_margin(mask, dims, r + 1);
            box_reference_op(mask, want, dims, r, MORPH_BOX_OPEN);
            memcpy(got, mask, nvox);
            morph_open(got, dims, r, 0.5, 0, DT_UINT8);
            if (memcmp(got, want, nvox) != 0)
                mismatch_open++;
        }
    MU_ASSERT("morph_close is a dilation then an erosion of the exact cube", mismatch_close == 0);
    MU_ASSERT("morph_open is an erosion then a dilation of the exact cube", mismatch_open == 0);
    free(mask);
    free(want);
    free(got);
}

int main(void)
{
    MU_RUN_TEST(test_matches_the_dense_version);
    MU_RUN_TEST(test_region_null_grows_everywhere);
    MU_RUN_TEST(test_edge_cases);
    MU_RUN_TEST(test_box_matches_the_definition);
    MU_RUN_TEST(test_box_in_place);
    MU_RUN_TEST(test_box_edge_cases);
    MU_RUN_TEST(test_morph_dilate_is_exact);
    MU_RUN_TEST(test_morph_close_open_are_exact);
    printf("%d tests run, %d failed\n", tests_run, tests_failed);
    return tests_failed ? 1 : 0;
}
