#include "minunit.h"
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "CAT_Cat12Utils.h"

/*
 * Lib/CAT_Cat12Utils.c ports CAT12's cat_vol_smooth3X, cat_vol_div,
 * cat_vol_morph, cat_vol_localstat, cat_vol_laplace3R and cat_vbdist.  The
 * voxel-for-voxel comparison with MATLAB lives in tests/cat12_utils_reference;
 * these tests check the properties that follow from the definitions and need
 * no MATLAB: what a filter does to a constant or linear field, how far a
 * dilation grows, which component survives, and what happens at the border.
 */

#define IDX(x, y, z, d) (((size_t)(z) * (d)[1] + (y)) * (d)[0] + (x))

static size_t
nvox(const int d[3])
{
    return (size_t)d[0] * d[1] * d[2];
}

static long
count(const unsigned char *v, size_t n)
{
    size_t i;
    long c = 0;

    for (i = 0; i < n; i++)
        c += v[i] != 0;
    return c;
}

/** \brief Solid box [x0,x1] x [y0,y1] x [z0,z1] set to 1. */
static void
box(unsigned char *v, const int d[3], int x0, int x1, int y0, int y1, int z0, int z1)
{
    int x, y, z;

    for (z = z0; z <= z1; z++)
        for (y = y0; y <= y1; y++)
            for (x = x0; x <= x1; x++)
                v[IDX(x, y, z, d)] = 1;
}

static void
test_smooth_keeps_constant(void)
{
    int d[3] = {23, 19, 17}, small[3] = {5, 9, 6};
    double s[7] = {0.2, 0.45, 0.75, 1.0, 1.7, 2.4, 4.0};
    float *v = (float *)malloc(nvox(d) * sizeof(float));
    size_t i;
    int k, bad = 0;

    for (k = 0; k < 7; k++)
    {
        for (i = 0; i < nvox(d); i++)
            v[i] = 3.25f;
        cat12_vol_smooth3X(v, d, s[k]);
        for (i = 0; i < nvox(d); i++)
            bad += fabsf(v[i] - 3.25f) > 1e-5f;
    }
    /* a dimension <= 6: the 5x5x5 kernel */
    for (i = 0; i < nvox(small); i++)
        v[i] = -1.5f;
    cat12_vol_smooth3X(v, small, 2.0);
    for (i = 0; i < nvox(small); i++)
        bad += fabsf(v[i] + 1.5f) > 1e-5f;
    MU_ASSERT("smooth3X keeps a constant field", bad == 0);
    free(v);
}

static void
test_smooth_spreads_and_cleans(void)
{
    int d[3] = {21, 21, 21};
    float *v = (float *)calloc(nvox(d), sizeof(float));
    size_t i, c = IDX(10, 10, 10, d);
    int bad = 0;

    /* a point spreads symmetrically and the peak drops */
    v[c] = 1.0f;
    cat12_vol_smooth3X(v, d, 0.8);
    MU_ASSERT("smooth3X lowers a point", v[c] < 1.0f && v[c] > 0.1f);
    MU_ASSERT("smooth3X is symmetric",
              v[IDX(9, 10, 10, d)] == v[IDX(11, 10, 10, d)] &&
                  v[IDX(10, 9, 10, d)] == v[IDX(10, 11, 10, d)] &&
                  v[IDX(10, 10, 9, d)] == v[IDX(10, 10, 11, d)]);
    MU_ASSERT("smooth3X with a 3x3x3 kernel reaches only the neighbours",
              v[IDX(12, 10, 10, d)] == 0.0f);

    /* s == 0 leaves the volume alone, but NaN and Inf become 0 */
    for (i = 0; i < nvox(d); i++)
        v[i] = 1.0f;
    v[c] = NAN;
    v[0] = INFINITY;
    cat12_vol_smooth3X(v, d, 0.0);
    MU_ASSERT("smooth3X with s = 0 zeroes NaN and Inf only",
              v[c] == 0.0f && v[0] == 0.0f && v[1] == 1.0f);
    v[c] = NAN;
    cat12_vol_smooth3X(v, d, 2.0);
    for (i = 0; i < nvox(d); i++)
        bad += !isfinite(v[i]);
    MU_ASSERT("smooth3X output is finite", bad == 0);
    free(v);
}

static void
test_dilate_box(void)
{
    int d[3] = {15, 16, 17};
    unsigned char *v = (unsigned char *)calloc(nvox(d), 1);
    double one[3] = {1, 1, 1}, aniso[3] = {1, 2, 4};

    v[IDX(7, 8, 8, d)] = 1;
    cat12_vol_morph_dilate(v, d, 1.0, NULL);
    MU_ASSERT("'d' 1 grows a voxel to 3x3x3", count(v, nvox(d)) == 27);

    memset(v, 0, nvox(d));
    v[IDX(7, 8, 8, d)] = 1;
    cat12_vol_morph_dilate(v, d, 1.4, one);
    MU_ASSERT("'d' 1.4 rounds to 3x3x3", count(v, nvox(d)) == 27);

    memset(v, 0, nvox(d));
    v[IDX(7, 8, 8, d)] = 1;
    cat12_vol_morph_dilate(v, d, 1.5, one);
    MU_ASSERT("'d' 1.5 rounds to 5x5x5", count(v, nvox(d)) == 125);

    memset(v, 0, nvox(d));
    v[IDX(7, 8, 8, d)] = 1;
    cat12_vol_morph_dilate(v, d, 2.0, aniso);
    MU_ASSERT("'d' 2 with voxels 1x2x4 gives 5x3x3",
              count(v, nvox(d)) == 45 && v[IDX(9, 9, 9, d)] && !v[IDX(7, 10, 8, d)]);

    memset(v, 0, nvox(d));
    v[0] = 1;
    cat12_vol_morph_dilate(v, d, 1.0, NULL);
    MU_ASSERT("'d' at the corner does not wrap", count(v, nvox(d)) == 8);

    memset(v, 0, nvox(d));
    v[IDX(7, 8, 8, d)] = 1;
    cat12_vol_morph_dilate(v, d, 0.4, NULL);
    MU_ASSERT("'d' with round(n) == 0 is a no-op", count(v, nvox(d)) == 1);
    free(v);
}

static void
test_lab(void)
{
    int d[3] = {20, 20, 20};
    unsigned char *v = (unsigned char *)calloc(nvox(d), 1);

    box(v, d, 2, 4, 2, 4, 2, 4);      /* 27 voxels */
    box(v, d, 10, 14, 10, 13, 10, 13); /* 80 voxels */
    v[IDX(6, 3, 3, d)] = 1;            /* 1 voxel */
    cat12_vol_morph_lab(v, d);
    MU_ASSERT("'l' keeps the largest component",
              count(v, nvox(d)) == 80 && v[IDX(12, 12, 12, d)] && !v[IDX(3, 3, 3, d)]);

    memset(v, 0, nvox(d));
    box(v, d, 12, 14, 12, 14, 12, 14); /* first voxel later in memory */
    box(v, d, 2, 4, 2, 4, 2, 4);
    cat12_vol_morph_lab(v, d);
    MU_ASSERT("'l' resolves a tie in memory order",
              count(v, nvox(d)) == 27 && v[IDX(3, 3, 3, d)]);

    memset(v, 0, nvox(d));
    box(v, d, 2, 4, 2, 4, 2, 4);
    box(v, d, 5, 7, 5, 7, 5, 7); /* touches the first only at an edge */
    cat12_vol_morph_lab(v, d);
    MU_ASSERT("'l' uses 6-connectivity", count(v, nvox(d)) == 27);

    memset(v, 0, nvox(d));
    cat12_vol_morph_lab(v, d);
    MU_ASSERT("'l' of an empty volume is empty", count(v, nvox(d)) == 0);

    /* [10 0.1]: components above 10% of the foreground */
    box(v, d, 0, 9, 0, 9, 0, 3);      /* 400 */
    box(v, d, 12, 19, 12, 19, 12, 13); /* 128 */
    box(v, d, 0, 2, 15, 17, 15, 17);  /* 27 */
    cat12_vol_morph_lab_multi(v, d, 10, 0.1);
    MU_ASSERT("'l' [10 0.1] keeps the two large components", count(v, nvox(d)) == 528);
    free(v);
}

static void
test_labdistclose_fills(void)
{
    int d[3] = {40, 40, 40}, x, y, z;
    unsigned char *v = (unsigned char *)calloc(nvox(d), 1);
    long n;

    /* hollow ball: radius 14, cavity radius 9 */
    for (z = 0; z < d[2]; z++)
        for (y = 0; y < d[1]; y++)
            for (x = 0; x < d[0]; x++)
            {
                double r2 = (x - 20) * (x - 20) + (y - 20) * (y - 20) + (z - 20) * (z - 20);
                v[IDX(x, y, z, d)] = r2 <= 14 * 14 && r2 > 9 * 9;
            }
    n = count(v, nvox(d));
    cat12_vol_morph_labdistclose(v, d, 4.8, NULL);
    MU_ASSERT("'ldc' fills the cavity of a hollow ball",
              v[IDX(20, 20, 20, d)] && count(v, nvox(d)) > n);
    MU_ASSERT("'ldc' does not grow into the background",
              !v[IDX(20, 20, 36, d)] && !v[0]);

    /* the vbdist path (n > 5) too */
    for (z = 0; z < d[2]; z++)
        for (y = 0; y < d[1]; y++)
            for (x = 0; x < d[0]; x++)
            {
                double r2 = (x - 20) * (x - 20) + (y - 20) * (y - 20) + (z - 20) * (z - 20);
                v[IDX(x, y, z, d)] = r2 <= 14 * 14 && r2 > 9 * 9;
            }
    cat12_vol_morph_labdistclose(v, d, 8.0, NULL);
    MU_ASSERT("'ldc' 8 fills the cavity", v[IDX(20, 20, 20, d)] && !v[0]);

    /* a notch narrower than 2n closes */
    memset(v, 0, nvox(d));
    box(v, d, 5, 34, 5, 34, 5, 34);
    for (z = 5; z <= 34; z++)
        for (x = 18; x <= 21; x++)
            for (y = 5; y <= 25; y++)
                v[IDX(x, y, z, d)] = 0;
    cat12_vol_morph_labdistclose(v, d, 4.8, NULL);
    MU_ASSERT("'ldc' closes a 4-voxel notch", v[IDX(20, 20, 20, d)]);
    free(v);
}

static void
test_labdistopen_and_open(void)
{
    int d[3] = {40, 30, 30};
    unsigned char *v = (unsigned char *)calloc(nvox(d), 1);
    double vx[3] = {1.5, 1.5, 1.5};

    /* two cubes joined by a one-voxel bridge, plus a speck */
    box(v, d, 3, 14, 5, 16, 5, 16);
    box(v, d, 22, 35, 5, 18, 5, 18);
    box(v, d, 15, 21, 10, 10, 10, 10);
    v[IDX(30, 26, 26, d)] = 1;
    cat12_vol_morph_labdistopen(v, d, 3.0, vx);
    MU_ASSERT("'ldo' cuts the bridge and keeps the larger cube",
              v[IDX(28, 12, 12, d)] && !v[IDX(8, 10, 10, d)] && !v[IDX(18, 10, 10, d)] &&
                  !v[IDX(30, 26, 26, d)]);

    /* 'o' 1: a slab of thickness 2 survives at the border, not inside */
    memset(v, 0, nvox(d));
    box(v, d, 0, 1, 0, 29, 0, 29);
    box(v, d, 20, 21, 0, 29, 0, 29);
    v[IDX(10, 15, 15, d)] = 1;
    cat12_vol_morph_open(v, d, 1.0, NULL);
    MU_ASSERT("'o' counts the outside as foreground",
              v[IDX(0, 15, 15, d)] && v[IDX(1, 15, 15, d)] && !v[IDX(20, 15, 15, d)] &&
                  !v[IDX(10, 15, 15, d)]);

    memset(v, 0, nvox(d));
    box(v, d, 10, 20, 10, 20, 10, 20);
    cat12_vol_morph_open(v, d, 1.0, NULL);
    MU_ASSERT("'o' keeps a cube", count(v, nvox(d)) == 11 * 11 * 11);
    free(v);
}

static void
test_vbdist(void)
{
    int d[3] = {12, 11, 10}, x, y, z, bad = 0, badi = 0;
    size_t n = (size_t)12 * 11 * 10, c = IDX(4, 5, 6, d), i;
    float *v = (float *)calloc(n, sizeof(float)), *D = (float *)malloc(n * sizeof(float));
    unsigned int *I = (unsigned int *)malloc(n * sizeof(unsigned int));
    unsigned char *m = (unsigned char *)malloc(n);
    double vx[3] = {1.0, 1.5, 2.25};

    v[c] = 1.0f;
    cat12_vbdist(v, NULL, d, vx, D, I);
    for (z = 0; z < d[2]; z++)
        for (y = 0; y < d[1]; y++)
            for (x = 0; x < d[0]; x++)
            {
                double e = sqrt(pow((x - 4) * 1.0, 2) + pow((y - 5) * 1.5, 2) + pow((z - 6) * 2.25, 2));
                bad += fabs(D[IDX(x, y, z, d)] - e) > 1e-5;
                badi += I[IDX(x, y, z, d)] != c;
            }
    MU_ASSERT("vbdist of one voxel is the Euclidean distance", bad == 0 && badi == 0);

    for (i = 0; i < n; i++)
        m[i] = (i % 12) < 8;
    cat12_vbdist(v, m, d, NULL, D, NULL);
    MU_ASSERT("vbdist leaves FLT_MAX outside the mask",
              D[IDX(10, 5, 6, d)] == FLT_MAX && D[IDX(6, 5, 6, d)] == 2.0f && D[c] == 0.0f);
    free(m);
    free(I);
    free(D);
    free(v);
}

static void
test_div(void)
{
    int d[3] = {16, 14, 12}, x, y, z, bad = 0;
    size_t n = (size_t)16 * 14 * 12, i;
    float *v = (float *)malloc(n * sizeof(float)), *o = (float *)malloc(n * sizeof(float));
    const float a = 0.01f;

    /* linear field: divergence 0 away from the border */
    for (z = 0; z < d[2]; z++)
        for (y = 0; y < d[1]; y++)
            for (x = 0; x < d[0]; x++)
                v[IDX(x, y, z, d)] = 1.0f + 0.1f * x + 0.05f * y - 0.02f * z;
    cat12_vol_div(v, d, o);
    for (z = 2; z < d[2] - 2; z++)
        for (y = 2; y < d[1] - 2; y++)
            for (x = 2; x < d[0] - 2; x++)
                bad += fabsf(o[IDX(x, y, z, d)]) > 1e-6f;
    MU_ASSERT("div of a linear field is 0 inside", bad == 0);

    /* a*x^2 + a*y^2: (2a + 2a)/1.5 */
    for (z = 0; z < d[2]; z++)
        for (y = 0; y < d[1]; y++)
            for (x = 0; x < d[0]; x++)
                v[IDX(x, y, z, d)] = 1.0f + a * x * x + a * y * y;
    cat12_vol_div(v, d, o);
    bad = 0;
    for (z = 2; z < d[2] - 2; z++)
        for (y = 2; y < d[1] - 2; y++)
            for (x = 2; x < d[0] - 2; x++)
                bad += fabsf(o[IDX(x, y, z, d)] - 4.0f * a / 1.5f) > 1e-5f;
    MU_ASSERT("div of a paraboloid is its Laplacian / 1.5", bad == 0);

    /* values below 1/3 are clamped: no divergence; NaN gives 0 */
    for (i = 0; i < n; i++)
        v[i] = -2.0f + 0.001f * (float)(i % 17);
    v[n / 2] = NAN;
    cat12_vol_div(v, d, o);
    bad = 0;
    for (i = 0; i < n; i++)
        bad += o[i] != 0.0f;
    MU_ASSERT("div ignores values below 1/3 and zeroes NaN", bad == 0);

    /* in place */
    for (z = 0; z < d[2]; z++)
        for (y = 0; y < d[1]; y++)
            for (x = 0; x < d[0]; x++)
                v[IDX(x, y, z, d)] = 1.0f + a * x * x;
    cat12_vol_div(v, d, o);
    cat12_vol_div(v, d, v);
    MU_ASSERT("div works in place", memcmp(v, o, n * sizeof(float)) == 0);
    free(o);
    free(v);
}

static void
test_localstat(void)
{
    int d[3] = {6, 5, 4};
    size_t n = 120, i;
    float *v = (float *)malloc(n * sizeof(float)), *o = (float *)malloc(n * sizeof(float));
    unsigned char *m = (unsigned char *)malloc(n);
    int bad = 0;

    for (i = 0; i < n; i++)
        v[i] = 2.0f;
    cat12_vol_localstat_mean(v, NULL, d, 1, o);
    for (i = 0; i < n; i++)
        bad += o[i] != 2.0f;
    MU_ASSERT("localstat keeps a constant field", bad == 0);

    /* zeros are neither averaged nor filled */
    v[IDX(2, 2, 2, d)] = 0.0f;
    v[IDX(3, 2, 2, d)] = 5.0f;
    cat12_vol_localstat_mean(v, NULL, d, 1, o);
    MU_ASSERT("localstat drops zeros",
              o[IDX(2, 2, 2, d)] == 0.0f && o[IDX(3, 2, 2, d)] == (5.0f + 2.0f * 5) / 6.0f &&
                  o[IDX(1, 2, 2, d)] == 2.0f);

    /* the mask limits both centres and neighbours */
    for (i = 0; i < n; i++)
    {
        v[i] = (float)(i % 7) + 1.0f;
        m[i] = i % 3 != 0;
    }
    cat12_vol_localstat_mean(v, m, d, 1, o);
    MU_ASSERT("localstat is 0 outside the mask", o[0] == 0.0f && o[3] == 0.0f);
    MU_ASSERT("localstat averages the face neighbours in the mask",
              o[IDX(1, 0, 0, d)] == (v[IDX(1, 0, 0, d)] + v[IDX(2, 0, 0, d)] + v[IDX(1, 1, 0, d)] +
                                     v[IDX(1, 0, 1, d)]) /
                                        4.0f);
    free(m);
    free(o);
    free(v);
}

static void
test_laplace(void)
{
    int d[3] = {21, 6, 6}, x, y, z, bad = 0;
    size_t n = (size_t)21 * 6 * 6;
    float *v = (float *)calloc(n, sizeof(float));
    unsigned char *m = (unsigned char *)calloc(n, 1);

    /* fixed planes x = 0 (0) and x = 20 (1), free inside: a linear ramp */
    for (z = 0; z < d[2]; z++)
        for (y = 0; y < d[1]; y++)
            for (x = 0; x < d[0]; x++)
            {
                v[IDX(x, y, z, d)] = (x == d[0] - 1) ? 1.0f : 0.0f;
                m[IDX(x, y, z, d)] = x > 0 && x < d[0] - 1;
            }
    cat12_vol_laplace3R(v, m, d, 1e-6);
    for (z = 0; z < d[2]; z++)
        for (y = 0; y < d[1]; y++)
            for (x = 0; x < d[0]; x++)
                bad += fabsf(v[IDX(x, y, z, d)] - x / 20.0f) > 2e-3f;
    MU_ASSERT("laplace3R relaxes to a linear ramp", bad == 0);
    free(m);
    free(v);
}

int
main(void)
{
    MU_RUN_TEST(test_smooth_keeps_constant);
    MU_RUN_TEST(test_smooth_spreads_and_cleans);
    MU_RUN_TEST(test_dilate_box);
    MU_RUN_TEST(test_lab);
    MU_RUN_TEST(test_labdistclose_fills);
    MU_RUN_TEST(test_labdistopen_and_open);
    MU_RUN_TEST(test_vbdist);
    MU_RUN_TEST(test_div);
    MU_RUN_TEST(test_localstat);
    MU_RUN_TEST(test_laplace);
    printf("Tests run: %d, failed: %d\n", tests_run, tests_failed);
    return tests_failed != 0;
}
