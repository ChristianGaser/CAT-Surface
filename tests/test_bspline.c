#include "minunit.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "CAT_Vol.h"

/*
 * bspline_sample() interpolates a volume with a quadratic B-spline at
 * positions normalised to [-1, 1], the way spline_resize.grid_sample does in
 * T1Prep.  The reference below is the same algorithm written the obvious way,
 * in double precision throughout: one line at a time for the prefilter, a
 * direct sum over the 27 taps for the samples.  The two have to agree to float
 * rounding, wherever the position lies, including outside the volume.
 */

static int DIMS[3] = {23, 19, 17};

/** \brief Prefilter of one line, as spline_resize's filter_last_dim. */
static void
ref_filter_line(double *c, int n)
{
    const double z = sqrt(8.0) - 3.0;
    const int horizon = (int)ceil(-30.0 / log(fabs(z)));
    double zk, acc;
    int i, k;

    if (n < 2)
        return;
    for (i = 0; i < n; i++)
        c[i] *= (1.0 - z) * (1.0 - 1.0 / z);
    if (horizon < n)
    {
        for (k = 0, acc = 0.0, zk = 1.0; k < horizon; k++, zk *= z)
            acc += c[k] * zk;
        c[0] = acc;
    }
    else
    {
        const double a = pow(z, n - 1);
        for (k = 1, acc = c[0], zk = 1.0; k < n - 1; k++)
        {
            zk *= z;
            acc += c[k] * (zk + a * a / zk);
        }
        c[0] = (acc + a * c[n - 1]) / (1.0 - a * a);
    }
    for (i = 1; i < n; i++)
        c[i] += z * c[i - 1];
    c[n - 1] = (z * c[n - 2] + c[n - 1]) * (z / (z * z - 1.0));
    for (i = n - 2; i >= 0; i--)
        c[i] = z * (c[i + 1] - c[i]);
}

/** \brief Prefilter of the volume, along z, then y, then x. */
static void
ref_prefilter(double *v, const int d[3])
{
    int n = d[0] > d[1] ? (d[0] > d[2] ? d[0] : d[2]) : (d[1] > d[2] ? d[1] : d[2]);
    double *line = (double *)malloc(sizeof(double) * n);
    int x, y, z;

#define AT(X, Y, Z) v[((size_t)(Z) * d[1] + (Y)) * d[0] + (X)]
    for (y = 0; y < d[1]; y++)
        for (x = 0; x < d[0]; x++)
        {
            for (z = 0; z < d[2]; z++)
                line[z] = AT(x, y, z);
            ref_filter_line(line, d[2]);
            for (z = 0; z < d[2]; z++)
                AT(x, y, z) = line[z];
        }
    for (z = 0; z < d[2]; z++)
        for (x = 0; x < d[0]; x++)
        {
            for (y = 0; y < d[1]; y++)
                line[y] = AT(x, y, z);
            ref_filter_line(line, d[1]);
            for (y = 0; y < d[1]; y++)
                AT(x, y, z) = line[y];
        }
    for (z = 0; z < d[2]; z++)
        for (y = 0; y < d[1]; y++)
        {
            for (x = 0; x < d[0]; x++)
                line[x] = AT(x, y, z);
            ref_filter_line(line, d[0]);
            for (x = 0; x < d[0]; x++)
                AT(x, y, z) = line[x];
        }
#undef AT
    free(line);
}

static double
ref_weight(double x)
{
    x = fabs(x);
    if (x < 0.5)
        return 0.75 - x * x;
    return (x < 1.5) ? 0.5 * (1.5 - x) * (1.5 - x) : 0.0;
}

static double
ref_coord(double g, int n, int align_corners)
{
    return align_corners ? (g + 1.0) * 0.5 * (n - 1) : ((g + 1.0) * n - 1.0) * 0.5;
}

/** \brief Direct sum over the 27 taps around one position; outside taps are 0. */
static double
ref_sample(const double *c, const int d[3], const float g[3], int align_corners)
{
    double u[3], sum = 0.0;
    int base[3], a, i, j, k;

    for (a = 0; a < 3; a++)
    {
        u[a] = ref_coord(g[a], d[a], align_corners);
        if (!isfinite(u[a]))
            return 0.0;
        base[a] = (int)floor(u[a] - 0.5);
    }
    for (k = base[2]; k < base[2] + 3; k++)
        for (j = base[1]; j < base[1] + 3; j++)
            for (i = base[0]; i < base[0] + 3; i++)
            {
                if (i < 0 || j < 0 || k < 0 || i >= d[0] || j >= d[1] || k >= d[2])
                    continue;
                sum += c[((size_t)k * d[1] + j) * d[0] + i] *
                       ref_weight(u[0] - i) * ref_weight(u[1] - j) * ref_weight(u[2] - k);
            }
    return sum;
}

static void
fill_random(float *v, size_t n, int seed)
{
    size_t i;
    srand(seed);
    for (i = 0; i < n; i++)
        v[i] = (float)(rand() % 1000) / 10.0f;
}

/** \brief The normalised coordinate of voxel centre i on an axis of n voxels. */
static float
centre(int i, int n)
{
    return (float)(2.0 * i / (n - 1) - 1.0);
}

static void
test_interpolates_the_voxel_centres(void)
{
    size_t nvox = (size_t)DIMS[0] * DIMS[1] * DIMS[2], p, n = 0;
    float *vol = (float *)malloc(sizeof(float) * nvox);
    float *grid = (float *)malloc(sizeof(float) * 3 * nvox);
    float *out = (float *)malloc(sizeof(float) * nvox);
    float *want = (float *)malloc(sizeof(float) * nvox);
    double maxerr = 0.0;
    int x, y, z;

    fill_random(vol, nvox, 1);
    /* interior centres: at the outermost voxels the zero padding of the taps
       takes away the mirrored neighbour, as it does in spline_resize */
    for (z = 1; z < DIMS[2] - 1; z++)
        for (y = 1; y < DIMS[1] - 1; y++)
            for (x = 1; x < DIMS[0] - 1; x++, n++)
            {
                grid[3 * n] = centre(x, DIMS[0]);
                grid[3 * n + 1] = centre(y, DIMS[1]);
                grid[3 * n + 2] = centre(z, DIMS[2]);
                want[n] = vol[((size_t)z * DIMS[1] + y) * DIMS[0] + x];
            }
    bspline_sample(vol, DIMS, grid, n, out, 1, 1, 0, 0.0f);
    for (p = 0; p < n; p++)
        maxerr = fmax(maxerr, fabs((double)out[p] - want[p]));
    MU_ASSERT("the samples reproduce the volume at its voxel centres", maxerr < 1e-3);
    free(vol);
    free(grid);
    free(out);
    free(want);
}

static void
test_reproduces_a_constant(void)
{
    size_t nvox = (size_t)DIMS[0] * DIMS[1] * DIMS[2], p, n = 5000;
    float *vol = (float *)malloc(sizeof(float) * nvox);
    float *grid = (float *)malloc(sizeof(float) * 3 * n);
    float *out = (float *)malloc(sizeof(float) * n);
    double maxerr = 0.0;
    int a;

    for (p = 0; p < nvox; p++)
        vol[p] = 7.25f;
    srand(2);
    /* anywhere at least 1.5 voxels inside, so that all 27 taps exist */
    for (p = 0; p < n; p++)
        for (a = 0; a < 3; a++)
        {
            double u = 1.5 + (DIMS[a] - 4.0) * rand() / (double)RAND_MAX;
            grid[3 * p + a] = (float)(2.0 * u / (DIMS[a] - 1) - 1.0);
        }
    bspline_sample(vol, DIMS, grid, n, out, 1, 1, 0, 0.0f);
    for (p = 0; p < n; p++)
        maxerr = fmax(maxerr, fabs(out[p] - 7.25));
    MU_ASSERT("a constant volume samples as the constant", maxerr < 1e-4);
    free(vol);
    free(grid);
    free(out);
}

static void
test_matches_the_reference(void)
{
    size_t nvox = (size_t)DIMS[0] * DIMS[1] * DIMS[2], p, n = 7919; /* not a multiple of the threads */
    float *vol = (float *)malloc(sizeof(float) * nvox);
    double *coef = (double *)malloc(sizeof(double) * nvox);
    float *grid = (float *)malloc(sizeof(float) * 3 * n);
    float *out = (float *)malloc(sizeof(float) * n);
    int align, pre, mask, a, mismatch = 0, runs = 0;

    fill_random(vol, nvox, 3);
    for (p = 0; p < nvox; p += 7)
        vol[p] = 0.0f; /* something for mask_value to find */
    srand(4);
    /* positions from well outside to well inside, on every axis */
    for (p = 0; p < n; p++)
        for (a = 0; a < 3; a++)
            grid[3 * p + a] = (float)(-1.3 + 2.6 * rand() / (double)RAND_MAX);

    for (align = 0; align <= 1; align++)
        for (pre = 0; pre <= 1; pre++)
            for (mask = 0; mask <= 1; mask++)
            {
                double maxerr = 0.0;

                for (p = 0; p < nvox; p++)
                    coef[p] = vol[p];
                if (pre)
                {
                    ref_prefilter(coef, DIMS);
                    if (mask)
                        for (p = 0; p < nvox; p++)
                            if (vol[p] == 0.0f)
                                coef[p] = 0.0;
                }
                bspline_sample(vol, DIMS, grid, n, out, align, pre, mask, 0.0f);
                for (p = 0; p < n; p++)
                    maxerr = fmax(maxerr, fabs(out[p] - ref_sample(coef, DIMS, grid + 3 * p, align)));
                runs++;
                if (!(maxerr < 1e-3)) /* the volume runs to 100 */
                {
                    mismatch++;
                    printf("  align_corners=%d prefilter=%d mask=%d: max |diff| %g\n",
                           align, pre, mask, maxerr);
                }
            }
    MU_ASSERT("bspline_sample matches the reference", mismatch == 0 && runs == 8);
    free(vol);
    free(coef);
    free(grid);
    free(out);
}

static void
test_edge_cases(void)
{
    size_t nvox = (size_t)DIMS[0] * DIMS[1] * DIMS[2];
    float *vol = (float *)malloc(sizeof(float) * nvox);
    float grid[12] = {NAN, 0.0f, 0.0f, 0.0f, INFINITY, 0.0f, 50.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f};
    float out[4] = {-1.0f, -1.0f, -1.0f, -1.0f};
    float line[5] = {1.0f, 2.0f, 4.0f, 8.0f, 16.0f};
    int thin[3] = {5, 1, 1};
    float at[3] = {0.0f, 0.0f, 0.0f};
    float res = -1.0f;

    fill_random(vol, nvox, 5);
    bspline_sample(vol, DIMS, grid, 0, out, 1, 1, 0, 0.0f);
    MU_ASSERT("no positions leave the output alone", out[0] == -1.0f);

    bspline_sample(vol, DIMS, grid, 4, out, 1, 1, 0, 0.0f);
    MU_ASSERT("a non-finite position samples as 0", out[0] == 0.0f && out[1] == 0.0f);
    MU_ASSERT("a position far outside samples as 0", out[2] == 0.0f);
    MU_ASSERT("the centre samples as something inside", out[3] != 0.0f && isfinite(out[3]));

    /* axes of a single voxel: no division by zero, a finite value */
    bspline_sample(line, thin, at, 1, &res, 1, 1, 0, 0.0f);
    MU_ASSERT("axes of one voxel give a finite sample", isfinite(res));
    free(vol);
}

int main(void)
{
    MU_RUN_TEST(test_interpolates_the_voxel_centres);
    MU_RUN_TEST(test_reproduces_a_constant);
    MU_RUN_TEST(test_matches_the_reference);
    MU_RUN_TEST(test_edge_cases);
    printf("%d tests run, %d failed\n", tests_run, tests_failed);
    return tests_failed ? 1 : 0;
}
