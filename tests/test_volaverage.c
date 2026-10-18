#include "minunit.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "CAT_VolumeStats.h"

/*
 * CAT_VolAverage takes every volume as one sample.
 *
 * It used to read one volume per FILE and divide by the number of files, so a
 * 4D file of 531 z-maps was a single sample of its first volume: the mean was
 * volume 0 and the standard deviation, 1/(n-1) with n = 1, NaN everywhere.
 * Where the samples are separate files nothing may change.
 */

#define NX 4
#define NY 3
#define NZ 2
#define NVOX (NX * NY * NZ)
#define NT 6

static const char *file_4d = "test_volaverage_4d.nii";
static const char *file_mid = "test_volaverage_mid.nii";
static char vol_names[NT][40];

/** \brief Voxel i of volume t: a bright "brain" with a little structure. */
static float
sample(int t, int i)
{
    return 50.0f + 3.0f * sinf(0.7f * (float)i + 1.3f * (float)t) + 0.25f * (float)t;
}

/** \brief Write nt volumes of nx*ny*nz float voxels as one NIfTI file (3D if nt is 1). */
static void
write_series(const char *file, int nx, int ny, int nz, int nt, const float *data)
{
    int dims[8] = {(nt > 1) ? 4 : 3, nx, ny, nz, nt, 1, 1, 1};
    nifti_image *nim = nifti_make_new_nim(dims, DT_FLOAT32, 1);

    memcpy(nim->data, data, sizeof(float) * (size_t)nx * ny * nz * nt);
    nifti_set_filenames(nim, file, 0, 1);
    nifti_image_write(nim);
    nifti_image_free(nim);
}

/**
 * \brief Write the series as one 4D file, as NT 3D files and as the 4D file
 *        file_mid with volumes 1 to 3; the last volume is shifted by outlier.
 */
static void
make_fixture(float outlier)
{
    float data[NT * NVOX];
    int t, i;

    for (t = 0; t < NT; t++)
        for (i = 0; i < NVOX; i++)
            data[t * NVOX + i] = sample(t, i) + ((t == NT - 1) ? outlier : 0.0f);

    write_series(file_4d, NX, NY, NZ, NT, data);
    write_series(file_mid, NX, NY, NZ, 3, data + NVOX);
    for (t = 0; t < NT; t++)
    {
        sprintf(vol_names[t], "test_volaverage_v%d.nii", t);
        write_series(vol_names[t], NX, NY, NZ, 1, data + t * NVOX);
    }
}

static void
remove_fixture(void)
{
    int t;

    remove(file_4d);
    remove(file_mid);
    for (t = 0; t < NT; t++)
        remove(vol_names[t]);
}

/** \brief Largest difference of two arrays is at most tol. */
static int
same_values(const double *a, const double *b, size_t n, double tol)
{
    size_t i;

    for (i = 0; i < n; i++)
        if (!(fabs(a[i] - b[i]) <= tol))
            return 0;
    return 1;
}

/** \brief Mean and std (n - 1) of the volumes first .. first + n - 1 by the two-pass formula. */
static void
reference(int first, int n, double *mean, double *sd)
{
    double s, d;
    int t, i;

    for (i = 0; i < NVOX; i++)
    {
        s = 0.0;
        for (t = first; t < first + n; t++)
            s += (double)sample(t, i);
        mean[i] = s / n;
        s = 0.0;
        for (t = first; t < first + n; t++)
        {
            d = (double)sample(t, i) - mean[i];
            s += d * d;
        }
        sd[i] = sqrt(s / (n - 1));
    }
}

/** \brief The volumes of a 4D file are the samples, as if they were separate files. */
static void
test_four_d_file_is_many_samples(void)
{
    CAT_VolAverageResult r4, r3;
    char *in4[1], *in3[NT];
    double mean[NVOX], sd[NVOX];
    int t, ok;

    make_fixture(0.0f);
    in4[0] = (char *)file_4d;
    for (t = 0; t < NT; t++)
        in3[t] = vol_names[t];

    ok = CAT_VolAverageScan(in4, 1, 0, &r4) == 0;
    MU_ASSERT("the 4D file is scanned", ok);
    if (ok)
    {
        MU_ASSERT("one sample per volume", r4.n_samples == NT && r4.n_volumes[0] == NT);
        MU_ASSERT("the matrix size is that of one volume",
                  r4.dims[0] == NX && r4.dims[1] == NY && r4.dims[2] == NZ &&
                  r4.n_voxels == NVOX);
        ok = CAT_VolAverageCompute(in4, 1, 1, 0, &r4) == 0;
        MU_ASSERT("mean and std of the 4D file", ok);
    }

    if (ok && CAT_VolAverageScan(in3, NT, 0, &r3) == 0 &&
        CAT_VolAverageCompute(in3, NT, 1, 0, &r3) == 0)
    {
        reference(0, NT, mean, sd);
        MU_ASSERT("the mean is taken over the volumes",
                  same_values(r4.mean, mean, NVOX, 1e-9));
        MU_ASSERT("the std is taken over the volumes",
                  same_values(r4.std, sd, NVOX, 1e-9));
        MU_ASSERT("the mean is that of the same volumes as separate files",
                  memcmp(r4.mean, r3.mean, sizeof(double) * NVOX) == 0);
        MU_ASSERT("the std is that of the same volumes as separate files",
                  memcmp(r4.std, r3.std, sizeof(double) * NVOX) == 0);
        CAT_VolAverageFree(&r3);
    }
    else if (ok)
        MU_ASSERT("the separate files are processed", 0);

    CAT_VolAverageFree(&r4);
    remove_fixture();
}

/** \brief 3D and 4D files can be mixed. */
static void
test_mixed_inputs(void)
{
    CAT_VolAverageResult mixed, sep;
    char *in_mixed[3], *in_sep[5];
    int t, ok;

    make_fixture(0.0f);
    in_mixed[0] = vol_names[0];
    in_mixed[1] = (char *)file_mid;
    in_mixed[2] = vol_names[4];
    for (t = 0; t < 5; t++)
        in_sep[t] = vol_names[t];

    ok = CAT_VolAverageScan(in_mixed, 3, 0, &mixed) == 0 &&
         CAT_VolAverageCompute(in_mixed, 3, 1, 0, &mixed) == 0;
    MU_ASSERT("the mixed inputs are processed", ok);
    if (ok)
    {
        MU_ASSERT("five samples in three files", mixed.n_samples == 5 &&
                  mixed.n_volumes[0] == 1 && mixed.n_volumes[1] == 3 &&
                  mixed.n_volumes[2] == 1);
        if (CAT_VolAverageScan(in_sep, 5, 0, &sep) == 0 &&
            CAT_VolAverageCompute(in_sep, 5, 1, 0, &sep) == 0)
        {
            MU_ASSERT("the mean equals that of five separate files",
                      memcmp(mixed.mean, sep.mean, sizeof(double) * NVOX) == 0);
            MU_ASSERT("the std equals that of five separate files",
                      memcmp(mixed.std, sep.std, sizeof(double) * NVOX) == 0);
            CAT_VolAverageFree(&sep);
        }
        else
            MU_ASSERT("the separate files are processed", 0);
    }

    CAT_VolAverageFree(&mixed);
    remove_fixture();
}

/** \brief A single sample has a mean but no standard deviation: say so, do not write NaN. */
static void
test_std_needs_two_volumes(void)
{
    CAT_VolAverageResult r;
    char *in[1];
    double v0[NVOX];
    int i;

    make_fixture(0.0f);
    in[0] = vol_names[0];
    for (i = 0; i < NVOX; i++)
        v0[i] = (double)sample(0, i);

    MU_ASSERT("a single volume is scanned", CAT_VolAverageScan(in, 1, 0, &r) == 0);
    MU_ASSERT("n = 1 has no standard deviation", CAT_VolAverageCompute(in, 1, 1, 0, &r) != 0);
    MU_ASSERT("n = 1 has a mean", CAT_VolAverageCompute(in, 1, 0, 0, &r) == 0 &&
              same_values(r.mean, v0, NVOX, 0.0));
    MU_ASSERT("and no std array", r.std == NULL);
    CAT_VolAverageFree(&r);
    remove_fixture();
}

/** \brief Inputs of a different matrix size are refused instead of read out of bounds. */
static void
test_matrix_size_mismatch_is_refused(void)
{
    CAT_VolAverageResult r;
    const char *other = "test_volaverage_other.nii";
    float data[(NX + 1) * NY * NZ];
    char *in[2];
    int i;

    make_fixture(0.0f);
    for (i = 0; i < (NX + 1) * NY * NZ; i++)
        data[i] = 1.0f;
    write_series(other, NX + 1, NY, NZ, 1, data);

    in[0] = vol_names[0];
    in[1] = (char *)other;
    MU_ASSERT("a different matrix size is refused", CAT_VolAverageScan(in, 2, 0, &r) != 0);
    MU_ASSERT("and nothing is left behind", r.header == NULL && r.n_volumes == NULL);

    remove(other);
    remove_fixture();
}

/**
 * \brief A voxel that is constant across the samples has std 0, not NaN.
 *
 * sum_x2 - n*mean^2 is 0 only up to rounding and about half of such voxels come
 * out slightly negative, which sqrt() turned into NaN.  100 volumes of 64
 * different constants.
 */
static void
test_constant_voxels_have_no_nan_std(void)
{
    enum { N = 64, VOLS = 100 };
    CAT_VolAverageResult r;
    const char *file = "test_volaverage_const.nii";
    float *data = (float *)malloc(sizeof(float) * N * VOLS);
    char *in[1];
    int t, i, bad = 0;

    for (t = 0; t < VOLS; t++)
        for (i = 0; i < N; i++)
            data[t * N + i] = 0.1f * (float)(i + 1);
    write_series(file, 4, 4, 4, VOLS, data);

    in[0] = (char *)file;
    if (CAT_VolAverageScan(in, 1, 0, &r) == 0 &&
        CAT_VolAverageCompute(in, 1, 1, 0, &r) == 0)
    {
        for (i = 0; i < N; i++)
            if (!(r.std[i] >= 0.0 && r.std[i] < 1e-5))
                bad++;
        MU_ASSERT("no NaN std where the samples are all equal", bad == 0);
        CAT_VolAverageFree(&r);
    }
    else
        MU_ASSERT("the constant series is processed", 0);

    free(data);
    remove(file);
}

/** \brief One z-score per volume, the same as for separate files; the outlier stands out. */
static void
test_zscores_per_volume(void)
{
    CAT_VolAverageResult r4, r3;
    char *in4[1], *in3[NT];
    double z4[NT], z3[NT], others = 0.0;
    int t, ok;

    make_fixture(20.0f);
    in4[0] = (char *)file_4d;
    for (t = 0; t < NT; t++)
        in3[t] = vol_names[t];

    ok = CAT_VolAverageScan(in4, 1, 0, &r4) == 0 &&
         CAT_VolAverageCompute(in4, 1, 1, 0, &r4) == 0 &&
         CAT_VolAverageZscore(in4, 1, &r4, 0, z4) == 0;
    MU_ASSERT("z-scores of the 4D file", ok);
    if (ok)
    {
        ok = CAT_VolAverageScan(in3, NT, 0, &r3) == 0 &&
             CAT_VolAverageCompute(in3, NT, 1, 0, &r3) == 0 &&
             CAT_VolAverageZscore(in3, NT, &r3, 0, z3) == 0;
        MU_ASSERT("z-scores of the separate files", ok);
        if (ok)
        {
            MU_ASSERT("one z-score per volume, in the order of the volumes",
                      same_values(z4, z3, NT, 1e-12));
            for (t = 0; t < NT - 1; t++)
                if (z4[t] > others)
                    others = z4[t];
            MU_ASSERT("the shifted volume has the largest z-score", z4[NT - 1] > 1.5 * others);
            CAT_VolAverageFree(&r3);
        }
    }

    CAT_VolAverageFree(&r4);
    remove_fixture();
}

/** \brief Z-scores need the standard deviation. */
static void
test_zscore_needs_std(void)
{
    CAT_VolAverageResult r;
    char *in[1];
    double z[NT];

    make_fixture(0.0f);
    in[0] = (char *)file_4d;
    if (CAT_VolAverageScan(in, 1, 0, &r) == 0 && CAT_VolAverageCompute(in, 1, 0, 0, &r) == 0)
        MU_ASSERT("no z-scores without the std", CAT_VolAverageZscore(in, 1, &r, 0, z) != 0);
    else
        MU_ASSERT("the mean is computed", 0);

    CAT_VolAverageFree(&r);
    remove_fixture();
}

/**
 * \brief NaN and Inf are read as 0 by the NIfTI library.
 *
 * A NaN-masked background (the z-maps of a voxel-wise model are NaN outside
 * its mask) is therefore 0 in the mean and in the std; the header documents it.
 */
static void
test_nonfinite_values_are_zero(void)
{
    CAT_VolAverageResult r;
    const char *file = "test_volaverage_nan.nii";
    float data[NT * NVOX];
    char *in[1];
    int t, i;

    for (t = 0; t < NT; t++)
        for (i = 0; i < NVOX; i++)
            data[t * NVOX + i] = sample(t, i);
    for (t = 0; t < NT; t++)
    {
        data[t * NVOX] = NAN;
        data[t * NVOX + 1] = INFINITY;
    }
    write_series(file, NX, NY, NZ, NT, data);

    in[0] = (char *)file;
    if (CAT_VolAverageScan(in, 1, 0, &r) == 0 && CAT_VolAverageCompute(in, 1, 1, 0, &r) == 0)
    {
        MU_ASSERT("NaN voxels have mean 0 and std 0", r.mean[0] == 0.0 && r.std[0] == 0.0);
        MU_ASSERT("Inf voxels have mean 0 and std 0", r.mean[1] == 0.0 && r.std[1] == 0.0);
        MU_ASSERT("the others are unaffected", r.mean[2] > 40.0 && r.std[2] > 0.0);
        CAT_VolAverageFree(&r);
    }
    else
        MU_ASSERT("the series with NaN is processed", 0);

    remove(file);
}

int main(void)
{
    MU_RUN_TEST(test_four_d_file_is_many_samples);
    MU_RUN_TEST(test_mixed_inputs);
    MU_RUN_TEST(test_std_needs_two_volumes);
    MU_RUN_TEST(test_matrix_size_mismatch_is_refused);
    MU_RUN_TEST(test_constant_voxels_have_no_nan_std);
    MU_RUN_TEST(test_zscores_per_volume);
    MU_RUN_TEST(test_zscore_needs_std);
    MU_RUN_TEST(test_nonfinite_values_are_zero);
    printf("%d tests run, %d failed\n", tests_run, tests_failed);
    return tests_failed ? 1 : 0;
}
