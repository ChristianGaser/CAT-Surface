#include "minunit.h"
#include <stdio.h>
#include "CAT_ThicknessQC.h"
#include "CAT_Vol.h"

#include <stdlib.h>
#include <string.h>
#include <math.h>

#define N 40
#define IDX(x, y, z) ((x) + (y) * N + (z) * N * N)

static int dims3[3] = {N, N, N};
static double vx1[3] = {1.0, 1.0, 1.0};

/* ------------------------------------------------------------------ */
/* defaults                                                            */
/* ------------------------------------------------------------------ */

static void test_qc_defaults(void)
{
    CAT_ThicknessQCOpts opts;

    CAT_ThicknessQCOptionsInit(&opts);

    /* cortex tops out near 4.5 mm; anything above that wants explaining */
    MU_ASSERT("threshold is an implausible thickness", opts.thresh >= 4.0);
    /* two banks back to back span at most ~5 mm, so half of that bounds a plate */
    MU_ASSERT("a plate radius is bounded by half a two-bank band",
              opts.plate_radius > 0.0 && opts.plate_radius <= 2.6);
    MU_ASSERT("connectivity is a legal value",
              opts.conn == 6 || opts.conn == 18 || opts.conn == 26);
}

/* ------------------------------------------------------------------ */
/* a plate and a ball are told apart by the inscribed radius            */
/* ------------------------------------------------------------------ */

static void test_qc_plate_vs_solid(void)
{
    const int nvox = N * N * N;
    float *gmt = (float *)malloc(sizeof(float) * nvox);
    float *cls = (float *)malloc(sizeof(float) * nvox);
    CAT_ThicknessQCOpts opts;
    CAT_ThicknessComponent *comps = NULL;
    int n_comps = 0;
    int x, y, z, i, n_plate = 0, n_solid = 0;
    double r_plate = -1.0, r_solid = -1.0;

    MU_ASSERT("alloc", gmt && cls);

    for (i = 0; i < nvox; i++)
        gmt[i] = 1.0f;

    /* a 4-voxel-thick plate: a glued sulcus seen end on.  Its extent along the
       sulcus is large, which is exactly why extent cannot be the criterion. */
    for (z = 4; z < 36; z++)
        for (y = 4; y < 36; y++)
            for (x = 6; x < 10; x++)
                gmt[IDX(x, y, z)] = 6.0f;

    /* a ball of radius 6: a subcortical mass.  Same label, same thickness, and
       it even holds fewer voxels than the plate. */
    for (z = 0; z < N; z++)
        for (y = 0; y < N; y++)
            for (x = 0; x < N; x++)
            {
                const double dx = x - 26.0, dy = y - 20.0, dz = z - 20.0;
                if (dx * dx + dy * dy + dz * dz <= 36.0)
                    gmt[IDX(x, y, z)] = 6.0f;
            }

    CAT_ThicknessQCOptionsInit(&opts);
    opts.thresh = 4.5;
    opts.plate_radius = 2.5;
    opts.min_volume = 20.0;

    MU_ASSERT("qc runs",
              CAT_VolThicknessQC(gmt, NULL, dims3, vx1, &opts, &comps, &n_comps,
                                 cls) == 0);
    MU_ASSERT("both structures were found", n_comps == 2);

    for (i = 0; i < n_comps; i++)
    {
        if (comps[i].shape == CAT_QC_PLATE)
        {
            n_plate++;
            r_plate = comps[i].max_radius;
        }
        else
        {
            n_solid++;
            r_solid = comps[i].max_radius;
        }
    }

    MU_ASSERT("the plate is classified as a plate", n_plate == 1);
    MU_ASSERT("the ball is classified as solid", n_solid == 1);

    /* a 4-voxel plate admits a sphere of radius 2, a radius-6 ball one of 6 */
    MU_ASSERT("the plate radius is half its thickness",
              r_plate > 1.5 && r_plate <= 2.5);
    MU_ASSERT("the ball radius is its own radius", r_solid > 4.0);

    /* the discriminator must not be size: here the plate is the larger object */
    {
        long n_plate_vox = 0, n_solid_vox = 0;
        for (i = 0; i < n_comps; i++)
        {
            if (comps[i].shape == CAT_QC_PLATE)
                n_plate_vox = comps[i].n_voxels;
            else
                n_solid_vox = comps[i].n_voxels;
        }
        MU_ASSERT("the plate outweighs the ball, so volume cannot be the test",
                  n_plate_vox > n_solid_vox);
    }

    /* the class map must agree with the report */
    MU_ASSERT("plate voxels are tagged 1",
              fabs(cls[IDX(7, 20, 20)] - (float)CAT_QC_PLATE) < 1e-6);
    MU_ASSERT("ball voxels are tagged 2",
              fabs(cls[IDX(26, 20, 20)] - (float)CAT_QC_SOLID) < 1e-6);
    MU_ASSERT("thin cortex is left untagged", cls[IDX(20, 2, 2)] == 0.0f);

    free(gmt);
    free(cls);
    free(comps);
}

/* ------------------------------------------------------------------ */
/* the label mask and the volume floor                                  */
/* ------------------------------------------------------------------ */

static void test_qc_mask_and_floor(void)
{
    const int nvox = N * N * N;
    float *gmt = (float *)malloc(sizeof(float) * nvox);
    float *lab = (float *)malloc(sizeof(float) * nvox);
    CAT_ThicknessQCOpts opts;
    CAT_ThicknessComponent *comps = NULL;
    int n_comps = 0;
    int x, y, z, i;

    MU_ASSERT("alloc", gmt && lab);

    for (i = 0; i < nvox; i++)
    {
        gmt[i] = 1.0f;
        lab[i] = 2.0f; /* GM everywhere */
    }

    /* one thick blob inside the ribbon, one outside it */
    for (z = 16; z < 24; z++)
        for (y = 16; y < 24; y++)
            for (x = 6; x < 14; x++)
                gmt[IDX(x, y, z)] = 6.0f;
    for (z = 16; z < 24; z++)
        for (y = 16; y < 24; y++)
            for (x = 26; x < 34; x++)
            {
                gmt[IDX(x, y, z)] = 6.0f;
                lab[IDX(x, y, z)] = 3.0f; /* WM: outside the cortical band */
            }

    CAT_ThicknessQCOptionsInit(&opts);
    opts.thresh = 4.5;

    MU_ASSERT("qc runs unmasked",
              CAT_VolThicknessQC(gmt, NULL, dims3, vx1, &opts, &comps, &n_comps,
                                 NULL) == 0);
    MU_ASSERT("without a label map both blobs are reported", n_comps == 2);
    free(comps);
    comps = NULL;

    MU_ASSERT("qc runs masked",
              CAT_VolThicknessQC(gmt, lab, dims3, vx1, &opts, &comps, &n_comps,
                                 NULL) == 0);
    MU_ASSERT("the label map confines the search to the ribbon", n_comps == 1);
    free(comps);
    comps = NULL;

    /* the floor drops anything smaller than it, measured as a volume */
    opts.min_volume = 1e9;
    MU_ASSERT("qc runs with a high floor",
              CAT_VolThicknessQC(gmt, NULL, dims3, vx1, &opts, &comps, &n_comps,
                                 NULL) == 0);
    MU_ASSERT("a floor above every component reports nothing", n_comps == 0);

    free(gmt);
    free(lab);
    free(comps);
}

/* ------------------------------------------------------------------ */
/* nothing implausible means nothing reported                           */
/* ------------------------------------------------------------------ */

static void test_qc_clean_map_is_silent(void)
{
    const int nvox = N * N * N;
    float *gmt = (float *)malloc(sizeof(float) * nvox);
    float *cls = (float *)malloc(sizeof(float) * nvox);
    CAT_ThicknessComponent *comps = NULL;
    int n_comps = -1;
    int i;

    MU_ASSERT("alloc", gmt && cls);

    /* a plausible ribbon everywhere: 2.5 mm */
    for (i = 0; i < nvox; i++)
        gmt[i] = 2.5f;

    MU_ASSERT("qc runs",
              CAT_VolThicknessQC(gmt, NULL, dims3, vx1, NULL, &comps, &n_comps,
                                 cls) == 0);
    MU_ASSERT("a plausible map yields no components", n_comps == 0);
    for (i = 0; i < nvox; i++)
        MU_ASSERT("and an empty class map", cls[i] == 0.0f);

    free(gmt);
    free(cls);
    free(comps);
}

int main(void)
{
    MU_RUN_TEST(test_qc_defaults);
    MU_RUN_TEST(test_qc_plate_vs_solid);
    MU_RUN_TEST(test_qc_mask_and_floor);
    MU_RUN_TEST(test_qc_clean_map_is_silent);
    printf("%d tests run, %d failed\n", tests_run, tests_failed);
    return tests_failed ? 1 : 0;
}
