#include "minunit.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "CAT_NiftiLib.h"
#include "CAT_BBReg.h"
#include "CAT_VolumeReg.h"

/*
 * Robust rigid registration on an analytic phantom.
 *
 * The phantom is a sum of Gaussian blobs of different size and contrast,
 * arranged without symmetry so that every rotation is identifiable.  Both
 * volumes are sampled from the formula, the moving one at the transformed
 * positions, so no interpolation of the test data enters the comparison.
 *
 * Both get noise, as every scan has.  Without it most residuals are those of
 * the flat background, exactly zero, and the MAD that scales the Tukey weights
 * collapses: every voxel that carries structure is then rejected as an
 * outlier.
 */

#define N 36
#define VOX 2.0

static const double blob[][5] = {
    /* x, y, z (mm), sigma (mm), amplitude */
    {0.0, 0.0, 0.0, 14.0, 60.0},
    {-12.0, 8.0, 4.0, 5.0, 80.0},
    {10.0, -6.0, 9.0, 4.0, -50.0},
    {6.0, 14.0, -8.0, 6.0, 70.0},
    {-9.0, -12.0, -10.0, 3.5, 90.0},
    {14.0, 4.0, -3.0, 3.0, 60.0},
    {-3.0, -2.0, 15.0, 4.5, -40.0},
};

/** \brief Phantom intensity at world position (x, y, z). */
static double
phantom(double x, double y, double z)
{
    double v = 10.0;
    size_t k;

    for (k = 0; k < sizeof(blob) / sizeof(blob[0]); k++)
    {
        double dx = x - blob[k][0], dy = y - blob[k][1], dz = z - blob[k][2];
        double s = blob[k][3];
        v += blob[k][4] * exp(-0.5 * (dx * dx + dy * dy + dz * dz) / (s * s));
    }
    return v;
}

/** \brief Header of an N^3 grid of VOX mm voxels centred on the origin. */
static nifti_image *
make_header(void)
{
    int dims[8] = {3, N, N, N, 1, 1, 1, 1};
    nifti_image *nii = nifti_make_new_nim(dims, DT_FLOAT32, 0);
    int i, j;

    for (i = 0; i < 4; i++)
        for (j = 0; j < 4; j++)
            nii->sto_xyz.m[i][j] = (i == j) ? (i < 3 ? VOX : 1.0f) : 0.0f;
    for (i = 0; i < 3; i++)
        nii->sto_xyz.m[i][3] = (float)(-0.5 * (N - 1) * VOX);
    nii->sto_ijk = nifti_mat44_inverse(nii->sto_xyz);
    nii->sform_code = NIFTI_XFORM_SCANNER_ANAT;
    nii->dx = nii->dy = nii->dz = VOX;
    return nii;
}

/** \brief Standard normal deviate from a fixed-seed generator (Box-Muller). */
static double
gauss(unsigned long *state)
{
    double u1, u2;

    *state = *state * 6364136223846793005UL + 1442695040888963407UL;
    u1 = ((*state >> 11) + 0.5) / 9007199254740992.0;
    *state = *state * 6364136223846793005UL + 1442695040888963407UL;
    u2 = ((*state >> 11) + 0.5) / 9007199254740992.0;
    return sqrt(-2.0 * log(u1)) * cos(6.283185307179586 * u2);
}

/**
 * \brief Sample the phantom on the grid of nii, times gain, plus noise.
 *
 * With m (fixed -> moving, row-major) the volume is the moving image of a
 * registration whose exact answer is m: voxel y holds phantom(m^-1 y).  With
 * m NULL it is the fixed image.  The noise has a standard deviation of 1
 * (the blobs reach 40 to 90) and its own seed per volume.
 */
static float *
sample(const nifti_image *nii, const double *m, double gain, unsigned long seed)
{
    float *vol = (float *)malloc(sizeof(float) * N * N * N);
    double inv[16];
    unsigned long state = seed;
    int x, y, z;

    if (m)
        CAT_BBReg_invert_matrix(m, inv);
    for (z = 0; z < N; z++)
        for (y = 0; y < N; y++)
            for (x = 0; x < N; x++)
            {
                double w[3], p[3];
                int r;
                for (r = 0; r < 3; r++)
                    w[r] = nii->sto_xyz.m[r][0] * x + nii->sto_xyz.m[r][1] * y +
                           nii->sto_xyz.m[r][2] * z + nii->sto_xyz.m[r][3];
                for (r = 0; r < 3; r++)
                    p[r] = m ? inv[4 * r] * w[0] + inv[4 * r + 1] * w[1] +
                                   inv[4 * r + 2] * w[2] + inv[4 * r + 3]
                             : w[r];
                vol[x + N * (y + N * z)] =
                    (float)(gain * phantom(p[0], p[1], p[2]) + gauss(&state));
            }
    return vol;
}

/** \brief Largest translation (mm) and rotation (rad) difference. */
static void
param_error(const CAT_RigidParams *a, const CAT_RigidParams *b,
            double *dt, double *dr)
{
    *dt = fmax(fabs(a->tx - b->tx), fmax(fabs(a->ty - b->ty), fabs(a->tz - b->tz)));
    *dr = fmax(fabs(a->rx - b->rx), fmax(fabs(a->ry - b->ry), fabs(a->rz - b->rz)));
}

static const CAT_RigidParams truth = {1.2, -0.8, 0.5, 0.035, -0.026, 0.017};

/**
 * \brief Registering prepared pyramids is the pairwise registration.
 *
 * Also checks that the known motion is recovered.
 */
static void
test_pyramids_match_pairwise(void)
{
    nifti_image *nii = make_header();
    int dims[3] = {N, N, N};
    double m[16], dt, dr;
    float *fixed, *moving;
    CAT_RigidParams p1 = {0, 0, 0, 0, 0, 0}, p2 = {0, 0, 0, 0, 0, 0};
    CAT_VolPyramid *fp, *mp;

    CAT_BBReg_params_to_matrix(&truth, m);
    fixed = sample(nii, NULL, 1.0, 1);
    moving = sample(nii, m, 1.0, 2);

    CAT_VolumeReg_register(fixed, nii, dims, moving, nii, dims, &p1,
                           3, 4.685, 20, 0);

    fp = CAT_VolumeReg_pyramid(fixed, nii, dims, 3);
    mp = CAT_VolumeReg_pyramid(moving, nii, dims, 3);
    MU_ASSERT("pyramids are built", fp != NULL && mp != NULL);
    if (fp && mp)
        CAT_VolumeReg_register_pyramids(fp, mp, &p2, 4.685, 20, 0);
    CAT_VolumeReg_pyramid_free(fp);
    CAT_VolumeReg_pyramid_free(mp);

    MU_ASSERT("pyramids: the same estimate as the pairwise call",
              memcmp(&p1, &p2, sizeof(p1)) == 0);

    param_error(&p1, &truth, &dt, &dr);
    if (dt >= 0.05 || dr >= 0.002)
        printf("  recovered t %.4f %.4f %.4f r %.5f %.5f %.5f\n",
               p1.tx, p1.ty, p1.tz, p1.rx, p1.ry, p1.rz);
    MU_ASSERT("the translation is recovered within 0.05 mm", dt < 0.05);
    MU_ASSERT("the rotation is recovered within 0.002 rad", dr < 0.002);

    free(fixed);
    free(moving);
    nifti_image_free(nii);
}

int main(void)
{
    MU_RUN_TEST(test_pyramids_match_pairwise);
    printf("%d tests run, %d failed\n", tests_run, tests_failed);
    return tests_failed ? 1 : 0;
}
