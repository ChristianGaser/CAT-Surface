#include "minunit.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "CAT_NiftiLib.h"
#include "CAT_BBReg.h"
#include "CAT_Realign.h"

/*
 * Motion correction of a 4D series in mcflirt's conventions.
 *
 * The conventions are checked against numbers mcflirt (FSL 6.0.7.23) wrote
 * for a synthetic series, since FSL is not available where the tests run.
 * The estimation is checked on an analytic phantom with known motion.  Its
 * accuracy is noise-limited at about 0.06 mm (0.1 mm with a separate, equally
 * noisy reference): mcflirt's RMS sphere of 80 mm is larger than the 72 mm
 * phantom and weights rotation errors heavily.  A failed registration is off
 * by more than a millimetre.
 */

/* ------------------------------------------------------------------------
 * mcflirt's own output for three volumes of a 60x72x60 series of 3 mm
 * voxels (radiological storage), registered with -reffile.  com is the
 * centre of mass of that reference, which CAT_VolRealignCentreOfMass()
 * reproduces; centre is the centre of its grid.
 * ------------------------------------------------------------------------ */
static const double fix_mat[3][16] = {
    {0.999996, 0.000832, -0.002572, 0.266024,
     -0.000831, 1.000000, 0.000257, 0.179145,
     0.002572, -0.000255, 0.999997, -0.065103,
     0, 0, 0, 1},
    {0.999998, 0.000416, -0.001935, 0.183538,
     -0.000414, 0.999999, 0.001474, 0.050212,
     0.001935, -0.001473, 0.999997, 0.083049,
     0, 0, 0, 1},
    {0.998756, -0.037452, -0.032911, 4.210960,
     0.038309, 0.998933, 0.025803, -6.327763,
     0.031909, -0.027032, 0.999125, -2.602680,
     0, 0, 0, 1}};
static const double fix_par[3][6] = {
    {0.000257089, 0.00257151, 0.000831878, 0.158934, 0.125018, 0.135275},
    {0.00147406, 0.00193471, 0.000416408, 0.0808293, 0.125251, 0.0979259},
    {0.0258201, 0.0329165, -0.0374808, -2.36614, -1.09306, -2.71078}};
static const double fix_abs[3] = {0.263863, 0.221221, 4.90578};
static const double fix_rel_1_0 = 0.108941; /* volume 6 against volume 5 */
static const double fix_com[3] = {88.527336, 106.044810, 75.825065};
static const double fix_centre[3] = {88.5, 106.5, 88.5};

/** \brief .par and RMS from mcflirt's matrices reproduce mcflirt's numbers. */
static void
test_mcflirt_conventions(void)
{
    static const double I[16] = {1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1};
    double par[6], dr = 0.0, dt = 0.0, da = 0.0;
    int k, i;

    for (k = 0; k < 3; k++)
    {
        CAT_VolRealignFslParams(fix_mat[k], fix_com, par);
        for (i = 0; i < 3; i++)
        {
            dr = fmax(dr, fabs(par[i] - fix_par[k][i]));
            dt = fmax(dt, fabs(par[3 + i] - fix_par[k][3 + i]));
        }
        da = fmax(da, fabs(CAT_VolRealignRms(fix_mat[k], I, fix_centre) - fix_abs[k]));
    }
    if (dr >= 2e-6 || dt >= 5e-4 || da >= 5e-4)
        printf("  conventions: rotation %.1e rad, translation %.1e mm, abs RMS %.1e mm\n",
               dr, dt, da);
    MU_ASSERT(".par rotations match mcflirt (R = Rx Ry Rz)", dr < 2e-6);
    MU_ASSERT(".par translations match mcflirt (about the centre of mass)", dt < 5e-4);
    MU_ASSERT("_abs.rms matches mcflirt (80 mm around the grid centre)", da < 5e-4);
    MU_ASSERT("_rel.rms matches mcflirt",
              fabs(CAT_VolRealignRms(fix_mat[1], fix_mat[0], fix_centre) - fix_rel_1_0) < 5e-4);
}

/* ------------------------------------------------------------------------
 * Analytic phantom: Gaussian blobs without symmetry on a grid of N^3 voxels
 * of VOX mm centred on the origin, plus noise (see tests/test_volumereg.c).
 * ------------------------------------------------------------------------ */
#define N 36
#define VOX 2.0
#define NT 6

static const double blob[][5] = {
    {0.0, 0.0, 0.0, 14.0, 60.0},
    {-12.0, 8.0, 4.0, 5.0, 80.0},
    {10.0, -6.0, 9.0, 4.0, -50.0},
    {6.0, 14.0, -8.0, 6.0, 70.0},
    {-9.0, -12.0, -10.0, 3.5, 90.0},
    {14.0, 4.0, -3.0, 3.0, 60.0},
    {-3.0, -2.0, 15.0, 4.5, -40.0},
};

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

/** \brief Header of the phantom grid; sign -1 stores x reversed (radiological). */
static nifti_image *
make_header(int sign)
{
    int dims[8] = {3, N, N, N, 1, 1, 1, 1};
    nifti_image *nii = nifti_make_new_nim(dims, DT_FLOAT32, 0);
    int i, j;

    for (i = 0; i < 4; i++)
        for (j = 0; j < 4; j++)
            nii->sto_xyz.m[i][j] = (i == j) ? (i < 3 ? VOX : 1.0f) : 0.0f;
    nii->sto_xyz.m[0][0] = (float)(sign * VOX);
    nii->sto_xyz.m[0][3] = (float)(sign * 0.5 * (N - 1) * VOX * -1.0);
    nii->sto_xyz.m[1][3] = (float)(-0.5 * (N - 1) * VOX);
    nii->sto_xyz.m[2][3] = (float)(-0.5 * (N - 1) * VOX);
    nii->sto_ijk = nifti_mat44_inverse(nii->sto_xyz);
    nii->sform_code = NIFTI_XFORM_SCANNER_ANAT;
    nii->dx = nii->dy = nii->dz = VOX;
    return nii;
}

/** \brief Sample gain * phantom(m^-1 y) + noise into vol; m NULL is the identity. */
static void
sample(const nifti_image *nii, const double *m, double gain, unsigned long seed,
       float *vol)
{
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
}

/* fixed -> moving motion of every volume: a drift with a jump at volume 4 */
static const CAT_RigidParams motion[NT] = {
    {0.3, -0.2, 0.1, 0.004, -0.003, 0.002},
    {0.5, -0.1, 0.2, 0.006, -0.002, 0.004},
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {0.6, 0.1, 0.4, 0.009, 0.001, 0.005},
    {1.8, -1.2, 1.0, 0.030, -0.020, 0.015},
    {1.9, -1.0, 1.1, 0.032, -0.019, 0.017},
};

/** \brief Largest RMS deviation (80 mm sphere at the origin) from the truth. */
static double
max_error(const CAT_VolRealignResult *res)
{
    static const double origin[3] = {0.0, 0.0, 0.0};
    double worst = 0.0;
    int t;

    for (t = 0; t < NT; t++)
    {
        double m[16], truth[16];
        CAT_BBReg_params_to_matrix(&motion[t], m);
        CAT_BBReg_invert_matrix(m, truth); /* volume -> reference */
        worst = fmax(worst, CAT_VolRealignRms(res->ras + 16 * t, truth, origin));
    }
    return worst;
}

/**
 * \brief A series registered to a brighter separate reference.
 *
 * The motion is recovered and the intensity scale reference = s * volume is
 * the gain, as with fMRIPrep's boldref from non-steady-state volumes.
 */
static void
test_series_to_reference(void)
{
    nifti_image *nii = make_header(-1);
    int dims[3] = {N, N, N}, t, ok = 1;
    size_t nvox = (size_t)N * N * N;
    float *series = (float *)malloc(sizeof(float) * nvox * NT);
    float *ref = (float *)malloc(sizeof(float) * nvox);
    CAT_VolRealignResult res;
    double err;

    for (t = 0; t < NT; t++)
    {
        double m[16];
        CAT_BBReg_params_to_matrix(&motion[t], m);
        sample(nii, m, 1.0, 10 + t, series + t * nvox);
    }
    sample(nii, NULL, 1.3, 99, ref);

    MU_ASSERT("estimation succeeds",
              CAT_VolRealignEstimate(series, nii, dims, NT, ref, nii, dims, NULL, &res) == 0);
    err = max_error(&res);
    if (err >= 0.2)
        printf("  separate reference: largest error %.3f mm\n", err);
    MU_ASSERT("separate reference: motion within 0.2 mm RMS", err < 0.2);
    for (t = 0; t < NT; t++)
        if (fabs(res.scale[t] - 1.3) > 0.02)
            ok = 0;
    MU_ASSERT("separate reference: the intensity scale is the gain", ok);
    MU_ASSERT("separate reference: no volume is the reference", res.ref_index == -1);

    CAT_VolRealignFree(&res);
    free(series);
    free(ref);
    nifti_image_free(nii);
}

/**
 * \brief A series registered to one of its volumes, stored neurologically.
 *
 * The reference volume is the identity exactly, as mcflirt writes it, and the
 * resampled reference volume is the volume itself.
 */
static void
test_series_to_own_volume(void)
{
    nifti_image *nii = make_header(1);
    int dims[3] = {N, N, N}, t, i, ok = 1, same = 1;
    size_t nvox = (size_t)N * N * N;
    float *series = (float *)malloc(sizeof(float) * nvox * NT);
    float *out = (float *)malloc(sizeof(float) * nvox * NT);
    CAT_VolRealignOpts opts;
    CAT_VolRealignResult res;
    double err;

    for (t = 0; t < NT; t++)
    {
        double m[16];
        CAT_BBReg_params_to_matrix(&motion[t], m);
        sample(nii, m, 1.0, 20 + t, series + t * nvox);
    }
    CAT_VolRealignOptionsInit(&opts);
    MU_ASSERT("default reference is volume nt/2", opts.refvol == -1);
    opts.refvol = 2; /* the volume without motion */

    MU_ASSERT("estimation succeeds",
              CAT_VolRealignEstimate(series, nii, dims, NT, NULL, NULL, NULL, &opts, &res) == 0);
    MU_ASSERT("own volume: volume 2 is the reference", res.ref_index == 2);
    for (i = 0; i < 16; i++)
        if (res.fsl[16 * 2 + i] != ((i % 5 == 0) ? 1.0 : 0.0))
            ok = 0;
    MU_ASSERT("own volume: its matrix is exactly the identity", ok);
    MU_ASSERT("own volume: its RMS and parameters are 0",
              res.rms_abs[2] == 0.0 && res.par[12] == 0.0 && res.par[17] == 0.0);
    err = max_error(&res);
    if (err >= 0.2)
        printf("  own volume: largest error %.3f mm\n", err);
    MU_ASSERT("own volume: motion within 0.2 mm RMS", err < 0.2);
    MU_ASSERT("own volume: nt - 1 relative RMS values average to the mean",
              fabs((res.rms_rel[0] + res.rms_rel[1] + res.rms_rel[2] + res.rms_rel[3] +
                    res.rms_rel[4]) / 5.0 - res.rms_rel_mean) < 1e-12);

    MU_ASSERT("resampling succeeds",
              CAT_VolRealignResample(series, nii, dims, NT, nii, dims, &res, out) == 0);
    for (i = 0; i < (int)nvox; i++)
        if (out[2 * nvox + i] != series[2 * nvox + i])
            same = 0;
    MU_ASSERT("own volume: the resampled reference is the volume itself", same);

    CAT_VolRealignFree(&res);
    free(series);
    free(out);
    nifti_image_free(nii);
}

/* ------------------------------------------------------------------------
 * A spherically symmetric object on a grid of NS^3 voxels of VOX mm whose
 * world origin is the corner voxel, 42 mm from the object, as a scanner
 * origin away from the head is.
 * ------------------------------------------------------------------------ */
#define NS 24
#define NS_T 4

/** \brief Gaussian blob of 100 on 10 at the grid centre, moved by dx mm in x, plus noise. */
static void
sample_sphere(double dx, unsigned long seed, float *vol)
{
    unsigned long state = seed;
    int x, y, z;

    for (z = 0; z < NS; z++)
        for (y = 0; y < NS; y++)
            for (x = 0; x < NS; x++)
            {
                double ux = x - 0.5 * NS - dx / VOX, uy = y - 0.5 * NS, uz = z - 0.5 * NS;
                vol[x + NS * (y + NS * z)] =
                    (float)(100.0 * exp(-(ux * ux + uy * uy + uz * uz) / 30.0) + 10.0 +
                            gauss(&state));
            }
}

/**
 * \brief A symmetric object is registered without running away.
 *
 * Rotations about the centre of the object leave the image unchanged, so the
 * data do not determine them.  The undamped Gauss-Newton step moved freely
 * along such directions: here the motion-free volumes, which start from the
 * estimate of the moved one, ended with their object 88 mm from the
 * reference's.  The rotation stays undetermined -- noise sets it, for mcflirt
 * as well -- so only where the object lands is checked.
 */
static void
test_symmetric_object(void)
{
    int dims8[8] = {3, NS, NS, NS, 1, 1, 1, 1}, dims[3] = {NS, NS, NS}, i, j, t;
    nifti_image *nii = nifti_make_new_nim(dims8, DT_FLOAT32, 0);
    size_t nvox = (size_t)NS * NS * NS;
    float *series = (float *)malloc(sizeof(float) * nvox * NS_T);
    float *ref = (float *)malloc(sizeof(float) * nvox);
    double centre = 0.5 * NS * VOX, worst = 0.0;
    CAT_VolRealignResult res;

    for (i = 0; i < 4; i++)
        for (j = 0; j < 4; j++)
            nii->sto_xyz.m[i][j] = (i == j) ? (i < 3 ? (float)VOX : 1.0f) : 0.0f;
    nii->sto_ijk = nifti_mat44_inverse(nii->sto_xyz);
    nii->sform_code = NIFTI_XFORM_SCANNER_ANAT;
    nii->dx = nii->dy = nii->dz = VOX;

    /* volume 2, the first registered, moves by one voxel */
    sample_sphere(0.0, 40, ref);
    for (t = 0; t < NS_T; t++)
        sample_sphere(t == 2 ? VOX : 0.0, 41 + t, series + t * nvox);

    MU_ASSERT("symmetric object: estimation succeeds",
              CAT_VolRealignEstimate(series, nii, dims, NS_T, ref, nii, dims, NULL, &res) == 0);

    /* volume -> reference has to map the object's centre onto the reference's */
    for (t = 0; t < NS_T; t++)
    {
        const double c[3] = {centre + (t == 2 ? VOX : 0.0), centre, centre};
        const double *m = res.ras + 16 * t;
        double d = 0.0;

        for (i = 0; i < 3; i++)
        {
            double e = m[4 * i] * c[0] + m[4 * i + 1] * c[1] + m[4 * i + 2] * c[2] +
                       m[4 * i + 3] - centre;
            d += e * e;
        }
        MU_ASSERT("symmetric object: the mapped centre is finite", isfinite(d));
        worst = fmax(worst, sqrt(d));
    }
    if (worst >= 0.5)
        printf("  symmetric object: lands %.3f mm from the reference's\n", worst);
    MU_ASSERT("symmetric object: every object lands within 0.5 mm", worst < 0.5);

    CAT_VolRealignFree(&res);
    free(series);
    free(ref);
    nifti_image_free(nii);
}

/** \brief The centre of mass subtracts the minimum and runs radiologically. */
static void
test_centre_of_mass(void)
{
    nifti_image *rad = make_header(-1), *neu = make_header(1);
    int dims[3] = {N, N, N};
    size_t nvox = (size_t)N * N * N, i;
    float *vol = (float *)calloc(nvox, sizeof(float));
    double com[3];

    /* one bright voxel at index (5, 7, 9) on a background of 100 */
    for (i = 0; i < nvox; i++)
        vol[i] = 100.0f;
    vol[5 + N * (7 + N * 9)] = 200.0f;

    CAT_VolRealignCentreOfMass(vol, rad, dims, com);
    MU_ASSERT("centre of mass: the background is subtracted",
              fabs(com[0] - 5 * VOX) < 1e-9 && fabs(com[1] - 7 * VOX) < 1e-9 &&
                  fabs(com[2] - 9 * VOX) < 1e-9);
    CAT_VolRealignCentreOfMass(vol, neu, dims, com);
    MU_ASSERT("centre of mass: x is reversed for neurological storage",
              fabs(com[0] - (N - 1 - 5) * VOX) < 1e-9 && fabs(com[1] - 7 * VOX) < 1e-9);

    free(vol);
    nifti_image_free(rad);
    nifti_image_free(neu);
}

/** \brief Count the lines of a text file; -1 if it cannot be read. */
static int
count_lines(const char *file)
{
    FILE *fp = fopen(file, "r");
    int n = 0, c;

    if (!fp)
        return -1;
    while ((c = fgetc(fp)) != EOF)
        if (c == '\n')
            n++;
    fclose(fp);
    return n;
}

/** \brief Output names append to the whole base, as mcflirt's do. */
static void
test_outputs(void)
{
    CAT_VolRealignResult res;
    char image[256];
    int t;

    CAT_VolRealignImageName("bold_mcf.nii.gz", image, sizeof(image));
    MU_ASSERT("a .nii.gz base is the image name", strcmp(image, "bold_mcf.nii.gz") == 0);
    CAT_VolRealignImageName("bold_mcf.nii", image, sizeof(image));
    MU_ASSERT("a .nii base is the image name", strcmp(image, "bold_mcf.nii") == 0);

    memset(&res, 0, sizeof(res));
    res.nt = 3;
    res.fsl = (double *)calloc(16 * 3, sizeof(double));
    res.par = (double *)calloc(6 * 3, sizeof(double));
    res.rms_abs = (double *)calloc(3, sizeof(double));
    res.rms_rel = (double *)calloc(2, sizeof(double));
    for (t = 0; t < 3; t++)
        res.fsl[16 * t] = res.fsl[16 * t + 5] = res.fsl[16 * t + 10] =
            res.fsl[16 * t + 15] = 1.0;

    MU_ASSERT("writing succeeds",
              CAT_VolRealignWrite("test_realign.nii.gz", &res, 1, 1, 1, 1) == 0);
    MU_ASSERT("MAT_0002 has four lines",
              count_lines("test_realign.nii.gz.mat/MAT_0002") == 4);
    MU_ASSERT(".par has one line per volume", count_lines("test_realign.nii.gz.par") == 3);
    MU_ASSERT("_abs.rms has nt lines", count_lines("test_realign.nii.gz_abs.rms") == 3);
    MU_ASSERT("_rel.rms has nt - 1 lines", count_lines("test_realign.nii.gz_rel.rms") == 2);
    MU_ASSERT("the means are written",
              count_lines("test_realign.nii.gz_abs_mean.rms") == 1 &&
                  count_lines("test_realign.nii.gz_rel_mean.rms") == 1);
    MU_ASSERT("an existing .mat directory is reused",
              CAT_VolRealignWrite("test_realign.nii.gz", &res, 1, 0, 0, 0) == 0);

    for (t = 0; t < 3; t++)
    {
        char name[64];
        snprintf(name, sizeof(name), "test_realign.nii.gz.mat/MAT_%04d", t);
        remove(name);
    }
    remove("test_realign.nii.gz.mat");
    remove("test_realign.nii.gz.par");
    remove("test_realign.nii.gz_abs.rms");
    remove("test_realign.nii.gz_rel.rms");
    remove("test_realign.nii.gz_abs_mean.rms");
    remove("test_realign.nii.gz_rel_mean.rms");
    CAT_VolRealignFree(&res);
}

int main(void)
{
    MU_RUN_TEST(test_mcflirt_conventions);
    MU_RUN_TEST(test_centre_of_mass);
    MU_RUN_TEST(test_series_to_reference);
    MU_RUN_TEST(test_series_to_own_volume);
    MU_RUN_TEST(test_symmetric_object);
    MU_RUN_TEST(test_outputs);
    printf("%d tests run, %d failed\n", tests_run, tests_failed);
    return tests_failed ? 1 : 0;
}
