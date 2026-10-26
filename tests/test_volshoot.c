#include "minunit.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "CAT_GeodesicShoot.h"

/*
 * Geodesic shooting on an analytic phantom.
 *
 * The phantom is a brain-like ellipsoid: a white matter core whose surface is
 * folded by a few sinusoids, wrapped in a grey matter shell of constant
 * thickness.  The five templates are the phantom's tissue maps smoothed from
 * coarse (Template_0) to fine (Template_4), as CAT's shooting templates are.
 * Tissue maps are evaluated from the formula at every position, so no
 * interpolation of test data enters the comparisons.
 */

#define NX 45
#define NY 53
#define NZ 45
#define VOX 1.5

static const size_t NV = (size_t)NX * NY * NZ;

/** \brief Smooth step from 1 (inside, d < 0) to 0 (outside), width w mm. */
static double
inside(double d, double w)
{
    return 0.5 * (1.0 - tanh(d / w));
}

/** \brief GM and WM probability at world position (x, y, z) in mm. */
static void
phantom(double x, double y, double z, double ventricle, double *gm, double *wm)
{
    double r = sqrt((x * x) / (25.0 * 25.0) + (y * y) / (31.0 * 31.0) + (z * z) / (25.0 * 25.0));
    double fold = 2.0 * sin(0.35 * x) * cos(0.3 * y) + 1.5 * sin(0.25 * z + 0.4 * y);
    double dwm = r * 25.0 - (16.0 + fold);  /* folded core */
    double dbr = r * 25.0 - 22.0;           /* brain surface */
    double pw = inside(dwm, 0.8), pb = inside(dbr, 0.8);
    double pv = 0.0;

    if (ventricle > 0.0)
    {
        double rv = sqrt(x * x / (ventricle * ventricle) + y * y / (1.6 * ventricle * 1.6 * ventricle) +
                         z * z / (ventricle * ventricle));
        pv = inside((rv - 1.0) * ventricle, 0.8);
    }
    *wm = pw * (1.0 - pv);
    *gm = (pb - pw > 0.0) ? pb - pw : 0.0;
}

/** \brief NIfTI-style matrix of the grid, centred on the origin. */
static void
grid_matrix(double M[16])
{
    memset(M, 0, sizeof(double) * 16);
    M[0] = M[5] = M[10] = VOX;
    M[3] = -VOX * (NX - 1) / 2.0;
    M[7] = -VOX * (NY - 1) / 2.0;
    M[11] = -VOX * (NZ - 1) / 2.0;
    M[15] = 1.0;
}

/** \brief Smooth displacement (mm) used to deform the input. */
static void
displacement(double x, double y, double z, double amp, double d[3])
{
    d[0] = amp * sin(0.08 * y + 0.3) * cos(0.07 * z);
    d[1] = amp * 0.8 * sin(0.07 * x) * cos(0.09 * z + 0.2);
    d[2] = amp * 0.6 * cos(0.08 * x + 0.1) * sin(0.06 * y);
}

/**
 * \brief Tissue maps on the grid, at the positions moved by -displacement,
 *        so that structures move by about +displacement.
 */
static void
make_maps(float *gm, float *wm, double amp, double ventricle)
{
    double M[16];
    int i, j, k;

    grid_matrix(M);
    for (k = 0; k < NZ; k++)
        for (j = 0; j < NY; j++)
            for (i = 0; i < NX; i++)
            {
                size_t v = i + (size_t)NX * (j + (size_t)NY * k);
                double x = M[0] * i + M[3], y = M[5] * j + M[7], z = M[10] * k + M[11];
                double d[3], g, w;
                displacement(x, y, z, amp, d);
                phantom(x - d[0], y - d[1], z - d[2], ventricle, &g, &w);
                gm[v] = (float)g;
                wm[v] = (float)w;
            }
}

/** \brief Box smoothing with radius r voxels along each axis, in place. */
static void
box_smooth(float *v, int r)
{
    float *t = (float *)malloc(sizeof(float) * NV);
    int d[3] = {NX, NY, NZ}, ax, i, j, k, o;

    for (ax = 0; ax < 3; ax++)
    {
        memcpy(t, v, sizeof(float) * NV);
        for (k = 0; k < NZ; k++)
            for (j = 0; j < NY; j++)
                for (i = 0; i < NX; i++)
                {
                    int c[3] = {i, j, k};
                    double s = 0.0;
                    int n = 0;
                    for (o = -r; o <= r; o++)
                    {
                        int cc[3] = {c[0], c[1], c[2]};
                        cc[ax] += o;
                        if (cc[ax] < 0 || cc[ax] >= d[ax]) continue;
                        s += t[cc[0] + (size_t)NX * (cc[1] + (size_t)NY * cc[2])];
                        n++;
                    }
                    v[i + (size_t)NX * (j + (size_t)NY * k)] = (float)(s / n);
                }
    }
    free(t);
}

/** \brief The five templates (GM and WM each), smoothed coarse to fine. */
static float **
make_templates(void)
{
    float **tpl = (float **)malloc(sizeof(float *) * CAT_SHOOT_NTEMPLATES);
    int l;

    for (l = 0; l < CAT_SHOOT_NTEMPLATES; l++)
    {
        tpl[l] = (float *)malloc(sizeof(float) * 2 * NV);
        make_maps(tpl[l], tpl[l] + NV, 0.0, 0.0);
        if (l < CAT_SHOOT_NTEMPLATES - 1)
        {
            box_smooth(tpl[l], CAT_SHOOT_NTEMPLATES - 1 - l);
            box_smooth(tpl[l] + NV, CAT_SHOOT_NTEMPLATES - 1 - l);
        }
    }
    return tpl;
}

/** \brief Trilinear value of vol at the 0-based voxel position (x, y, z). */
static double
trilinear(const float *vol, double x, double y, double z)
{
    int i = (int)floor(x), j = (int)floor(y), k = (int)floor(z), a, b, c;
    double fx = x - i, fy = y - j, fz = z - k, s = 0.0;

    for (c = 0; c < 2; c++)
        for (b = 0; b < 2; b++)
            for (a = 0; a < 2; a++)
            {
                int ii = i + a, jj = j + b, kk = k + c;
                double w = (a ? fx : 1 - fx) * (b ? fy : 1 - fy) * (c ? fz : 1 - fz);
                if (ii < 0) ii = 0;
                if (jj < 0) jj = 0;
                if (kk < 0) kk = 0;
                if (ii >= NX) ii = NX - 1;
                if (jj >= NY) jj = NY - 1;
                if (kk >= NZ) kk = NZ - 1;
                s += w * vol[ii + (size_t)NX * (jj + (size_t)NY * kk)];
            }
    return s;
}

/** \brief Mean squared GM+WM difference between the template and the pulled maps. */
static double
mismatch(const float *gm, const float *wm, const float *tpl, const float *y)
{
    double s = 0.0;
    size_t v;

    for (v = 0; v < NV; v++)
    {
        double px = y ? y[v] : (double)(v % NX);
        double py = y ? y[NV + v] : (double)((v / NX) % NY);
        double pz = y ? y[2 * NV + v] : (double)(v / ((size_t)NX * NY));
        double dg = trilinear(gm, px, py, pz) - tpl[v];
        double dw = trilinear(wm, px, py, pz) - tpl[NV + v];
        s += dg * dg + dw * dw;
    }
    return s / NV;
}

static int
run(const float *gm, const float *wm, float **tpl, double regstr, int threads,
    float *y, float *jac, CAT_VolShootInfo *info)
{
    const float *cls[2] = {gm, wm};
    const float *t[CAT_SHOOT_NTEMPLATES];
    CAT_VolShootOpts o;
    double M[16];
    int dims[3] = {NX, NY, NZ}, l;

    for (l = 0; l < CAT_SHOOT_NTEMPLATES; l++) t[l] = tpl[l];
    grid_matrix(M);
    CAT_VolShootOptionsInit(&o);
    o.regstr = regstr;
    o.n_threads = threads;
    return CAT_VolShoot(cls, 2, dims, M, NULL, t, dims, M, 1.5, NULL, &o, y, NULL, jac, info);
}

static int
all_finite(const float *v, size_t n)
{
    size_t i;
    for (i = 0; i < n; i++)
        if (!isfinite(v[i])) return 0;
    return 1;
}

/* The finest template registered to the series moves little.  It is not the
 * identity: the coarse levels match the sharp maps to smoothed templates. */
static void
test_identity(void)
{
    float **tpl = make_templates();
    float *y = (float *)malloc(sizeof(float) * 3 * NV), *jac = (float *)malloc(sizeof(float) * NV);
    CAT_VolShootInfo info;
    double mx = 0.0, mean = 0.0, jmin = HUGE_VAL, jmax = -HUGE_VAL, M[16];
    size_t v, nb = 0;
    int ret, l;

    ret = run(tpl[4], tpl[4] + NV, tpl, 0.5, 2, y, jac, &info);
    MU_ASSERT("identity: shooting succeeds", ret == 0);
    grid_matrix(M);
    for (v = 0; v < NV; v++)
    {
        int i = (int)(v % NX), j = (int)((v / NX) % NY), k = (int)(v / ((size_t)NX * NY));
        double d = sqrt((y[v] - i) * (y[v] - i) + (y[NV + v] - j) * (y[NV + v] - j) +
                        (y[2 * NV + v] - k) * (y[2 * NV + v] - k)), g, w;
        phantom(M[0] * i + M[3], M[5] * j + M[7], M[10] * k + M[11], 0.0, &g, &w);
        if (g + w < 0.5) continue;
        mean += d;
        nb++;
        if (d > mx) mx = d;
        if (jac[v] < jmin) jmin = jac[v];
        if (jac[v] > jmax) jmax = jac[v];
    }
    mean /= nb;
    if (mean >= 0.25 || mx >= 1.0 || jmin <= 0.8 || jmax >= 1.25)
        printf("  identity: |y - x| in the brain: mean %.3f, max %.3f voxels; det %.3f..%.3f\n",
               mean, mx, jmin, jmax);
    MU_ASSERT("identity: y stays near the identity in the brain", mean < 0.25 && mx < 1.0);
    MU_ASSERT("identity: Jacobians near 1 in the brain", jmin > 0.8 && jmax < 1.25);

    for (l = 0; l < CAT_SHOOT_NTEMPLATES; l++) free(tpl[l]);
    free(tpl);
    free(y);
    free(jac);
}

/* A known smooth deformation is undone: the mismatch drops and the field
 * follows the displacement. */
static void
test_recovers_deformation(void)
{
    float **tpl = make_templates();
    float *gm = (float *)malloc(sizeof(float) * NV), *wm = (float *)malloc(sizeof(float) * NV);
    float *y = (float *)malloc(sizeof(float) * 3 * NV), *jac = (float *)malloc(sizeof(float) * NV);
    float *y4 = (float *)malloc(sizeof(float) * 3 * NV);
    CAT_VolShootInfo info;
    double M[16], before, after, err = 0.0;
    size_t v, nb = 0;
    int ret, l;

    make_maps(gm, wm, 3.0, 0.0);
    before = mismatch(gm, wm, tpl[4], NULL);
    ret = run(gm, wm, tpl, 0.5, 1, y, jac, &info);
    MU_ASSERT("deformation: shooting succeeds", ret == 0);
    MU_ASSERT("deformation: finite field", all_finite(y, 3 * NV));
    after = mismatch(gm, wm, tpl[4], y);
    if (after >= 0.3 * before)
        printf("  deformation: mismatch %.5f -> %.5f\n", before, after);
    MU_ASSERT("deformation: the mismatch drops below 30%", after < 0.3 * before);
    MU_ASSERT("deformation: positive Jacobians", info.det_min > 0.0);

    /* inside the brain y - x follows the displacement */
    grid_matrix(M);
    for (v = 0; v < NV; v++)
    {
        int i = (int)(v % NX), j = (int)((v / NX) % NY), k = (int)(v / ((size_t)NX * NY));
        double x = M[0] * i + M[3], yy = M[5] * j + M[7], z = M[10] * k + M[11], d[3], g, w;
        phantom(x, yy, z, 0.0, &g, &w);
        if (g + w < 0.5) continue;
        displacement(x, yy, z, 3.0, d);
        err += fabs((y[v] - i) * VOX - d[0]) + fabs((y[NV + v] - j) * VOX - d[1]) +
               fabs((y[2 * NV + v] - k) * VOX - d[2]);
        nb++;
    }
    err /= 3.0 * nb;
    if (err >= 0.6)
        printf("  deformation: mean |error| %.3f mm\n", err);
    MU_ASSERT("deformation: displacement recovered within 0.6 mm on average", err < 0.6);

    /* the same result for every number of threads */
    ret = run(gm, wm, tpl, 0.5, 4, y4, jac, &info);
    MU_ASSERT("threads: shooting succeeds", ret == 0);
    MU_ASSERT("threads: 1 and 4 threads give identical fields",
              memcmp(y, y4, sizeof(float) * 3 * NV) == 0);

    /* SPM's default schedule converges as well */
    ret = run(gm, wm, tpl, 4, 2, y4, jac, &info);
    MU_ASSERT("regstr 4: shooting succeeds", ret == 0);
    after = mismatch(gm, wm, tpl[4], y4);
    MU_ASSERT("regstr 4: the mismatch drops below 30%", after < 0.3 * before);

    for (l = 0; l < CAT_SHOOT_NTEMPLATES; l++) free(tpl[l]);
    free(tpl);
    free(gm);
    free(wm);
    free(y);
    free(y4);
    free(jac);
}

/* Very large ventricles get a class of their own and the run stays finite. */
static void
test_hydro(void)
{
    float **tpl = make_templates();
    float *gm = (float *)malloc(sizeof(float) * NV), *wm = (float *)malloc(sizeof(float) * NV);
    float *y = (float *)malloc(sizeof(float) * 3 * NV), *jac = (float *)malloc(sizeof(float) * NV);
    CAT_VolShootInfo info;
    int ret, l;

    make_maps(gm, wm, 0.0, 11.0);
    ret = run(gm, wm, tpl, 0.5, 2, y, jac, &info);
    MU_ASSERT("hydro: shooting succeeds", ret == 0);
    if (info.hydro_levels == 0)
        printf("  hydro: the ventricle class was not used\n");
    MU_ASSERT("hydro: the ventricle class is used", info.hydro_levels > 0);
    MU_ASSERT("hydro: finite field", all_finite(y, 3 * NV) && all_finite(jac, NV));

    for (l = 0; l < CAT_SHOOT_NTEMPLATES; l++) free(tpl[l]);
    free(tpl);
    free(gm);
    free(wm);
    free(y);
    free(jac);
}

/* Invalid input is refused, not crashed on. */
static void
test_invalid(void)
{
    CAT_VolShootOpts o;
    float y[3];
    int dims[3] = {NX, NY, NZ};
    double M[16];

    grid_matrix(M);
    CAT_VolShootOptionsInit(&o);
    MU_ASSERT("defaults: regstr 0.5", o.regstr == 0.5);
    MU_ASSERT("defaults: 64 iterations", o.nits == 64);
    MU_ASSERT("invalid: no maps", CAT_VolShoot(NULL, 2, dims, M, NULL, NULL, dims, M, 1.5,
                                               NULL, &o, y, NULL, NULL, NULL) == 1);
    {
        float a = 0.0f;
        const float *cls[2] = {&a, &a}, *t[CAT_SHOOT_NTEMPLATES] = {&a, &a, &a, &a, &a};
        int nits;
        for (nits = 1; nits < 3; nits++)
        {
            o.nits = nits;
            MU_ASSERT("invalid: fewer than three iterations",
                      CAT_VolShoot(cls, 2, dims, M, NULL, t, dims, M, 1.5, NULL, &o, y, NULL,
                                   NULL, NULL) == 1);
        }
        o.nits = 64;
    }
    o.regstr = 2.0;
    {
        float a = 0.0f;
        const float *cls[2] = {&a, &a}, *t[CAT_SHOOT_NTEMPLATES] = {&a, &a, &a, &a, &a};
        MU_ASSERT("invalid: regstr 2 is not supported",
                  CAT_VolShoot(cls, 2, dims, M, NULL, t, dims, M, 1.5, NULL, &o, y, NULL,
                               NULL, NULL) == 1);
    }
}

static void
test_invalid_iterations(void)
{
    float **tpl = make_templates();
    float *y = (float *)malloc(sizeof(float) * 3 * NV);
    const float *cls[2] = {tpl[4], tpl[4] + NV};
    const float *t[CAT_SHOOT_NTEMPLATES];
    CAT_VolShootOpts o;
    CAT_VolShootInfo info;
    double M[16];
    int nits, l;

    for (l = 0; l < CAT_SHOOT_NTEMPLATES; l++) t[l] = tpl[l];
    grid_matrix(M);
    CAT_VolShootOptionsInit(&o);
    for (nits = 1; nits < 3; nits++)
    {
        o.nits = nits;
        MU_ASSERT("shooting rejects fewer than three iterations",
                  CAT_VolShoot(cls, 2, (int[]){NX, NY, NZ}, M, NULL, t, (int[]){NX, NY, NZ},
                               M, 1.5, NULL, &o, y, NULL, NULL, &info) == 1);
    }

    free(y);
    for (l = 0; l < CAT_SHOOT_NTEMPLATES; l++) free(tpl[l]);
    free(tpl);
}

int main(void)
{
    MU_RUN_TEST(test_invalid);
    MU_RUN_TEST(test_invalid_iterations);
    MU_RUN_TEST(test_identity);
    MU_RUN_TEST(test_recovers_deformation);
    MU_RUN_TEST(test_hydro);
    printf("%d tests run, %d failed\n", tests_run, tests_failed);
    return tests_failed ? 1 : 0;
}
