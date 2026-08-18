#include "minunit.h"
#include <stdio.h>
#include "CAT_Sheetness.h"
#include "CAT_SulcusRepair.h"
#include "CAT_Amap.h"
#include "CAT_Vol.h"

#include <stdlib.h>
#include <string.h>
#include <math.h>

/* ICM() is internal to CAT_Amap.c and deliberately not part of the public
   header; declared here so the no-op property of the anisotropic prior can be
   checked directly. */
extern void ICM(unsigned char *prob, unsigned char *label, int n_classes, int *dims,
                double beta, int iterations, double *voxelsize, int verbose,
                const double *class_weights, const CAT_MrfAnisoField *aniso);

#define N 24
#define IDX(x, y, z) ((x) + (y) * N + (z) * N * N)

static int dims3[3] = {N, N, N};
static double vx1[3] = {1.0, 1.0, 1.0};

/* ------------------------------------------------------------------ */
/* eigen decomposition                                                 */
/* ------------------------------------------------------------------ */

static void test_eigen_diagonal(void)
{
    /* diag(1, -5, 2): sorted by magnitude -> 1, 2, -5 */
    double a[6] = {1.0, 0.0, 0.0, -5.0, 0.0, 2.0};
    double ev[3], n[3];

    CAT_EigenSym3(a, ev, n);

    MU_ASSERT("|l1| smallest", fabs(ev[0] - 1.0) < 1e-9);
    MU_ASSERT("|l2| middle", fabs(ev[1] - 2.0) < 1e-9);
    MU_ASSERT("|l3| largest keeps its sign", fabs(ev[2] + 5.0) < 1e-9);
    MU_ASSERT("normal points along y", fabs(fabs(n[1]) - 1.0) < 1e-6);
}

static void test_eigen_general(void)
{
    /* symmetric matrix with known trace and determinant */
    double a[6] = {4.0, 1.0, 1.0, 4.0, 1.0, 4.0};
    double ev[3], n[3];
    double trace, det;

    CAT_EigenSym3(a, ev, n);

    trace = ev[0] + ev[1] + ev[2];
    det = ev[0] * ev[1] * ev[2];

    /* trace = 12, det = 4*(16-1) - 1*(4-1) + 1*(1-4) = 60 - 3 - 3 = 54 */
    MU_ASSERT("trace preserved", fabs(trace - 12.0) < 1e-6);
    MU_ASSERT("determinant preserved", fabs(det - 54.0) < 1e-5);
    MU_ASSERT("eigenvector is unit length",
              fabs(sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]) - 1.0) < 1e-6);
}

/* ------------------------------------------------------------------ */
/* sheetness separates a sheet from a blob                             */
/* ------------------------------------------------------------------ */

static void test_sheetness_sheet_vs_blob(void)
{
    float *vol = (float *)malloc(sizeof(float) * N * N * N);
    float *S = (float *)malloc(sizeof(float) * N * N * N);
    float *nrm = (float *)malloc(sizeof(float) * 3 * N * N * N);
    CAT_SheetnessOpts opts;
    int x, y, z;
    double s_sheet, s_blob, nx_abs;

    MU_ASSERT("alloc", vol && S && nrm);

    /* bright background, one dark plane at x = 8, one dark blob at (16,16,16) */
    for (z = 0; z < N; z++)
        for (y = 0; y < N; y++)
            for (x = 0; x < N; x++)
                vol[IDX(x, y, z)] = 2.0f;

    for (z = 0; z < N; z++)
        for (y = 0; y < N; y++)
            vol[IDX(8, y, z)] = 1.0f;

    for (z = 15; z <= 17; z++)
        for (y = 15; y <= 17; y++)
            for (x = 15; x <= 17; x++)
                vol[IDX(x, y, z)] = 1.0f;

    CAT_SheetnessOptionsInit(&opts);
    opts.polarity = -1; /* dark sheets */
    opts.sigma_min = 0.5;
    opts.sigma_max = 1.0;
    opts.n_scales = 2;

    MU_ASSERT("sheetness runs",
              CAT_VolSheetness(vol, S, nrm, NULL, dims3, vx1, &opts) == 0);

    s_sheet = S[IDX(8, 12, 12)];
    s_blob = S[IDX(16, 16, 16)];

    MU_ASSERT("thin dark sheet responds", s_sheet > 0.5);
    MU_ASSERT("dark blob is rejected", s_blob < 0.2);
    MU_ASSERT("sheet outranks blob", s_sheet > 3.0 * s_blob + 0.2);

    /* the sheet normal must point across the plane, i.e. along x */
    nx_abs = fabs((double)nrm[3 * IDX(8, 12, 12) + 0]);
    MU_ASSERT("normal is across the sheet", nx_abs > 0.95);

    /* the same structure is invisible to the opposite polarity */
    opts.polarity = 1;
    MU_ASSERT("sheetness runs",
              CAT_VolSheetness(vol, S, nrm, NULL, dims3, vx1, &opts) == 0);
    MU_ASSERT("bright-sheet filter ignores a dark sheet", S[IDX(8, 12, 12)] < 0.05);

    free(vol);
    free(S);
    free(nrm);
}

/* ------------------------------------------------------------------ */
/* the oriented median keeps what the isotropic one destroys           */
/* ------------------------------------------------------------------ */

static void test_oriented_median_preserves_sheet(void)
{
    const int nvox = N * N * N;
    float *iso = (float *)malloc(sizeof(float) * nvox);
    float *ori = (float *)malloc(sizeof(float) * nvox);
    float *flat = (float *)malloc(sizeof(float) * nvox);
    float *S = (float *)malloc(sizeof(float) * nvox);
    float *nrm = (float *)malloc(sizeof(float) * 3 * nvox);
    int x, y, z, i;

    MU_ASSERT("alloc", iso && ori && flat && S && nrm);

    for (z = 0; z < N; z++)
        for (y = 0; y < N; y++)
            for (x = 0; x < N; x++)
                iso[IDX(x, y, z)] = 2.0f;
    for (z = 0; z < N; z++)
        for (y = 0; y < N; y++)
            iso[IDX(12, y, z)] = 1.0f;

    memcpy(ori, iso, sizeof(float) * nvox);
    memcpy(flat, iso, sizeof(float) * nvox);

    /* an ideal sheet: sheetness 1 with a normal along x */
    for (i = 0; i < nvox; i++)
    {
        S[i] = 1.0f;
        nrm[3 * i + 0] = 1.0f;
        nrm[3 * i + 1] = 0.0f;
        nrm[3 * i + 2] = 0.0f;
    }

    /* isotropic: nine dark values among twenty-seven, the median is bright */
    MU_ASSERT("isotropic runs",
              CAT_VolOrientedMedian(iso, NULL, NULL, NULL, dims3, 1) == 0);
    MU_ASSERT("isotropic median erases the sheet", iso[IDX(12, 12, 12)] > 1.9f);

    /* oriented: neighbours across the sheet are excluded, the median stays dark */
    MU_ASSERT("oriented runs",
              CAT_VolOrientedMedian(ori, S, nrm, NULL, dims3, 1) == 0);
    MU_ASSERT("oriented median keeps the sheet", ori[IDX(12, 12, 12)] < 1.1f);

    /* zero sheetness must reproduce the isotropic result exactly */
    for (i = 0; i < nvox; i++)
        S[i] = 0.0f;
    MU_ASSERT("oriented runs",
              CAT_VolOrientedMedian(flat, S, nrm, NULL, dims3, 1) == 0);
    for (i = 0; i < nvox; i++)
        MU_ASSERT("zero sheetness is exactly the isotropic filter",
                  flat[i] == iso[i]);

    free(iso);
    free(ori);
    free(flat);
    free(S);
    free(nrm);
}

/* ------------------------------------------------------------------ */
/* intensity normalization onto the label axis                         */
/* ------------------------------------------------------------------ */

static void test_normalize_to_label(void)
{
    const int nvox = N * N * N;
    float *t1 = (float *)malloc(sizeof(float) * nvox);
    float *lab = (float *)malloc(sizeof(float) * nvox);
    float *out = (float *)malloc(sizeof(float) * nvox);
    int i;

    MU_ASSERT("alloc", t1 && lab && out);

    /* three equally sized slabs at arbitrary intensities 40 / 90 / 140 */
    for (i = 0; i < nvox; i++)
    {
        int z = i / (N * N);
        if (z < 8)
        {
            lab[i] = 1.0f;
            t1[i] = 40.0f;
        }
        else if (z < 16)
        {
            lab[i] = 2.0f;
            t1[i] = 90.0f;
        }
        else
        {
            lab[i] = 3.0f;
            t1[i] = 140.0f;
        }
    }

    MU_ASSERT("normalization runs",
              CAT_VolNormalizeToLabel(t1, lab, out, dims3) == 0);

    MU_ASSERT("CSF maps to 1", fabs(out[IDX(4, 4, 4)] - 1.0f) < 1e-4);
    MU_ASSERT("GM maps to 2", fabs(out[IDX(4, 4, 12)] - 2.0f) < 1e-4);
    MU_ASSERT("WM maps to 3", fabs(out[IDX(4, 4, 20)] - 3.0f) < 1e-4);

    free(t1);
    free(lab);
    free(out);
}

/* ------------------------------------------------------------------ */
/* glued sulcus: the label map lost the CSF, the intensity did not     */
/* ------------------------------------------------------------------ */

static void test_recover_sulcal_csf(void)
{
    const int nvox = N * N * N;
    float *t1 = (float *)malloc(sizeof(float) * nvox);
    float *lab = (float *)malloc(sizeof(float) * nvox);
    CAT_SulcusRepairOpts opts;
    int x, y, z;
    float before, after, elsewhere_before, elsewhere_after;

    MU_ASSERT("alloc", t1 && lab);

    /* WM slabs at x < 6 and x > 17, GM in between, and a dark CSF sheet at
       x = 12 that the label map does not know about (it says pure GM) */
    for (z = 0; z < N; z++)
        for (y = 0; y < N; y++)
            for (x = 0; x < N; x++)
            {
                const int i = IDX(x, y, z);
                if (x < 6 || x > 17)
                {
                    lab[i] = 3.0f;
                    t1[i] = 140.0f;
                }
                else
                {
                    lab[i] = 2.0f;
                    t1[i] = 90.0f;
                }
            }

    /* a few real CSF voxels so the class levels can be estimated at all */
    for (z = 0; z < N; z++)
        for (y = 0; y < N; y++)
        {
            lab[IDX(0, y, z)] = 1.0f;
            t1[IDX(0, y, z)] = 40.0f;
        }

    /* the sulcus: dark in the intensity, still labelled GM */
    for (z = 0; z < N; z++)
        for (y = 0; y < N; y++)
            t1[IDX(12, y, z)] = 45.0f;

    elsewhere_before = lab[IDX(9, 12, 12)];
    before = lab[IDX(12, 12, 12)];

    CAT_SulcusRepairOptionsInit(&opts);
    opts.sheet_sigma_min = 0.5;
    opts.sheet_sigma_max = 1.0;
    opts.sheet_n_scales = 2;

    MU_ASSERT("recovery runs",
              CAT_VolRecoverSulcalCSF(t1, lab, NULL, dims3, vx1, &opts) == 0);

    after = lab[IDX(12, 12, 12)];
    elsewhere_after = lab[IDX(9, 12, 12)];

    MU_ASSERT("the glued sulcus was opened", after < before - 0.4f);
    MU_ASSERT("it was not opened past what the intensity supports", after >= 1.0f);
    MU_ASSERT("ordinary GM was left alone",
              fabs(elsewhere_after - elsewhere_before) < 1e-5);

    free(t1);
    free(lab);
}

/* ------------------------------------------------------------------ */
/* broken gyrus: a gap in a WM blade, and a decoy with no second bank  */
/* ------------------------------------------------------------------ */

static void test_reconnect_gyri(void)
{
    const int nvox = N * N * N;
    float *t1 = (float *)malloc(sizeof(float) * nvox);
    float *lab = (float *)malloc(sizeof(float) * nvox);
    CAT_SulcusRepairOpts opts;
    int x, y, z;
    float gap_after, decoy_before, decoy_after;

    MU_ASSERT("alloc", t1 && lab);

    /* GM everywhere, with CSF and WM present so the levels can be estimated */
    for (z = 0; z < N; z++)
        for (y = 0; y < N; y++)
            for (x = 0; x < N; x++)
            {
                const int i = IDX(x, y, z);
                lab[i] = 2.0f;
                t1[i] = 90.0f;
            }
    for (z = 0; z < N; z++)
        for (y = 0; y < N; y++)
        {
            lab[IDX(0, y, z)] = 1.0f;
            t1[IDX(0, y, z)] = 40.0f;
        }

    /* a WM blade one voxel thick at x = 12, interrupted at z = 12 */
    for (z = 0; z < N; z++)
        for (y = 0; y < N; y++)
        {
            lab[IDX(12, y, z)] = 3.0f;
            t1[IDX(12, y, z)] = 140.0f;
        }
    for (y = 0; y < N; y++)
    {
        lab[IDX(12, y, 12)] = 2.0f; /* the missegmentation */
        t1[IDX(12, y, 12)] = 138.0f; /* but the intensity is still WM-bright */
    }

    /* decoy: a lone WM-bright voxel with no WM facing it on two sides */
    lab[IDX(5, 5, 5)] = 2.0f;
    t1[IDX(5, 5, 5)] = 140.0f;
    decoy_before = lab[IDX(5, 5, 5)];

    CAT_SulcusRepairOptionsInit(&opts);
    opts.sheet_sigma_min = 0.5;
    opts.sheet_sigma_max = 1.0;
    opts.sheet_n_scales = 2;

    MU_ASSERT("reconnection runs",
              CAT_VolReconnectGyri(t1, lab, NULL, dims3, vx1, &opts) == 0);

    gap_after = lab[IDX(12, 12, 12)];
    decoy_after = lab[IDX(5, 5, 5)];

    MU_ASSERT("the blade gap was bridged", gap_after > 2.4f);
    MU_ASSERT("the isolated bright voxel was not touched",
              fabs(decoy_after - decoy_before) < 1e-5);

    free(t1);
    free(lab);
}

/* ------------------------------------------------------------------ */
/* the anisotropic MRF is a no-op where nothing sheet-like was found   */
/* ------------------------------------------------------------------ */

static void test_aniso_mrf_is_noop_without_sheets(void)
{
    const int nvox = N * N * N;
    const int n_classes = 3;
    unsigned char *prob = (unsigned char *)malloc((size_t)nvox * n_classes);
    unsigned char *lab_ref = (unsigned char *)malloc(nvox);
    unsigned char *lab_beta = (unsigned char *)malloc(nvox);
    unsigned char *lab_potts = (unsigned char *)malloc(nvox);
    float *S = (float *)malloc(sizeof(float) * nvox);
    float *nrm = (float *)malloc(sizeof(float) * 3 * nvox);
    CAT_MrfAnisoField aniso;
    double cw[3] = {1.0, 1.0, 1.0};
    int i, c;

    MU_ASSERT("alloc", prob && lab_ref && lab_beta && lab_potts && S && nrm);

    /* a noisy but deterministic three-class field */
    for (i = 0; i < nvox; i++)
    {
        int z = i / (N * N);
        unsigned char base = (unsigned char)(z < 8 ? 1 : (z < 16 ? 2 : 3));
        if ((i % 37) == 0)
            base = (unsigned char)(base % 3 + 1);
        lab_ref[i] = base;
        for (c = 0; c < n_classes; c++)
            prob[i + c * nvox] = (unsigned char)((c + 1 == base) ? 200 : 28);
        S[i] = 0.0f;
        nrm[3 * i + 0] = 1.0f;
        nrm[3 * i + 1] = 0.0f;
        nrm[3 * i + 2] = 0.0f;
    }
    memcpy(lab_beta, lab_ref, nvox);
    memcpy(lab_potts, lab_ref, nvox);

    ICM(prob, lab_ref, n_classes, dims3, 0.3, 5, vx1, 0, cw, NULL);

    aniso.mode = MRF_ANISO_BETA;
    aniso.sheetness = S;
    aniso.normal = nrm;
    aniso.strength = 1.0;
    aniso.sigma = 0.5;
    ICM(prob, lab_beta, n_classes, dims3, 0.3, 5, vx1, 0, cw, &aniso);

    aniso.mode = MRF_ANISO_POTTS;
    ICM(prob, lab_potts, n_classes, dims3, 0.3, 5, vx1, 0, cw, &aniso);

    for (i = 0; i < nvox; i++)
    {
        MU_ASSERT("local-beta mode is a no-op at zero sheetness",
                  lab_beta[i] == lab_ref[i]);
        MU_ASSERT("anisotropic-Potts mode is a no-op at zero sheetness",
                  lab_potts[i] == lab_ref[i]);
    }

    free(prob);
    free(lab_ref);
    free(lab_beta);
    free(lab_potts);
    free(S);
    free(nrm);
}

/* ------------------------------------------------------------------ */
/* a sheet-aligned MRF stops the prior from eating a one-voxel sheet   */
/* ------------------------------------------------------------------ */

static void test_aniso_mrf_protects_thin_sheet(void)
{
    const int nvox = N * N * N;
    const int n_classes = 2;
    unsigned char *prob = (unsigned char *)malloc((size_t)nvox * n_classes);
    unsigned char *lab_iso = (unsigned char *)malloc(nvox);
    unsigned char *lab_ani = (unsigned char *)malloc(nvox);
    float *S = (float *)malloc(sizeof(float) * nvox);
    float *nrm = (float *)malloc(sizeof(float) * 3 * nvox);
    CAT_MrfAnisoField aniso;
    double cw[2] = {1.0, 1.0};
    int x, y, z, c;
    int surv_iso = 0, surv_ani = 0, total = 0;

    MU_ASSERT("alloc", prob && lab_iso && lab_ani && S && nrm);

    /* class 2 everywhere, with a one-voxel-thick class-1 sheet at x = 12.
       The data term only mildly prefers class 1 on the sheet, so an isotropic
       Potts prior -- where 18 of the 26 neighbours disagree -- overrules it,
       which is the shrinking bias in its simplest form. */
    for (z = 0; z < N; z++)
        for (y = 0; y < N; y++)
            for (x = 0; x < N; x++)
            {
                const int idx = IDX(x, y, z);
                const int on_sheet = (x == 12);

                lab_iso[idx] = (unsigned char)(on_sheet ? 1 : 2);
                for (c = 0; c < n_classes; c++)
                    prob[idx + c * nvox] = 40;
                prob[idx + (on_sheet ? 0 : 1) * nvox] = 90;

                S[idx] = on_sheet ? 1.0f : 0.0f;
                nrm[3 * idx + 0] = 1.0f;
                nrm[3 * idx + 1] = 0.0f;
                nrm[3 * idx + 2] = 0.0f;
            }
    memcpy(lab_ani, lab_iso, nvox);

    ICM(prob, lab_iso, n_classes, dims3, 0.6, 10, vx1, 0, cw, NULL);

    aniso.mode = MRF_ANISO_POTTS;
    aniso.sheetness = S;
    aniso.normal = nrm;
    aniso.strength = 1.0;
    aniso.sigma = 0.5;
    ICM(prob, lab_ani, n_classes, dims3, 0.6, 10, vx1, 0, cw, &aniso);

    for (z = 1; z < N - 1; z++)
        for (y = 1; y < N - 1; y++)
        {
            const int idx = IDX(12, y, z);
            total++;
            if (lab_iso[idx] == 1)
                surv_iso++;
            if (lab_ani[idx] == 1)
                surv_ani++;
        }

    MU_ASSERT("the isotropic prior eats most of the sheet",
              surv_iso < total / 2);
    MU_ASSERT("the anisotropic prior keeps it",
              surv_ani > (9 * total) / 10);

    free(prob);
    free(lab_iso);
    free(lab_ani);
    free(S);
    free(nrm);
}

int main(void)
{
    MU_RUN_TEST(test_eigen_diagonal);
    MU_RUN_TEST(test_eigen_general);
    MU_RUN_TEST(test_sheetness_sheet_vs_blob);
    MU_RUN_TEST(test_oriented_median_preserves_sheet);
    MU_RUN_TEST(test_normalize_to_label);
    MU_RUN_TEST(test_recover_sulcal_csf);
    MU_RUN_TEST(test_reconnect_gyri);
    MU_RUN_TEST(test_aniso_mrf_is_noop_without_sheets);
    MU_RUN_TEST(test_aniso_mrf_protects_thin_sheet);
    printf("%d tests run, %d failed\n", tests_run, tests_failed);
    return tests_failed ? 1 : 0;
}
