#include "minunit.h"
#include <stdio.h>
#include "CAT_Sheetness.h"
#include "CAT_SulcusRepair.h"
#include "CAT_Vol.h"
#include "CAT_Math.h"

#include <stdlib.h>
#include <string.h>
#include <math.h>

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
/* the gain scales the response without breaking the zero invariant    */
/* ------------------------------------------------------------------ */

static void test_sheetness_gain(void)
{
    float *vol = (float *)malloc(sizeof(float) * N * N * N);
    float *S1 = (float *)malloc(sizeof(float) * N * N * N);
    float *S2 = (float *)malloc(sizeof(float) * N * N * N);
    CAT_SheetnessOpts opts;
    int x, y, z, i;
    int n_zero = 0, n_above_half = 0;

    MU_ASSERT("alloc", vol && S1 && S2);

    /* one dark plane at x = 8, as in the sheet-vs-blob phantom */
    for (z = 0; z < N; z++)
        for (y = 0; y < N; y++)
            for (x = 0; x < N; x++)
                vol[IDX(x, y, z)] = 2.0f;
    for (z = 0; z < N; z++)
        for (y = 0; y < N; y++)
            vol[IDX(8, y, z)] = 1.0f;

    CAT_SheetnessOptionsInit(&opts);
    MU_ASSERT("gain defaults to 1", opts.gain == 1.0);

    opts.polarity = -1;
    opts.sigma_min = 0.5;
    opts.sigma_max = 1.0;
    opts.n_scales = 2;

    /* a strongly attenuated response, so that doubling it stays below the
       clamp and the linearity can be checked exactly */
    opts.gain = 0.25;
    MU_ASSERT("sheetness runs",
              CAT_VolSheetness(vol, S1, NULL, NULL, dims3, vx1, &opts) == 0);
    opts.gain = 0.5;
    MU_ASSERT("sheetness runs",
              CAT_VolSheetness(vol, S2, NULL, NULL, dims3, vx1, &opts) == 0);

    for (i = 0; i < N * N * N; i++)
    {
        /* linear in the gain wherever the clamp does not bite */
        if (S2[i] < 1.0f)
            MU_ASSERT("response is linear in the gain",
                      fabs((double)S2[i] - 2.0 * (double)S1[i]) < 1e-5);
        /* and a zero response stays zero at any gain, which is what keeps
           every oriented operator identical to its isotropic counterpart
           away from thin structures */
        if (S1[i] == 0.0f)
        {
            MU_ASSERT("zero response stays zero", S2[i] == 0.0f);
            n_zero++;
        }
    }
    MU_ASSERT("the phantom has voxels with no response", n_zero > 0);

    /* the point of a gain above 1: lift a response that is too weak to
       clear the hard 0.5 threshold of the oriented median */
    opts.gain = 0.25;
    MU_ASSERT("sheetness runs",
              CAT_VolSheetness(vol, S1, NULL, NULL, dims3, vx1, &opts) == 0);
    for (i = 0; i < N * N * N; i++)
        MU_ASSERT("attenuated response cannot reach the median threshold",
                  S1[i] <= 0.5f);

    opts.gain = 4.0;
    MU_ASSERT("sheetness runs",
              CAT_VolSheetness(vol, S2, NULL, NULL, dims3, vx1, &opts) == 0);
    for (i = 0; i < N * N * N; i++)
    {
        MU_ASSERT("gain never pushes the response out of [0,1]",
                  S2[i] >= 0.0f && S2[i] <= 1.0f);
        if (S2[i] > 0.5f)
            n_above_half++;
    }
    MU_ASSERT("a gain above 1 lifts the sheet over the threshold",
              n_above_half > 0);

    free(vol);
    free(S1);
    free(S2);
}

/* ------------------------------------------------------------------ */
/* the admission cutoff drops the three neighbour classes in turn      */
/* ------------------------------------------------------------------ */

static void test_oriented_median_cutoff(void)
{
    const int nvox = N * N * N;
    float *vol = (float *)malloc(sizeof(float) * nvox);
    float *out = (float *)malloc(sizeof(float) * nvox);
    float *S = (float *)malloc(sizeof(float) * nvox);
    float *nrm = (float *)malloc(sizeof(float) * 3 * nvox);
    int x, y, z, i;

    MU_ASSERT("alloc", vol && out && S && nrm);

    /* a one-voxel-thick sheet at x = 12 with a known constant sheetness and a
       normal along x, so the admitted set follows the arithmetic exactly */
    for (z = 0; z < N; z++)
        for (y = 0; y < N; y++)
            for (x = 0; x < N; x++)
                vol[IDX(x, y, z)] = 2.0f;
    for (z = 0; z < N; z++)
        for (y = 0; y < N; y++)
            vol[IDX(12, y, z)] = 1.0f;

    for (i = 0; i < nvox; i++)
    {
        S[i] = 0.4f;
        nrm[3 * i + 0] = 1.0f;
        nrm[3 * i + 1] = 0.0f;
        nrm[3 * i + 2] = 0.0f;
    }

    /* s = 0.4, so the 6 face neighbours drop out from cutoff 0.4 downwards,
       the 12 edge neighbours from 0.2 and the 8 corners from 0.1333.
       A one-voxel sheet needs the edge neighbours gone, i.e. cutoff <= s/2. */

    /* above s: nothing is rejected, so this is the isotropic median */
    memcpy(out, vol, sizeof(float) * nvox);
    MU_ASSERT("oriented median runs",
              CAT_VolOrientedMedian(out, S, nrm, NULL, dims3, 0.5, 1) == 0);
    MU_ASSERT("a cutoff above the sheetness erases the sheet",
              out[IDX(12, 12, 12)] > 1.9f);

    /* between s and s/2: only the faces are gone, still not enough */
    memcpy(out, vol, sizeof(float) * nvox);
    MU_ASSERT("oriented median runs",
              CAT_VolOrientedMedian(out, S, nrm, NULL, dims3, 0.3, 1) == 0);
    MU_ASSERT("dropping the face neighbours alone does not keep the sheet",
              out[IDX(12, 12, 12)] > 1.9f);

    /* at s/2 the edge neighbours go and the 9 in-plane offsets become half
       the admitted set, which is where the sheet starts to survive */
    memcpy(out, vol, sizeof(float) * nvox);
    MU_ASSERT("oriented median runs",
              CAT_VolOrientedMedian(out, S, nrm, NULL, dims3, 0.2, 1) == 0);
    MU_ASSERT("the sheet survives from cutoff = s/2 downwards",
              out[IDX(12, 12, 12)] < 1.1f);

    /* The default is calibrated against the normalized response, where p99.9
       is 1 by construction and the p99 level of the reference scan sits at
       0.20/0.33 = 0.61.  Preservation must begin at or below that, which is
       exactly what 2 * 0.30 = 0.60 encodes. */
    for (i = 0; i < nvox; i++)
        S[i] = 0.61f;
    memcpy(out, vol, sizeof(float) * nvox);
    MU_ASSERT("oriented median runs",
              CAT_VolOrientedMedian(out, S, nrm, NULL, dims3, 0.0, 1) == 0);
    MU_ASSERT("the default cutoff keeps a sheet at the normalized p99 level",
              out[IDX(12, 12, 12)] < 1.1f);
    /* The default has to sit inside the anchored scale -- a cutoff at or above
       0.5 would demand the top 0.1% of the response before a thin sheet is
       preserved at all, which is the saturation failure the anchor removed. */
    MU_ASSERT("the default is a fraction of the anchor",
              CAT_ORIENTED_MEDIAN_CUTOFF > 0.0 &&
              CAT_ORIENTED_MEDIAN_CUTOFF < 0.5);

    /* and it must still leave a clearly sub-threshold response isotropic:
       preservation begins at 2 * cutoff, so anything below that is untouched */
    for (i = 0; i < nvox; i++)
        S[i] = (float)(CAT_ORIENTED_MEDIAN_CUTOFF * 0.9);
    memcpy(out, vol, sizeof(float) * nvox);
    MU_ASSERT("oriented median runs",
              CAT_VolOrientedMedian(out, S, nrm, NULL, dims3, 0.0, 1) == 0);
    MU_ASSERT("a response below the cutoff is not protected by the default",
              out[IDX(12, 12, 12)] > 1.9f);

    /* and the zero-sheetness invariant survives every cutoff */
    for (i = 0; i < nvox; i++)
        S[i] = 0.0f;
    {
        float *iso = (float *)malloc(sizeof(float) * nvox);
        MU_ASSERT("alloc", iso != NULL);
        /* the reference is the same routine with no field at all, so that the
           comparison isolates the cutoff and not the boundary handling */
        memcpy(iso, vol, sizeof(float) * nvox);
        MU_ASSERT("oriented median runs",
                  CAT_VolOrientedMedian(iso, NULL, NULL, NULL, dims3, 0.0, 1) == 0);
        for (i = 0; i < 4; i++)
        {
            const double cut[4] = {0.5, 0.25, 0.1, 0.01};
            memcpy(out, vol, sizeof(float) * nvox);
            MU_ASSERT("oriented median runs",
                      CAT_VolOrientedMedian(out, S, nrm, NULL, dims3, cut[i], 1) == 0);
            MU_ASSERT("zero sheetness is the isotropic median at any cutoff",
                      memcmp(out, iso, sizeof(float) * nvox) == 0);
        }
        free(iso);
    }

    free(vol);
    free(out);
    free(S);
    free(nrm);
}

/* ------------------------------------------------------------------ */
/* the repair ramp reaches full strength at twice the threshold        */
/* ------------------------------------------------------------------ */

static void test_repair_ramp_reaches_full_strength(void)
{
    CAT_SulcusRepairOpts opts;

    CAT_SulcusRepairOptionsInit(&opts);

    /* The blend weight is strength * (s - thresh) / thresh, clamped to
       [0, strength].  Anchoring the ramp at a sheetness of 1 instead assumes a
       response that saturates; a real one does not, which left the correction
       at a few percent of its nominal strength.  The relation mirrors the
       oriented median, where a thin sheet is protected from twice the cutoff. */
    /* The response is anchored to its own p99.9 = 1 (CAT_SheetnessOpts::
       normalize), so these thresholds are read as fractions of that anchor and
       mean the same thing on every dataset.  Full strength is reached at twice
       the threshold, which must land at or below the anchor -- otherwise the
       correction would need the top 0.1% of the response merely to reach its
       nominal strength, which is the failure the anchor was introduced to
       remove. */
    MU_ASSERT("the CSF threshold is a fraction of the anchor",
              opts.csf_thresh > 0.0 && opts.csf_thresh <= 0.5);
    MU_ASSERT("the WM threshold is a fraction of the anchor",
              opts.wm_thresh > 0.0 && opts.wm_thresh <= 0.5);

    MU_ASSERT("full CSF strength is reached at or below the anchor",
              2.0 * opts.csf_thresh <= 1.0);
    MU_ASSERT("full WM strength is reached at or below the anchor",
              2.0 * opts.wm_thresh <= 1.0);

    /* and it has to sit at the p99 level of the reference response, which the
       anchor puts at 0.20/0.33 = 0.61 -- the same point the oriented median's
       default cutoff encodes */
    MU_ASSERT("the CSF ramp tops out around the reference p99",
              2.0 * opts.csf_thresh > 0.5 && 2.0 * opts.csf_thresh <= 0.7);

    /* and the ramp must never exceed the nominal strength: the steeper slope
       makes (s - thresh) / thresh pass 1 long before the sheetness does */
    {
        const double s_hi = 1.0;
        double w = opts.csf_strength * (s_hi - opts.csf_thresh) / opts.csf_thresh;
        if (w > opts.csf_strength)
            w = opts.csf_strength;
        MU_ASSERT("the blend weight is clamped to csf_strength, not to 1",
                  w <= opts.csf_strength + 1e-12);
        MU_ASSERT("csf_strength still leaves the label short of a full swap",
                  opts.csf_strength < 1.0);
    }
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
              CAT_VolOrientedMedian(iso, NULL, NULL, NULL, dims3, 0.0, 1) == 0);
    MU_ASSERT("isotropic median erases the sheet", iso[IDX(12, 12, 12)] > 1.9f);

    /* oriented: neighbours across the sheet are excluded, the median stays dark */
    MU_ASSERT("oriented runs",
              CAT_VolOrientedMedian(ori, S, nrm, NULL, dims3, 0.0, 1) == 0);
    MU_ASSERT("oriented median keeps the sheet", ori[IDX(12, 12, 12)] < 1.1f);

    /* zero sheetness must reproduce the isotropic result exactly */
    for (i = 0; i < nvox; i++)
        S[i] = 0.0f;
    MU_ASSERT("oriented runs",
              CAT_VolOrientedMedian(flat, S, nrm, NULL, dims3, 0.0, 1) == 0);
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

static void test_strengthen_wm_blades(void)
{
    const int nvox = N * N * N;
    float *t1 = (float *)malloc(sizeof(float) * nvox);
    float *lab = (float *)malloc(sizeof(float) * nvox);
    CAT_SulcusRepairOpts opts;
    int x, y, z;
    float gap_after, tip_after, decoy_before, decoy_after;

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

    /* a WM blade one voxel thick at x = 12, running from z = 0 to z = 17 and
       interrupted at z = 12; beyond z = 17 it is the blade *tip*, still bright
       in the intensity but labelled GM -- the gyral-crown case */
    for (z = 0; z <= 17; z++)
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
    for (z = 18; z <= 19; z++)
        for (y = 0; y < N; y++)
        {
            lab[IDX(12, y, z)] = 2.0f;   /* the classifier lost the tip ... */
            t1[IDX(12, y, z)] = 136.0f;  /* ... though it is still nearly WM */
        }

    /* decoy: a lone WM-bright voxel far from any blade */
    lab[IDX(5, 5, 5)] = 2.0f;
    t1[IDX(5, 5, 5)] = 140.0f;
    decoy_before = lab[IDX(5, 5, 5)];

    CAT_SulcusRepairOptionsInit(&opts);
    opts.sheet_sigma_min = 0.5;
    opts.sheet_sigma_max = 1.0;
    opts.sheet_n_scales = 2;

    MU_ASSERT("strengthening runs",
              CAT_VolStrengthenWmBlades(t1, lab, NULL, dims3, vx1, &opts) == 0);

    gap_after = lab[IDX(12, 12, 12)];
    tip_after = lab[IDX(12, 12, 18)];
    decoy_after = lab[IDX(5, 5, 5)];

    MU_ASSERT("the blade gap was closed", gap_after > 2.4f);

    /* the point of the rewrite: a blade *tip* has WM only behind it, so the
       old two-sided test could never fire here, which is exactly where blades
       break in real data */
    MU_ASSERT("the lost blade tip was strengthened", tip_after > 2.4f);

    MU_ASSERT("the isolated bright voxel was not touched",
              fabs(decoy_after - decoy_before) < 1e-5);
    MU_ASSERT("nothing was raised past what the intensity supports",
              tip_after <= 3.0f && gap_after <= 3.0f);

    free(t1);
    free(lab);
}

/* ------------------------------------------------------------------ */
/* buried sulci in a PPM: open the valley, never cut the ridge         */
/* ------------------------------------------------------------------ */

static void test_open_ppm_sulci(void)
{
    const int nvox = N * N * N;
    float *ppm = (float *)malloc(sizeof(float) * nvox);
    float *before = (float *)malloc(sizeof(float) * nvox);
    CAT_PpmSulciOpts opts;
    int x, y, z;

    MU_ASSERT("alloc", ppm && before);

    /* Two structures, both one voxel thick, both above the isovalue.
     *
     *   x = 8   a buried sulcus: a valley floor at 0.60 between banks at 0.95.
     *           Crossing a sulcus the PPM dips, and here the dip stops short of
     *           0.5, so marching cubes would bridge the two banks.
     *   x = 16  a gyral blade: a ridge at 0.95 between flanks at 0.60.  The
     *           same amplitudes, the opposite sign of curvature -- if the
     *           correction keyed on "thin sheet" alone it would cut this open. */
    for (z = 0; z < N; z++)
        for (y = 0; y < N; y++)
            for (x = 0; x < N; x++)
                ppm[IDX(x, y, z)] = 0.95f;

    for (z = 0; z < N; z++)
        for (y = 0; y < N; y++)
        {
            ppm[IDX(8, y, z)] = 0.60f;   /* valley floor */

            ppm[IDX(15, y, z)] = 0.60f;  /* ridge flank  */
            ppm[IDX(16, y, z)] = 0.95f;  /* ridge crest  */
            ppm[IDX(17, y, z)] = 0.60f;  /* ridge flank  */
        }

    memcpy(before, ppm, sizeof(float) * nvox);

    CAT_PpmSulciOptionsInit(&opts);
    opts.sigma_min = 0.5;
    opts.sigma_max = 1.0;
    opts.n_scales = 2;

    MU_ASSERT("open runs",
              CAT_VolOpenPpmSulci(ppm, NULL, NULL, dims3, vx1, 0.5, &opts) == 0);

    MU_ASSERT("buried sulcus is pushed below the isovalue",
              ppm[IDX(8, 12, 12)] < 0.5f);
    MU_ASSERT("the gyral crest is left alone",
              ppm[IDX(16, 12, 12)] == before[IDX(16, 12, 12)]);
    MU_ASSERT("solid values well above the isovalue are left alone",
              ppm[IDX(4, 12, 12)] == before[IDX(4, 12, 12)]);

    /* strictly one-sided: nothing may ever be raised */
    {
        int i, raised = 0;
        for (i = 0; i < nvox; i++)
            if (ppm[i] > before[i])
                raised++;
        MU_ASSERT("no value is ever raised", raised == 0);
    }

    /* strength 0 must reproduce the input exactly */
    memcpy(ppm, before, sizeof(float) * nvox);
    opts.strength = 0.0;
    MU_ASSERT("open runs",
              CAT_VolOpenPpmSulci(ppm, NULL, NULL, dims3, vx1, 0.5, &opts) == 0);
    {
        int i, changed = 0;
        for (i = 0; i < nvox; i++)
            if (ppm[i] != before[i])
                changed++;
        MU_ASSERT("strength 0 is an exact no-op", changed == 0);
    }

    free(ppm);
    free(before);
}

/* ------------------------------------------------------------------ */
/* the sulcus guard: protect a nearby sulcus, but only where there is  */
/* one -- a guard that vetoes every blade is no better than disabling  */
/* the strengthening altogether                                        */
/* ------------------------------------------------------------------ */

static void build_blade_phantom(float *t1, float *lab, int with_sulcus)
{
    int x, y, z;

    for (z = 0; z < N; z++)
        for (y = 0; y < N; y++)
            for (x = 0; x < N; x++)
            {
                const int i = IDX(x, y, z);
                lab[i] = 2.0f;
                t1[i] = 90.0f;
            }

    /* some real CSF so the class levels can be estimated at all */
    for (z = 0; z < N; z++)
        for (y = 0; y < N; y++)
        {
            lab[IDX(0, y, z)] = 1.0f;
            t1[IDX(0, y, z)] = 40.0f;
        }

    /* a one-voxel WM blade at x = 12, interrupted at z = 12 */
    for (z = 0; z < N; z++)
        for (y = 0; y < N; y++)
        {
            lab[IDX(12, y, z)] = 3.0f;
            t1[IDX(12, y, z)] = 140.0f;
        }
    for (y = 0; y < N; y++)
    {
        lab[IDX(12, y, 12)] = 2.0f;
        t1[IDX(12, y, 12)] = 138.0f;
    }

    /* the occipital case: a sulcal floor one voxel from the blade tip */
    if (with_sulcus)
        for (y = 0; y < N; y++)
            t1[IDX(13, y, 12)] = 45.0f;
}

static void test_wm_sulcus_guard(void)
{
    const int nvox = N * N * N;
    float *t1 = (float *)malloc(sizeof(float) * nvox);
    float *lab = (float *)malloc(sizeof(float) * nvox);
    CAT_SulcusRepairOpts opts;
    float gap_guarded, gap_unguarded, gap_no_sulcus;

    MU_ASSERT("alloc", t1 && lab);

    CAT_SulcusRepairOptionsInit(&opts);
    opts.sheet_sigma_min = 0.5;
    opts.sheet_sigma_max = 1.0;
    opts.sheet_n_scales = 2;

    /* sulcus adjacent, guard off: the blade is raised and buries the sulcus */
    build_blade_phantom(t1, lab, 1);
    opts.wm_sulcus_guard = 0.0;
    MU_ASSERT("strengthening runs",
              CAT_VolStrengthenWmBlades(t1, lab, NULL, dims3, vx1, &opts) == 0);
    gap_unguarded = lab[IDX(12, 12, 12)];

    /* sulcus adjacent, guard on: the blade is left alone */
    build_blade_phantom(t1, lab, 1);
    opts.wm_sulcus_guard = 1.0;
    MU_ASSERT("strengthening runs",
              CAT_VolStrengthenWmBlades(t1, lab, NULL, dims3, vx1, &opts) == 0);
    gap_guarded = lab[IDX(12, 12, 12)];

    /* no sulcus anywhere near, guard on: the blade is still repaired */
    build_blade_phantom(t1, lab, 0);
    opts.wm_sulcus_guard = 1.0;
    MU_ASSERT("strengthening runs",
              CAT_VolStrengthenWmBlades(t1, lab, NULL, dims3, vx1, &opts) == 0);
    gap_no_sulcus = lab[IDX(12, 12, 12)];

    MU_ASSERT("without the guard the blade buries the sulcus",
              gap_unguarded > 2.4f);
    MU_ASSERT("the guard stops that", gap_guarded < gap_unguarded - 0.4f);
    MU_ASSERT("the guard does not block a blade with no sulcus nearby",
              gap_no_sulcus > 2.4f);

    free(t1);
    free(lab);
}

/* ------------------------------------------------------------------ */
/* percentile normalization makes the response scale data-independent  */
/* ------------------------------------------------------------------ */

static void test_sheetness_normalize(void)
{
    const int nvox = N * N * N;
    float *vol = (float *)malloc(sizeof(float) * nvox);
    float *raw = (float *)malloc(sizeof(float) * nvox);
    float *norm = (float *)malloc(sizeof(float) * nvox);
    float *scaled = (float *)malloc(sizeof(float) * nvox);
    CAT_SheetnessOpts opts;
    double pct[2] = {50.0, 99.9}, val[2];
    int x, y, z, i, order_kept = 1, zeros_kept = 1;

    MU_ASSERT("alloc", vol && raw && norm && scaled);

    /* a sheet plus a much stronger step, so the automatic noise scale is set
       by the step and the sheet response lands well below 1 */
    for (z = 0; z < N; z++)
        for (y = 0; y < N; y++)
            for (x = 0; x < N; x++)
                vol[IDX(x, y, z)] = (x < 4) ? 0.0f : 2.0f;
    for (z = 0; z < N; z++)
        for (y = 0; y < N; y++)
            vol[IDX(14, y, z)] = 1.7f;

    CAT_SheetnessOptionsInit(&opts);
    opts.polarity = -1;
    opts.sigma_min = 0.5;
    opts.sigma_max = 1.0;
    opts.n_scales = 2;

    opts.normalize = 0.0;
    MU_ASSERT("raw runs",
              CAT_VolSheetness(vol, raw, NULL, NULL, dims3, vx1, &opts) == 0);

    opts.normalize = 1.0;
    MU_ASSERT("normalized runs",
              CAT_VolSheetness(vol, norm, NULL, NULL, dims3, vx1, &opts) == 0);

    memcpy(scaled, norm, sizeof(float) * nvox);
    get_prctile(scaled, nvox, val, pct, 1, DT_FLOAT32);
    MU_ASSERT("p99.9 is moved onto the anchor", fabs(val[1] - 1.0) < 0.02);

    /* the scaling must be a single positive factor: same zero set, same order */
    for (i = 0; i < nvox; i++)
    {
        if ((raw[i] == 0.0f) != (norm[i] == 0.0f))
            zeros_kept = 0;
        if (i > 0 && raw[i] > raw[i - 1] && norm[i] < norm[i - 1])
            order_kept = 0;
    }
    MU_ASSERT("the zero set is unchanged", zeros_kept);
    MU_ASSERT("the ranking is unchanged", order_kept);

    /* a second dataset differing only by a global intensity factor has to give
       the same normalized response -- that is the whole point of the anchor */
    for (i = 0; i < nvox; i++)
        vol[i] *= 7.0f;
    MU_ASSERT("normalized runs",
              CAT_VolSheetness(vol, scaled, NULL, NULL, dims3, vx1, &opts) == 0);
    {
        double maxdiff = 0.0;
        for (i = 0; i < nvox; i++)
            if (fabs((double)scaled[i] - (double)norm[i]) > maxdiff)
                maxdiff = fabs((double)scaled[i] - (double)norm[i]);
        MU_ASSERT("a rescaled image gives the same normalized response",
                  maxdiff < 1e-3);
    }

    free(vol);
    free(raw);
    free(norm);
    free(scaled);
}

/* ------------------------------------------------------------------ */
/* the skeleton collapses a wide response onto its ridge               */
/* ------------------------------------------------------------------ */

static void test_sheetness_skeleton(void)
{
    const int nvox = N * N * N;
    float *vol = (float *)malloc(sizeof(float) * nvox);
    float *wide = (float *)malloc(sizeof(float) * nvox);
    float *thin = (float *)malloc(sizeof(float) * nvox);
    CAT_SheetnessOpts opts;
    int x, y, z, i, n_wide = 0, n_thin = 0, subset = 1;

    MU_ASSERT("alloc", vol && wide && thin);

    /* one dark sheet, and a scale large enough that the response spreads to
       either side of it -- the situation a big sigma_max creates */
    for (z = 0; z < N; z++)
        for (y = 0; y < N; y++)
            for (x = 0; x < N; x++)
                vol[IDX(x, y, z)] = 2.0f;
    for (z = 0; z < N; z++)
        for (y = 0; y < N; y++)
            vol[IDX(12, y, z)] = 1.0f;

    CAT_SheetnessOptionsInit(&opts);
    opts.polarity = -1;
    opts.sigma_min = 1.0;
    opts.sigma_max = 2.0;
    opts.n_scales = 2;
    opts.normalize = 0.0;   /* compare raw values, not two different anchors */

    opts.skeletonize = 0;
    MU_ASSERT("wide runs",
              CAT_VolSheetness(vol, wide, NULL, NULL, dims3, vx1, &opts) == 0);
    opts.skeletonize = 1;
    MU_ASSERT("skeleton runs",
              CAT_VolSheetness(vol, thin, NULL, NULL, dims3, vx1, &opts) == 0);

    /* the ridge keeps its value exactly: thinning must not rescale anything */
    MU_ASSERT("the value on the ridge is untouched",
              thin[IDX(12, 12, 12)] == wide[IDX(12, 12, 12)]);

    /* across the sheet only the ridge survives */
    MU_ASSERT("the shoulders responded before thinning",
              wide[IDX(11, 12, 12)] > 0.0f && wide[IDX(13, 12, 12)] > 0.0f);
    MU_ASSERT("and are gone after it",
              thin[IDX(11, 12, 12)] == 0.0f && thin[IDX(13, 12, 12)] == 0.0f);

    /* globally it can only remove, never add or raise */
    for (i = 0; i < nvox; i++)
    {
        if (wide[i] > 0.0f)
            n_wide++;
        if (thin[i] > 0.0f)
            n_thin++;
        if (thin[i] != 0.0f && thin[i] != wide[i])
            subset = 0;
    }
    MU_ASSERT("the skeleton is a subset of the response", subset);
    MU_ASSERT("and a strictly smaller one", n_thin < n_wide);
    MU_ASSERT("without emptying it", n_thin > 0);

    free(vol);
    free(wide);
    free(thin);
}

/* ------------------------------------------------------------------ */
/* the sulcal barrier: bound the CSF distance by the collision midline */
/* ------------------------------------------------------------------ */

static void build_bank_phantom(float *lab, int gap_is_csf)
{
    int x, y, z;

    for (z = 0; z < N; z++)
        for (y = 0; y < N; y++)
            for (x = 0; x < N; x++)
            {
                const int i = IDX(x, y, z);
                /* two WM cores with a grey shell, facing each other */
                if (x >= 4 && x < 8)
                    lab[i] = 3.0f;          /* WM bank A */
                else if (x >= 16 && x < 20)
                    lab[i] = 3.0f;          /* WM bank B */
                else if (x >= 8 && x < 16)
                    lab[i] = 2.0f;          /* GM: both banks and the gap */
                else
                    lab[i] = 1.0f;          /* CSF outside */
            }

    /* the sulcus itself: present in one phantom, lost in the other */
    if (gap_is_csf)
        for (z = 0; z < N; z++)
            for (y = 0; y < N; y++)
                for (x = 11; x < 13; x++)
                    lab[IDX(x, y, z)] = 1.0f;
}

static void test_sulcal_barrier(void)
{
    const int nvox = N * N * N;
    float *lab = (float *)malloc(sizeof(float) * nvox);
    float *dwm = (float *)malloc(sizeof(float) * nvox);
    float *medial = (float *)malloc(sizeof(float) * nvox);
    unsigned char *band = (unsigned char *)malloc(nvox);
    long n_ok, n_fused;
    int i;

    MU_ASSERT("alloc", lab && dwm && medial && band);

    /* A properly segmented sulcus splits the band in two, so no front ever
       meets another and there is nothing to find -- which is exactly why the
       barrier is free to be left on. */
    build_bank_phantom(lab, 1);
    for (i = 0; i < nvox; i++)
        dwm[i] = (lab[i] > GWM) ? 1.0f : 0.0f;
    euclidean_distance(dwm, NULL, dims3, NULL, 0);
    for (i = 0; i < nvox; i++)
        band[i] = (lab[i] > CGM && lab[i] < GWM) ? 1 : 0;
    n_ok = CAT_VolSulcalMedialSet(dwm, band, medial, dims3, 0.0, 1.0);

    /* Fuse the banks and the collision appears, with no intensity involved. */
    build_bank_phantom(lab, 0);
    for (i = 0; i < nvox; i++)
        dwm[i] = (lab[i] > GWM) ? 1.0f : 0.0f;
    euclidean_distance(dwm, NULL, dims3, NULL, 0);
    for (i = 0; i < nvox; i++)
        band[i] = (lab[i] > CGM && lab[i] < GWM) ? 1 : 0;
    n_fused = CAT_VolSulcalMedialSet(dwm, band, medial, dims3, 0.0, 1.0);

    MU_ASSERT("a fused sulcus produces a collision", n_fused > 0);
    MU_ASSERT("a segmented one produces far fewer", n_ok < n_fused / 4);

    /* and the collision sits on the midline, not against either bank */
    MU_ASSERT("the medial set is at the midline",
              medial[IDX(12, 12, 12)] > 0.0f);
    MU_ASSERT("and not at the WM boundary",
              medial[IDX(9, 12, 12)] == 0.0f &&
              medial[IDX(15, 12, 12)] == 0.0f);

    /* it must never appear outside the band it was given */
    for (i = 0; i < nvox; i++)
        if (medial[i] > 0.0f && !band[i])
        {
            MU_ASSERT("the medial set stays inside the mask", 0);
            break;
        }

    /* The dual: run the fronts inward from the pial boundary instead, and the
       same routine becomes a detector for lost white-matter blades.  Nothing
       about it changes except which distance map it is handed. */
    {
        float *dcsf = (float *)malloc(sizeof(float) * nvox);
        long g_ok, g_lost;

        MU_ASSERT("alloc", dcsf != NULL);

        /* a blade between two sulci, present then lost */
        build_bank_phantom(lab, 1);
        for (i = 0; i < nvox; i++)
            if (lab[i] > CGM && lab[i] < GWM && (i % N) >= 11 && (i % N) < 13)
                lab[i] = 3.0f; /* a blade where the sulcus was */

        for (i = 0; i < nvox; i++)
            dcsf[i] = (lab[i] < CGM) ? 1.0f : 0.0f;
        euclidean_distance(dcsf, NULL, dims3, NULL, 0);
        for (i = 0; i < nvox; i++)
            band[i] = (lab[i] > CGM && lab[i] < GWM) ? 1 : 0;
        g_ok = CAT_VolSulcalMedialSet(dcsf, band, medial, dims3, 0.0, 1.0);

        build_bank_phantom(lab, 0); /* blade lost: the gyrus is all grey */
        for (i = 0; i < nvox; i++)
            dcsf[i] = (lab[i] < CGM) ? 1.0f : 0.0f;
        euclidean_distance(dcsf, NULL, dims3, NULL, 0);
        for (i = 0; i < nvox; i++)
            band[i] = (lab[i] > CGM && lab[i] < GWM) ? 1 : 0;
        g_lost = CAT_VolSulcalMedialSet(dcsf, band, medial, dims3, 0.0, 1.0);

        MU_ASSERT("a lost blade produces a collision", g_lost > 0);
        MU_ASSERT("a surviving blade splits the band, so far fewer",
                  g_ok < g_lost / 4);

        free(dcsf);
    }

    free(lab);
    free(dwm);
    free(medial);
    free(band);
}

int main(void)
{
    MU_RUN_TEST(test_eigen_diagonal);
    MU_RUN_TEST(test_eigen_general);
    MU_RUN_TEST(test_sheetness_sheet_vs_blob);
    MU_RUN_TEST(test_sheetness_gain);
    MU_RUN_TEST(test_sheetness_normalize);
    MU_RUN_TEST(test_sheetness_skeleton);
    MU_RUN_TEST(test_sulcal_barrier);
    MU_RUN_TEST(test_oriented_median_preserves_sheet);
    MU_RUN_TEST(test_oriented_median_cutoff);
    MU_RUN_TEST(test_normalize_to_label);
    MU_RUN_TEST(test_recover_sulcal_csf);
    MU_RUN_TEST(test_repair_ramp_reaches_full_strength);
    MU_RUN_TEST(test_strengthen_wm_blades);
    MU_RUN_TEST(test_open_ppm_sulci);
    MU_RUN_TEST(test_wm_sulcus_guard);
    printf("%d tests run, %d failed\n", tests_run, tests_failed);
    return tests_failed ? 1 : 0;
}
