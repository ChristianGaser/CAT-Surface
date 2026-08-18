/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h>

#include "CAT_SulcusRepair.h"
#include "CAT_Sheetness.h"
#include "CAT_Vol.h"
#include "CAT_Math.h"

/**
 * \brief Fill a CAT_SulcusRepairOpts with defaults tuned for 0.5 mm data.
 *
 * \param opts (out) option block to initialize; NULL is ignored
 */
void CAT_SulcusRepairOptionsInit(CAT_SulcusRepairOpts *opts)
{
    if (!opts)
        return;

    opts->sheet_sigma_min = 0.3;
    opts->sheet_sigma_max = 1.0;
    opts->sheet_n_scales = 3;
    opts->sheet_strength = 1.0;

    opts->csf_min_dist = 1.5;
    opts->csf_min_wmdist = 0.75;
    opts->csf_thresh = 0.1;
    opts->csf_strength = 0.8;

    opts->wm_thresh = 0.1;
    opts->wm_strength = 0.8;
    opts->wm_min_int = 2.1;
    opts->wm_max_gap = 3;

    opts->band_min_dist = 1.5;
    opts->band_window = 4;
    opts->band_strength = 0.7;

    opts->verbose = 0;
}

/* median of a value list, using the library helper on a double copy */
static double median_of(double *buf, int n)
{
    if (n <= 0)
        return 0.0;
    return get_median_double(buf, n, 0);
}

/**
 * \brief Rescale an intensity image onto the 1..3 axis of a PVE label map.
 *
 * The class levels are medians rather than means so that a few misclassified
 * voxels -- vessels inside the WM sample, dura inside the CSF sample -- cannot
 * drag a level away. Voxels are subsampled by two in each direction, which is
 * ample for a robust level estimate and keeps the temporary buffers small.
 *
 * \param t1        (in)  intensity image, arbitrary units
 * \param label     (in)  PVE label image in [0..3]
 * \param t1n       (out) float[nvox] rescaled image on the 1..3 axis
 * \param dims      (in)  {nx, ny, nz}
 * \return 0 on success, -1 on invalid arguments, -3 if a class is unpopulated
 */
int CAT_VolNormalizeToLabel(const float *t1, const float *label, float *t1n,
                            int dims[3])
{
    const int nx = dims ? dims[0] : 0;
    const int ny = dims ? dims[1] : 0;
    const int nz = dims ? dims[2] : 0;
    const int xy = nx * ny;
    const int nvox = xy * nz;
    double *bc = NULL, *bg = NULL, *bw = NULL;
    double icsf, igm, iwm;
    int nc = 0, ng = 0, nw = 0;
    int x, y, z, i;
    int cap;

    if (!t1 || !label || !t1n || !dims || nvox <= 0)
        return -1;

    cap = (nx / 2 + 1) * (ny / 2 + 1) * (nz / 2 + 1);
    bc = (double *)malloc(sizeof(double) * (size_t)cap);
    bg = (double *)malloc(sizeof(double) * (size_t)cap);
    bw = (double *)malloc(sizeof(double) * (size_t)cap);
    if (!bc || !bg || !bw)
    {
        free(bc);
        free(bg);
        free(bw);
        return -2;
    }

    for (z = 0; z < nz; z += 2)
        for (y = 0; y < ny; y += 2)
            for (x = 0; x < nx; x += 2)
            {
                const int idx = x + y * nx + z * xy;
                const float l = label[idx];
                if (l > 2.9f && nw < cap)
                    bw[nw++] = (double)t1[idx];
                else if (l > 1.9f && l < 2.1f && ng < cap)
                    bg[ng++] = (double)t1[idx];
                else if (l > 0.9f && l < 1.1f && nc < cap)
                    bc[nc++] = (double)t1[idx];
            }

    if (nw < 10 || ng < 10 || nc < 10)
    {
        free(bc);
        free(bg);
        free(bw);
        return -3;
    }

    icsf = median_of(bc, nc);
    igm = median_of(bg, ng);
    iwm = median_of(bw, nw);

    free(bc);
    free(bg);
    free(bw);

    if (!(igm > icsf) || !(iwm > igm))
        return -3;

    for (i = 0; i < nvox; i++)
    {
        double v = (double)t1[i];
        double out;

        if (v <= igm)
            out = 1.0 + (v - icsf) / (igm - icsf);
        else
            out = 2.0 + (v - igm) / (iwm - igm);

        if (out < 0.0)
            out = 0.0;
        if (out > 3.0)
            out = 3.0;
        t1n[i] = (float)out;
    }

    return 0;
}

/**
 * \brief Recover sulcal CSF that the classifier missed (opens glued sulci).
 *
 * The gate is the conjunction of three independent pieces of evidence, which
 * is what keeps the operation from behaving like yet another filter:
 *
 *   - the label map reports no CSF within csf_min_dist (so this really is a
 *     place where the classifier committed to solid GM),
 *   - the intensity image has a dark-sheet response above csf_thresh (so
 *     there really is a thin dark structure there, not a blob or noise),
 *   - the voxel is at least csf_min_wmdist away from WM (ACE's T > 1 guard,
 *     so thin cortex is never eaten from the inside).
 *
 * Where all three hold the label is blended towards min(label, intensity), so
 * the correction can only lower a label and never below what the intensity
 * itself supports.
 *
 * \param t1        (in)     bias-corrected intensity image
 * \param label     (in/out) PVE label image in [0..3], corrected in place
 * \param sheetness (out)    float[nvox] dark-sheet response, or NULL to skip
 * \param dims      (in)     {nx, ny, nz}
 * \param voxelsize (in)     voxel spacing in mm
 * \param opts      (in)     parameters; NULL selects the defaults
 * \return 0 on success, negative on error
 */
int CAT_VolRecoverSulcalCSF(const float *t1, float *label, float *sheetness,
                            int dims[3], double voxelsize[3],
                            const CAT_SulcusRepairOpts *opts)
{
    CAT_SulcusRepairOpts defaults;
    CAT_SheetnessOpts sopts;
    const int nvox = dims ? dims[0] * dims[1] * dims[2] : 0;
    float *t1n = NULL, *S = NULL, *dcsf = NULL, *dwm = NULL;
    unsigned char *roi = NULL;
    int own_S = 0;
    int i, rc = 0;
    long n_changed = 0;

    if (!t1 || !label || !dims || !voxelsize || nvox <= 0)
        return -1;

    if (!opts)
    {
        CAT_SulcusRepairOptionsInit(&defaults);
        opts = &defaults;
    }

    t1n = (float *)malloc(sizeof(float) * (size_t)nvox);
    dcsf = (float *)malloc(sizeof(float) * (size_t)nvox);
    dwm = (float *)malloc(sizeof(float) * (size_t)nvox);
    roi = (unsigned char *)malloc(sizeof(unsigned char) * (size_t)nvox);
    if (sheetness)
        S = sheetness;
    else
    {
        S = (float *)malloc(sizeof(float) * (size_t)nvox);
        own_S = 1;
    }

    if (!t1n || !dcsf || !dwm || !roi || !S)
    {
        rc = -2;
        goto cleanup;
    }

    rc = CAT_VolNormalizeToLabel(t1, label, t1n, dims);
    if (rc != 0)
        goto cleanup;

    /* dark sheets are only interesting inside tissue */
    for (i = 0; i < nvox; i++)
        roi[i] = (label[i] > CSF) ? 1 : 0;

    CAT_SheetnessOptionsInit(&sopts);
    sopts.sigma_min = opts->sheet_sigma_min;
    sopts.sigma_max = opts->sheet_sigma_max;
    sopts.n_scales = opts->sheet_n_scales;
    sopts.gain = opts->sheet_strength;
    sopts.polarity = -1; /* dark sheet: a valley, l3 > 0 */
    sopts.verbose = opts->verbose;

    if (opts->verbose)
        fprintf(stderr, "Sulcal CSF recovery: dark-sheet filter on intensity.\n");

    rc = CAT_VolSheetness(t1n, S, NULL, roi, dims, voxelsize, &sopts);
    if (rc != 0)
        goto cleanup;

    /* distance to the CSF the classifier did find, and to WM, both in mm */
    for (i = 0; i < nvox; i++)
        dcsf[i] = (label[i] < CGM) ? 1.0f : 0.0f;
    euclidean_distance(dcsf, NULL, dims, voxelsize, 0);

    for (i = 0; i < nvox; i++)
        dwm[i] = (label[i] > GWM) ? 1.0f : 0.0f;
    euclidean_distance(dwm, NULL, dims, voxelsize, 0);

    for (i = 0; i < nvox; i++)
    {
        double w, target, s;

        if (!(label[i] > CGM && label[i] < GWM))
            continue;
        if (dcsf[i] <= (float)opts->csf_min_dist)
            continue;
        if (dwm[i] <= (float)opts->csf_min_wmdist)
            continue;

        s = (double)S[i];
        if (s <= opts->csf_thresh)
            continue;

        /* Ramp the strength in above the threshold so there is no hard edge.
           The ramp reaches full strength at twice the threshold, mirroring the
           oriented median, where a thin sheet is protected from twice its
           cutoff.  Anchoring it at a sheetness of 1 instead -- as this did
           originally -- assumes a response that saturates, and a real one does
           not: on a 0.5 mm MPRAGE the dark-sheet map reaches 0.56, so a voxel
           at the threshold-clearing 0.35 was blended by 0.06 rather than by
           csf_strength, which made the whole correction invisible. */
        if (opts->csf_thresh > 0.0)
            w = opts->csf_strength * (s - opts->csf_thresh) / opts->csf_thresh;
        else
            w = opts->csf_strength;
        if (w < 0.0)
            w = 0.0;
        if (w > opts->csf_strength)
            w = opts->csf_strength;

        /* one-sided: only ever towards the intensity-implied value */
        target = (double)t1n[i];
        if (target >= (double)label[i])
            continue;

        label[i] = (float)((1.0 - w) * (double)label[i] + w * target);
        if (label[i] < (float)CSF)
            label[i] = (float)CSF;
        n_changed++;
    }

    if (opts->verbose)
        fprintf(stderr, "  lowered %ld voxels towards CSF.\n", n_changed);

cleanup:
    free(t1n);
    free(dcsf);
    free(dwm);
    free(roi);
    if (own_S)
        free(S);
    return rc;
}

/**
 * \brief Strengthen thin WM blades the classifier under-labelled.
 *
 * The structures this protects are the fine white-matter fingers reaching into
 * the gyral crowns. They are one to two voxels across at their far end, so
 * partial volume pulls their intensity down towards GM, and a classifier that
 * resolves the trunk of the blade correctly still tends to lose its tip. What
 * is lost is not a hole inside the blade but its last millimetre, and that is
 * enough to corrupt the WM distance map and with it the thickness and the
 * central surface along the whole gyrus.
 *
 * This step therefore asks whether a voxel *continues* a blade, not whether it
 * sits in a gap inside one. The candidate set is "labelled GM, bright-sheet-like,
 * and brighter than the label admits", and a candidate is accepted when it is
 * geodesically connected to existing WM *through the candidate set itself*
 * (downcut_float with everything else marked as an inert barrier) and lies no
 * further than wm_max_gap voxels from that WM. Growing through the candidate set
 * is what makes the repair follow the blade instead of dilating isotropically,
 * and it is also what supplies the connectivity evidence: a bright speck sitting
 * alone in GM is never reached, while the tip of a real blade is reached from the
 * trunk behind it.
 *
 * An earlier version required WM on two *opposite* sides within wm_max_gap, the
 * dual of ACE's shock test. That condition is unsatisfiable at the end of a
 * blade -- there is WM behind the tip and grey matter in front of it -- so the
 * step could only ever repair interior gaps and never the crowns, which is where
 * blades actually break. On a 0.5 mm MPRAGE it rejected 86.8% of otherwise
 * eligible voxels.
 *
 * What keeps this from firing inside a sulcus is not that test but the polarity
 * guard of the sheetness filter: a sulcus is a *dark* sheet, a valley, and the
 * bright-sheet filter used here (polarity +1) does not respond to it at all.
 * The intensity floor wm_min_int and the one-sided blend do the rest -- a label
 * is only ever raised, and never past what the intensity itself supports.
 *
 * \param t1        (in)     bias-corrected intensity image
 * \param label     (in/out) PVE label image in [0..3], corrected in place
 * \param sheetness (out)    float[nvox] bright-sheet response, or NULL to skip
 * \param dims      (in)     {nx, ny, nz}
 * \param voxelsize (in)     voxel spacing in mm
 * \param opts      (in)     parameters; NULL selects the defaults
 * \return 0 on success, negative on error
 */
int CAT_VolStrengthenWmBlades(const float *t1, float *label, float *sheetness,
                              int dims[3], double voxelsize[3],
                              const CAT_SulcusRepairOpts *opts)
{
    CAT_SulcusRepairOpts defaults;
    CAT_SheetnessOpts sopts;
    const int nx = dims ? dims[0] : 0;
    const int ny = dims ? dims[1] : 0;
    const int nz = dims ? dims[2] : 0;
    const int nvox = nx * ny * nz;
    float *t1n = NULL, *S = NULL, *grow = NULL, *dwm = NULL;
    unsigned char *roi = NULL, *cand = NULL;
    double unit[3] = {1.0, 1.0, 1.0};
    int own_S = 0;
    int i, rc = 0;
    long n_cand = 0, n_changed = 0;
    double dd[2];

    if (!t1 || !label || !dims || !voxelsize || nvox <= 0)
        return -1;
    if (nx < 3 || ny < 3 || nz < 3)
        return -1;

    if (!opts)
    {
        CAT_SulcusRepairOptionsInit(&defaults);
        opts = &defaults;
    }

    t1n = (float *)malloc(sizeof(float) * (size_t)nvox);
    grow = (float *)malloc(sizeof(float) * (size_t)nvox);
    dwm = (float *)malloc(sizeof(float) * (size_t)nvox);
    roi = (unsigned char *)malloc(sizeof(unsigned char) * (size_t)nvox);
    cand = (unsigned char *)malloc(sizeof(unsigned char) * (size_t)nvox);
    if (sheetness)
        S = sheetness;
    else
    {
        S = (float *)malloc(sizeof(float) * (size_t)nvox);
        own_S = 1;
    }

    if (!t1n || !grow || !dwm || !roi || !cand || !S)
    {
        rc = -2;
        goto cleanup;
    }

    rc = CAT_VolNormalizeToLabel(t1, label, t1n, dims);
    if (rc != 0)
        goto cleanup;

    for (i = 0; i < nvox; i++)
        roi[i] = (label[i] > CGM) ? 1 : 0;

    CAT_SheetnessOptionsInit(&sopts);
    sopts.sigma_min = opts->sheet_sigma_min;
    sopts.sigma_max = opts->sheet_sigma_max;
    sopts.n_scales = opts->sheet_n_scales;
    sopts.gain = opts->sheet_strength;
    sopts.polarity = 1; /* bright sheet: a ridge, l3 < 0 */
    sopts.verbose = opts->verbose;

    if (opts->verbose)
        fprintf(stderr, "WM blade strengthening: bright-sheet filter on intensity.\n");

    rc = CAT_VolSheetness(t1n, S, NULL, roi, dims, voxelsize, &sopts);
    if (rc != 0)
        goto cleanup;

    /* Distance to existing WM, in voxels rather than mm: wm_max_gap is a
       neighbourhood reach, so it should not change meaning with the sampling. */
    for (i = 0; i < nvox; i++)
        dwm[i] = (label[i] > GWM) ? 1.0f : 0.0f;
    euclidean_distance(dwm, NULL, dims, unit, 0);

    /* candidate = GM-labelled, sheet-like, brighter than the label admits, and
       close enough to WM that it can plausibly be part of the same blade */
    for (i = 0; i < nvox; i++)
    {
        cand[i] = 0;

        if (!(label[i] > CGM && label[i] <= GWM))
            continue;
        if ((double)S[i] <= opts->wm_thresh)
            continue;
        if ((double)t1n[i] < opts->wm_min_int)
            continue;
        /* one-sided from the outset: nothing to do where the intensity does
           not already ask for more WM than the label has */
        if ((double)t1n[i] <= (double)label[i])
            continue;
        if ((double)dwm[i] > (double)opts->wm_max_gap)
            continue;

        cand[i] = 1;
        n_cand++;
    }

    if (opts->verbose)
        fprintf(stderr, "  %ld blade candidates.\n", n_cand);

    /* geodesic growth from WM, confined to the candidate set:
       WM = seed (1), candidates = free (0), everything else = barrier */
    for (i = 0; i < nvox; i++)
    {
        if (label[i] > (float)GWM)
            grow[i] = 1.0f;
        else if (cand[i])
            grow[i] = 0.0f;
        else
            grow[i] = -FLT_MAX;
    }

    dd[0] = 1.0; /* distance weight */
    dd[1] = 0.0; /* no intensity cost: the candidate set is the constraint */
    downcut_float(grow, t1n, NULL, dims, 0.5, voxelsize, dd);

    for (i = 0; i < nvox; i++)
    {
        double w, target;

        if (!cand[i] || grow[i] <= 0.0f)
            continue;

        /* same ramp as the CSF step: full strength at twice the threshold */
        if (opts->wm_thresh > 0.0)
            w = opts->wm_strength * ((double)S[i] - opts->wm_thresh) /
                opts->wm_thresh;
        else
            w = opts->wm_strength;
        if (w < 0.0)
            w = 0.0;
        if (w > opts->wm_strength)
            w = opts->wm_strength;

        /* one-sided: only ever towards the intensity-implied value */
        target = (double)t1n[i];
        if (target <= (double)label[i])
            continue;

        label[i] = (float)((1.0 - w) * (double)label[i] + w * target);
        if (label[i] > (float)WM)
            label[i] = (float)WM;
        n_changed++;
    }

    if (opts->verbose)
        fprintf(stderr, "  raised %ld voxels towards WM.\n", n_changed);

cleanup:
    free(t1n);
    free(grow);
    free(dwm);
    free(roi);
    free(cand);
    if (own_S)
        free(S);
    return rc;
}

/**
 * \brief Re-estimate partial volume in a narrow band from local class means.
 *
 * The band is the set of voxels where the label map claims no CSF within
 * band_min_dist. Inside it, a voxel is only refitted when the intensity has a
 * genuine local minimum along the sheet normal -- checking across the normal
 * rather than in all directions is what distinguishes a sulcal dip from a
 * local darkening caused by noise or by residual bias.
 *
 * The two bracketing class levels are re-estimated inside a local window, so
 * the refit is immune to whatever residual bias the global classifier left
 * behind; when the window contains too little CSF to estimate a level -- which
 * is precisely the failing case -- the global level is used instead. The
 * result is applied one-sided (it can only lower a label), so this step can
 * open a sulcus but never close one.
 *
 * \param t1        (in)     bias-corrected intensity image
 * \param label     (in/out) PVE label image in [0..3], corrected in place
 * \param dims      (in)     {nx, ny, nz}
 * \param voxelsize (in)     voxel spacing in mm
 * \param opts      (in)     parameters; NULL selects the defaults
 * \return 0 on success, negative on error
 */
int CAT_VolRefinePveNarrowBand(const float *t1, float *label,
                               int dims[3], double voxelsize[3],
                               const CAT_SulcusRepairOpts *opts)
{
    CAT_SulcusRepairOpts defaults;
    CAT_SheetnessOpts sopts;
    const int nx = dims ? dims[0] : 0;
    const int ny = dims ? dims[1] : 0;
    const int nz = dims ? dims[2] : 0;
    const int xy = nx * ny;
    const int nvox = xy * nz;
    float *S = NULL, *nrm = NULL, *dcsf = NULL, *t1n = NULL, *out = NULL;
    unsigned char *roi = NULL;
    double *wbuf = NULL;
    double gcsf, ggm;
    int x, y, z, i, rc = 0;
    int win, cap;
    long n_changed = 0;

    if (!t1 || !label || !dims || !voxelsize || nvox <= 0)
        return -1;
    if (nx < 3 || ny < 3 || nz < 3)
        return -1;

    if (!opts)
    {
        CAT_SulcusRepairOptionsInit(&defaults);
        opts = &defaults;
    }

    win = opts->band_window > 0 ? opts->band_window : 4;
    cap = (2 * win + 1) * (2 * win + 1) * (2 * win + 1);

    S = (float *)malloc(sizeof(float) * (size_t)nvox);
    nrm = (float *)malloc(sizeof(float) * 3 * (size_t)nvox);
    dcsf = (float *)malloc(sizeof(float) * (size_t)nvox);
    t1n = (float *)malloc(sizeof(float) * (size_t)nvox);
    out = (float *)malloc(sizeof(float) * (size_t)nvox);
    roi = (unsigned char *)malloc(sizeof(unsigned char) * (size_t)nvox);
    wbuf = (double *)malloc(sizeof(double) * (size_t)cap);

    if (!S || !nrm || !dcsf || !t1n || !out || !roi || !wbuf)
    {
        rc = -2;
        goto cleanup;
    }

    rc = CAT_VolNormalizeToLabel(t1, label, t1n, dims);
    if (rc != 0)
        goto cleanup;

    /* global class levels on the raw intensity, used as the fallback */
    {
        int nc = 0, ng = 0;
        double *bc = (double *)malloc(sizeof(double) * (size_t)(nvox / 8 + 1));
        double *bg = (double *)malloc(sizeof(double) * (size_t)(nvox / 8 + 1));
        if (!bc || !bg)
        {
            free(bc);
            free(bg);
            rc = -2;
            goto cleanup;
        }
        for (z = 0; z < nz; z += 2)
            for (y = 0; y < ny; y += 2)
                for (x = 0; x < nx; x += 2)
                {
                    const int idx = x + y * nx + z * xy;
                    if (label[idx] > 0.9f && label[idx] < 1.1f)
                        bc[nc++] = (double)t1[idx];
                    else if (label[idx] > 1.9f && label[idx] < 2.1f)
                        bg[ng++] = (double)t1[idx];
                }
        if (nc < 10 || ng < 10)
        {
            free(bc);
            free(bg);
            rc = -3;
            goto cleanup;
        }
        gcsf = median_of(bc, nc);
        ggm = median_of(bg, ng);
        free(bc);
        free(bg);
    }

    for (i = 0; i < nvox; i++)
        roi[i] = (label[i] > CSF) ? 1 : 0;

    CAT_SheetnessOptionsInit(&sopts);
    sopts.sigma_min = opts->sheet_sigma_min;
    sopts.sigma_max = opts->sheet_sigma_max;
    sopts.n_scales = opts->sheet_n_scales;
    sopts.gain = opts->sheet_strength;
    sopts.polarity = -1;
    sopts.verbose = opts->verbose;

    if (opts->verbose)
        fprintf(stderr, "Narrow-band PVE refit: locating dips across the sheet normal.\n");

    rc = CAT_VolSheetness(t1n, S, nrm, roi, dims, voxelsize, &sopts);
    if (rc != 0)
        goto cleanup;

    for (i = 0; i < nvox; i++)
        dcsf[i] = (label[i] < CGM) ? 1.0f : 0.0f;
    euclidean_distance(dcsf, NULL, dims, voxelsize, 0);

    memcpy(out, label, sizeof(float) * (size_t)nvox);

    for (z = win; z < nz - win; z++)
        for (y = win; y < ny - win; y++)
            for (x = win; x < nx - win; x++)
            {
                const int idx = x + y * nx + z * xy;
                double n0, n1, n2, lcsf, lgm, frac, refit, w;
                int px, py, pz, mx, my, mz;
                int wx, wy, wz, nc = 0, ng = 0;

                if (!(label[idx] > CGM && label[idx] < GWM))
                    continue;
                if (dcsf[idx] <= (float)opts->band_min_dist)
                    continue;

                n0 = (double)nrm[3 * idx + 0];
                n1 = (double)nrm[3 * idx + 1];
                n2 = (double)nrm[3 * idx + 2];
                if (n0 == 0.0 && n1 == 0.0 && n2 == 0.0)
                    continue;

                /* nearest voxel one step along +/- the sheet normal */
                px = x + (int)floor(n0 + 0.5);
                py = y + (int)floor(n1 + 0.5);
                pz = z + (int)floor(n2 + 0.5);
                mx = x - (int)floor(n0 + 0.5);
                my = y - (int)floor(n1 + 0.5);
                mz = z - (int)floor(n2 + 0.5);

                if (px == x && py == y && pz == z)
                    continue;

                /* a genuine dip: darker than both sides across the sheet */
                if (!(t1[idx] < t1[px + py * nx + pz * xy] &&
                      t1[idx] < t1[mx + my * nx + mz * xy]))
                    continue;

                /* local class levels; fall back to global where a class is
                   absent, which is exactly the situation we are repairing */
                for (wz = z - win; wz <= z + win; wz++)
                    for (wy = y - win; wy <= y + win; wy++)
                        for (wx = x - win; wx <= x + win; wx++)
                        {
                            const int w_idx = wx + wy * nx + wz * xy;
                            if (label[w_idx] > 0.9f && label[w_idx] < 1.1f)
                                wbuf[nc++] = (double)t1[w_idx];
                        }
                lcsf = (nc >= 10) ? median_of(wbuf, nc) : gcsf;

                for (wz = z - win; wz <= z + win; wz++)
                    for (wy = y - win; wy <= y + win; wy++)
                        for (wx = x - win; wx <= x + win; wx++)
                        {
                            const int w_idx = wx + wy * nx + wz * xy;
                            if (label[w_idx] > 1.9f && label[w_idx] < 2.1f)
                                wbuf[ng++] = (double)t1[w_idx];
                        }
                lgm = (ng >= 10) ? median_of(wbuf, ng) : ggm;

                if (!(lgm > lcsf))
                    continue;

                frac = ((double)t1[idx] - lcsf) / (lgm - lcsf);
                if (frac < 0.0)
                    frac = 0.0;
                if (frac > 1.0)
                    frac = 1.0;
                refit = 1.0 + frac;

                if (refit >= (double)label[idx])
                    continue; /* one-sided: never re-glue */

                w = opts->band_strength;
                if (w < 0.0)
                    w = 0.0;
                if (w > 1.0)
                    w = 1.0;

                out[idx] = (float)((1.0 - w) * (double)label[idx] + w * refit);
                n_changed++;
            }

    memcpy(label, out, sizeof(float) * (size_t)nvox);

    if (opts->verbose)
        fprintf(stderr, "  refitted %ld voxels.\n", n_changed);

cleanup:
    free(S);
    free(nrm);
    free(dcsf);
    free(t1n);
    free(out);
    free(roi);
    free(wbuf);
    return rc;
}
