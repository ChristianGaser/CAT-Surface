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
    opts->sheet_normalize = CAT_SHEETNESS_NORMALIZE;
    opts->sheet_skeleton = 0;
    opts->sheet_strength = 1.0;

    opts->csf_min_dist = 1.5;
    opts->csf_min_wmdist = 0.75;
    opts->csf_thresh = 0.3;
    opts->csf_strength = 0.8;

    opts->wm_thresh = 0.3;
    opts->wm_strength = 0.8;
    opts->wm_min_int = 2.1;
    opts->wm_max_gap = 3;
    opts->wm_sulcus_guard = 1.0;

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
    sopts.normalize = opts->sheet_normalize;
    sopts.skeletonize = opts->sheet_skeleton;
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
    float *t1n = NULL, *S = NULL, *grow = NULL, *dwm = NULL, *S_dark = NULL;
    unsigned char *roi = NULL, *cand = NULL;
    double unit[3] = {1.0, 1.0, 1.0};
    int own_S = 0;
    int i, rc = 0;
    long n_cand = 0, n_changed = 0, n_vetoed = 0;
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
    if (opts->wm_sulcus_guard > 0.0 && opts->csf_thresh > 0.0)
        S_dark = (float *)malloc(sizeof(float) * (size_t)nvox);
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
    sopts.normalize = opts->sheet_normalize;
    sopts.skeletonize = opts->sheet_skeleton;
    sopts.polarity = 1; /* bright sheet: a ridge, l3 < 0 */
    sopts.verbose = opts->verbose;

    if (opts->verbose)
        fprintf(stderr, "WM blade strengthening: bright-sheet filter on intensity.\n");

    rc = CAT_VolSheetness(t1n, S, NULL, roi, dims, voxelsize, &sopts);
    if (rc != 0)
        goto cleanup;

    /* Second pass with the opposite polarity: where does a sulcus run?
     *
     * The bright ridge a blade tip sits on is real, so the polarity guard has
     * nothing to reject -- what must not happen is that filling the tip also
     * fills the dark sheet one voxel behind it.  Where cortex is thin, as in
     * the occipital lobe, the two are that close, which is why blade
     * strengthening buries sulci there and nowhere else in particular.
     *
     * Dilating the response by one voxel is what turns "this voxel is a
     * sulcus" into "this voxel touches a sulcus", which is the question that
     * actually matters when deciding whether to raise it. */
    if (S_dark)
    {
        CAT_SheetnessOpts dopts = sopts;
        float *dark_nrm = (float *)malloc(sizeof(float) * 3 * (size_t)nvox);
        int x, y, z;

        if (!dark_nrm)
        {
            rc = -2;
            goto cleanup;
        }

        dopts.polarity = -1;

        if (opts->verbose)
            fprintf(stderr, "  sulcus guard: dark-sheet filter on intensity.\n");

        rc = CAT_VolSheetness(t1n, S_dark, dark_nrm, roi, dims, voxelsize, &dopts);
        if (rc != 0)
        {
            free(dark_nrm);
            goto cleanup;
        }

        /* A polarity of -1 alone is not enough here.  The *shoulders* of a
           bright ridge curve upward into it, so their second derivative across
           the ridge is positive and the dark-sheet filter answers on them --
           which would make the guard veto every blade, including the ones with
           no sulcus anywhere near.  A response is therefore only kept where the
           intensity confirms it: the voxel has to be darker than pure grey
           matter and a genuine local minimum along its own normal, exactly the
           test CAT_VolRecoverSulcalCSF uses to find a sulcal floor. */
        for (z = 0; z < nz; z++)
            for (y = 0; y < ny; y++)
                for (x = 0; x < nx; x++)
                {
                    const int idx = x + y * nx + z * nx * ny;
                    double n0, n1, n2;
                    int px, py, pz, mx, my, mz;

                    if (S_dark[idx] <= 0.0f)
                        continue;

                    if ((double)t1n[idx] >= GM)
                    {
                        S_dark[idx] = 0.0f;
                        continue;
                    }

                    n0 = (double)dark_nrm[3 * idx + 0];
                    n1 = (double)dark_nrm[3 * idx + 1];
                    n2 = (double)dark_nrm[3 * idx + 2];

                    px = x + (int)floor(n0 + 0.5);
                    py = y + (int)floor(n1 + 0.5);
                    pz = z + (int)floor(n2 + 0.5);
                    mx = x - (int)floor(n0 + 0.5);
                    my = y - (int)floor(n1 + 0.5);
                    mz = z - (int)floor(n2 + 0.5);

                    if ((px == x && py == y && pz == z) ||
                        px < 0 || px >= nx || py < 0 || py >= ny ||
                        pz < 0 || pz >= nz ||
                        mx < 0 || mx >= nx || my < 0 || my >= ny ||
                        mz < 0 || mz >= nz)
                    {
                        S_dark[idx] = 0.0f;
                        continue;
                    }

                    if (!(t1n[idx] < t1n[px + py * nx + pz * nx * ny] &&
                          t1n[idx] < t1n[mx + my * nx + mz * nx * ny]))
                        S_dark[idx] = 0.0f;
                }

        free(dark_nrm);

        /* Dilating by one voxel is what turns "this voxel is a sulcal floor"
           into "this voxel touches one", which is the question that actually
           matters when deciding whether a blade tip may be raised. */
        grey_dilate(S_dark, dims, 1, DT_FLOAT32);
    }

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

        /* veto where a sulcus runs alongside: same ramp shape as everywhere
           else in this module, reaching a full veto at twice csf_thresh */
        if (S_dark)
        {
            double g = ((double)S_dark[i] - opts->csf_thresh) / opts->csf_thresh;
            if (g < 0.0)
                g = 0.0;
            if (g > 1.0)
                g = 1.0;
            w *= (1.0 - opts->wm_sulcus_guard * g);
            if (w <= 0.0)
            {
                n_vetoed++;
                continue;
            }
        }

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
    {
        fprintf(stderr, "  raised %ld voxels towards WM.\n", n_changed);
        if (S_dark)
            fprintf(stderr, "  sulcus guard vetoed %ld voxels.\n", n_vetoed);
    }

cleanup:
    free(t1n);
    free(grow);
    free(dwm);
    free(roi);
    free(cand);
    free(S_dark);
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
    sopts.normalize = opts->sheet_normalize;
    sopts.skeletonize = opts->sheet_skeleton;
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

/**
 * \brief Fill a CAT_PpmSulciOpts with defaults for a 0..1 PPM at 0.5 mm.
 *
 * The band of 0.25 keeps the operation to values in (0.5, 0.75] for the usual
 * isovalue: a buried sulcus sits just above the threshold, whereas solid white
 * matter is close to 1 and must never be touched.  The margin of 0.05 puts an
 * opened valley just far enough below the isovalue for marching cubes to
 * separate the banks without carving a wide gap.
 *
 * \param opts (out) option block to initialize; NULL is ignored
 */
void CAT_PpmSulciOptionsInit(CAT_PpmSulciOpts *opts)
{
    if (!opts)
        return;

    opts->sigma_min = 0.3;
    opts->sigma_max = 3.0;
    opts->n_scales = 3;
    opts->sheet_normalize = CAT_SHEETNESS_NORMALIZE;
    opts->sheet_skeleton = 1;
    opts->sheet_strength = 10.0;
    opts->thresh = 0.3;
    opts->band = 0.25;
    opts->margin = 0.05;
    opts->strength = 1.0;
    opts->offset = 0.0;
    opts->cutoff = 0.0;
    opts->verbose = 0;
}

/**
 * \brief Push buried sulcal valleys in a PPM below the isovalue.
 *
 * See the header for the three gating conditions.  The local-minimum test is
 * the important one: it is evaluated along the sheet normal rather than over
 * the whole neighbourhood, because that is what distinguishes a sulcal valley
 * floor from an ordinary dark voxel, and it is what makes the operation
 * structurally unable to fire on a gyral crown -- crossing a blade the PPM has
 * a maximum along the normal, never a minimum.
 *
 * \param ppm       (in/out) percentage position map in [0,1], corrected in place
 * \param sheetness (in)     precomputed dark-sheet response, or NULL to compute
 * \param normal    (in)     precomputed unit sheet normals, or NULL to compute
 * \param dims      (in)     {nx, ny, nz}
 * \param voxelsize (in)     voxel spacing in mm
 * \param isovalue  (in)     threshold the surface will be extracted at
 * \param opts      (in)     parameters; NULL selects the defaults
 * \return 0 on success, negative on error
 */
int CAT_VolOpenPpmSulci(float *ppm, const float *sheetness, const float *normal,
                        int dims[3], double voxelsize[3], double isovalue,
                        const CAT_PpmSulciOpts *opts)
{
    CAT_PpmSulciOpts defaults;
    CAT_SheetnessOpts sopts;
    const int nx = dims ? dims[0] : 0;
    const int ny = dims ? dims[1] : 0;
    const int nz = dims ? dims[2] : 0;
    const int xy = nx * ny;
    const int nvox = xy * nz;
    const float *S = sheetness;
    const float *nrm = normal;
    float *own_S = NULL, *own_nrm = NULL;
    int x, y, z, rc = 0;
    long n_changed = 0;

    if (!ppm || !dims || !voxelsize || nvox <= 0)
        return -1;
    if (nx < 3 || ny < 3 || nz < 3)
        return -1;

    if (!opts)
    {
        CAT_PpmSulciOptionsInit(&defaults);
        opts = &defaults;
    }

    if (!S || !nrm)
    {
        own_S = (float *)malloc(sizeof(float) * (size_t)nvox);
        own_nrm = (float *)malloc(sizeof(float) * 3 * (size_t)nvox);
        if (!own_S || !own_nrm)
        {
            free(own_S);
            free(own_nrm);
            return -2;
        }

        CAT_SheetnessOptionsInit(&sopts);
        sopts.sigma_min = opts->sigma_min;
        sopts.sigma_max = opts->sigma_max;
        sopts.n_scales = opts->n_scales;
        sopts.gain = opts->sheet_strength;
        sopts.normalize = opts->sheet_normalize;
        sopts.skeletonize = opts->sheet_skeleton;
        sopts.polarity = -1; /* a sulcus is a valley in the PPM */
        sopts.verbose = opts->verbose;

        if (opts->verbose)
            fprintf(stderr, "Open buried sulci: dark-sheet filter on the PPM.\n");

        rc = CAT_VolSheetness(ppm, own_S, own_nrm, NULL, dims, voxelsize, &sopts);
        if (rc != 0)
        {
            free(own_S);
            free(own_nrm);
            return rc;
        }
        S = own_S;
        nrm = own_nrm;
    }

    for (z = 1; z < nz - 1; z++)
        for (y = 1; y < ny - 1; y++)
            for (x = 1; x < nx - 1; x++)
            {
                const int idx = x + y * nx + z * xy;
                double s, w, target, n0, n1, n2;
                int px, py, pz, mx, my, mz;

                /* only marginal values: above the isovalue but not solid WM */
                if (!((double)ppm[idx] > isovalue &&
                      (double)ppm[idx] < isovalue + opts->band))
                    continue;

                s = (double)S[idx];
                if (s <= opts->thresh)
                    continue;

                n0 = (double)nrm[3 * idx + 0];
                n1 = (double)nrm[3 * idx + 1];
                n2 = (double)nrm[3 * idx + 2];
                if (n0 == 0.0 && n1 == 0.0 && n2 == 0.0)
                    continue;

                px = x + (int)floor(n0 + 0.5);
                py = y + (int)floor(n1 + 0.5);
                pz = z + (int)floor(n2 + 0.5);
                mx = x - (int)floor(n0 + 0.5);
                my = y - (int)floor(n1 + 0.5);
                mz = z - (int)floor(n2 + 0.5);

                if (px == x && py == y && pz == z)
                    continue;
                if (px < 0 || px >= nx || py < 0 || py >= ny || pz < 0 || pz >= nz ||
                    mx < 0 || mx >= nx || my < 0 || my >= ny || mz < 0 || mz >= nz)
                    continue;

                /* a valley floor, not a ridge: this is what a gyral crown can
                   never satisfy, so the test doubles as the gyrus guard */
                if (!(ppm[idx] < ppm[px + py * nx + pz * xy] &&
                      ppm[idx] < ppm[mx + my * nx + mz * xy]))
                    continue;

                w = opts->strength * (s - opts->thresh) / (1.0 - opts->thresh);
                if (w < 0.0)
                    w = 0.0;
                if (w > 1.0)
                    w = 1.0;

                target = isovalue - opts->margin;
                if (target < 0.0)
                    target = 0.0;

                /* one-sided: only ever lower a value */
                if (target >= (double)ppm[idx])
                    continue;

                ppm[idx] = (float)((1.0 - w) * (double)ppm[idx] + w * target);
                n_changed++;
            }

    if (opts->verbose)
        fprintf(stderr, "  opened %ld voxels below the isovalue.\n", n_changed);

    free(own_S);
    free(own_nrm);
    return 0;
}
