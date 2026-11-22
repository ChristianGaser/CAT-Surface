#include "minunit.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "CAT_UnifiedSeg.h"

/*
 * Unified segmentation on an analytic head phantom.
 *
 * The phantom is a head of concentric shells: a white matter core whose
 * surface is folded by a few sinusoids, grey matter, CSF, bone, soft tissue
 * and air.  The five templates are its class maps smoothed from coarse
 * (Template_0) to fine (Template_4).  The image shows the same head at
 * another position, zoom and rotation and with a smooth deformation, with one
 * intensity per class, a smooth bias field and noise.  Class maps and the
 * image are evaluated from the formula, so the truth is known at every voxel.
 */

#define TX 45
#define TY 53
#define TZ 45
#define TVOX 4.0

#define IX 50
#define IY 58
#define IZ 50
#define IVOX 3.6

#define NCLS CAT_USEG_NCLASSES

static const size_t NT = (size_t)TX * TY * TZ;
static const size_t NI = (size_t)IX * IY * IZ;

/* intensity of GM, WM, CSF, bone, soft tissue, air; soft tissue varies by
 * +-30 around its value, as fat and muscle do */
static const double MU[NCLS] = {65.0, 100.0, 30.0, 12.0, 150.0, 3.0};

/* voxels deep in white matter that can be made brighter: every third one
 * along each axis within HOT_RADIUS mm (template) of the centre of the head,
 * too few to widen the Gaussian of white matter by much; HOT_Z is how many SD
 * of white matter they are raised by */
#define HOT_RADIUS 12.0
#define HOT_Z 5.6

/** \brief Smooth step from 1 (inside, d < 0) to 0 (outside), width w mm. */
static double
inside(double d, double w)
{
    return 0.5 * (1.0 - tanh(d / w));
}

/** \brief Class probabilities at template position (x, y, z) in mm. */
static void
phantom(double x, double y, double z, double p[NCLS])
{
    double r = 67.0 * sqrt((x * x) / (67.0 * 67.0) + (y * y) / (80.0 * 80.0) + (z * z) / (67.0 * 67.0));
    double fold = 4.0 * sin(0.13 * x) * cos(0.11 * y) + 2.7 * sin(0.09 * z + 0.15 * y);
    double w = 1.9;
    double in0 = inside(r - (25.0 + fold), w); /* white matter */
    double in1 = inside(r - 38.0, w);          /* brain */
    double in2 = inside(r - 46.0, w);          /* inside the skull */
    double in3 = inside(r - 55.0, w);          /* skull */
    double in4 = inside(r - 65.0, w);          /* head */

    p[1] = in0;
    p[0] = (in1 - in0 > 0.0) ? in1 - in0 : 0.0;
    p[2] = (in2 - in1 > 0.0) ? in2 - in1 : 0.0;
    p[3] = (in3 - in2 > 0.0) ? in3 - in2 : 0.0;
    p[4] = (in4 - in3 > 0.0) ? in4 - in3 : 0.0;
    p[5] = 1.0 - in4;
}

/** \brief NIfTI-style matrix of a grid centred on the origin. */
static void
grid_matrix(int nx, int ny, int nz, double vox, double M[16])
{
    memset(M, 0, sizeof(double) * 16);
    M[0] = M[5] = M[10] = vox;
    M[3] = -vox * (nx - 1) / 2.0;
    M[7] = -vox * (ny - 1) / 2.0;
    M[11] = -vox * (nz - 1) / 2.0;
    M[15] = 1.0;
}

/**
 * \brief The world mapping from the image to the template: a zoom of 1.1, a
 *        rotation of 0.1 rad about z and a shift.
 */
static void
true_affine(double A[16])
{
    double c = cos(0.1), s = sin(0.1), zoom = 1.1;

    memset(A, 0, sizeof(double) * 16);
    A[0] = zoom * c;
    A[1] = -zoom * s;
    A[4] = zoom * s;
    A[5] = zoom * c;
    A[10] = zoom;
    A[3] = 6.0;
    A[7] = -4.0;
    A[11] = 5.0;
    A[15] = 1.0;
}

/** \brief Template position of the image position s: affine, then a smooth deformation. */
static void
to_template(const double s[3], double t[3])
{
    double A[16], a[3];
    int r;

    true_affine(A);
    for (r = 0; r < 3; r++) a[r] = A[r * 4] * s[0] + A[r * 4 + 1] * s[1] + A[r * 4 + 2] * s[2] + A[r * 4 + 3];
    t[0] = a[0] + 3.0 * sin(0.03 * a[1] + 0.3) * cos(0.026 * a[2]);
    t[1] = a[1] + 2.4 * sin(0.026 * a[0]) * cos(0.034 * a[2] + 0.2);
    t[2] = a[2] + 1.8 * cos(0.03 * a[0] + 0.1) * sin(0.022 * a[1]);
}

/** \brief Log of the bias field at image position s. */
static double
log_bias(const double s[3])
{
    return 0.25 * sin(s[0] / 45.0 + 0.3) + 0.2 * cos(s[1] / 50.0) - 0.15 * sin(s[2] / 40.0);
}

/** \brief Box smoothing with radius r voxels along each axis, in place. */
static void
box_smooth(float *v, int r)
{
    float *t = (float *)malloc(sizeof(float) * NT);
    int d[3] = {TX, TY, TZ}, ax, i, j, k, o;

    for (ax = 0; ax < 3; ax++)
    {
        memcpy(t, v, sizeof(float) * NT);
        for (k = 0; k < TZ; k++)
            for (j = 0; j < TY; j++)
                for (i = 0; i < TX; i++)
                {
                    int c[3] = {i, j, k};
                    double s = 0.0;
                    int n = 0;
                    for (o = -r; o <= r; o++)
                    {
                        int cc[3] = {c[0], c[1], c[2]};
                        cc[ax] += o;
                        if (cc[ax] < 0 || cc[ax] >= d[ax]) continue;
                        s += t[cc[0] + (size_t)TX * (cc[1] + (size_t)TY * cc[2])];
                        n++;
                    }
                    v[i + (size_t)TX * (j + (size_t)TY * k)] = (float)(s / n);
                }
    }
    free(t);
}

/** \brief The five templates of six classes, smoothed coarse to fine. */
static float **
make_templates(void)
{
    float **tpl = (float **)malloc(sizeof(float *) * CAT_SHOOT_NTEMPLATES);
    double M[16], p[NCLS];
    int l, i, j, k, c;

    grid_matrix(TX, TY, TZ, TVOX, M);
    for (l = 0; l < CAT_SHOOT_NTEMPLATES; l++)
    {
        tpl[l] = (float *)malloc(sizeof(float) * NCLS * NT);
        for (k = 0; k < TZ; k++)
            for (j = 0; j < TY; j++)
                for (i = 0; i < TX; i++)
                {
                    size_t v = i + (size_t)TX * (j + (size_t)TY * k);
                    phantom(M[0] * i + M[3], M[5] * j + M[7], M[10] * k + M[11], p);
                    /* no class falls below 1e-4, as in CAT's templates */
                    for (c = 0; c < NCLS; c++)
                        tpl[l][c * NT + v] = (float)((p[c] + 1e-4) / (1.0 + NCLS * 1e-4));
                }
        if (l < CAT_SHOOT_NTEMPLATES - 1)
            for (c = 0; c < NCLS; c++) box_smooth(tpl[l] + c * NT, (CAT_SHOOT_NTEMPLATES - l) / 2);
    }
    return tpl;
}

static void
free_templates(float **tpl)
{
    int l;
    for (l = 0; l < CAT_SHOOT_NTEMPLATES; l++) free(tpl[l]);
    free(tpl);
}

/** \brief A normal deviate from a linear congruential generator. */
static double
lcg_normal(unsigned long long *state)
{
    double u1, u2;
    *state = *state * 6364136223846793005ULL + 1442695040888963407ULL;
    u1 = ((double)((*state >> 11) & 0xFFFFFFFFFFFFFULL) + 0.5) / 4503599627370496.0;
    *state = *state * 6364136223846793005ULL + 1442695040888963407ULL;
    u2 = ((double)((*state >> 11) & 0xFFFFFFFFFFFFFULL) + 0.5) / 4503599627370496.0;
    return sqrt(-2.0 * log(u1)) * cos(6.283185307179586 * u2);
}

/** \brief Whether voxel (i, j, k) at template position t is one of the bright ones. */
static int
is_hot(int i, int j, int k, const double t[3])
{
    return i % 3 == 0 && j % 3 == 0 && k % 3 == 0 &&
           sqrt(t[0] * t[0] + t[1] * t[1] + t[2] * t[2]) < HOT_RADIUS;
}

static double hot_factor = 1.0; /* brightness of these voxels */

/**
 * \brief The image with noise of SD sd and the bias field, and the class of
 *        every voxel (truth).
 */
static void
make_image(float *img, unsigned char *truth, double sd, int with_bias)
{
    double M[16], p[NCLS], s[3], t[3];
    unsigned long long state = 12345ULL;
    int i, j, k, c;

    grid_matrix(IX, IY, IZ, IVOX, M);
    for (k = 0; k < IZ; k++)
        for (j = 0; j < IY; j++)
            for (i = 0; i < IX; i++)
            {
                size_t v = i + (size_t)IX * (j + (size_t)IY * k);
                double val = 0.0;
                int best = 0;
                s[0] = M[0] * i + M[3];
                s[1] = M[5] * j + M[7];
                s[2] = M[10] * k + M[11];
                to_template(s, t);
                phantom(t[0], t[1], t[2], p);
                for (c = 0; c < NCLS; c++)
                {
                    double mu = MU[c];
                    if (c == 4) mu += 30.0 * sin(0.07 * t[1] + 0.05 * t[2]);
                    val += p[c] * mu;
                    if (p[c] > p[best]) best = c;
                }
                if (is_hot(i, j, k, t)) val *= hot_factor;
                if (with_bias) val *= exp(log_bias(s));
                val += sd * lcg_normal(&state);
                img[v] = (float)((val > 0.2) ? val : 0.2);
                if (truth) truth[v] = (unsigned char)best;
            }
}

/** \brief Dice coefficient of class c between the largest posterior and the truth. */
static double
dice(const unsigned char *prob, const unsigned char *truth, int c)
{
    size_t v, a = 0, b = 0, ab = 0;
    int k;

    for (v = 0; v < NI; v++)
    {
        int best = 0;
        for (k = 1; k < NCLS; k++)
            if (prob[k * NI + v] > prob[best * NI + v]) best = k;
        a += (best == c);
        b += (truth[v] == c);
        ab += (best == c && truth[v] == c);
    }
    return (a + b) ? 2.0 * (double)ab / (double)(a + b) : 0.0;
}

/** \brief Options for the phantom: its voxels are 4 mm, one Gaussian per class. */
static void
phantom_options(CAT_UnifiedSegOpts *o)
{
    int k;

    CAT_UnifiedSegOptionsInit(o);
    for (k = 0; k < NCLS; k++) o->ngaus[k] = 1;
    o->shoot.lowres = 6.5;
    o->samp = 7.0;
}

static int
run(const float *img, float **tpl, const CAT_UnifiedSegOpts *o, unsigned char *prob,
    float *label, float *bias, CAT_UnifiedSegResult *res)
{
    const float *t[CAT_SHOOT_NTEMPLATES];
    const int dims[3] = {IX, IY, IZ}, tdims[3] = {TX, TY, TZ};
    double M[16], Mt[16];
    int l;

    for (l = 0; l < CAT_SHOOT_NTEMPLATES; l++) t[l] = tpl[l];
    grid_matrix(IX, IY, IZ, IVOX, M);
    grid_matrix(TX, TY, TZ, TVOX, Mt);
    return CAT_UnifiedSeg(img, dims, M, t, tdims, Mt, NULL, NULL, o, prob, label, NULL, NULL,
                          bias, NULL, NULL, NULL, res);
}

/** \brief The smaller variance of the Gaussians of GM and WM (one each). */
static double
tissue_variance(const CAT_UnifiedSegResult *res)
{
    double v = HUGE_VAL;
    int g;

    for (g = 0; g < res->n_gauss; g++)
        if (res->lkp[g] < 2 && res->vr[g] < v) v = res->vr[g];
    return v;
}

/* The classes, the affine matrix, the bias field and the label map of the
 * phantom are recovered, with both bases of the bias field. */
static void
test_phantom(void)
{
    float **tpl = make_templates();
    float *img = (float *)malloc(sizeof(float) * NI);
    float *label = (float *)malloc(sizeof(float) * NI);
    float *bias = (float *)malloc(sizeof(float) * NI);
    unsigned char *prob = (unsigned char *)malloc(NCLS * NI);
    unsigned char *truth = (unsigned char *)malloc(NI);
    CAT_UnifiedSegOpts o;
    CAT_UnifiedSegResult res;
    double A[16], M[16], s[3], t[3], d, sxy = 0, sx = 0, sy = 0, sxx = 0, syy = 0, n = 0, lw = 0, lg = 0;
    size_t v, nw = 0, ng = 0;
    int model, i, j, k, ret;

    make_image(img, truth, 2.0, 1);
    grid_matrix(IX, IY, IZ, IVOX, M);
    for (model = CAT_USEG_BIAS_SPLINE; model >= CAT_USEG_BIAS_DCT; model--)
    {
        phantom_options(&o);
        o.bias_model = model;
        ret = run(img, tpl, &o, prob, label, bias, &res);
        MU_ASSERT("segmentation of the phantom succeeds", ret == 0);
        if (ret) continue;
        MU_ASSERT("grey matter is recovered", dice(prob, truth, 0) > 0.95);
        MU_ASSERT("white matter is recovered", dice(prob, truth, 1) > 0.97);
        MU_ASSERT("CSF is recovered", dice(prob, truth, 2) > 0.90);
        MU_ASSERT("bone is recovered", dice(prob, truth, 3) > 0.75);
        MU_ASSERT("soft tissue is recovered", dice(prob, truth, 4) > 0.70);
        MU_ASSERT("air is recovered", dice(prob, truth, 5) > 0.97);

        /* the affine matrix puts the centre of the head where it belongs */
        true_affine(A);
        s[0] = s[1] = s[2] = 0.0;
        to_template(s, t);
        d = 0.0;
        for (i = 0; i < 3; i++)
        {
            double e = res.affine[i * 4 + 3] - t[i];
            d += e * e;
        }
        /* (what the affine matrix leaves, the deformation takes) */
        MU_ASSERT("the affine matrix finds the position of the head", sqrt(d) < 8.0);
        d = pow(fabs(res.affine[0] * (res.affine[5] * res.affine[10] - res.affine[6] * res.affine[9]) -
                     res.affine[1] * (res.affine[4] * res.affine[10] - res.affine[6] * res.affine[8]) +
                     res.affine[2] * (res.affine[4] * res.affine[9] - res.affine[5] * res.affine[8])),
                1.0 / 3.0);
        MU_ASSERT("the affine matrix finds the zoom of the head", fabs(d - 1.1) < 0.05);

        /* the bias field is the inverse of the one applied, up to a factor */
        sxy = sx = sy = sxx = syy = n = 0.0;
        for (k = 0; k < IZ; k++)
            for (j = 0; j < IY; j++)
                for (i = 0; i < IX; i++)
                {
                    double x, y;
                    v = i + (size_t)IX * (j + (size_t)IY * k);
                    if (truth[v] > 2) continue; /* the brain */
                    s[0] = M[0] * i + M[3];
                    s[1] = M[5] * j + M[7];
                    s[2] = M[10] * k + M[11];
                    x = -log_bias(s);
                    y = log((double)bias[v]);
                    sx += x;
                    sy += y;
                    sxx += x * x;
                    syy += y * y;
                    sxy += x * y;
                    n += 1.0;
                }
        d = (sxy - sx * sy / n) / sqrt((sxx - sx * sx / n) * (syy - sy * sy / n));
        MU_ASSERT("the bias field is recovered", d > 0.95);

        /* the label is 3 in white and 2 in grey matter */
        lw = lg = 0.0;
        nw = ng = 0;
        for (k = 0; k < IZ; k++)
            for (j = 0; j < IY; j++)
                for (i = 0; i < IX; i++)
                {
                    double p[NCLS];
                    v = i + (size_t)IX * (j + (size_t)IY * k);
                    s[0] = M[0] * i + M[3];
                    s[1] = M[5] * j + M[7];
                    s[2] = M[10] * k + M[11];
                    to_template(s, t);
                    phantom(t[0], t[1], t[2], p);
                    if (p[1] > 0.99)
                    {
                        lw += label[v];
                        nw++;
                    }
                    if (p[0] > 0.99)
                    {
                        lg += label[v];
                        ng++;
                    }
                }
        MU_ASSERT("the label of pure white matter is 3", nw > 100 && fabs(lw / (double)nw - 3.0) < 0.1);
        MU_ASSERT("the label of pure grey matter is 2", ng > 100 && fabs(lg / (double)ng - 2.0) < 0.15);
        MU_ASSERT("partial volumes were estimated", res.pve_done);
        MU_ASSERT("the bias model is reported", res.bias_model == model && res.bias_coef != NULL);
        MU_ASSERT("the volumes are positive", res.volume[0] > 0 && res.volume[1] > 0);
        CAT_UnifiedSegResultFree(&res);
    }

    free(img);
    free(label);
    free(bias);
    free(prob);
    free(truth);
    free_templates(tpl);
}

/* The result does not depend on the number of threads, and voxels without
 * data keep their prior without disturbing the rest. */
static void
test_threads_and_missing(void)
{
    float **tpl = make_templates();
    float *img = (float *)malloc(sizeof(float) * NI);
    unsigned char *p1 = (unsigned char *)malloc(NCLS * NI);
    unsigned char *p2 = (unsigned char *)malloc(NCLS * NI);
    unsigned char *truth = (unsigned char *)malloc(NI);
    CAT_UnifiedSegOpts o;
    CAT_UnifiedSegResult r1, r2;
    size_t v, bad = 0;
    int i, j, k, c, ret1, ret2;

    make_image(img, truth, 2.0, 1);
    phantom_options(&o);
    o.n_threads = 1;
    ret1 = run(img, tpl, &o, p1, NULL, NULL, &r1);
    o.n_threads = 3;
    ret2 = run(img, tpl, &o, p2, NULL, NULL, &r2);
    MU_ASSERT("both thread counts succeed", ret1 == 0 && ret2 == 0);
    if (!ret1 && !ret2)
    {
        MU_ASSERT("the class maps do not depend on the number of threads",
                  memcmp(p1, p2, NCLS * NI) == 0);
        MU_ASSERT("the parameters do not depend on the number of threads",
                  memcmp(r1.mn, r2.mn, sizeof(r1.mn)) == 0 && r1.ll == r2.ll);
    }
    CAT_UnifiedSegResultFree(&r1);
    CAT_UnifiedSegResultFree(&r2);

    /* a block without data in the brain */
    for (k = 20; k < 26; k++)
        for (j = 24; j < 30; j++)
            for (i = 20; i < 26; i++) img[i + (size_t)IX * (j + (size_t)IY * k)] = 0.0f;
    o.n_threads = 0;
    ret1 = run(img, tpl, &o, p2, NULL, NULL, &r1);
    MU_ASSERT("an image with missing voxels is segmented", ret1 == 0);
    if (!ret1)
    {
        for (v = 0; v < NI; v++)
        {
            int s = 0;
            for (c = 0; c < NCLS; c++) s += p2[c * NI + v];
            if (s < 252 || s > 258) bad++;
        }
        MU_ASSERT("the posteriors sum to one everywhere", bad == 0);
        MU_ASSERT("white matter is recovered next to missing voxels", dice(p2, truth, 1) > 0.9);
        /* the block lies in white matter: its prior says so */
        v = 23 + (size_t)IX * (27 + (size_t)IY * 23);
        MU_ASSERT("a voxel without data gets its prior", p2[1 * NI + v] > 128);
    }
    CAT_UnifiedSegResultFree(&r1);

    free(img);
    free(p1);
    free(p2);
    free(truth);
    free_templates(tpl);
}

/* Without noise the Gaussian of white matter degenerates; noise added for
 * the estimation alone gives it a width, and the segmentation is computed
 * from the image without it. */
static void
test_noise(void)
{
    float **tpl = make_templates();
    float *img = (float *)malloc(sizeof(float) * NI);
    float *label = (float *)malloc(sizeof(float) * NI);
    unsigned char *prob = (unsigned char *)malloc(NCLS * NI);
    unsigned char *truth = (unsigned char *)malloc(NI);
    CAT_UnifiedSegOpts o;
    CAT_UnifiedSegResult res;
    double sd_plain = 0.0, sd_noise = 0.0, est_plain = 0.0, add, excess;
    int ret, g;

    make_image(img, truth, 0.0, 0);
    phantom_options(&o);
    ret = run(img, tpl, &o, prob, NULL, NULL, &res);
    MU_ASSERT("the image without noise is segmented", ret == 0);
    if (!ret)
    {
        MU_ASSERT("no noise is added unless asked for", res.noise_added == 0.0);
        est_plain = res.noise_estimated;
        for (g = 0; g < res.n_gauss; g++)
            if (res.lkp[g] == 1) sd_plain = sqrt(res.vr[g]);
    }
    CAT_UnifiedSegResultFree(&res);

    o.noise = 0.05;
    ret = run(img, tpl, &o, prob, label, NULL, &res);
    MU_ASSERT("the image is segmented with added noise", ret == 0);
    if (!ret)
    {
        /* the noise of the image and the added noise make up the level */
        MU_ASSERT("noise is added up to the level asked for",
                  fabs(sqrt(res.noise_added * res.noise_added + est_plain * est_plain) - 0.05) < 1e-3);
        for (g = 0; g < res.n_gauss; g++)
            if (res.lkp[g] == 1) sd_noise = sqrt(res.vr[g]);
        MU_ASSERT("the added noise widens the Gaussian of white matter", sd_noise > 1.5 * sd_plain);
        MU_ASSERT("grey matter is recovered with added noise", dice(prob, truth, 0) > 0.9);
        MU_ASSERT("white matter is recovered with added noise", dice(prob, truth, 1) > 0.95);
        /* the segmented image has no added noise: the pure tissues of the
         * partial volume model are narrower by it than the Gaussians.  The
         * Gaussians live in bias-corrected intensities and the difference is
         * taken over the brain voxels of the lattice, so the two agree only
         * roughly (the ratio is 0.91 here). */
        add = res.noise_added * (res.intensity_range[1] - res.intensity_range[0]);
        excess = tissue_variance(&res) - res.pve_sd * res.pve_sd;
        MU_ASSERT("partial volumes are estimated", res.pve_done);
        MU_ASSERT("the pure tissues lose the added noise", fabs(excess / (add * add) - 1.0) < 0.25);
    }
    CAT_UnifiedSegResultFree(&res);

    /* an image that is noisy already gets none */
    make_image(img, truth, 6.0, 0);
    o.noise = 0.01;
    ret = run(img, tpl, &o, prob, NULL, NULL, &res);
    MU_ASSERT("a noisy image gets no added noise",
              ret == 0 && res.noise_added == 0.0 && res.noise_estimated > est_plain);
    CAT_UnifiedSegResultFree(&res);

    free(img);
    free(label);
    free(prob);
    free(truth);
    free_templates(tpl);
}

/* White matter that is brighter than its Gaussian allows -- a gyrus next to
 * a coil element, steeper than the bias field can follow -- must stay brain.
 * With a floor of 1e-4 under the priors of the classes outside the brain, as
 * SPM and CAT's templates have it, the broad Gaussian of soft tissue takes
 * such voxels from about 4.5 SD above white matter and the label map gets a
 * hole: at HOT_Z all bright voxels of this phantom were lost (mean label
 * 0.02). Without the floor they stay brain up to about 7 SD. */
static void
test_bright_white_matter(void)
{
    float **tpl = make_templates();
    float *img = (float *)malloc(sizeof(float) * NI);
    float *label = (float *)malloc(sizeof(float) * NI);
    unsigned char *prob = (unsigned char *)malloc(NCLS * NI);
    CAT_UnifiedSegOpts o;
    CAT_UnifiedSegResult res;
    double M[16], s[3], t[3], lab = 0.0, cv = 0.0;
    size_t v, nhot = 0, lost = 0;
    int i, j, k, g, ret;

    grid_matrix(IX, IY, IZ, IVOX, M);
    phantom_options(&o);
    o.samp = IVOX; /* every voxel, so that a bright one counts the same wherever it lies */

    /* the spread of white matter without bright voxels sets their brightness */
    make_image(img, NULL, 2.0, 1);
    ret = run(img, tpl, &o, prob, NULL, NULL, &res);
    MU_ASSERT("the image without bright voxels is segmented", ret == 0);
    if (!ret)
        for (g = 0; g < res.n_gauss; g++)
            if (res.lkp[g] == 1) cv = sqrt(res.vr[g]) / res.mn[g];
    CAT_UnifiedSegResultFree(&res);
    MU_ASSERT("white matter has a spread", cv > 0.0);

    hot_factor = 1.0 + HOT_Z * cv;
    make_image(img, NULL, 2.0, 1);
    hot_factor = 1.0;
    ret = run(img, tpl, &o, prob, label, NULL, &res);
    MU_ASSERT("the image with bright voxels is segmented", ret == 0);
    if (!ret)
    {
        for (k = 0; k < IZ; k++)
            for (j = 0; j < IY; j++)
                for (i = 0; i < IX; i++)
                {
                    v = i + (size_t)IX * (j + (size_t)IY * k);
                    s[0] = M[0] * i + M[3];
                    s[1] = M[5] * j + M[7];
                    s[2] = M[10] * k + M[11];
                    to_template(s, t);
                    if (!is_hot(i, j, k, t)) continue;
                    nhot++;
                    lab += label[v];
                    if (prob[0 * NI + v] + prob[1 * NI + v] + prob[2 * NI + v] < 128) lost++;
                }
        MU_ASSERT("there are bright voxels", nhot >= 3);
        MU_ASSERT("white matter that is too bright stays brain", lost == 0);
        MU_ASSERT("and keeps the label of white matter", nhot && lab / (double)nhot > 2.9);
    }
    CAT_UnifiedSegResultFree(&res);

    free(img);
    free(label);
    free(prob);
    free_templates(tpl);
}

/* With the SANLM filter the parameters still come from the image as it came:
 * they are those of a run without the filter, and only the segmentation
 * changes. */
static void
test_filter_keeps_estimation(void)
{
    float **tpl = make_templates();
    float *img = (float *)malloc(sizeof(float) * NI);
    unsigned char *p1 = (unsigned char *)malloc(NCLS * NI);
    unsigned char *p2 = (unsigned char *)malloc(NCLS * NI);
    unsigned char *truth = (unsigned char *)malloc(NI);
    float *label = (float *)malloc(sizeof(float) * NI);
    CAT_UnifiedSegOpts o;
    CAT_UnifiedSegResult r1, r2;
    int ret1, ret2;

    make_image(img, truth, 4.0, 1);
    phantom_options(&o);
    ret1 = run(img, tpl, &o, p1, label, NULL, &r1);
    o.sanlm = 1;
    ret2 = run(img, tpl, &o, p2, label, NULL, &r2);
    MU_ASSERT("the runs with and without the filter succeed", ret1 == 0 && ret2 == 0);
    if (!ret1 && !ret2)
    {
        MU_ASSERT("the filter does not change the Gaussians",
                  memcmp(r1.mn, r2.mn, sizeof(r1.mn)) == 0 && memcmp(r1.vr, r2.vr, sizeof(r1.vr)) == 0);
        MU_ASSERT("the filter does not change the affine matrix",
                  memcmp(r1.affine, r2.affine, sizeof(r1.affine)) == 0);
        MU_ASSERT("the filter changes the segmentation", memcmp(p1, p2, NCLS * NI) != 0);
        MU_ASSERT("white matter is recovered from the filtered image", dice(p2, truth, 1) > 0.97);
        /* the pure tissues of the partial volume model lie in the order and
         * about the spacing of the phantom: (65 - 30) / (100 - 30) */
        MU_ASSERT("the pure tissues are ordered CSF, GM, WM",
                  r2.pve_mean[0] < r2.pve_mean[1] && r2.pve_mean[1] < r2.pve_mean[2]);
        MU_ASSERT("grey matter lies half way between CSF and white matter",
                  fabs((r2.pve_mean[1] - r2.pve_mean[0]) / (r2.pve_mean[2] - r2.pve_mean[0]) - 0.5) < 0.12);
        /* a pure tissue is as wide as the narrowest tissue Gaussian, less
         * what the filter took out of the image that is segmented */
        MU_ASSERT("without the filter a pure tissue has the width of its Gaussian",
                  fabs(r1.pve_sd * r1.pve_sd / tissue_variance(&r1) - 1.0) < 1e-6);
        MU_ASSERT("with the filter it is narrower", r2.pve_sd < 0.97 * r1.pve_sd);
    }
    CAT_UnifiedSegResultFree(&r1);
    CAT_UnifiedSegResultFree(&r2);

    free(img);
    free(p1);
    free(p2);
    free(truth);
    free(label);
    free_templates(tpl);
}

/* The SANLM filter, SPM's number of Gaussians, the MRF and the JSON file. */
static void
test_options(void)
{
    float **tpl = make_templates();
    float *img = (float *)malloc(sizeof(float) * NI);
    unsigned char *prob = (unsigned char *)malloc(NCLS * NI);
    unsigned char *truth = (unsigned char *)malloc(NI);
    CAT_UnifiedSegOpts o;
    CAT_UnifiedSegResult res;
    char line[512];
    FILE *fp;
    int ret, found = 0;

    make_image(img, truth, 4.0, 1);
    CAT_UnifiedSegOptionsInit(&o);
    o.shoot.lowres = 6.5;
    o.samp = 7.0;
    o.sanlm = 1;
    o.mrf = 1.0;
    ret = run(img, tpl, &o, prob, NULL, NULL, &res);
    MU_ASSERT("filter, default Gaussians and MRF run", ret == 0);
    if (!ret)
    {
        MU_ASSERT("the default has 13 Gaussians", res.n_gauss == 13);
        MU_ASSERT("white matter is recovered with the filter", dice(prob, truth, 1) > 0.9);
        MU_ASSERT("the JSON file is written",
                  CAT_UnifiedSegWriteJson("test_unifiedseg_tmp.json", &res, &o, "image", "tpl") == 0);
        fp = fopen("test_unifiedseg_tmp.json", "r");
        while (fp && fgets(line, sizeof(line), fp))
            if (strstr(line, "\"mn\"")) found = 1;
        if (fp) fclose(fp);
        remove("test_unifiedseg_tmp.json");
        MU_ASSERT("the JSON file holds the means", found);
    }
    CAT_UnifiedSegResultFree(&res);

    free(img);
    free(prob);
    free(truth);
    free_templates(tpl);
}

static void
test_invalid(void)
{
    float **tpl = make_templates();
    float *img = (float *)calloc(NI, sizeof(float));
    unsigned char *prob = (unsigned char *)malloc(NCLS * NI);
    const float *t[CAT_SHOOT_NTEMPLATES];
    const int dims[3] = {IX, IY, IZ}, tdims[3] = {TX, TY, TZ};
    CAT_UnifiedSegOpts o;
    double M[16], Mt[16];
    int l;

    for (l = 0; l < CAT_SHOOT_NTEMPLATES; l++) t[l] = tpl[l];
    grid_matrix(IX, IY, IZ, IVOX, M);
    grid_matrix(TX, TY, TZ, TVOX, Mt);
    phantom_options(&o);
    MU_ASSERT("an image without data is refused",
              CAT_UnifiedSeg(img, dims, M, t, tdims, Mt, NULL, NULL, &o, prob, NULL, NULL, NULL, NULL,
                             NULL, NULL, NULL, NULL) == 3);
    MU_ASSERT("a missing image is refused",
              CAT_UnifiedSeg(NULL, dims, M, t, tdims, Mt, NULL, NULL, &o, prob, NULL, NULL, NULL, NULL,
                             NULL, NULL, NULL, NULL) == 1);
    t[2] = NULL;
    MU_ASSERT("a missing template is refused",
              CAT_UnifiedSeg(img, dims, M, t, tdims, Mt, NULL, NULL, &o, prob, NULL, NULL, NULL, NULL,
                             NULL, NULL, NULL, NULL) == 1);
    t[2] = tpl[2];
    o.ngaus[0] = 0;
    MU_ASSERT("a class without a Gaussian is refused",
              CAT_UnifiedSeg(img, dims, M, t, tdims, Mt, NULL, NULL, &o, prob, NULL, NULL, NULL, NULL,
                             NULL, NULL, NULL, NULL) == 1);
    phantom_options(&o);
    o.affine_mode = CAT_USEG_AFFINE_FIXED;
    MU_ASSERT("a fixed affine matrix must be given",
              CAT_UnifiedSeg(img, dims, M, t, tdims, Mt, NULL, NULL, &o, prob, NULL, NULL, NULL, NULL,
                             NULL, NULL, NULL, NULL) == 1);

    free(img);
    free(prob);
    free_templates(tpl);
}

int main(void)
{
    MU_RUN_TEST(test_invalid);
    MU_RUN_TEST(test_phantom);
    MU_RUN_TEST(test_threads_and_missing);
    MU_RUN_TEST(test_noise);
    MU_RUN_TEST(test_bright_white_matter);
    MU_RUN_TEST(test_filter_keeps_estimation);
    MU_RUN_TEST(test_options);
    printf("%d tests run, %d failed\n", tests_run, tests_failed);
    return tests_failed ? 1 : 0;
}
