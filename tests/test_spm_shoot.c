#include "minunit.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "spm_shoot.h"

/*
 * Properties of the SPM26 shooting primitives in 3rdparty/spm_shoot that hold
 * without MATLAB.  Their agreement with the mex files themselves is checked
 * by tests/spm_shoot_reference (make_reference.m, compare_reference.c).
 *
 * The grid has odd sizes that are not multiples of 2 or 3, so that the
 * colours of the Gauss-Seidel sweeps meet across the periodic boundary.
 */

#define NX 23
#define NY 19
#define NZ 17
#define M ((size_t)NX*NY*NZ)

static mwSize dm[3] = {NX, NY, NZ};
static const double vx[3] = {1.5, 1.2, 1.8};

/* --------------------------------------------------------------------- */

static unsigned long long rng_state = 20261026ULL;

/** \brief Uniform deviate in [-1, 1) from a fixed-seed generator. */
static double urand(void)
{
    rng_state = rng_state*6364136223846793005ULL + 1442695040888963407ULL;
    return (double)(rng_state >> 11)/(double)(1ULL << 53)*2.0 - 1.0;
}

/** \brief nc volumes of periodic noise, box-smoothed w times, max |f| = amp. */
static float *smooth_field(int nc, int w, double amp)
{
    float  *f = (float *)malloc(nc*M*sizeof(float));
    double *t = (double *)malloc(M*sizeof(double)), *u = (double *)malloc(M*sizeof(double));
    size_t  i, j, k, c, n;
    double  mx = 0.0;

    for (c = 0; c < (size_t)nc; c++)
    {
        for (i = 0; i < M; i++) t[i] = urand();
        for (n = 0; n < (size_t)w; n++)
        {
            for (k = 0; k < NZ; k++)
                for (j = 0; j < NY; j++)
                    for (i = 0; i < NX; i++)
                    {
                        size_t ip = (i + 1) % NX, im = (i + NX - 1) % NX;
                        size_t jp = (j + 1) % NY, jm = (j + NY - 1) % NY;
                        size_t kp = (k + 1) % NZ, km = (k + NZ - 1) % NZ;
                        u[i + NX*(j + NY*k)] = (t[i + NX*(j + NY*k)]*2.0
                            + t[ip + NX*(j + NY*k)] + t[im + NX*(j + NY*k)]
                            + t[i + NX*(jp + NY*k)] + t[i + NX*(jm + NY*k)]
                            + t[i + NX*(j + NY*kp)] + t[i + NX*(j + NY*km)])/8.0;
                    }
            memcpy(t, u, M*sizeof(double));
        }
        for (i = 0; i < M; i++) f[c*M + i] = (float)t[i];
    }
    for (i = 0; i < nc*M; i++) if (fabs(f[i]) > mx) mx = fabs(f[i]);
    for (i = 0; i < nc*M; i++) f[i] = (float)(f[i]*amp/mx);
    free(t);
    free(u);
    return f;
}

/** \brief Identity deformation (1-based voxel coordinates). */
static float *identity(void)
{
    float *y = (float *)malloc(3*M*sizeof(float));
    size_t i, j, k, o = 0;
    for (k = 0; k < NZ; k++)
        for (j = 0; j < NY; j++)
            for (i = 0; i < NX; i++, o++)
            {
                y[o]       = (float)(i + 1);
                y[o + M]   = (float)(j + 1);
                y[o + 2*M] = (float)(k + 1);
            }
    return y;
}

/** \brief Field of identity matrices (9 components, element (r,c) at r+3c). */
static float *identity_jac(void)
{
    float *J = (float *)calloc(9*M, sizeof(float));
    size_t i;
    for (i = 0; i < M; i++) J[i] = J[i + 4*M] = J[i + 8*M] = 1.0f;
    return J;
}

/** \brief Regularisation of CAT12's shooting, voxel sizes vx, schedule s. */
static void params(double s, const double rp[5], double prm[8], double param[8])
{
    int k;
    for (k = 0; k < 3; k++) prm[k] = vx[k];
    for (k = 0; k < 5; k++) prm[k + 3] = rp[k]*s*vx[0]*vx[1]*vx[2];
    spm_shoot_param(prm, param);
}

static const double rp_cat[5] = {1e-4, 0.001, 0.2, 0.05, 0.2};

static double maxabsdiff(const float *a, const float *b, size_t n)
{
    double d = 0.0;
    size_t i;
    for (i = 0; i < n; i++)
        if (fabs((double)a[i] - (double)b[i]) > d) d = fabs((double)a[i] - (double)b[i]);
    return d;
}

static double rel_l2(const float *a, const float *b, size_t n)
{
    double s = 0.0, r = 0.0;
    size_t i;
    for (i = 0; i < n; i++)
    {
        s += ((double)a[i] - (double)b[i])*((double)a[i] - (double)b[i]);
        r += (double)b[i]*(double)b[i];
    }
    return sqrt(s/r);
}

/* --------------------------------------------------------------------- */

/** \brief Pulling through the identity returns the input exactly. */
static void test_pull_identity(void)
{
    mwSize dm0[4] = {NX, NY, NZ, 2};
    float *f = smooth_field(2, 1, 5.0), *y = identity();
    float *g = (float *)malloc(2*M*sizeof(float));

    spm_pullc(dm0, M, 2, y, f, g);
    MU_ASSERT("pullc through the identity returns the input", memcmp(f, g, 2*M*sizeof(float)) == 0);
    spm_pull(dm0, M, 2, y, f, g);
    MU_ASSERT("pull through the identity returns the input", memcmp(f, g, 2*M*sizeof(float)) == 0);
    free(f);
    free(y);
    free(g);
}

/** \brief pushc is the adjoint of pullc: <pull f, g> = <f, push g>. */
static void test_push_adjoint(void)
{
    mwSize dm0[4] = {NX, NY, NZ, 1};
    float *v = smooth_field(3, 3, 2.5), *y = (float *)malloc(3*M*sizeof(float));
    float *f = smooth_field(1, 0, 1.0), *g = smooth_field(1, 0, 1.0);
    float *pf = (float *)malloc(M*sizeof(float)), *pg = (float *)calloc(M, sizeof(float));
    double a = 0.0, b = 0.0, s = 0.0;
    size_t i;

    spm_smalldef(dm, 1.0f, v, y);   /* up to 2.5 voxels, wrapping at the edges */
    spm_pullc(dm0, M, 1, y, f, pf);
    spm_pushc(dm0, M, 1, y, g, pg, NULL);
    for (i = 0; i < M; i++)
    {
        a += (double)pf[i]*(double)g[i];
        b += (double)f[i]*(double)pg[i];
        s += fabs((double)pf[i]*(double)g[i]);
    }
    if (fabs(a - b) > 1e-5*s)
        printf("  adjoint: <pull f,g> = %.9g, <f,push g> = %.9g\n", a, b);
    MU_ASSERT("pushc is the adjoint of pullc", fabs(a - b) <= 1e-5*s);
    free(v);
    free(y);
    free(f);
    free(g);
    free(pf);
    free(pg);
}

/** \brief A zero velocity gives the identity and identity Jacobians (det 1). */
static void test_smalldef_zero(void)
{
    float *v = (float *)calloc(3*M, sizeof(float)), *y = (float *)malloc(3*M*sizeof(float));
    float *J = (float *)malloc(9*M*sizeof(float)), *d = (float *)malloc(M*sizeof(float));
    float *id = identity(), *Jid = identity_jac();
    size_t i;
    int    ok = 1;

    spm_smalldef_jac1(dm, 0.125f, v, y, J);
    MU_ASSERT("smalldef_jac1 of zero is the identity", memcmp(y, id, 3*M*sizeof(float)) == 0);
    MU_ASSERT("smalldef_jac1 of zero has identity Jacobians", maxabsdiff(J, Jid, 9*M) == 0.0);
    spm_smalldef_jac(dm, -0.125f, v, y, J);
    MU_ASSERT("smalldef_jac of zero has identity Jacobians", maxabsdiff(J, Jid, 9*M) == 0.0);
    spm_determinant(dm, J, d);
    for (i = 0; i < M; i++) ok &= (d[i] == 1.0f);
    MU_ASSERT("determinant of identity Jacobians is 1", ok);
    spm_def2det(dm, id, d, -1);
    for (ok = 1, i = 0; i < M; i++) ok &= (d[i] == 1.0f);
    MU_ASSERT("def2det of the identity is 1", ok);
    free(v);
    free(y);
    free(J);
    free(d);
    free(id);
    free(Jid);
}

/** \brief A zero right-hand side is already converged for cgs3. */
static void test_cgs_zero_rhs(void)
{
    double prm[8], param[8];
    float *A = (float *)calloc(6*M, sizeof(float));
    float *b = (float *)calloc(3*M, sizeof(float));
    float *x = (float *)calloc(3*M, sizeof(float));
    float *r = (float *)calloc(3*M, sizeof(float));
    float *p = (float *)calloc(3*M, sizeof(float));
    float *Ap = (float *)calloc(3*M, sizeof(float));
    size_t i;
    int ok = 1;

    params(4.0, rp_cat, prm, param);
    spm_cgs3(dm, A, b, param, 1e-4, 10, x, r, p, Ap);
    for (i = 0; i < 3*M; i++)
        ok &= isfinite(x[i]) && x[i] == 0.0f;
    MU_ASSERT("cgs3 keeps the zero solution for a zero right-hand side", ok);

    free(A);
    free(b);
    free(x);
    free(r);
    free(p);
    free(Ap);
}

/** \brief Composing with the identity on either side changes nothing. */
static void test_composition_identity(void)
{
    float *v = smooth_field(3, 3, 3.0), *y = (float *)malloc(3*M*sizeof(float));
    float *J = (float *)malloc(9*M*sizeof(float)), *id = identity(), *Jid = identity_jac();
    float *C = (float *)malloc(3*M*sizeof(float)), *JC = (float *)malloc(9*M*sizeof(float));

    spm_smalldef_jac1(dm, 1.0f, v, y, J);
    spm_composition_jacobian(dm, M, y, J, id, Jid, C, JC);       /* y o id */
    MU_ASSERT("y o id = y", maxabsdiff(C, y, 3*M) < 1e-5);
    MU_ASSERT("J(id) * I = J", maxabsdiff(JC, J, 9*M) < 1e-6);
    spm_composition_jacobian(dm, M, id, Jid, y, J, C, JC);       /* id o y */
    MU_ASSERT("id o y = y", maxabsdiff(C, y, 3*M) < 1e-5);
    MU_ASSERT("I * J = J", maxabsdiff(JC, J, 9*M) < 1e-6);
    spm_composition(dm, M, y, id, C);
    MU_ASSERT("composition without Jacobians: y o id = y", maxabsdiff(C, y, 3*M) < 1e-5);
    free(v);
    free(y);
    free(J);
    free(id);
    free(Jid);
    free(C);
    free(JC);
}

/** \brief fmg3 with a zero Hessian inverts vel2mom. */
static void test_fmg_inverts_vel2mom(void)
{
    double prm[8], param[8];
    float *u = smooth_field(3, 3, 1.0), *b = (float *)malloc(3*M*sizeof(float));
    float *H = (float *)calloc(6*M, sizeof(float)), *w = (float *)calloc(3*M, sizeof(float));
    float *scratch;
    double e1, e2;
    int    k;

    params(4.0, rp_cat, prm, param);
    spm_vel2mom(dm, u, param, b);

    /* CAT12's setting: 3 cycles of full multigrid, 2 relaxation sweeps */
    scratch = (float *)calloc(spm_fmg3_scratchsize(dm, 1), sizeof(float));
    spm_fmg3(dm, H, b, param, 3, 2, w, scratch);
    e1 = rel_l2(w, u, 3*M);
    /* further V-cycles from that estimate converge */
    for (k = 0; k < 20; k++)
        spm_fmg3(dm, H, b, param, 2, 2, w, scratch);
    e2 = rel_l2(w, u, 3*M);
    if (e1 > 0.1 || e2 > 1e-4) printf("  fmg: relative error %g, after V-cycles %g\n", e1, e2);
    /* one pass of full multigrid is only an approximation (here 6%) */
    MU_ASSERT("fmg3 (3 cycles, 2 sweeps) recovers u within 10%", e1 < 0.1);
    MU_ASSERT("V-cycles converge to u within 1e-4", e2 < 1e-4);

    /* the same without a Hessian */
    memset(w, 0, 3*M*sizeof(float));
    free(scratch);
    scratch = (float *)calloc(spm_fmg3_scratchsize(dm, 0), sizeof(float));
    spm_fmg3(dm, NULL, b, param, 3, 2, w, scratch);
    MU_ASSERT("fmg3 without a Hessian recovers u within 10%", rel_l2(w, u, 3*M) < 0.1);

    free(u);
    free(b);
    free(H);
    free(w);
    free(scratch);
}

/** \brief The kernel is the impulse response of vel2mom (the layout the
    Green's function of the driver relies on). */
static void test_kernel_is_impulse_response(void)
{
    static const double rp_be[5] = {1e-4, 0.001, 0.2, 0.0, 0.0};
    double prm[8], param[8], d, kmax;
    float *K = (float *)calloc(9*M, sizeof(float)), *delta = (float *)calloc(3*M, sizeof(float));
    float *g = (float *)malloc(3*M*sizeof(float));
    int    r, c;
    size_t i;

    /* linear elasticity: 9 blocks, (r,c) at K + (r+3c)*M */
    params(4.0, rp_cat, prm, param);
    spm_kernel(dm, param, K);
    for (c = 0, d = kmax = 0.0; c < 3; c++)
    {
        memset(delta, 0, 3*M*sizeof(float));
        delta[c*M] = 1.0f;
        spm_vel2mom(dm, delta, param, g);
        for (r = 0; r < 3; r++)
            for (i = 0; i < M; i++)
            {
                double k = K[(r + 3*c)*M + i];
                if (fabs(k) > kmax) kmax = fabs(k);
                if (fabs(g[r*M + i] - k) > d) d = fabs(g[r*M + i] - k);
            }
    }
    MU_ASSERT("kernel blocks are the impulse responses of vel2mom", d < 1e-5*kmax);

    /* bending energy: one scalar kernel, scaled by 1/s[c]^2 per component */
    params(4.0, rp_be, prm, param);
    memset(K, 0, M*sizeof(float));
    spm_kernel(dm, param, K);
    for (c = 0, d = kmax = 0.0; c < 3; c++)
    {
        memset(delta, 0, 3*M*sizeof(float));
        delta[c*M] = 1.0f;
        spm_vel2mom(dm, delta, param, g);
        for (i = 0; i < M; i++)
        {
            double k = K[i]/(param[c]*param[c]);
            if (fabs(k) > kmax) kmax = fabs(k);
            if (fabs(g[c*M + i] - k) > d) d = fabs(g[c*M + i] - k);
        }
    }
    MU_ASSERT("scalar kernel is the impulse response of vel2mom", d < 1e-5*kmax);
    free(K);
    free(delta);
    free(g);
}

/** \brief B-spline coefficients interpolate the volume at the grid points. */
static void test_bspline_interpolates(void)
{
    static const int deg2[3] = {2, 2, 2}, wr[3] = {1, 1, 1};
    static const int deg3[3] = {3, 3, 3}, mir[3] = {0, 0, 0};
    float  *f = smooth_field(1, 2, 3.0);
    double *c = (double *)malloc(M*sizeof(double)), *x = (double *)malloc(3*M*sizeof(double));
    double *s = (double *)malloc(M*sizeof(double)), *s2 = (double *)malloc(M*sizeof(double));
    double *d0 = (double *)malloc(M*sizeof(double)), *d1 = (double *)malloc(M*sizeof(double));
    double *d2 = (double *)malloc(M*sizeof(double)), e;
    double  p[3] = {NAN, 2.0, 2.0}, q[3] = {2.0, 0.9, 2.0}, fp[1], gp[1];
    size_t  i, j, k, o = 0;

    for (k = 0; k < NZ; k++)
        for (j = 0; j < NY; j++)
            for (i = 0; i < NX; i++, o++)
            {
                x[o] = (double)(i + 1); x[o + M] = (double)(j + 1); x[o + 2*M] = (double)(k + 1);
            }
    MU_ASSERT("bsplinc3 succeeds", spm_bsplinc3(f, dm, deg2, wr, c) == 0);
    spm_bsplins3(c, dm, M, x, x + M, x + 2*M, deg2, wr, s, NULL, NULL, NULL);
    for (e = 0.0, i = 0; i < M; i++) if (fabs(s[i] - f[i]) > e) e = fabs(s[i] - f[i]);
    MU_ASSERT("degree 2, periodic: interpolates the grid values", e < 1e-10);
    spm_bsplins3(c, dm, M, x, x + M, x + 2*M, deg2, wr, s2, d0, d1, d2);
    MU_ASSERT("values with and without derivatives agree", memcmp(s, s2, M*sizeof(double)) == 0);

    MU_ASSERT("bsplinc3 (degree 3, mirrored) succeeds", spm_bsplinc3(f, dm, deg3, mir, c) == 0);
    spm_bsplins3(c, dm, M, x, x + M, x + 2*M, deg3, mir, s, NULL, NULL, NULL);
    for (e = 0.0, i = 0; i < M; i++) if (fabs(s[i] - f[i]) > e) e = fabs(s[i] - f[i]);
    MU_ASSERT("degree 3, mirrored: interpolates the grid values", e < 1e-10);

    /* out of bounds: NaN outside a mirrored dimension, 0 for NaN in a
       periodic one (as the arm64 mex) */
    spm_bsplins3(c, dm, 1, &q[0], &q[1], &q[2], deg3, mir, fp, gp, NULL, NULL);
    MU_ASSERT("mirrored: outside [1-0.05, n+0.05] gives NaN", isnan(fp[0]) && gp[0] == 0.0);
    MU_ASSERT("bsplinc3 rejects degree 8", spm_bsplinc3(f, dm, (int[]){8, 2, 2}, wr, c) == -1);
    spm_shoot_clear_error();
    spm_bsplinc3(f, dm, deg2, wr, c);
    spm_bsplins3(c, dm, 1, &p[0], &p[1], &p[2], deg2, wr, fp, gp, NULL, NULL);
    MU_ASSERT("periodic: a NaN coordinate gives 0", fp[0] == 0.0 && gp[0] == 0.0);
    free(f);
    free(c);
    free(x);
    free(s);
    free(s2);
    free(d0);
    free(d1);
    free(d2);
}

/** \brief Results with one and with four threads are bit-identical (and,
    except for trapprox, identical to the serial SPM code). */
static void test_threads_identical(void)
{
    static const double rp_me[5] = {1e-4, 0.001, 0.0, 0.0, 0.0};
    static const double rp_le[5] = {1e-4, 0.0, 0.0, 0.05, 0.2};
    const double *rps[3] = {rp_cat, rp_me, rp_le};
    static const int deg2[3] = {2, 2, 2}, wr[3] = {1, 1, 1};
    mwSize dm0[4] = {NX, NY, NZ, 3};
    double prm[8], param[8], tr[3];
    float *v = smooth_field(3, 3, 2.0), *H = smooth_field(6, 2, 0.3), *b = smooth_field(3, 1, 0.1);
    float *f = smooth_field(3, 1, 1.0), *y = (float *)malloc(3*M*sizeof(float));
    float *J = (float *)malloc(9*M*sizeof(float));
    float *out[3], *scratch = (float *)malloc(spm_fmg3_scratchsize(dm, 1)*sizeof(float));
    double *c[3], *x = (double *)malloc(3*M*sizeof(double)), *bs[3];
    int    mode, nt[3] = {0, 1, 4}, k, same, same0;
    size_t i, off, n = 44*M;

    for (i = 0; i < 3*M; i++) x[i] = 1.0 + (double)(i % M)*0.37 + (double)(i/M)*0.11;
    for (i = 0; i < M; i++)       /* diagonally dominant, so positive definite */
    {
        H[i]       = 0.61f + fabsf(H[i]);
        H[i + M]   = 0.61f + fabsf(H[i + M]);
        H[i + 2*M] = 0.61f + fabsf(H[i + 2*M]);
    }
    spm_smalldef_jac1(dm, 1.0f, v, y, J);

    for (mode = 0; mode < 3; mode++)
    {
        float *o;
        spm_shoot_set_threads(nt[mode]);
        spm_shoot_set_grain(0.0);   /* even this small grid on all threads */
        out[mode] = o = (float *)calloc(n, sizeof(float));
        c[mode]   = (double *)malloc(M*sizeof(double));
        bs[mode]  = (double *)malloc(4*M*sizeof(double));
        off = 0;

        params(4.0, rp_cat, prm, param);
        spm_vel2mom(dm, v, param, o + off);                              off += 3*M;
        for (k = 0; k < 3; k++)     /* relax_all, relax_me, relax_le */
        {
            params(4.0, rps[k], prm, param);
            spm_shoot_fmg_param(prm, param);
            memset(scratch, 0, spm_fmg3_scratchsize(dm, 1)*sizeof(float));
            spm_fmg3(dm, H, b, param, 3, 2, o + off, scratch);           off += 3*M;
        }
        tr[mode] = spm_trapprox(dm, H, param);
        spm_pullc(dm0, M, 3, y, f, o + off);                             off += 3*M;
        spm_pushc_grads(dm, dm, y, J, v, o + off);                       off += 3*M;
        spm_composition_jacobian(dm, M, y, J, y, J, o + off, o + off + 3*M); off += 12*M;
        spm_smalldef_jac1(dm, 0.5f, v, o + off, o + off + 3*M);          off += 12*M;
        spm_determinant(dm, J, o + off);                                 off += M;
        spm_def2det(dm, y, o + off, -1);                                 off += M;
        if (off != n) printf("  output layout is wrong\n");
        spm_bsplinc3(f, dm, deg2, wr, c[mode]);
        spm_bsplins3(c[mode], dm, M, x, x + M, x + 2*M, deg2, wr,
                     bs[mode], bs[mode] + M, bs[mode] + 2*M, bs[mode] + 3*M);
    }
    spm_shoot_set_threads(0);
    spm_shoot_set_grain(-1.0);

    same  = memcmp(out[1], out[2], n*sizeof(float)) == 0 &&
            memcmp(c[1], c[2], M*sizeof(double)) == 0 &&
            memcmp(bs[1], bs[2], 4*M*sizeof(double)) == 0 && tr[1] == tr[2];
    same0 = memcmp(out[0], out[2], n*sizeof(float)) == 0 &&
            memcmp(c[0], c[2], M*sizeof(double)) == 0 &&
            memcmp(bs[0], bs[2], 4*M*sizeof(double)) == 0;
    MU_ASSERT("1 and 4 threads give bit-identical results", same);
    MU_ASSERT("4 threads give the serial SPM results", same0);
    MU_ASSERT("trapprox with threads is within 1e-12 of the serial sum",
              fabs(tr[2] - tr[0]) <= 1e-12*fabs(tr[0]));
    for (mode = 0; mode < 3; mode++)
    {
        free(out[mode]);
        free(c[mode]);
        free(bs[mode]);
    }
    free(v);
    free(H);
    free(b);
    free(f);
    free(y);
    free(J);
    free(scratch);
    free(x);
}

/** \brief Errors are reported instead of aborting. */
static void test_errors(void)
{
    double persp[16] = {1,0,0,0.1, 0,1,0,0, 0,0,1,0, 0,0,0,1};
    float  Mf[4][3];

    spm_shoot_clear_error();
    MU_ASSERT("no error to begin with", spm_shoot_get_error() == NULL);
    spm_set_bound(7);
    MU_ASSERT("an undefined boundary condition is reported", spm_shoot_get_error() != NULL);
    spm_set_bound(0);
    MU_ASSERT("circulant boundary restored", spm_get_bound() == 0);
    spm_shoot_clear_error();
    MU_ASSERT("a perspective matrix is rejected", spm_shoot_mat(persp, Mf) == -1);
    MU_ASSERT("... and reported", spm_shoot_get_error() != NULL);
    spm_shoot_clear_error();
}

int main(void)
{
    spm_set_bound(0);
    MU_RUN_TEST(test_pull_identity);
    MU_RUN_TEST(test_push_adjoint);
    MU_RUN_TEST(test_smalldef_zero);
    MU_RUN_TEST(test_cgs_zero_rhs);
    MU_RUN_TEST(test_composition_identity);
    MU_RUN_TEST(test_fmg_inverts_vel2mom);
    MU_RUN_TEST(test_kernel_is_impulse_response);
    MU_RUN_TEST(test_bspline_interpolates);
    MU_RUN_TEST(test_threads_identical);
    MU_RUN_TEST(test_errors);
    printf("%d tests run, %d failed\n", tests_run, tests_failed);
    return tests_failed ? 1 : 0;
}
