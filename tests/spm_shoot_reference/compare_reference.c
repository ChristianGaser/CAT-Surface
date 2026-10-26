/*
 * Compare the CAT-Surface port of the SPM26 shooting primitives
 * (3rdparty/spm_shoot) with the SPM26 mex files, using the fixtures written
 * by make_reference.m, for several thread settings; optionally time them.
 *
 * Usage: compare_reference DIR [--modes LIST] [--case small|real]
 *                              [--time N] [--tmodes LIST]
 *   DIR      fixtures of make_reference.m (manifest.txt, *.bin)
 *   --modes  thread modes to compare, default 0,1,8,8g: "n" runs with
 *            spm_shoot_set_threads(n), "ng" also sets the grain to 0, so
 *            that even tiny grids (the coarse multigrid levels, the small
 *            fixture) are split among all n threads
 *   --case   only the small or only the real-size fixture
 *   --time N time the main primitives on the real-size fixture instead:
 *            best of N runs for every mode of --tmodes (default 0,1,2,4,8)
 *
 * For every output it prints the number of elements that differ from the
 * mex (NaN equals NaN), the largest absolute difference, the largest
 * difference relative to max|reference|, the number of NaN mismatches, and
 * the number of elements that differ from mode 0 and from mode 1.
 *
 * Built by build.sh (no autotools, no libCAT).
 *
 * Copyright (C) 2026 Christian Gaser, University of Jena.
 * Licence: GPL-2.0-or-later.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "spm_shoot.h"

/* ---------------------------------------------------------------------- */
/* Fixtures                                                                */

typedef struct
{
    char   name[96];
    int    bytes;       /* 4: f32, 8: f64 */
    int    ndim;
    size_t dim[8];
    size_t n;
} entry;

static entry      ent[1024];
static int        nent = 0;
static const char *refdir;

static void read_manifest(void)
{
    char  path[4096], type[8];
    FILE *fp;
    int   k;

    snprintf(path, sizeof(path), "%s/manifest.txt", refdir);
    fp = fopen(path, "r");
    if (!fp)
    {
        fprintf(stderr, "cannot open %s\n", path);
        exit(2);
    }
    while (nent < 1024 && fscanf(fp, "%95s %7s %d", ent[nent].name, type, &ent[nent].ndim) == 3)
    {
        ent[nent].bytes = strcmp(type, "f64") ? 4 : 8;
        ent[nent].n = 1;
        for (k = 0; k < ent[nent].ndim && k < 8; k++)
        {
            if (fscanf(fp, "%zu", &ent[nent].dim[k]) != 1) exit(2);
            ent[nent].n *= ent[nent].dim[k];
        }
        nent++;
    }
    fclose(fp);
}

static const entry *find(const char *name)
{
    int i;
    for (i = 0; i < nent; i++)
        if (!strcmp(ent[i].name, name)) return &ent[i];
    return NULL;
}

/* Load array <pf>_<s>; bytes must match (4 = float, 8 = double) */
static void *load(const char *pf, const char *s, int bytes, size_t *n)
{
    char         name[128], path[4096];
    const entry *e;
    FILE        *fp;
    void        *buf;

    snprintf(name, sizeof(name), "%s_%s", pf, s);
    e = find(name);
    if (!e || e->bytes != bytes)
    {
        fprintf(stderr, "fixture %s missing or of wrong type\n", name);
        exit(2);
    }
    snprintf(path, sizeof(path), "%s/%s.bin", refdir, name);
    fp  = fopen(path, "rb");
    buf = malloc(e->n*(size_t)bytes);
    if (!fp || !buf || fread(buf, (size_t)bytes, e->n, fp) != e->n)
    {
        fprintf(stderr, "cannot read %s\n", path);
        exit(2);
    }
    fclose(fp);
    if (n) *n = e->n;
    return buf;
}

static float  *loadf(const char *pf, const char *s) { return (float *)load(pf, s, 4, NULL); }
static double *loadd(const char *pf, const char *s) { return (double *)load(pf, s, 8, NULL); }

static void dims(const char *pf, mwSize dm[3])
{
    char name[128];
    const entry *e;
    snprintf(name, sizeof(name), "%s_v", pf);
    e = find(name);
    if (!e)
    {
        fprintf(stderr, "no %s\n", name);
        exit(2);
    }
    dm[0] = e->dim[0];
    dm[1] = e->dim[1];
    dm[2] = e->dim[2];
}

static void *zalloc(size_t n, size_t sz)
{
    void *p = calloc(n ? n : 1, sz);
    if (!p)
    {
        fprintf(stderr, "out of memory\n");
        exit(2);
    }
    return p;
}

/* ---------------------------------------------------------------------- */
/* Tests: each fills up to 4 outputs and names their reference arrays     */

typedef struct
{
    char    ref[96];
    void   *data;
    int     bytes;
    size_t  n;
} out_t;

static void set_out(out_t *o, const char *pf, const char *ref, void *data, int bytes, size_t n)
{
    snprintf(o->ref, sizeof(o->ref), "%s_%s", pf, ref);
    o->data  = data;
    o->bytes = bytes;
    o->n     = n;
}

static int t_vel2mom(const char *pf, const char *var, out_t *o)
{
    mwSize dm[3];
    double param[8], *prm;
    float  *v, *g;
    char   s[64];
    size_t m;

    dims(pf, dm);
    m = dm[0]*dm[1]*dm[2];
    snprintf(s, sizeof(s), "prm_%s", var);
    prm = loadd(pf, s);
    spm_shoot_param(prm, param);
    v = loadf(pf, "v");
    g = (float *)zalloc(3*m, sizeof(float));
    spm_vel2mom(dm, v, param, g);
    snprintf(s, sizeof(s), "vel2mom_%s", var);
    set_out(&o[0], pf, s, g, 4, 3*m);
    free(prm);
    free(v);
    return 1;
}

static int t_fmg(const char *pf, const char *var, out_t *o)
{
    mwSize dm[3];
    double *prm;
    float  *A, *b, *u, *u0 = NULL;
    char   s[64], r[64];
    size_t m;
    int    cyc = 3, nit = 2;

    dims(pf, dm);
    m = dm[0]*dm[1]*dm[2];
    if (!strcmp(var, "vcyc") || !strcmp(var, "u0"))
    {
        prm = loadd(pf, "prm_s4");
        if (!strcmp(var, "vcyc"))
            u0 = loadf(pf, "fmg_u0");
        else
            cyc = nit = 1;
        snprintf(r, sizeof(r), "fmg_%s", var);
    }
    else
    {
        snprintf(s, sizeof(s), "prm_%s", var);
        prm = loadd(pf, s);
        snprintf(r, sizeof(r), "fmg_%s", var);
    }
    A = loadf(pf, "A");
    b = loadf(pf, "b");
    u = (float *)zalloc(3*m, sizeof(float));
    if (spm_shoot_fmg(dm, A, b, prm, cyc, nit, u0, u))
    {
        fprintf(stderr, "spm_shoot_fmg failed: %s\n", spm_shoot_get_error());
        exit(2);
    }
    set_out(&o[0], pf, r, u, 4, 3*m);
    free(prm);
    free(A);
    free(b);
    free(u0);
    return 1;
}

static int t_mom2vel(const char *pf, const char *var, out_t *o)
{
    mwSize dm[3];
    double param[8], *prm;
    float  *b, *u, *scratch;
    char   s[64];
    size_t m;

    dims(pf, dm);
    m = dm[0]*dm[1]*dm[2];
    snprintf(s, sizeof(s), "prm_%s", strcmp(var, "all") ? var : "s4");
    prm = loadd(pf, s);
    spm_shoot_param(prm, param);    /* no extra penalty in 'mom2vel' */
    b = loadf(pf, "b");
    u = (float *)zalloc(3*m, sizeof(float));
    scratch = (float *)zalloc(spm_fmg3_scratchsize(dm, 0), sizeof(float));
    spm_fmg3(dm, NULL, b, param, 3, 2, u, scratch);
    snprintf(s, sizeof(s), "mom2vel_%s", var);
    set_out(&o[0], pf, s, u, 4, 3*m);
    free(prm);
    free(b);
    free(scratch);
    return 1;
}

static int t_trapprox(const char *pf, const char *var, out_t *o)
{
    mwSize dm[3];
    double param[8], *prm, *t;
    float  *A;
    char   s[64];

    dims(pf, dm);
    snprintf(s, sizeof(s), "prm_%s", var);
    prm = loadd(pf, s);
    spm_shoot_param(prm, param);
    A = loadf(pf, "A");
    t = (double *)zalloc(2, sizeof(double));
    t[0] = spm_trapprox(dm, A, param);
    t[1] = (double)(dm[0]*dm[1]*dm[2]*3) - t[0];
    snprintf(s, sizeof(s), "trapprox_%s", var);
    set_out(&o[0], pf, s, t, 8, 2);
    free(prm);
    free(A);
    return 1;
}

static int t_kernel(const char *pf, const char *var, out_t *o)
{
    mwSize dm[3];
    double param[8], *prm;
    float  *F;
    size_t m, nF;

    dims(pf, dm);
    m   = dm[0]*dm[1]*dm[2];
    prm = loadd(pf, strcmp(var, "be") ? "prm_s4" : "prm_be");
    spm_shoot_param(prm, param);
    nF  = (param[6] == 0.0 && param[7] == 0.0) ? m : 9*m;
    F   = (float *)zalloc(nF, sizeof(float));
    spm_kernel(dm, param, F);
    set_out(&o[0], pf, strcmp(var, "be") ? "kernel" : "kernel_be", F, 4, nF);
    free(prm);
    return 1;
}

static int t_smalldef(const char *pf, const char *var, out_t *o)
{
    mwSize dm[3];
    float  *v, *y, *J = NULL;
    size_t m;

    dims(pf, dm);
    m = dm[0]*dm[1]*dm[2];
    v = loadf(pf, "v");
    y = (float *)zalloc(3*m, sizeof(float));
    if (!strcmp(var, "1"))
    {
        spm_smalldef(dm, (float)0.37, v, y);
        set_out(&o[0], pf, "smalldef1_y", y, 4, 3*m);
        free(v);
        return 1;
    }
    J = (float *)zalloc(9*m, sizeof(float));
    spm_smalldef_jac1(dm, (float)(strcmp(var, "m") ? 1.0/8.0 : -1.0/8.0), v, y, J);
    set_out(&o[0], pf, strcmp(var, "m") ? "smalldef_y" : "smalldefm_y", y, 4, 3*m);
    set_out(&o[1], pf, strcmp(var, "m") ? "smalldef_J" : "smalldefm_J", J, 4, 9*m);
    free(v);
    return 2;
}

static int t_comp(const char *pf, const char *var, out_t *o)
{
    mwSize dm[3];
    float  *dp, *dJ, *phi, *Jphi, *C, *JC;
    size_t m;
    int    nan = !strcmp(var, "nan");

    dims(pf, dm);
    m  = dm[0]*dm[1]*dm[2];
    C  = (float *)zalloc(3*m, sizeof(float));
    if (!strcmp(var, "2"))
    {   /* spm_diffeo('comp', phi, dpm) */
        phi = loadf(pf, "phi");
        dp  = loadf(pf, "smalldefm_y");
        spm_composition(dm, m, phi, dp, C);
        set_out(&o[0], pf, "comp2_y", C, 4, 3*m);
        free(phi);
        free(dp);
        return 1;
    }
    /* spm_diffeo('comp', dp, phi, dJ, Jphi) */
    dp   = loadf(pf, "smalldef_y");
    dJ   = loadf(pf, "smalldef_J");
    phi  = loadf(pf, nan ? "phin" : "phi");
    Jphi = loadf(pf, "Jphi");
    JC   = (float *)zalloc(9*m, sizeof(float));
    spm_composition_jacobian(dm, m, dp, dJ, phi, Jphi, C, JC);
    set_out(&o[0], pf, nan ? "comp_nan_y" : "comp_y", C, 4, 3*m);
    set_out(&o[1], pf, nan ? "comp_nan_J" : "comp_J", JC, 4, 9*m);
    free(dp);
    free(dJ);
    free(phi);
    free(Jphi);
    return 2;
}

static int t_pushg(const char *pf, const char *var, out_t *o)
{
    mwSize dm[3];
    float  *m0, *phi, *Jphi = NULL, *po;
    size_t m;

    dims(pf, dm);
    m   = dm[0]*dm[1]*dm[2];
    m0  = loadf(pf, "m0");
    phi = loadf(pf, strcmp(var, "nan") ? "phi" : "phin");
    if (strcmp(var, "noJ")) Jphi = loadf(pf, "Jphi");
    po  = (float *)zalloc(3*m, sizeof(float));
    spm_pushc_grads(dm, dm, phi, Jphi, m0, po);
    set_out(&o[0], pf, !strcmp(var, "noJ") ? "pushg_noJ" : (!strcmp(var, "nan") ? "pushg_nan" : "pushg"),
            po, 4, 3*m);
    free(m0);
    free(phi);
    free(Jphi);
    return 1;
}

static int t_pull(const char *pf, const char *var, out_t *o)
{
    mwSize dm[3], dm0[4];
    float  *f, *phi, *f1;
    size_t m;

    dims(pf, dm);
    m   = dm[0]*dm[1]*dm[2];
    dm0[0] = dm[0]; dm0[1] = dm[1]; dm0[2] = dm[2]; dm0[3] = 3;
    f   = loadf(pf, "f");
    phi = loadf(pf, strcmp(var, "nan") ? "phi" : "phin");
    f1  = (float *)zalloc(3*m, sizeof(float));
    if (!strcmp(var, "pull"))
        spm_pull(dm0, m, 3, phi, f, f1);
    else
        spm_pullc(dm0, m, 3, phi, f, f1);
    set_out(&o[0], pf, !strcmp(var, "pull") ? "pull" : (!strcmp(var, "nan") ? "pullc_nan" : "pullc"),
            f1, 4, 3*m);
    free(f);
    free(phi);
    return 1;
}

static int t_push(const char *pf, const char *var, out_t *o)
{
    mwSize dm[3], dm0[4];
    float  *f1, *phi, *f0, *c0;
    size_t m;
    int    circ = !strcmp(var, "pushc");

    dims(pf, dm);
    m   = dm[0]*dm[1]*dm[2];
    dm0[0] = dm[0]; dm0[1] = dm[1]; dm0[2] = dm[2]; dm0[3] = 3;
    f1  = loadf(pf, "pullc");
    phi = loadf(pf, "phi");
    f0  = (float *)zalloc(3*m, sizeof(float));
    c0  = (float *)zalloc(m, sizeof(float));
    if (circ)
        spm_pushc(dm0, m, 3, phi, f1, f0, c0);
    else
        spm_push(dm0, m, 3, phi, f1, f0, c0);
    set_out(&o[0], pf, circ ? "pushc_f" : "push_f", f0, 4, 3*m);
    set_out(&o[1], pf, circ ? "pushc_c" : "push_c", c0, 4, m);
    free(f1);
    free(phi);
    return 2;
}

static int t_det(const char *pf, const char *var, out_t *o)
{
    mwSize dm[3];
    float  *J, *d;
    size_t m;

    (void)var;
    dims(pf, dm);
    m = dm[0]*dm[1]*dm[2];
    J = loadf(pf, "Jphi");
    d = (float *)zalloc(m, sizeof(float));
    spm_determinant(dm, J, d);
    set_out(&o[0], pf, "det", d, 4, m);
    free(J);
    return 1;
}

static int t_def2det(const char *pf, const char *var, out_t *o)
{
    mwSize dm[3];
    float  *y, *d;
    size_t m;
    int    jac = !strcmp(var, "jac");

    dims(pf, dm);
    m = dm[0]*dm[1]*dm[2];
    y = loadf(pf, "phi");
    d = (float *)zalloc(jac ? 9*m : m, sizeof(float));
    if (jac)
        spm_def2jac(dm, y, d, -1);
    else
        spm_def2det(dm, y, d, -1);
    set_out(&o[0], pf, jac ? "def2jac" : "def2det", d, 4, jac ? 9*m : m);
    free(y);
    return 1;
}

static int t_invdef(const char *pf, const char *var, out_t *o)
{
    mwSize dm[3], di[3];
    float  *y, *iy, M1[4][3], M2[4][3];
    double *M, *d, eye[16] = {1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1};
    int    id = !strcmp(var, "id");
    size_t n;

    dims(pf, dm);
    y = loadf(pf, "phi");
    if (id)
    {
        di[0] = dm[0]; di[1] = dm[1]; di[2] = dm[2];
        spm_shoot_mat(eye, M1);
        spm_shoot_mat(eye, M2);
    }
    else
    {
        d = loadd(pf, "invdef_dim");
        di[0] = (mwSize)d[0]; di[1] = (mwSize)d[1]; di[2] = (mwSize)d[2];
        free(d);
        M = loadd(pf, "invdef_M1");
        spm_shoot_mat(M, M1);
        free(M);
        M = loadd(pf, "invdef_M2");
        spm_shoot_mat(M, M2);
        free(M);
    }
    n  = 3*di[0]*di[1]*di[2];
    iy = (float *)zalloc(n, sizeof(float));
    spm_invdef(dm, y, di, iy, M1, M2);
    set_out(&o[0], pf, id ? "invdef_id" : "invdef", iy, 4, n);
    free(y);
    return 1;
}

static int t_multiscale(const char *pf, const char *var, out_t *o)
{
    mwSize dm[3], na[3], nc[3];
    float  *f, *c, *b;
    size_t nout;
    int    k;

    dims(pf, dm);
    for (k = 0; k < 3; k++)
    {
        na[k] = dm[k];
        nc[k] = (dm[k] + 1)/2;
    }
    if (!strcmp(var, "up"))
    {   /* resize the restricted volume back up to dm */
        f = loadf(pf, "restrict");
        for (k = 0; k < 3; k++)
        {
            na[k] = (dm[k] + 1)/2;
            nc[k] = dm[k];
        }
    }
    else
        f = loadf(pf, "f");     /* its first volume */
    nout = nc[0]*nc[1]*nc[2];
    c = (float *)zalloc(nout, sizeof(float));
    b = (float *)zalloc(4*nc[0]*nc[1] + na[0]*nc[1], sizeof(float));
    if (!strcmp(var, "restrict"))
        spm_restrict_vol(na, f, nc, c, b);
    else
        spm_resize_vol(na, f, nc, c, b);
    set_out(&o[0], pf, !strcmp(var, "restrict") ? "restrict" : (!strcmp(var, "up") ? "resize_up" : "resize_dn"),
            c, 4, nout);
    free(f);
    free(b);
    return 1;
}

static int t_bsplinc(const char *pf, const char *var, out_t *o)
{
    mwSize  dm[3];
    float  *g;
    double *c;
    int     deg[3] = {2, 2, 2}, wr[3] = {1, 1, 1}, sts;
    size_t  m;
    const char *ref = "bs_c";

    dims(pf, dm);
    m = dm[0]*dm[1]*dm[2];
    if (!strcmp(var, "m"))
    {
        g = loadf(pf, "bs_gn");
        deg[0] = deg[1] = deg[2] = 3;
        wr[0] = wr[1] = wr[2] = 0;
        ref = "bs_cm";
    }
    else if (!strcmp(var, "x"))
    {
        g = loadf(pf, "bs_gn");
        deg[0] = 4; deg[1] = 1; deg[2] = 2;
        wr[0] = 1; wr[1] = 0; wr[2] = 1;
        ref = "bs_cx";
    }
    else
        g = loadf(pf, "bs_lg");
    c = (double *)zalloc(m, sizeof(double));
    sts = spm_bsplinc3(g, dm, deg, wr, c);
    if (sts)
    {
        fprintf(stderr, "spm_bsplinc3 failed (%d): %s\n", sts, spm_shoot_get_error());
        exit(2);
    }
    set_out(&o[0], pf, ref, c, 8, m);
    free(g);
    return 1;
}

static int t_bsplins(const char *pf, const char *var, out_t *o)
{
    mwSize  dm[3];
    double *c, *x, *f, *d0 = NULL, *d1 = NULL, *d2 = NULL;
    int     deg[3] = {2, 2, 2}, wr[3] = {1, 1, 1};
    size_t  n, i, j, k, l;
    const char *pre = "bs_xr";
    int     id = !strcmp(var, "id"), fonly = !strcmp(var, "fonly");
    char    s[64];

    dims(pf, dm);
    if (!strcmp(var, "m"))
    {
        c = loadd(pf, "bs_cm");
        deg[0] = deg[1] = deg[2] = 3;
        wr[0] = wr[1] = wr[2] = 0;
        pre = "bs_m";
    }
    else if (!strcmp(var, "x"))
    {
        c = loadd(pf, "bs_cx");
        deg[0] = 4; deg[1] = 1; deg[2] = 2;
        wr[0] = 1; wr[1] = 0; wr[2] = 1;
        pre = "bs_x";
    }
    else
        c = loadd(pf, "bs_c");

    if (id)
    {   /* identity grid, 1-based */
        n = dm[0]*dm[1]*dm[2];
        x = (double *)zalloc(3*n, sizeof(double));
        for (k = 0, l = 0; k < dm[2]; k++)
            for (j = 0; j < dm[1]; j++)
                for (i = 0; i < dm[0]; i++, l++)
                {
                    x[l]       = (double)(i + 1);
                    x[l + n]   = (double)(j + 1);
                    x[l + 2*n] = (double)(k + 1);
                }
        pre = "bs_id";
    }
    else
        x = (double *)load(pf, "bs_xr", 8, &n), n /= 3;

    f = (double *)zalloc(n, sizeof(double));
    if (!fonly)
    {
        d0 = (double *)zalloc(n, sizeof(double));
        d1 = (double *)zalloc(n, sizeof(double));
        d2 = (double *)zalloc(n, sizeof(double));
    }
    spm_bsplins3(c, dm, n, x, x + n, x + 2*n, deg, wr, f, d0, d1, d2);
    if (fonly)
    {
        set_out(&o[0], pf, "bs_xr_fonly", f, 8, n);
        free(c);
        free(x);
        return 1;
    }
    snprintf(s, sizeof(s), "%s_f", pre);  set_out(&o[0], pf, s, f, 8, n);
    snprintf(s, sizeof(s), "%s_d0", pre); set_out(&o[1], pf, s, d0, 8, n);
    snprintf(s, sizeof(s), "%s_d1", pre); set_out(&o[2], pf, s, d1, 8, n);
    snprintf(s, sizeof(s), "%s_d2", pre); set_out(&o[3], pf, s, d2, 8, n);
    free(c);
    free(x);
    return 4;
}

typedef int (*test_fn)(const char *pf, const char *var, out_t *o);
typedef struct
{
    const char *label;
    test_fn     fn;
    const char *var;
    int         small_only;
} test_t;

static const test_t tests[] = {
    {"vel2mom s32",           t_vel2mom,    "s32",      0},
    {"vel2mom s4",            t_vel2mom,    "s4",       0},
    {"vel2mom s1",            t_vel2mom,    "s1",       0},
    {"vel2mom me",            t_vel2mom,    "me",       1},
    {"vel2mom be",            t_vel2mom,    "be",       1},
    {"vel2mom le",            t_vel2mom,    "le",       1},
    {"kernel le",             t_kernel,     "le",       0},
    {"kernel be",             t_kernel,     "be",       1},
    {"trapprox s32",          t_trapprox,   "s32",      0},
    {"trapprox s4",           t_trapprox,   "s4",       0},
    {"trapprox s1",           t_trapprox,   "s1",       0},
    {"fmg s32",               t_fmg,        "s32",      0},
    {"fmg s4",                t_fmg,        "s4",       0},
    {"fmg s1",                t_fmg,        "s1",       0},
    {"fmg cyc1 nit1",         t_fmg,        "u0",       0},
    {"fmg V-cycles",          t_fmg,        "vcyc",     0},
    {"fmg me",                t_fmg,        "me",       1},
    {"fmg be",                t_fmg,        "be",       1},
    {"fmg le",                t_fmg,        "le",       1},
    {"mom2vel me",            t_mom2vel,    "me",       1},
    {"mom2vel be",            t_mom2vel,    "be",       1},
    {"mom2vel le",            t_mom2vel,    "le",       1},
    {"mom2vel all",           t_mom2vel,    "all",      1},
    {"restrict",              t_multiscale, "restrict", 0},
    {"resize down",           t_multiscale, "dn",       0},
    {"resize up",             t_multiscale, "up",       0},
    {"smalldef_jac1 1/8",     t_smalldef,   "p",        0},
    {"smalldef_jac1 -1/8",    t_smalldef,   "m",        0},
    {"smalldef 0.37",         t_smalldef,   "1",        0},
    {"composition_jacobian",  t_comp,       "j",        0},
    {"composition",           t_comp,       "2",        0},
    {"composition NaN",       t_comp,       "nan",      1},
    {"pushc_grads",           t_pushg,      "J",        0},
    {"pushc_grads no J",      t_pushg,      "noJ",      0},
    {"pushc_grads NaN",       t_pushg,      "nan",      1},
    {"pullc",                 t_pull,       "pullc",    0},
    {"pull",                  t_pull,       "pull",     0},
    {"pullc NaN",             t_pull,       "nan",      1},
    {"pushc",                 t_push,       "pushc",    0},
    {"push",                  t_push,       "push",     0},
    {"determinant",           t_det,        "",         0},
    {"def2det",               t_def2det,    "det",      0},
    {"def2jac",               t_def2det,    "jac",      1},
    {"invdef",                t_invdef,     "M",        0},
    {"invdef identity",       t_invdef,     "id",       0},
    {"bsplinc [2 2 2 1 1 1]", t_bsplinc,    "",         0},
    {"bsplinc [3 3 3 0 0 0]", t_bsplinc,    "m",        1},
    {"bsplinc [4 1 2 1 0 1]", t_bsplinc,    "x",        1},
    {"bsplins identity grid", t_bsplins,    "id",       0},
    {"bsplins points",        t_bsplins,    "xr",       0},
    {"bsplins points f only", t_bsplins,    "fonly",    0},
    {"bsplins [3 3 3 0 0 0]", t_bsplins,    "m",        1},
    {"bsplins [4 1 2 1 0 1]", t_bsplins,    "x",        1},
};

/* ---------------------------------------------------------------------- */
/* Comparison                                                              */

typedef struct
{
    size_t ndiff, nnan;
    double maxabs, refmax;
} cmp_t;

static double val(const void *p, int bytes, size_t i)
{
    return (bytes == 4) ? (double)((const float *)p)[i] : ((const double *)p)[i];
}

/* a == b counts NaN == NaN and -0 == +0 as equal */
static cmp_t compare(const void *a, const void *b, int bytes, size_t n)
{
    cmp_t  c = {0, 0, 0.0, 0.0};
    size_t i;

    for (i = 0; i < n; i++)
    {
        double x = val(a, bytes, i), y = val(b, bytes, i), d;
        if (isfinite(y) && fabs(y) > c.refmax) c.refmax = fabs(y);
        if (isnan(x) && isnan(y)) continue;
        if (x == y) continue;
        c.ndiff++;
        if (isnan(x) || isnan(y))
        {
            c.nnan++;
            continue;
        }
        d = fabs(x - y);
        if (d > c.maxabs) c.maxabs = d;
    }
    return c;
}

typedef struct
{
    int  n;
    int  grain0;
    char label[16];
} thrmode;

static int parse_modes(const char *s, thrmode *md, int maxn)
{
    int  k = 0;
    char buf[256], *tok;

    snprintf(buf, sizeof(buf), "%s", s);
    for (tok = strtok(buf, ","); tok && k < maxn; tok = strtok(NULL, ","))
    {
        md[k].n      = atoi(tok);
        md[k].grain0 = (strchr(tok, 'g') != NULL);
        snprintf(md[k].label, sizeof(md[k].label), "%s", tok);
        k++;
    }
    return k;
}

static void set_mode(const thrmode *md)
{
    spm_shoot_set_threads(md->n);
    spm_shoot_set_grain(md->grain0 ? 0.0 : -1.0);   /* -1: default */
}

static int run_compare(const char *only, const thrmode *md, int nmodes)
{
    const char *cases[2] = {"small", "real"};
    int  ic, it, im, k, bad = 0;

    printf("%-22s %-13s %-4s %9s %9s %10s %10s %6s | %9s %10s | %9s\n",
           "function", "output", "mode", "n", "ndiff", "max abs", "max rel", "NaN",
           "ndiff/m0", "maxabs/m0", "ndiff/m1");
    for (ic = 0; ic < 2; ic++)
    {
        const char *pf = cases[ic];
        if (only && strcmp(only, pf)) continue;
        printf("---- %s\n", pf);
        for (it = 0; it < (int)(sizeof(tests)/sizeof(tests[0])); it++)
        {
            out_t base0[4], base1[4];
            int   nb0 = 0, nb1 = 0;

            if (tests[it].small_only && strcmp(pf, "small")) continue;
            for (im = 0; im < nmodes; im++)
            {
                out_t o[4];
                int   no;

                set_mode(&md[im]);
                spm_shoot_clear_error();
                no = tests[it].fn(pf, tests[it].var, o);
                for (k = 0; k < no; k++)
                {
                    const char *oname = o[k].ref + strlen(pf) + 1;
                    void  *ref = load(pf, oname, o[k].bytes, NULL);
                    cmp_t  c = compare(o[k].data, ref, o[k].bytes, o[k].n);
                    printf("%-22s %-13s %-4s %9zu %9zu %10.3g %10.3g %6zu |",
                           k ? "" : tests[it].label, oname, md[im].label, o[k].n, c.ndiff,
                           c.maxabs, c.refmax > 0 ? c.maxabs/c.refmax : 0.0, c.nnan);
                    if (md[im].n == 0 && nb0 == k)
                    {
                        base0[k] = o[k];
                        nb0++;
                        printf(" %9s %10s |", "-", "-");
                    }
                    else if (nb0 > k)
                    {
                        cmp_t d = compare(o[k].data, base0[k].data, o[k].bytes, o[k].n);
                        printf(" %9zu %10.3g |", d.ndiff, d.maxabs);
                    }
                    else
                        printf(" %9s %10s |", "", "");
                    if (md[im].n == 1 && !md[im].grain0 && nb1 == k)
                    {
                        base1[k] = o[k];
                        nb1++;
                        printf(" %9s\n", "-");
                    }
                    else if (md[im].n >= 1 && nb1 > k)
                    {
                        cmp_t d = compare(o[k].data, base1[k].data, o[k].bytes, o[k].n);
                        printf(" %9zu\n", d.ndiff);
                        if (d.ndiff) bad++;
                    }
                    else
                        printf(" %9s\n", "");
                    free(ref);
                    if (!((md[im].n == 0 && base0[k].data == o[k].data) ||
                          (md[im].n == 1 && !md[im].grain0 && nb1 > k && base1[k].data == o[k].data)))
                        free(o[k].data);
                }
                if (spm_shoot_get_error())
                    printf("    error: %s\n", spm_shoot_get_error());
            }
            for (k = 0; k < nb0; k++) free(base0[k].data);
            for (k = 0; k < nb1; k++) free(base1[k].data);
        }
    }
    printf("%d outputs of threaded modes differ from mode 1\n", bad);
    return bad;
}

/* ---------------------------------------------------------------------- */
/* Timing on the real-size fixture                                         */

static double now(void)
{
    struct timespec t;
    clock_gettime(CLOCK_MONOTONIC, &t);
    return (double)t.tv_sec + 1e-9*(double)t.tv_nsec;
}

static void run_time(int reps, const thrmode *md, int nmodes)
{
    const char *pf = "real";
    mwSize  dm[3], dm0[4], di[3];
    size_t  m, n, i, j, k, l;
    double *prm, param[8], *Md, *dd, *coef, *xg, *bf, *bd0, *bd1, *bd2, t0, best;
    float  *v, *A, *b, *u, *g, *f, *phi, *Jphi, *f1, *m0, *po, *dp, *dJ, *C, *JC, *y, *J, *d, *iy, *lg;
    float   M1[4][3], M2[4][3];
    int     im, it, r, deg[3] = {2, 2, 2}, wr[3] = {1, 1, 1};
    const char *names[] = {"vel2mom", "fmg (cyc 3, nit 2)", "pullc (3 volumes)", "pushc_grads",
                           "composition_jacobian", "smalldef_jac1", "determinant", "def2det",
                           "invdef", "bsplinc3", "bsplins3 (+derivatives)"};
    const int ntimed = (int)(sizeof(names)/sizeof(names[0]));

    dims(pf, dm);
    m   = dm[0]*dm[1]*dm[2];
    dm0[0] = dm[0]; dm0[1] = dm[1]; dm0[2] = dm[2]; dm0[3] = 3;
    prm = loadd(pf, "prm_s4");
    spm_shoot_param(prm, param);
    v    = loadf(pf, "v");
    A    = loadf(pf, "A");
    b    = loadf(pf, "b");
    f    = loadf(pf, "f");
    phi  = loadf(pf, "phi");
    Jphi = loadf(pf, "Jphi");
    m0   = loadf(pf, "m0");
    dp   = loadf(pf, "smalldef_y");
    dJ   = loadf(pf, "smalldef_J");
    lg   = loadf(pf, "bs_lg");
    coef = loadd(pf, "bs_c");
    dd   = loadd(pf, "invdef_dim");
    di[0] = (mwSize)dd[0]; di[1] = (mwSize)dd[1]; di[2] = (mwSize)dd[2];
    Md = loadd(pf, "invdef_M1"); spm_shoot_mat(Md, M1); free(Md);
    Md = loadd(pf, "invdef_M2"); spm_shoot_mat(Md, M2); free(Md);
    g    = (float *)zalloc(3*m, sizeof(float));
    u    = (float *)zalloc(3*m, sizeof(float));
    f1   = (float *)zalloc(3*m, sizeof(float));
    po   = (float *)zalloc(3*m, sizeof(float));
    C    = (float *)zalloc(3*m, sizeof(float));
    JC   = (float *)zalloc(9*m, sizeof(float));
    y    = (float *)zalloc(3*m, sizeof(float));
    J    = (float *)zalloc(9*m, sizeof(float));
    d    = (float *)zalloc(m, sizeof(float));
    iy   = (float *)zalloc(3*di[0]*di[1]*di[2], sizeof(float));
    bf   = (double *)zalloc(m, sizeof(double));
    bd0  = (double *)zalloc(m, sizeof(double));
    bd1  = (double *)zalloc(m, sizeof(double));
    bd2  = (double *)zalloc(m, sizeof(double));
    n    = m;
    xg   = (double *)zalloc(3*n, sizeof(double));
    for (k = 0, l = 0; k < dm[2]; k++)
        for (j = 0; j < dm[1]; j++)
            for (i = 0; i < dm[0]; i++, l++)
            {
                xg[l] = (double)(i+1); xg[l+n] = (double)(j+1); xg[l+2*n] = (double)(k+1);
            }

    printf("%dx%dx%d, best of %d runs, seconds\n%-26s", (int)dm[0], (int)dm[1], (int)dm[2], reps, "function");
    for (im = 0; im < nmodes; im++) printf(" %9s", md[im].label);
    printf("\n");
    for (it = 0; it < ntimed; it++)
    {
        printf("%-26s", names[it]);
        for (im = 0; im < nmodes; im++)
        {
            set_mode(&md[im]);
            best = 1e30;
            for (r = 0; r < reps; r++)
            {
                if (it == 3) memset(po, 0, 3*m*sizeof(float));
                t0 = now();
                switch (it)
                {
                    case 0: spm_vel2mom(dm, v, param, g); break;
                    case 1: spm_shoot_fmg(dm, A, b, prm, 3, 2, NULL, u); break;
                    case 2: spm_pullc(dm0, m, 3, phi, f, f1); break;
                    case 3: spm_pushc_grads(dm, dm, phi, Jphi, m0, po); break;
                    case 4: spm_composition_jacobian(dm, m, dp, dJ, phi, Jphi, C, JC); break;
                    case 5: spm_smalldef_jac1(dm, 0.125f, v, y, J); break;
                    case 6: spm_determinant(dm, Jphi, d); break;
                    case 7: spm_def2det(dm, phi, d, -1); break;
                    case 8: spm_invdef(dm, phi, di, iy, M1, M2); break;
                    case 9: spm_bsplinc3(lg, dm, deg, wr, bf); break;
                    case 10: spm_bsplins3(coef, dm, n, xg, xg+n, xg+2*n, deg, wr, bf, bd0, bd1, bd2); break;
                }
                t0 = now() - t0;
                if (t0 < best) best = t0;
            }
            printf(" %9.4f", best);
            fflush(stdout);
        }
        printf("\n");
    }
}

int main(int argc, char **argv)
{
    thrmode      md[16], tmd[16];
    int         nmodes, ntmodes, i, reps = 0;
    const char *only = NULL;

    if (argc < 2)
    {
        fprintf(stderr, "Usage: %s DIR [--modes 0,1,8,8g] [--case small|real] "
                        "[--time N] [--tmodes 0,1,2,4,8]\n", argv[0]);
        return 1;
    }
    refdir  = argv[1];
    nmodes  = parse_modes("0,1,8,8g", md, 16);
    ntmodes = parse_modes("0,1,2,4,8", tmd, 16);
    for (i = 2; i < argc; i++)
    {
        if (!strcmp(argv[i], "--modes") && i + 1 < argc)
            nmodes = parse_modes(argv[++i], md, 16);
        else if (!strcmp(argv[i], "--tmodes") && i + 1 < argc)
            ntmodes = parse_modes(argv[++i], tmd, 16);
        else if (!strcmp(argv[i], "--case") && i + 1 < argc)
            only = argv[++i];
        else if (!strcmp(argv[i], "--time") && i + 1 < argc)
            reps = atoi(argv[++i]);
    }
    read_manifest();
    if (reps > 0)
    {
        run_time(reps, tmd, ntmodes);
        return 0;
    }
    return run_compare(only, md, nmodes) ? 1 : 0;
}
