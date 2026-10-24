/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Metric relaxation of a spherical map, as FreeSurfer's mris_sphere does it:
 * the distances between neighbouring vertices on the sphere are pulled towards
 * their distances on the folded surface, with the gradient averaged over a
 * shrinking number of neighbour iterations so that area moves at a large
 * scale first.
 *
 * Copyright Christian Gaser, University of Jena.
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "CAT_SphereRelax.h"
#include "CAT_Smooth.h"

/* averaging scales of the gradient, coarse to fine (FreeSurfer's MRISunfold) */
static const int relax_averages[] = {1024, 256, 64, 16, 4, 1, 0};

/* slope of the nonlinear area term (FreeSurfer's NEG_AREA_K) */
#define NLAREA_K 10.0

/* the first step of a scale moves no vertex further than this share of the
   radius; momentum and backtracking take it from there */
#define FIRST_STEP 0.005
#define MOMENTUM 0.9
#define MAX_BACKTRACK 12

/* Laplacian iterations that may remove the folds left after a scale */
#define UNFOLD_ITER 200

/* One level of the averaging hierarchy: a graph whose vertices are aggregates
   of the vertices of the level below (level 0 is the mesh itself). */
struct ml_level
{
    int n;          /* vertices of this level */
    int *off, *adj; /* their neighbours, CSR; NULL at level 0 (uses nbr) */
    int *agg;       /* vertex -> aggregate at the next level; NULL at the coarsest */
    double *inv;    /* 1 / members of each vertex of this level (level > 0) */
    double scale;   /* mesh vertices per vertex of this level */
    double *buf;    /* 6 n doubles of scratch: two 3-vectors per vertex */
};

typedef struct
{
    int n, n_pairs, n_faces;
    int *pi, *pj;      /* vertex pairs within n_rings, pi < pj */
    double *target;    /* their distance on the surface, scaled to the sphere */
    int *f;            /* three vertices per triangle */
    double *a0;        /* area of each triangle on the surface, scaled */
    int *n_nbr, **nbr; /* 1-ring neighbours, for the averaging */
    double l_nlarea;
    int n_levels;      /* levels of the averaging hierarchy, 0 = the mesh */
    struct ml_level *lv;
} relax_ctx;

/**
 * \brief Set the default options of the metric relaxation.
 *
 * \param opts (out) options to initialise
 */
void CAT_SphereRelaxOptionsInit(CAT_SphereRelaxOpts *opts)
{
    opts->n_rings = 2;
    opts->max_iter = 30;
    opts->l_nlarea = 1000.0;
    opts->tol = 1e-4;
    opts->verbose = 0;
}

static double tri_area(const double *x, int a, int b, int c)
{
    double u[3], v[3], w[3];
    int k;

    for (k = 0; k < 3; k++)
    {
        u[k] = x[3 * b + k] - x[3 * a + k];
        v[k] = x[3 * c + k] - x[3 * a + k];
    }
    w[0] = u[1] * v[2] - u[2] * v[1];
    w[1] = u[2] * v[0] - u[0] * v[2];
    w[2] = u[0] * v[1] - u[1] * v[0];
    return 0.5 * sqrt(w[0] * w[0] + w[1] * w[1] + w[2] * w[2]);
}

/* the vertex pairs within n_rings rings, collected ring by ring from every
   vertex; returns -1 if out of memory */
static int build_pairs(relax_ctx *c, const double *v, int n_rings, double scale)
{
    int *mark = (int *)malloc(sizeof(int) * c->n);
    int *queue = (int *)malloc(sizeof(int) * c->n);
    int cap = 16 * c->n, i;

    c->pi = (int *)malloc(sizeof(int) * cap);
    c->pj = (int *)malloc(sizeof(int) * cap);
    c->n_pairs = 0;
    if (!mark || !queue || !c->pi || !c->pj)
    {
        free(mark);
        free(queue);
        return -1;
    }
    for (i = 0; i < c->n; i++)
        mark[i] = -1;

    for (i = 0; i < c->n; i++)
    {
        int head = 0, tail = 0, ring, q;

        queue[tail++] = i;
        mark[i] = i;
        for (ring = 0; ring < n_rings; ring++)
        {
            int end = tail;
            for (; head < end; head++)
            {
                int u = queue[head], k;
                for (k = 0; k < c->n_nbr[u]; k++)
                {
                    int w = c->nbr[u][k];
                    if (mark[w] == i)
                        continue;
                    mark[w] = i;
                    queue[tail++] = w;
                }
            }
        }
        for (q = 1; q < tail; q++)
        {
            if (queue[q] < i)
                continue;
            if (c->n_pairs == cap)
            {
                int *a, *b;
                cap *= 2;
                a = (int *)realloc(c->pi, sizeof(int) * cap);
                b = a ? (int *)realloc(c->pj, sizeof(int) * cap) : NULL;
                if (a)
                    c->pi = a;
                if (!a || !b)
                {
                    free(mark);
                    free(queue);
                    return -1;
                }
                c->pj = b;
            }
            c->pi[c->n_pairs] = i;
            c->pj[c->n_pairs] = queue[q];
            c->n_pairs++;
        }
    }
    free(mark);
    free(queue);

    c->target = (double *)malloc(sizeof(double) * c->n_pairs);
    if (!c->target)
        return -1;
    for (i = 0; i < c->n_pairs; i++)
    {
        const double *a = v + 3 * c->pi[i], *b = v + 3 * c->pj[i];
        double dx = a[0] - b[0], dy = a[1] - b[1], dz = a[2] - b[2];
        c->target[i] = scale * sqrt(dx * dx + dy * dy + dz * dz);
    }
    return 0;
}

/* log(1 + exp(z)) without overflow */
static double softplus(double z)
{
    return (z > 30.0) ? z : log1p(exp(z));
}

/**
 * \brief Energy of a configuration and, if g is not NULL, its gradient.
 *
 * E = sum over pairs (d - target)^2 + l_nlarea * sum over triangles
 * a0 * softplus(-k q) / k, with q the signed area on the sphere over a0.  The
 * area term is FreeSurfer's nonlinear area term: negligible for a triangle of
 * normal size (q = 1 gives 4.5e-5 of a0) and a strong push outwards for one
 * that collapses or flips.
 */
static double energy_grad(const relax_ctx *c, const double *x, double *g)
{
    double E = 0.0;
    int p, t;

    if (g)
        memset(g, 0, sizeof(double) * 3 * c->n);
    for (p = 0; p < c->n_pairs; p++)
    {
        const int i = c->pi[p], j = c->pj[p];
        double dx = x[3 * i] - x[3 * j], dy = x[3 * i + 1] - x[3 * j + 1];
        double dz = x[3 * i + 2] - x[3 * j + 2];
        double d = sqrt(dx * dx + dy * dy + dz * dz), r = d - c->target[p];

        E += r * r;
        if (g && d > 1e-12)
        {
            double w = 2.0 * r / d;
            g[3 * i] += w * dx;
            g[3 * i + 1] += w * dy;
            g[3 * i + 2] += w * dz;
            g[3 * j] -= w * dx;
            g[3 * j + 1] -= w * dy;
            g[3 * j + 2] -= w * dz;
        }
    }
    if (c->l_nlarea <= 0.0)
        return E;

    for (t = 0; t < c->n_faces; t++)
    {
        const int i0 = c->f[3 * t], i1 = c->f[3 * t + 1], i2 = c->f[3 * t + 2];
        const double *x0 = x + 3 * i0, *x1 = x + 3 * i1, *x2 = x + 3 * i2;
        double nrm[3], u[3], v[3], cr[3], len, a, q;
        int k;

        for (k = 0; k < 3; k++)
        {
            nrm[k] = x0[k] + x1[k] + x2[k];
            u[k] = x1[k] - x0[k];
            v[k] = x2[k] - x0[k];
        }
        len = sqrt(nrm[0] * nrm[0] + nrm[1] * nrm[1] + nrm[2] * nrm[2]);
        if (len < 1e-12)
            continue;
        for (k = 0; k < 3; k++)
            nrm[k] /= len;
        cr[0] = u[1] * v[2] - u[2] * v[1];
        cr[1] = u[2] * v[0] - u[0] * v[2];
        cr[2] = u[0] * v[1] - u[1] * v[0];
        a = 0.5 * (cr[0] * nrm[0] + cr[1] * nrm[1] + cr[2] * nrm[2]);
        q = a / c->a0[t];
        E += c->l_nlarea * c->a0[t] * softplus(-NLAREA_K * q) / NLAREA_K;

        if (g && NLAREA_K * q < 50.0)
        {
            /* dE/da, times the gradient of the signed area of each corner,
               0.5 * n x (opposite edge) */
            double s = -c->l_nlarea / (1.0 + exp(NLAREA_K * q));
            const double *xs[3] = {x0, x1, x2};
            const int is[3] = {i0, i1, i2};
            for (k = 0; k < 3; k++)
            {
                const double *pa = xs[(k + 1) % 3], *pb = xs[(k + 2) % 3];
                double e[3];
                e[0] = pb[0] - pa[0];
                e[1] = pb[1] - pa[1];
                e[2] = pb[2] - pa[2];
                g[3 * is[k]] += s * 0.5 * (nrm[1] * e[2] - nrm[2] * e[1]);
                g[3 * is[k] + 1] += s * 0.5 * (nrm[2] * e[0] - nrm[0] * e[2]);
                g[3 * is[k] + 2] += s * 0.5 * (nrm[0] * e[1] - nrm[1] * e[0]);
            }
        }
    }
    return E;
}

/* triangles whose orientation is flipped against the outward direction */
static int count_folds(const relax_ctx *c, const double *x)
{
    int t, folds = 0;

    for (t = 0; t < c->n_faces; t++)
    {
        const double *x0 = x + 3 * c->f[3 * t], *x1 = x + 3 * c->f[3 * t + 1];
        const double *x2 = x + 3 * c->f[3 * t + 2];
        double u[3], v[3], cr[3];
        int k;

        for (k = 0; k < 3; k++)
        {
            u[k] = x1[k] - x0[k];
            v[k] = x2[k] - x0[k];
        }
        cr[0] = u[1] * v[2] - u[2] * v[1];
        cr[1] = u[2] * v[0] - u[0] * v[2];
        cr[2] = u[0] * v[1] - u[1] * v[0];
        if (cr[0] * (x0[0] + x1[0] + x2[0]) + cr[1] * (x0[1] + x1[1] + x2[1]) +
            cr[2] * (x0[2] + x1[2] + x2[2]) < 0.0)
            folds++;
    }
    return folds;
}

/* remove the radial component of g at every vertex */
static void tangent(const double *x, double *g, int n)
{
    int i;

    for (i = 0; i < n; i++)
    {
        const double *p = x + 3 * i;
        double *q = g + 3 * i;
        double r2 = p[0] * p[0] + p[1] * p[1] + p[2] * p[2], d;

        if (r2 < 1e-24)
            continue;
        d = (q[0] * p[0] + q[1] * p[1] + q[2] * p[2]) / r2;
        q[0] -= d * p[0];
        q[1] -= d * p[1];
        q[2] -= d * p[2];
    }
}

/* One pass of 1-ring averaging at level l: dst = mean over the vertex and its
   neighbours (FreeSurfer's MRISaverageGradients). */
static void average_pass(const relax_ctx *c, int l, const double *src, double *dst)
{
    const struct ml_level *lv = &c->lv[l];
    int i, k;

    for (i = 0; i < lv->n; i++)
    {
        double s0 = src[3 * i], s1 = src[3 * i + 1], s2 = src[3 * i + 2];
        int deg = l ? lv->off[i + 1] - lv->off[i] : c->n_nbr[i];
        const int *nb = l ? lv->adj + lv->off[i] : c->nbr[i];
        for (k = 0; k < deg; k++)
        {
            const double *q = src + 3 * nb[k];
            s0 += q[0];
            s1 += q[1];
            s2 += q[2];
        }
        dst[3 * i] = s0 / (1.0 + deg);
        dst[3 * i + 1] = s1 / (1.0 + deg);
        dst[3 * i + 2] = s2 / (1.0 + deg);
    }
}

/* n passes at level l on x, using tmp; the result ends in x */
static void average_passes(const relax_ctx *c, int l, double *x, double *tmp, int n)
{
    double *src = x, *dst = tmp, *swap;
    int pass;

    for (pass = 0; pass < n; pass++)
    {
        average_pass(c, l, src, dst);
        swap = src;
        src = dst;
        dst = swap;
    }
    if (src != x)
        memcpy(x, src, sizeof(double) * 3 * c->lv[l].n);
}

/* Average the gradient as n_avg passes of 1-ring averaging do (FreeSurfer's
   MRISaverageGradients).  Up to ML_EXACT passes run as they are.  Beyond, the
   gradient goes up the hierarchy to the coarsest level that still takes
   ML_MIN_PASSES passes -- smoothed ML_PRE passes at every level before its
   aggregates are averaged, since the gradient is almost all high frequency
   (1024 passes keep 2.6% of its norm) and would otherwise alias into the
   coarse levels -- is averaged there for the diffusion time that is left,
   scaled down by the level's mesh vertices per vertex, and comes back down,
   each level smoothed ML_POST passes.  On a T1Prep hemisphere this reproduces
   1024, 256 and 64 passes of the real gradient with r 0.987, 0.991 and 0.998,
   and the relaxation takes a few seconds instead of 46. */
#define ML_EXACT 64
#define ML_MIN_PASSES 8
#define ML_PRE 32
#define ML_POST 8

static void average_gradient(const relax_ctx *c, double *g, double *tmp, int n_avg)
{
    double spent = 0.0;
    int l, L = 0, i, k, n_coarse;

    if (n_avg <= ML_EXACT || c->n_levels < 2)
    {
        average_passes(c, 0, g, tmp, n_avg);
        return;
    }
    while (L + 1 < c->n_levels && n_avg / c->lv[L + 1].scale >= ML_MIN_PASSES)
        L++;

    /* up: smooth, then level l+1 holds the means of its aggregates */
    memcpy(c->lv[0].buf, g, sizeof(double) * 3 * c->n);
    for (l = 0; l < L; l++)
    {
        const struct ml_level *f = &c->lv[l], *q = &c->lv[l + 1];
        average_passes(c, l, f->buf, f->buf + 3 * f->n, ML_PRE);
        spent += (ML_PRE + ML_POST) * f->scale;
        memset(q->buf, 0, sizeof(double) * 3 * q->n);
        for (i = 0; i < f->n; i++)
            for (k = 0; k < 3; k++)
                q->buf[3 * f->agg[i] + k] += f->buf[3 * i + k];
        for (i = 0; i < q->n; i++)
            for (k = 0; k < 3; k++)
                q->buf[3 * i + k] *= q->inv[i];
    }
    n_coarse = (int)((n_avg - spent) / c->lv[L].scale + 0.5);
    average_passes(c, L, c->lv[L].buf, c->lv[L].buf + 3 * c->lv[L].n,
                   n_coarse > 1 ? n_coarse : 1);

    /* down: every vertex takes its aggregate's value, then is smoothed */
    for (l = L - 1; l >= 0; l--)
    {
        const struct ml_level *f = &c->lv[l], *q = &c->lv[l + 1];
        for (i = 0; i < f->n; i++)
            for (k = 0; k < 3; k++)
                f->buf[3 * i + k] = q->buf[3 * f->agg[i] + k];
        average_passes(c, l, f->buf, f->buf + 3 * f->n, ML_POST);
    }
    memcpy(g, c->lv[0].buf, sizeof(double) * 3 * c->n);
}

/* The averaging hierarchy: an unassigned vertex and its unassigned neighbours
   form an aggregate, aggregates are adjacent when members are, repeated until
   fewer than ML_COARSEST vertices are left (168k -> 46k -> 13k -> 3.5k -> ...
   on a T1Prep hemisphere).  Returns -1 if out of memory. */
#define ML_COARSEST 500
#define ML_MAX_LEVELS 16

static int build_hierarchy(relax_ctx *c)
{
    int l;

    c->lv = (struct ml_level *)calloc(ML_MAX_LEVELS, sizeof(struct ml_level));
    if (!c->lv)
        return -1;
    c->lv[0].n = c->n;
    c->lv[0].scale = 1.0;
    c->n_levels = 1;
    for (l = 0; l + 1 < ML_MAX_LEVELS && c->lv[l].n > ML_COARSEST; l++)
    {
        struct ml_level *f = &c->lv[l], *q = &c->lv[l + 1];
        int i, k, a, nq = 0, *count, *start, *members, *stamp;

        f->agg = (int *)malloc(sizeof(int) * f->n);
        if (!f->agg)
            return -1;
        for (i = 0; i < f->n; i++)
            f->agg[i] = -1;
        for (i = 0; i < f->n; i++)
        {
            int deg = l ? f->off[i + 1] - f->off[i] : c->n_nbr[i];
            const int *nb = l ? f->adj + f->off[i] : c->nbr[i];
            if (f->agg[i] >= 0)
                continue;
            f->agg[i] = nq;
            for (k = 0; k < deg; k++)
                if (f->agg[nb[k]] < 0)
                    f->agg[nb[k]] = nq;
            nq++;
        }

        /* members of every aggregate, then its neighbours without repeats */
        count = (int *)calloc(nq + 1, sizeof(int));
        start = (int *)malloc(sizeof(int) * (nq + 1));
        members = (int *)malloc(sizeof(int) * f->n);
        stamp = (int *)malloc(sizeof(int) * nq);
        q->off = (int *)malloc(sizeof(int) * (nq + 1));
        q->inv = (double *)malloc(sizeof(double) * nq);
        if (!count || !start || !members || !stamp || !q->off || !q->inv)
        {
            free(count); free(start); free(members); free(stamp);
            return -1;
        }
        for (i = 0; i < f->n; i++)
            count[f->agg[i]]++;
        start[0] = 0;
        for (a = 0; a < nq; a++)
        {
            start[a + 1] = start[a] + count[a];
            q->inv[a] = 1.0 / count[a];
        }
        for (a = 0; a < nq; a++)
            count[a] = 0;
        for (i = 0; i < f->n; i++)
            members[start[f->agg[i]] + count[f->agg[i]]++] = i;

        {
            int pass, total = 0;
            for (pass = 0; pass < 2; pass++)
            {
                for (a = 0; a < nq; a++)
                    stamp[a] = -1;
                total = 0;
                for (a = 0; a < nq; a++)
                {
                    int m;
                    if (pass)
                        q->off[a] = total;
                    for (m = start[a]; m < start[a + 1]; m++)
                    {
                        int fi = members[m];
                        int deg = l ? f->off[fi + 1] - f->off[fi] : c->n_nbr[fi];
                        const int *nb = l ? f->adj + f->off[fi] : c->nbr[fi];
                        for (k = 0; k < deg; k++)
                        {
                            int b = f->agg[nb[k]];
                            if (b == a || stamp[b] == a)
                                continue;
                            stamp[b] = a;
                            if (pass)
                                q->adj[total] = b;
                            total++;
                        }
                    }
                }
                if (!pass)
                {
                    q->adj = (int *)malloc(sizeof(int) * (total > 0 ? total : 1));
                    if (!q->adj)
                    {
                        free(count); free(start); free(members); free(stamp);
                        return -1;
                    }
                }
            }
            q->off[nq] = total;
        }
        free(count);
        free(start);
        free(members);
        free(stamp);
        q->n = nq;
        q->scale = (double)c->n / nq;
        c->n_levels++;
    }
    for (l = 0; l < c->n_levels; l++)
    {
        c->lv[l].buf = (double *)malloc(sizeof(double) * 6 * c->lv[l].n);
        if (!c->lv[l].buf)
            return -1;
    }
    return 0;
}

/* every vertex onto the sphere of the given radius */
static void project(double *x, int n, double radius)
{
    int i;

    for (i = 0; i < n; i++)
    {
        double *p = x + 3 * i;
        double r = sqrt(p[0] * p[0] + p[1] * p[1] + p[2] * p[2]);
        if (r > 1e-12)
        {
            p[0] *= radius / r;
            p[1] *= radius / r;
            p[2] *= radius / r;
        }
    }
}

/* Laplacian step on the vertices of folded triangles and their 1-ring until no
   more than target triangles are folded (or max_iter): a fold that a step of
   the descent made is pushed back by the area term, and this removes what is
   left of it.  Returns the number of folded triangles. */
static int unfold(const relax_ctx *c, double *x, double *tmp, double radius,
                  int target, int max_iter)
{
    char *mark = (char *)malloc(c->n);
    int it, t, i, k, folds = count_folds(c, x);

    if (!mark)
        return folds;
    for (it = 0; it < max_iter && folds > target; it++)
    {
        memset(mark, 0, c->n);
        for (t = 0; t < c->n_faces; t++)
        {
            const int *tri = c->f + 3 * t;
            const double *x0 = x + 3 * tri[0], *x1 = x + 3 * tri[1], *x2 = x + 3 * tri[2];
            double u[3], v[3], cr[3];

            for (k = 0; k < 3; k++)
            {
                u[k] = x1[k] - x0[k];
                v[k] = x2[k] - x0[k];
            }
            cr[0] = u[1] * v[2] - u[2] * v[1];
            cr[1] = u[2] * v[0] - u[0] * v[2];
            cr[2] = u[0] * v[1] - u[1] * v[0];
            if (cr[0] * (x0[0] + x1[0] + x2[0]) + cr[1] * (x0[1] + x1[1] + x2[1]) +
                cr[2] * (x0[2] + x1[2] + x2[2]) < 0.0)
                mark[tri[0]] = mark[tri[1]] = mark[tri[2]] = 1;
        }
        for (i = 0; i < c->n; i++)
            if (mark[i] == 1)
                for (k = 0; k < c->n_nbr[i]; k++)
                    if (!mark[c->nbr[i][k]])
                        mark[c->nbr[i][k]] = 2;
        for (i = 0; i < c->n; i++)
            if (mark[i] && c->n_nbr[i] > 0)
            {
                double s3[3] = {0.0, 0.0, 0.0};
                for (k = 0; k < c->n_nbr[i]; k++)
                {
                    const double *q = x + 3 * c->nbr[i][k];
                    s3[0] += q[0];
                    s3[1] += q[1];
                    s3[2] += q[2];
                }
                for (k = 0; k < 3; k++)
                    tmp[3 * i + k] = s3[k] / c->n_nbr[i];
            }
        for (i = 0; i < c->n; i++)
            if (mark[i] && c->n_nbr[i] > 0)
            {
                double r = sqrt(tmp[3 * i] * tmp[3 * i] + tmp[3 * i + 1] * tmp[3 * i + 1] +
                                tmp[3 * i + 2] * tmp[3 * i + 2]);
                for (k = 0; k < 3; k++)
                    x[3 * i + k] = (r > 1e-12) ? tmp[3 * i + k] * radius / r : x[3 * i + k];
            }
        folds = count_folds(c, x);
    }
    free(mark);
    return folds;
}

/* median of n values by Hoare selection on a copy; -1 if out of memory */
static double median_of(const double *values, int n)
{
    double *a = (double *)malloc(sizeof(double) * n), m;
    int lo = 0, hi = n - 1, k = n / 2;

    if (!a)
        return -1.0;
    memcpy(a, values, sizeof(double) * n);
    while (lo < hi)
    {
        double pivot = a[lo + (hi - lo) / 2];
        int i = lo, j = hi;
        while (i <= j)
        {
            while (a[i] < pivot)
                i++;
            while (a[j] > pivot)
                j--;
            if (i <= j)
            {
                double t = a[i];
                a[i] = a[j];
                a[j] = t;
                i++;
                j--;
            }
        }
        if (k <= j)
            hi = j;
        else if (k >= i)
            lo = i;
        else
            break;
    }
    m = a[k];
    free(a);
    return m;
}

static void free_ctx(relax_ctx *c)
{
    int l;

    if (c->lv)
    {
        for (l = 0; l < ML_MAX_LEVELS; l++)
        {
            free(c->lv[l].off);
            free(c->lv[l].adj);
            free(c->lv[l].agg);
            free(c->lv[l].inv);
            free(c->lv[l].buf);
        }
        free(c->lv);
    }
    free(c->pi);
    free(c->pj);
    free(c->target);
    free(c->f);
    free(c->a0);
    if (c->n_nbr)
        delete_all_polygon_point_neighbours(c->n_nbr, c->nbr);
}

/**
 * \brief Reduce the metric distortion of a spherical map of a surface.
 *
 * Gradient descent with momentum on E (see energy_grad()), on the sphere: the
 * gradient is made tangent, averaged over relax_averages[] neighbour
 * iterations, and a step is accepted when it lowers the energy, else it is
 * halved.  A scale ends when no step is accepted, after opts->max_iter steps,
 * or when the energy falls by less than opts->tol; folds left after it are
 * removed by a local Laplacian (unfold()).  As in FreeSurfer, a step may fold a
 * triangle -- the area term makes that expensive and pushes it back -- since
 * refusing every such step stalled the descent: inflated spheres of T1Prep's
 * meshes carry ~200 triangles crushed below 1% of their area, which flip under
 * any step (energy -3% instead of -55%).  The original distances are Euclidean distances on the folded
 * surface between vertices up to opts->n_rings apart, as FreeSurfer's distance
 * term uses them, scaled by sqrt(area of the sphere / area of the surface).
 *
 * \param surface (in)     the surface the sphere maps, same mesh as sphere
 * \param sphere  (in/out) spherical map of surface, centred on the origin
 *                         on return
 * \param opts    (in)     options; NULL for CAT_SphereRelaxOptionsInit()
 * \return 0 on success, -1 if the meshes differ, are not triangulated, or
 *         memory runs out
 */
int CAT_SphereRelaxMetric(polygons_struct *surface, polygons_struct *sphere,
                          const CAT_SphereRelaxOpts *opts)
{
    CAT_SphereRelaxOpts defaults;
    relax_ctx c;
    double *v = NULL, *x = NULL, *xn = NULL, *g = NULL, *gt = NULL, *tmp = NULL;
    double *delta = NULL, radius = 0.0, centre[3], lo[3], hi[3];
    double area_v = 0.0, area_x = 0.0, scale, E, median_a0;
    int i, k, s, folds0, folds, ret = -1;

    if (!opts)
    {
        CAT_SphereRelaxOptionsInit(&defaults);
        opts = &defaults;
    }
    memset(&c, 0, sizeof(c));
    if (surface->n_points != sphere->n_points || surface->n_items != sphere->n_items)
    {
        fprintf(stderr, "CAT_SphereRelaxMetric: surface and sphere differ\n");
        return -1;
    }
    for (i = 0; i < sphere->n_items; i++)
        if (GET_OBJECT_SIZE(*sphere, i) != 3)
        {
            fprintf(stderr, "CAT_SphereRelaxMetric: mesh is not triangulated\n");
            return -1;
        }

    c.n = sphere->n_points;
    c.n_faces = sphere->n_items;
    c.l_nlarea = opts->l_nlarea;
    v = (double *)malloc(sizeof(double) * 3 * c.n);
    x = (double *)malloc(sizeof(double) * 3 * c.n);
    xn = (double *)malloc(sizeof(double) * 3 * c.n);
    g = (double *)malloc(sizeof(double) * 3 * c.n);
    gt = (double *)malloc(sizeof(double) * 3 * c.n);
    tmp = (double *)malloc(sizeof(double) * 3 * c.n);
    delta = (double *)malloc(sizeof(double) * 3 * c.n);
    c.f = (int *)malloc(sizeof(int) * 3 * c.n_faces);
    c.a0 = (double *)malloc(sizeof(double) * c.n_faces);
    if (!v || !x || !xn || !g || !gt || !tmp || !delta || !c.f || !c.a0)
        goto done;

    /* the centre of the sphere is the centre of its bounding box; the vertex
       centroid is not, wherever the vertices are spread unevenly -- as on
       every distorted map -- and re-projecting about it would move them */
    for (k = 0; k < 3; k++)
    {
        lo[k] = 1e30;
        hi[k] = -1e30;
    }
    for (i = 0; i < c.n; i++)
        for (k = 0; k < 3; k++)
        {
            v[3 * i + k] = Point_coord(surface->points[i], k);
            x[3 * i + k] = Point_coord(sphere->points[i], k);
            lo[k] = fmin(lo[k], x[3 * i + k]);
            hi[k] = fmax(hi[k], x[3 * i + k]);
        }
    for (k = 0; k < 3; k++)
        centre[k] = 0.5 * (lo[k] + hi[k]);
    for (i = 0; i < c.n; i++)
    {
        for (k = 0; k < 3; k++)
            x[3 * i + k] -= centre[k];
        radius += sqrt(x[3 * i] * x[3 * i] + x[3 * i + 1] * x[3 * i + 1] +
                       x[3 * i + 2] * x[3 * i + 2]) / c.n;
    }
    project(x, c.n, radius);

    for (i = 0; i < c.n_faces; i++)
        for (k = 0; k < 3; k++)
            c.f[3 * i + k] = sphere->indices[POINT_INDEX(sphere->end_indices, i, k)];
    for (i = 0; i < c.n_faces; i++)
    {
        area_v += tri_area(v, c.f[3 * i], c.f[3 * i + 1], c.f[3 * i + 2]);
        area_x += tri_area(x, c.f[3 * i], c.f[3 * i + 1], c.f[3 * i + 2]);
    }
    if (area_v <= 0.0 || area_x <= 0.0)
        goto done;
    scale = sqrt(area_x / area_v);

    /* scaled surface area of every triangle, floored so that a degenerate
       triangle of the surface does not make the area term blow up */
    for (i = 0; i < c.n_faces; i++)
        c.a0[i] = scale * scale * tri_area(v, c.f[3 * i], c.f[3 * i + 1], c.f[3 * i + 2]);
    median_a0 = median_of(c.a0, c.n_faces);
    if (median_a0 < 0.0)
        goto done;
    for (i = 0; i < c.n_faces; i++)
        if (c.a0[i] < 1e-3 * median_a0)
            c.a0[i] = 1e-3 * median_a0;

    get_all_polygon_point_neighbours(sphere, &c.n_nbr, &c.nbr);
    if (build_pairs(&c, v, opts->n_rings < 1 ? 1 : opts->n_rings, scale) != 0 ||
        build_hierarchy(&c) != 0)
        goto done;

    folds0 = count_folds(&c, x);
    E = energy_grad(&c, x, g);
    if (opts->verbose)
        printf("Sphere relaxation: %d pairs, scale %.3f, %d folded triangles, "
               "energy/pair %.4f\n", c.n_pairs, scale, folds0, E / c.n_pairs);

    for (s = 0; s < (int)(sizeof(relax_averages) / sizeof(relax_averages[0])); s++)
    {
        double alpha = -1.0;
        int it;

        memset(delta, 0, sizeof(double) * 3 * c.n);
        for (it = 0; it < opts->max_iter; it++)
        {
            double En = E, gmax = 0.0;
            int bt, accepted = 0;

            memcpy(gt, g, sizeof(double) * 3 * c.n);
            tangent(x, gt, c.n);
            average_gradient(&c, gt, tmp, relax_averages[s]);
            tangent(x, gt, c.n);
            if (alpha < 0.0)
            {
                for (i = 0; i < 3 * c.n; i++)
                    if (fabs(gt[i]) > gmax)
                        gmax = fabs(gt[i]);
                if (gmax <= 0.0)
                    break;
                alpha = FIRST_STEP * radius / gmax;
            }
            for (bt = 0; bt < MAX_BACKTRACK; bt++)
            {
                for (i = 0; i < 3 * c.n; i++)
                    xn[i] = x[i] + MOMENTUM * delta[i] - alpha * gt[i];
                project(xn, c.n, radius);
                En = energy_grad(&c, xn, NULL);
                if (En < E)
                {
                    accepted = 1;
                    break;
                }
                alpha *= 0.5;
                for (i = 0; i < 3 * c.n; i++)
                    delta[i] *= 0.5;
            }
            if (!accepted)
                break;
            for (i = 0; i < 3 * c.n; i++)
            {
                delta[i] = xn[i] - x[i];
                x[i] = xn[i];
            }
            energy_grad(&c, x, g);
            alpha *= 1.2;
            if ((E - En) < opts->tol * E)
            {
                E = En;
                it++;
                break;
            }
            E = En;
        }
        folds = count_folds(&c, x);
        if (folds > folds0)
        {
            folds = unfold(&c, x, tmp, radius, folds0, UNFOLD_ITER);
            E = energy_grad(&c, x, g);
        }
        if (opts->verbose)
            printf("  averages %4d: %2d steps, energy/pair %.4f, %d folded triangles\n",
                   relax_averages[s], it, E / c.n_pairs, folds);
    }

    for (i = 0; i < c.n; i++)
        fill_Point(sphere->points[i], x[3 * i], x[3 * i + 1], x[3 * i + 2]);
    compute_polygon_normals(sphere);
    ret = 0;

done:
    if (ret != 0)
        fprintf(stderr, "CAT_SphereRelaxMetric: out of memory or empty mesh\n");
    free_ctx(&c);
    free(v);
    free(x);
    free(xn);
    free(g);
    free(gt);
    free(tmp);
    free(delta);
    return ret;
}
