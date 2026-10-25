/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Removal of the folds of a spherical map.  Every cluster of folded triangles
 * is re-embedded on a disk of the mesh around it, with the rim of the disk
 * fixed, by minimising the barrier energy of Garanzha et al. (2021,
 * "Foldover-free maps in 50 lines of code"): it untangles a map of a disk from
 * any start and keeps every triangle close to its shape on the surface.
 *
 * Copyright Christian Gaser, University of Jena.
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "CAT_SphereUntangle.h"
#include "CAT_Smooth.h"
#include "CAT_Math.h"

/* a triangle is bad when its signed height (twice its area over its longest
   edge) is below this share of the radius: folded, or so thin that writing the
   sphere in single precision could flip it */
#define H_MIN 1e-6

/* rings of the mesh grown around a cluster, tried in turn */
static const int patch_rings[] = {2, 4, 8, 16, 32};
#define N_PATCH_RINGS ((int)(sizeof(patch_rings) / sizeof(patch_rings[0])))

/* weight of the area term against the shape term of the energy */
#define THETA 0.5

/* reference triangles below this share of the median area, or below this
   shape quality (1 for an equilateral triangle), are blended towards an
   equilateral triangle */
#define REF_FLOOR 0.1
#define SHAPE_MIN 0.2

/* L-BFGS memory and iterations per value of eps; values of eps per patch */
#define LBFGS_M 8
#define LBFGS_ITER 200
#define MAX_STAGES 100

/* stages without fewer folded triangles after which a patch is given up */
#define ABORT_STAGES 8

/* rounds over the clusters that are left: an accepted patch leaves every
   triangle it touches good and nothing outside it moves, so a later round
   only helps where patches overlap */
#define MAX_ROUNDS 3

/* the rim of a patch must lie within this cosine of its centre (about 72
   degrees), so that its gnomonic projection stays well conditioned */
#define MIN_COS 0.3

typedef struct
{
    int n, n_faces;
    int *f;              /* three vertices per triangle */
    double *x;           /* the map, centred on the origin */
    double radius, hmin;
    int orient;          /* +1: triangles counter-clockwise seen from outside */
    double *rinv;        /* inverse of the reference edge matrix, 4 per triangle */
    double *aref;        /* reference area per triangle */
    int *n_nbr, **nbr;   /* 1-ring neighbours */
    int *vf_off, *vf;    /* triangles of every vertex, CSR */
    int *ff;             /* the triangles across the three edges of a triangle */
    int *vmark, *bmark, *fmark, *fseen, *nmark, *next, stamp;
    int *loc;            /* vertex -> its number in the current patch */
    char *moved;         /* vertices moved by an accepted patch */
} untangle_ctx;

/* a disk of the mesh: free vertices S, fixed rim B, triangles F */
typedef struct
{
    int nS, nB, nF, nL, fstamp;
    int *S, *B, *F, *comp, *L;  /* L: the rim in boundary order */
    int *lt;             /* the triangles in patch vertex numbers */
    double *rinv, *aref; /* their reference, gathered */
    double *pfix;        /* gnomonic coordinates of the rim */
    double *z;           /* gnomonic coordinates of the free vertices */
    double *saved;       /* their position before the patch is applied */
    double c[3], e1[3], e2[3];
} patch;

typedef struct
{
    int cap;
    double *g, *d, *zn, *gn, *s, *y, *rho, *alpha;
} lbfgs_buf;

static double signed_height(const untangle_ctx *c, const double *x, int t)
{
    const double *a = x + 3 * c->f[3 * t], *b = x + 3 * c->f[3 * t + 1];
    const double *d = x + 3 * c->f[3 * t + 2];
    double u[3], v[3], w[3], m[3], e, emax, nm, area;
    int k;

    for (k = 0; k < 3; k++)
    {
        u[k] = b[k] - a[k];
        v[k] = d[k] - a[k];
        m[k] = a[k] + b[k] + d[k];
    }
    w[0] = u[1] * v[2] - u[2] * v[1];
    w[1] = u[2] * v[0] - u[0] * v[2];
    w[2] = u[0] * v[1] - u[1] * v[0];
    nm = sqrt(m[0] * m[0] + m[1] * m[1] + m[2] * m[2]);
    if (nm <= 0.0)
        return 0.0;
    area = 0.5 * c->orient * (w[0] * m[0] + w[1] * m[1] + w[2] * m[2]) / nm;
    emax = u[0] * u[0] + u[1] * u[1] + u[2] * u[2];
    e = v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
    if (e > emax)
        emax = e;
    e = (d[0] - b[0]) * (d[0] - b[0]) + (d[1] - b[1]) * (d[1] - b[1]) +
        (d[2] - b[2]) * (d[2] - b[2]);
    if (e > emax)
        emax = e;
    return emax > 0.0 ? 2.0 * area / sqrt(emax) : 0.0;
}

/* the reference of every triangle: its shape on the surface, scaled to the
   sphere, in a frame of its own; degenerate triangles are blended towards an
   equilateral one so that they do not dictate slivers.  Returns -1 if out of
   memory, -2 if most of the surface is degenerate. */
static int build_reference(untangle_ctx *c, const double *v, double scale)
{
    double *r = (double *)malloc(sizeof(double) * 3 * c->n_faces);
    double amin;
    int t, k;

    if (!r)
        return -1;
    for (t = 0; t < c->n_faces; t++)
    {
        const double *p0 = v + 3 * c->f[3 * t], *p1 = v + 3 * c->f[3 * t + 1];
        const double *p2 = v + 3 * c->f[3 * t + 2];
        double a[3], b[3], la, bx = 0.0, by = 0.0;

        for (k = 0; k < 3; k++)
        {
            a[k] = scale * (p1[k] - p0[k]);
            b[k] = scale * (p2[k] - p0[k]);
        }
        la = sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
        if (la > 0.0)
        {
            bx = (a[0] * b[0] + a[1] * b[1] + a[2] * b[2]) / la;
            by = b[0] * b[0] + b[1] * b[1] + b[2] * b[2] - bx * bx;
            by = by > 0.0 ? sqrt(by) : 0.0;
        }
        r[3 * t] = la;
        r[3 * t + 1] = bx;
        r[3 * t + 2] = by;
        c->aref[t] = 0.5 * la * by;
    }
    amin = REF_FLOOR * get_median_double(c->aref, c->n_faces, 0);
    if (!(amin > 0.0))
    {
        free(r);
        return -2;
    }
    for (t = 0; t < c->n_faces; t++)
    {
        double la = r[3 * t], bx = r[3 * t + 1], by = r[3 * t + 2];
        double area = c->aref[t], l2, w, side, det;

        /* shape quality: 4 sqrt(3) area / sum of the squared edges */
        l2 = la * la + bx * bx + by * by + (bx - la) * (bx - la) + by * by;
        w = l2 > 0.0 ? 1.0 - 4.0 * sqrt(3.0) * area / l2 / SHAPE_MIN : 1.0;
        if (w < 0.0)
            w = 0.0;
        if (area < amin)
            w = 1.0;
        side = sqrt(4.0 * (area > amin ? area : amin) / sqrt(3.0));
        la = (1.0 - w) * la + w * side;
        bx = (1.0 - w) * bx + w * 0.5 * side;
        by = (1.0 - w) * by + w * 0.5 * sqrt(3.0) * side;
        /* the edge matrix [[la, bx], [0, by]] has the edges as columns */
        det = la * by;
        c->rinv[4 * t] = 1.0 / la;
        c->rinv[4 * t + 1] = -bx / det;
        c->rinv[4 * t + 2] = 0.0;
        c->rinv[4 * t + 3] = 1.0 / by;
        c->aref[t] = 0.5 * det;
    }
    free(r);
    return 0;
}

/* triangles of every vertex and the neighbour across every edge */
static int build_topology(untangle_ctx *c, polygons_struct *sphere)
{
    int *fill, i, t, k;

    get_all_polygon_point_neighbours(sphere, &c->n_nbr, &c->nbr);
    c->vf_off = (int *)calloc(c->n + 1, sizeof(int));
    c->vf = (int *)malloc(sizeof(int) * 3 * c->n_faces);
    c->ff = (int *)malloc(sizeof(int) * 3 * c->n_faces);
    fill = (int *)malloc(sizeof(int) * c->n);
    if (!c->n_nbr || !c->vf_off || !c->vf || !c->ff || !fill)
    {
        free(fill);
        return -1;
    }
    for (i = 0; i < 3 * c->n_faces; i++)
        c->vf_off[c->f[i] + 1]++;
    for (i = 0; i < c->n; i++)
    {
        c->vf_off[i + 1] += c->vf_off[i];
        fill[i] = c->vf_off[i];
    }
    for (t = 0; t < c->n_faces; t++)
        for (k = 0; k < 3; k++)
            c->vf[fill[c->f[3 * t + k]]++] = t;
    free(fill);

    for (t = 0; t < c->n_faces; t++)
        for (k = 0; k < 3; k++)
        {
            int a = c->f[3 * t + k], b = c->f[3 * t + (k + 1) % 3], j;

            c->ff[3 * t + k] = -1;
            for (j = c->vf_off[a]; j < c->vf_off[a + 1]; j++)
            {
                int s = c->vf[j];
                if (s != t && (c->f[3 * s] == b || c->f[3 * s + 1] == b ||
                               c->f[3 * s + 2] == b))
                {
                    c->ff[3 * t + k] = s;
                    break;
                }
            }
        }
    return 0;
}

static int find_root(int *parent, int i)
{
    while (parent[i] != i)
    {
        parent[i] = parent[parent[i]];
        i = parent[i];
    }
    return i;
}

/* the clusters of bad triangles that share a vertex: the vertices of cluster
   k in seeds[off[k] .. off[k + 1]); returns the number of clusters and the
   number of bad triangles in *n_bad */
static int find_clusters(untangle_ctx *c, int *parent, int *cid, int *off,
                         int *seeds, int *n_bad)
{
    int t, k, i, n_cl = 0, mark;

    for (i = 0; i < c->n; i++)
    {
        parent[i] = i;
        cid[i] = -1;
    }
    *n_bad = 0;
    mark = ++c->stamp;
    for (t = 0; t < c->n_faces; t++)
    {
        if (signed_height(c, c->x, t) >= c->hmin)
            continue;
        (*n_bad)++;
        for (k = 0; k < 3; k++)
            c->vmark[c->f[3 * t + k]] = mark;
        for (k = 1; k < 3; k++)
        {
            int a = find_root(parent, c->f[3 * t]);
            int b = find_root(parent, c->f[3 * t + k]);
            if (a != b)
                parent[b] = a;
        }
    }
    if (*n_bad == 0)
        return 0;

    /* number the roots, then give every vertex the number of its root */
    for (i = 0; i < c->n; i++)
        if (c->vmark[i] == mark)
        {
            int r = find_root(parent, i);
            if (cid[r] < 0)
                cid[r] = n_cl++;
        }
    for (i = 0; i < c->n; i++)
        if (c->vmark[i] == mark)
            cid[i] = cid[find_root(parent, i)];

    memset(off, 0, sizeof(int) * (n_cl + 1));
    for (i = 0; i < c->n; i++)
        if (c->vmark[i] == mark)
            off[cid[i] + 1]++;
    for (k = 0; k < n_cl; k++)
        off[k + 1] += off[k];
    /* parent[] is not needed any more: it serves as the fill pointer */
    for (k = 0; k < n_cl; k++)
        parent[k] = off[k];
    for (i = 0; i < c->n; i++)
        if (c->vmark[i] == mark)
            seeds[parent[cid[i]]++] = i;
    return n_cl;
}

/* the triangles that touch a free vertex */
static void collect_faces(untangle_ctx *c, patch *p)
{
    int i, j;

    p->fstamp = ++c->stamp;
    p->nF = 0;
    for (i = 0; i < p->nS; i++)
        for (j = c->vf_off[p->S[i]]; j < c->vf_off[p->S[i] + 1]; j++)
        {
            int t = c->vf[j];
            if (c->fmark[t] != p->fstamp)
            {
                c->fmark[t] = p->fstamp;
                p->F[p->nF++] = t;
            }
        }
}

/* Grows the free vertices k rings around the seed, adds the holes that the
   patch encloses (the triangles outside it that are reached from its rim and
   form a small component across edges) and collects its triangles and rim.
   Returns the Euler characteristic of the patch (1 for a disk), or -1 if out
   of memory. */
static int grow_patch(untangle_ctx *c, const int *seed, int n_seed, int k,
                      patch *p)
{
    int vs = ++c->stamp, head = 0, ring, i, j, top, n_comp, outside, cap;
    int *cstart, *csize, *cbig, bs, n_be, filled = 0;

    p->nS = 0;
    for (i = 0; i < n_seed; i++)
        if (c->vmark[seed[i]] != vs)
        {
            c->vmark[seed[i]] = vs;
            p->S[p->nS++] = seed[i];
        }
    for (ring = 0; ring < k; ring++)
    {
        int end = p->nS;
        for (; head < end; head++)
        {
            int u = p->S[head];
            for (j = 0; j < c->n_nbr[u]; j++)
            {
                int w = c->nbr[u][j];
                if (c->vmark[w] != vs)
                {
                    c->vmark[w] = vs;
                    p->S[p->nS++] = w;
                }
            }
        }
    }
    collect_faces(c, p);

    /* the components outside the patch, flooded from the triangles across
       its rim; a component that grows beyond cap is the outside */
    cstart = (int *)malloc(sizeof(int) * (3 * p->nF + 1));
    csize = (int *)malloc(sizeof(int) * (3 * p->nF + 1));
    cbig = (int *)malloc(sizeof(int) * (3 * p->nF + 1));
    if (!cstart || !csize || !cbig)
    {
        free(cstart);
        free(csize);
        free(cbig);
        return -1;
    }
    cap = 4 * p->nF + 2000;
    top = 0;
    n_comp = 0;
    {
        int seen = ++c->stamp;
        for (i = 0; i < p->nF; i++)
            for (j = 0; j < 3; j++)
            {
                int s = c->ff[3 * p->F[i] + j], q, big = 0;

                if (s < 0 || c->fmark[s] == p->fstamp || c->fseen[s] == seen)
                    continue;
                cstart[n_comp] = top;
                c->fseen[s] = seen;
                p->comp[top++] = s;
                for (q = cstart[n_comp]; q < top && !big; q++)
                {
                    int t = p->comp[q], e;
                    for (e = 0; e < 3; e++)
                    {
                        int r = c->ff[3 * t + e];
                        if (r < 0 || c->fmark[r] == p->fstamp || c->fseen[r] == seen)
                            continue;
                        c->fseen[r] = seen;
                        p->comp[top++] = r;
                        if (top - cstart[n_comp] > cap)
                        {
                            big = 1;
                            break;
                        }
                    }
                }
                csize[n_comp] = top - cstart[n_comp];
                cbig[n_comp] = big;
                n_comp++;
            }
    }
    outside = -1;
    for (i = 0; i < n_comp; i++)
        if (outside < 0 || (cbig[i] && !cbig[outside]) ||
            (cbig[i] == cbig[outside] && csize[i] > csize[outside]))
            outside = i;
    for (i = 0; i < n_comp; i++)
    {
        int q;
        if (i == outside || cbig[i])
            continue;
        for (q = cstart[i]; q < cstart[i] + csize[i]; q++)
            for (j = 0; j < 3; j++)
            {
                int v = c->f[3 * p->comp[q] + j];
                if (c->vmark[v] != vs)
                {
                    c->vmark[v] = vs;
                    p->S[p->nS++] = v;
                    filled = 1;
                }
            }
    }
    free(cstart);
    free(csize);
    free(cbig);
    if (filled)
        collect_faces(c, p);

    /* a vertex left out of the free set whose triangles all belong to the
       patch lies inside it -- the filling above only sees triangles outside
       -- and would be held where it is; free it, its triangles are in the
       patch already */
    {
        int checked = ++c->stamp, q;

        for (i = 0; i < p->nF; i++)
            for (j = 0; j < 3; j++)
            {
                int v = c->f[3 * p->F[i] + j], inside = 1;

                if (c->vmark[v] == vs || c->bmark[v] == checked)
                    continue;
                c->bmark[v] = checked;
                for (q = c->vf_off[v]; q < c->vf_off[v + 1]; q++)
                    if (c->fmark[c->vf[q]] != p->fstamp)
                    {
                        inside = 0;
                        break;
                    }
                if (inside)
                {
                    c->vmark[v] = vs;
                    p->S[p->nS++] = v;
                }
            }
    }

    /* the rim and the boundary edges: chi = V - E + F, E = (3 F + edges on
       the boundary) / 2 */
    bs = ++c->stamp;
    p->nB = 0;
    n_be = 0;
    for (i = 0; i < p->nF; i++)
        for (j = 0; j < 3; j++)
        {
            int v = c->f[3 * p->F[i] + j], s = c->ff[3 * p->F[i] + j];
            if (c->vmark[v] != vs && c->bmark[v] != bs)
            {
                c->bmark[v] = bs;
                p->B[p->nB++] = v;
            }
            if (s < 0 || c->fmark[s] != p->fstamp)
                n_be++;
        }
    return p->nS + p->nB - (3 * p->nF + n_be) / 2 + p->nF;
}

static void gnomonic(const patch *p, const double *x, double radius, double *q)
{
    double r = sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]);
    double d = (x[0] * p->c[0] + x[1] * p->c[1] + x[2] * p->c[2]) / r;
    double s;

    if (d < 1e-3)
        d = 1e-3;
    s = radius / (r * d);
    q[0] = s * (x[0] * p->e1[0] + x[1] * p->e1[1] + x[2] * p->e1[2]);
    q[1] = s * (x[0] * p->e2[0] + x[1] * p->e2[1] + x[2] * p->e2[2]);
}

/* the gnomonic plane of the patch, centred on its rim, and the patch in its
   coordinates; -1 if the rim spreads too far for the projection */
static int frame_patch(untangle_ctx *c, patch *p)
{
    double nrm, axis[3] = {0.0, 0.0, 0.0};
    int i, k, imin = 0;

    for (k = 0; k < 3; k++)
        p->c[k] = 0.0;
    for (i = 0; i < p->nB; i++)
    {
        const double *x = c->x + 3 * p->B[i];
        double r = sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]);
        for (k = 0; k < 3; k++)
            p->c[k] += x[k] / r;
    }
    nrm = sqrt(p->c[0] * p->c[0] + p->c[1] * p->c[1] + p->c[2] * p->c[2]);
    if (nrm <= 0.0)
        return -1;
    for (k = 0; k < 3; k++)
        p->c[k] /= nrm;
    for (i = 0; i < p->nB; i++)
    {
        const double *x = c->x + 3 * p->B[i];
        double r = sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]);
        if ((x[0] * p->c[0] + x[1] * p->c[1] + x[2] * p->c[2]) / r < MIN_COS)
            return -1;
    }

    /* e1 orthogonal to c, e2 = c x e1 (negated for clockwise triangles) */
    for (k = 1; k < 3; k++)
        if (fabs(p->c[k]) < fabs(p->c[imin]))
            imin = k;
    axis[imin] = 1.0;
    p->e1[0] = p->c[1] * axis[2] - p->c[2] * axis[1];
    p->e1[1] = p->c[2] * axis[0] - p->c[0] * axis[2];
    p->e1[2] = p->c[0] * axis[1] - p->c[1] * axis[0];
    nrm = sqrt(p->e1[0] * p->e1[0] + p->e1[1] * p->e1[1] + p->e1[2] * p->e1[2]);
    for (k = 0; k < 3; k++)
        p->e1[k] /= nrm;
    p->e2[0] = c->orient * (p->c[1] * p->e1[2] - p->c[2] * p->e1[1]);
    p->e2[1] = c->orient * (p->c[2] * p->e1[0] - p->c[0] * p->e1[2]);
    p->e2[2] = c->orient * (p->c[0] * p->e1[1] - p->c[1] * p->e1[0]);

    for (i = 0; i < p->nS; i++)
    {
        c->loc[p->S[i]] = i;
        gnomonic(p, c->x + 3 * p->S[i], c->radius, p->z + 2 * i);
    }
    for (i = 0; i < p->nB; i++)
    {
        c->loc[p->B[i]] = p->nS + i;
        gnomonic(p, c->x + 3 * p->B[i], c->radius, p->pfix + 2 * i);
    }
    for (i = 0; i < p->nF; i++)
    {
        int t = p->F[i];
        for (k = 0; k < 3; k++)
            p->lt[3 * i + k] = c->loc[c->f[3 * t + k]];
        for (k = 0; k < 4; k++)
            p->rinv[4 * i + k] = c->rinv[4 * t + k];
        p->aref[i] = c->aref[t];
    }
    return 0;
}

static double dot(const double *a, const double *b, int n)
{
    double s = 0.0;
    int i;

    for (i = 0; i < n; i++)
        s += a[i] * b[i];
    return s;
}

/* The rim of the patch in the order of its boundary edges (each runs along a
   patch triangle, so the patch lies to its left); -1 unless it is one loop. */
static int rim_loop(untangle_ctx *c, patch *p)
{
    int mark = ++c->stamp, i, j, n_edges = 0, start = -1, v;

    for (i = 0; i < p->nF; i++)
        for (j = 0; j < 3; j++)
        {
            int t = p->F[i], s = c->ff[3 * t + j];
            int a = c->f[3 * t + j], b = c->f[3 * t + (j + 1) % 3];

            if (s >= 0 && c->fmark[s] == p->fstamp)
                continue;
            if (c->nmark[a] == mark)
                return -1; /* the boundary passes a vertex twice */
            c->nmark[a] = mark;
            c->next[a] = b;
            start = a;
            n_edges++;
        }
    if (start < 0)
        return -1;
    p->nL = 0;
    v = start;
    do
    {
        if (p->nL == n_edges || c->nmark[v] != mark)
            return -1;
        p->L[p->nL++] = v;
        v = c->next[v];
    } while (v != start);
    return p->nL == n_edges ? 0 : -1;
}

static double orient2d(const double *a, const double *b, const double *d)
{
    return (b[0] - a[0]) * (d[1] - a[1]) - (b[1] - a[1]) * (d[0] - a[0]);
}

/* The area inside the rim if, in the gnomonic plane, the rim is a simple
   polygon around the patch (counter-clockwise, no two edges touching but
   neighbours) -- only then can the patch be embedded without folds -- and 0
   otherwise. */
static double rim_simple(const untangle_ctx *c, const patch *p)
{
    double area = 0.0;
    int i, j, n = p->nL;

    for (i = 0; i < n; i++)
    {
        const double *a = p->pfix + 2 * (c->loc[p->L[i]] - p->nS);
        const double *b = p->pfix + 2 * (c->loc[p->L[(i + 1) % n]] - p->nS);
        area += a[0] * b[1] - b[0] * a[1];
    }
    if (!(area > 0.0))
        return 0.0;
    for (i = 0; i < n; i++)
    {
        const double *a = p->pfix + 2 * (c->loc[p->L[i]] - p->nS);
        const double *b = p->pfix + 2 * (c->loc[p->L[(i + 1) % n]] - p->nS);

        for (j = i + 2; j < n; j++)
        {
            const double *d = p->pfix + 2 * (c->loc[p->L[j]] - p->nS);
            const double *e = p->pfix + 2 * (c->loc[p->L[(j + 1) % n]] - p->nS);

            if (i == 0 && j == n - 1)
                continue; /* neighbours across the start of the loop */
            if (orient2d(a, b, d) * orient2d(a, b, e) <= 0.0 &&
                orient2d(d, e, a) * orient2d(d, e, b) <= 0.0)
                return 0.0;
        }
    }
    return 0.5 * area;
}

/* Scales the reference of the patch to the area inside its rim, so that det J
   is 1 on average: the free vertices cannot cover more than that, and where
   the map is compressed by orders of magnitude eps would otherwise exceed
   every determinant and the barrier would no longer tell a fold from a
   thin triangle. */
static void normalise_patch(patch *p, double rim_area)
{
    double sum = 0.0, s;
    int i, k;

    for (i = 0; i < p->nF; i++)
        sum += p->aref[i];
    s = sqrt(rim_area / sum);
    for (i = 0; i < p->nF; i++)
    {
        for (k = 0; k < 4; k++)
            p->rinv[4 * i + k] /= s;
        p->aref[i] *= s * s;
    }
}

/* y = L z for the graph Laplacian of the free vertices, the rim held at 0 */
static void laplacian(const untangle_ctx *c, const patch *p, const double *z, double *y)
{
    int i, j;

    for (i = 0; i < p->nS; i++)
    {
        int v = p->S[i];
        double sx = c->n_nbr[v] * z[2 * i], sy = c->n_nbr[v] * z[2 * i + 1];

        for (j = 0; j < c->n_nbr[v]; j++)
        {
            int l = c->loc[c->nbr[v][j]];
            if (l < p->nS)
            {
                sx -= z[2 * l];
                sy -= z[2 * l + 1];
            }
        }
        y[2 * i] = sx;
        y[2 * i + 1] = sy;
    }
}

/* Tutte's embedding of the patch: every free vertex at the mean of its
   neighbours, the rim fixed -- fold-free for a convex rim and nearly so for
   any simple one, which makes it the start of the minimisation.  Solved by
   conjugate gradients from the current positions; r, q and Aq hold 2 nS. */
static void tutte_start(const untangle_ctx *c, patch *p, double *r, double *q,
                        double *Aq)
{
    int n = 2 * p->nS, i, j, it;
    double rr, bb = 0.0;

    /* r = b - L z, b = the sum over the fixed neighbours */
    laplacian(c, p, p->z, Aq);
    for (i = 0; i < p->nS; i++)
    {
        int v = p->S[i];
        double bx = 0.0, by = 0.0;

        for (j = 0; j < c->n_nbr[v]; j++)
        {
            int l = c->loc[c->nbr[v][j]];
            if (l >= p->nS)
            {
                bx += p->pfix[2 * (l - p->nS)];
                by += p->pfix[2 * (l - p->nS) + 1];
            }
        }
        r[2 * i] = bx - Aq[2 * i];
        r[2 * i + 1] = by - Aq[2 * i + 1];
        bb += bx * bx + by * by;
    }
    memcpy(q, r, sizeof(double) * n);
    rr = dot(r, r, n);
    for (it = 0; it < 5 * n && rr > 1e-24 * bb; it++)
    {
        double alpha, beta, rr_new;

        laplacian(c, p, q, Aq);
        alpha = rr / dot(q, Aq, n);
        for (i = 0; i < n; i++)
        {
            p->z[i] += alpha * q[i];
            r[i] -= alpha * Aq[i];
        }
        rr_new = dot(r, r, n);
        beta = rr_new / rr;
        rr = rr_new;
        for (i = 0; i < n; i++)
            q[i] = r[i] + beta * q[i];
    }
}

static double chi(double D, double eps)
{
    double r = sqrt(eps * eps + D * D);
    return D >= 0.0 ? 0.5 * (D + r) : 0.5 * eps * eps / (r - D);
}

static const double *coord(const patch *p, const double *z, int i)
{
    return i < p->nS ? z + 2 * i : p->pfix + 2 * (i - p->nS);
}

/* Jacobian of triangle i from its reference to the patch; returns its
   determinant */
static double jacobian(const patch *p, const double *z, int i, double J[4])
{
    const double *a = coord(p, z, p->lt[3 * i]), *b = coord(p, z, p->lt[3 * i + 1]);
    const double *d = coord(p, z, p->lt[3 * i + 2]), *R = p->rinv + 4 * i;
    double P00 = b[0] - a[0], P10 = b[1] - a[1], P01 = d[0] - a[0], P11 = d[1] - a[1];

    J[0] = P00 * R[0] + P01 * R[2];
    J[1] = P00 * R[1] + P01 * R[3];
    J[2] = P10 * R[0] + P11 * R[2];
    J[3] = P10 * R[1] + P11 * R[3];
    return J[0] * J[3] - J[1] * J[2];
}

static double min_det(const patch *p, const double *z)
{
    double J[4], m = HUGE_VAL;
    int i;

    for (i = 0; i < p->nF; i++)
    {
        double D = jacobian(p, z, i, J);
        if (D < m)
            m = D;
    }
    return m;
}

static int count_neg(const patch *p, const double *z)
{
    double J[4];
    int i, n = 0;

    for (i = 0; i < p->nF; i++)
        n += jacobian(p, z, i, J) <= 0.0;
    return n;
}

static void add_grad(const patch *p, double *g, int i, double gx, double gy)
{
    if (i < p->nS)
    {
        g[2 * i] += gx;
        g[2 * i + 1] += gy;
    }
}

/* E = sum_t A_ref [(1 - THETA) |J|^2 + THETA (det J^2 + 1)] / chi(det J, eps)
   and its gradient with respect to the free coordinates (if g); +inf where
   chi underflows */
static double patch_energy(const patch *p, const double *z, double eps, double *g)
{
    double E = 0.0, J[4];
    int i;

    if (g)
        memset(g, 0, sizeof(double) * 2 * p->nS);
    for (i = 0; i < p->nF; i++)
    {
        double D = jacobian(p, z, i, J), r = sqrt(eps * eps + D * D);
        double ch = chi(D, eps), F, w = p->aref[i];
        const double *R = p->rinv + 4 * i;

        if (!(ch > 1e-300))
            return HUGE_VAL;
        F = (1.0 - THETA) * (J[0] * J[0] + J[1] * J[1] + J[2] * J[2] + J[3] * J[3]) +
            THETA * (D * D + 1.0);
        E += w * F / ch;
        if (g)
        {
            /* dE/dJ = a (2 (1 - THETA) J + 2 THETA D cof J) - b cof J,
               with dchi/dD = chi / r; dE/dP = dE/dJ R^-T */
            double a = w / ch, b = w * F / (ch * r);
            double c1 = 2.0 * (1.0 - THETA) * a, c2 = 2.0 * THETA * a * D - b;
            double G00 = c1 * J[0] + c2 * J[3], G01 = c1 * J[1] - c2 * J[2];
            double G10 = c1 * J[2] - c2 * J[1], G11 = c1 * J[3] + c2 * J[0];
            double dP00 = G00 * R[0] + G01 * R[1], dP01 = G00 * R[2] + G01 * R[3];
            double dP10 = G10 * R[0] + G11 * R[1], dP11 = G10 * R[2] + G11 * R[3];

            add_grad(p, g, p->lt[3 * i + 1], dP00, dP10);
            add_grad(p, g, p->lt[3 * i + 2], dP01, dP11);
            add_grad(p, g, p->lt[3 * i], -dP00 - dP01, -dP10 - dP11);
        }
    }
    return E;
}

static int lbfgs_reserve(lbfgs_buf *b, int n)
{
    if (n <= b->cap)
        return 0;
    free(b->g);
    b->g = (double *)malloc(sizeof(double) * (size_t)n * (4 + 2 * LBFGS_M));
    if (!b->g)
    {
        b->cap = 0;
        return -1;
    }
    b->d = b->g + n;
    b->zn = b->d + n;
    b->gn = b->zn + n;
    b->s = b->gn + n;
    b->y = b->s + (size_t)n * LBFGS_M;
    b->cap = n;
    return 0;
}

/* L-BFGS with a backtracking (Armijo) line search; a step along the plain
   gradient moves no coordinate further than step0.  Returns the energy at
   the minimum found, z holds its coordinates. */
static double minimise(const patch *p, double *z, double eps, double step0,
                       lbfgs_buf *b)
{
    int n = 2 * p->nS, it, i, k, m = 0, newest = -1;
    double E, rho[LBFGS_M], alpha[LBFGS_M];

    E = patch_energy(p, z, eps, b->g);
    if (!isfinite(E))
        return E;
    for (it = 0; it < LBFGS_ITER; it++)
    {
        double gd, t = 1.0, En = HUGE_VAL, sy, gmax = 0.0;
        int ls;

        /* two-loop recursion: d = H g */
        memcpy(b->d, b->g, sizeof(double) * n);
        for (k = 0; k < m; k++)
        {
            int j = (newest - k + LBFGS_M) % LBFGS_M;
            alpha[j] = rho[j] * dot(b->s + (size_t)j * n, b->d, n);
            for (i = 0; i < n; i++)
                b->d[i] -= alpha[j] * b->y[(size_t)j * n + i];
        }
        if (m > 0)
        {
            const double *yj = b->y + (size_t)newest * n;
            double gam = dot(b->s + (size_t)newest * n, yj, n) / dot(yj, yj, n);
            for (i = 0; i < n; i++)
                b->d[i] *= gam;
        }
        for (k = m - 1; k >= 0; k--)
        {
            int j = (newest - k + LBFGS_M) % LBFGS_M;
            double beta = rho[j] * dot(b->y + (size_t)j * n, b->d, n);
            for (i = 0; i < n; i++)
                b->d[i] += (alpha[j] - beta) * b->s[(size_t)j * n + i];
        }
        gd = dot(b->g, b->d, n);
        if (m == 0 || !(gd > 0.0))
        {
            /* no curvature yet, or not a descent direction: the gradient */
            m = 0;
            for (i = 0; i < n; i++)
                if (fabs(b->g[i]) > gmax)
                    gmax = fabs(b->g[i]);
            if (gmax <= 0.0)
                break;
            for (i = 0; i < n; i++)
                b->d[i] = b->g[i] * step0 / gmax;
            gd = dot(b->g, b->d, n);
        }

        /* d is the ascent direction: step along -d */
        for (ls = 0; ls < 40; ls++)
        {
            for (i = 0; i < n; i++)
                b->zn[i] = z[i] - t * b->d[i];
            En = patch_energy(p, b->zn, eps, b->gn);
            if (isfinite(En) && En <= E - 1e-4 * t * gd)
                break;
            t *= 0.5;
        }
        if (ls == 40)
            break;

        newest = (newest + 1) % LBFGS_M;
        for (i = 0; i < n; i++)
        {
            b->s[(size_t)newest * n + i] = b->zn[i] - z[i];
            b->y[(size_t)newest * n + i] = b->gn[i] - b->g[i];
        }
        sy = dot(b->s + (size_t)newest * n, b->y + (size_t)newest * n, n);
        if (sy > 0.0)
        {
            rho[newest] = 1.0 / sy;
            if (m < LBFGS_M)
                m++;
        }
        else
        {
            /* no positive curvature along the step: forget the history */
            newest = (newest - 1 + LBFGS_M) % LBFGS_M;
            m = 0;
        }
        memcpy(z, b->zn, sizeof(double) * n);
        memcpy(b->g, b->gn, sizeof(double) * n);
        if (E - En <= 1e-10 * fabs(E))
        {
            E = En;
            break;
        }
        E = En;
    }
    return E;
}

/* Lifts the free vertices from the plane onto the sphere and keeps them if
   every triangle of the patch is good; restores them otherwise. */
static int apply_patch(untangle_ctx *c, patch *p)
{
    int i, k;

    for (i = 0; i < p->nS; i++)
    {
        double *x = c->x + 3 * p->S[i], q[3], nrm;

        for (k = 0; k < 3; k++)
        {
            p->saved[3 * i + k] = x[k];
            q[k] = p->c[k] + (p->z[2 * i] * p->e1[k] + p->z[2 * i + 1] * p->e2[k]) / c->radius;
        }
        nrm = sqrt(q[0] * q[0] + q[1] * q[1] + q[2] * q[2]);
        for (k = 0; k < 3; k++)
            x[k] = c->radius * q[k] / nrm;
    }
    for (i = 0; i < p->nF; i++)
        if (signed_height(c, c->x, p->F[i]) < c->hmin)
            break;
    if (i == p->nF)
        return 1;
    for (i = 0; i < p->nS; i++)
        for (k = 0; k < 3; k++)
            c->x[3 * p->S[i] + k] = p->saved[3 * i + k];
    return 0;
}

/* Minimises the energy of the patch while eps is lowered as in Garanzha et
   al. (2021): after each minimisation eps is set so that chi of the worst
   triangle drops by the relative decrease of the energy (at least 10%), until
   no triangle is folded.  Returns 1 and applies the patch if it ends up good. */
static int untangle_patch(untangle_ctx *c, patch *p, lbfgs_buf *b, int *stages)
{
    double mindet = min_det(p, p->z), eps, E0, E1, sigma, mu, step0 = 0.0;
    int stage, i, n_neg, best_neg = p->nF + 1, since_best = 0;

    for (i = 0; i < p->nF; i++)
        step0 += p->aref[i];
    step0 = 0.1 * sqrt(step0 / p->nF);
    /* the start rule of the paper assumes determinants of order 1; a start
       without folds, whose smallest determinant can be far below that where
       the map is compressed, starts with a barrier that is already hard */
    eps = mindet > 0.0 ? 0.1 * mindet : sqrt(1e-12 + 0.04 * mindet * mindet);
    for (stage = 0; stage < MAX_STAGES; stage++)
    {
        E0 = patch_energy(p, p->z, eps, NULL);
        E1 = minimise(p, p->z, eps, step0, b);
        if (!isfinite(E0) || !isfinite(E1))
            break;
        mindet = min_det(p, p->z);
        /* a patch whose folds stop going away has no fold-free embedding */
        n_neg = count_neg(p, p->z);
        if (n_neg < best_neg)
        {
            best_neg = n_neg;
            since_best = 0;
        }
        else if (n_neg > 0 && ++since_best >= ABORT_STAGES)
            break;
        sigma = 1.0 - E1 / E0;
        if (sigma < 0.1)
            sigma = 0.1;
        mu = (1.0 - sigma) * chi(mindet, eps);
        if (mindet < mu)
            eps = 2.0 * sqrt(mu * (mu - mindet));
        else
            eps = 0.1 * mindet < 1e-10 ? 0.1 * mindet : 1e-10;
        if (mindet > 0.0 && E0 - E1 <= 1e-3 * E1 && apply_patch(c, p))
        {
            *stages = stage + 1;
            return 1;
        }
    }
    *stages = stage;
    return mindet > 0.0 && apply_patch(c, p);
}

static void free_ctx(untangle_ctx *c, patch *p, lbfgs_buf *b)
{
    free(c->f);
    free(c->x);
    free(c->rinv);
    free(c->aref);
    if (c->n_nbr)
        delete_all_polygon_point_neighbours(c->n_nbr, c->nbr);
    free(c->vf_off);
    free(c->vf);
    free(c->ff);
    free(c->vmark);
    free(c->bmark);
    free(c->fmark);
    free(c->fseen);
    free(c->nmark);
    free(c->next);
    free(c->loc);
    free(c->moved);
    free(p->S);
    free(p->B);
    free(p->F);
    free(p->comp);
    free(p->L);
    free(p->lt);
    free(p->rinv);
    free(p->aref);
    free(p->pfix);
    free(p->z);
    free(p->saved);
    free(b->g);
}

/**
 * \brief Remove the folded (negative-area) triangles of a spherical map.
 *
 * Bad triangles (signed height below 1e-6 of the radius) that share a vertex
 * form a cluster; every cluster is re-embedded on a disk of the mesh around
 * it, grown by 2, 4, 8, 16 or 32 rings until it works, with the vertices just
 * outside the disk fixed.  In the gnomonic plane of the disk the barrier
 * energy of Garanzha et al. (2021) -- the shape and area of every triangle
 * against its shape on \p surface, divided by a regularised determinant -- is
 * minimised by L-BFGS while the regularisation is lowered until no triangle
 * is folded.  A disk is accepted only if all its triangles end up good.
 * A map without bad triangles is left untouched, bit for bit.
 *
 * \param surface (in)     the surface the sphere maps, same mesh as \p sphere
 * \param sphere  (in/out) spherical map of \p surface
 * \param verbose (in)     print the clusters and what was done with them
 * \return the number of bad triangles left (0: none), or -1 if the meshes
 *         differ, are not triangulated, or memory runs out
 */
int CAT_SphereRemoveFolds(polygons_struct *surface, polygons_struct *sphere,
                          int verbose)
{
    untangle_ctx c;
    patch p;
    lbfgs_buf b;
    double *v = NULL, centre[3], lo[3], hi[3], area_v = 0.0, area_x = 0.0;
    int *parent = NULL, *cid = NULL, *off = NULL, *seeds = NULL;
    int i, k, t, n_bad, n_bad0, n_cl, round, pos = 0, ret = -1;
    const char *why = "out of memory";

    memset(&c, 0, sizeof(c));
    memset(&p, 0, sizeof(p));
    memset(&b, 0, sizeof(b));
    if (surface->n_points != sphere->n_points || surface->n_items != sphere->n_items)
    {
        fprintf(stderr, "CAT_SphereRemoveFolds: surface and sphere differ\n");
        return -1;
    }
    for (i = 0; i < sphere->n_items; i++)
        if (GET_OBJECT_SIZE(*sphere, i) != 3)
        {
            fprintf(stderr, "CAT_SphereRemoveFolds: mesh is not triangulated\n");
            return -1;
        }

    c.n = sphere->n_points;
    c.n_faces = sphere->n_items;
    c.f = (int *)malloc(sizeof(int) * 3 * c.n_faces);
    c.x = (double *)malloc(sizeof(double) * 3 * c.n);
    c.vmark = (int *)calloc(c.n, sizeof(int));
    c.bmark = (int *)calloc(c.n, sizeof(int));
    c.fmark = (int *)calloc(c.n_faces, sizeof(int));
    c.fseen = (int *)calloc(c.n_faces, sizeof(int));
    c.nmark = (int *)calloc(c.n, sizeof(int));
    c.next = (int *)malloc(sizeof(int) * c.n);
    parent = (int *)malloc(sizeof(int) * c.n);
    cid = (int *)malloc(sizeof(int) * c.n);
    off = (int *)malloc(sizeof(int) * (c.n + 1));
    seeds = (int *)malloc(sizeof(int) * c.n);
    if (!c.f || !c.x || !c.vmark || !c.bmark || !c.fmark || !c.fseen ||
        !c.nmark || !c.next || !parent || !cid || !off || !seeds)
        goto done;

    /* centre of the bounding box, as the relaxation takes it */
    for (k = 0; k < 3; k++)
    {
        lo[k] = 1e30;
        hi[k] = -1e30;
    }
    for (i = 0; i < c.n; i++)
        for (k = 0; k < 3; k++)
        {
            double xk = Point_coord(sphere->points[i], k);
            lo[k] = xk < lo[k] ? xk : lo[k];
            hi[k] = xk > hi[k] ? xk : hi[k];
        }
    for (k = 0; k < 3; k++)
        centre[k] = 0.5 * (lo[k] + hi[k]);
    for (i = 0; i < c.n; i++)
    {
        for (k = 0; k < 3; k++)
            c.x[3 * i + k] = Point_coord(sphere->points[i], k) - centre[k];
        c.radius += sqrt(c.x[3 * i] * c.x[3 * i] + c.x[3 * i + 1] * c.x[3 * i + 1] +
                         c.x[3 * i + 2] * c.x[3 * i + 2]) / c.n;
    }
    c.hmin = H_MIN * c.radius;
    for (t = 0; t < c.n_faces; t++)
        for (k = 0; k < 3; k++)
            c.f[3 * t + k] = sphere->indices[POINT_INDEX(sphere->end_indices, t, k)];

    /* the orientation of most triangles is the right one */
    c.orient = 1;
    for (t = 0; t < c.n_faces; t++)
        pos += signed_height(&c, c.x, t) > 0.0;
    if (2 * pos < c.n_faces)
        c.orient = -1;

    n_cl = find_clusters(&c, parent, cid, off, seeds, &n_bad);
    n_bad0 = n_bad;
    if (verbose)
        printf("Fold removal: %d bad triangles in %d clusters\n", n_bad, n_cl);
    if (n_bad == 0)
    {
        ret = 0;
        goto done;
    }

    v = (double *)malloc(sizeof(double) * 3 * c.n);
    c.rinv = (double *)malloc(sizeof(double) * 4 * c.n_faces);
    c.aref = (double *)malloc(sizeof(double) * c.n_faces);
    c.loc = (int *)malloc(sizeof(int) * c.n);
    c.moved = (char *)calloc(c.n, 1);
    p.S = (int *)malloc(sizeof(int) * c.n);
    p.B = (int *)malloc(sizeof(int) * c.n);
    p.F = (int *)malloc(sizeof(int) * c.n_faces);
    p.comp = (int *)malloc(sizeof(int) * c.n_faces);
    p.L = (int *)malloc(sizeof(int) * c.n);
    p.lt = (int *)malloc(sizeof(int) * 3 * c.n_faces);
    p.rinv = (double *)malloc(sizeof(double) * 4 * c.n_faces);
    p.aref = (double *)malloc(sizeof(double) * c.n_faces);
    p.pfix = (double *)malloc(sizeof(double) * 2 * c.n);
    p.z = (double *)malloc(sizeof(double) * 2 * c.n);
    p.saved = (double *)malloc(sizeof(double) * 3 * c.n);
    if (!v || !c.rinv || !c.aref || !c.loc || !c.moved || !p.S || !p.B || !p.F ||
        !p.comp || !p.L || !p.lt || !p.rinv || !p.aref || !p.pfix || !p.z || !p.saved)
        goto done;

    for (i = 0; i < c.n; i++)
        for (k = 0; k < 3; k++)
            v[3 * i + k] = Point_coord(surface->points[i], k);
    for (t = 0; t < c.n_faces; t++)
    {
        const double *a = v + 3 * c.f[3 * t], *bb = v + 3 * c.f[3 * t + 1];
        const double *d = v + 3 * c.f[3 * t + 2];
        double u[3], w[3], cr[3];

        for (k = 0; k < 3; k++)
        {
            u[k] = bb[k] - a[k];
            w[k] = d[k] - a[k];
        }
        cr[0] = u[1] * w[2] - u[2] * w[1];
        cr[1] = u[2] * w[0] - u[0] * w[2];
        cr[2] = u[0] * w[1] - u[1] * w[0];
        area_v += 0.5 * sqrt(cr[0] * cr[0] + cr[1] * cr[1] + cr[2] * cr[2]);
    }
    area_x = 4.0 * PI * c.radius * c.radius;
    k = area_v > 0.0 ? build_reference(&c, v, sqrt(area_x / area_v)) : -2;
    if (k != 0)
    {
        if (k == -2)
            why = "the surface is degenerate";
        goto done;
    }
    if (build_topology(&c, sphere) != 0)
        goto done;

    for (round = 0; round < MAX_ROUNDS && n_cl > 0; round++)
    {
        int cleared = 0;

        for (k = 0; k < n_cl; k++)
        {
            int r, stages = 0, chi_patch = 0, ok = 0;
            double rim_area;

            for (r = 0; r < N_PATCH_RINGS; r++)
            {
                chi_patch = grow_patch(&c, seeds + off[k], off[k + 1] - off[k],
                                       patch_rings[r], &p);
                if (chi_patch < 0)
                    goto done;
                if (chi_patch != 1)
                    continue;
                for (i = 0; i < p.nF; i++)
                    if (signed_height(&c, c.x, p.F[i]) < c.hmin)
                        break;
                if (i == p.nF)
                {
                    /* an earlier patch has resolved this cluster */
                    ok = 1;
                    break;
                }
                if (rim_loop(&c, &p) != 0)
                {
                    if (verbose > 2)
                        printf("    %d rings: rim is not one loop\n", patch_rings[r]);
                    continue;
                }
                if (frame_patch(&c, &p) != 0)
                {
                    if (verbose > 2)
                        printf("    %d rings: patch too large\n", patch_rings[r]);
                    break;
                }
                rim_area = rim_simple(&c, &p);
                if (rim_area <= 0.0)
                {
                    if (verbose > 2)
                        printf("    %d rings: rim not simple\n", patch_rings[r]);
                    continue;
                }
                normalise_patch(&p, rim_area);
                if (lbfgs_reserve(&b, 2 * p.nS) != 0)
                    goto done;
                tutte_start(&c, &p, b.d, b.zn, b.gn);
                if (verbose > 2)
                    printf("    %d rings: %d free, %d triangles, %d folded after Tutte\n",
                           patch_rings[r], p.nS, p.nF, count_neg(&p, p.z));
                if (untangle_patch(&c, &p, &b, &stages))
                {
                    for (i = 0; i < p.nS; i++)
                        c.moved[p.S[i]] = 1;
                    ok = 1;
                    break;
                }
            }
            cleared += ok;
            if (verbose > 1)
                printf("  cluster %d (%d vertices): %s, %d rings, %d free vertices, "
                       "%d stages\n", k, off[k + 1] - off[k], ok ? "removed" : "left",
                       patch_rings[r < N_PATCH_RINGS ? r : N_PATCH_RINGS - 1],
                       p.nS, stages);
        }
        n_cl = find_clusters(&c, parent, cid, off, seeds, &n_bad);
        if (verbose)
            printf("  round %d: %d clusters removed, %d bad triangles left\n",
                   round + 1, cleared, n_bad);
        if (cleared == 0)
            break;
    }

    for (i = 0; i < c.n; i++)
        if (c.moved[i])
            fill_Point(sphere->points[i], c.x[3 * i] + centre[0],
                       c.x[3 * i + 1] + centre[1], c.x[3 * i + 2] + centre[2]);
    if (verbose)
        printf("Fold removal: %d of %d bad triangles removed\n", n_bad0 - n_bad, n_bad0);
    ret = n_bad;

done:
    if (ret < 0)
        fprintf(stderr, "CAT_SphereRemoveFolds: %s\n", why);
    free(v);
    free(parent);
    free(cid);
    free(off);
    free(seeds);
    free_ctx(&c, &p, &b);
    return ret;
}
