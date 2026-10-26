/*
 * John Ashburner
 * Copyright (C) 2001-2022 Wellcome Centre for Human Neuroimaging
 * Modified for CAT-Surface: spm_bsplinc.c and spm_bsplins.c for volumes in
 * memory instead of the mex gateways, pthreads (see README.md).
 */

#include "shoot_namespace.h"
#include <math.h>
#include <stdlib.h>
#include "shoot_types.h"
#include "bsplines.h"
#include "shoot_parallel.h"
#include "spm_shoot.h"

/***************************************************************************************
Deconvolve the B-spline basis functions from the image volume (vol_coeffs of
spm_bsplinc.c).  Every 1D line is filtered on its own, so the lines are
split among the threads.
*/
typedef struct
{
    const mwSize *dim;
    double       *c;
    void        (*splinc)(IMAGE_DTYPE c[], int m, double p[], int np);
    double        p[4];
    int           np;
} splinc_args;

/* Deconvolve along the fastest dimension (X), slices k=lo..hi-1 */
static void splinc_x(void *ctx, size_t lo, size_t hi)
{
    splinc_args *a = (splinc_args *)ctx;
    const mwSize *dim = a->dim;
    double *cp;
    mwSize j, k;

    for(k=lo; k<hi; k++)
    {
        for(j=0; j<dim[1]; j++)
        {
            cp = &a->c[dim[0]*(j+dim[1]*k)];
            a->splinc(cp, (int)dim[0], a->p, a->np);
        }
    }
}

/* Deconvolve along the middle dimension (Y), slices k=lo..hi-1 */
static void splinc_y(void *ctx, size_t lo, size_t hi)
{
    splinc_args *a = (splinc_args *)ctx;
    const mwSize *dim = a->dim;
    double *cp, f[10240];
    mwSize i, j, k, n = dim[0];

    for(k=lo; k<hi; k++)
    {
        for(i=0;i<dim[0];i++)
        {
            cp = &a->c[i+dim[0]*dim[1]*k];
            for(j=0; j<dim[1]; j++, cp+=n)
                f[j] = *cp;
            a->splinc(f, (int)dim[1], a->p, a->np);
            cp = &a->c[i+dim[0]*dim[1]*k];
            for(j=0; j<dim[1]; j++, cp+=n)
                *cp = f[j];
        }
    }
}

/* Deconvolve along the slowest dimension (Z), rows j=lo..hi-1 */
static void splinc_z(void *ctx, size_t lo, size_t hi)
{
    splinc_args *a = (splinc_args *)ctx;
    const mwSize *dim = a->dim;
    double *cp, f[10240];
    mwSize i, j, k, n = dim[0]*dim[1];

    for(j=lo; j<hi; j++)
    {
        for(i=0;i<dim[0];i++)
        {
            cp = &a->c[i+dim[0]*j];
            for(k=0; k<dim[2]; k++, cp+=n)
                f[k] = *cp;
            a->splinc(f, (int)dim[2], a->p, a->np);
            cp = &a->c[i+dim[0]*j];
            for(k=0; k<dim[2]; k++, cp+=n)
                *cp = f[k];
        }
    }
}

int spm_bsplinc3(const float *vol, const mwSize dim[3], const int deg[3],
                 const int wrapflags[3], double *coef)
{
    splinc_args a;
    mwSize i, n;
    int k;

    for(k=0; k<3; k++)
        if (deg[k]<0 || deg[k]>7)
        {
            mexErrMsgTxt("Bad spline degree.");
            return(-1);
        }

    /* Check that dimensions don't exceed size of f */
    if (dim[1]>10240 || dim[2]>10240)
    {
        mexErrMsgTxt("Problem with deconvolution.");
        return(1);
    }

    /* Do a straight copy: the mex reads the volume through resample()
       (nearest neighbour), which returns 1.0*value + 0.0 */
    n = dim[0]*dim[1]*dim[2];
    for(i=0; i<n; i++)
    {
        coef[i] = 1.0*(double)vol[i] + 0.0;

        /* Not sure how best to handle NaNs */
        if (!mxIsFinite(coef[i])) coef[i] = 0.0;
    }

    a.dim = dim;
    a.c   = coef;

    /* Deconvolve along the fastest dimension (X) */
    if (deg[0]>1 && dim[0]>1)
    {
        if (get_poles(deg[0], &a.np, a.p))
        {
            mexErrMsgTxt("Problem with deconvolution.");
            return(1);
        }
        a.splinc = wrapflags[0] ? splinc_wrap : splinc_mirror;
        spm_shoot_parallel(dim[2], 20.0*(double)(dim[0]*dim[1]), splinc_x, &a);
    }

    /* Deconvolve along the middle dimension (Y) */
    if (deg[1]>1 && dim[1]>1)
    {
        if (get_poles(deg[1], &a.np, a.p))
        {
            mexErrMsgTxt("Problem with deconvolution.");
            return(1);
        }
        a.splinc = wrapflags[1] ? splinc_wrap : splinc_mirror;
        spm_shoot_parallel(dim[2], 25.0*(double)(dim[0]*dim[1]), splinc_y, &a);
    }

    /* Deconvolve along the slowest dimension (Z) */
    if (deg[2]>1 && dim[2]>1)
    {
        if (get_poles(deg[2], &a.np, a.p))
        {
            mexErrMsgTxt("Problem with deconvolution.");
            return(1);
        }
        a.splinc = wrapflags[2] ? splinc_wrap : splinc_mirror;
        spm_shoot_parallel(dim[1], 25.0*(double)(dim[0]*dim[2]), splinc_z, &a);
    }
    return(0);
}

/***************************************************************************************
Resample points and, optionally, their derivatives (fun and dfun of
spm_bsplins.c).  The points are split among the threads.
*/
#define TINY 5e-2

/* Coordinates beyond this make the float-to-int conversion in weights()
   undefined; the arm64 mex then returns 0 (all weights vanish) */
#define BIG 1073741824.0

typedef struct
{
    const double *c;
    int           m0, m1, m2;
    const double *x0, *x1, *x2;
    int           d[3];
    int           cond;
    int         (*bnd[3])();
    double       *f, *df0, *df1, *df2;
} bsplins_args;

static void bsplins_range(void *ctx, size_t lo, size_t hi)
{
    bsplins_args *a = (bsplins_args *)ctx;
    double NaN = mxGetNaN();
    int deriv = (a->df0 != 0 || a->df1 != 0 || a->df2 != 0);
    int cond = a->cond, m0 = a->m0, m1 = a->m1, m2 = a->m2;
    size_t j;

    for(j=lo; j<hi; j++)
    {
        double x0 = a->x0[j], x1 = a->x1[j], x2 = a->x2[j];
        double f = NaN, g0 = 0.0, g1 = 0.0, g2 = 0.0;

        if (((cond&1) | (x0>=1-TINY && x0<=m0+TINY)) &&
            ((cond&2) | (x1>=1-TINY && x1<=m1+TINY)) &&
            ((cond&4) | (x2>=1-TINY && x2<=m2+TINY)))
        {
            if (((cond&1) && !(fabs(x0)<BIG)) ||
                ((cond&2) && !(fabs(x1)<BIG)) ||
                ((cond&4) && !(fabs(x2)<BIG)))
                f = 0.0;
            else if (deriv)
                f = dsample3((double *)a->c, m0,m1,m2, x0-1,x1-1,x2-1, (int *)a->d,
                             &g0,&g1,&g2, a->bnd);
            else
                f = sample3((double *)a->c, m0,m1,m2, x0-1,x1-1,x2-1, (int *)a->d, a->bnd);
        }
        a->f[j] = f;
        if (a->df0) a->df0[j] = g0;
        if (a->df1) a->df1[j] = g1;
        if (a->df2) a->df2[j] = g2;
    }
}

void spm_bsplins3(const double *coef, const mwSize dim[3], mwSize n,
                  const double *x0, const double *x1, const double *x2,
                  const int deg[3], const int wrapflags[3],
                  double *f, double *d0, double *d1, double *d2)
{
    bsplins_args a;
    int k;

    a.c    = coef;
    a.m0   = (int)dim[0];
    a.m1   = (int)dim[1];
    a.m2   = (int)dim[2];
    a.x0   = x0;
    a.x1   = x1;
    a.x2   = x2;
    a.cond = 0;
    for(k=0; k<3; k++)
    {
        a.d[k]   = deg[k];
        a.bnd[k] = mirror;
        if (wrapflags[k])
        {
            a.bnd[k] = wrap;
            a.cond  += 1<<k;
        }
    }
    a.f   = f;
    a.df0 = d0;
    a.df1 = d1;
    a.df2 = d2;
    spm_shoot_parallel(n, (d0!=0 || d1!=0 || d2!=0) ? 120.0 : 60.0, bsplins_range, &a);
}
