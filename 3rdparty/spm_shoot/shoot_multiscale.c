/*
 * John Ashburner
 * Copyright (C) 2011-2022 Wellcome Centre for Human Neuroimaging
 * Modified for CAT-Surface: no mex, spm_ symbol prefix, pthreads (see README.md).
 */

#include "shoot_namespace.h"

#include <math.h>
#include <stdlib.h>
#include "shoot_types.h"
#include "shoot_boundary.h"
#include "shoot_parallel.h"

/* 2nd degree B-spline basis */
static double wt2(double x)
{
    x = fabs(x);
    if (x < 0.5)
        return(0.75 - x*x);
    if (x < 1.5)
    {
        x = 1.5 - x;
        return(0.5*x*x);
    }
    return(0.0);
}

/* Linear interpolation (1st degree B-spline) basis */
static double wt1(double x)
{
    x = fabs(x);
    if (x < 1.0)
        return(1.0-x);
    return(0.0);
}

/* Note that restriction uses linear interpolation, whereas prolongation
   uses 2nd degree B-spline. For "bending energy", the restriction and
   prolongation operators need to be at least this degree.  See Numerical
   Recipes for further details. */
static void restrict_plane(mwSize na[], float *a,  mwSize nc[], float *c, float *b)
{
    mwSignedIndex i, j, o;
    double loc, s;
    float *ap, *bp, *cp;
    /* a - na[0]*na[1]
     * c - nc[0]*nc[1]
     * b - na[0]*nc[1]
     */
    if (na[1]==1)
    {
        for(j=0; j<(mwSignedIndex)nc[1]*na[0]; j++)
            b[j] = a[j];
    }
    /*
    else if (!(na[1]%2))
    {
        for(j=0; j<nc[1]; j++)
        {
            mwSignedIndex oo,om,op;
            o   = 2*j;
            oo  = o*na[0];
            om  = bound(o-1,na[1])*na[0];
            op  = bound(o+1,na[1])*na[0];
            for(ap=a, bp=b+j*na[0], cp=ap+na[0]; ap<cp; ap++, bp++)
                *bp = 0.25*(ap[om]+ap[op])+0.5*ap[oo];
        }
    }
    */
    else
    {
        s = (double)na[1]/(double)nc[1];
        for(j=0; j<(mwSignedIndex)nc[1]; j++)
        {
            double w0, w1, w2, w3;
            mwSignedIndex o0, o1, o2, o3;
            /* loc = s*j; */
            loc = ((double)j+0.5)*s-0.5;
            o   = (mwSignedIndex)floor(loc);
            o0  = bound(o-1,na[1])*na[0];
            o1  = bound(o  ,na[1])*na[0];
            o2  = bound(o+1,na[1])*na[0];
            o3  = bound(o+2,na[1])*na[0];
            w0  = wt1(((double)(o-1)-loc)/2.0)/2.0;
            w1  = wt1(((double)(o  )-loc)/2.0)/2.0;
            w2  = wt1(((double)(o+1)-loc)/2.0)/2.0;
            w3  = wt1(((double)(o+2)-loc)/2.0)/2.0;
            for(ap=a, bp=b+j*na[0], cp=ap+na[0]; ap<cp; ap++, bp++)
                *bp = (float)(w0*ap[o0]+w1*ap[o1]+w2*ap[o2]+w3*ap[o3]);
        }
    }

    if (na[0]==1)
    {
        for(j=0; j<(mwSignedIndex)nc[0]*nc[1]; j++)
            c[j] = b[j];
    }
    /*
    else if (!(na[0]%2))
    {
        for(i=0; i<nc[0]; i++)
        {
            mwSignedIndex om,op;
            o   = 2*i;
            om  = bound(o-1,na[0]);
            op  = bound(o+1,na[0]);
            for(bp=b, cp=c+i, ap=bp+na[0]*nc[1]; bp<ap; bp+=na[0], cp+=nc[0])
                *cp = 0.25*(bp[om]+bp[op])+0.5*bp[o];
        }
    }
    */
    else
    {
        s = (double)na[0]/(double)nc[0];
        for(i=0; i<(mwSignedIndex)nc[0]; i++)
        {
            double w0, w1, w2, w3;
            mwSignedIndex o0, o1, o2, o3;
            /* loc = s*i; */
            loc = ((double)i+0.5)*s-0.5;
            o   = (mwSignedIndex)floor(loc);
            o0  = bound(o-1,na[0]); 
            o1  = bound(o  ,na[0]); 
            o2  = bound(o+1,na[0]); 
            o3  = bound(o+2,na[0]); 
            w0  = wt1(((double)(o-1)-loc)/2.0)/2.0;
            w1  = wt1(((double)(o  )-loc)/2.0)/2.0;
            w2  = wt1(((double)(o+1)-loc)/2.0)/2.0;
            w3  = wt1(((double)(o+2)-loc)/2.0)/2.0;
            for(bp=b, cp=c+i, ap=bp+na[0]*nc[1]; bp<ap; bp+=na[0], cp+=nc[0])
                *cp = (float)(w0*bp[o0]+w1*bp[o1]+w2*bp[o2]+w3*bp[o3]);
        }
    }
}

/* Output slices k0..k1-1; the planes in b are filled at the first one */
static void restrict_vol_k(mwSize na[], float *a, mwSize nc[], float *c, float *b, mwSignedIndex k0, mwSignedIndex k1)
{
    mwSignedIndex j, k, o=-999999, oo;
    mwSize m;
    double loc, s;
    float *bp, *cp, *pl[4];

    if (na[2]==1)
    {
        restrict_plane(na,a,nc,c,b);
        return;
    }
    m     = nc[0]*nc[1];
    pl[0] = b;
    pl[1] = b + m;
    pl[2] = b + m*2;
    pl[3] = b + m*3;
    bp    = b + m*4;

    /* This would work if loc = s*k and would be faster
    if (!(na[2]%2))
    {
        restrict_plane(na, a+na[0]*na[1]*2,nc,pl[2],bp);
        for(k=0; k<nc[2]; k++)
        {
            mwSignedIndex om, op;
            float *tp;
            o   = 2*k;
            om  = bound(o-1,na[2]);
            op  = bound(o+1,na[2]);

            tp    = pl[0];
            pl[0] = pl[2];
            pl[2] = tp;
            restrict_plane(na, a+na[0]*na[1]*o ,nc,pl[1],bp);
            restrict_plane(na, a+na[0]*na[1]*op,nc,pl[2],bp);

            cp  = c+nc[0]*nc[1]*k;
            for(j=0; j<nc[0]*nc[1]; j++)
            {
                cp[j] = 0.25*(pl[0][j]+pl[2][j])+0.5*pl[1][j];
            }
        }
    }
    else
    */
    {
        s      = (double)na[2]/(double)nc[2];
        for(k=k0; k<k1; k++)
        {
            double w0, w1, w2, w3;
            mwSignedIndex o0, o1, o2, o3;
            /* loc = s*k; */
            loc = ((double)k+0.5)*s-0.5;
            oo  = o;
            o   = (mwSignedIndex)floor(loc);

            o0  = bound(o-1,na[2]);
            o1  = bound(o  ,na[2]);
            o2  = bound(o+1,na[2]);
            o3  = bound(o+2,na[2]);
            w0  = wt1(((double)(o-1)-loc)/2.0)/2.0;
            w1  = wt1(((double)(o  )-loc)/2.0)/2.0;
            w2  = wt1(((double)(o+1)-loc)/2.0)/2.0;
            w3  = wt1(((double)(o+2)-loc)/2.0)/2.0;

            if (o==oo)
            {   /* do nothing */
            }
            else if (o==oo+1)
            {   /* Shift by 1 */
                float *tp;
                tp    = pl[0];
                pl[0] = pl[1];
                pl[1] = pl[2];
                pl[2] = pl[3];
                pl[3] = tp;
                restrict_plane(na, a+na[0]*na[1]*o3,nc,pl[3],bp);
            }
            else if (o==oo+2)
            {   /* Shift by 2 */
                float *tp;
                tp    = pl[0];
                pl[0] = pl[2];
                pl[2] = tp;
                tp    = pl[1];
                pl[1] = pl[3];
                pl[3] = tp;
                restrict_plane(na, a+na[0]*na[1]*o2,nc,pl[2],bp);
                restrict_plane(na, a+na[0]*na[1]*o3,nc,pl[3],bp);
            }
            else if (o==oo+3)
            {   /* Shift by 2 */
                float *tp;
                tp    = pl[0];
                pl[0] = pl[3];
                pl[3] = tp;
                restrict_plane(na, a+na[0]*na[1]*o1,nc,pl[1],bp);
                restrict_plane(na, a+na[0]*na[1]*o2,nc,pl[2],bp);
                restrict_plane(na, a+na[0]*na[1]*o3,nc,pl[3],bp);
            }
            else
            {   /* Read everything */
                restrict_plane(na, a+na[0]*na[1]*o0,nc,pl[0],bp);
                restrict_plane(na, a+na[0]*na[1]*o1,nc,pl[1],bp);
                restrict_plane(na, a+na[0]*na[1]*o2,nc,pl[2],bp);
                restrict_plane(na, a+na[0]*na[1]*o3,nc,pl[3],bp);
            }
            cp  = c+nc[0]*nc[1]*k;
            for(j=0; j<(mwSignedIndex)nc[0]*nc[1]; j++)
            {
                cp[j] = (float)(w0*pl[0][j]+w1*pl[1][j]+w2*pl[2][j]+w3*pl[3][j]);
            }
        }
    }
}

static void resized_plane(mwSize na[], float *a,  mwSize nc[], float *c, float *b)
{
    mwSignedIndex i, j, o,oc,om,op;
    double loc, s, w, wm, wp;
    float *ap, *bp, *cp;
    /* a - na[0]*na[1]
     * c - nc[0]*nc[1]
     * b - na[0]*nc[1]
     */

    s = (double)na[1]/(double)nc[1];
    for(j=0; j<(mwSignedIndex)nc[1]; j++)
    {
        /* loc = j*s; */
        loc = ((double)j+0.5)*s-0.5;
        o   = (mwSignedIndex)floor(loc+0.5);
        oc  = bound(o  ,na[1])*na[0];
        om  = bound(o-1,na[1])*na[0];
        op  = bound(o+1,na[1])*na[0];
        w   = wt2((double) o   -loc);
        wp  = wt2((double)(o+1)-loc);
        wm  = wt2((double)(o-1)-loc);
        for(ap=a, bp=b+j*na[0], cp=ap+na[0]; ap<cp; ap++, bp++)
            *bp = (float)(wm*ap[om]+w*ap[oc]+wp*ap[op]);
    }
    s = (double)na[0]/(double)nc[0];
    for(i=0; i<(mwSignedIndex)nc[0]; i++)
    {
        /* loc = i*s; */
        loc = ((double)i+0.5)*s-0.5;
        o   = (mwSignedIndex)floor(loc+0.5);
        oc  = bound(o  ,na[0]);
        om  = bound(o-1,na[0]);
        op  = bound(o+1,na[0]);
        w   = wt2((double) o   -loc);
        wp  = wt2((double)(o+1)-loc);
        wm  = wt2((double)(o-1)-loc);
        for(bp=b, cp=c+i, ap=bp+na[0]*nc[1]; bp<ap; bp+=na[0], cp+=nc[0])
            *cp = (float)(wm*bp[om]+w*bp[oc]+wp*bp[op]);
    }
}

/* Output slices k0..k1-1; the planes in b are filled at the first one */
static void resize_vol_k(mwSize na[], float *a, mwSize nc[], float *c, float *b, mwSignedIndex k0, mwSignedIndex k1)
{
    mwSignedIndex j, k, o=-999999,oc,om,op, oo;
    mwSize m;
    double loc, s, w, wm, wp;
    float *bp, *cp, *pl[3];

    if (na[2]==1 && nc[2]==1)
    {
        resized_plane(na,a,nc,c,b);
        return;
    }

    m     = nc[0]*nc[1];
    pl[0] = b;
    pl[1] = b + m;
    pl[2] = b + m*2;
    bp    = b + m*3;

    for(k=k0; k<k1; k++)
    {
        s   = (double)na[2]/(double)nc[2];
        /* loc = k*s; */
        loc = ((double)k+0.5)*s-0.5;
        oo  = o;
        o   = (mwSignedIndex)floor(loc+0.5);
        oc  = bound(o  ,na[2]);
        om  = bound(o-1,na[2]);
        op  = bound(o+1,na[2]);

        if (o==oo)
        {   /* do nothing */
        }
        else if (o==oo+1)
        {   /* Shift by 1 */
            float *tp;
            tp    = pl[0];
            pl[0] = pl[1];
            pl[1] = pl[2];
            pl[2] = tp;
            resized_plane(na, a+na[0]*na[1]*op,nc,pl[2],bp);
        }
        else if (o==oo+2)
        {   /* Shift by 2 */
            float *tp;
            tp    = pl[0];
            pl[0] = pl[2];
            pl[2] = tp;
            resized_plane(na, a+na[0]*na[1]*oc,nc,pl[1],bp);
            resized_plane(na, a+na[0]*na[1]*op,nc,pl[2],bp);
        }
        else
        {   /* Read everything */
            resized_plane(na, a+na[0]*na[1]*om,nc,pl[0],bp);
            resized_plane(na, a+na[0]*na[1]*oc,nc,pl[1],bp);
            resized_plane(na, a+na[0]*na[1]*op,nc,pl[2],bp);
        }
        w   = wt2((double) o   -loc);
        wp  = wt2((double)(o+1)-loc);
        wm  = wt2((double)(o-1)-loc);
        cp  = c+nc[0]*nc[1]*k;
        for(j=0; j<(mwSignedIndex)nc[0]*nc[1]; j++)
        {
            cp[j] = (float)(wm*pl[0][j]+w*pl[1][j]+wp*pl[2][j]);
        }
    }
}

/* Threads: each takes a block of output slices with its own planes, which
   hold the same values as in the serial loop, so the result is the same. */
typedef struct
{
    mwSize        *na, *nc;
    float         *a, *c, *b;
    mwSignedIndex  k0, k1;
    int            up;
} transfer_slot;

static void *transfer_worker(void *arg)
{
    transfer_slot *s = (transfer_slot *)arg;
    if (s->up)
        resize_vol_k(s->na, s->a, s->nc, s->c, s->b, s->k0, s->k1);
    else
        restrict_vol_k(s->na, s->a, s->nc, s->c, s->b, s->k0, s->k1);
    return NULL;
}

static void transfer(mwSize na[], float *a, mwSize nc[], float *c, float *b, int up)
{
    transfer_slot slot[64];
    float *buf = 0;
    mwSize bsz = 4*nc[0]*nc[1] + na[0]*nc[1]; /* the size of b */
    int t, nt = 0;

    if (nc[2]>1 && na[2]>1)
        nt = spm_shoot_nworkers(nc[2], (up ? 0.8 : 1.2)*(double)(na[0]*na[1]));
    if (nt>64) nt = 64;
    if (nt>1) buf = (float *)malloc(bsz*(mwSize)nt*sizeof(float));
    if (buf==0)
    {
        if (up)
            resize_vol_k(na, a, nc, c, b, 0, (mwSignedIndex)nc[2]);
        else
            restrict_vol_k(na, a, nc, c, b, 0, (mwSignedIndex)nc[2]);
        return;
    }
    for(t=0; t<nt; t++)
    {
        slot[t].na = na;
        slot[t].nc = nc;
        slot[t].a  = a;
        slot[t].c  = c;
        slot[t].b  = buf + bsz*(mwSize)t;
        slot[t].k0 = (mwSignedIndex)(nc[2]*(mwSize)t/(mwSize)nt);
        slot[t].k1 = (mwSignedIndex)(nc[2]*(mwSize)(t+1)/(mwSize)nt);
        slot[t].up = up;
    }
    cat_parallel_run(nt, transfer_worker, slot, sizeof(transfer_slot));
    free(buf);
}

void restrict_vol(mwSize na[], float *a, mwSize nc[], float *c, float *b)
{
    transfer(na, a, nc, c, b, 0);
}

void resize_vol(mwSize na[], float *a, mwSize nc[], float *c, float *b)
{
    transfer(na, a, nc, c, b, 1);
}
