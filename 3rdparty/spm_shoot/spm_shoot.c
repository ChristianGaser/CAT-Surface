/*
 * Error reporting and spm_diffeo gateway equivalents for the CAT-Surface
 * port of the SPM26 shooting primitives.  See spm_shoot.h.
 *
 * spm_shoot_fmg_param, spm_shoot_fmg and spm_shoot_mat reproduce what
 * spm_diffeo.c (John Ashburner, Copyright (C) 2011-2022 Wellcome Centre for
 * Human Neuroimaging) does around the computational code.
 *
 * Copyright (C) 2026 Christian Gaser, University of Jena.
 * Licence: GPL-2.0-or-later, as the SPM code it accompanies.
 */

#include "shoot_namespace.h"
#include <stdlib.h>
#include <string.h>
#include "shoot_types.h"
#include "shoot_boundary.h"
#include "shoot_regularisers.h"
#include "shoot_optim3d.h"
#include "shoot_multiscale.h"
#include "shoot_diffeo3d.h"
#include "shoot_expm3.h"
#include "shoot_parallel.h"
/* Last, so that its prototypes are checked against the ones above */
#include "spm_shoot.h"

#if !defined(_WIN32) && !defined(_WIN64)
#include <pthread.h>
/* invdef's threads can raise errors */
static pthread_mutex_t spm_shoot_err_mtx = PTHREAD_MUTEX_INITIALIZER;
#define ERR_LOCK   pthread_mutex_lock(&spm_shoot_err_mtx)
#define ERR_UNLOCK pthread_mutex_unlock(&spm_shoot_err_mtx)
#else
#define ERR_LOCK
#define ERR_UNLOCK
#endif

static const char *spm_shoot_err = 0;

void spm_shoot_error(const char *msg)
{
    ERR_LOCK;
    if (spm_shoot_err == 0) spm_shoot_err = msg;
    ERR_UNLOCK;
}

const char *spm_shoot_get_error(void)
{
    const char *msg;
    ERR_LOCK;
    msg = spm_shoot_err;
    ERR_UNLOCK;
    return msg;
}

void spm_shoot_clear_error(void)
{
    ERR_LOCK;
    spm_shoot_err = 0;
    ERR_UNLOCK;
}

void spm_shoot_param(const double prm[8], double param[8])
{
    param[0] = 1/prm[0];
    param[1] = 1/prm[1];
    param[2] = 1/prm[2];
    param[3] = prm[3];
    param[4] = prm[4];
    param[5] = prm[5];
    param[6] = prm[6];
    param[7] = prm[7];
}

/* fmg3_mexFunction of spm_diffeo.c */
#define eps 1.1921e-7
void spm_shoot_fmg_param(const double prm[8], double param[8])
{
    spm_shoot_param(prm, param);
    {
     /* Penalise absolute displacements slightly in case supplied Hessian is too small.
        Extra penalty based on sum of values in centre of difference operator, scaled
        by some arbitrary multiple of eps.
      */
        double v0   = param[0]*param[0],
               v1   = param[1]*param[1],
               v2   = param[2]*param[2],
               lam1 = param[4], lam2 = param[5],
               mu   = param[6], lam  = param[7],
               w000, wx000, wy000, wz000;
        w000      =  lam2*(6*(v0*v0+v1*v1+v2*v2) +8*(v0*v1+v0*v2+v1*v2)) +lam1*2*(v0+v1+v2);
        wx000     =  2*mu*(2*v0+v1+v2)/v0+2*lam + w000/v0;
        wy000     =  2*mu*(v0+2*v1+v2)/v1+2*lam + w000/v1;
        wz000     =  2*mu*(v0+v1+2*v2)/v2+2*lam + w000/v2;
        param[3] += (wx000 + wy000 + wz000)*eps*0.0001;
    }
}
#undef eps

int spm_shoot_fmg(const mwSize dm[3], const float *H, const float *g,
                  const double prm[8], int cyc, int nit, const float *v0,
                  float *v)
{
    mwSize n[3], m, i;
    double param[8];
    float *scratch;

    n[0] = dm[0];
    n[1] = dm[1];
    n[2] = dm[2];
    m    = 3*dm[0]*dm[1]*dm[2];
    spm_shoot_fmg_param(prm, param);

    /* the mex callocs the scratch space and copies (or zeroes) the start */
    scratch = (float *)calloc(fmg3_scratchsize(n, H!=0), sizeof(float));
    if (scratch == 0)
    {
        mexErrMsgTxt("Out of memory.");
        return -1;
    }
    if (v0 != 0)
    {
        if (v0 != v) memcpy(v, v0, m*sizeof(float));
    }
    else
        for(i=0; i<m; i++) v[i] = 0.0f;

    /* fmg3 reads but does not write H and g */
    fmg3(n, (float *)H, (float *)g, param, cyc, nit, v, scratch);
    free(scratch);
    return 0;
}

/* get_mat of spm_diffeo.c */
int spm_shoot_mat(const double M[16], float Mf[4][3])
{
    int i, j;

    for(i=0; i<3; i++)
        for(j=0; j<4; j++)
            Mf[j][i] = (float)M[i+4*j];

    if (M[3+4*0] != 0.0 || M[3+4*1] != 0.0 || M[3+4*2] != 0.0 || M[3+4*3] != 1.0)
    {
        mexErrMsgTxt("No perspective projections allowed.");
        return -1;
    }
    return 0;
}
