/*
 * pthreads support for the CAT-Surface port of the SPM26 shooting primitives;
 * it replaces SPM's OpenMP pragmas.
 *
 * Copyright (C) 2026 Christian Gaser, University of Jena.
 * Licence: GPL-2.0-or-later, as the SPM code it accompanies.
 */

#ifndef SHOOT_PARALLEL_H
#define SHOOT_PARALLEL_H

#include <stddef.h>

/*
 * Provided by libCAT (defined in Lib/CAT_Vol.c, declared in Include/CAT_Vol.h).
 * CAT_Vol.h is not included because it defines macros (GM, WM, CSF, index)
 * that clash with this code.  cat_parallel_run() starts the workers with
 * pthread_create and joins them; on Windows it runs them one after another.
 */
extern void cat_parallel_run(int nthreads, void *(*worker)(void *), void *args,
                             size_t argsz);
extern int  cat_parallel_nthreads(int n_units);

/* Thread control, see spm_shoot.h */
extern void   spm_shoot_set_threads(int n);
extern int    spm_shoot_get_threads(void);
extern void   spm_shoot_set_grain(double ns);
extern double spm_shoot_get_grain(void);

/*
 * Parallel for over the units [0,n): fun(ctx,lo,hi) is called on contiguous
 * blocks, each by one worker.  cost is the approximate run time of one unit
 * in nanoseconds; blocks are not made smaller than the grain.  In serial
 * mode (threads 0), or when there is too little work, fun(ctx,0,n) runs in
 * the calling thread, which is the SPM loop unchanged.
 */
typedef void (*spm_shoot_range_fun)(void *ctx, size_t lo, size_t hi);
extern int  spm_shoot_nworkers(size_t n, double cost);
extern void spm_shoot_parallel(size_t n, double cost, spm_shoot_range_fun fun,
                               void *ctx);

/*
 * Teams for the multi-colour Gauss-Seidel sweeps of shoot_regularisers.c.
 * In each sweep the slices k = k0, k0+st, ... < n of one colour are claimed
 * by the team's threads; a thread processes the slices it claims in
 * ascending order and spm_shoot_sweep_end() waits until every slice of the
 * sweep is done.  Slices of one colour are further apart than the stencil
 * reaches, except the last and the first, which can meet across the
 * periodic boundary: the thread that claims the first slice also processes
 * the last one, after it.  Every slice therefore sees exactly the values it
 * sees in SPM's serial loop, so the result does not depend on the number of
 * threads or on their timing.  Claiming is dynamic, so the team finishes
 * even if fewer threads run than were requested.
 */
typedef struct spm_shoot_team spm_shoot_team;
typedef struct
{
    spm_shoot_team *team;
    int             sweep;  /* sweep of the claims below, -1 before the first */
    int             head;   /* this thread processed the first slice */
    long            q;      /* slice index (within the colour) claimed last */
    long            ndone;  /* slices processed in this sweep */
} spm_shoot_thr;

extern int  spm_shoot_team_size(long nk, double cost);
extern void spm_shoot_team_run(int nthreads, int nsweeps,
                               void (*fun)(void *ctx, spm_shoot_thr *thr),
                               void *ctx);
extern int  spm_shoot_claim(spm_shoot_thr *thr, int sweep, long k, long k0,
                            long st, long n, long reach);
extern void spm_shoot_sweep_end(spm_shoot_thr *thr, int sweep, long k0,
                                long st, long n);

#endif /* SHOOT_PARALLEL_H */
