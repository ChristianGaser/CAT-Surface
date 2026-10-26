/*
 * pthreads support for the CAT-Surface port of the SPM26 shooting primitives;
 * it replaces SPM's OpenMP pragmas.  See shoot_parallel.h.
 *
 * Copyright (C) 2026 Christian Gaser, University of Jena.
 * Licence: GPL-2.0-or-later, as the SPM code it accompanies.
 */

#include <stdlib.h>
#include "shoot_parallel.h"

#if !defined(_WIN32) && !defined(_WIN64)
#include <pthread.h>
#define SPM_SHOOT_PTHREADS 1
#endif

/* Upper limit for the number of threads */
#define SPM_SHOOT_MAX_THREADS 64

/* Default grain: minimum work per worker (ns).  Starting and joining a
   thread costs 10-20 us and a sweep barrier 2-25 us on an M2 Ultra, so a
   worker needs about 0.1 ms of work to pay for itself. */
#define SPM_SHOOT_GRAIN 1.0e5

static int    spm_shoot_threads = 0;
static double spm_shoot_grain   = SPM_SHOOT_GRAIN;

void spm_shoot_set_threads(int n)
{
    if (n < 0) n = 0;
    if (n > SPM_SHOOT_MAX_THREADS) n = SPM_SHOOT_MAX_THREADS;
    spm_shoot_threads = n;
}

int spm_shoot_get_threads(void)
{
    return spm_shoot_threads;
}

void spm_shoot_set_grain(double ns)
{
    spm_shoot_grain = (ns >= 0.0) ? ns : SPM_SHOOT_GRAIN;
}

double spm_shoot_get_grain(void)
{
    return spm_shoot_grain;
}

/* Number of workers for n units of cost ns each: 0 in serial mode, else
   1..threads, limited by n and by the grain. */
int spm_shoot_nworkers(size_t n, double cost)
{
    double nw;

    if (spm_shoot_threads <= 0) return 0;
#ifndef SPM_SHOOT_PTHREADS
    /* cat_parallel_run() would run the blocks one after another */
    (void)n; (void)cost;
    return 1;
#endif
    if (spm_shoot_grain > 0.0)
        nw = (double)n*cost/spm_shoot_grain;
    else
        nw = (double)n;
    if (nw > (double)spm_shoot_threads) nw = (double)spm_shoot_threads;
    if (nw > (double)n) nw = (double)n;
    return (nw < 1.0) ? 1 : (int)nw;
}

/***************************************************************************/
/* Parallel for */

typedef struct
{
    spm_shoot_range_fun fun;
    void               *ctx;
    size_t              lo, hi;
} range_slot;

static void *range_worker(void *arg)
{
    range_slot *s = (range_slot *)arg;
    if (s->lo < s->hi) s->fun(s->ctx, s->lo, s->hi);
    return NULL;
}

void spm_shoot_parallel(size_t n, double cost, spm_shoot_range_fun fun, void *ctx)
{
    range_slot slot[SPM_SHOOT_MAX_THREADS];
    int        t, nw;

    nw = spm_shoot_nworkers(n, cost);
    if (nw <= 1)
    {
        if (n > 0) fun(ctx, 0, n);
        return;
    }
    for (t = 0; t < nw; t++)
    {
        slot[t].fun = fun;
        slot[t].ctx = ctx;
        slot[t].lo  = n*(size_t)t/(size_t)nw;
        slot[t].hi  = n*(size_t)(t+1)/(size_t)nw;
    }
    cat_parallel_run(nw, range_worker, slot, sizeof(range_slot));
}

/***************************************************************************/
/* Sweep teams */

struct spm_shoot_team
{
#ifdef SPM_SHOOT_PTHREADS
    pthread_mutex_t mtx;
    pthread_cond_t  cnd;
#endif
    long *next;     /* per sweep: next slice index to hand out */
    long *done;     /* per sweep: slices finished */
};

typedef struct
{
    void          (*fun)(void *ctx, spm_shoot_thr *thr);
    void           *ctx;
    spm_shoot_thr   thr;
} team_slot;

#ifdef SPM_SHOOT_PTHREADS
static void *team_worker(void *arg)
{
    team_slot *s = (team_slot *)arg;
    s->fun(s->ctx, &s->thr);
    return NULL;
}
#endif

/* Number of threads for a team whose sweeps have nk slices of cost ns each:
   0 (serial code) unless at least two threads get work.  A team starts its
   threads once per call and then only meets at the end of each sweep, which
   costs about a quarter of starting them, hence a quarter of the grain per
   sweep (measured with fmg3 on 113x137x113: 0.37 s with the full grain,
   0.33 s with a quarter, 8 threads). */
int spm_shoot_team_size(long nk, double cost)
{
#ifdef SPM_SHOOT_PTHREADS
    int nt;
    if (nk < 2) return 0;
    nt = spm_shoot_nworkers((size_t)nk, 4.0*cost);
    return (nt >= 2) ? nt : 0;
#else
    (void)nk; (void)cost;
    return 0; /* cat_parallel_run() is sequential without pthreads */
#endif
}

void spm_shoot_team_run(int nthreads, int nsweeps,
                        void (*fun)(void *ctx, spm_shoot_thr *thr), void *ctx)
{
#ifdef SPM_SHOOT_PTHREADS
    spm_shoot_team tm;
    team_slot      slot[SPM_SHOOT_MAX_THREADS];
    long          *buf;
    int            t;

    if (nthreads > SPM_SHOOT_MAX_THREADS) nthreads = SPM_SHOOT_MAX_THREADS;
    buf = (nthreads >= 2 && nsweeps > 0) ? (long *)calloc((size_t)nsweeps*2, sizeof(long)) : NULL;
    if (buf == NULL)
    {
        fun(ctx, NULL);
        return;
    }
    if (pthread_mutex_init(&tm.mtx, NULL) != 0)
    {
        free(buf);
        fun(ctx, NULL);
        return;
    }
    if (pthread_cond_init(&tm.cnd, NULL) != 0)
    {
        pthread_mutex_destroy(&tm.mtx);
        free(buf);
        fun(ctx, NULL);
        return;
    }
    tm.next = buf;
    tm.done = buf + nsweeps;
    for (t = 0; t < nthreads; t++)
    {
        slot[t].fun       = fun;
        slot[t].ctx       = ctx;
        slot[t].thr.team  = &tm;
        slot[t].thr.sweep = -1;
        slot[t].thr.head  = 0;
        slot[t].thr.q     = -1;
        slot[t].thr.ndone = 0;
    }
    cat_parallel_run(nthreads, team_worker, slot, sizeof(team_slot));
    pthread_cond_destroy(&tm.cnd);
    pthread_mutex_destroy(&tm.mtx);
    free(buf);
#else
    (void)nthreads; (void)nsweeps;
    fun(ctx, NULL);
#endif
}

#ifdef SPM_SHOOT_PTHREADS
static long nclass(long k0, long st, long n)
{
    return (n > k0) ? (n - k0 + st - 1)/st : 0;
}

static long take(spm_shoot_team *tm, int sweep)
{
    long q;
    pthread_mutex_lock(&tm->mtx);
    q = tm->next[sweep]++;
    pthread_mutex_unlock(&tm->mtx);
    return q;
}
#endif

/* Does this thread process slice k of the colour k0, k0+st, ... < n in this
   sweep?  reach is how far the stencil reaches along the slice dimension. */
int spm_shoot_claim(spm_shoot_thr *thr, int sweep, long k, long k0, long st,
                    long n, long reach)
{
#ifdef SPM_SHOOT_PTHREADS
    long nk, q, neff;
    int  wraps;

    nk    = nclass(k0, st, n);
    q     = (k - k0)/st;
    wraps = (nk >= 2) && (n - (nk - 1)*st <= reach); /* last slice meets the first */
    neff  = nk - wraps;                              /* slices handed out */

    if (thr->sweep != sweep)
    {
        thr->sweep = sweep;
        thr->head  = 0;
        thr->q     = -1;
        thr->ndone = 0;
    }
    if (q < neff)
    {
        /* A claim is only taken here, at a position not beyond what has
           been handed out, so it is never behind the loop. */
        if (thr->q < q) thr->q = take(thr->team, sweep);
        if (thr->q != q) return 0;
        if (q == 0) thr->head = 1;
        thr->ndone++;
        return 1;
    }
    if (wraps && q == nk - 1 && thr->head)
    {
        thr->ndone++;
        return 1;
    }
    return 0;
#else
    (void)thr; (void)sweep; (void)k; (void)k0; (void)st; (void)n; (void)reach;
    return 1;
#endif
}

/* Report the slices this thread finished and wait for the rest of the sweep */
void spm_shoot_sweep_end(spm_shoot_thr *thr, int sweep, long k0, long st, long n)
{
#ifdef SPM_SHOOT_PTHREADS
    spm_shoot_team *tm = thr->team;
    long nk, mine;

    nk   = nclass(k0, st, n);
    mine = (thr->sweep == sweep) ? thr->ndone : 0;
    pthread_mutex_lock(&tm->mtx);
    tm->done[sweep] += mine;
    if (tm->done[sweep] >= nk)
        pthread_cond_broadcast(&tm->cnd);
    else
        while (tm->done[sweep] < nk)
            pthread_cond_wait(&tm->cnd, &tm->mtx);
    pthread_mutex_unlock(&tm->mtx);
#else
    (void)thr; (void)sweep; (void)k0; (void)st; (void)n;
#endif
}
