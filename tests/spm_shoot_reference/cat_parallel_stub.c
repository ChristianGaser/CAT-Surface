/*
 * cat_parallel_run() and cat_parallel_nthreads() as in Lib/CAT_Vol.c, for
 * building the spm_shoot tests by hand without libCAT (see build.sh).  Do not
 * link it together with libCAT, which defines both.
 *
 * Copyright (C) 2026 Christian Gaser, University of Jena.
 */

#include <stdlib.h>
#if !defined(_WIN32) && !defined(_WIN64)
#include <pthread.h>
#endif

#ifndef MAX_NTHREADS
#define MAX_NTHREADS 4
#endif

void cat_parallel_run(int nthreads, void *(*worker)(void *), void *args, size_t argsz)
{
    int t;
#if defined(_WIN32) || defined(_WIN64)
    for (t = 0; t < nthreads; ++t)
        worker((char *)args + (size_t)t * argsz);
#else
    pthread_t *threads = (pthread_t *)calloc((size_t)nthreads, sizeof(pthread_t));
    if (!threads)
    {
        for (t = 0; t < nthreads; ++t)
            worker((char *)args + (size_t)t * argsz);
        return;
    }
    for (t = 0; t < nthreads; ++t)
        pthread_create(&threads[t], NULL, worker, (char *)args + (size_t)t * argsz);
    for (t = 0; t < nthreads; ++t)
        pthread_join(threads[t], NULL);
    free(threads);
#endif
}

int cat_parallel_nthreads(int n_units)
{
    int n = (n_units < MAX_NTHREADS) ? n_units : MAX_NTHREADS;
    return (n < 1) ? 1 : n;
}
