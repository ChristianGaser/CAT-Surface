/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Compare the C ports in Lib/CAT_Cat12Utils.c with the CAT12 MATLAB originals.
 *
 *   compare_reference DIR
 *
 * DIR holds manifest.txt and the raw files written by make_reference.m.  For
 * every case the C port runs on the same input and the output is compared
 * with MATLAB's: the largest absolute difference, the number of differing
 * voxels (binary outputs and the exact ports must give 0) and, for float
 * outputs, the largest absolute value of the reference.  A summary per
 * function follows.  The exit status is non-zero if a port that should be
 * exact differs or a smoothing differs by more than 1e-5 of the data range.
 *
 * Copyright Christian Gaser, University of Jena.
 */

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "CAT_Cat12Utils.h"
#include "CAT_Vol.h"

#define MAXF 32

/** Summary of one function. */
typedef struct
{
    char func[MAXF];
    int ncase;
    int ndiffcase;
    long ndiff;
    double maxabs;
    double maxrel;
} summary;

static summary sums[MAXF];
static int nsums = 0;

static summary *
get_summary(const char *func)
{
    int i;

    for (i = 0; i < nsums; i++)
        if (strcmp(sums[i].func, func) == 0)
            return &sums[i];
    memset(&sums[nsums], 0, sizeof(summary));
    strncpy(sums[nsums].func, func, MAXF - 1);
    return &sums[nsums++];
}

static void *
read_raw(const char *dir, const char *name, size_t n, size_t sz)
{
    char path[4096];
    FILE *f;
    void *p = malloc(n * sz);

    snprintf(path, sizeof(path), "%s/%s", dir, name);
    f = fopen(path, "rb");
    if (!f || !p || fread(p, sz, n, f) != n)
    {
        fprintf(stderr, "cannot read %s\n", path);
        exit(EXIT_FAILURE);
    }
    fclose(f);
    return p;
}

/** \brief Differences of two float volumes; NaN equals NaN, Inf equals Inf. */
static void
cmp_f32(const float *a, const float *b, size_t n, long *ndiff, double *maxabs, double *maxref)
{
    size_t i;

    *ndiff = 0;
    *maxabs = 0.0;
    *maxref = 0.0;
    for (i = 0; i < n; i++)
    {
        const int an = isnan(a[i]), bn = isnan(b[i]);
        if (!bn && !isinf(b[i]) && fabs(b[i]) > *maxref && b[i] != FLT_MAX)
            *maxref = fabs(b[i]);
        if (an && bn)
            continue;
        if (an != bn || ((isinf(a[i]) || isinf(b[i])) && a[i] != b[i]))
        {
            (*ndiff)++;
            *maxabs = INFINITY;
            continue;
        }
        if (a[i] != b[i])
        {
            double d = fabs((double)a[i] - (double)b[i]);
            (*ndiff)++;
            if (d > *maxabs)
                *maxabs = d;
        }
    }
}

static long
cmp_u8(const unsigned char *a, const unsigned char *b, size_t n)
{
    size_t i;
    long ndiff = 0;

    for (i = 0; i < n; i++)
        ndiff += (a[i] != 0) != (b[i] != 0);
    return ndiff;
}

int
main(int argc, char **argv)
{
    char path[4096], line[8192];
    FILE *mf;
    int status = 0, i;

    if (argc < 2)
    {
        fprintf(stderr, "Usage: %s DIR\n", argv[0]);
        return 1;
    }
    snprintf(path, sizeof(path), "%s/manifest.txt", argv[1]);
    mf = fopen(path, "r");
    if (!mf)
    {
        fprintf(stderr, "cannot open %s\n", path);
        return 1;
    }

    printf("%-34s %-9s %-13s %-28s %12s %9s %10s %7s\n", "case", "function", "dims",
           "parameters", "max|diff|", "#diff", "max|ref|", "ms");
    while (fgets(line, sizeof(line), mf))
    {
        char name[512], func[MAXF], fin[512], fmask[512], fout[512], fout2[512] = "";
        int d[3], nf;
        double p[4], vx[3], maxabs = 0.0, maxref = 0.0;
        long ndiff = 0;
        size_t n;
        unsigned char *mask = NULL;
        clock_t t0;
        double ms;
        int exact = 1, isfloat = 1;
        summary *s;

        nf = sscanf(line, "%511s %31s %d %d %d %lf %lf %lf %lf %511s %511s %511s %511s", name,
                    func, &d[0], &d[1], &d[2], &p[0], &p[1], &p[2], &p[3], fin, fmask, fout,
                    fout2);
        if (nf < 12)
            continue;
        n = (size_t)d[0] * d[1] * d[2];
        if (strcmp(fmask, "-") != 0)
            mask = (unsigned char *)read_raw(argv[1], fmask, n, 1);
        vx[0] = p[1];
        vx[1] = p[2];
        vx[2] = p[3];

        t0 = clock();
        if (strcmp(func, "smooth3X") == 0 || strcmp(func, "div") == 0 ||
            strcmp(func, "localstat") == 0 || strcmp(func, "laplace") == 0)
        {
            float *in = (float *)read_raw(argv[1], fin, n, 4);
            float *ref = (float *)read_raw(argv[1], fout, n, 4);
            float *out = (float *)malloc(n * sizeof(float));
            size_t k;

            if (strcmp(func, "smooth3X") == 0)
            {
                memcpy(out, in, n * sizeof(float));
                cat12_vol_smooth3X(out, d, p[0]);
                exact = 0;
            }
            else if (strcmp(func, "div") == 0)
                cat12_vol_div(in, d, out);
            else if (strcmp(func, "localstat") == 0)
                cat12_vol_localstat_mean(in, mask, d, (int)p[0], out);
            else
            {
                /* cat_vol_laplace3R(u/10, M, TH)*10 in single precision */
                for (k = 0; k < n; k++)
                    out[k] = (float)((double)in[k] / 10.0);
                cat12_vol_laplace3R(out, mask, d, p[0]);
                for (k = 0; k < n; k++)
                    out[k] = (float)((double)out[k] * 10.0);
            }
            ms = 1000.0 * (double)(clock() - t0) / CLOCKS_PER_SEC;
            cmp_f32(out, ref, n, &ndiff, &maxabs, &maxref);

            if (strcmp(func, "laplace") == 0)
            {
                /* libCAT's laplace3R() with TH in double, for the record */
                long nd2;
                double ma2, mr2;
                int dd[3] = {d[0], d[1], d[2]};
                for (k = 0; k < n; k++)
                    out[k] = (float)((double)in[k] / 10.0);
                laplace3R(out, mask, dd, p[0]);
                for (k = 0; k < n; k++)
                    out[k] = (float)((double)out[k] * 10.0);
                cmp_f32(out, ref, n, &nd2, &ma2, &mr2);
                printf("    (libCAT laplace3R with TH in double: max|diff| %g, #diff %ld)\n", ma2,
                       nd2);
            }
            free(out);
            free(ref);
            free(in);
        }
        else if (strcmp(func, "vbdist") == 0)
        {
            float *in = (float *)read_raw(argv[1], fin, n, 4);
            float *ref = (float *)read_raw(argv[1], fout, n, 4);
            unsigned int *refi = (unsigned int *)read_raw(argv[1], fout2, n, 4);
            float *D = (float *)malloc(n * sizeof(float));
            unsigned int *I = (unsigned int *)malloc(n * sizeof(unsigned int));
            long ndi = 0;
            size_t k;

            cat12_vbdist(in, mask, d, p, D, I);
            ms = 1000.0 * (double)(clock() - t0) / CLOCKS_PER_SEC;
            cmp_f32(D, ref, n, &ndiff, &maxabs, &maxref);
            for (k = 0; k < n; k++)
                ndi += (I[k] + 1 != refi[k]);
            if (ndi)
                printf("    (nearest index differs at %ld voxels)\n", ndi);
            ndiff += ndi;
            free(I);
            free(D);
            free(refi);
            free(ref);
            free(in);
        }
        else
        {
            unsigned char *vol = (unsigned char *)read_raw(argv[1], fin, n, 1);
            unsigned char *ref = (unsigned char *)read_raw(argv[1], fout, n, 1);

            isfloat = 0;
            if (strcmp(func, "lab") == 0)
                cat12_vol_morph_lab(vol, d);
            else if (strcmp(func, "labmulti") == 0)
                cat12_vol_morph_lab_multi(vol, d, p[0], p[1]);
            else if (strcmp(func, "dilate") == 0)
                cat12_vol_morph_dilate(vol, d, p[0], vx);
            else if (strcmp(func, "open") == 0)
                cat12_vol_morph_open(vol, d, p[0], vx);
            else if (strcmp(func, "ldc") == 0)
                cat12_vol_morph_labdistclose(vol, d, p[0], vx);
            else if (strcmp(func, "ldo") == 0)
                cat12_vol_morph_labdistopen(vol, d, p[0], vx);
            else
            {
                fprintf(stderr, "unknown function %s\n", func);
                return 1;
            }
            ms = 1000.0 * (double)(clock() - t0) / CLOCKS_PER_SEC;
            ndiff = cmp_u8(vol, ref, n);
            maxabs = ndiff ? 1.0 : 0.0;
            free(ref);
            free(vol);
        }
        free(mask);

        {
            char dims[32], par[64];
            snprintf(dims, sizeof(dims), "%dx%dx%d", d[0], d[1], d[2]);
            if (strcmp(func, "lab") == 0 || strcmp(func, "div") == 0)
                par[0] = 0;
            else if (strcmp(func, "smooth3X") == 0 || strcmp(func, "localstat") == 0 ||
                     strcmp(func, "laplace") == 0)
                snprintf(par, sizeof(par), "%.4g", p[0]);
            else if (strcmp(func, "labmulti") == 0)
                snprintf(par, sizeof(par), "[%g %g]", p[0], p[1]);
            else if (strcmp(func, "vbdist") == 0)
                snprintf(par, sizeof(par), "vx [%g %g %g]", p[0], p[1], p[2]);
            else
                snprintf(par, sizeof(par), "%.4g vx [%g %g %g]", p[0], p[1], p[2], p[3]);
            if (isfloat)
                printf("%-34s %-9s %-13s %-28s %12.4g %9ld %10.4g %7.1f\n", name, func, dims, par,
                       maxabs, ndiff, maxref, ms);
            else
                printf("%-34s %-9s %-13s %-28s %12s %9ld %10s %7.1f\n", name, func, dims, par, "",
                       ndiff, "", ms);
        }

        s = get_summary(func);
        s->ncase++;
        s->ndiff += ndiff;
        if (ndiff)
            s->ndiffcase++;
        if (maxabs > s->maxabs)
            s->maxabs = maxabs;
        if (maxref > 0 && maxabs / maxref > s->maxrel)
            s->maxrel = maxabs / maxref;
        if (exact ? ndiff != 0 : maxabs > 1e-5 * maxref)
            status = 1;
    }
    fclose(mf);

    printf("\n%-10s %6s %14s %12s %14s %14s\n", "function", "cases", "cases w/ diff", "#diff",
           "max|diff|", "max rel diff");
    for (i = 0; i < nsums; i++)
        printf("%-10s %6d %14d %12ld %14.4g %14.4g\n", sums[i].func, sums[i].ncase,
               sums[i].ndiffcase, sums[i].ndiff, sums[i].maxabs, sums[i].maxrel);
    printf("\n%s\n", status ? "FAILED" : "OK");
    return status;
}
