/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "ParseArgv.h"
#include "CAT_NiftiLib.h"
#include "CAT_Calc.h"
#include "CAT_ArgvCheck.h"
#include "CAT_VolumeStats.h"

char *expression = NULL;
char *datatype   = NULL;
int   rescale    = 0;
int   verbose    = 0;

extern int ParseArgv(int *argcPtr, char **argv, ArgvInfo *argTable, int flags);

static ArgvInfo argTable[] = {
    {"-expression", ARGV_STRING, (char *) 1, (char *) &expression,
         "Formula to evaluate at every voxel (required)."},
    {"-dt", ARGV_STRING, (char *) 1, (char *) &datatype,
         "Output data type: uint8, int8, uint16, int16, uint32, int32,\n\
         int64, float32 (default) or float64."},
    {"-rescale", ARGV_CONSTANT, (char *) 1, (char *) &rescale,
         "Rescale the result to fill the range of an integer output type,\n\
         storing the scale factor in the NIfTI header.  Without this the\n\
         computed values are written verbatim (and rounded, for an integer\n\
         type), which is what a calculator should normally do.  It matters\n\
         when the result is fractional but the output type is not."},
    {"-v", ARGV_CONSTANT, (char *) 1, (char *) &verbose,
         "Be verbose."},
    {NULL, ARGV_END, NULL, NULL, NULL}
};

/* Names accepted by -dt, mapped to the NIfTI codes write_nifti_double() takes */
static const struct {
    const char *name;
    int         code;
} datatypes[] = {
    {"uint8",   DT_UINT8},
    {"int8",    DT_INT8},
    {"uint16",  DT_UINT16},
    {"int16",   DT_INT16},
    {"uint32",  DT_UINT32},
    {"int32",   DT_INT32},
    {"int64",   DT_INT64},
    {"float32", DT_FLOAT32},
    {"float64", DT_FLOAT64},
    {NULL,      0}
};

static int
datatype_from_name(const char *name)
{
    int i;

    for (i = 0; datatypes[i].name != NULL; i++) {
        if (strcmp(datatypes[i].name, name) == 0)
            return datatypes[i].code;
    }
    return -1;
}

static void usage(char *executable)
{
    char *usage_str = "\n\
Usage: %s -expression \"<formula>\" [options] <in1.nii> [<in2.nii> ...] <out.nii>\n\
\n\
    Voxel-wise image calculator, in the spirit of SPM's spm_imcalc.  The\n\
    formula is evaluated once per voxel over the input images, which must all\n\
    share a grid, and the result is written to <out.nii>.\n\
\n\
    The inputs are addressed as i1, i2, ... in the order given on the command\n\
    line.  A 4D file counts as one input per volume: its volumes take the\n\
    next numbers in turn, so with a single 4D file i3 is its third volume,\n\
    and 3D and 4D files may be mixed.  The result is always one 3D image.\n\
    Numeric literals and the constants pi, inf and nan are available,\n\
    as are the operators\n\
\n\
        + - * / ^        arithmetic ('^' is right associative)\n\
        < <= > >= == ~=  comparison, yielding 1.0 or 0.0\n\
        & |  ! ~         logical and, or, not\n\
\n\
    and the element-wise functions abs, sqrt, exp, log, log2, log10, sin,\n\
    cos, tan, asin, acos, atan, floor, ceil, round, sign, isnan, isinf,\n\
    isfinite, min, max, pow, atan2 and mod.  Because everything here is\n\
    element-wise already, the MATLAB spellings .* ./ .^ are accepted as\n\
    synonyms of * / ^, so an expression can be pasted over from spm_imcalc\n\
    unchanged.  Comparisons returning 1.0 or 0.0 make \"(i1>0.5).*i2\" the\n\
    idiomatic way to mask one image with another.\n\
\n\
Matrix mode:\n\
\n\
    The variable X stands for the whole vector of input values at the current\n\
    voxel, over all files and all volumes of 4D files, which is what makes a\n\
    formula independent of how many images are passed.  Three reductions\n\
    consume it:\n\
\n\
        mean(X)     mean across images\n\
        median(X)   median across images\n\
        std(X)      sample standard deviation across images (normalised n-1)\n\
\n\
    X is an ordinary value, not a token that only these three accept, so a\n\
    reduction may take a computed argument: std(X-i1) is the spread about the\n\
    first scan, and mean(X>0.5) is the fraction of scans above a threshold.\n\
\n\
    The reductions skip non-finite values rather than propagating them, as\n\
    get_mean_double() and get_std_double() do elsewhere in the library, so a\n\
    NaN arising inside the formula -- from log() of a negative, or 0/0 --\n\
    drops out of the reduction instead of blanking the result.  Note that the\n\
    NIfTI reader replaces every non-finite value in an input file with zero\n\
    before it gets here, so a NaN stored in a file arrives as a zero and is\n\
    counted as one.  A voxel at which nothing finite remains yields NaN, as\n\
    does std() of fewer than two finite values.\n\
\n\
    All inputs are held in memory as doubles, so a long file list or a long\n\
    4D series costs 8 bytes per voxel per volume.\n\
\n\
Options:\n\
    -expression <str>  Formula to evaluate (required).\n\
    -dt <str>          Output data type (default: float32).\n\
    -rescale           Scale an integer result to fill the type's range.\n\
    -v                 Be verbose.\n\
\n\
Examples:\n\
    %s -expression \"i1 - i2\" scan1.nii scan2.nii difference.nii\n\
    %s -expression \"(i1>0.5).*i2\" mask.nii t1.nii masked.nii\n\
    %s -expression \"i1./mean(X)\" -v sub*.nii ratio.nii\n\
    %s -expression \"median(X)\" sub*.nii median.nii\n\
    %s -expression \"std(X)\" sub*.nii sd.nii\n\
    %s -expression \"median(X)\" zmaps4d.nii median.nii\n\
    %s -expression \"mean(abs(X)>1.96)\" zmaps4d.nii share_extreme.nii\n\
    %s -expression \"mean(X)>0.5\" -dt uint8 sub*.nii majority_mask.nii\n\n";

    fprintf(stderr, usage_str, executable, executable, executable, executable,
            executable, executable, executable, executable, executable);
}

/* Main program */
int
main(int argc, char *argv[])
{
    char **infiles, *outfile, err[256];
    int i, k, n_files, n_img, dt;
    double **images = NULL, **buffers = NULL, *result = NULL, slope;
    nifti_image *hdr;
    CAT_VolAverageResult scan;
    CAT_CalcExpr *expr = NULL;

    if (ParseArgv(&argc, argv, argTable, 0) || (argc < 3) ||
        (expression == NULL)) {
        usage(argv[0]);
        (void) fprintf(stderr, "     %s -help\n\n", argv[0]);
        exit(EXIT_FAILURE);
    }

    /* a leftover that looks like an option is a typo or a removed flag;
       ParseArgv() keeps it, so it would be read as a file name below */
    if (cat_check_unknown_options(argc, argv, argv[0]))
    {
        usage(argv[0]);
        (void) fprintf(stderr, "     %s -help\n\n", argv[0]);
        exit(EXIT_FAILURE);
    }

    n_files = argc - 2;
    infiles = &argv[1];
    outfile = argv[argc - 1];

    dt = (datatype != NULL) ? datatype_from_name(datatype) : DT_FLOAT32;
    if (dt < 0) {
        fprintf(stderr, "Unknown data type '%s'.\n", datatype);
        exit(EXIT_FAILURE);
    }

    /* The headers only: every volume is one input, so a 4D file contributes
       one per volume, and the formula has to be checked against that number.
       Files of another matrix size or with more than four dimensions are
       refused here. */
    if (CAT_VolAverageScan(infiles, n_files, 0, &scan))
        exit(EXIT_FAILURE);
    n_img = scan.n_samples;

    /* No reslicing here: an expression evaluated across grids that do not
       line up is meaningless, so refuse a different voxel size rather than
       guess. */
    for (i = 1; i < n_files; i++) {
        hdr = nifti_image_read(infiles[i], 0);
        k = (hdr != NULL) && equal_image_dimensions(scan.header, hdr);
        if (hdr != NULL)
            nifti_image_free(hdr);
        if (!k) {
            fprintf(stderr, "Error: %s does not share the grid of %s.\n",
                    infiles[i], infiles[0]);
            goto fail;
        }
    }

    /* Parse before reading any data: a typo in the formula should not cost
       the time it takes to load a long file list. */
    expr = CAT_CalcParse(expression, n_img, err, sizeof(err));
    if (expr == NULL) {
        fprintf(stderr, "Error in expression \"%s\":\n  %s.\n", expression, err);
        goto fail;
    }

    if (verbose) {
        fprintf(stdout, "Expression: %s\n", expression);
        fprintf(stdout, "Inputs:     %d image%s%s\n", n_img,
                n_img == 1 ? "" : "s",
                CAT_CalcUsesMatrix(expr) ? " (matrix mode)" : "");
        for (i = 0, k = 1; i < n_files; k += scan.n_volumes[i], i++) {
            if (scan.n_volumes[i] > 1)
                fprintf(stdout, "  i%d-i%d %s (%d volumes)\n", k,
                        k + scan.n_volumes[i] - 1, infiles[i], scan.n_volumes[i]);
            else
                fprintf(stdout, "  i%-3d %s\n", k, infiles[i]);
        }
    }

    if (CAT_VolAverageRead(infiles, n_files, &scan, &images, &buffers))
        goto fail;

    result = (double *) malloc(sizeof(double) * scan.n_voxels);
    if (result == NULL) {
        fprintf(stderr, "Memory allocation error.\n");
        goto fail;
    }

    if (CAT_CalcApply(expr, images, n_img, scan.n_voxels, result) != 0)
        goto fail;

    /* slope 0.0 asks write_nifti_double() to fit the data to the output type
       and record the factor; 1.0 writes the values as they are. */
    slope = rescale ? 0.0 : 1.0;

    if (verbose)
        fprintf(stdout, "Output:     %s\n", outfile);

    if (!write_nifti_double(outfile, result, dt, slope, scan.dims, scan.voxelsize,
                            scan.header))
        goto fail;

    CAT_VolAverageReadFree(images, buffers, n_files);
    free(result);
    CAT_CalcFree(expr);
    CAT_VolAverageFree(&scan);

    return EXIT_SUCCESS;

fail:
    CAT_VolAverageReadFree(images, buffers, n_files);
    free(result);
    CAT_CalcFree(expr);
    CAT_VolAverageFree(&scan);

    return EXIT_FAILURE;
}
