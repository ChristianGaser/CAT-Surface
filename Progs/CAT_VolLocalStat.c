/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 *
 */

#include <float.h>
#include <stdlib.h>
#include <string.h>

#if !defined(_WIN32) && !defined(_WIN64)
#include <libgen.h>
#endif

#include "ParseArgv.h"
#include "CAT_NiftiLib.h"
#include "CAT_Vol.h"
#include "CAT_Sheetness.h"

/* Defaults */
int    stat_func          = F_MEAN;
int    dist               = 1;
int    iters              = 1;
int    use_euclidean_dist = 0;
int    oriented           = 0;
char  *guide_file         = NULL;
double sheet_sigma_min    = 0.3;
double sheet_sigma_max    = 1.0;
int    sheet_scales       = 3;
int    sheet_polarity     = 0;
double sheet_strength     = 1.0;
int    verbose            = 0;

static
ArgvInfo argTable[] = {
    {"-stat", ARGV_INT, (char *) 1, (char *) &stat_func,
         "Statistic function: 0=mean, 1=min, 2=max, 3=std, 7=median, 12=close, 13=open (default: 0=mean)."},
    {"-dist", ARGV_INT, (char *) 1, (char *) &dist,
         "Search distance from voxel center in voxels (1..10, default: 1)."},
    {"-iter", ARGV_INT, (char *) 1, (char *) &iters,
         "Number of iterations (default: 1)."},
    {"-euclid", ARGV_CONSTANT, (char *) 1, (char *) &use_euclidean_dist,
         "Use Euclidean distance instead of block distance (default: block)."},
    {"-oriented", ARGV_CONSTANT, (char *) 1, (char *) &oriented,
         "Run the median over a sheetness-oriented neighbourhood instead of an\n\
         isotropic one (-stat 7 only).  An isotropic median penalizes boundary\n\
         area, so it removes thin structures whichever side of the label\n\
         boundary they lie on: the same filter that opens a glued sulcus closes\n\
         a cerebellar fissure, and tuning trades one against the other.  The\n\
         oriented variant admits only neighbours lying in the plane of the local\n\
         sheet, so it averages along a thin structure and never across it.\n\
         Where no sheet is detected every neighbour is admitted and the result\n\
         is identical to the isotropic median.  Uses a fixed 3x3x3 neighbourhood,\n\
         so -dist and -euclid are ignored; -iter still applies."},
    {"-guide", ARGV_STRING, (char *) 1, (char *) &guide_file,
         "Volume the orientation is estimated from (default: the input itself).\n\
         Pass the intensity image when filtering a label map -- a T1 still shows\n\
         the dip through a tight sulcus long after the label map committed to\n\
         pure GM.  Must have the same dimensions as the input."},
    {"-sheet-sigma-min", ARGV_FLOAT, (char *) 1, (char *) &sheet_sigma_min,
         "Smallest scale of the sheetness filter in mm (default: 0.3)."},
    {"-sheet-sigma-max", ARGV_FLOAT, (char *) 1, (char *) &sheet_sigma_max,
         "Largest scale of the sheetness filter in mm (default: 1.0)."},
    {"-sheet-scales", ARGV_INT, (char *) 1, (char *) &sheet_scales,
         "Number of log-spaced sheetness scales (default: 3)."},
    {"-sheet-polarity", ARGV_INT, (char *) 1, (char *) &sheet_polarity,
         "1 = bright sheets, -1 = dark sheets, 0 = either (default)."},
    {"-sheet-strength", ARGV_FLOAT, (char *) 1, (char *) &sheet_strength,
         "How far the filter may deviate from isotropic, 0..1 (default: 1.0).\n\
         0 reproduces the isotropic filter exactly."},
    {"-v", ARGV_CONSTANT, (char *) 1, (char *) &verbose,
         "Be verbose."},
    {NULL, ARGV_END, NULL, NULL, NULL}
};

void usage(char *executable)
{
    char *usage_str = "\n\
Usage: %s [options] <input.nii> [<output.nii>]\n\
\n\
    Apply local statistics to a 3D NIfTI volume. For each voxel the chosen\n\
    statistic is computed over a neighbourhood defined by the search distance\n\
    and, optionally, restricted to an Euclidean sphere.\n\
\n\
    Available statistics (-stat):\n\
        0  mean     (default)\n\
        1  min      (erosion)\n\
        2  max      (dilation)\n\
        3  std\n\
        7  median\n\
       12  close    (grey closing: max then min)\n\
       13  open     (grey opening: min then max)\n\
\n\
Options:\n\
    -stat   <int>    Statistic function (default: 0=mean).\n\
    -dist   <int>    Search distance 1..10 (default: 1).\n\
    -iter   <int>    Number of iterations (default: 1).\n\
    -euclid          Use Euclidean distance instead of block distance.\n\
    -oriented        Sheetness-oriented median (-stat 7 only); see below.\n\
    -guide  <file>   Volume the orientation is estimated from.\n\
    -v               Be verbose.\n\
\n\
    With -oriented the median runs over a neighbourhood oriented by a Hessian\n\
    sheetness filter (see CAT_VolSheetness): a neighbour at offset d is admitted\n\
    only when 1 - s*(dhat.n)^2 > 0.5, with s the local sheetness and n the local\n\
    sheet normal.  At s = 0 every neighbour is admitted and the operator is the\n\
    plain isotropic median; at s = 1 only offsets within 45 degrees of the sheet\n\
    plane survive, so the filter can no longer close a thin structure.\n\
\n\
Example:\n\
    %s -stat 7 -dist 2 -iter 3 input.nii output.nii\n\
    %s -stat 7 -oriented -guide t1_corr.nii label.nii label_filtered.nii\n\n";

    fprintf(stderr, usage_str, executable, executable, executable);
}

/* Return a short name for the selected statistic (used for auto-naming). */
static const char *
stat_name(int func)
{
    switch (func) {
    case F_MEAN:   return "mean";
    case F_MIN:    return "min";
    case F_MAX:    return "max";
    case F_STD:    return "std";
    case F_MEDIAN: return oriented ? "orimedian" : "median";
    case F_CLOSE:  return "close";
    case F_OPEN:   return "open";
    default:       return "stat";
    }
}

/* Main program */
int
main(int argc, char *argv[])
{
    char *infile, outfile[1024];
    int dims[3];
    float *input;
    double voxelsize[3];
    nifti_image *nii_ptr;
    float *guide = NULL, *sheetness = NULL, *normal = NULL;
    nifti_image *guide_ptr = NULL;

    /* Parse arguments */
    if (ParseArgv(&argc, argv, argTable, 0) || (argc < 2)) {
        usage(argv[0]);
        (void) fprintf(stderr, "     %s -help\n\n", argv[0]);
        exit(EXIT_FAILURE);
    }

    infile = argv[1];

    if (oriented && stat_func != F_MEDIAN) {
        fprintf(stderr, "Error: -oriented is only implemented for -stat 7 (median).\n");
        exit(EXIT_FAILURE);
    }

    /* Validate stat_func */
    if (stat_func != F_MEAN && stat_func != F_MIN && stat_func != F_MAX &&
        stat_func != F_STD  && stat_func != F_MEDIAN &&
        stat_func != F_CLOSE && stat_func != F_OPEN) {
        fprintf(stderr, "Error: unsupported statistic function %d.\n", stat_func);
        exit(EXIT_FAILURE);
    }

    if (verbose) {
        fprintf(stdout, "Applying local %s filter (dist=%d, iter=%d, %s) to %s\n",
                stat_name(stat_func), dist, iters,
                use_euclidean_dist ? "euclidean" : "block", infile);
    }

    /* Read input NIfTI */
    nii_ptr = read_nifti_float(infile, &input, 0);
    if (nii_ptr == NULL) {
        fprintf(stderr, "Error reading %s.\n", infile);
        return EXIT_FAILURE;
    }

    voxelsize[0] = nii_ptr->dx;
    voxelsize[1] = nii_ptr->dy;
    voxelsize[2] = nii_ptr->dz;
    dims[0] = nii_ptr->nx;
    dims[1] = nii_ptr->ny;
    dims[2] = nii_ptr->nz;

    /* Estimate the orientation field the oriented median needs.  It is
       recomputed here rather than read from a file because a 3-vector per
       voxel needs a 4-D image, and the Hessian pass is cheap next to the
       cost of moving such a field through disk. */
    if (oriented) {
        int nvox = dims[0] * dims[1] * dims[2];
        CAT_SheetnessOpts sopts;
        int i;

        if (guide_file) {
            guide_ptr = read_nifti_float(guide_file, &guide, 0);
            if (guide_ptr == NULL) {
                fprintf(stderr, "Error reading %s.\n", guide_file);
                return EXIT_FAILURE;
            }
            if (guide_ptr->nx != dims[0] || guide_ptr->ny != dims[1] ||
                guide_ptr->nz != dims[2]) {
                fprintf(stderr, "Guide volume must have the same dimensions as the input.\n");
                return EXIT_FAILURE;
            }
        }

        sheetness = (float *) malloc(sizeof(float) * (size_t) nvox);
        normal = (float *) malloc(sizeof(float) * 3 * (size_t) nvox);
        if (!sheetness || !normal) {
            fprintf(stderr, "Memory allocation error\n");
            exit(EXIT_FAILURE);
        }

        CAT_SheetnessOptionsInit(&sopts);
        sopts.sigma_min = sheet_sigma_min;
        sopts.sigma_max = sheet_sigma_max;
        sopts.n_scales  = sheet_scales;
        sopts.polarity  = sheet_polarity;
        sopts.verbose   = verbose;

        if (CAT_VolSheetness(guide ? guide : input, sheetness, normal, NULL,
                             dims, voxelsize, &sopts) != 0) {
            fprintf(stderr, "Sheetness estimation failed.\n");
            exit(EXIT_FAILURE);
        }

        if (sheet_strength >= 0.0 && sheet_strength < 1.0)
            for (i = 0; i < nvox; i++)
                sheetness[i] *= (float) sheet_strength;
    }

    /* Apply local statistic */
    if (oriented) {
        CAT_VolOrientedMedian(input, sheetness, normal, NULL, dims, iters);
    } else if (stat_func == F_CLOSE) {
        /* Grey closing: dilation (max) followed by erosion (min) */
        localstat3(input, NULL, dims, dist, F_MAX, iters, use_euclidean_dist, DT_FLOAT32);
        localstat3(input, NULL, dims, dist, F_MIN, iters, use_euclidean_dist, DT_FLOAT32);
    } else if (stat_func == F_OPEN) {
        /* Grey opening: erosion (min) followed by dilation (max) */
        localstat3(input, NULL, dims, dist, F_MIN, iters, use_euclidean_dist, DT_FLOAT32);
        localstat3(input, NULL, dims, dist, F_MAX, iters, use_euclidean_dist, DT_FLOAT32);
    } else {
        localstat3(input, NULL, dims, dist, stat_func, iters, use_euclidean_dist, DT_FLOAT32);
    }

    /* Build output filename */
    if (argc == 3) {
        (void) sprintf(outfile, "%s", argv[2]);
    } else {
#if !defined(_WIN32) && !defined(_WIN64)
        (void) sprintf(outfile, "%s/%s_%s", dirname(infile),
                       stat_name(stat_func), basename(infile));
#else
        usage(argv[0]);
        return EXIT_FAILURE;
#endif
    }

    /* Write result */
    if (!write_nifti_float(outfile, input, DT_FLOAT32, 1.0, dims, voxelsize, nii_ptr))
        exit(EXIT_FAILURE);

    free(input);
    if (guide)
        free(guide);
    if (sheetness)
        free(sheetness);
    if (normal)
        free(normal);

    return EXIT_SUCCESS;
}
