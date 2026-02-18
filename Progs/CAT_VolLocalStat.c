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

/* Defaults */
int    stat_func          = F_MEAN;
int    dist               = 1;
int    iters              = 1;
int    use_euclidean_dist = 0;
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
    -v               Be verbose.\n\
\n\
Example:\n\
    %s -stat 7 -dist 2 -iter 3 input.nii output.nii\n\n";

    fprintf(stderr, usage_str, executable, executable);
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
    case F_MEDIAN: return "median";
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

    /* Parse arguments */
    if (ParseArgv(&argc, argv, argTable, 0) || (argc < 2)) {
        usage(argv[0]);
        (void) fprintf(stderr, "     %s -help\n\n", argv[0]);
        exit(EXIT_FAILURE);
    }

    infile = argv[1];

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

    /* Apply local statistic */
    if (stat_func == F_CLOSE) {
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

    return EXIT_SUCCESS;
}
