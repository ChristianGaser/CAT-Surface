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
double sheet_cutoff       = 0.0;   /* 0 selects CAT_ORIENTED_MEDIAN_CUTOFF */
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
         "Overall gain on the sheetness before it is used (default: 1.0).\n\
         0 reproduces the isotropic filter exactly.  Values above 1 amplify a\n\
         response too weak to matter: the oriented median admits every neighbour\n\
         unless the sheetness exceeds 0.5, so a map peaking at 0.5 leaves it\n\
         bit-identical to the isotropic median.  Reaching 0.5 is only where it\n\
         starts to differ, though: a one-voxel sheet has 9 of 27 neighbours in\n\
         its own plane, and the 8 edge neighbours that would outvote them drop\n\
         out only at a sheetness of 1, so the median returns the sheet value\n\
         only once the response saturates.  Check the map with CAT_VolSheetness\n\
         and set this to about 1/max, not 0.6/max.  The noise floor is amplified\n\
         along with the sheets, and the strongest responses -- not necessarily\n\
         the sulci -- come up first, so watch the map rather than guessing.\n\
         Lowering -c in CAT_VolSheetness is the principled alternative."},
    {"-sheet-cutoff", ARGV_FLOAT, (char *) 1, (char *) &sheet_cutoff,
         "Admission cutoff of the oriented median (default: 0.10).  A neighbour\n\
         at offset d is admitted when s*(dhat.n)^2 < cutoff, so the three\n\
         neighbour classes drop out in turn: the 6 face neighbours at s = cutoff,\n\
         the 12 edge neighbours at s = 2*cutoff, the 8 corners at s = 3*cutoff.\n\
         The 9 offsets lying in the sheet plane are always admitted.  A\n\
         one-voxel-thick sheet is preserved from s = 2*cutoff upwards, where\n\
         those 9 first make up half the admitted set, so the cutoff is half the\n\
         sheetness at which the filter starts to protect a thin structure.  Set\n\
         it from the response CAT_VolSheetness actually gives on your data: a\n\
         map peaking near 0.5 wants roughly 0.25, one peaking near 0.3 wants\n\
         0.15.  Larger values confine the effect to a thinner rim around the\n\
         strongest responses; 0.5 was the historical value and needed a\n\
         saturated response before it did anything at all."},
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
    -sheet-strength <float>  Gain on the sheetness (default: 1.0); see below.\n\
    -sheet-cutoff   <float>  Admission cutoff (default: 0.10); see below.\n\
    -v               Be verbose.\n\
\n\
    With -oriented the median runs over a neighbourhood oriented by a Hessian\n\
    sheetness filter (see CAT_VolSheetness): a neighbour at offset d is admitted\n\
    only when s*(dhat.n)^2 < cutoff, with s the local sheetness, n the local\n\
    sheet normal and cutoff set by -sheet-cutoff.  At s = 0 every neighbour is\n\
    admitted and the operator is the plain isotropic median.\n\
\n\
    The 9 offsets lying in the sheet plane are always admitted; the other three\n\
    classes drop out at s = cutoff (6 faces), s = 2*cutoff (12 edges) and\n\
    s = 3*cutoff (8 corners).  A one-voxel-thick sheet puts the sheet value on\n\
    the 9 in-plane offsets only, so it survives the median from s = 2*cutoff\n\
    upwards, where those 9 first make up half the admitted set.\n\
\n\
    So the two numbers to match are the cutoff and the response your data\n\
    actually produces: write the map with CAT_VolSheetness, and pick a cutoff of\n\
    about half the sheetness your sulci reach.  If the effect is confined to a\n\
    thin rim around the strongest responses, the cutoff is too high for the map;\n\
    lower it, or raise -sheet-strength, or widen the scale range.\n\
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
        sopts.gain      = sheet_strength;
        sopts.verbose   = verbose;

        if (CAT_VolSheetness(guide ? guide : input, sheetness, normal, NULL,
                             dims, voxelsize, &sopts) != 0) {
            fprintf(stderr, "Sheetness estimation failed.\n");
            exit(EXIT_FAILURE);
        }
    }

    /* Apply local statistic */
    if (oriented) {
        CAT_VolOrientedMedian(input, sheetness, normal, NULL, dims, sheet_cutoff,
                              iters);
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
