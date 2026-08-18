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

/* Defaults mirror CAT_SheetnessOptionsInit() */
double sigma_min = 0.3;
double sigma_max = 1.0;
int    n_scales  = 3;
double alpha     = 0.5;
double beta      = 0.5;
double c_noise   = -1.0;
double strength  = 1.0;
int    polarity  = 0;
int    verbose   = 0;

static
ArgvInfo argTable[] = {
    {"-sigma-min", ARGV_FLOAT, (char *) 1, (char *) &sigma_min,
         "Smallest Gaussian scale in mm (default: 0.3)."},
    {"-sigma-max", ARGV_FLOAT, (char *) 1, (char *) &sigma_max,
         "Largest Gaussian scale in mm (default: 1.0).  The range should\n\
         bracket the structures being looked for.  A sulcal CSF sheet or a\n\
         gyral WM blade that survives at all is one to three voxels thick at\n\
         0.5 mm; larger scales start responding to the cortical ribbon itself."},
    {"-scales", ARGV_INT, (char *) 1, (char *) &n_scales,
         "Number of log-spaced scales between the two sigmas (default: 3)."},
    {"-polarity", ARGV_INT, (char *) 1, (char *) &polarity,
         "Which sign of the dominant curvature to accept: 1 = bright sheets\n\
         (ridges, e.g. gyral WM blades), -1 = dark sheets (valleys, e.g.\n\
         sulcal CSF), 0 = either (default).  The guard is what stops a\n\
         WM-blade detector from responding to sulci and vice versa."},
    {"-alpha", ARGV_FLOAT, (char *) 1, (char *) &alpha,
         "Plate-vs-tube sensitivity (default: 0.5).  Larger values accept\n\
         more tube-like structures as sheets."},
    {"-beta", ARGV_FLOAT, (char *) 1, (char *) &beta,
         "Blob-vs-plate sensitivity (default: 0.5)."},
    {"-c", ARGV_FLOAT, (char *) 1, (char *) &c_noise,
         "Structure-vs-noise sensitivity.  Negative (default) selects the\n\
         automatic estimate, which makes the filter independent of the\n\
         intensity units of the input.  The automatic value is half the largest\n\
         Hessian norm in the volume, so on a whole head it is set by the\n\
         scalp/air step and the cortical response collapses; lowering it is the\n\
         principled way to bring the response back up."},
    {"-strength", ARGV_FLOAT, (char *) 1, (char *) &strength,
         "Overall gain on the response (default: 1.0).  The blunt alternative to\n\
         -c when the intensity units are unknown: the map is multiplied by this\n\
         and clamped to [0,1].  Values above 1 amplify a response that is too\n\
         weak to reach the thresholds the consumers gate on -- notably the hard\n\
         0.5 of the oriented median, below which it is exactly the isotropic\n\
         median.  0 reproduces the isotropic filters exactly.  The same knob is\n\
         called -sheet-strength in the tools that consume the field."},
    {"-v", ARGV_CONSTANT, (char *) 1, (char *) &verbose,
         "Be verbose."},
    {NULL, ARGV_END, NULL, NULL, NULL}
};

static void usage(char *executable)
{
    char *usage_str = "\n\
Usage: %s [options] <input.nii> [<output.nii>]\n\
\n\
    Multi-scale Hessian sheetness (plate) filter.\n\
\n\
    Isotropic regularizers -- a local median, a Potts MRF, total variation --\n\
    penalize boundary area.  A thin structure has an extreme area-to-volume\n\
    ratio, so deleting it is always the cheaper labelling.  That shrinking bias\n\
    is why one and the same median filter opens glued sulci in one place and\n\
    closes cerebellar fissures in another: to an isotropic prior the two are\n\
    the same object and only the sign of the label differs.\n\
\n\
    A sheetness filter escapes this by using a *shape* prior instead of a\n\
    smoothness prior.  From the eigenvalues of the Hessian, |l1| <= |l2| <=\n\
    |l3|, the local shape is read off directly:\n\
\n\
        sheet (plate):  |l1| ~ |l2| ~ 0,  |l3| large\n\
        tube  (line):   |l1| ~ 0,         |l2| ~ |l3| large\n\
        blob:           |l1| ~ |l2| ~ |l3| large\n\
\n\
    It therefore keeps thin sheets and ignores blobs, and it shrinks nothing,\n\
    because it makes no statement about boundary length.  Sulcal CSF is a dark\n\
    sheet in a T1 (-polarity -1), a gyral white-matter blade is a bright sheet\n\
    (-polarity 1); the same operator finds both.\n\
\n\
    The response is written as a float map in [0,1].  Its main use is to\n\
    check the scale range and the polarity on a new protocol before enabling\n\
    any of the options that consume it -- CAT_VolLocalStat -oriented,\n\
    CAT_VolThicknessPbt -oriented-filter and CAT_VolSulcusRepair, all of\n\
    which estimate the field\n\
    themselves rather than reading it from a file.  The accompanying sheet\n\
    normals are not written: a 3-vector per voxel needs a 4-D image, which\n\
    write_nifti_float() does not produce, and recomputing the field is cheap.\n\
\n\
Options:\n\
    -sigma-min <float>  Smallest scale in mm (default: 0.3).\n\
    -sigma-max <float>  Largest scale in mm (default: 1.0).\n\
    -scales    <int>    Number of log-spaced scales (default: 3).\n\
    -polarity  <int>    1 = bright sheets, -1 = dark sheets, 0 = either.\n\
    -alpha     <float>  Plate-vs-tube sensitivity (default: 0.5).\n\
    -beta      <float>  Blob-vs-plate sensitivity (default: 0.5).\n\
    -c         <float>  Noise sensitivity; negative selects automatic.\n\
    -strength  <float>  Overall gain on the response (default: 1.0).\n\
    -v                  Be verbose.\n\
\n\
Example:\n\
    %s -polarity -1 -v t1_corr.nii sheetness.nii\n\
    %s -stat 7 -oriented -guide t1_corr.nii label.nii label_filtered.nii\n\
\n\
References:\n\
    Sato et al., Med Image Anal 2(2):143-168, 1998.\n\
    Frangi et al., MICCAI 1998, LNCS 1496:130-137.\n\
    Descoteaux et al., Med Image Anal 10(4):638-651, 2006.\n\n";

    fprintf(stderr, usage_str, executable, executable, "CAT_VolLocalStat");
}

/* Main program */
int
main(int argc, char *argv[])
{
    char *infile, outfile[1024];
    int dims[3], nvox, rc;
    float *input, *sheetness;
    double voxelsize[3];
    nifti_image *nii_ptr;
    CAT_SheetnessOpts opts;

    if (ParseArgv(&argc, argv, argTable, 0) || (argc < 2)) {
        usage(argv[0]);
        (void) fprintf(stderr, "     %s -help\n\n", argv[0]);
        exit(EXIT_FAILURE);
    }

    infile = argv[1];

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
    nvox = dims[0] * dims[1] * dims[2];

    sheetness = (float *) malloc(sizeof(float) * (size_t) nvox);
    if (!sheetness) {
        fprintf(stderr, "Memory allocation error\n");
        exit(EXIT_FAILURE);
    }

    CAT_SheetnessOptionsInit(&opts);
    opts.sigma_min = sigma_min;
    opts.sigma_max = sigma_max;
    opts.n_scales  = n_scales;
    opts.alpha     = alpha;
    opts.beta      = beta;
    opts.c         = c_noise;
    opts.gain      = strength;
    opts.polarity  = polarity;
    opts.verbose   = verbose;

    if (verbose)
        fprintf(stdout, "Sheetness of %s: %d scales in [%g, %g] mm, polarity %d, gain %g.\n",
                infile, n_scales, sigma_min, sigma_max, polarity, strength);

    rc = CAT_VolSheetness(input, sheetness, NULL, NULL, dims, voxelsize, &opts);
    if (rc != 0) {
        fprintf(stderr, "Sheetness estimation failed (%d).\n", rc);
        exit(EXIT_FAILURE);
    }

    if (argc == 3) {
        (void) sprintf(outfile, "%s", argv[2]);
    } else {
#if !defined(_WIN32) && !defined(_WIN64)
        (void) sprintf(outfile, "%s/sheet_%s", dirname(infile), basename(infile));
#else
        usage(argv[0]);
        return EXIT_FAILURE;
#endif
    }

    if (!write_nifti_float(outfile, sheetness, DT_FLOAT32, 1.0, dims, voxelsize, nii_ptr))
        exit(EXIT_FAILURE);

    free(input);
    free(sheetness);

    return EXIT_SUCCESS;
}
