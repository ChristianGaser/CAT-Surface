/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 *
 */

#include <float.h>
#include <stdlib.h>

#if !defined(_WIN32) && !defined(_WIN64)
#include <libgen.h>
#endif

#include "ParseArgv.h"
#include "CAT_NiftiLib.h"
#include "CAT_Vol.h"
#include "CAT_Sheetness.h"

double fwhm = 8.0;
int use_mask = 0;
int oriented = 0;
int iters = 1;
double sheet_sigma = 0.5;
char *guide_file = NULL;
double sheet_sigma_min = 0.3;
double sheet_sigma_max = 1.0;
int sheet_scales = 3;
int sheet_polarity = 0;
double sheet_strength = 1.0;
int verbose = 0;

static
ArgvInfo argTable[] = {
    {"-fwhm", ARGV_FLOAT, (char *) 1, (char *) &fwhm, 
             "FWHM in mm."},
    {"-mask", ARGV_CONSTANT, (char *) 1, (char *) &use_mask,
             "Use masked smoothing (default no masking)."},
    {"-oriented", ARGV_CONSTANT, (char *) 1, (char *) &oriented,
             "Diffuse along thin sheets instead of across them (coherence-enhancing\n\
             smoothing).  Each 3x3x3 neighbour is weighted by\n\
             (1-s) + s*exp(-(dhat.n)^2 / 2 sigma^2), with s the local sheetness\n\
             from a Hessian plate filter and n the local sheet normal, so\n\
             tangential neighbours keep full weight while neighbours across the\n\
             sheet are suppressed.  At s = 0 the weights reduce to plain inverse\n\
             distance and the operator is ordinary local averaging.  This is an\n\
             iterated local kernel, not a Gaussian: -fwhm does not apply and the\n\
             amount of smoothing is set by -iter."},
    {"-iter", ARGV_INT, (char *) 1, (char *) &iters,
             "Number of passes for -oriented (default: 1)."},
    {"-sheet-sigma", ARGV_FLOAT, (char *) 1, (char *) &sheet_sigma,
             "Angular width of the anisotropy for -oriented (default: 0.5).\n\
             Smaller values confine the diffusion more tightly to the sheet plane."},
    {"-guide", ARGV_STRING, (char *) 1, (char *) &guide_file,
             "Volume the orientation is estimated from (default: the input itself).\n\
             Pass the intensity image when smoothing a label or probability map.\n\
             Must have the same dimensions as the input."},
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
             0 reproduces isotropic local averaging exactly."},
    {"-v", ARGV_CONSTANT, (char *) 1, (char *) &verbose,
             "Be verbose."},
     {NULL, ARGV_END, NULL, NULL, NULL}
};


/* Main program */
int
main(int argc, char *argv[])
{
    char *infile, outfile[1024];
    int i, dims[3];
    float *input;
    double voxelsize[3], s[3];
    nifti_image *nii_ptr;
    float *guide = NULL, *sheetness = NULL, *normal = NULL;
    nifti_image *guide_ptr = NULL;

    /* Get arguments */
#if !defined(_WIN32) && !defined(_WIN64)
    if (ParseArgv(&argc, argv, argTable, 0) ||(argc < 2)) {
#else
    if (ParseArgv(&argc, argv, argTable, 0) ||(argc < 3)) {
#endif
        (void) fprintf(stderr, 
        "\nUsage: %s [options] in.nii [out.nii]\n", argv[0]);
        (void) fprintf(stderr, 
        "     %s -help\n\n", argv[0]);
        exit(EXIT_FAILURE);
    }
    
    infile = argv[1];
    
    if (verbose) {
        if (oriented)
            fprintf(stdout,"Filtering %s with %d oriented pass(es).\n", infile, iters);
        else
            fprintf(stdout,"Filtering %s with FWHM of %gmm.\n", infile, fwhm);
    }
    
    /* read first image to get image parameters */
    nii_ptr = read_nifti_float(infile, &input, 0);
    if (nii_ptr == NULL) {
        fprintf(stderr,"Error reading %s.\n", infile);
        return(EXIT_FAILURE);
    }

    /* only allow isotropic filtering */
    for (i=0; i<3; i++) s[i] = fwhm;
    
    voxelsize[0] = nii_ptr->dx;
    voxelsize[1] = nii_ptr->dy;
    voxelsize[2] = nii_ptr->dz;
    dims[0] = nii_ptr->nx;
    dims[1] = nii_ptr->ny;
    dims[2] = nii_ptr->nz;
    
    if (oriented) {
        /* The orientation field is recomputed here rather than read from a
           file: a 3-vector per voxel needs a 4-D image, and the Hessian pass
           is cheap next to the cost of moving such a field through disk. */
        int nvox = dims[0] * dims[1] * dims[2];
        CAT_SheetnessOpts sopts;

        if (guide_file) {
            guide_ptr = read_nifti_float(guide_file, &guide, 0);
            if (guide_ptr == NULL) {
                fprintf(stderr,"Error reading %s.\n", guide_file);
                return(EXIT_FAILURE);
            }
            if (guide_ptr->nx != dims[0] || guide_ptr->ny != dims[1] ||
                guide_ptr->nz != dims[2]) {
                fprintf(stderr,"Guide volume must have the same dimensions as the input.\n");
                return(EXIT_FAILURE);
            }
        }

        sheetness = (float *) malloc(sizeof(float) * (size_t) nvox);
        normal = (float *) malloc(sizeof(float) * 3 * (size_t) nvox);
        if (!sheetness || !normal) {
            fprintf(stderr,"Memory allocation error\n");
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
            fprintf(stderr,"Sheetness estimation failed.\n");
            exit(EXIT_FAILURE);
        }

        if (sheet_strength >= 0.0 && sheet_strength < 1.0)
            for (i = 0; i < nvox; i++)
                sheetness[i] *= (float) sheet_strength;

        CAT_VolOrientedSmooth(input, sheetness, normal, NULL, dims,
                              sheet_sigma, iters);
    } else
        smooth3(input, dims, voxelsize, s, use_mask, DT_FLOAT32);
    
    /* if not defined use original name as basename for output */
    if (argc == 3)
            (void) sprintf(outfile, "%s", argv[2]); 
    else {
        #if !defined(_WIN32) && !defined(_WIN64)
            if (oriented)
                (void) sprintf(outfile, "%s/os%d%s", dirname(infile), iters, basename(infile)); 
            else
                (void) sprintf(outfile, "%s/s%g%s", dirname(infile), fwhm, basename(infile)); 
        #endif
    }

    /* write data using same data type and rescale */
    if (!write_nifti_float(outfile, input, nii_ptr->datatype, 0.0, dims, voxelsize, nii_ptr)) 
        exit(EXIT_FAILURE);

    if (guide)
        free(guide);
    if (sheetness)
        free(sheetness);
    if (normal)
        free(normal);
    
    return(EXIT_SUCCESS);

}
