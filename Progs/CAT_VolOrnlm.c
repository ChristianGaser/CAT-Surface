/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 *
 */

#include <float.h>
#include <stdlib.h>
#if !defined(_WIN32) && !defined(_WIN64)
#include <libgen.h>
#endif

#include "ParseArgv.h"
#include "CAT_NiftiLib.h"
#include "CAT_Nlm.h"

double h_ornlm = 0.05;
double sigma_ornlm = -1.0;

static ArgvInfo argTable[] = {
    {"-h", ARGV_FLOAT, (char *) 1, (char *) &h_ornlm,
         "Specifies h as a smoothing parameter controlling the decay of the exponential function."},

    {"-sigma", ARGV_FLOAT, (char *) 1, (char *) &sigma_ornlm,
         "Specifies sigma for Rician noise correction (default: equal to h)."},
         
     {NULL, ARGV_END, NULL, NULL, NULL}
};

void usage(char *executable)
{
    char *usage_str = "\n\
Usage: %s [options] <input.nii> [<output_denoised.nii>]\n\
\n\
    This program applies a 3D Optimized Blockwise Non-Local Means (ORNLM) filter to\n\
    reduce noise in medical images, particularly MRI scans.\n\
\n\
Options:\n\
    -h <float>               Smoothing parameter (default: 0.05)\n\
    -sigma <float>           Sigma for Rician noise correction. If not specified,\n\
                             it defaults to the value of h.\n\
\n\
Example:\n\
    %s -h 0.05 input.nii output_denoised.nii\n\n";
    
    fprintf(stderr, usage_str, executable, executable);
}

/* Main program */

int main(int argc, char *argv[])
{
    char *infile, outfile[1024];
    int dims[3];
    float *input;
    double voxelsize[3];
    nifti_image *nii_ptr;
    
    if (ParseArgv(&argc, argv, argTable, 0) ||(argc < 2)) {
        usage(argv[0]);
        (void) fprintf(stderr, "     %s -help\n\n", argv[0]);
         exit(EXIT_FAILURE);
    }
    
    infile  = argv[1];

    /* if not defined use original name as basename for output */
    if(argc == 3)
        (void) sprintf(outfile, "%s", argv[2]); 
    else {
        #if !defined(_WIN32) && !defined(_WIN64)
            (void) sprintf(outfile, "%s/ornlm_%s", dirname(infile), basename(infile)); 
        #else
            usage(argv[0]);
            return( 1 );
        #endif
    }
    
    /* read first image to get image parameters */
    nii_ptr = read_nifti_float(infile, &input, 0);
    if (nii_ptr == NULL) {
        fprintf(stderr,"Error reading %s.\n", infile);
        return(EXIT_FAILURE);
    }

    voxelsize[0] = nii_ptr->dx;
    voxelsize[1] = nii_ptr->dy;
    voxelsize[2] = nii_ptr->dz;
    dims[0] = nii_ptr->nx;
    dims[1] = nii_ptr->ny;
    dims[2] = nii_ptr->nz;
    
    int v = 3, f = 1;
    float noise_sigma = (float)h_ornlm;
    if (sigma_ornlm >= 0.0) noise_sigma = (float)sigma_ornlm;
    
    ornlm(input, v, f, (float)h_ornlm, noise_sigma, dims);

    if (!write_nifti_float( outfile, input, DT_FLOAT32, 1.0, dims, voxelsize, nii_ptr)) 
        exit(EXIT_FAILURE);

    free(input);
    
    return(EXIT_SUCCESS);
}
