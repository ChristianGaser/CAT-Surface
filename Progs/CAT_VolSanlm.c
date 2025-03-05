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

int is_rician = 0;

static
ArgvInfo argTable[] = {
    {"-rician", ARGV_CONSTANT, (char *) 1, (char *) &is_rician,
         "Use Rician noise estimation. MRIs can have Gaussian or Rician distributed noise with uniform or nonuniform variance across the image. If SNR is high enough (>3) noise can be well approximated by Gaussian noise in the foreground. However, for SENSE reconstruction or DTI data a Rician distribution is expected. Please note that the Rician noise estimation is sensitive for large signals in the neighbourhood and can lead to artefacts (e.g. cortex can be affected by very high values in the scalp or in blood vessels."},
     {NULL, ARGV_END, NULL, NULL, NULL}
};

void usage(char *executable)
{
    char *usage_str = "\n\
Usage: %s [options] <input.nii> [<output_denoised.nii>]\n\
\n\
    This program applies a 3D Spatial-Adaptive Non-Local Means (SA-NLM) filter to\n\
    reduce noise in medical images, particularly MRI scans. It extends the Non-Local Means\n\
    filter by incorporating spatial adaptation to improve denoising while preserving\n\
    structural details.\n\
\n\
    The filtering process includes:\n\
\n\
    1. **Noise Estimation:**\n\
       - Determines noise characteristics in the image.\n\
       - Supports both Gaussian and Rician noise models (see -rician option).\n\
\n\
    2. **Patch-Based Similarity Computation:**\n\
       - Compares intensity patches across the image to find similar structures.\n\
       - Assigns adaptive weights based on local noise variance.\n\
\n\
    3. **Denoising with Spatial Adaptation:**\n\
       - Enhances edge preservation and reduces artifacts.\n\
       - Adapts filtering strength depending on anatomical features.\n\
\n\
Options:\n\
    -rician                  Use Rician noise estimation. Recommended for MRI data \n\
                             with non-Gaussian noise distribution.\n\
\n\
Example:\n\
    %s -rician input.nii output_denoised.nii\n\n";
    
    fprintf(stderr, usage_str, executable, executable);
}

void anlm(float* ima, int v, int f, int is_rician, const int* dims);

/* Main program */

int main(int argc, char *argv[])
{
    char *infile, outfile[1024];
    int i, j, dims[3];
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
            (void) sprintf(outfile, "%s/n%s", dirname(infile), basename(infile)); 
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
    
    
    anlm(input, 3, 1, is_rician, dims);

    if (!write_nifti_float( outfile, input, DT_FLOAT32, 1.0, dims, voxelsize, nii_ptr)) 
        exit(EXIT_FAILURE);

    free(input);
    
    return(EXIT_SUCCESS);

}