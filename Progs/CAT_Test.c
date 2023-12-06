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
#if !defined(_WIN32)
#include <libgen.h>
#endif

#include "ParseArgv.h"
#include "CAT_NiftiLib.h"
#include "CAT_Vol.h"

int rician = 0;
static
ArgvInfo argTable[] = {
        {"-rician", ARGV_CONSTANT, (char *) 1, (char *) &rician,
                 "Use Rician noise estimation. MRIs can have Gaussian or Rician distributed noise with uniform or nonuniform variance across the image. If SNR is high enough (>3) noise can be well approximated by Gaussian noise in the foreground. However, for SENSE reconstruction or DTI data a Rician distribution is expected. Please note that the Rician noise estimation is sensitive for large signals in the neighbourhood and can lead to artefacts (e.g. cortex can be affected by very high values in the scalp or in blood vessels."},
         {NULL, ARGV_END, NULL, NULL, NULL}
};

/* Main program */

int main(int argc, char *argv[])
{
        char *infile, outfile[1024];
        int i, j, dims[3];
        float *input;
        double separations[3];
        nifti_image *nii_ptr;
        
        if (ParseArgv(&argc, argv, argTable, 0) ||(argc < 2)) {
                 (void) fprintf(stderr, "\nUsage: %s [options] in.nii [out.nii]\nSpatial adaptive non-local means denoising filter.\n\n", argv[0]);
                 (void) fprintf(stderr, "         %s -help\n\n", argv[0]);
         exit(EXIT_FAILURE);
        }
        
        infile  = argv[1];

        /* if not defined use original name as basename for output */
        if(argc == 3)
                (void) sprintf(outfile, "%s", argv[2]); 
        else {
                #if !defined(_WIN32)
                        (void) sprintf(outfile, "%s/n%s", dirname(infile), basename(infile)); 
                #else
                        fprintf(stderr,"\nUsage: %s input.nii output.nii\n\n\
                            Spatial adaptive non-local means denoising filter.\n\n", argv[0]);
                        return( 1 );
                #endif
        }
        
        /* read first image to get image parameters */
        nii_ptr = read_nifti_float(infile, &input, 0);
        if(nii_ptr == NULL) {
                fprintf(stderr,"Error reading %s.\n", infile);
                return(EXIT_FAILURE);
        }

        separations[0] = nii_ptr->dx;
        separations[1] = nii_ptr->dy;
        separations[2] = nii_ptr->dz;
        dims[0] = nii_ptr->nx;
        dims[1] = nii_ptr->ny;
        dims[2] = nii_ptr->nz;
        
        
        //anlm(input, 3, 1, rician, dims);

        if (!write_nifti_float( outfile, input, DT_FLOAT32, 1.0, dims, separations, nii_ptr)) 
                exit(EXIT_FAILURE);

        free(input);
        
        return(EXIT_SUCCESS);

}