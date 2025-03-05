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

#include "ParseArgv.h"
#include "CAT_NiftiLib.h"
#include "CAT_Vol.h"


/* Main program */
int main(int argc, char *argv[])
{
    char *infile, *outfile;
    int i, j, dims[3];
    float *input;
    double voxelsize[3];
    nifti_image *nii_ptr;
    
    if (argc < 2) {
         (void) fprintf(stderr, "\nUsage: %s in.nii out.nii]\n\n", argv[0]);
         (void) fprintf(stderr, "     %s -help\n\n", argv[0]);
     exit(EXIT_FAILURE);
    }
    
    infile  = argv[1];
    outfile = argv[2];
    
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

    vol_approx(input, dims, voxelsize);

    if (!write_nifti_float(outfile, input, DT_FLOAT32, 1.0, dims, voxelsize, nii_ptr)) 
        exit(EXIT_FAILURE);

    free(input);
    
    return(EXIT_SUCCESS);

}