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
    float *input_float;
    double voxelsize[3];
    nifti_image *nii_ptr;
    
    if (argc < 2) {
         (void) fprintf(stderr, "\nUsage: %s in.nii out.nii\n\n", argv[0]);
         (void) fprintf(stderr, "     %s -help\n\n", argv[0]);
     exit(EXIT_FAILURE);
    }
    
    infile  = argv[1];
    outfile = argv[2];
    
    /* read first image to get image parameters */
    nii_ptr = read_nifti_float(infile, &input_float, 0);
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
    int nvol = dims[0]*dims[1]*dims[2];

        float *vol_float = (float *)malloc(nvol * sizeof(float));
        unsigned char *mask = (unsigned char *)malloc(nvol * sizeof(unsigned char));
        for (i = 0; i < nvol; i++) {
            vol_float[i] = input_float[i];
            mask[i] = input_float[i] != 0;
        }
        
        median3(vol_float, mask, dims, 5, DT_FLOAT32);
        for (i = 0; i < nvol; i++)
            input_float[i] = (vol_float[i] > 0.5 || input_float[i] > 0.5) ? MAX(vol_float[i], input_float[i]) : MIN(vol_float[i], input_float[i]);
        free(vol_float);
        free(mask);

    if (!write_nifti_float(outfile, input_float, DT_FLOAT32, 1.0, dims, voxelsize, nii_ptr)) 
        exit(EXIT_FAILURE);

    free(input_float);
    
    return(EXIT_SUCCESS);

}