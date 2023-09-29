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

int verbose = 0;

extern int ParseArgv(int *argcPtr, char **argv, ArgvInfo *argTable, int flags);

static ArgvInfo argTable[] = {
        {"-v", ARGV_CONSTANT, (char *) 1, (char *) &verbose,
                  "Be verbose."},
         {NULL, ARGV_END, NULL, NULL, NULL}
};

int main(int argc, char *argv[])
{
        char *infile, outfile[1024];
        int i, j, dims[3];
        float F, *input, *src, *dist_CSF, *dist_WM;
        unsigned int *mask;
        double separations[3];
        nifti_image *src_ptr;
        
        if (ParseArgv(&argc, argv, argTable, 0) ||(argc < 2)) {
                 (void) fprintf(stderr, "\nUsage: %s [options] in.nii [out.nii]\nProjection-based thickness estimation\n", argv[0]);
                 (void) fprintf(stderr, "         %s -help\n\n", argv[0]);
         exit(EXIT_FAILURE);
        }
        
        infile  = argv[1];

        /* if not defined use original name as basename for output */
        if(argc == 3)
                (void) sprintf(outfile, "%s", argv[2]); 
        else {
                #if !defined(_WIN32)
                        (void) sprintf(outfile, "%s/ppi_%s", dirname(infile), basename(infile)); 
                #else
                        fprintf(stderr,"\nUsage: %s input.nii output.nii\n
                Projection-based thickness estimation.\n\n", argv[0]);
                        return( 1 );
                #endif
        }
        
        /* read source image */
        src_ptr = read_nifti_float(infile, &src, 0);
        if(src_ptr == NULL) {
                fprintf(stderr,"Error reading %s.\n", infile);
                return(EXIT_FAILURE);
        }

        /* get dimensions and voxel size */
        separations[0] = src_ptr->dx;
        separations[1] = src_ptr->dy;
        separations[2] = src_ptr->dz;
        dims[0] = src_ptr->nx;
        dims[1] = src_ptr->ny;
        dims[2] = src_ptr->nz;
        
        mask     = (unsigned int *)malloc(sizeof(unsigned int)*src_ptr->nvox);
        input    = (float *)malloc(sizeof(float)*src_ptr->nvox);
        dist_CSF = (float *)malloc(sizeof(float)*src_ptr->nvox);
        dist_WM  = (float *)malloc(sizeof(float)*src_ptr->nvox);
        
        /* check for memory faults */
        if ((input == NULL) || (mask == NULL) || (dist_CSF == NULL) || (dist_WM == NULL)) {
                fprintf(stderr,"Memory allocation error\n");
                exit(EXIT_FAILURE);
        }
        
        /* prepare map outside CSF and mask to obtain distance map */
        for (i = 0; i < src_ptr->nvox; i++) {
                input[i] = (src[i] > 1.5) ? 1.0 : 0.0;
                mask[i]  = (src[i] < 3.0) ? 1 : 0;
        }        

        /* obtain CSF distance map */
        vbdist(input, mask, dims, separations);
        memcpy(dist_CSF, input, sizeof(double)*src_ptr->nvox);
        
        /* prepare map outside WM and mask to obtain distance map */
        for (i = 0; i < src_ptr->nvox; i++) {
                input[i] = (src[i] < 2.5) ? 1.0 : 0.0;
                mask[i]  = (src[i] > 1.0) ? 1 : 0;
        }        

        /* obtain WM distance map */
        vbdist(input, mask, dims, separations);
        memcpy(dist_WM, input, sizeof(double)*src_ptr->nvox);
        

        /* 

    % estimate thickness with PBT approach
    % CSF and WM distance maps (based on a binary boundary without PVE)
    dist_CSF = cat_vbdist(single(Yp0 > 1.5), Yp0 < 3); 
    dist_WM = cat_vbdist(single(Yp0 < 2.5), Yp0 > 1); 

    % projection-based thickness mapping
    [Ygmt,Ypp] = cat_vol_pbtp(round(Yp0), dist_WM, dist_CSF); 

    % Because dist_CSF and dist_WM measure a grid-based distance (defined as the 
    % center of a voxel), we have to correct 1 voxel, and finally correct
    % for the _isotropic_ size of our voxel-grid.
    Ygmt = (Ygmt - 1) / mean(vx_vol); 
            
        */        
        
        if (!write_nifti_float(outfile, dist_WM, DT_FLOAT32, 1.0, dims, separations, src_ptr)) 
                exit(EXIT_FAILURE);

        free(mask);
        free(input);
        free(dist_CSF);
        free(dist_WM);
        
        return(EXIT_SUCCESS);

}