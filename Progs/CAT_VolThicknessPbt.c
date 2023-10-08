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
    char *infile, out_GMT[1024], out_RPM[1024], out_CSFD[1024], out_WMD[1024];
    int i, j, dims[3];
    float *input, *src, *dist_CSF, *dist_WM, *GMT, *RPM;
    float mean_vx_size;
    unsigned int *mask;
    double separations[3];
    nifti_image *src_ptr, *out_ptr;
    
    if (ParseArgv(&argc, argv, argTable, 0) ||(argc < 2)) {
         (void) fprintf(stderr, "\nUsage: %s [options] in.nii [out.nii]\nProjection-based thickness estimation\n", argv[0]);
         (void) fprintf(stderr, "     %s -help\n\n", argv[0]);
     exit(EXIT_FAILURE);
    }
    
    infile  = argv[1];

    /* if not defined use original name as basename for output */
    if(argc == 4) {
        (void) sprintf(out_GMT, "%s", argv[2]); 
        (void) sprintf(out_RPM, "%s", argv[3]); 
    } else {
        #if !defined(_WIN32)
            (void) sprintf(out_GMT, "%s/gmt_%s", dirname(infile), basename(infile)); 
            (void) sprintf(out_RPM, "%s/ppi_%s", dirname(infile), basename(infile)); 
        #else
            fprintf(stderr,"\nUsage: %s input.nii GMT.nii RPM.nii\n
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
    
    mask   = (unsigned int *)malloc(sizeof(unsigned int)*src_ptr->nvox);
    input  = (float *)malloc(sizeof(float)*src_ptr->nvox);
    dist_CSF = (float *)malloc(sizeof(float)*src_ptr->nvox);
    dist_WM  = (float *)malloc(sizeof(float)*src_ptr->nvox);

    GMT = (float *)malloc(sizeof(float)*src_ptr->nvox);
    RPM = (float *)malloc(sizeof(float)*src_ptr->nvox);
    
    /* check for memory faults */
    if ((input == NULL) || (mask == NULL) || (dist_CSF == NULL) ||
           (dist_WM == NULL) || (GMT == NULL) || (RPM == NULL)) {
        fprintf(stderr,"Memory allocation error\n");
        exit(EXIT_FAILURE);
    }
    
    /* prepare map outside CSF and mask to obtain distance map */
    for (i = 0; i < src_ptr->nvox; i++) {
        input[i] = (src[i] < 1.5) ? 1.0f : 0.0f;
        mask[i]  = (src[i] < 3.0) ? 1   : 0;
        mask[i] = 1;
    }    

    /* obtain CSF distance map */
    vbdist(input, mask, dims, separations);
    for (i = 0; i < src_ptr->nvox; i++)
        dist_CSF[i] = input[i] + 0.0;

    out_ptr = nifti_copy_nim_info(src_ptr);
    (void) sprintf(out_CSFD, "%s/dist_CSF_%s", dirname(infile), basename(infile)); 
    if (!write_nifti_float(out_CSFD, dist_CSF, DT_FLOAT32, 1.0, dims, separations, out_ptr)) 
        exit(EXIT_FAILURE);
            
    /* prepare map outside WM and mask to obtain distance map */
    for (i = 0; i < src_ptr->nvox; i++) {
        input[i] = (src[i] > 2.5) ? 1.0f : 0.0f;
        mask[i]  = (src[i] > 1.0) ? 1   : 0;
        mask[i] = 1;
    }    

    /* obtain WM distance map */
    vbdist(input, mask, dims, separations);
    for (i = 0; i < src_ptr->nvox; i++)
        dist_WM[i] = input[i] + 0.0;

    (void) sprintf(out_WMD, "%s/dist_WM_%s", dirname(infile), basename(infile)); 
    if (!write_nifti_float(out_WMD, dist_WM, DT_FLOAT32, 1.0, dims, separations, out_ptr)) 
        exit(EXIT_FAILURE);
    
    projection_based_thickness(src, dist_WM, dist_CSF, GMT, RPM, dims, separations); 

    /* Because dist_CSF and dist_WM measure a grid-based distance (defined as the 
       center of a voxel), we have to correct by 1 voxel, and finally correct
       for the isotropic size of our voxel-grid.
    */
    mean_vx_size = (separations[0]+separations[1]+separations[2])/3.0;
    for (i = 0; i < src_ptr->nvox; i++)
        GMT[i] = (GMT[i] - 1.0)*mean_vx_size;
       
    if (!write_nifti_float(out_GMT, GMT, DT_FLOAT32, 1.0, dims, separations, out_ptr)) 
        exit(EXIT_FAILURE);
    if (!write_nifti_float(out_RPM, RPM, DT_FLOAT32, 1.0, dims, separations, out_ptr)) 
        exit(EXIT_FAILURE);

    free(mask);
    free(input);
    free(dist_CSF);
    free(dist_WM);
    free(GMT);
    free(RPM);    
    
    return(EXIT_SUCCESS);

}