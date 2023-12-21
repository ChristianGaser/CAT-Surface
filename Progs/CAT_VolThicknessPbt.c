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
int no_minimum_thickness = 0;
int n_avgs = 4;

extern int ParseArgv(int *argcPtr, char **argv, ArgvInfo *argTable, int flags);

static ArgvInfo argTable[] = {
    {"-v", ARGV_CONSTANT, (char *) 1, (char *) &verbose,
          "Be verbose."},
    {"-n_avgs", ARGV_INT, (char *) 1, (char *) &n_avgs,
          "Number of averages for estimating WM and CSF distance by shifting the border between GM/WM and GM/CSF to obtain a less noisy measure."},
    {"-no_min_thickness", ARGV_CONSTANT, (char *) 1, (char *) &no_minimum_thickness,
          "Don't use two thickness measures from sulci and gyri and estimate minimum, but use the simpler approach based on sulci only."},
     {NULL, ARGV_END, NULL, NULL, NULL}
};

int main(int argc, char *argv[])
{
    char *infile, out_GMT[1024], out_PPM[1024];
    int i, j, dims[3];
    float *input, *src, *dist_CSF, *dist_WM, *GMT, *GMT2, *PPM, *PPM_filtered;
    float mean_vx_size;
    unsigned int *mask;
    double separations[3], add_value;
    nifti_image *src_ptr, *out_ptr;
    
    if (ParseArgv(&argc, argv, argTable, 0) ||(argc < 2)) {
         (void) fprintf(stderr, "\nUsage: %s [options] in.nii [GMT.nii PPM.nii]\n", argv[0]);
         (void) fprintf(stderr, "         Projection-based thickness estimation based on PVE label image, where background\n");
         (void) fprintf(stderr, "         has value 1, gray matter 2, and white matter 3 and there are two additional PVE classes (CGM/GWM).\n");
         (void) fprintf(stderr, "         A voxel-wise thickness and percentage position map are saved.\n");
         (void) fprintf(stderr, "         If no output names are defined the input name will be prepended by gmt_ and ppm_\n");
         (void) fprintf(stderr, "     %s -help\n\n", argv[0]);
     exit(EXIT_FAILURE);
    }
    
    infile  = argv[1];

    /* if not defined use original name as basename for output */
    if(argc == 4) {
        (void) sprintf(out_GMT, "%s", argv[2]); 
        (void) sprintf(out_PPM, "%s", argv[3]); 
    } else {
        #if !defined(_WIN32)
            (void) sprintf(out_GMT, "%s/gmt_%s", dirname(infile), basename(infile)); 
            (void) sprintf(out_PPM, "%s/ppm_%s", dirname(infile), basename(infile)); 
        #else
            fprintf(stderr,"\nUsage: %s input.nii GMT.nii PPM.nii\n\n", argv[0]);
            return( 1 );
        #endif
    }
    
    /* read source image */
    src_ptr = read_nifti_float(infile, &src, 0);
    if(src_ptr == NULL) {
        fprintf(stderr,"Error reading %s.\n", infile);
        return(EXIT_FAILURE);
    }

    out_ptr = nifti_copy_nim_info(src_ptr);

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

    GMT  = (float *)malloc(sizeof(float)*src_ptr->nvox);
    PPM  = (float *)malloc(sizeof(float)*src_ptr->nvox);
    PPM_filtered = (float *)malloc(sizeof(float)*src_ptr->nvox);
    
    /* check for memory faults */
    if ((input == NULL) || (mask == NULL) || (dist_CSF == NULL) ||
           (dist_WM == NULL) || (GMT == NULL) || (PPM == NULL) || (PPM_filtered == NULL)) {
        fprintf(stderr,"Memory allocation error\n");
        exit(EXIT_FAILURE);
    }
    
    /* initialize distances */
    for (i = 0; i < src_ptr->nvox; i++) {
        dist_CSF[i] = 0.0;
        dist_WM[i]  = 0.0;
    }
    
    for (j = 0; j < n_avgs; j++) {
        /* estimate value for shifting the border to obtain a less noisy measure by averaging distances */
        add_value = ((double)j + 1.0) / ((double)n_avgs + 1.0) - 0.5;
        /* prepare map outside CSF and mask to obtain distance map for CSF */
        for (i = 0; i < src_ptr->nvox; i++) {
            input[i] = (src[i] < (CGM+add_value)) ? 1.0f : 0.0f;
            mask[i]  = (src[i] < WM) ? 1 : 0;
        }    
    
        /* obtain CSF distance map */
        if (verbose) fprintf(stderr,"Estimate CSF distance map.\n");
        vbdist(input, mask, dims, separations);
        for (i = 0; i < src_ptr->nvox; i++)
            dist_CSF[i] += input[i];
                
        /* prepare map outside WM and mask to obtain distance map for WN */
        for (i = 0; i < src_ptr->nvox; i++) {
            input[i] = (src[i] > (GWM+add_value)) ? 1.0f : 0.0f;
            mask[i]  = (src[i] > CSF) ? 1 : 0;
        }    
    
        /* obtain WM distance map */
        if (verbose) fprintf(stderr,"Estimate WM distance map.\n");
        vbdist(input, mask, dims, separations);
        for (i = 0; i < src_ptr->nvox; i++)
            dist_WM[i] += input[i];
    }
    
    /* estimate average */
    if (n_avgs > 1) {
        for (i = 0; i < src_ptr->nvox; i++) {
            dist_CSF[i] /= (float) n_avgs;
            dist_WM[i]  /= (float) n_avgs;
        }
    }

    if (verbose) fprintf(stderr,"Estimate thickness map.\n");
    /* first reconstruct sulci */
    projection_based_thickness(src, dist_WM, dist_CSF, GMT, dims, separations); 

    /* only use reconstruction of sulci */
    if (no_minimum_thickness) {
        /* use minimum to reduce issues with meninges */
        for (i = 0; i < src_ptr->nvox; i++)
            GMT[i]  = MIN(GMT[i],  dist_WM[i]+dist_CSF[i]);
    } else { /* use both reconstruction of sulci as well as gyri and use minimum of both */
      
        /* we need the inverse of src: 4 - src */
        for (i = 0; i < src_ptr->nvox; i++)
            input[i] = 4.0 - src[i];
            
        GMT2 = (float *)malloc(sizeof(float)*src_ptr->nvox);
        if (GMT2 == NULL) {
            fprintf(stderr,"Memory allocation error\n");
            exit(EXIT_FAILURE);
        }
        
        /* then reconstruct gyri by using the inverse of src and switching the WM and CSF distance */
        projection_based_thickness(input, dist_CSF, dist_WM, GMT2, dims, separations); 
    
        /* use minimum for each measure to reduce issues with meninges */
        for (i = 0; i < src_ptr->nvox; i++) {
            GMT[i]  = MIN(GMT[i],  dist_WM[i]+dist_CSF[i]);
            GMT2[i] = MIN(GMT2[i], dist_WM[i]+dist_CSF[i]);
        }
        
        /* finally use minimum of both thickness measures */
        for (i = 0; i < src_ptr->nvox; i++)
            GMT[i] = MIN(GMT[i], GMT2[i]);
            
        free(GMT2);
    }
        
    /* use masked smoothing for thickness map
       this should be probably replaced by simple_approx */
    if (verbose) fprintf(stderr,"Correct thickness map.\n");
    double s[] = {1.0, 1.0, 1.0};
    //smooth_float(GMT, dims, separations, s, 1);

    /* init PPM */
    for (i = 0; i < src_ptr->nvox; i++)
        PPM[i] = (src[i] >= GWM) ? 1.0f : 0.0f;

    /* Estimate percentage position map (PPM)
       We first create a corrected CSF distance map with reconstructed sulci.
       If gyri were reconstructed too than also the dist_WM have to be
       corrected to avoid underestimation of the position map with surfaces 
       running to close to the WM. */
    if (verbose) fprintf(stderr,"Correct percentage position map.\n");
    for (i = 0; i < src_ptr->nvox; i++) {
        if ((src[i] >= (CGM+add_value)) && (src[i] < (GWM+add_value)) && GMT[i] > 1e-15) {
            PPM[i] = MIN(dist_CSF[i], GMT[i]-dist_WM[i]) / GMT[i];
        }
        if (PPM[i] < 0.0) PPM[i] = 0.0;
    }
    
    /* finally minimize outliers in the PPM using median-filter */
    for (i = 0; i < src_ptr->nvox; i++)
        PPM_filtered[i] = PPM[i];
    median3_float(PPM_filtered, dims);
    
    /* protect values in sulci and only replace other areas with median-filtered values */
    for (i = 0; i < src_ptr->nvox; i++)
        if (PPM[i] > 0.25)
            PPM[i] = PPM_filtered[i];

    /* we have to correct for the isotropic size of our voxel-grid */
    mean_vx_size = (separations[0]+separations[1]+separations[2])/3.0;
    for (i = 0; i < src_ptr->nvox; i++) 
        GMT[i] *= mean_vx_size;
    
    /* apply (masked) smoothing */
    if (verbose) fprintf(stderr,"Final correction\n");
    s[0] = s[1] = s[2] = 0.9;
    //smooth_float(GMT, dims, separations, s, 0);
    //smooth_float(PPM, dims, separations, s, 0);
    
    /* save GMT and PPM image */
    if (!write_nifti_float(out_GMT, GMT, DT_FLOAT32, 1.0, dims, separations, out_ptr)) 
        exit(EXIT_FAILURE);
    if (!write_nifti_float(out_PPM, PPM, DT_FLOAT32, 1.0, dims, separations, out_ptr)) 
        exit(EXIT_FAILURE);

    free(dist_CSF);
    free(dist_WM);
    free(GMT);
    free(PPM);    
    free(PPM_filtered);    
    free(mask);
    free(input);

    return(EXIT_SUCCESS);

}