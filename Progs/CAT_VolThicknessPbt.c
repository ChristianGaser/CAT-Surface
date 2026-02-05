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

#include <bicpl.h>
#include <ParseArgv.h>
#include "CAT_NiftiLib.h"
#include "CAT_Vol.h"
#include "CAT_Math.h"

int fast = 0;
int verbose = 0;
int n_avgs = 2;
int n_median_filter = 2;
double range = 0.45;
double downsample = 0.0;
double fill_thresh = 0.5;
double correct_voxelsize = 0.5;

static ArgvInfo argTable[] = {
  {"-verbose", ARGV_CONSTANT, (char *) 1, (char *) &verbose,
    "Enable verbose mode. Provides detailed output during processing for debugging\n\
    and monitoring."},

  {"-fast", ARGV_CONSTANT, (char *) 1, (char *) &fast,
    "Enable fast mode in order to get a very quick and rougher estimate of thickness only."},

  {"-n-avgs", ARGV_INT, (char *) 1, (char *) &n_avgs,
    "Specify the number of averages for distance estimation. Used for averaging\n\
    the distances in White Matter (WM) and Cerebrospinal Fluid (CSF) to obtain a\n\
    less noisy measure. A higher number results in smoother but potentially less\n\
    accurate measures."},

  {"-fill-holes", ARGV_FLOAT, (char *) 1, (char *) &fill_thresh,
    "Fill remaining holes in the PPM image using the defined threshold.\n\
    To maximize the filling-effect, this threshold should be the same as used for the\n\
    subsequent Marching Cubes approach (e.g. 0.5). Set to '0' to disable filling."},

  {"-range", ARGV_FLOAT, (char *) 1, (char *) &range,
    "Extend range for masking of euclidean distance estimation. A slight increase\n\
    of range (i.e 0.3) helps in obtaining a more stable distance estimation."},

  {"-downsample", ARGV_FLOAT, (char *) 1, (char *) &downsample,
    "Downsample PPM and GMT image to defined resolution since we do not need that 0.5mm\n\
    spatial resolution for the subsequent steps. Set to '0' to disable downsampling."},

  {"-median-filter", ARGV_INT, (char *) TRUE, (char *) &n_median_filter,
    "Specify the number of iterations to apply a median filter to areas\n\
     where the gradient of the thresholded image indicates larger clusters.\n\
     These clusters may point to potential topology artifacts and regions\n\
     with high local variations. This process helps to smooth these areas, \n\
     improving the quality of the surface reconstruction in subsequent steps."},
  
  {"-correct-voxelsize", ARGV_FLOAT, (char *) 1, (char *) &correct_voxelsize,
    "Amount of thickness correction for voxel-size, since we observed a systematic \n\
     shift to smaller thickness values of half voxel-size."},

  {NULL, ARGV_END, NULL, NULL, NULL}
};

private void usage(char *executable)
{
    char *usage_str = "\n\
Usage: %s [options] <input.nii> output_GMT.nii output_PPM.nii [output_WMD.nii output_CSD.nii GMT1.nii GMT2.nii]\n\
\n\
    This program performs projection-based cortical thickness estimation and\n\
    percentage position mapping (PPM) from a given PVE label image described in:\n\
        Dahnke R, Yotter RA, Gaser C.\n\
        Cortical thickness and central surface estimation.\n\
        Neuroimage. 2013 Jan 15;65:336-48.\n\
    \n\
    The process involves:\n\
    \n\
    1. **Distance Estimation:**\n\
       - Estimating distances in White Matter (WM) and Cerebrospinal Fluid (CSF)\n\
         by shifting the border between Gray Matter (GM)/WM and GM/CSF.\n\
       - Averaging distances over a specified number of iterations (n-avgs) to\n\
         obtain a less noisy measure.\n\
    \n\
    2. **Thickness Estimation:**\n\
       - Reconstructing sulci and optionally gyri to estimate cortical thickness.\n\
    \n\
    3. **PPM Calculation:**\n\
       - Estimating the percentage position map (PPM), representing the relative\n\
         position within the cortical ribbon.\n\
    \n\
    4. **Final Correction and Smoothing:**\n\
       - Correcting for isotropic voxel size and applying masked smoothing\n\
         to the thickness map based on the specified Full Width Half Maximum (FWHM).\n\
\n\
Options:\n\
    -verbose                   Enable verbose mode for detailed output during processing.\n\
    -n-avgs <int>              Set the number of averages for distance estimation.\n\
    -fill-holes <float>        Define the threshold to fill holes in the PPM image.\n\
    -downsample <float>        Downsample PPM and GMT image to defined resolution.\n\
    -correct-voxelsize <float> Amount of correction of thickness by voxel-size.\n\
\n\
Example:\n\
    %s -verbose -n-avgs 4 input.nii gmt_output.nii ppm_output.nii\n\n";

    fprintf(stderr,"%s\n %s\n",usage_str, executable);
}

int main(int argc, char *argv[])
{
    char out_GMT[1024], out_PPM[1024], out_CSD[1024], out_WMD[1024];
    char out_GMT1[1024], out_GMT2[1024];
    int i, j, dims[3], dims_reduced[3], replace = 0;
    float *src;
    float sum_dist, abs_dist;
    double voxelsize[3], voxelsize_reduced[3], samp[3], s[3], slope, add_value;
    double threshold[2], prctile[2];
    
    initialize_argument_processing(argc, argv);

    if (ParseArgv(&argc, argv, argTable, 0) ||(argc < 2)) {
        usage(argv[0]);
        fprintf(stderr, "     %s -help\n\n", argv[0]);
        exit(EXIT_FAILURE);
    }
    
    char *infile = argv[1];

    /* Determine output filenames based on input filename or command-line arguments */
    if (argc == 8) {
        (void) sprintf(out_GMT1, "%s", argv[6]); 
        (void) sprintf(out_GMT2, "%s", argv[7]);
    }
    if (argc >= 6) {
        (void) sprintf(out_WMD, "%s", argv[4]); 
        (void) sprintf(out_CSD, "%s", argv[5]);
    }
    if (argc >= 4) {
        (void) sprintf(out_GMT, "%s", argv[2]); 
        (void) sprintf(out_PPM, "%s", argv[3]); 
    } else {
        #if !defined(_WIN32) && !defined(_WIN64)
            (void) sprintf(out_GMT, "%s/gmt_%s", dirname(infile), basename(infile)); 
            (void) sprintf(out_PPM, "%s/ppm_%s", dirname(infile), basename(infile)); 
        #else
            fprintf(stderr,"\nUsage: %s input.nii GMT.nii PPM.nii\n\n", argv[0]);
            return(1);
        #endif
    }

    /* read source image */
    nifti_image *src_ptr = read_nifti_float(infile, &src, 0);
    if (!src_ptr) {
        fprintf(stderr,"Error reading %s.\n", infile);
        return(EXIT_FAILURE);
    }
    
    /* Number of voxels */
    int nvox = src_ptr->nvox;

    /* Prepare output NIfTI images */
    nifti_image *out_ptr = nifti_copy_nim_info(src_ptr);
    nifti_image *out_ptr_reduced = nifti_copy_nim_info(src_ptr);

    /* Retrieve dimensions and voxel size from source image */
    voxelsize[0] = src_ptr->dx;
    voxelsize[1] = src_ptr->dy;
    voxelsize[2] = src_ptr->dz;
    dims[0] = src_ptr->nx;
    dims[1] = src_ptr->ny;
    dims[2] = src_ptr->nz;
    float mean_vx_size = (voxelsize[0]+voxelsize[1]+voxelsize[2])/3.0;

    unsigned char *mask = (unsigned char *)malloc(sizeof(unsigned char)*nvox);
    float *input = (float *)malloc(sizeof(float)*nvox);
    float *dist_CSF = (float *)malloc(sizeof(float)*nvox);
    float *dist_WM = (float *)malloc(sizeof(float)*nvox);
    float *GMT = (float *)malloc(sizeof(float)*nvox);
    float *GMT1 = (float *)malloc(sizeof(float)*nvox);
    float *GMT2 = (float *)malloc(sizeof(float)*nvox);
    float *PPM = (float *)malloc(sizeof(float)*nvox);
    float *gyrus_mask = (float *)malloc(sizeof(float)*nvox);
    
    /* check for memory faults */
    if (!input || !mask || !dist_CSF || !dist_WM || !GMT || !GMT2 || !PPM) {
        fprintf(stderr,"Memory allocation error\n");
        exit(EXIT_FAILURE);
    }
    
    /* Initialize distances for CSF and WM */
    for (i = 0; i < nvox; i++) {
        dist_CSF[i] = 0.0;
        dist_WM[i]  = 0.0;
    }
    
    /* Change some defaults for fast option */
    if (fast) {
        n_avgs /= 2;
        n_median_filter = 0;
        fill_thresh = 0.0;
        downsample = 0.0;
    }
    
    /* Ensure that n_avgs is at least 1 */
    n_avgs = (n_avgs < 1) ? 1 : n_avgs;

    /* Median-filtering of input with euclidean distance helps a bit */
    localstat3(src, NULL, dims, 1, F_MEDIAN, 1, 1, DT_FLOAT32);
    
    /* Process each average for distance estimation */
    for (j = 0; j < n_avgs; j++) {
        
        /* estimate value for shifting the border to obtain a less noisy measure by averaging distances */
        add_value = ((double)j + 1.0) / ((double)n_avgs + 1.0) - 0.5;
        
        /* prepare map outside CSF and mask to obtain distance map for CSF */
        for (i = 0; i < nvox; i++) {
            input[i] = (src[i] < (CGM + add_value)) ? 1.0f : 0.0f;
            mask[i]  = (src[i] < (GWM + range)) ? 1 : 0;
        }    
    
        /* obtain CSF distance map */
        if (verbose && (j == 0)) fprintf(stderr,"Estimate CSF distance map.\n");
        euclidean_distance(input, mask, dims, NULL, replace);
        for (i = 0; i < nvox; i++)
            dist_CSF[i] += input[i];
                
        /* prepare map outside WM and mask to obtain distance map for WN */
        for (i = 0; i < nvox; i++) {
            input[i] = (src[i] > (GWM + add_value)) ? 1.0f : 0.0f;
            mask[i]  = (src[i] > (CGM - range)) ? 1 : 0;
        }    
    
        /* obtain WM distance map */
        if (verbose && (j == 0)) fprintf(stderr,"Estimate WM distance map.\n");
        euclidean_distance(input, mask, dims, NULL, replace);
        for (i = 0; i < nvox; i++)
            dist_WM[i] += input[i];
    }
        
    /* Calculate average distances if n_avgs > 1 */
    if (n_avgs > 1) {
        for (i = 0; i < nvox; i++) {
            dist_CSF[i] /= (float) n_avgs;
            dist_WM[i]  /= (float) n_avgs;
        }
    }
        
    /* Estimate cortical thickness (first using sulci measures */
    if (verbose) fprintf(stderr,"Estimate thickness map.\n");
    for (i = 0; i < nvox; i++) input[i] = roundf(src[i]);

    projection_based_thickness(input, dist_WM, dist_CSF, GMT1, dims, voxelsize);
    
    /* use both reconstruction of sulci as well as gyri and use minimum of both */
    /* we need the inverse of src: 4 - src */
    for (i = 0; i < nvox; i++) input[i] = roundf(4.0 - src[i]);

    /* Then reconstruct gyri by using the inverse of src and switching the WM and CSF distance */
    projection_based_thickness(input, dist_CSF, dist_WM, GMT2, dims, voxelsize);
    
    /* Use minimum/maximum to reduce issues with meninges */
    for (i = 0; i < nvox; i++) {
        sum_dist = dist_WM[i] + dist_CSF[i];
        GMT1[i] = MIN(sum_dist, MAX(0.0, GMT1[i] - 0.125*(GMT1[i]  < sum_dist)));

        /* Limit GMT2 to thick regions */
        GMT2[i] = (GMT2[i] > 0.0) * MAX(GMT2[i], 1.75/mean_vx_size);

        GMT2[i] = MIN(sum_dist, MAX(0.0, GMT2[i] - 0.125*(GMT2[i]  < sum_dist)));
    }
        
    /* Approximate thickness values outside GM */
    if (verbose) fprintf(stderr,"Fill values using Euclidean distance approach\n");

    euclidean_distance(GMT1, NULL, dims, NULL, 1);
    euclidean_distance(GMT2, NULL, dims, NULL, 1);
    
    /* Apply iterative median filter with 4 iterations on a 2 mm resampled image */
    median_subsample3(GMT1, dims, voxelsize, 4, 2.0, DT_FLOAT32);
    median_subsample3(GMT2, dims, voxelsize, 4, 2.0, DT_FLOAT32);

    /* Use minimum of thickness measures */
    for (i = 0; i < nvox; i++) GMT[i] = MIN(GMT1[i], GMT2[i]);
   
    for (i = 0; i < nvox; i++)
        mask[i] = (GMT[i] > 1.0) ? 1 : 0;

    median3(GMT, mask, dims, 3, DT_FLOAT32);

    /* Re-estimate CSF distance using corrected GM thickness */
    for (i = 0; i < nvox; i++) {
        if ((src[i] > CGM) && (src[i] < GWM) && (GMT[i] > 1e-15))
            dist_CSF[i] = MIN(dist_CSF[i], GMT[i] - dist_WM[i]);
    }
    clip_data(dist_CSF, nvox, 0.0, 1E15, DT_FLOAT32);

    for (i = 0; i < nvox; i++)
        mask[i] = (src[i] > CGM && src[i] < GWM) ? 0 : 1;
    
    /* Initialize and estimate percentage position map (PPM) */
    for (i = 0; i < nvox; i++)
        PPM[i] = (src[i] >= GWM) ? 1.0f : 0.0f;

    /* Estimate percentage position map (PPM)
       We use a smooth map that indicates gyri and sulci to mix the separate
       estimates of PPM for gyri and sulci */
    if (verbose) fprintf(stderr,"Estimate percentage position map.\n");
    smooth_gyri_mask(src, gyrus_mask, dims, voxelsize, CGM, 8.0);
    for (i = 0; i < nvox; i++) {
        if ((src[i] > CGM) && (src[i] < GWM) && (GMT[i] > 1e-15)) {
            float PPM_sulci = MAX(0.0, GMT[i] - dist_WM[i]) / GMT[i];
            float PPM_gyri  = MAX(0.0, GMT1[i] - dist_WM[i]) / GMT1[i];
            PPM[i] = gyrus_mask[i]*PPM_sulci + (1.0 - gyrus_mask[i])*PPM_gyri;
        }
    }
    clip_data(PPM, nvox, 0.0, 1.0, DT_FLOAT32);
    
    /* Fill remaining holes with ones that may otherwise cause topology artefacts */
    if (fill_thresh > 0.0)
        fill_holes(PPM, dims, fill_thresh, 1.0, DT_FLOAT32);

    /* 2x Median-filtering of PPM with euclidean distance */
    if (!fast) localstat3(PPM, NULL, dims, 1, F_MEDIAN, 2, 1, DT_FLOAT32);

    /* Apply voxel size correction */
    for (i = 0; i < nvox; i++) {
        GMT[i] += correct_voxelsize;
        GMT[i] *= mean_vx_size;
    }

    /* Preprocessing with Median Filter:
       - Apply an iterative median filter to areas where the gradient of
         the thresholded image indicates larger clusters.
       - Use a weighted average of the original and median filtered images.
       - Weighting is estimated using gradient of the input image and.
         morphological operations to find larger clusters */
    if (n_median_filter) {
        float *vol_smoothed = (float *)malloc(sizeof(float)*nvox);
        
        for (i = 0; i < nvox; i++) vol_smoothed[i] = PPM[i];

        s[0] = s[1] = s[2] = 4.0;
        smooth3(vol_smoothed, dims, voxelsize, s, 0, DT_FLOAT32);

        for (i = 0; i < nvox; i++) vol_smoothed[i] = PPM[i] - vol_smoothed[i];

        prctile[0] = 0.1; prctile[1] = 99.0;
        get_prctile(vol_smoothed, dims[0]*dims[1]*dims[2], threshold, prctile, 1, DT_FLOAT32);  

        /* Threshold the difference image */
        for (i = 0; i < nvox; i++)
            vol_smoothed[i] = ((vol_smoothed[i] > threshold[1]) && (GMT[i] > 1.5)) ? 1.0 : 0.0;
            
        /* Apply morphological operations to find and slightly increase 
           regions with larger clusters */
        morph_close(vol_smoothed, dims, 1, 0.5, DT_FLOAT32);
        morph_open(vol_smoothed, dims, 1, 0.0, 0, DT_FLOAT32);
        morph_dilate(vol_smoothed, dims, 3, 0.0, DT_FLOAT32);
        
        /* Smooth the gradient image to later apply weighted average */
        s[0] = s[1] = s[2] = 3.0;
        smooth3(vol_smoothed, dims, voxelsize, s, 0, DT_FLOAT32);

        for (i = 0; i < nvox; i++)
            input[i] = PPM[i];
            
        /* Apply iterative median filter with Euclidean distance */
        localstat3(input, NULL, dims, 1, F_MEDIAN, n_median_filter, 1, DT_FLOAT32);

        /* Calculate weighted average of filtered and original input */
        for (i = 0; i < nvox; i++)
            PPM[i] = (1.0 - vol_smoothed[i])*PPM[i] + vol_smoothed[i]*input[i];
            
        free(vol_smoothed);
    }
    
    clip_data(PPM, nvox, 0.0, 1.0, DT_FLOAT32);
    
    /* Downsample images */
    if (downsample > 0.0) {

        for (i = 0; i<3; i++) {
            s[i] = 1.2;
            voxelsize_reduced[i] = downsample;
            samp[i] = voxelsize_reduced[i]/voxelsize[i];
        }
    
        /*  Define grid dimensions */
        for (i = 0; i<3; i++) 
            dims_reduced[i] = (int) ceil((dims[i]-1)/((double) samp[i]))+1;
    
        out_ptr_reduced->nvox = dims_reduced[0] * dims_reduced[1] * dims_reduced[2];
        
        /* Correct affine matrix */
        for (i = 0; i < 3; i++) {
            for (j = 0; j < 3; j++) {
                out_ptr_reduced->sto_xyz.m[i][j] = out_ptr->sto_xyz.m[i][j]*samp[i];
            }
        }
        out_ptr_reduced->sto_ijk = nifti_mat44_inverse(out_ptr_reduced->sto_xyz);
    
        smooth3(GMT, dims, voxelsize, s, 0, DT_FLOAT32);
        smooth3(PPM, dims, voxelsize, s, 0, DT_FLOAT32);

        /* Save GMT and PPM image */
        slope = 1.0;
        float *vol_reduced = (float *)malloc(sizeof(float)*out_ptr_reduced->nvox);

        subsample3(GMT, vol_reduced, dims, dims_reduced, DT_FLOAT32);
        if (!write_nifti_float(out_GMT, vol_reduced, DT_FLOAT32, slope, dims_reduced, voxelsize_reduced, out_ptr_reduced)) 
            exit(EXIT_FAILURE);
    
        subsample3(PPM, vol_reduced, dims, dims_reduced, DT_FLOAT32);
        if (!write_nifti_float(out_PPM, vol_reduced, DT_FLOAT32, slope, dims_reduced, voxelsize_reduced, out_ptr_reduced)) 
            exit(EXIT_FAILURE);

        free(vol_reduced);
    } else {
        if (!write_nifti_float(out_GMT, GMT, DT_FLOAT32, slope, dims, voxelsize, out_ptr)) 
            exit(EXIT_FAILURE);
    
        if (!write_nifti_float(out_PPM, PPM, DT_FLOAT32, slope, dims, voxelsize, out_ptr)) 
            exit(EXIT_FAILURE);
    }

    if (argc >= 6) {
        if (!write_nifti_float(out_CSD, dist_CSF, DT_FLOAT32, slope, dims, voxelsize, out_ptr)) 
            exit(EXIT_FAILURE);
        if (!write_nifti_float(out_WMD, dist_WM, DT_FLOAT32, slope, dims, voxelsize, out_ptr)) 
            exit(EXIT_FAILURE);
    }

    if (argc == 8) {
        for (i = 0; i < nvox; i++)
            mask[i] = (GMT1[i] > 1.0) ? 1 : 0;
        median3(GMT1, mask, dims, 3, DT_FLOAT32);
        median3(GMT2, mask, dims, 3, DT_FLOAT32);
        for (i = 0; i < nvox; i++) {
            GMT1[i] *= mean_vx_size;
            GMT2[i] *= mean_vx_size;
        }
        if (!write_nifti_float(out_GMT1, GMT1, DT_FLOAT32, slope, dims, voxelsize, out_ptr)) 
            exit(EXIT_FAILURE);
        if (!write_nifti_float(out_GMT2, GMT2, DT_FLOAT32, slope, dims, voxelsize, out_ptr)) 
            exit(EXIT_FAILURE);
    }
            
    free(mask);
    free(dist_CSF);
    free(dist_WM);
    free(GMT);
    free(GMT1);
    free(GMT2);
    free(PPM);
    free(input);
    free(gyrus_mask);

    return(EXIT_SUCCESS);

}
