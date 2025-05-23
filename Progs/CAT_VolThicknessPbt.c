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

int verbose = 0;
int n_avgs = 8;
int n_median_filter = 2;
double sharpening = 0.0;
double downsample = 0.0;
double fwhm = 0.0;
double fill_thresh = 0.5;

static ArgvInfo argTable[] = {
  {"-verbose", ARGV_CONSTANT, (char *) 1, (char *) &verbose,
    "Enable verbose mode. Provides detailed output during processing for debugging\n\
    and monitoring."},

  {"-n-avgs", ARGV_INT, (char *) 1, (char *) &n_avgs,
    "Specify the number of averages for distance estimation. Used for averaging\n\
    the distances in White Matter (WM) and Cerebrospinal Fluid (CSF) to obtain a\n\
    less noisy measure. A higher number results in smoother but potentially less\n\
    accurate measures."},

  {"-fwhm", ARGV_FLOAT, (char *) 1, (char *) &fwhm,
    "Set the Full Width Half Maximum (FWHM) value for final thickness smoothing.\n\
    This value determines the extent of smoothing applied, using a mask to prevent\n\
    smearing values outside the Gray Matter (GM) areas."},

  {"-fill-holes", ARGV_FLOAT, (char *) 1, (char *) &fill_thresh,
    "Fill remaining holes in the PPM image using the defined threshold.\n\
    To maximize the filling-effect, this threshold should be the same as used for the\n\
    subsequent Marching Cubes approach (e.g. 0.5). Set to '0' to disable filling."},

  {"-downsample", ARGV_FLOAT, (char *) 1, (char *) &downsample,
    "Downsample PPM and GMT image to defined resolution since we do not need that 0.5mm\n\
    spatial resolution for the subsequent steps. Set to '0' to disable downsampling."},

  {"-median-filter", ARGV_INT, (char *) TRUE, (char *) &n_median_filter,
    "Specify the number of iterations to apply a median filter to areas\n\
     where the gradient of the thresholded image indicates larger clusters.\n\
     These clusters may point to potential topology artifacts and regions\n\
     with high local variations. This process helps to smooth these areas, \n\
     improving the quality of the surface reconstruction in subsequent steps."},
  
  {"-sharpen", ARGV_FLOAT, (char *) 1, (char *) &sharpening,
    "Amount of sharpening the PPM map by adding the difference between the unsmoothed and \n\
     smoothed PPM map. Set to '0' to disable sharpening"},

  {NULL, ARGV_END, NULL, NULL, NULL}
};

private void usage(char *executable)
{
    char *usage_str = "\n\
Usage: %s [options] <input.nii> output_GMT.nii output_PPM.nii [output_WMD.nii output_CSD.nii]\n\
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
    -verbose                Enable verbose mode for detailed output during processing.\n\
    -n-avgs <int>           Set the number of averages for distance estimation.\n\
    -fwhm <float>           Define FWHM for final thickness smoothing.\n\
    -fill-holes <float>     Define the threshold to fill holes in the PPM image.\n\
    -downsample <float>     Downsample PPM and GMT image to defined resolution.\n\
    -sharpen <float>        Amount of sharpening the PPM map.\n\
\n\
Example:\n\
    %s -verbose -n-avgs 4 -fwhm 2 input.nii gmt_output.nii ppm_output.nii\n\n";

    fprintf(stderr,"%s\n %s\n",usage_str, executable);
}

int main(int argc, char *argv[])
{
    char *infile, out_GMT[1024], out_PPM[1024], out_CSD[1024], out_WMD[1024];
    int i, j, dims[3], dims_reduced[3], replace = 0;
    float *input, *src, *dist_CSF, *dist_WM, *GMT, *GMT1, *GMT2, *PPM, *PPM0, *vol_reduced;
    float *vol_smoothed, dist_CSF_val, dist_WM_val, mean_vx_size;
    float sum_dist, abs_dist;
    unsigned char *mask;
    double voxelsize[3], voxelsize_reduced[3], samp[3], s[3], slope, add_value;
    double threshold[2], prctile[2];
    nifti_image *src_ptr, *out_ptr, *out_ptr_reduced;
    
    initialize_argument_processing(argc, argv);

    if (ParseArgv(&argc, argv, argTable, 0) ||(argc < 2)) {
        usage(argv[0]);
        fprintf(stderr, "     %s -help\n\n", argv[0]);
        exit(EXIT_FAILURE);
    }
    
    infile  = argv[1];

    /* Determine output filenames based on input filename or command-line arguments */
    if (argc == 6) {
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
    src_ptr = read_nifti_float(infile, &src, 0);
    if (!src_ptr) {
        fprintf(stderr,"Error reading %s.\n", infile);
        return(EXIT_FAILURE);
    }

    /* Prepare output NIfTI images */
    out_ptr = nifti_copy_nim_info(src_ptr);
    out_ptr_reduced = nifti_copy_nim_info(src_ptr);

    /* Retrieve dimensions and voxel size from source image */
    voxelsize[0] = src_ptr->dx;
    voxelsize[1] = src_ptr->dy;
    voxelsize[2] = src_ptr->dz;
    dims[0] = src_ptr->nx;
    dims[1] = src_ptr->ny;
    dims[2] = src_ptr->nz;
    mean_vx_size = (voxelsize[0]+voxelsize[1]+voxelsize[2])/3.0;

    mask = (unsigned char *)malloc(sizeof(unsigned char)*src_ptr->nvox);
    input = (float *)malloc(sizeof(float)*src_ptr->nvox);
    dist_CSF = (float *)malloc(sizeof(float)*src_ptr->nvox);
    dist_WM = (float *)malloc(sizeof(float)*src_ptr->nvox);
    GMT = (float *)malloc(sizeof(float)*src_ptr->nvox);
    GMT1 = (float *)malloc(sizeof(float)*src_ptr->nvox);
    GMT2 = (float *)malloc(sizeof(float)*src_ptr->nvox);
    PPM = (float *)malloc(sizeof(float)*src_ptr->nvox);
    
    /* check for memory faults */
    if (!input || !mask || !dist_CSF || !dist_WM || !GMT || !GMT2 || !PPM) {
        fprintf(stderr,"Memory allocation error\n");
        exit(EXIT_FAILURE);
    }
    
    /* Initialize distances for CSF and WM */
    for (i = 0; i < src_ptr->nvox; i++) {
        dist_CSF[i] = 0.0;
        dist_WM[i]  = 0.0;
    }
    
    /* Ensure that n_avgs is at least 1 */
    n_avgs = (n_avgs < 1) ? 1 : n_avgs;

    /* Median-filtering of input with euclidean distance helps a bit */
    localstat3(src, NULL, dims, 1, F_MEDIAN, 1, 1, DT_FLOAT32);

    double range = 0.3; /* Value of 0.3 was worst */
    
    /* Process each average for distance estimation */
    for (j = 0; j < n_avgs; j++) {
        
        /* estimate value for shifting the border to obtain a less noisy measure by averaging distances */
        add_value = ((double)j + 1.0) / ((double)n_avgs + 1.0) - 0.5;
        
        /* prepare map outside CSF and mask to obtain distance map for CSF */
        for (i = 0; i < src_ptr->nvox; i++) {
            input[i] = (src[i] < (CGM + add_value)) ? 1.0f : 0.0f;
            mask[i]  = (src[i] < GWM + range) ? 1 : 0;
        }    
    
        /* obtain CSF distance map */
        if (verbose && (j == 0)) fprintf(stderr,"Estimate CSF distance map.\n");
        euclidean_distance(input, mask, dims, NULL, replace);
        for (i = 0; i < src_ptr->nvox; i++)
            dist_CSF[i] += input[i];
                
        /* prepare map outside WM and mask to obtain distance map for WN */
        for (i = 0; i < src_ptr->nvox; i++) {
            input[i] = (src[i] > (GWM + add_value)) ? 1.0f : 0.0f;
            mask[i]  = (src[i] > CGM - range) ? 1 : 0;
        }    
    
        /* obtain WM distance map */
        if (verbose && (j == 0)) fprintf(stderr,"Estimate WM distance map.\n");
        euclidean_distance(input, mask, dims, NULL, replace);
        for (i = 0; i < src_ptr->nvox; i++)
            dist_WM[i] += input[i];
    }
        
    /* Calculate average distances if n_avgs > 1 */
    if (n_avgs > 1) {
        for (i = 0; i < src_ptr->nvox; i++) {
            dist_CSF[i] /= (float) n_avgs;
            dist_WM[i]  /= (float) n_avgs;
        }
    }
    
    int gyrus_recon2 = 0; // not working yet!!!
    
    /* Estimate cortical thickness (first using sulci measures */
    if (verbose) fprintf(stderr,"Estimate thickness map.\n");
    for (i = 0; i < src_ptr->nvox; i++) input[i] = roundf(src[i]);

    projection_based_thickness(input, dist_WM, dist_CSF, GMT1, dims, voxelsize);
    
    if (gyrus_recon2) {
        vol_approx(GMT1, dims, voxelsize);
        for (i = 0; i < src_ptr->nvox; i++) {
            dist_WM[i] = MIN(dist_WM[i], GMT1[i]);
            dist_CSF[i] = MIN(dist_CSF[i], GMT1[i]);
            mask[i] = (dist_WM[i] > 0.0) & (dist_CSF[i] > 0.0);
            if (mask[i] > 0) 
                dist_CSF[i] = MIN(dist_CSF[i], GMT1[i] - dist_WM[i] );
        }
    }

    /* use both reconstruction of sulci as well as gyri and use minimum of both */
    /* we need the inverse of src: 4 - src */
    for (i = 0; i < src_ptr->nvox; i++) input[i] = roundf(4.0 - src[i]);

    /* Then reconstruct gyri by using the inverse of src and switching the WM and CSF distance */
    projection_based_thickness(input, dist_CSF, dist_WM, GMT2, dims, voxelsize);
    
    if (gyrus_recon2) {
        vol_approx(GMT2, dims, voxelsize);
        memcpy(GMT, GMT2, src_ptr->nvox*sizeof(float));
        median3(GMT, NULL, dims, 3, DT_FLOAT32);
        for (i = 0; i < src_ptr->nvox; i++) {
            if ((mask[i] > 0.0) & (dist_CSF[i] > GMT[i])) /* GMT is here the median filtered GMT2 */ 
                dist_WM[i] = MIN(dist_WM[i], GMT2[i] - dist_CSF[i] );
        }
    }

    /* Correct distance values by half of a voxel */
    double dist_corr = 0.0; // 0.5 was not working better, which I used before
    for (i = 0; i < src_ptr->nvox; i++) {
        dist_CSF[i] -= dist_corr;
        dist_WM[i]  -= dist_corr;
        if (dist_CSF[i] < 0.0) dist_CSF[i] = 0.0;
        if (dist_WM[i]  < 0.0) dist_WM[i]  = 0.0;
    }

    if (gyrus_recon2 == 0) {
        /* Use minimum/maximum to reduce issues with meninges */
        for (i = 0; i < src_ptr->nvox; i++) {
            sum_dist = dist_WM[i] + dist_CSF[i];
            GMT1[i] = MIN(sum_dist, MAX(0.0, GMT1[i] - 0.5*(GMT1[i]  < sum_dist)));
        }
    
        /* Limit GMT2 to thick regions */
        for (i = 0; i < src_ptr->nvox; i++)
            GMT2[i] = (GMT2[i] > 0.0) * MAX(GMT2[i], 1.75/mean_vx_size);
    
        for (i = 0; i < src_ptr->nvox; i++) {
            sum_dist = dist_WM[i] + dist_CSF[i];
            GMT2[i] = MIN(sum_dist, MAX(0.0, GMT2[i] - 0.5*(GMT2[i]  < sum_dist)));
        }
            
        /* Use minimum of thickness measures */
        for (i = 0; i < src_ptr->nvox; i++)
            GMT[i] = MIN(GMT1[i], GMT2[i]);
       
        for (i = 0; i < src_ptr->nvox; i++)
            mask[i] = (GMT[i] > 1.0) ? 1 : 0;
    
        median3(GMT, mask, dims, 3, DT_FLOAT32);
        //median3(dist_WM, NULL, dims, 1, DT_FLOAT32);
        //median3(dist_CSF, NULL, dims, 1, DT_FLOAT32);
    
        /* Re-estimate CSF distance using corrected GM thickness */
        for (i = 0; i < src_ptr->nvox; i++) {
            if ((src[i] > CGM) && (src[i] < GWM) && (GMT[i] > 1e-15))
                dist_CSF[i] = MIN(dist_CSF[i], GMT[i] - dist_WM[i]);
            if (dist_CSF[i] < 0.0) dist_CSF[i] = 0.0;
        }
    }

    for (i = 0; i < src_ptr->nvox; i++)
        mask[i] = (src[i] > CGM && src[i] < GWM) ? 0 : 1;

    /* Approximate thickness values outside GM */
    if (fwhm >= 0.0) {
        if (verbose) fprintf(stderr,"Fill values using Euclidean distance approach\n");
        euclidean_distance(GMT, mask, dims, NULL, 1);
        laplace3R(GMT, mask, dims, 0.1);
    } else {
        if (verbose) fprintf(stderr,"Fill values using Approximation approach\n");
        vol_approx(GMT, dims, voxelsize);
    }
    
    /* Apply final smoothing */
    if (fwhm > 0.0) {
        if (verbose) fprintf(stderr,"Final correction\n");
        s[0] = s[1] = s[2] = fwhm;
        smooth3(GMT, dims, voxelsize, s, 1, DT_FLOAT32);
    }
    
    /* Initialize and estimate percentage position map (PPM) */
    for (i = 0; i < src_ptr->nvox; i++)
        PPM[i] = (src[i] >= GWM) ? 1.0f : 0.0f;

    /* Estimate percentage position map (PPM)
       We first create a corrected CSF distance map with reconstructed sulci.
       If gyri were reconstructed too than also the dist_WM have to be
       corrected to avoid underestimation of the position map with surfaces 
       running to close to the WM. */
    if (verbose) fprintf(stderr,"Estimate percentage position map.\n");
    for (i = 0; i < src_ptr->nvox; i++) {
        if ((src[i] > CGM) && (src[i] < GWM) && (GMT[i] > 1e-15))
            PPM[i] = MAX(0.0, dist_CSF[i] / GMT[i]);
        if (PPM[i] < 0.0) PPM[i] = 0.0;
        if (PPM[i] > 1.0) PPM[i] = 1.0;
        if ((PPM[i] < 0.9) && (src[i] > GWM)) PPM[i] = 1.0;
    }
    
    /* Fill remaining holes with ones that may otherwise cause topology artefacts */
    if (fill_thresh > 0.0)
        fill_holes(PPM, dims, fill_thresh, 1.0, DT_FLOAT32);

    /* 2x Median-filtering of PPM with euclidean distance */
    localstat3(PPM, NULL, dims, 1, F_MEDIAN, 2, 1, DT_FLOAT32);

    /* Apply isotropic voxel size correction */
    for (i = 0; i < src_ptr->nvox; i++) 
        GMT[i] *= mean_vx_size;

    /* Preprocessing with Median Filter:
       - Apply an iterative median filter to areas where the gradient of
         the thresholded image indicates larger clusters.
       - Use a weighted average of the original and median filtered images.
       - Weighting is estimated using gradient of the input image and.
         morphological operations to find larger clusters */
    if (n_median_filter) {
        vol_smoothed = (float *)malloc(sizeof(float)*src_ptr->nvox);
        
        for (i = 0; i < src_ptr->nvox; i++) vol_smoothed[i] = PPM[i];

        s[0] = s[1] = s[2] = 4.0;
        smooth3(vol_smoothed, dims, voxelsize, s, 0, DT_FLOAT32);

        for (i = 0; i < src_ptr->nvox; i++) vol_smoothed[i] = PPM[i] - vol_smoothed[i];

        prctile[0] = 0.1; prctile[1] = 99.0;
        get_prctile(vol_smoothed, dims[0]*dims[1]*dims[2], threshold, prctile, 1, DT_FLOAT32);  

        /* Threshold the difference image */
        for (i = 0; i < src_ptr->nvox; i++)
            vol_smoothed[i] = ((vol_smoothed[i] > threshold[1]) && (GMT[i] > 1.5)) ? 1.0 : 0.0;
            
        /* Apply morphological operations to find and slightly increase 
           regions with larger clusters */
        morph_close(vol_smoothed, dims, 1, 0.5, DT_FLOAT32);
        morph_open(vol_smoothed, dims, 1, 0.0, DT_FLOAT32);
        morph_dilate(vol_smoothed, dims, 3, 0.0, DT_FLOAT32);
        
        /* Smooth the gradient image to later apply weighted average */
        s[0] = s[1] = s[2] = 3.0;
        smooth3(vol_smoothed, dims, voxelsize, s, 0, DT_FLOAT32);

        for (i = 0; i < src_ptr->nvox; i++)
            input[i] = PPM[i];
            
        /* Apply iterative median filter with Euclidean distance */
        localstat3(input, NULL, dims, 1, F_MEDIAN, n_median_filter, 1, DT_FLOAT32);

        /* Calculate weighted average of filtered and original input */
        for (i = 0; i < src_ptr->nvox; i++)
            PPM[i] = (1.0 - vol_smoothed[i])*PPM[i] + vol_smoothed[i]*input[i];
            
        free(vol_smoothed);
    }
    
    /* Apply sharpening by subtracting smoothed PPM map */
    if (sharpening > 0.0) {
        vol_smoothed = (float *)malloc(sizeof(float)*src_ptr->nvox);
        PPM0 = (float *)malloc(sizeof(float)*src_ptr->nvox);

        for (i = 0; i < src_ptr->nvox; i++) PPM0[i] = PPM[i];
        for (i = 0; i < src_ptr->nvox; i++) vol_smoothed[i] = PPM[i];

        j = 1;
        
        /* Do 3 steps with j = [1 2 4] */
        while (j < 5) {
            /* Smoothing PPM map */
            s[0] = s[1] = s[2] = 4.0*j;
            smooth3(vol_smoothed, dims, voxelsize, s, 0, DT_FLOAT32);
            
            /* Use iterative levels of smoothing */
            for (i = 0; i < src_ptr->nvox; i++)
                PPM[i] += sharpening*j*(PPM[i] - vol_smoothed[i] + 0.001*j);
            j *= 2;
        }
        
        /* Correct thickness values */
        for (i = 0; i < src_ptr->nvox; i++) vol_smoothed[i] = PPM[i] - PPM0[i];
        median3(vol_smoothed, NULL, dims, 1, DT_FLOAT32);
        for (i = 0; i < src_ptr->nvox; i++) GMT[i] += 2*vol_smoothed[i];
        
        free(vol_smoothed);
        free(PPM0);
    }

    for (i = 0; i < src_ptr->nvox; i++) {
        if (PPM[i] < 0.0) PPM[i] = 0.0;
        if (PPM[i] > 1.0) PPM[i] = 1.0;
    }
    
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
        //smooth3(PPM, dims, voxelsize, s, 0, DT_FLOAT32);

        /* Save GMT and PPM image */
        slope = 1.0;
        vol_reduced = (float *)malloc(sizeof(float)*out_ptr_reduced->nvox);

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

    if (argc == 6) {
        if (!write_nifti_float(out_CSD, dist_CSF, DT_FLOAT32, slope, dims, voxelsize, out_ptr)) 
            exit(EXIT_FAILURE);
        if (!write_nifti_float(out_WMD, dist_WM, DT_FLOAT32, slope, dims, voxelsize, out_ptr)) 
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

    return(EXIT_SUCCESS);

}