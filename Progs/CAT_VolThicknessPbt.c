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

#include <bicpl.h>
#include <ParseArgv.h>
#include "CAT_NiftiLib.h"
#include "CAT_Vol.h"

int verbose = 0;
int n_avgs = 4;
int thin_cortex = 1;

double fwhm = 3.0;
double min_thickness = 0.5;
double max_thickness = 5.0;

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

  {"-min-thickness", ARGV_FLOAT, (char *) 1, (char *) &min_thickness,
    "Set the minimum thickness that is expected. Values below that minimum thickness\n\
    are set to zero and will be approximated by using the replace option in the vbdist approach."},

  {"-max-thickness", ARGV_FLOAT, (char *) 1, (char *) &max_thickness,
    "Set the maximum thickness that is expected. Values exceeding that maximum thickness\n\
    are set to that value."},

  {"-no-thin-cortex", ARGV_CONSTANT, (char *) 0, (char *) &thin_cortex,
    "Disable the correction for the typical underestimation of GM thickness in data\n\
    where the border between GM and WM is not correctly estimated by the PVE segmentation."},

  {NULL, ARGV_END, NULL, NULL, NULL}
};

private void usage(char *executable)
{
    char *usage_str = "\n\
Usage: %s [options] <input.nii> [output_GMT.nii output_PPM.nii output_WMD.nii output_CSD.nii]\n\
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
       - Using minimum thickness measures from sulci and gyri for more accurate\n\
         representation.\n\
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
    -verbose           Enable verbose mode for detailed output during processing.\n\
    -n-avgs <int>      Set the number of averages for distance estimation.\n\
    -fwhm <float>      Define FWHM for final thickness smoothing.\n\
    -no-min-thickness  Use a simpler thickness estimation approach based on sulci only.\n\
    -no-thin-cortex    Do not slightly shift border between GM/WM for thinner cortices.\n\
    -min-thickness     Set the minimum thickness that is expected.\n\
    -max-thickness     Set the maximum thickness that is expected.\n\
\n\
Example:\n\
    %s -verbose -n-avgs 4 -fwhm 2.5 input.nii gmt_output.nii ppm_output.nii\n\n";

    fprintf(stderr,"%s\n %s\n",usage_str, executable);
}

int main(int argc, char *argv[])
{
    char *infile, out_GMT[1024], out_PPM[1024], out_CSD[1024], out_WMD[1024];
    int i, j, dims[3], replace = 0;
    float *input, *src, *dist_CSF, *dist_WM, *GMT, *GMT2, *PPM;
    float dist_CSF_val, dist_WM_val, mean_vx_size;
    float mean_GMT, abs_dist;
    unsigned char *mask;
    double voxelsize[3], slope, add_value;
    nifti_image *src_ptr, *out_ptr;
    
    if (ParseArgv(&argc, argv, argTable, 0) ||(argc < 2)) {
        usage(argv[0]);
        fprintf(stderr, "     %s -help\n\n", argv[0]);
        exit(EXIT_FAILURE);
    }
    
    initialize_argument_processing(argc, argv);

    infile  = argv[1];

    /* Determine output filenames based on input filename or command-line arguments */
    if(argc == 6) {
        (void) sprintf(out_WMD, "%s", argv[4]); 
        (void) sprintf(out_CSD, "%s", argv[5]);
    }
    if(argc >= 4) {
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
    if (!src_ptr) {
        fprintf(stderr,"Error reading %s.\n", infile);
        return(EXIT_FAILURE);
    }

    /* Prepare output NIfTI images */
    out_ptr = nifti_copy_nim_info(src_ptr);

    /* Retrieve dimensions and voxel size from source image */
    voxelsize[0] = src_ptr->dx;
    voxelsize[1] = src_ptr->dy;
    voxelsize[2] = src_ptr->dz;
    dims[0] = src_ptr->nx;
    dims[1] = src_ptr->ny;
    dims[2] = src_ptr->nz;
    mean_vx_size = (voxelsize[0]+voxelsize[1]+voxelsize[2])/3.0;

    mask     = (unsigned char *)malloc(sizeof(unsigned char)*src_ptr->nvox);
    input    = (float *)malloc(sizeof(float)*src_ptr->nvox);
    dist_CSF = (float *)malloc(sizeof(float)*src_ptr->nvox);
    dist_WM  = (float *)malloc(sizeof(float)*src_ptr->nvox);

    GMT = (float *)malloc(sizeof(float)*src_ptr->nvox);
    PPM = (float *)malloc(sizeof(float)*src_ptr->nvox);
    
    /* check for memory faults */
    if (!input || !mask || !dist_CSF || !dist_WM || !GMT || !PPM) {
        fprintf(stderr,"Memory allocation error\n");
        exit(EXIT_FAILURE);
    }
    
    /* Initialize distances for CSF and WM */
    for (i = 0; i < src_ptr->nvox; i++) {
        dist_CSF[i] = 0.0;
        dist_WM[i]  = 0.0;
    }
    
    /* Process each average for distance estimation */
    for (j = 0; j < n_avgs; j++) {
        
        /* estimate value for shifting the border to obtain a less noisy measure by averaging distances */
        add_value = ((double)j + 1.0) / ((double)n_avgs + 1.0) - 0.5;        
        
        /* prepare map outside CSF and mask to obtain distance map for CSF */
        for (i = 0; i < src_ptr->nvox; i++) {
            input[i] = (src[i] <= (CGM+add_value)) ? 1.0f : 0.0f;
            mask[i]  = (src[i] <= GWM + 0.1*thin_cortex) ? 1 : 0;
        }    
    
        /* obtain CSF distance map */
        if (verbose && (j == 0)) fprintf(stderr,"Estimate CSF distance map.\n");
        vbdist(input, mask, dims, NULL, replace);
        for (i = 0; i < src_ptr->nvox; i++)
            dist_CSF[i] += input[i];
                
        /* prepare map outside WM and mask to obtain distance map for WN */
        for (i = 0; i < src_ptr->nvox; i++) {
            input[i] = (src[i] >= (GWM+add_value)) ? 1.0f : 0.0f;
            mask[i]  = (src[i] >= CGM) ? 1 : 0;
        }    
    
        /* obtain WM distance map */
        if (verbose && (j == 0)) fprintf(stderr,"Estimate WM distance map.\n");
        vbdist(input, mask, dims, NULL, replace);
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

    /* fill small holes that cause topology artefacts */
    fill_holes(dist_WM, 1E-3, dims, DT_FLOAT32);
    fill_holes(dist_CSF,1E-3, dims, DT_FLOAT32);
    
    /* Estimate cortical thickness (first using sulci measures */
    if (verbose) fprintf(stderr,"Estimate thickness map.\n");
    for (i = 0; i < src_ptr->nvox; i++) input[i] = src[i];

    projection_based_thickness(input, dist_WM, dist_CSF, GMT, dims, voxelsize);

    /* use minimum/maximum to reduce issues with meninges */
    for (i = 0; i < src_ptr->nvox; i++)
        GMT[i] = MIN(dist_WM[i]+dist_CSF[i], GMT[i]);
    
    /* get mean thickness */
    mean_GMT = get_mean_float(GMT, src_ptr->nvox);
    
    /* use both reconstruction of sulci as well as gyri and use minimum of both */
    /* we need the inverse of src: 4 - src */
    for (i = 0; i < src_ptr->nvox; i++)
        input[i] = (4.0 - src[i]);

    GMT2 = (float *)malloc(sizeof(float)*src_ptr->nvox);
    if (!GMT2) {
        fprintf(stderr,"Memory allocation error\n");
        exit(EXIT_FAILURE);
    }

    /* then reconstruct gyri by using the inverse of src and switching the WM and CSF distance */
    projection_based_thickness(input, dist_CSF, dist_WM, GMT2, dims, voxelsize);

    for (i = 0; i < src_ptr->nvox; i++)
        GMT2[i] = MIN(dist_WM[i]+dist_CSF[i], GMT2[i]);

    /* get overall mean thickness */
    mean_GMT = mean_GMT/2.0 + get_mean_float(GMT2, src_ptr->nvox)/2.0;

    /* use average of thickness measures */
    for (i = 0; i < src_ptr->nvox; i++)      
        GMT[i] = (GMT[i] + GMT2[i])/2;

    free(GMT2);
    free(input);
   
    median3(GMT, NULL, dims, DT_FLOAT32);

    /* Re-estimate CSF distance using corrected GM thickness */
    for (i = 0; i < src_ptr->nvox; i++)
        dist_CSF[i] = GMT[i] - dist_WM[i];
    
    /* Approximate thickness values outside GM or below minimum thickness */
    for (i = 0; i < src_ptr->nvox; i++) {
        GMT[i]  = ((GMT[i]*mean_vx_size) < min_thickness) ? 0.0 : GMT[i];
        GMT[i]  = ((GMT[i]*mean_vx_size) > max_thickness) ? max_thickness : GMT[i];
        mask[i] = (GMT[i] == 0) ? 1 : 0;
    }
    vbdist(GMT, mask, dims, NULL, 1);

    /* Apply final smoothing */
    if (fwhm > 0.0) {
        if (verbose) fprintf(stderr,"Final correction\n");
        double s[3] = {fwhm, fwhm, fwhm};
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
            PPM[i] = MIN(dist_CSF[i], (GMT[i]-dist_WM[i])) / GMT[i];
        if (PPM[i] < 0.0) PPM[i] = 0.0;
        if (PPM[i] > 1.0) PPM[i] = 1.0;
    }
    
    /* fill small holes that cause topology artefacts */
    fill_holes(PPM, 1E-3, dims, DT_FLOAT32);

    /* Apply isotropic voxel size correction */
    for (i = 0; i < src_ptr->nvox; i++) 
        GMT[i] *= mean_vx_size;
    
    /* save GMT and PPM image */
    slope = 1.0;
    if (!write_nifti_float(out_GMT, GMT, DT_FLOAT32, slope, dims, voxelsize, out_ptr)) 
        exit(EXIT_FAILURE);
    if (!write_nifti_float(out_PPM, PPM, DT_FLOAT32, slope, dims, voxelsize, out_ptr)) 
        exit(EXIT_FAILURE);

    if(argc == 6) {
        if (!write_nifti_float(out_CSD, dist_CSF, DT_FLOAT32, slope, dims, voxelsize, out_ptr)) 
            exit(EXIT_FAILURE);
        if (!write_nifti_float(out_WMD, dist_WM, DT_FLOAT32, slope, dims, voxelsize, out_ptr)) 
            exit(EXIT_FAILURE);
    }
            
    free(mask);
    free(dist_CSF);
    free(dist_WM);
    free(GMT);
    free(PPM);

    return(EXIT_SUCCESS);

}