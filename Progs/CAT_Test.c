/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 *
 */

#include <float.h>
#include <math.h>

#include <bicpl.h>
#include <ParseArgv.h>
#include "CAT_Amap.h"
#include "CAT_NiftiLib.h"
#include "CAT_Vol.h"

char *label_filename;
int iters_amap = 50;
int subsample = 96;
int iters_ICM = 50;
int pve = 1;
int write_seg[3] = {0, 1, 0};
int write_label = 1;
int write_corr = 1;
int write_bias = 0;
int debug = 0;
double weight_LAS = 0.5;
double weight_MRF = 0.0;
double bias_fwhm = 10.0;

static ArgvInfo argTable[] = {
    {"-label", ARGV_STRING, (char *) 1, (char *) &label_filename, 
         "File containing segmentation labels for initialization."},
         
    {"-iters", ARGV_INT, (char *) 1, (char *) &iters_amap,
         "Specifies the number of iterations for the Amap approach to terminate."},
         
    {"-sub", ARGV_INT, (char *) 1, (char *) &subsample,
         "Defines the subsampling factor for Amap approach, which will be scaled\n\
         internally by voxel size."},
         
    {"-iters-icm", ARGV_INT, (char *) 1, (char *) &iters_ICM,
         "Sets the number of iterations for the Iterative Conditional Mode (ICM)\n\
         algorithm."},
         
    {"-mrf", ARGV_FLOAT, (char *) 1, (char *) &weight_MRF,
         "Determines the weight of the Markov Random Field (MRF) prior, a value\n\
         between 0 and 1."},
         
    {"-las", ARGV_FLOAT, (char *) 1, (char *) &weight_LAS,
         "Determines the weight of the local adaptive segmentation (LAS), a value\n\
         between 0 and 1. Only used if bias correction is applied."},
         
    {"-bias-fwhm", ARGV_FLOAT, (char *) 1, (char *) &bias_fwhm,
         "Specifies the Full Width Half Maximum (FWHM) value for the bias correction\n\
         smoothing kernel."},
         
    {"-pve", ARGV_INT, (char *) 1, (char *) &pve,
         "Option to use Partial Volume Estimation with 5 classes (1) or not (0).\n\
         Default setting is 1."},
         
    {"-write-seg", ARGV_INT, (char *) 3, (char *) &write_seg,
         "Option to write segmentation results as separate images. Requires three integers\n\
         indicating whether to save each tissue class (CSF/GM/WM) with '1' for yes."},
         
    {"-write-label", ARGV_CONSTANT, (char *) 1, (char *) &write_label,
         "Enable writing the label image. This is the default setting."},
         
    {"-nowrite-label", ARGV_CONSTANT, (char *) 0, (char *) &write_label,
         "Disable writing the label image."},
         
    {"-write-corr", ARGV_CONSTANT, (char *) 1, (char *) &write_corr,
         "Enable writing the nu-corrected image. This is the default setting."},
         
    {"-write-bias", ARGV_CONSTANT, (char *) 1, (char *) &write_bias,
         "Enable writing the nu-correction (bias field)."},
         
    {"-nowrite-corr", ARGV_CONSTANT, (char *) 0, (char *) &write_corr,
         "Disable writing the nu-corrected image."},
         
    {"-debug", ARGV_CONSTANT, (char *) 1, (char *) &debug,
         "Enable debug mode to print additional debug information."},
         
    {NULL, ARGV_END, NULL, NULL, NULL}
};


private void
usage(
    char *executable)
{
    char *usage_str = "\n\
        CAT_VolAmap: Segmentation with adaptive MAP\n\
         usage: CAT_VolAmap [options] -label label.nii in.nii [out.nii]\n\n";
    fprintf(stderr,"%s\n %s\n",usage_str, executable);
}

int
main(int argc, char *argv[])
{
    nifti_image   *src_ptr, *label_ptr;
    int       n_classes;
    char      *input_filename, *output_filename, *basename, *extension;
    int       i, j, dims[3], n_pure_classes = 3;;
    int       x, y, z, z_area, y_dims;
    char      *arg_string, buffer[1024];
    unsigned char *label, *prob;
    float     *src, *buffer_vol, *biasfield;
    double    slope, offset, val, max_vol, min_vol, voxelsize[3];    
    char *label_arr[] = {"CSF", "GM", "WM"};

    /* Get arguments */
    if (ParseArgv(&argc, argv, argTable, 0) || (argc < 2)) {
        usage(argv[0]);
        fprintf(stderr, "     %s -help\n\n", argv[0]);
        exit(EXIT_FAILURE);
    }
    
    initialize_argument_processing(argc, argv);

    if (!get_string_argument(NULL, &input_filename)) {
        usage(argv[0]);
        fprintf(stderr, "     %s -help\n\n", argv[0]);
        return(1);
    }

    /* if not defined use original name as basename for output */
    if (!get_string_argument(NULL, &output_filename))
        output_filename = argv[1];
        
    /* get basename */
    basename = nifti_makebasename(output_filename);

    /* deal with extension */
    extension = nifti_find_file_extension(output_filename);
    
    /* if no valid extension was found use .nii */
    if (!extension) {
        fprintf(stdout,"Use .nii as extension for %s.\n",output_filename);
        strcpy(extension, ".nii");
    }

    if (iters_amap == 0) {
        fprintf(stderr,"To estimate segmentation you need at least one iteration.\n");
        exit(EXIT_FAILURE);
    }
    
    if (pve) n_classes = 5;
    else     n_classes = 3;
    
    /* read data */
    src_ptr = read_nifti_float(input_filename, &src, 0);
    
    if (!src_ptr) {
        fprintf(stderr,"Error reading %s.\n",input_filename);
        exit(EXIT_FAILURE);
    }

    /* read label and check for same size */
    if (!label_filename) {
        fprintf(stderr,"Label image has to be defined\n");
        exit(EXIT_FAILURE);
    }
            
    label = (unsigned char *)malloc(sizeof(unsigned char)*src_ptr->nvox);
    prob  = (unsigned char *)malloc(sizeof(unsigned char)*src_ptr->nvox*n_classes);
    biasfield = (float *)malloc(sizeof(float)*src_ptr->nvox);
    
    if (!label || !prob || !biasfield) {
        fprintf(stderr,"Memory allocation error\n");
        exit(EXIT_FAILURE);
    }

    /* read volume */
    label_ptr = read_nifti_float(label_filename, &buffer_vol, 0);
    if (!label_ptr) {
        fprintf(stderr,"Error reading %s.\n", label_filename);
        free(label);
        free(prob);
        return(EXIT_FAILURE);
    }
    
    /* check size */
    if (!equal_image_dimensions(src_ptr, label_ptr)) {     
        fprintf(stderr,"Label and source image have different size\n");
        free(label);
        free(prob);
        exit(EXIT_FAILURE);
    }
    
    /* scale label image to a range 0..3 */
    for (i = 0; i < src_ptr->nvox; i++) 
        label[i] = (unsigned char) round(buffer_vol[i]);

    double mean[n_classes], mu[n_pure_classes], var[n_pure_classes];
    for (i = 0; i < n_pure_classes; i++)
        mu[i] = 0;

    /* get min/max */
    min_vol = FLT_MAX; max_vol = -FLT_MAX;
    for (i = 0; i < src_ptr->nvox; i++) {
        min_vol = MIN((double)src[i], min_vol);
        max_vol = MAX((double)src[i], max_vol);
    }

    if (debug)
        fprintf(stdout,"Intensity range: %3.2f - %3.2f\n", min_vol, max_vol);

    /* correct images with values < 0 */
    if (min_vol < 0) {
        for (i = 0; i < src_ptr->nvox; i++)
            src[i] = src[i] - (float)min_vol;
    }

    voxelsize[0] = src_ptr->dx;
    voxelsize[1] = src_ptr->dy;
    voxelsize[2] = src_ptr->dz;
    dims[0] = src_ptr->nx;
    dims[1] = src_ptr->ny;
    dims[2] = src_ptr->nz;

    /* apply bias correction first for GM+WM and subsequently for WM only with less
     * smoothing to emphasize subcortical structures */
    if (bias_fwhm > 0.0) {
        fprintf(stdout,"Bias correction\n");
        correct_bias(src, biasfield, label, dims, voxelsize, bias_fwhm, weight_LAS, 0);
    }

    Amap(src, label, prob, mean, n_pure_classes, iters_amap, subsample, dims, pve, weight_MRF, voxelsize, iters_ICM, offset, bias_fwhm);

    /* PVE */
    if (pve) {
        fprintf(stdout,"Calculate Partial Volume Estimate.\n");
        Pve5(src, prob, label, mean, dims);
    }
        
    /* write nu-corrected volume */
    if (write_corr) {

        slope = 0.0;
        sprintf(buffer, "%s_corr%s",basename,extension);
        if (!write_nifti_float(buffer, src, DT_UINT16, slope, 
                        dims, voxelsize, src_ptr))
            exit(EXIT_FAILURE);
    }

    /* write bias field */
    if (write_bias) {

        slope = 0.0;
        sprintf(buffer, "%s_bias%s",basename,extension);
        if (!write_nifti_float(buffer, biasfield, DT_FLOAT, slope, 
                        dims, voxelsize, src_ptr))
            exit(EXIT_FAILURE);
    }

    /* write labeled volume */
    if (write_label) {

        /* different ranges for pve */
        if (pve) slope = 3.0/255.0;
        else slope = 1.0;

        for (i = 0; i < src_ptr->nvox; i++)
            src[i] = (float)label[i];

        sprintf(buffer, "%s_seg%s",basename,extension);
        if (!write_nifti_float(buffer, src, DT_UINT8, slope, 
                        dims, voxelsize, src_ptr))
            exit(EXIT_FAILURE);
    }
    
    /* write PVE segmentations for each class */
    if (write_seg[0] || write_seg[1] || write_seg[2]) {
        
        slope = 1.0/255.0;

        for (j = 0; j < n_pure_classes; j++) {
            if (write_seg[j]) {

                sprintf(buffer, "%s_label-%s_probseg%s",basename,label_arr[j],extension); 
                
                for (i = 0; i < src_ptr->nvox; i++)
                    src[i] = prob[i+(j*src_ptr->nvox)];
                
                if (!write_nifti_float(buffer, src, DT_UINT8, slope, 
                                dims, voxelsize, src_ptr))
                    exit(EXIT_FAILURE);
            }
        }        
    }
    
    free(prob);
    free(label);
    free(biasfield);
    
    return(EXIT_SUCCESS);
}