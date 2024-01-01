/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 *
 */

#include <float.h>
#include <math.h>

#include "ParseArgv.h"
#include "CAT_NiftiLib.h"
#include "CAT_Vol.h"

typedef  enum  { PDw, T1w, T2w, Unknown } Modality ;
                 
char *label_filename = NULL;

static ArgvInfo argTable[] = {
    {"-label", ARGV_STRING, (char *) 1, (char *) &label_filename, 
         "Segmentation label for initialization."},
     {NULL, ARGV_END, NULL, NULL, NULL}
};


static int usage(void)
{
    static const char msg[] = {
         "CAT_Test: Test\n"
         "usage: CAT_Test [options] -label label.nii in.nii [out.nii] []\n"
    };
    fprintf(stderr, "%s", msg);
    exit(EXIT_FAILURE);
}

int
main(int argc, char **argv)
{
    /* NIFTI stuff */
    nifti_image   *src_ptr, *label_ptr;
    char      *input_filename, *output_filename, *basename, *extension;
    int       i, j, count, dims[3], nvox;
    int       x, y, z, z_area, y_dims, modality;
    char      *arg_string, buffer[1024];
    unsigned char *label = NULL;
    float     *src, *buffer_vol;
    double    mu[3], std, min_mu, slope, voxelsize[3];
    
    /* Get arguments */
    if (ParseArgv(&argc, argv, argTable, 0) || (argc < 2)) {
        (void) fprintf(stderr, "\nUsage: %s [options] -label label.nii in.nii [out.nii]\n", argv[0]);
        (void) fprintf(stderr, "     %s -help\n\n", argv[0]);
        exit(EXIT_FAILURE);
    }
    
    input_filename  = argv[1];

    /* if not defined use original name as basename for output */
    if (argc > 2)
        output_filename = argv[2];
    else  output_filename = argv[1];
    
    /* get basename */
    basename = nifti_makebasename(output_filename);

    /* deal with extension */
    extension = nifti_find_file_extension(output_filename);
    
    /* if no valid extension was found use .nii */
    if (extension == NULL) {
        fprintf(stdout,"Use .nii as extension for %s.\n",output_filename);
        strcpy(extension, ".nii");
    }

    /* read data */
    src_ptr = read_nifti_float(input_filename, &src, 0);
    
    if (src_ptr == NULL) {
        fprintf(stderr,"Error reading %s.\n",input_filename);
        exit(EXIT_FAILURE);
    }

    nvox = src_ptr->nvox;
    label = (unsigned char *)malloc(sizeof(unsigned char)*nvox);
    
    if (label == NULL) {
        fprintf(stderr,"Memory allocation error\n");
        exit(EXIT_FAILURE);
    }

    /* read label and check for same size */
    if (label_filename == NULL) {
        fprintf(stderr,"Label image has to be defined\n");
        exit(EXIT_FAILURE);
    }
            
    /* read volume */
    label_ptr = read_nifti_float(label_filename, &buffer_vol, 0);
    if (label_ptr == NULL) {
        fprintf(stderr,"Error reading %s.\n", label_filename);
        return(EXIT_FAILURE);
    }
    
    /* check size */ 
    if (!equal_image_dimensions(src_ptr, label_ptr)) {     
        fprintf(stderr,"Label and source image have different size\n");
        exit(EXIT_FAILURE);
    }
    
    /* estimate mean for each label class */
    for (j = 0; j < 3; j++) {
        for (i = 0; i < nvox; i++)
            label[i] = (unsigned char) (round(buffer_vol[i]) == (j+1));
        mu[j] = get_masked_mean_array_float(src, nvox, label);
    }

    /* estimate image contrast (T1w/T2w/PDw) */
    if (mu[0] < mu[1] && mu[1] < mu[2]) {
        modality = 1; // T1w
    } else if (mu[0] > mu[1] && mu[1] > mu[2]) {
        modality = 2; // T2w
    } else if (mu[0] < mu[2] && mu[1] < mu[2]) {
        modality = 0; // PDw - WM is maximum
    } else {
        modality = 3; // PDw - WM is minimum or other conditions
    }
    
    /* get mean for WM normalized by std */
    for (i = 0; i < nvox; i++)
        label[i] = (unsigned char) (round(buffer_vol[i]) == 3);
    min_mu = MIN(fabs(mu[2] - mu[1]), fabs(mu[1] - mu[0]));
    std = get_masked_std_array_float(src, nvox, label) / min_mu;
    fprintf(stderr, "%g\n", std);
    
    for (i = 0; i < nvox; i++) {    
        if (label[i] == 3)
            src[i] -= mu[2];
        else src[i] = 0.0;
    }

    voxelsize[0] = src_ptr->dx;
    voxelsize[1] = src_ptr->dy;
    voxelsize[2] = src_ptr->dz;
    dims[0] = src_ptr->nx;
    dims[1] = src_ptr->ny;
    dims[2] = src_ptr->nz;
        
    /* write nu-corrected volume */
    if (1) {

        slope = 2.0/65535.0;
        sprintf(buffer, "%s_corr%s",basename,extension);
        if (!write_nifti_float(buffer, src, DT_INT16, slope, 
                        dims, voxelsize, src_ptr))
            exit(EXIT_FAILURE);
    }
    
    free(label);
    
    return(EXIT_SUCCESS);
}