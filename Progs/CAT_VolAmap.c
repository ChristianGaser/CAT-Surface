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
#include "CAT_Amap.h"
#include "CAT_NiftiLib.h"

char *label_filename = NULL;
int n_pure_classes = 3;
int iters_amap = 200;
int subsample = 32;
int iters_ICM = 50;
int pve = 5;
int write_seg[3] = {0, 1, 0};
int write_label = 1;
int debug = 0;
double weight_MRF = 0.15;

static ArgvInfo argTable[] = {
        {"-label", ARGV_STRING, (char *) 1, (char *) &label_filename, 
                         "Segmentation label for initialization."},
        {"-iters", ARGV_INT, (char *) 1, (char *) &iters_amap,
                         "Number of iterations to end."},
        {"-sub", ARGV_INT, (char *) 1, (char *) &subsample,
                         "Subsampling for Amap approach."},
        {"-iters_icm", ARGV_INT, (char *) 1, (char *) &iters_ICM,
                         "Number of iterations for Iterative Conditional Mode (ICM)."},
        {"-mrf", ARGV_FLOAT, (char *) 1, (char *) &weight_MRF,
                         "Weight of MRF prior (0..1)."},
        {"-pve", ARGV_INT, (char *) 1, (char *) &pve,
                         "Use Partial Volume Estimation with 5 classes (5), 6 classes (6) or do not use PVE (0)."},
        {"-write_seg", ARGV_INT, (char *) 3, (char *) &write_seg,
                         "Write fuzzy segmentations as separate images. Three numbers should be given, while a '1' indicates that this tissue class should be saved. Order is CSF/GM/WM."},
        {"-write_label", ARGV_CONSTANT, (char *) 1, (char *) &write_label,
                         "Write label image (default)."},
        {"-nowrite_label", ARGV_CONSTANT, (char *) 0, (char *) &write_label,
                         "Do not write label image."},
        {"-debug", ARGV_CONSTANT, (char *) 1, (char *) &debug,
                         "Print debug information."},
         {NULL, ARGV_END, NULL, NULL, NULL}
};


static int usage(void)
{
        static const char msg[] = {
                        "CAT_VolAmap: Segmentation with adaptive MAP\n"
                        "usage: CAT_VolAmap [options] in.nii [out.nii]\n"
        };
        fprintf(stderr, "%s", msg);
        exit(EXIT_FAILURE);
}

int
main( int argc, char **argv )
{
        /* NIFTI stuff */
        nifti_image   *src_ptr, *label_ptr;
        int           n_classes;
        char          *input_filename, *output_filename, *basename, *extension;
        int           i, j, dims[3];
        int           x, y, z, z_area, y_dims;
        char          *arg_string, buffer[1024];
        unsigned char *label, *prob;
        float         *src, *buffer_vol;
        double        slope;
        double        offset, val, max_vol, min_vol, voxelsize[3];

        /* Get arguments */
        if (ParseArgv(&argc, argv, argTable, 0) || (argc < 2)) {
                (void) fprintf(stderr, "\nUsage: %s [options] in.nii [out.nii]\n", argv[0]);
                (void) fprintf(stderr, "         %s -help\n\n", argv[0]);
                exit(EXIT_FAILURE);
        }
        
        input_filename  = argv[1];

        /* if not defined use original name as basename for output */
        if (argc == 3)
                output_filename = argv[2];
        else    output_filename = argv[1];
        
        /* get basename */
        basename = nifti_makebasename(output_filename);

        /* deal with extension */
        extension = nifti_find_file_extension(output_filename);
        
        /* if no valid extension was found use .nii */
        if (extension == NULL) {
                fprintf(stdout,"Use .nii as extension for %s.\n",output_filename);
                strcpy( extension, ".nii");
        }

        if (iters_amap == 0) {
                fprintf(stdout,"To estimate segmentation you need at least one iteration.\n");
                exit(EXIT_FAILURE);
        }
        
        if ((pve != 0) && (pve != 5) && (pve != 6)) {
                fprintf(stderr,"Value for pve can be either 0 (no PVE), 5 or 6 (5 or 6 classes).\n");
                exit(EXIT_FAILURE);
        }

        switch(pve) {
        case 0:
                n_classes = 3;
                break;
        case 5:
                n_classes = 5;
                break;
        case 6:
                n_classes = 6;
                break;
        }
        
        /* read data */
        src_ptr = read_nifti_float(input_filename, &src, 0);
        
        if (src_ptr == NULL) {
                fprintf(stderr,"Error reading %s.\n",input_filename);
                exit(EXIT_FAILURE);
        }

        label = (unsigned char *)malloc(sizeof(unsigned char)*src_ptr->nvox);
        prob    = (unsigned char *)malloc(sizeof(unsigned char)*src_ptr->nvox*n_classes);
        
        if ((label == NULL) || (prob == NULL)) {
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
        
        /* scale label image to a range 0..3 */
        for (i = 0; i < label_ptr->nvox; i++) 
                label[i] = (unsigned char) round(buffer_vol[i]);


        double mean[n_classes], mu[n_pure_classes], var[n_pure_classes];
        for (i = 0; i < n_pure_classes; i++)
                mu[i] = 0;

        /* get min/max */
        min_vol =        FLT_MAX; max_vol = -FLT_MAX;
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

        /* add offset to ensure that CSF values are much larger than background noise */
        offset = 0.2*max_vol;    
        for (i = 0; i < src_ptr->nvox; i++)
                if (label[i] > 0) src[i] += (float)offset;

        voxelsize[0] = src_ptr->dx;
        voxelsize[1] = src_ptr->dy;
        voxelsize[2] = src_ptr->dz;
        dims[0] = src_ptr->nx;
        dims[1] = src_ptr->ny;
        dims[2] = src_ptr->nz;

        Amap( src, label, prob, mean, n_pure_classes, iters_amap, subsample, dims, pve, weight_MRF, voxelsize, iters_ICM, offset);

        /* PVE */
        if (pve) {
                fprintf(stdout,"Calculate Partial Volume Estimate.\n");
                if (pve==6)
                        Pve6(src, prob, label, mean, dims);
                else
                        Pve5(src, prob, label, mean, dims);
        }
                
        /* write labeled volume */
        if (write_label) {

                /* different ranges for pve */
                if (pve) slope = 3.0/255.0;
                else slope = 1.0;

                for (i = 0; i < src_ptr->nvox; i++)
                        src[i] = (float)label[i];

                (void) sprintf( buffer, "%s_seg%s",basename,extension); 

                if (!write_nifti_float(buffer, src, DT_UINT8, slope, 
                                                dims, voxelsize, src_ptr))
                        exit(EXIT_FAILURE);
                
        }
        
        /* write fuzzy segmentations for each class */
        if (write_seg[0] || write_seg[1] || write_seg[2]) {
                
                slope = 1.0/255.0;

                for (j = 0; j < n_pure_classes; j++) {
                        if (write_seg[j]) {
                                (void) sprintf( buffer, "%s_prob%d%s",basename,j,extension); 
                                
                                for (i = 0; i < src_ptr->nvox; i++)
                                        src[i] = prob[i+(j*src_ptr->nvox)];
                                
                                if (!write_nifti_float( buffer, src, DT_UINT8, slope, 
                                                                dims, voxelsize, src_ptr))
                                        exit(EXIT_FAILURE);
                        
                        }
                }                
        }
        
        free(src);
        free(prob);
        free(label);
        
        return(EXIT_SUCCESS);
}