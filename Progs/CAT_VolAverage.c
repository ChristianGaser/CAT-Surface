/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 *
 */

#include "ParseArgv.h"
#include "CAT_NiftiLib.h"

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include "CAT_SafeAlloc.h"

int verbose = 0;

extern int ParseArgv(int *argcPtr, char **argv, ArgvInfo *argTable, int flags);

char *std_filename = NULL;
char *zscore_filename = NULL;
char *zscore_txt_filename = NULL;

static ArgvInfo argTable[] = {
    {"-std", ARGV_STRING, (char *) 1, (char *) &std_filename, 
          "Write standard deviation."},
    {"-zscore", ARGV_STRING, (char *) 1, (char *) &zscore_filename, 
          "Write filenames and z-scores in csv-file."},
    {"-zscore-txt", ARGV_STRING, (char *) 1, (char *) &zscore_txt_filename, 
          "Write z-scores in txt-file."},
    {"-v", ARGV_CONSTANT, (char *) 1, (char *) &verbose,
          "Be verbose."},
     {NULL, ARGV_END, NULL, NULL, NULL}
};

/* Main program */

int main(int argc, char *argv[])
{
    char **infiles, *outfile;
    int i, j, nfiles, dims[3], n_voxels;
    double *avg, *sum_squares, *input, voxelsize[3], *zscore;
    nifti_image *nii_ptr, *nii_ptr2;
    FILE *fid;

    /* Get arguments */
    if (ParseArgv(&argc, argv, argTable, 0) || (argc < 2)) {
        (void) fprintf(stderr, 
        "\nUsage: %s [-std std-out.nii] [-zscore zscore.csv] in1.nii ... inx.nii out.nii\n", argv[0]);
        (void) fprintf(stderr, 
        "     %s -help\n\n", argv[0]);
        exit(EXIT_FAILURE);
    }
    
    nfiles = argc - 2;
    infiles = &argv[1];
    outfile = argv[argc-1];

    /* Make sure that we have something to process */
    if (nfiles == 0) {
        (void) fprintf(stderr, "Error: No input files specified\n");
        exit(EXIT_FAILURE);
    }
    
    /* read first image to get image parameters */
    nii_ptr = read_nifti_double(infiles[0], &input, 0);
    if(nii_ptr == NULL) {
        fprintf(stderr,"Error reading %s.\n", infiles[0]);
        return(EXIT_FAILURE);
    }
    fprintf(stdout,"%3d: %s\n",1, infiles[0]);

    voxelsize[0] = nii_ptr->dx;
    voxelsize[1] = nii_ptr->dy;
    voxelsize[2] = nii_ptr->dz;
    dims[0] = nii_ptr->nx;
    dims[1] = nii_ptr->ny;
    dims[2] = nii_ptr->nz;

    /* prepare average image */
    avg = (double *)malloc(sizeof(double)*nii_ptr->nvox);
    for (j=0; j<nii_ptr->nvox; j++) 
        avg[j] = 0.0;
        
    /* prepare sum of squares image */
    if ((std_filename != NULL) || (zscore_filename != NULL) || (zscore_txt_filename != NULL)) {
        sum_squares = (double *)malloc(sizeof(double)*nii_ptr->nvox);
        for (j=0; j<nii_ptr->nvox; j++) 
            sum_squares[j] = 0.0;
    }        
    
    /* prepare zscore */
    if ((zscore_filename != NULL) || (zscore_txt_filename != NULL))
        zscore = (double *)malloc(sizeof(double)*nfiles);
        
    free(input);
    
    /* read images and check for image parameters */
    for (i=0; i<nfiles; i++) {
        fprintf(stdout,"run1 %5d/%d:\t%s\n",i+1, nfiles, infiles[i]);
        nii_ptr2 = read_nifti_double(infiles[i], &input, 0);
        if(nii_ptr2 == NULL) {
            fprintf(stderr,"Error reading %s.\n", infiles[i]);
            return(EXIT_FAILURE);
        }
        
        /* check for dimensions */
        if (verbose)
            equal_image_dimensions(nii_ptr,nii_ptr2);
        
        /* calculate average */
        for (j=0; j<nii_ptr->nvox; j++) 
            avg[j] += (input[j]/(double)nfiles);

        /* calculate sum of squares */
        if ((std_filename != NULL) || (zscore_filename != NULL) || (zscore_txt_filename != NULL)) {
            for (j=0; j<nii_ptr->nvox; j++) 
                sum_squares[j] += input[j]*input[j];
        }
        free(input);
    }

    if ((std_filename != NULL) || (zscore_filename != NULL) || (zscore_txt_filename != NULL)) {
        for (j=0; j<nii_ptr->nvox; j++) 
            sum_squares[j] = sqrt(1.0/((double)nfiles-1.0)*(sum_squares[j] - (double)nfiles*avg[j]*avg[j]));
            
        if (std_filename != NULL) {
            if (!write_nifti_double( std_filename, sum_squares, DT_FLOAT32, 1.0, dims, voxelsize, nii_ptr)) 
                exit(EXIT_FAILURE);
        }
        
        if ((zscore_filename != NULL) || (zscore_txt_filename != NULL)) {
            /* compute global mean threshold (to exclude background) */
            double mean_nz = 0.0, mean_ab = 0.0, global_thresh, thresh8, z;
            int cnt_nz = 0, cnt_ab = 0;
            for (j = 0; j < nii_ptr->nvox; j++) {
                if (avg[j] != 0.0) { mean_nz += avg[j]; cnt_nz++; }
            }
            if (cnt_nz > 0) mean_nz /= (double)cnt_nz;
            thresh8 = mean_nz / 8.0;
            for (j = 0; j < nii_ptr->nvox; j++) {
                if (avg[j] > thresh8) { mean_ab += avg[j]; cnt_ab++; }
            }
            if (cnt_ab > 0) mean_ab /= (double)cnt_ab;
            global_thresh = 0.25 * mean_ab;

            for (i=0; i<nfiles; i++) {
                fprintf(stdout,"run2 %5d/%d:\t%s\n",i+1, nfiles, infiles[i]);
                nii_ptr = read_nifti_double(infiles[i], &input, 0);
                if(nii_ptr == NULL) {
                    fprintf(stderr,"Error reading %s.\n", infiles[i]);
                    return(EXIT_FAILURE);
                }

                zscore[i] = 0.0;
                n_voxels = 0;
                for (j=0; j<nii_ptr->nvox; j++) {
                    if ((sum_squares[j]==0.0) || (isnan(sum_squares[j]))
                        || avg[j] <= global_thresh)
                        continue;
                    z = (input[j] - avg[j]) / sum_squares[j];
                    zscore[i] += z * z * z * z;
                    n_voxels++;
                }
                if (n_voxels > 0)
                    zscore[i] = pow(zscore[i] / (double)n_voxels, 0.25);
                free(input);
            }
            if (zscore_filename != NULL) {
                fid = SAFE_FOPEN(zscore_filename, "w");
                fprintf(fid,"filename,z-score\n");
                for (i=0; i<nfiles; i++)
                    fprintf(fid,"%s,%g\n",infiles[i],zscore[i]);
                fclose(fid);
            }
            if (zscore_txt_filename != NULL) {
                fid = SAFE_FOPEN(zscore_txt_filename, "w");
                for (i=0; i<nfiles; i++)
                    fprintf(fid,"%g\n",zscore[i]);
                fclose(fid);
            }

            free(zscore);
        }        
        free(sum_squares);
    }

    if (!write_nifti_double( outfile, avg, DT_FLOAT32, 1.0, dims, voxelsize, nii_ptr)) 
        exit(EXIT_FAILURE);

    free(avg);
    
    return(EXIT_SUCCESS);

}
