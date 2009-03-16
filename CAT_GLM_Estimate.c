/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 *
 */

#include <volume_io/internal_volume_io.h>
#include <bicpl.h>
#include <time_stamp.h>

#include "CAT_Pinv.h"

#define VERBOSE 0
#define EPS 1e-15
#define DEBUG 1

/* Definitions for accessing information on each dimension */
#define NUM_DIMENSIONS 3        /* Number of dimensions of volume */
#define MAX_FILES 1500          /* Number of possible files */

/* Types */
typedef struct {
   long nslices;                     /* Number of slices */
   long nrows;                       /* Number of rows in inputdata */
   long ncolumns;                    /* Number of columns in inputdata */
   double maximum;                   /* Volume maximum */
   double minimum;                   /* Volume minimum */
   double step[NUM_DIMENSIONS];    /* Step sizes for dimensions */
   double start[NUM_DIMENSIONS];    /* Start positions for dimensions */
   char dimension_names[NUM_DIMENSIONS][MAX_NC_NAME]; /* Dimension names */
} Volume_Info;

void
usage(char *executable)
{
        char *usage_str = "\n\
NAME\n\
    CAT_GLM_Estimate - estimation of a General Linear Model (GLM)\n\n\
SYNOPSIS\n\
    CAT_GLM_Estimate file1_grp1 file2_grp1 ... + file1_grp2 file2_grp2 ... : covariate_file \n\n\
DESCRIPTION\n\
    CAT_GLM_Estimate estimates the beta parameters for the GLM and\n\
    writes the resulting parameters and the t-values for each column\n\
    of the design matrix to the current working directory.  The files\n\
    are named beta_xxxx and T_xxxx, where xxxx are numbered according\n\
    to the corresponding column of the design matrix. Additionally the\n\
    residual standard deviation is saved in a file ResSD_xxxx.\n\
\n\
    The input files can be either in minc-format (using the extension\n\
    .mnc) or in ascii-format (using .txt as extension) and the respective\n\
    result files are saved in the same format and with the same extension.\n\
    The basic model at each element of the input files is of the form\n\
    Y = X*B + e, for data Y, design matrix X, (unknown) parameters B, and\n\
    residual errors e. The errors are assumed to have normal distribution.\n\
\n\
    Covariates for regression (correlation) or AnCova models can be defined\n\
    using ':' as the delimiter.  Each covariate should be defined as a\n\
    seperate file.\n\
\n\
    The following files are written:\n\
    ResMS       - residual mean square\n\
    beta_xxxx   - parameter estimates, numbered according to the\n\
                                  corresponding column of the design matrix.\n\
    T_xxxx      - T-values, numbered according to the corresponding\n\
                                  column of the design matrix.  These values\n\
                                  are the parameter estimates divided by\n\
                                  residual standard deviation.\n\
    ResSD_xxxx  - estimated residual standard deviation, numbered according\n\
                                  to the corresponding column of the design\n\
                                  matrix.\n\n\
    These files can be used with glm_mat to calculate different contrasts\n\
    of factor levels.\n";

        fprintf(stderr, usage_str, executable);
}

void
normalize_vector(Real *data, int length)
{
        Real mean;
        int i;
    
        mean = 0.0;
        for (i = 0; i < length; i++)
                mean += data[i];
        mean = mean / length;

        for (i = 0; i < length; i++)
                data[i] -= mean;

}

/* ----------------------------------------------------------------------- */
/*  estimation of GLM */
/* ----------------------------------------------------------------------- */
int
estimate(char **infiles, char *arg_string, int argc)
{
        char                 *outfile;
        FILE                 *fp;
        char                 buffer[1024];
        int                  i, j, k, n_subj, n_vals, n_tmp, prev_n_vals;
        int                  n_beta, rank, erdf, counter;
        int                  in_icvid[MAX_FILES], out_beta_icvid[MAX_FILES];
        int                  out_T_icvid[MAX_FILES], out_res_icvid[MAX_FILES];
        int                  out_resMS_icvid, *idx;
        int                  islice, n_dims, n_cov;
        double               *outdata;
        Real                 **vals, *tmpvals, *data, **indata;
        Real                 **G, *v, **inv_G, **transp_G, **pinv_GG;
        Real                 **beta, **estimates, **resSD, sum, *result;
        progress_struct      progress;
        Volume_Info          vol_info[MAX_FILES];

        if (filename_extension_matches(infiles[0], "mnc"))
                n_dims = 3;
        else    n_dims = 1;

        counter = 0;
        n_beta = 1;
        n_cov = 0;
    
        ALLOC(idx, argc-1);

        for (i = 0; i < argc-1; i++) {
                if (equal_strings(infiles[i], "+")) {
                        n_beta++;
                } else if (equal_strings(infiles[i], ":")) {
                        n_beta++;
                        n_cov++;
                        /* skip file containing covariates by incr. i */
                        i++;
                } else {
                        idx[counter] = i;
                        /* check for files */
                        if (!file_exists(infiles[i])) {
                                fprintf(stderr, "\nFile %s not found.\n",
                                        infiles[i]);
                                return(1);
                        }
                        /* read input icvid for minc files */
                        if (n_dims == 3)
                                in_icvid[counter] = get_volume_info(infiles[i],
                                                                   &vol_info[i],
                                                                   NC_DOUBLE);
                        counter++;
                }
        }
    
        n_subj = argc - n_cov - n_beta;
    
        ALLOC2D(G, n_subj, n_beta);
        ALLOC2D(inv_G, n_beta, n_subj);

        /* initialize design matrix G to zero */
        for (j = 0; j < n_beta; j++) {
                for (i = 0; i < n_subj; i++)
                        G[i][j] = 0;     
        }

        /* read data and build design matrix */
        if (n_dims == 1)
                printf("Read files:");
        fflush(stdout);
    
        counter = 0;
        for (i = 0, j = 0; i < argc-1; i++) {
                /* define covariates */
                if (equal_strings(infiles[i], ":")) {
                        /* count columns and files */
                        i++; j++;
                        if (input_values_any_format(infiles[i], &n_tmp,
                                                 &tmpvals) != OK) {
                                fprintf(stderr, "\nError reading file %s\n",
                                        infiles[i]);
                                return(1);
                        }
                        if (n_subj != n_tmp) {
                                fprintf(stderr, "\nError: Number of files differs from number of rows in design matrix.\n");
                                return(1);
                        }
                        /* normalize covariates to zero mean */
                        normalize_vector(tmpvals, n_tmp);

                        /* build design matrix G */
                        for (k = 0; k < n_subj; k++) {
                                G[k][j] = tmpvals[k];
                        }
                } else if (equal_strings(infiles[i], "+")) {
                        /* define groups */
                        j++;
                } else if (n_dims == 1) {
                        if (input_values_any_format(infiles[i], &n_vals,
                                                 &tmpvals) != OK) {
                                fprintf(stderr, "\nError reading file %s\n",
                                       infiles[i]);
                                return(1);
                        }
                        if (i == 0) {
                                ALLOC2D(vals, n_subj, n_vals);
                        } else {
                                if (prev_n_vals != n_vals) {
                                        fprintf(stderr, "\nError: Wrong number of values in %s\n",infiles[i]);
                                        return(1);
                                }
                        }
                        for (k = 0; k < n_vals; k++)
                                vals[counter][k] = tmpvals[k];
                        prev_n_vals = n_vals;
                        G[counter][j] = 1;
                        counter++;
                } else if (n_dims == 3) {
                        G[counter][j] = 1;
                        counter++;
                }

                if (n_dims == 1) {
                        printf(".");
                        fflush(stdout);
                }
        }

        if (n_dims == 1)
                printf("\n");

        if (n_dims == 3)
                n_vals = vol_info[0].nrows * vol_info[0].ncolumns;
    
        /* compute pseudo inverse from design matrix */        
        rank = pinv(n_subj, n_beta, G, inv_G);

        /* effective residual d.f. */        
        erdf = n_subj - rank;

        printf("Effective residual d.f.: %d\n",erdf);

        /* check estimability */
        if (erdf < 0) {
                fprintf(stderr, "This design is unestimable! (df=%d).\n",erdf);
                return(1);
        }

        if (erdf == 0) {
                fprintf(stderr, "This design has no res! (df=0).\n");
                return(1);
        }

        /* open log file */
        if ((fp = fopen("glm.log", "w")) == 0) {
                fprintf(stderr, "Couldn't open file glm.log.\n");
                return(1);
        }

        fprintf(fp, "[df]\n%d\n", erdf);
        fprintf(fp, "\n[design matrix]\n");

        for (i = 0; i < n_subj; i++) {
                if (VERBOSE)
                        printf("%s\t", infiles[idx[i]]);
                fprintf(fp, "%s\t", infiles[idx[i]]);          
                for (j = 0; j < n_beta; j++) {
                        if (VERBOSE)
                                printf("%6.3f ", G[i][j]);
                        fprintf(fp, "%6.3f ", G[i][j]);
                }
                if (VERBOSE)
                        printf("\n");
                fprintf(fp, "\n");
        }
    
        fprintf(fp, "\n[history]\n%s\n", arg_string);
        fclose(fp);

        ALLOC(data, n_vals);
        ALLOC(outdata, n_vals);
        ALLOC(v, n_vals);
        ALLOC2D(resSD, n_beta, n_vals);
        ALLOC2D(indata, n_subj, n_vals);
        ALLOC2D(estimates, n_subj, n_vals);
        ALLOC2D(beta, n_beta, n_vals);
        ALLOC2D(pinv_GG, n_beta, n_subj);
        ALLOC2D(transp_G, n_beta, n_subj);

        /* calculate pinv(G'*G) */
        transpose(n_subj, n_beta, G, transp_G);
        matrix_multiply(n_beta, n_subj, n_beta, transp_G, G, pinv_GG);
        (void) pinv(n_beta, n_beta, pinv_GG, pinv_GG);

        /* ---------------------------------------------------------------- */
        /*  estimation for minc-files */
        /* ---------------------------------------------------------------- */
        if (n_dims == 3) {
                /* Check vol info about dimensions, voxel size and origin */
                for (i = 1; i < n_subj; i++) {
                        if (vol_info[idx[i]].nslices != vol_info[0].nslices ||
                            vol_info[idx[i]].nrows != vol_info[0].nrows ||
                            vol_info[idx[i]].ncolumns != vol_info[0].ncolumns) {
                                fprintf(stderr,
                                        "Dimensions of file %s differ.\n",
                                        infiles[idx[i]]);
                                return(1);
                        }
                        if (vol_info[idx[i]].step[0] != vol_info[0].step[0] ||
                            vol_info[idx[i]].step[1] != vol_info[0].step[1] ||
                            vol_info[idx[i]].step[2] != vol_info[0].step[2]) {
                                fprintf(stderr,
                                        "Voxel size of file %s differs.\n",
                                        infiles[idx[i]]);
                                return(1);
                        }
                        if (vol_info[idx[i]].start[0] != vol_info[0].start[0] ||
                            vol_info[idx[i]].start[1] != vol_info[0].start[1] ||
                            vol_info[idx[i]].start[2] != vol_info[0].start[2]) {
                                fprintf(stderr,
                                        "Origin of file %s differs.\n",
                                        infiles[idx[i]]);
                                return(1);
                        }
                }

                /* prepare output for resMS */
                outfile = create_string("ResMS.mnc");
                out_resMS_icvid = save_volume_info(in_icvid[0], outfile,
                                                   arg_string, &vol_info[0],
                                                   NC_FLOAT);

                /* prepare output files for beta-, T-, and ResSD-images */
                for (j = 0; j < n_beta; j++) {
                        outfile = create_string("beta");
                        sprintf(buffer, "_%04d.mnc", j+1);
                        concat_to_string(&outfile, buffer);
                        out_beta_icvid[j] = save_volume_info(in_icvid[0],
                                                             outfile,
                                                             arg_string, 
                                                             &vol_info[0],
                                                             NC_FLOAT);
                        outfile = create_string("T");
                        sprintf(buffer, "_%04d.mnc", j+1);
                        concat_to_string(&outfile, buffer);
                        out_T_icvid[j] = save_volume_info(in_icvid[0], outfile,
                                                          arg_string, 
                                                          &vol_info[0],
                                                          NC_FLOAT);
                        outfile = create_string("ResSD");
                        sprintf(buffer, "_%04d.mnc", j+1);
                        concat_to_string(&outfile, buffer);
                        out_res_icvid[j] = save_volume_info(in_icvid[0],
                                                            outfile, arg_string,
                                                            &vol_info[0],
                                                            NC_FLOAT);
                }

                initialize_progress_report(&progress, FALSE,
                                           vol_info[0].nslices, "Estimate GLM");

                /* Loop through slices */
                for (islice = 0; islice < vol_info[0].nslices; islice++) {
                        /* get input data */    
                        for (i = 0; i < n_subj; i++) {
                                get_volume_slice_double(in_icvid[i],
                                                        &vol_info[idx[i]],
                                                        islice, data);
                                for (k = 0; k < n_vals; k++)
                                        indata[i][k] = data[k];
                        }

                        /* calculate beta by multiplying pseudo inverse */
                        /* from design matrix with transposed data */    
                        matrix_multiply(n_beta, n_subj, n_vals, inv_G,
                                        indata, beta);

                        /* calculate fitted data: estimates = G*beta */
                        matrix_multiply(n_subj, n_beta, n_vals, G, beta,
                                        estimates);

                        /* calculate estimated residual standard deviation: */
                        /* ResSD = sqrt(sum(res.^2)/erdf) */
                        for (k = 0; k < n_vals; k++) {
                                sum = 0.0;
                                for (i = 0; i < n_subj; i++) {
                                        estimates[i][k] = indata[i][k] -
                                                          estimates[i][k];
                                        sum += estimates[i][k] *
                                               estimates[i][k];
                                }
                                v[k] = sum / (Real) erdf;
                        }

                        for (k = 0; k < n_vals; k++)
                                outdata[k] = (double) v[k];
                        save_volume_slice_double(out_resMS_icvid, &vol_info[0],
                                                 islice, outdata);

                        for (j = 0; j < n_beta; j++) {
                                for (k = 0; k < n_vals; k++)
                                        resSD[j][k] = sqrt(v[k] *
                                                           pinv_GG[j][j]);
                        }

                        /* save slices */
                        for (j = 0; j < n_beta; j++) {
                                for (k = 0; k < n_vals; k++)
                                        outdata[k] = (double) resSD[j][k];

                                save_volume_slice_double(out_res_icvid[j],
                                                         &vol_info[0], islice,
                                                         outdata);
                                for (k = 0; k < n_vals; k++)
                                        outdata[k] = (double) beta[j][k];

                                save_volume_slice_double(out_beta_icvid[j],
                                                         &vol_info[0], islice,
                                                         outdata);
                                for (k = 0; k < n_vals; k++)
                                        outdata[k] = (double) beta[j][k] /
                                                     (resSD[j][k] + EPS);

                                save_volume_slice_double(out_T_icvid[j],
                                                         &vol_info[0], islice,
                                                         outdata);
                        }
                        update_progress_report(&progress, islice + 1);    
                }  // islice++
                terminate_progress_report(&progress);

                /* Free up data and close files */
                close_volume(out_resMS_icvid);
                for (i = 0; i < n_subj; i++)
                        close_volume(in_icvid[i]);

                for (j = 0; j < n_beta; j++) {
                        close_volume(out_beta_icvid[j]);
                        close_volume(out_T_icvid[j]);
                }
        }  // if n_dims == 3

        /* ---------------------------------------------------------------- */
        /*  estimation for ascii-files */
        /* ---------------------------------------------------------------- */

        if (n_dims == 1) {
                /* multiply pseudo inverse from design matrix with */
                /* transposed data */    
                for (k = 0; k < n_vals; k++) {
                        for (j = 0; j < n_beta; j++) {
                                sum = 0.0;
                                for (i = 0; i < n_subj; i++)
                                sum += inv_G[j][i] * vals[i][k];        
                                beta[j][k] = sum;
                        }
                }

                /* write betas for each column of design matrix */
                for (j = 0; j < n_beta; j++) {
                        outfile = create_string("beta");
                        sprintf(buffer, "_%04d.txt", j+1);
                        concat_to_string(&outfile, buffer);
    
                        output_values_any_format(outfile, n_vals, beta);
                }

                /* calculate fitted data: estimates = G*beta */
                matrix_multiply(n_subj, n_beta, n_vals, G, beta, estimates);

                /* calculate estimated squared residual standard deviation: */
                /* v = sum((vals - estimates).^2)/erdf */
                for (k = 0; k < n_vals; k++) {
                        sum = 0.0;
                        for (i = 0; i < n_subj; i++) {
                                estimates[i][k] = vals[i][k] - estimates[i][k];
                                sum += estimates[i][k] * estimates[i][k];
                        }
                        v[k] = sum / erdf;
                }

                outfile = create_string("ResMS.txt");    
                output_values_any_format(outfile, n_vals, v);

                /* write beta and beta/ResSD for each column of design matrix */
                for (j = 0; j < n_beta; j++) {
                        /* ResSD values*/
                        outfile = create_string("ResSD");
                        sprintf(buffer, "_%04d.txt", j+1);
                        concat_to_string(&outfile, buffer);
                        
                        ALLOC(result, n_vals);
                            
                        for (k = 0; k < n_vals; k++) {
                                result[k] = sqrt(v[k] * pinv_GG[j][j]);
                        }
                        output_values_any_format(outfile, n_vals, result);

                        /* T-values */
                        outfile = create_string("T");
                        sprintf(buffer, "_%04d.txt", j+1);
                        concat_to_string(&outfile, buffer);

                        for (k = 0; k < n_vals; k++) {
                                result[k] = beta[j][k] /
                                         (sqrt(v[k] * pinv_GG[j][j]) + EPS);
                        }
                        output_values_any_format(outfile, n_vals, result);
                        FREE(result);
                }
        } // if n_dims == 1
    
        FREE(data);
        FREE(outdata);
        FREE(v);
        FREE(idx);
        FREE2D(resSD);
        FREE2D(indata);
        FREE2D(estimates);
        FREE2D(transp_G);
        FREE2D(pinv_GG);
        FREE2D(beta);
   
        return(0);
}

int
main(int argc, char *argv[])
{
        char *arg_string, **infiles;

        /* Save time stamp and args */
        arg_string = time_stamp(argc, argv);
    
        initialize_argument_processing(argc, argv);

        infiles = &argv[1];

        /* get filenames */
        if (argc < 2) {
                usage(argv[0]);
                return(1);
        }
        estimate(infiles, arg_string, argc);
}
