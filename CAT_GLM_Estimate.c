/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 *
 */

#include <bicpl.h>
#include <time_stamp.h>

#include "CAT_Pinv.h"

#define VERBOSE 0
#define EPS 1e-15
#define DEBUG 1

/* Definitions for accessing information on each dimension */
#define MAX_FILES 1500          /* Number of possible files */

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
        int                  n_beta, rank, erdf, counter, *idx;
        int                  n_cov;
        Real                 **vals, *tmpvals, *data, **indata;
        Real                 **G, *v, **inv_G, **transp_G, **pinv_GG;
        Real                 **beta, *beta0, **estimates, **resSD, sum, *result;
        progress_struct      progress;

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
                                exit(EXIT_FAILURE);
                        }
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
                                exit(EXIT_FAILURE);
                        }
                        if (n_subj != n_tmp) {
                                fprintf(stderr, "\nError: Number of files differs from number of rows in design matrix.\n");
                                exit(EXIT_FAILURE);
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
                } else {
                        if (input_values_any_format(infiles[i], &n_vals,
                                                 &tmpvals) != OK) {
                                fprintf(stderr, "\nError reading file %s\n",
                                       infiles[i]);
                                exit(EXIT_FAILURE);
                        }
                        if (i == 0) {
                                ALLOC2D(vals, n_subj, n_vals);
                        } else {
                                if (prev_n_vals != n_vals) {
                                        fprintf(stderr, "\nError: Wrong number of values in %s\n",infiles[i]);
                                        exit(EXIT_FAILURE);
                                }
                        }
                        for (k = 0; k < n_vals; k++)
                                vals[counter][k] = tmpvals[k];
                        prev_n_vals = n_vals;
                        G[counter][j] = 1;
                        counter++;
                }

                printf(".");
                fflush(stdout);
        }

        printf("\n");
    
        /* compute pseudo inverse from design matrix */        
        rank = pinv(n_subj, n_beta, G, inv_G);

        /* effective residual d.f. */        
        erdf = n_subj - rank;

        printf("Effective residual d.f.: %d\n",erdf);

        /* check estimability */
        if (erdf < 0) {
                fprintf(stderr, "This design is unestimable! (df=%d).\n",erdf);
                exit(EXIT_FAILURE);
        }

        if (erdf == 0) {
                fprintf(stderr, "This design has no res! (df=0).\n");
                exit(EXIT_FAILURE);
        }

        /* open log file */
        if ((fp = fopen("glm.log", "w")) == 0) {
                fprintf(stderr, "Couldn't open file glm.log.\n");
                exit(EXIT_FAILURE);
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
        ALLOC(v, n_vals);
        ALLOC(beta0, n_vals);
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
        /*  estimation for ascii-files */
        /* ---------------------------------------------------------------- */

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
                for (k = 0; k < n_vals; k++) beta0[k] = beta[j][k];
    
                output_values_any_format(outfile, n_vals, beta0);
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
                        result[k] = beta[j][k] / (sqrt(v[k] * pinv_GG[j][j]) + EPS);
                }
                output_values_any_format(outfile, n_vals, result);
                FREE(result);
        }
    
        FREE(data);
        FREE(v);
        FREE(idx);
        FREE(beta0);
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
                exit(EXIT_FAILURE);
        }
        estimate(infiles, arg_string, argc);
        
        return(EXIT_SUCCESS);
}
