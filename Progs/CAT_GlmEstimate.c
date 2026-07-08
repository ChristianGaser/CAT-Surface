/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 *
 */

#include <bicpl.h>

#include "CAT_Math.h"
#include "CAT_SurfaceIO.h"
#include "CAT_NiftiLib.h"

#define VERBOSE 0
#define EPS 1e-15

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
    The input files can either be surface data in freesurfer-format,\n\
    gifti-format (using .gii) or ascii-format (using .txt), or NIfTI\n\
    volumes (using .nii or .nii.gz).  Surface results are always saved\n\
    in gifti-format (.gii) and volume results in NIfTI-format (.nii).\n\
    All input files must be of the same kind and have identical\n\
    dimensions.\n\
    The basic model at each element of the input files is of the form\n\
    Y = X*B + e, for data Y, design matrix X, (unknown) parameters B, and\n\
    residual errors e. The errors are assumed to have normal distribution.\n\
\n\
    Covariates for regression (correlation) or AnCova models can be defined\n\
    using ':' as the delimiter.  Each covariate should be defined as a\n\
    seperate file.\n\
\n\
    The following files are written:\n\
    ResMS   - residual mean square\n\
    beta_xxxx   - parameter estimates, numbered according to the\n\
                  corresponding column of the design matrix.\n\
    T_xxxx  - T-values, numbered according to the corresponding\n\
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


/* ----------------------------------------------------------------------- */
/*  helpers to abstract surface vs. volume I/O */
/* ----------------------------------------------------------------------- */

/* Return non-zero if the filename refers to a NIfTI volume (.nii/.nii.gz). */
static int
filename_is_volume(char *filename)
{
    return filename_extension_matches(filename, "nii") ||
           filename_extension_matches(filename, "gz");
}

/* Write a result vector either as gifti surface data or NIfTI volume. */
static void
write_result_values(char *outfile, int is_volume, int n_vals, double *data,
                     int dim[], double vox[], nifti_image *nii_ptr)
{
    if (is_volume) {
        if (!write_nifti_double(outfile, data, DT_FLOAT32, 1.0, dim, vox,
                                nii_ptr)) {
            fprintf(stderr, "\nError writing file %s\n", outfile);
            exit(EXIT_FAILURE);
        }
    } else {
        output_values_any_format(outfile, n_vals, data, TYPE_DOUBLE);
    }
}


/* ----------------------------------------------------------------------- */
/*  estimation of GLM */
/* ----------------------------------------------------------------------- */
int
estimate(char **infiles, int argc)
{
    char         *outfile, *ext;
    FILE         *fp;
    char         buffer[1024];
    int          i, j, k, n_subj, n_vals, n_tmp, prev_n_vals;
    int          n_beta, rank, erdf, counter, *idx;
    int          n_cov, is_volume, dim[3];
    double       **vals, *tmpvals, *data, **indata;
    double       **G, *v, **inv_G, **transp_G, **pinv_GG;
    double       **beta, *beta0, **estimates, **resSD, sum, *result;
    double       vox[3];
    nifti_image  *nii_ptr;
    progress_struct progress;

    counter = 0;
    n_beta  = 1;
    n_cov   = 0;
    is_volume = 0;
    nii_ptr = NULL;

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

    /* decide surface vs. volume mode from the first data file; all data */
    /* files are expected to be of the same kind */
    if (counter > 0)
        is_volume = filename_is_volume(infiles[idx[0]]);

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
            normalize_double(tmpvals, n_tmp);

            /* build design matrix G */
            for (k = 0; k < n_subj; k++) {
                G[k][j] = tmpvals[k];
            }
        } else if (equal_strings(infiles[i], "+")) {
            /* define groups */
            j++;
        } else {
            if (is_volume) {
                /* header only: read_nifti_double still fills tmpvals with */
                /* the voxel values and frees the raw nii_ptr->data buffer */
                nifti_image *cur_ptr = read_nifti_double(infiles[i],
                                                         &tmpvals, 0);
                if (cur_ptr == NULL) {
                    fprintf(stderr, "\nError reading file %s\n",
                           infiles[i]);
                    exit(EXIT_FAILURE);
                }
                n_vals = (int) cur_ptr->nvox;
                /* keep the first volume header as output reference */
                if (nii_ptr == NULL) {
                    nii_ptr = cur_ptr;
                    dim[0] = nii_ptr->nx;
                    dim[1] = nii_ptr->ny;
                    dim[2] = nii_ptr->nz;
                    vox[0] = nii_ptr->dx;
                    vox[1] = nii_ptr->dy;
                    vox[2] = nii_ptr->dz;
                } else {
                    /* raw data already freed inside read_nifti_double */
                    cur_ptr->data = NULL;
                    nifti_image_free(cur_ptr);
                }
            } else {
                if (input_values_any_format(infiles[i], &n_vals,
                             &tmpvals) != OK) {
                    fprintf(stderr, "\nError reading file %s\n",
                           infiles[i]);
                    exit(EXIT_FAILURE);
                }
            }

            if (counter == 0) {
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

            /* release per-file input buffer (surface data uses bicpl ALLOC, */
            /* volume data uses malloc inside read_nifti_double) */
            if (is_volume)
                free(tmpvals);
            else
                FREE(tmpvals);

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
    
    fprintf(fp, "\n[history]\n");
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

    /* output extension: gifti for surfaces, NIfTI for volumes */
    ext = is_volume ? "nii" : "gii";

    /* write betas for each column of design matrix */
    for (j = 0; j < n_beta; j++) {
        outfile = create_string("beta");
        sprintf(buffer, "_%04d.%s", j+1, ext);
        concat_to_string(&outfile, buffer);
        for (k = 0; k < n_vals; k++) beta0[k] = beta[j][k];

        write_result_values(outfile, is_volume, n_vals, beta0, dim, vox,
                             nii_ptr);
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

    sprintf(buffer, "ResMS.%s", ext);
    outfile = create_string(buffer);
    write_result_values(outfile, is_volume, n_vals, v, dim, vox, nii_ptr);

    /* write beta and beta/ResSD for each column of design matrix */
    for (j = 0; j < n_beta; j++) {
        /* ResSD values*/
        outfile = create_string("ResSD");
        sprintf(buffer, "_%04d.%s", j+1, ext);
        concat_to_string(&outfile, buffer);

        ALLOC(result, n_vals);

        for (k = 0; k < n_vals; k++) {
            result[k] = sqrt(v[k] * pinv_GG[j][j]);
        }
        write_result_values(outfile, is_volume, n_vals, result, dim, vox,
                             nii_ptr);

        /* T-values */
        outfile = create_string("T");
        sprintf(buffer, "_%04d.%s", j+1, ext);
        concat_to_string(&outfile, buffer);

        for (k = 0; k < n_vals; k++) {
            result[k] = beta[j][k] / (sqrt(v[k] * pinv_GG[j][j]) + EPS);
        }
        write_result_values(outfile, is_volume, n_vals, result, dim, vox,
                             nii_ptr);
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
    FREE2D(vals);

    if (nii_ptr != NULL) {
        /* data was freed by read_nifti_double / write_nifti_double */
        nii_ptr->data = NULL;
        nifti_image_free(nii_ptr);
    }

    return(0);
}

int
main(int argc, char *argv[])
{
    char **infiles;

    initialize_argument_processing(argc, argv);

    infiles = &argv[1];

    /* get filenames */
    if (argc < 2) {
        usage(argv[0]);
        exit(EXIT_FAILURE);
    }
    estimate(infiles, argc);
    
    return(EXIT_SUCCESS);
}
