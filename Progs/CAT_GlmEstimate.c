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
#include "CAT_GlmFormula.h"

#define VERBOSE 0

void
usage(char *executable)
{
    char *usage_str = "\n\
NAME\n\
    CAT_GLM_Estimate - estimation of a General Linear Model (GLM)\n\n\
SYNOPSIS\n\
    CAT_GLM_Estimate file1_grp1 file2_grp1 ... + file1_grp2 file2_grp2 ... : covariate_file \n\
    CAT_GLM_Estimate scan1 scan2 ... scanX -formula \"~ group.csv + age.csv\"\n\n\
DESCRIPTION\n\
    CAT_GLM_Estimate estimates the beta parameters for the GLM and\n\
    writes the resulting parameter estimates for each column of the\n\
    design matrix to the current working directory.  The files are\n\
    named beta_xxxx, where xxxx are numbered according to the\n\
    corresponding column of the design matrix.  Additionally the\n\
    residual mean square is saved in a file ResMS (as used by SPM12).\n\
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
    The design matrix can be defined in two ways:\n\
\n\
    1) Positionally: groups are separated by '+' and each covariate is\n\
       introduced with ':' followed by a file of one value per scan.\n\
\n\
    2) With an R-style model formula via -formula.  The scans are listed\n\
       first, then a formula whose right-hand side names one file per\n\
       variable (one value per scan).  Numeric files become covariates,\n\
       text files become factors; factor(file) forces a factor.  The\n\
       operators '+', ':' and '*' are supported, and a leading '0' or\n\
       '-1' drops the intercept.  Examples:\n\
         -formula \"~ group.csv\"                    group comparison\n\
         -formula \"~ age.csv\"                      regression\n\
         -formula \"~ group.csv * age.csv\"          group x covariate\n\
         -formula \"~ group.csv * sex.csv * site.csv\"  three-way ANOVA\n\
\n\
    Formula limitations: '*' and ':' cannot be mixed within a single\n\
    term (e.g. a:b*c); only treatment contrasts (with intercept) and\n\
    cell-means coding (leading '0'/'-1') are supported - there is no\n\
    poly(), spline, I(), nesting ('/') or '%%in%%'.  Each variable file\n\
    must hold exactly one value per scan and is matched to the scans by\n\
    row order.\n\
\n\
    The following files are written:\n\
    beta_xxxx   - parameter estimates, numbered according to the\n\
                  corresponding column of the design matrix.\n\
    ResMS   - residual mean square (as used by SPM12).\n\n\
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
/*  shared estimation core                                                 */
/*                                                                         */
/*  Given the ordered list of scan files and a ready-made design matrix G  */
/*  (n_subj x n_beta), estimate the betas and residual mean square and     */
/*  write beta_xxxx and ResMS.  The design matrix is not freed here.       */
/* ----------------------------------------------------------------------- */
static int
estimate_core(char **scan_files, int n_subj, double **G, int n_beta,
              char **colnames)
{
    char         *outfile, *ext, buffer[1024];
    FILE         *fp;
    int          i, j, k, n_vals = 0, prev_n_vals = 0, rank, erdf;
    int          is_volume, dim[3];
    double       vox[3], sum;
    double       **vals = NULL, *tmpvals, *v, *beta0;
    double       **inv_G, **beta, **estimates;
    nifti_image  *nii_ptr = NULL;

    is_volume = filename_is_volume(scan_files[0]);

    /* pseudo inverse of the design matrix and effective residual d.f. */
    ALLOC2D(inv_G, n_beta, n_subj);
    rank = pinv(n_subj, n_beta, G, inv_G);
    erdf = n_subj - rank;
    printf("Effective residual d.f.: %d\n", erdf);
    if (erdf < 0) {
        fprintf(stderr, "This design is unestimable! (df=%d).\n", erdf);
        FREE2D(inv_G);
        return 0;
    }
    if (erdf == 0) {
        fprintf(stderr, "This design has no res! (df=0).\n");
        FREE2D(inv_G);
        return 0;
    }

    /* log the design matrix (with column labels when available) */
    if ((fp = fopen("glm.log", "w")) == NULL) {
        fprintf(stderr, "Couldn't open file glm.log.\n");
        FREE2D(inv_G);
        return 0;
    }
    fprintf(fp, "[df]\n%d\n", erdf);
    fprintf(fp, "\n[design matrix]\n");
    if (colnames != NULL) {
        fprintf(fp, "%-40s", "");
        for (j = 0; j < n_beta; j++)
            fprintf(fp, " %14s", colnames[j]);
        fprintf(fp, "\n");
    }
    for (i = 0; i < n_subj; i++) {
        fprintf(fp, "%-40s", scan_files[i]);
        for (j = 0; j < n_beta; j++)
            fprintf(fp, " %14.4f", G[i][j]);
        fprintf(fp, "\n");
    }
    fprintf(fp, "\n[history]\n");
    fclose(fp);

    /* read the scan data (surface or volume) */
    printf("Read files:");
    fflush(stdout);
    for (i = 0; i < n_subj; i++) {
        if (is_volume) {
            /* header only: read_nifti_double still fills tmpvals with the */
            /* voxel values and frees the raw nii_ptr->data buffer */
            nifti_image *cur_ptr = read_nifti_double(scan_files[i],
                                                     &tmpvals, 0);
            if (cur_ptr == NULL) {
                fprintf(stderr, "\nError reading file %s\n", scan_files[i]);
                return 0;
            }
            n_vals = (int) cur_ptr->nvox;
            if (nii_ptr == NULL) {
                nii_ptr = cur_ptr;
                dim[0] = nii_ptr->nx;
                dim[1] = nii_ptr->ny;
                dim[2] = nii_ptr->nz;
                vox[0] = nii_ptr->dx;
                vox[1] = nii_ptr->dy;
                vox[2] = nii_ptr->dz;
            } else {
                cur_ptr->data = NULL;
                nifti_image_free(cur_ptr);
            }
        } else {
            if (input_values_any_format(scan_files[i], &n_vals,
                                        &tmpvals) != OK) {
                fprintf(stderr, "\nError reading file %s\n", scan_files[i]);
                return 0;
            }
        }

        if (i == 0) {
            ALLOC2D(vals, n_subj, n_vals);
        } else if (prev_n_vals != n_vals) {
            fprintf(stderr, "\nError: Wrong number of values in %s\n",
                    scan_files[i]);
            return 0;
        }
        for (k = 0; k < n_vals; k++)
            vals[i][k] = tmpvals[k];
        prev_n_vals = n_vals;

        /* surface data uses bicpl ALLOC, volume data uses malloc */
        if (is_volume)
            free(tmpvals);
        else
            FREE(tmpvals);

        printf(".");
        fflush(stdout);
    }
    printf("\n");

    ALLOC(v, n_vals);
    ALLOC(beta0, n_vals);
    ALLOC2D(estimates, n_subj, n_vals);
    ALLOC2D(beta, n_beta, n_vals);

    /* beta = pinv(G) * Y */
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

    /* write betas for each column of the design matrix */
    for (j = 0; j < n_beta; j++) {
        outfile = create_string("beta");
        sprintf(buffer, "_%04d.%s", j+1, ext);
        concat_to_string(&outfile, buffer);
        for (k = 0; k < n_vals; k++) beta0[k] = beta[j][k];
        write_result_values(outfile, is_volume, n_vals, beta0, dim, vox,
                             nii_ptr);
    }

    /* fitted data and residual mean square v = sum((Y - G*beta).^2)/erdf */
    matrix_multiply(n_subj, n_beta, n_vals, G, beta, estimates);
    for (k = 0; k < n_vals; k++) {
        sum = 0.0;
        for (i = 0; i < n_subj; i++) {
            estimates[i][k] = vals[i][k] - estimates[i][k];
            sum += estimates[i][k] * estimates[i][k];
        }
        v[k] = sum / erdf;
    }

    /* residual mean square (ResMS), as used by SPM12 */
    sprintf(buffer, "ResMS.%s", ext);
    outfile = create_string(buffer);
    write_result_values(outfile, is_volume, n_vals, v, dim, vox, nii_ptr);

    FREE(v);
    FREE(beta0);
    FREE2D(estimates);
    FREE2D(beta);
    FREE2D(vals);
    FREE2D(inv_G);

    if (nii_ptr != NULL) {
        /* data was freed by read_nifti_double / write_nifti_double */
        nii_ptr->data = NULL;
        nifti_image_free(nii_ptr);
    }

    return 1;
}


/* ----------------------------------------------------------------------- */
/*  legacy interface: groups via '+' and covariates via ':'                */
/* ----------------------------------------------------------------------- */
int
estimate(char **infiles, int argc)
{
    char         **scan_files;
    int          i, j, k, n_subj, n_tmp, n_beta, n_cov, counter, *idx, ret;
    double       **G, *tmpvals;

    counter = 0;
    n_beta  = 1;
    n_cov   = 0;

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
            if (!file_exists(infiles[i])) {
                fprintf(stderr, "\nFile %s not found.\n", infiles[i]);
                exit(EXIT_FAILURE);
            }
            counter++;
        }
    }

    n_subj = argc - n_cov - n_beta;

    ALLOC2D(G, n_subj, n_beta);
    for (j = 0; j < n_beta; j++)
        for (i = 0; i < n_subj; i++)
            G[i][j] = 0;

    /* build the design matrix from the group/covariate layout */
    counter = 0;
    for (i = 0, j = 0; i < argc-1; i++) {
        if (equal_strings(infiles[i], ":")) {
            i++; j++;
            if (input_values_any_format(infiles[i], &n_tmp, &tmpvals) != OK) {
                fprintf(stderr, "\nError reading file %s\n", infiles[i]);
                exit(EXIT_FAILURE);
            }
            if (n_subj != n_tmp) {
                fprintf(stderr, "\nError: Number of files differs from number of rows in design matrix.\n");
                exit(EXIT_FAILURE);
            }
            /* normalize covariates to zero mean */
            normalize_double(tmpvals, n_tmp);
            for (k = 0; k < n_subj; k++)
                G[k][j] = tmpvals[k];
            FREE(tmpvals);
        } else if (equal_strings(infiles[i], "+")) {
            j++;
        } else {
            G[counter][j] = 1;
            counter++;
        }
    }

    /* collect the ordered scan file list and estimate */
    ALLOC(scan_files, n_subj);
    for (i = 0; i < n_subj; i++)
        scan_files[i] = infiles[idx[i]];

    ret = estimate_core(scan_files, n_subj, G, n_beta, NULL);

    FREE(scan_files);
    FREE(idx);
    FREE2D(G);

    return ret ? 0 : 1;
}


/* ----------------------------------------------------------------------- */
/*  formula interface: R-style model formula                               */
/* ----------------------------------------------------------------------- */
static int
estimate_formula(char **scan_files, int n_subj, const char *formula)
{
    GlmDesign design;
    int       i, ret;

    for (i = 0; i < n_subj; i++) {
        if (!file_exists(scan_files[i])) {
            fprintf(stderr, "\nFile %s not found.\n", scan_files[i]);
            return 1;
        }
    }

    if (!glm_build_design(formula, n_subj, &design)) {
        fprintf(stderr, "Error building design matrix from formula.\n");
        return 1;
    }
    printf("Design: %d observations x %d columns\n",
           design.n_obs, design.n_beta);

    ret = estimate_core(scan_files, n_subj, design.G, design.n_beta,
                        design.colnames);

    glm_free_design(&design);
    return ret ? 0 : 1;
}


int
main(int argc, char *argv[])
{
    char **infiles, **scan_files, *formula = NULL;
    int  i, a, n_scan, rc;

    initialize_argument_processing(argc, argv);

    if (argc < 2) {
        usage(argv[0]);
        exit(EXIT_FAILURE);
    }

    /* look for the -formula option */
    for (i = 1; i < argc; i++) {
        if (equal_strings(argv[i], "-formula")) {
            if (i+1 >= argc) {
                usage(argv[0]);
                exit(EXIT_FAILURE);
            }
            formula = argv[i+1];
            break;
        }
    }

    if (formula != NULL) {
        /* remaining positional arguments are the scan files */
        ALLOC(scan_files, argc);
        n_scan = 0;
        for (a = 1; a < argc; a++) {
            if (equal_strings(argv[a], "-formula")) {
                a++;                 /* skip the formula string */
                continue;
            }
            scan_files[n_scan++] = argv[a];
        }
        if (n_scan == 0) {
            fprintf(stderr, "No scan files specified.\n");
            exit(EXIT_FAILURE);
        }
        rc = estimate_formula(scan_files, n_scan, formula);
        FREE(scan_files);
        return rc ? EXIT_FAILURE : EXIT_SUCCESS;
    }

    infiles = &argv[1];
    estimate(infiles, argc);

    return EXIT_SUCCESS;
}
