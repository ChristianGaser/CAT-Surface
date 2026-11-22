/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 */

/* Unified segmentation with geodesic shooting: the six classes of a shooting
 * template series (GM, WM, CSF, bone, soft tissue, air) from one image.  This
 * is a thin CLI wrapper around CAT_UnifiedSeg().
 *
 * Input  : the image and the five shooting templates.
 * Output : whatever is asked for by an output option: class maps, a label map
 *          with partial volumes, the bias-corrected image, the deformation
 *          and the parameters as JSON.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "ParseArgv.h"
#include "CAT_NiftiLib.h"
#include "CAT_ArgvCheck.h"
#include "CAT_UnifiedSeg.h"

static CAT_UnifiedSegOpts opts; /* defaults set in main() before parsing */

static char *prob_file = NULL;
static char *label_file = NULL;
static char *pve_file = NULL;
static char *corrected_file = NULL;
static char *bias_file = NULL;
static char *y_file = NULL;
static char *iy_file = NULL;
static char *jac_file = NULL;
static char *params_file = NULL;
static char *affine_file = NULL;
static char *mask_file = NULL;
static char *affreg_name = NULL;
static char *bias_model_name = NULL;
static int affine_fixed = 0;
static int affine_refine = 0;
static int no_pve = 0;
static int verbose2 = 0;

static ArgvInfo argTable[] = {
    {"-prob", ARGV_STRING, (char *)1, (char *)&prob_file,
     "Write the six class maps (GM, WM, CSF, bone, soft tissue, air) as one 4-D file."},
    {"-label", ARGV_STRING, (char *)1, (char *)&label_file,
     "Write the label map: 0 background, 1 CSF, 2 GM, 3 WM, partial volumes in between."},
    {"-pve", ARGV_STRING, (char *)1, (char *)&pve_file,
     "Write the fractions of the six classes as one 4-D file: CSF, GM and WM from\n\t\t\tthe partial volume model, bone, soft tissue and air their posteriors\n\t\t\t(with -pve-head from a partial volume model of the head)."},
    {"-corrected", ARGV_STRING, (char *)1, (char *)&corrected_file,
     "Write the bias-corrected image."},
    {"-bias", ARGV_STRING, (char *)1, (char *)&bias_file,
     "Write the bias field (corrected = image * field)."},
    {"-y", ARGV_STRING, (char *)1, (char *)&y_file,
     "Write the deformation in SPM's y_ format: on the template grid, the image\n\t\t\tposition in mm of every template voxel."},
    {"-iy", ARGV_STRING, (char *)1, (char *)&iy_file,
     "Write the inverse deformation: on the image grid, the template position in mm."},
    {"-jac", ARGV_STRING, (char *)1, (char *)&jac_file,
     "Write the Jacobian determinant of the deformation (template grid)."},
    {"-params", ARGV_STRING, (char *)1, (char *)&params_file,
     "Write the estimated parameters and a summary of the run as JSON."},
    {"-sanlm", ARGV_CONSTANT, (char *)1, (char *)&opts.sanlm,
     "Segment the image denoised with the SANLM filter; the parameters are still\n\t\t\testimated from the image as it came."},
    {"-sanlm-strength", ARGV_FLOAT, (char *)1, (char *)&opts.sanlm_strength,
     "Strength of the SANLM filter."},
    {"-rician", ARGV_CONSTANT, (char *)1, (char *)&opts.sanlm_rician,
     "Rician noise model of the SANLM filter."},
    {"-noise", ARGV_FLOAT, (char *)1, (char *)&opts.noise,
     "Add Gaussian noise to the image the parameters are estimated from, until its\n\t\t\tnoise SD reaches this fraction of the intensity range (e.g. 0.02).  The\n\t\t\tsegmentation itself is computed without it."},
    {"-noise-seed", ARGV_INT, (char *)1, (char *)&opts.noise_seed,
     "Seed of the added noise."},
    {"-ngaus", ARGV_INT, (char *)CAT_USEG_NCLASSES, (char *)opts.ngaus,
     "Gaussians per class (six values)."},
    {"-samp", ARGV_FLOAT, (char *)1, (char *)&opts.samp,
     "Distance in mm between the voxels the parameters are estimated from."},
    {"-bias-model", ARGV_STRING, (char *)1, (char *)&bias_model_name,
     "Basis functions of the bias field: spline (cubic B-splines) or dct (as SPM)."},
    {"-biasfwhm", ARGV_FLOAT, (char *)1, (char *)&opts.biasfwhm,
     "Smoothness of the bias field as a FWHM in mm; 0: no bias correction."},
    {"-biasreg", ARGV_FLOAT, (char *)1, (char *)&opts.biasreg,
     "Regularisation of the bias field."},
    {"-wp-reg", ARGV_FLOAT, (char *)1, (char *)&opts.wp_reg,
     "Holds the class weights near 1/6, as a multiple of the prior mass of a class;\n\t\t\t0: free as in SPM, large: fixed."},
    {"-tol", ARGV_FLOAT, (char *)1, (char *)&opts.tol,
     "Stopping criterion of the EM steps (log-likelihood change per voxel)."},
    {"-affreg", ARGV_STRING, (char *)1, (char *)&affreg_name,
     "Regularisation of the affine registration: mni, eastern, subj, rigid or none."},
    {"-affine", ARGV_STRING, (char *)1, (char *)&affine_file,
     "Text file with a 4x4 image-to-template world matrix: one more starting point\n\t\t\tof the affine registration."},
    {"-affine-refine", ARGV_CONSTANT, (char *)1, (char *)&affine_refine,
     "Refine the matrix of -affine only, without the other starting points."},
    {"-affine-fixed", ARGV_CONSTANT, (char *)1, (char *)&affine_fixed,
     "Use the matrix of -affine as it is."},
    {"-mask", ARGV_STRING, (char *)1, (char *)&mask_file,
     "Mask on the image grid: voxels where it is zero carry no data (lesions)."},
    {"-mrf", ARGV_FLOAT, (char *)1, (char *)&opts.mrf,
     "Strength of the Markov random field cleanup of the class maps (SPM: 1)."},
    {"-cleanup", ARGV_INT, (char *)1, (char *)&opts.cleanup,
     "SPM's cleanup of the brain classes: 0 off, 1 light, 2 thorough."},
    {"-no-pve", ARGV_CONSTANT, (char *)1, (char *)&no_pve,
     "No partial volume model: the label is the posterior mean."},
    {"-pve-head", ARGV_CONSTANT, (char *)1, (char *)&opts.pve_head,
     "Partial volume model of the head as well: the fractions of bone, soft tissue\n\t\t\tand air from mixtures of soft tissue with bone and with air, instead of\n\t\t\ttheir posteriors."},
    {"-regstr", ARGV_FLOAT, (char *)1, (char *)&opts.shoot.regstr,
     "Shooting: 1 runs every scheduled iteration, (0,1) leaves a level early."},
    {"-reg-scale", ARGV_FLOAT, (char *)1, (char *)&opts.shoot.reg_scale,
     "Shooting: factor on the regularisation; smaller values deform more."},
    {"-lowres", ARGV_FLOAT, (char *)1, (char *)&opts.shoot.lowres,
     "Shooting: coarsest registration resolution in mm."},
    {"-nits", ARGV_INT, (char *)1, (char *)&opts.shoot.nits,
     "Shooting: maximum number of Gauss-Newton iterations."},
    {"-threads", ARGV_INT, (char *)1, (char *)&opts.n_threads,
     "Number of worker threads (default: the library's)."},
    {"-verbose", ARGV_CONSTANT, (char *)1, (char *)&opts.verbose,
     "Print the stages and the estimated parameters."},
    {"-debug", ARGV_CONSTANT, (char *)1, (char *)&verbose2,
     "Print every iteration."},
    {NULL, ARGV_END, NULL, NULL, NULL}};

static void
usage(const char *executable)
{
    CAT_UnifiedSegOpts d;

    CAT_UnifiedSegOptionsInit(&d);
    fprintf(stderr,
            "\nUsage: %s [options] <image> <Template_0> [<Template_1> .. <Template_4>]\n\n"
            "Unified segmentation with geodesic shooting.  The image is classified into\n"
            "the six classes of a shooting template series -- GM, WM, CSF, bone, soft\n"
            "tissue, air -- with the model of SPM's unified segmentation: a mixture of\n"
            "Gaussians per class, a smooth bias field, and the templates as priors, which\n"
            "are deformed onto the image by geodesic shooting on a grid of increasing\n"
            "resolution:\n\n"
            "    %s -prob c_sub.nii -label p0_sub.nii -params sub_seg8.json \\\n"
            "        sub.nii Template_0_GS.nii\n\n"
            "A single template whose name contains _0_ names the series: the other four\n"
            "replace that with _1_ .. _4_.  Every template holds the six classes.\n"
            "Voxels that are zero or not finite carry no data and get their prior, so\n"
            "skull-stripped and defaced images can be segmented.\n\n"
            "Outputs (at least one):\n"
            "  -prob <file>       the six class maps, 4-D\n"
            "  -label <file>      label map with partial volumes (0, CSF 1, GM 2, WM 3)\n"
            "  -pve <file>        fractions of the six classes, 4-D (bone, soft tissue\n"
            "                     and air: their posteriors unless -pve-head)\n"
            "  -corrected <file>  bias-corrected image\n"
            "  -bias <file>       bias field\n"
            "  -y <file>          deformation, SPM's y_ format (template grid)\n"
            "  -iy <file>         inverse deformation (image grid)\n"
            "  -jac <file>        Jacobian determinants (template grid)\n"
            "  -params <file>     estimated parameters and summary, JSON\n\n"
            "Options (defaults in brackets):\n"
            "  -sanlm             segment the SANLM-filtered image (parameters are\n"
            "                     estimated from the image as it came)\n"
            "  -sanlm-strength <f>  strength of the filter [%g]\n"
            "  -rician            Rician noise model of the filter\n"
            "  -noise <f>         estimate the parameters from the image with noise\n"
            "                     added up to this fraction of its intensity range;\n"
            "                     the segmentation is computed without it [%g]\n"
            "  -noise-seed <n>    seed of the added noise [%u]\n"
            "  -ngaus <6 values>  Gaussians per class [%d %d %d %d %d %d]\n"
            "  -samp <mm>         distance of the voxels used for estimation [%g]\n"
            "  -bias-model <name> spline or dct (SPM's basis) [%s]\n"
            "  -biasfwhm <mm>     smoothness of the bias field; 0: none [%g]\n"
            "  -biasreg <f>       regularisation of the bias field [%g]\n"
            "  -wp-reg <f>        holds the class weights near 1/6; 0: free as SPM [%g]\n"
            "  -tol <f>           stopping criterion of the EM steps [%g]\n"
            "  -affreg <name>     mni, eastern, subj, rigid or none [mni]\n"
            "  -affine <file>     4x4 image-to-template world matrix, a starting point\n"
            "  -affine-refine     refine that matrix only\n"
            "  -affine-fixed      use that matrix as it is\n"
            "  -mask <file>       zero marks voxels without data\n"
            "  -mrf <f>           Markov random field cleanup of the class maps [%g]\n"
            "  -cleanup <n>       SPM's cleanup of the brain classes: 0 off, 1 light,\n"
            "                     2 thorough [%d]\n"
            "  -no-pve            label map without the partial volume model\n"
            "  -pve-head          partial volume model of bone, soft tissue and air too\n"
            "  -regstr <v>        shooting: 1 all iterations, (0,1) early level switch [%g]\n"
            "  -reg-scale <f>     shooting: factor on the regularisation [%g]\n"
            "  -lowres <mm>       shooting: coarsest resolution [%g]\n"
            "  -nits <n>          shooting: maximum Gauss-Newton iterations [%d]\n"
            "  -threads <n>       worker threads\n"
            "  -verbose           print stages and parameters\n"
            "  -debug             print every iteration\n\n",
            executable, executable, d.sanlm_strength, d.noise, d.noise_seed, d.ngaus[0],
            d.ngaus[1], d.ngaus[2], d.ngaus[3], d.ngaus[4], d.ngaus[5], d.samp,
            (d.bias_model == CAT_USEG_BIAS_DCT) ? "dct" : "spline", d.biasfwhm,
            d.biasreg, d.wp_reg, d.tol, d.mrf, d.cleanup, d.shoot.regstr, d.shoot.reg_scale,
            d.shoot.lowres, d.shoot.nits);
}

/* sto_xyz of a header as a row-major double matrix */
static void
header_matrix(const nifti_image *nii, double M[16])
{
    int r, c;
    for (r = 0; r < 4; r++)
        for (c = 0; c < 4; c++) M[4 * r + c] = (double)nii->sto_xyz.m[r][c];
}

static int
read_matrix(const char *file, double M[16])
{
    FILE *fp = fopen(file, "r");
    int i;
    if (!fp) return 1;
    for (i = 0; i < 16; i++)
        if (fscanf(fp, "%lf", &M[i]) != 1)
        {
            fclose(fp);
            return 1;
        }
    fclose(fp);
    return 0;
}

static float *
read_map(const char *file, int dims[3], int *nvol, nifti_image **hdr)
{
    float *data = NULL;
    nifti_image *nii = read_nifti_float(file, &data, 0);
    if (!nii)
    {
        fprintf(stderr, "Error reading %s.\n", file);
        exit(EXIT_FAILURE);
    }
    dims[0] = nii->nx;
    dims[1] = nii->ny;
    dims[2] = nii->nz;
    if (nvol)
        *nvol = (nii->nt > 1 ? nii->nt : 1) * (nii->nu > 1 ? nii->nu : 1) *
                (nii->nv > 1 ? nii->nv : 1) * (nii->nw > 1 ? nii->nw : 1);
    if (hdr) *hdr = nii;
    else nifti_image_free(nii);
    return data;
}

/* Template_0 -> Template_k by replacing the last "_0_" */
static char *
template_name(const char *first, int k)
{
    const char *p = NULL, *q = first;
    char *name;
    while ((q = strstr(q, "_0_")) != NULL) p = q++;
    if (!p) return NULL;
    name = (char *)malloc(strlen(first) + 1);
    if (!name) return NULL;
    strcpy(name, first);
    name[p - first + 1] = (char)('0' + k);
    return name;
}

/* Write nvol volumes of bytes with the scaling 1/255 on the grid of hdr. */
static int
write_bytes(const char *file, unsigned char *data, const int dims[3], int nvol,
            const nifti_image *hdr)
{
    nifti_image *nim = nifti_copy_nim_info(hdr);

    if (!nim) return 1;
    nim->ndim = nim->dim[0] = (nvol > 1) ? 4 : 3;
    nim->nx = nim->dim[1] = dims[0];
    nim->ny = nim->dim[2] = dims[1];
    nim->nz = nim->dim[3] = dims[2];
    nim->nt = nim->dim[4] = nvol;
    nim->nu = nim->dim[5] = 1;
    nim->nv = nim->dim[6] = 1;
    nim->nw = nim->dim[7] = 1;
    nim->nvox = (size_t)dims[0] * dims[1] * dims[2] * nvol;
    nim->datatype = DT_UINT8;
    nim->nbyper = 1;
    nim->scl_slope = 1.0f / 255.0f;
    nim->scl_inter = 0.0f;
    nim->cal_min = 0.0f;
    nim->cal_max = 1.0f;
    if (nifti_set_filenames(nim, file, 0, 1))
    {
        nifti_image_free(nim);
        return 1;
    }
    nim->data = data;
    nifti_image_write(nim);
    nim->data = NULL; /* the caller's */
    nifti_image_free(nim);
    return 0;
}

/* Write a deformation in SPM's 5-D layout (x, y, z, 1, 3) on the grid of
 * grid: the field holds voxel coordinates (0-based) of the other grid, which
 * vmat turns into mm.  Undefined voxels (NaN) stay NaN. */
static int
write_deformation(const char *file, const float *vox, const int dims[3],
                  const double vmat[16], const nifti_image *grid)
{
    nifti_image *nim = nifti_copy_nim_info(grid);
    size_t n = (size_t)dims[0] * dims[1] * dims[2], i;
    float *mm;
    int j;

    if (!nim) return 1;
    mm = (float *)malloc(sizeof(float) * 3 * n);
    if (!mm)
    {
        nifti_image_free(nim);
        return 1;
    }
    for (i = 0; i < n; i++)
    {
        double a = vox[i], b = vox[n + i], c = vox[2 * n + i];
        for (j = 0; j < 3; j++)
            mm[j * n + i] = (float)(vmat[4 * j] * a + vmat[4 * j + 1] * b +
                                    vmat[4 * j + 2] * c + vmat[4 * j + 3]);
    }
    nim->ndim = nim->dim[0] = 5;
    nim->nx = nim->dim[1] = dims[0];
    nim->ny = nim->dim[2] = dims[1];
    nim->nz = nim->dim[3] = dims[2];
    nim->nt = nim->dim[4] = 1;
    nim->nu = nim->dim[5] = 3;
    nim->nv = nim->dim[6] = 1;
    nim->nw = nim->dim[7] = 1;
    nim->nvox = 3 * n;
    nim->datatype = DT_FLOAT32;
    nim->nbyper = 4;
    nim->scl_slope = 1.0f;
    nim->scl_inter = 0.0f;
    nim->intent_code = NIFTI_INTENT_VECTOR;
    if (nifti_set_filenames(nim, file, 0, 1))
    {
        free(mm);
        nifti_image_free(nim);
        return 1;
    }
    nim->data = mm;
    nifti_image_write(nim);
    nifti_image_free(nim); /* frees mm as nim->data */
    return 0;
}

static int
parse_affreg(const char *name)
{
    static const char *names[] = {"none", "mni", "eastern", "subj", "rigid"};
    int i;
    for (i = 0; i < 5; i++)
        if (!strcmp(name, names[i])) return i;
    return -1;
}

int
main(int argc, char *argv[])
{
    const float *tpl[CAT_SHOOT_NTEMPLATES];
    float *src, *label = NULL, *corrected = NULL, *bias = NULL, *y = NULL, *iy = NULL, *jac = NULL;
    unsigned char *prob = NULL, *pve = NULL, *mask = NULL;
    nifti_image *hdr_in = NULL, *hdr_tpl = NULL;
    double mat[16], tmat[16], Aff[16], vox[3], tvox[3];
    int dims[3], d2[3], tdims[3], k, nvol, ret;
    CAT_UnifiedSegResult res;
    size_t n, nt, i;

    CAT_UnifiedSegOptionsInit(&opts);
    if (ParseArgv(&argc, argv, argTable, 0) ||
        cat_check_unknown_options(argc, argv, argv[0]))
    {
        usage(argv[0]);
        exit(EXIT_FAILURE);
    }
    if (!(argc == 3 || argc == 2 + CAT_SHOOT_NTEMPLATES) ||
        (!prob_file && !label_file && !pve_file && !corrected_file && !bias_file &&
         !y_file && !iy_file && !jac_file && !params_file))
    {
        usage(argv[0]);
        exit(EXIT_FAILURE);
    }
    if (verbose2) opts.verbose = 2;
    if (no_pve) opts.pve = 0;
    if (opts.cleanup < 0 || opts.cleanup > 2)
    {
        fprintf(stderr, "Error: -cleanup must be 0, 1 or 2.\n");
        exit(EXIT_FAILURE);
    }
    if (affreg_name && (opts.affreg = parse_affreg(affreg_name)) < 0)
    {
        fprintf(stderr, "Error: unknown -affreg %s.\n", affreg_name);
        exit(EXIT_FAILURE);
    }
    if (bias_model_name)
    {
        if (!strcmp(bias_model_name, "dct")) opts.bias_model = CAT_USEG_BIAS_DCT;
        else if (!strcmp(bias_model_name, "spline")) opts.bias_model = CAT_USEG_BIAS_SPLINE;
        else
        {
            fprintf(stderr, "Error: unknown -bias-model %s.\n", bias_model_name);
            exit(EXIT_FAILURE);
        }
    }
    if ((affine_fixed || affine_refine) && !affine_file)
    {
        fprintf(stderr, "Error: -affine-fixed and -affine-refine need -affine.\n");
        exit(EXIT_FAILURE);
    }
    opts.affine_mode = affine_fixed ? CAT_USEG_AFFINE_FIXED
                                    : (affine_refine ? CAT_USEG_AFFINE_REFINE : CAT_USEG_AFFINE_ESTIMATE);

    /* the image */
    src = read_map(argv[1], dims, NULL, &hdr_in);
    header_matrix(hdr_in, mat);
    vox[0] = hdr_in->dx;
    vox[1] = hdr_in->dy;
    vox[2] = hdr_in->dz;
    n = (size_t)dims[0] * dims[1] * dims[2];

    /* the templates */
    for (k = 0; k < CAT_SHOOT_NTEMPLATES; k++)
    {
        nifti_image *h;
        char *tn = (argc == 3) ? template_name(argv[2], k) : NULL;
        const char *name = (argc == 3) ? tn : argv[2 + k];
        if (!name)
        {
            fprintf(stderr, "Error: %s does not name a series (no _0_); give all %d templates.\n",
                    argv[2], CAT_SHOOT_NTEMPLATES);
            exit(EXIT_FAILURE);
        }
        tpl[k] = read_map(name, d2, &nvol, &h);
        if (nvol != CAT_USEG_NCLASSES)
        {
            fprintf(stderr, "Error: %s holds %d volumes, not the %d classes.\n", name, nvol,
                    CAT_USEG_NCLASSES);
            exit(EXIT_FAILURE);
        }
        if (k == 0)
        {
            memcpy(tdims, d2, sizeof(tdims));
            header_matrix(h, tmat);
            hdr_tpl = h;
        }
        else
        {
            if (d2[0] != tdims[0] || d2[1] != tdims[1] || d2[2] != tdims[2])
            {
                fprintf(stderr, "Error: the templates differ in size.\n");
                exit(EXIT_FAILURE);
            }
            nifti_image_free(h);
        }
        if (opts.verbose) printf("level %d: template %s\n", k, name);
        free(tn);
    }
    nt = (size_t)tdims[0] * tdims[1] * tdims[2];
    tvox[0] = hdr_tpl->dx;
    tvox[1] = hdr_tpl->dy;
    tvox[2] = hdr_tpl->dz;

    if (affine_file && read_matrix(affine_file, Aff))
    {
        fprintf(stderr, "Error reading the 4x4 matrix in %s.\n", affine_file);
        exit(EXIT_FAILURE);
    }
    if (mask_file)
    {
        float *m = read_map(mask_file, d2, NULL, NULL);
        if (d2[0] != dims[0] || d2[1] != dims[1] || d2[2] != dims[2])
        {
            fprintf(stderr, "Error: %s is not on the image grid.\n", mask_file);
            exit(EXIT_FAILURE);
        }
        mask = (unsigned char *)malloc(n);
        if (!mask)
        {
            fprintf(stderr, "Error: out of memory.\n");
            exit(EXIT_FAILURE);
        }
        for (i = 0; i < n; i++) mask[i] = (m[i] != 0.0f && isfinite(m[i]));
        free(m);
    }

    if (prob_file) prob = (unsigned char *)malloc(CAT_USEG_NCLASSES * n);
    if (label_file) label = (float *)malloc(sizeof(float) * n);
    if (pve_file) pve = (unsigned char *)malloc(CAT_USEG_NCLASSES * n);
    if (corrected_file) corrected = (float *)malloc(sizeof(float) * n);
    if (bias_file) bias = (float *)malloc(sizeof(float) * n);
    if (y_file) y = (float *)malloc(sizeof(float) * 3 * nt);
    if (iy_file) iy = (float *)malloc(sizeof(float) * 3 * n);
    if (jac_file) jac = (float *)malloc(sizeof(float) * nt);
    if ((prob_file && !prob) || (label_file && !label) || (pve_file && !pve) ||
        (corrected_file && !corrected) || (bias_file && !bias) || (y_file && !y) ||
        (iy_file && !iy) || (jac_file && !jac))
    {
        fprintf(stderr, "Error: out of memory.\n");
        exit(EXIT_FAILURE);
    }

    ret = CAT_UnifiedSeg(src, dims, mat, tpl, tdims, tmat, mask, affine_file ? Aff : NULL,
                         &opts, prob, label, pve, corrected, bias, y, iy, jac, &res);
    if (ret)
    {
        fprintf(stderr, "Error: the segmentation failed (%s).\n",
                (ret == 2) ? "out of memory" : ((ret == 3) ? "estimation failed" : "invalid input"));
        exit(EXIT_FAILURE);
    }

    if (prob_file && write_bytes(prob_file, prob, dims, CAT_USEG_NCLASSES, hdr_in))
        exit(EXIT_FAILURE);
    if (pve_file && write_bytes(pve_file, pve, dims, CAT_USEG_NCLASSES, hdr_in)) exit(EXIT_FAILURE);
    if (label_file && !write_nifti_float(label_file, label, DT_FLOAT32, 1.0, dims, vox, hdr_in))
        exit(EXIT_FAILURE);
    if (corrected_file &&
        !write_nifti_float(corrected_file, corrected, DT_FLOAT32, 1.0, dims, vox, hdr_in))
        exit(EXIT_FAILURE);
    if (bias_file && !write_nifti_float(bias_file, bias, DT_FLOAT32, 1.0, dims, vox, hdr_in))
        exit(EXIT_FAILURE);
    if (y_file && write_deformation(y_file, y, tdims, mat, hdr_tpl)) exit(EXIT_FAILURE);
    if (iy_file && write_deformation(iy_file, iy, dims, tmat, hdr_in)) exit(EXIT_FAILURE);
    if (jac_file && !write_nifti_float(jac_file, jac, DT_FLOAT32, 1.0, tdims, tvox, hdr_tpl))
        exit(EXIT_FAILURE);
    if (params_file && CAT_UnifiedSegWriteJson(params_file, &res, &opts, argv[1], argv[2]))
    {
        fprintf(stderr, "Error writing %s.\n", params_file);
        exit(EXIT_FAILURE);
    }
    CAT_UnifiedSegResultFree(&res);
    free(src);
    for (k = 0; k < CAT_SHOOT_NTEMPLATES; k++) free((void *)tpl[k]);
    free(mask);
    free(prob);
    free(pve);
    free(label);
    free(corrected);
    free(bias);
    free(y);
    free(iy);
    free(jac);
    nifti_image_free(hdr_in);
    nifti_image_free(hdr_tpl);
    return EXIT_SUCCESS;
}
