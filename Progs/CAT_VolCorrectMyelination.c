/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 */

/**
 * \file CAT_VolCorrectMyelination.c
 * \brief CLI tool to correct PVE labels for myelination bias using the raw T1w.
 *
 * Usage:
 *   CAT_VolCorrectMyelination [options] pve_label.nii t1w.nii output.nii
 *
 * This tool reads a PVE label image and the corresponding raw T1w image,
 * applies the myelination correction from CAT_VolMyelinCorrection, and
 * writes the corrected PVE labels.
 *
 * Intended to be called before PBT thickness estimation.
 */

#include <stdlib.h>
#include <stdio.h>

#include <bicpl.h>
#include <ParseArgv.h>

#include "CAT_NiftiLib.h"
#include "CAT_VolMyelinCorrection.h"

/* Default option values */
double erosion_mm = 3.0;
double k_intensity = 1.0;
double grad_percentile = 25.0;
double dist_mm = 2.0;
double max_correction = 0.5;
double min_cluster_mm3 = 5.0;
double gm_grad_pct = 50.0;
double max_gm_grad_dist = 3.0;
int n_median_filter = 1;
int no_wm_corr = 0;
int csf_corr = 0;
int verbose = 0;

static ArgvInfo argTable[] = {
    {"-erosion", ARGV_FLOAT, (char *)TRUE, (char *)&erosion_mm,
     "Erosion radius (mm) defining the deep WM core (default: 3.0)."},
    {"-k", ARGV_FLOAT, (char *)TRUE, (char *)&k_intensity,
     "Intensity threshold: flag if T1w < WM_mean - k * WM_std (default: 1.0)."},
    {"-grad-pct", ARGV_FLOAT, (char *)TRUE, (char *)&grad_percentile,
     "Percentile of boundary gradient below which gradient is 'weak' (default: 25.0)."},
    {"-dist", ARGV_FLOAT, (char *)TRUE, (char *)&dist_mm,
     "Minimum distance (mm) from deep WM core to flag (default: 2.0)."},
    {"-max-corr", ARGV_FLOAT, (char *)TRUE, (char *)&max_correction,
     "Maximum PVE shift toward GM (default: 0.5)."},
    {"-min-cluster", ARGV_FLOAT, (char *)TRUE, (char *)&min_cluster_mm3,
     "Minimum cluster volume in mm^3 to keep (default: 5.0)."},
    {"-gm-grad-pct", ARGV_FLOAT, (char *)TRUE, (char *)&gm_grad_pct,
     "Percentile of GM gradient for double-gradient detection (default: 50.0)."},
    {"-max-gm-dist", ARGV_FLOAT, (char *)TRUE, (char *)&max_gm_grad_dist,
     "Max distance (mm) from high-gradient GM for myelination detection (default: 3.0). "
     "Set to 0 to disable."},
    {"-median", ARGV_INT, (char *)TRUE, (char *)&n_median_filter,
     "Iterations of 3x3x3 median filter on correction field (default: 1)."},
    {"-no-wm-corr", ARGV_CONSTANT, (char *)TRUE, (char *)&no_wm_corr,
     "Disable WM/GM (myelination) boundary correction."},
    {"-csf-corr", ARGV_CONSTANT, (char *)TRUE, (char *)&csf_corr,
     "Disable GM/CSF (pial) boundary correction."},
    {"-verbose", ARGV_CONSTANT, (char *)TRUE, (char *)&verbose,
     "Enable verbose output."},
    {NULL, ARGV_END, NULL, NULL, NULL}};

static void
usage(const char *executable)
{
    fprintf(stderr,
            "\nUsage: %s [options] pve_label.nii t1w.nii output.nii\n\n"
            "Correct PVE tissue labels for increased cortical myelination\n"
            "using the raw T1w intensity image.\n\n"
            "This tool identifies voxels near the GM/WM and GM/CSF boundaries\n"
            "that were likely misclassified due to myelination or PVE effects,\n"
            "and shifts their PVE values toward GM.\n\n"
            "Both corrections are enabled by default. Use -no-wm-corr or\n"
            "-csf-corr to disable/enable either side individually.\n\n"
            "The corrected PVE map can then be used for thickness estimation\n"
            "(e.g., CAT_VolThicknessPbt) with reduced boundary bias.\n\n",
            executable);
}

int main(int argc, char *argv[])
{
    char *pve_file, *t1w_file, *out_file;
    float *pve_data, *t1w_data;
    nifti_image *pve_nii, *t1w_nii;
    int dims[3];
    double voxelsize[3];

    /* Parse optional flags */
    if (ParseArgv(&argc, argv, argTable, 0))
    {
        usage(argv[0]);
        exit(EXIT_FAILURE);
    }

    initialize_argument_processing(argc, argv);

    if (!get_string_argument(NULL, &pve_file) ||
        !get_string_argument(NULL, &t1w_file) ||
        !get_string_argument(NULL, &out_file))
    {
        usage(argv[0]);
        exit(EXIT_FAILURE);
    }

    /* Load PVE label volume */
    pve_nii = read_nifti_float(pve_file, &pve_data, 0);
    if (!pve_nii)
    {
        fprintf(stderr, "Error reading PVE label volume: %s\n", pve_file);
        return EXIT_FAILURE;
    }

    /* Load T1w intensity volume */
    t1w_nii = read_nifti_float(t1w_file, &t1w_data, 0);
    if (!t1w_nii)
    {
        fprintf(stderr, "Error reading T1w volume: %s\n", t1w_file);
        return EXIT_FAILURE;
    }

    /* Verify dimensions match */
    if (pve_nii->nx != t1w_nii->nx ||
        pve_nii->ny != t1w_nii->ny ||
        pve_nii->nz != t1w_nii->nz)
    {
        fprintf(stderr, "Error: PVE and T1w volumes have different dimensions.\n");
        return EXIT_FAILURE;
    }

    dims[0] = pve_nii->nx;
    dims[1] = pve_nii->ny;
    dims[2] = pve_nii->nz;
    voxelsize[0] = pve_nii->dx;
    voxelsize[1] = pve_nii->dy;
    voxelsize[2] = pve_nii->dz;

    /* Set up options */
    CAT_MyelinCorrOptions opts;
    CAT_MyelinCorrOptionsInit(&opts);
    opts.erosion_mm = erosion_mm;
    opts.k_intensity = k_intensity;
    opts.grad_percentile = grad_percentile;
    opts.dist_mm = dist_mm;
    opts.max_correction = max_correction;
    opts.n_median_filter = n_median_filter;
    opts.min_cluster_mm3 = min_cluster_mm3;
    opts.gm_grad_pct = gm_grad_pct;
    opts.max_gm_grad_dist = max_gm_grad_dist;
    opts.correct_wm = !no_wm_corr;
    opts.correct_csf = csf_corr;
    opts.verbose = verbose;

    /* Apply correction */
    if (CAT_VolCorrectMyelination(pve_data, t1w_data, dims, voxelsize, &opts) != 0)
    {
        fprintf(stderr, "Error during myelination correction.\n");
        return EXIT_FAILURE;
    }

    /* Write corrected PVE */
    nifti_image *out_nii = nifti_copy_nim_info(pve_nii);
    if (!write_nifti_float(out_file, pve_data, DT_FLOAT32, 1.0,
                           dims, voxelsize, out_nii))
    {
        fprintf(stderr, "Error writing output: %s\n", out_file);
        return EXIT_FAILURE;
    }

    if (verbose)
        fprintf(stderr, "Corrected PVE written to %s\n", out_file);

    free(pve_data);
    free(t1w_data);

    return EXIT_SUCCESS;
}
