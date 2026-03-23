/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 *
 */

#include <float.h>
#include <stdlib.h>

#if !defined(_WIN32) && !defined(_WIN64)
#include <libgen.h>
#endif

#include <bicpl.h>
#include <ParseArgv.h>
#include "CAT_NiftiLib.h"
#include "CAT_Vol.h"
#include "CAT_VolPbt.h"
#include "CAT_Math.h"

int fast = 0;
int verbose = 0;
int n_avgs = 2;
int n_median_filter = 2;
int median_subsample = 2;
int blood_vessel_correction = 1;
double range = 0.45;
double downsample = 0.0;
double fill_thresh = 0.5;
double correct_voxelsize = 0.5;

static ArgvInfo argTable[] = {
    {"-verbose", ARGV_CONSTANT, (char *)1, (char *)&verbose,
     "Enable verbose mode. Provides detailed output during processing for debugging\n\
    and monitoring."},

    {"-fast", ARGV_CONSTANT, (char *)1, (char *)&fast,
     "Enable fast mode in order to get a very quick and rougher estimate of thickness only."},

    {"-no-blood-vessel-correction", ARGV_CONSTANT, (char *)0, (char *)&blood_vessel_correction,
     "Enable blood-vessel correction before thickness estimation (0 disables, >0 enables)."},

    {"-no-bvc", ARGV_CONSTANT, (char *)0, (char *)&blood_vessel_correction,
     "Enable blood-vessel correction before thickness estimation (0 disables, >0 enables)."},

    {"-n-avgs", ARGV_INT, (char *)1, (char *)&n_avgs,
     "Specify the number of averages for distance estimation. Used for averaging\n\
    the distances in White Matter (WM) and Cerebrospinal Fluid (CSF) to obtain a\n\
    less noisy measure. A higher number results in smoother but potentially less\n\
    accurate measures."},

    {"-fill-holes", ARGV_FLOAT, (char *)1, (char *)&fill_thresh,
     "Fill remaining holes in the PPM image using the defined threshold.\n\
    To maximize the filling-effect, this threshold should be the same as used for the\n\
    subsequent Marching Cubes approach (e.g. 0.5). Set to '0' to disable filling."},

    {"-range", ARGV_FLOAT, (char *)1, (char *)&range,
     "Extend range for masking of euclidean distance estimation. A slight increase\n\
    of range (i.e 0.3) helps in obtaining a more stable distance estimation."},

    {"-downsample", ARGV_FLOAT, (char *)1, (char *)&downsample,
     "Downsample PPM and GMT image to defined resolution since we do not need that 0.5mm\n\
    spatial resolution for the subsequent steps. Set to '0' to disable downsampling."},

    {"-median-filter", ARGV_INT, (char *)TRUE, (char *)&n_median_filter,
     "Specify the number of iterations for weighted local median filtering of the\n\
        final PPM image. The filter is not applied uniformly: CAT first estimates a\n\
        topology-artifact weight map from the positive residual PPM - smooth(PPM),\n\
        keeps only high-residual voxels in sufficiently thick cortex, regularizes this\n\
        mask by close/open/dilate and smoothing, and then blends original PPM with the\n\
        locally median-filtered PPM. Higher weights mean stronger median-filter\n\
        influence. Set to 0 to disable this cleanup."},

    {"-median-subsample", ARGV_INT, (char *)TRUE, (char *)&median_subsample,
     "Specify the size of subsampling for the median filter to smooth local\n\
     thickness values"},

    {"-correct-voxelsize", ARGV_FLOAT, (char *)1, (char *)&correct_voxelsize,
     "Amount of thickness correction for voxel-size, since we observed a systematic \n\
     shift to smaller thickness values of half voxel-size."},

    {NULL, ARGV_END, NULL, NULL, NULL}};

private void usage(char *executable)
{
    char *usage_str = "\n\
Usage: %s [options] <input.nii> output_GMT.nii output_PPM.nii [output_WMD.nii output_CSD.nii GMT1.nii GMT2.nii]\n\
\n\
    This program performs projection-based cortical thickness estimation and\n\
    percentage position mapping (PPM) from a given PVE label image described in:\n\
        Dahnke R, Yotter RA, Gaser C.\n\
        Cortical thickness and central surface estimation.\n\
        Neuroimage. 2013 Jan 15;65:336-48.\n\
    \n\
    The process involves:\n\
    \n\
    1. **Distance Estimation:**\n\
       - Estimating distances in White Matter (WM) and Cerebrospinal Fluid (CSF)\n\
         by shifting the border between Gray Matter (GM)/WM and GM/CSF.\n\
       - Averaging distances over a specified number of iterations (n-avgs) to\n\
         obtain a less noisy measure.\n\
    \n\
    2. **Thickness Estimation:**\n\
       - Reconstructing sulci and optionally gyri to estimate cortical thickness.\n\
    \n\
    3. **PPM Calculation:**\n\
       - Estimating the percentage position map (PPM), representing the relative\n\
         position within the cortical ribbon.\n\
    \n\
    4. **Final Correction and Smoothing:**\n\
       - Correcting for isotropic voxel size and applying masked smoothing\n\
         to the thickness map based on the specified Full Width Half Maximum (FWHM).\n\
\n\
Options:\n\
    -verbose                   Enable verbose mode for detailed output during processing.\n\
    -range <float>             Extend range for masking of euclidean distance estimation.\n\
    -no-bvc                    Disable blood-vessel correction.\n\
    -fast                      Enable fast mode in order to get a very quick and rougher estimate of thickness only.\n\
    -n-avgs <int>              Set the number of averages for distance estimation.\n\
    -fill-holes <float>        Define the threshold to fill holes in the PPM image.\n\
    -downsample <float>        Downsample PPM and GMT image to defined resolution.\n\
    -median-filter <int>       Iterations for weighted local PPM median filtering; higher values increase the cleanup where the topology-artifact weight map is high.\n\
    -correct-voxelsize <float> Amount of correction of thickness by voxel-size.\n\
\n\
Example:\n\
    %s -verbose -n-avgs 4 input.nii gmt_output.nii ppm_output.nii\n\n";

    fprintf(stderr, "%s\n %s\n", usage_str, executable);
}

int main(int argc, char *argv[])
{
    char out_GMT[1024], out_PPM[1024], out_CSD[1024], out_WMD[1024];
    int i, j, dims[3], dims_reduced[3];
    float *src;
    double voxelsize[3], voxelsize_reduced[3], samp[3], s[3], slope;

    initialize_argument_processing(argc, argv);

    if (ParseArgv(&argc, argv, argTable, 0) || (argc < 2))
    {
        usage(argv[0]);
        fprintf(stderr, "     %s -help\n\n", argv[0]);
        exit(EXIT_FAILURE);
    }

    char *infile = argv[1];

    /* Determine output filenames based on input filename or command-line arguments */
    if (argc >= 6)
    {
        (void)sprintf(out_WMD, "%s", argv[4]);
        (void)sprintf(out_CSD, "%s", argv[5]);
    }
    if (argc >= 4)
    {
        (void)sprintf(out_GMT, "%s", argv[2]);
        (void)sprintf(out_PPM, "%s", argv[3]);
    }
    else
    {
#if !defined(_WIN32) && !defined(_WIN64)
        (void)sprintf(out_GMT, "%s/gmt_%s", dirname(infile), basename(infile));
        (void)sprintf(out_PPM, "%s/ppm_%s", dirname(infile), basename(infile));
#else
        fprintf(stderr, "\nUsage: %s input.nii GMT.nii PPM.nii\n\n", argv[0]);
        return (1);
#endif
    }

    /* read source image */
    nifti_image *src_ptr = read_nifti_float(infile, &src, 0);
    if (!src_ptr)
    {
        fprintf(stderr, "Error reading %s.\n", infile);
        return (EXIT_FAILURE);
    }

    /* Number of voxels */
    int nvox = src_ptr->nvox;

    /* Prepare output NIfTI images */
    nifti_image *out_ptr = nifti_copy_nim_info(src_ptr);

    /* Retrieve dimensions and voxel size from source image */
    voxelsize[0] = src_ptr->dx;
    voxelsize[1] = src_ptr->dy;
    voxelsize[2] = src_ptr->dz;
    dims[0] = src_ptr->nx;
    dims[1] = src_ptr->ny;
    dims[2] = src_ptr->nz;
    float *dist_CSF = (float *)malloc(sizeof(float) * nvox);
    float *dist_WM = (float *)malloc(sizeof(float) * nvox);
    float *GMT = (float *)malloc(sizeof(float) * nvox);
    float *PPM = (float *)malloc(sizeof(float) * nvox);

    /* check for memory faults */
    if (!dist_CSF || !dist_WM || !GMT || !PPM)
    {
        fprintf(stderr, "Memory allocation error\n");
        exit(EXIT_FAILURE);
    }

    /* Change some defaults for fast option */
    if (fast)
    {
        n_avgs /= 2;
        n_median_filter = 0;
        fill_thresh = 0.0;
        downsample = 0.0;
    }

    /* Ensure that n_avgs is at least 1 */
    n_avgs = (n_avgs < 1) ? 1 : n_avgs;

    /* Optional blood-vessel correction before any other operation */
    if (blood_vessel_correction > 0.0)
    {
        if (verbose)
            fprintf(stderr, "Apply blood-vessel correction on input PVE map.\n");
        blood_vessel_correction_pve_float(src, dims, voxelsize);
    }

    CAT_PbtOptions opts;
    CAT_PbtOptionsInit(&opts);
    opts.n_avgs = n_avgs;
    opts.n_median_filter = n_median_filter;
    opts.range = range;
    opts.fill_thresh = fill_thresh;
    opts.correct_voxelsize = correct_voxelsize;
    opts.fast = fast;
    opts.verbose = verbose;
    opts.median_subsample = median_subsample;

    if (CAT_VolComputePbt(src, GMT, PPM, dist_CSF, dist_WM, dims, voxelsize, &opts) != 0)
    {
        fprintf(stderr, "Error computing projection-based thickness.\n");
        exit(EXIT_FAILURE);
    }

    /* Downsample images */
    slope = 1.0;
    if (downsample > 0.0)
    {

        nifti_image *out_ptr_reduced = nifti_copy_nim_info(src_ptr);

        for (i = 0; i < 3; i++)
        {
            s[i] = 1.2;
            voxelsize_reduced[i] = downsample;
            samp[i] = voxelsize_reduced[i] / voxelsize[i];
        }

        /*  Define grid dimensions */
        for (i = 0; i < 3; i++)
            dims_reduced[i] = (int)ceil((dims[i] - 1) / ((double)samp[i])) + 1;

        out_ptr_reduced->nvox = dims_reduced[0] * dims_reduced[1] * dims_reduced[2];

        /* Correct affine matrix */
        for (i = 0; i < 3; i++)
        {
            for (j = 0; j < 3; j++)
            {
                out_ptr_reduced->sto_xyz.m[i][j] = out_ptr->sto_xyz.m[i][j] * samp[i];
            }
        }
        out_ptr_reduced->sto_ijk = nifti_mat44_inverse(out_ptr_reduced->sto_xyz);

        smooth3(GMT, dims, voxelsize, s, 0, DT_FLOAT32);
        smooth3(PPM, dims, voxelsize, s, 0, DT_FLOAT32);

        /* Save GMT and PPM image */
        float *vol_reduced = (float *)malloc(sizeof(float) * out_ptr_reduced->nvox);

        subsample3(GMT, vol_reduced, dims, dims_reduced, DT_FLOAT32);
        if (!write_nifti_float(out_GMT, vol_reduced, DT_FLOAT32, slope, dims_reduced, voxelsize_reduced, out_ptr_reduced))
            exit(EXIT_FAILURE);

        subsample3(PPM, vol_reduced, dims, dims_reduced, DT_FLOAT32);
        if (!write_nifti_float(out_PPM, vol_reduced, DT_FLOAT32, slope, dims_reduced, voxelsize_reduced, out_ptr_reduced))
            exit(EXIT_FAILURE);

        free(vol_reduced);
    }
    else
    {
        if (!write_nifti_float(out_GMT, GMT, DT_FLOAT32, slope, dims, voxelsize, out_ptr))
            exit(EXIT_FAILURE);

        if (!write_nifti_float(out_PPM, PPM, DT_FLOAT32, slope, dims, voxelsize, out_ptr))
            exit(EXIT_FAILURE);
    }

    if (argc >= 6)
    {
        if (!write_nifti_float(out_CSD, dist_CSF, DT_FLOAT32, slope, dims, voxelsize, out_ptr))
            exit(EXIT_FAILURE);
        if (!write_nifti_float(out_WMD, dist_WM, DT_FLOAT32, slope, dims, voxelsize, out_ptr))
            exit(EXIT_FAILURE);
    }

    free(dist_CSF);
    free(dist_WM);
    free(GMT);
    free(PPM);
    free(src);

    return (EXIT_SUCCESS);
}
