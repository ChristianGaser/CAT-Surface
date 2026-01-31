/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 *
 * CAT_VolLayerSmooth - Anisotropic smoothing along cortical layers
 *
 * This tool performs Laplace-guided anisotropic smoothing that smooths
 * data along cortical layers while preserving the laminar structure.
 *
 */

#include <float.h>
#include <stdlib.h>
#include <string.h>

#include <bicpl.h>
#include <ParseArgv.h>
#include "CAT_NiftiLib.h"
#include "CAT_Vol.h"
#include "CAT_LayerSmooth.h"

/* Argument defaults */
double fwhm = 3.0;        /* FWHM smoothing kernel in mm */
double extend = 0.0;      /* Extension outside cortical band in mm */
int verbose = 0;
int debug = 0;            /* Debug mode: save intermediate outputs */

static ArgvInfo argTable[] = {
    {"-fwhm", ARGV_FLOAT, (char *)1, (char *)&fwhm,
     "Smoothing kernel size in mm (FWHM). Default: 3.0 mm.\n\t\t   This controls the spatial extent of smoothing along the cortical surface."},
    
    {"-extend", ARGV_FLOAT, (char *)1, (char *)&extend,
     "Extension outside cortical band in mm. Default: 0.0 (strict GM only).\n\t\t   Setting this > 0 allows smoothing to include voxels slightly outside\n\t\t   the strict GM band (e.g., partial volume voxels in CSF or WM)."},
    
    {"-verbose", ARGV_CONSTANT, (char *)1, (char *)&verbose,
     "Enable verbose mode for detailed output during processing."},
    
    {"-debug", ARGV_CONSTANT, (char *)1, (char *)&debug,
     "Debug mode: save intermediate outputs (depth field, WM/CSF distances).\n\t\t   Files are saved with prefixes: depth_, distWM_, distCSF_"},
    
    {NULL, ARGV_END, NULL, NULL, NULL}
};

static void usage(char *executable)
{
    char *usage_str = "\n\
Usage: %s [options] <pve_label.nii> <input.nii> <output.nii>\n\
\n\
    This program performs anisotropic smoothing along cortical layers.\n\
    It smooths data within iso-depth contours (parallel to cortical surface)\n\
    while preserving the laminar structure and avoiding smoothing across\n\
    GM/WM or GM/CSF boundaries.\n\
    \n\
    The algorithm:\n\
    1. Computes cortical depth field from the PVE segmentation:\n\
       phi = WM_distance / (WM_distance + CSF_distance)\n\
       This gives values from 0 (at WM surface) to 1 (at pial surface).\n\
    \n\
    2. For each cortical voxel, applies Gaussian smoothing weighted by:\n\
       - Euclidean distance (standard spatial Gaussian kernel)\n\
       - Depth similarity: exp(-(delta_phi)^2 / (2*sigma_depth^2))\n\
    \n\
    3. Neighbors with different cortical depths contribute less to smoothing,\n\
       effectively smoothing along layers rather than across them.\n\
    \n\
    This approach is similar to LAYNII's LN_LAYER_SMOOTH but uses a continuous\n\
    depth field instead of discrete layers for smoother transitions.\n\
\n\
    Required inputs:\n\
        pve_label.nii  - PVE segmentation image (1=CSF, 2=GM, 3=WM with partial volumes)\n\
        input.nii      - Image to be smoothed (same dimensions as PVE label)\n\
        output.nii     - Output smoothed image\n\
\n\
Options:\n";

    fprintf(stderr, usage_str, executable);
}

int main(int argc, char *argv[])
{
    char *pve_file, *input_file, *output_file;
    int dims[3], dims_seg[3];
    float *seg, *data;
    float *depth = NULL, *dist_WM = NULL, *dist_CSF = NULL;
    double voxelsize[3];
    nifti_image *seg_ptr, *data_ptr;
    
    /* Initialize argument processing */
    initialize_argument_processing(argc, argv);
    
    /* Parse arguments */
    if (ParseArgv(&argc, argv, argTable, 0) || (argc < 4)) {
        usage(argv[0]);
        fprintf(stderr, "     %s -help\n\n", argv[0]);
        exit(EXIT_FAILURE);
    }
    
    /* Get filenames */
    pve_file = argv[1];
    input_file = argv[2];
    output_file = argv[3];
    
    /* Validate FWHM */
    if (fwhm <= 0.0) {
        fprintf(stderr, "Error: FWHM must be positive (got %g)\n", fwhm);
        exit(EXIT_FAILURE);
    }
    
    /* Read PVE segmentation */
    if (verbose) fprintf(stderr, "Reading PVE segmentation: %s\n", pve_file);
    seg_ptr = read_nifti_float(pve_file, &seg, 0);
    if (!seg_ptr) {
        fprintf(stderr, "Error reading PVE file: %s\n", pve_file);
        exit(EXIT_FAILURE);
    }
    
    /* Get dimensions from segmentation */
    dims_seg[0] = seg_ptr->nx;
    dims_seg[1] = seg_ptr->ny;
    dims_seg[2] = seg_ptr->nz;
    
    /* Read input image */
    if (verbose) fprintf(stderr, "Reading input image: %s\n", input_file);
    data_ptr = read_nifti_float(input_file, &data, 0);
    if (!data_ptr) {
        fprintf(stderr, "Error reading input file: %s\n", input_file);
        exit(EXIT_FAILURE);
    }
    
    /* Get dimensions and voxel size from input */
    dims[0] = data_ptr->nx;
    dims[1] = data_ptr->ny;
    dims[2] = data_ptr->nz;
    voxelsize[0] = data_ptr->dx;
    voxelsize[1] = data_ptr->dy;
    voxelsize[2] = data_ptr->dz;
    
    /* Check dimension consistency */
    if (dims[0] != dims_seg[0] || dims[1] != dims_seg[1] || dims[2] != dims_seg[2]) {
        fprintf(stderr, "Error: Dimension mismatch between PVE (%dx%dx%d) and input (%dx%dx%d)\n",
                dims_seg[0], dims_seg[1], dims_seg[2], dims[0], dims[1], dims[2]);
        exit(EXIT_FAILURE);
    }
    
    if (verbose) {
        fprintf(stderr, "Image dimensions: %d x %d x %d\n", dims[0], dims[1], dims[2]);
        fprintf(stderr, "Voxel size: %.3f x %.3f x %.3f mm\n", voxelsize[0], voxelsize[1], voxelsize[2]);
        fprintf(stderr, "Smoothing FWHM: %.2f mm\n", fwhm);
        if (extend > 0.0) {
            fprintf(stderr, "Extension outside GM: %.2f mm\n", extend);
        }
        if (debug) {
            fprintf(stderr, "Debug mode: will save intermediate outputs\n");
        }
    }
    
    /* In debug mode, allocate arrays for intermediate outputs */
    if (debug) {
        int nvox = dims[0] * dims[1] * dims[2];
        depth = (float *)malloc(sizeof(float) * nvox);
        dist_WM = (float *)malloc(sizeof(float) * nvox);
        dist_CSF = (float *)malloc(sizeof(float) * nvox);
        
        if (!depth || !dist_WM || !dist_CSF) {
            fprintf(stderr, "Error: Memory allocation failed for debug arrays\n");
            exit(EXIT_FAILURE);
        }
        
        /* Compute depth field and distances */
        if (verbose) fprintf(stderr, "Computing cortical depth field...\n");
        compute_cortical_depth(seg, depth, dist_WM, dist_CSF, dims, voxelsize, extend);
        
        /* Save debug outputs */
        char debug_filename[1024];
        char *base = strrchr(output_file, '/');
        char *dir_part = "";
        char base_copy[1024];
        
        if (base) {
            strncpy(base_copy, output_file, base - output_file + 1);
            base_copy[base - output_file + 1] = '\0';
            dir_part = base_copy;
            base++;
        } else {
            base = output_file;
        }
        
        /* Save depth field */
        snprintf(debug_filename, sizeof(debug_filename), "%sdepth_%s", dir_part, base);
        if (verbose) fprintf(stderr, "Saving depth field: %s\n", debug_filename);
        if (!write_nifti_float(debug_filename, depth, DT_FLOAT32, 1.0, dims, voxelsize, data_ptr)) {
            fprintf(stderr, "Error writing depth field: %s\n", debug_filename);
        }
        
        /* Save WM distance */
        snprintf(debug_filename, sizeof(debug_filename), "%sdistWM_%s", dir_part, base);
        if (verbose) fprintf(stderr, "Saving WM distance: %s\n", debug_filename);
        if (!write_nifti_float(debug_filename, dist_WM, DT_FLOAT32, 1.0, dims, voxelsize, data_ptr)) {
            fprintf(stderr, "Error writing WM distance: %s\n", debug_filename);
        }
        
        /* Save CSF distance */
        snprintf(debug_filename, sizeof(debug_filename), "%sdistCSF_%s", dir_part, base);
        if (verbose) fprintf(stderr, "Saving CSF distance: %s\n", debug_filename);
        if (!write_nifti_float(debug_filename, dist_CSF, DT_FLOAT32, 1.0, dims, voxelsize, data_ptr)) {
            fprintf(stderr, "Error writing CSF distance: %s\n", debug_filename);
        }
        
        free(depth);
        free(dist_WM);
        free(dist_CSF);
    return(1);
    }
    
    /* Perform layer-guided smoothing */
    if (verbose) fprintf(stderr, "Performing layer-guided anisotropic smoothing...\n");
    smooth_within_cortex_float(data, seg, dims, voxelsize, fwhm, extend);
    
    /* Write output */
    if (verbose) fprintf(stderr, "Writing output: %s\n", output_file);
    if (!write_nifti_float(output_file, data, DT_FLOAT32, 1.0, dims, voxelsize, data_ptr)) {
        fprintf(stderr, "Error writing output file: %s\n", output_file);
        exit(EXIT_FAILURE);
    }
    
    /* Cleanup */
    free(seg);
    free(data);
    
    /* Clear data pointers before freeing nifti structs to avoid double-free */
    seg_ptr->data = NULL;
    data_ptr->data = NULL;
    
    nifti_image_free(seg_ptr);
    nifti_image_free(data_ptr);
    
    if (verbose) fprintf(stderr, "Done.\n");
    
    return EXIT_SUCCESS;
}
