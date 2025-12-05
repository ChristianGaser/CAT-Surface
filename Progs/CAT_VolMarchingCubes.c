/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 *
 */

#include <ParseArgv.h>
#include "CAT_MarchingCubes.h"
#include "CAT_SurfaceIO.h"

/* argument defaults */
char *label_filename = NULL;
double min_threshold = 0.5;
double pre_fwhm = 2.0;
double dist_morph = FLT_MAX;
double strength_gyri_mask = 0.1;
int iter_laplacian = 50;
int n_median_filter = 2;
int verbose = 0;
int n_iter = 10;
int fast = 0;

/* the argument table */
static ArgvInfo argTable[] = {
  {"-label", ARGV_STRING, (char *) NULL, (char *) &label_filename, 
    "File containing segmentation labels for creating smooth mask of gyral\n\
    and sulcal areas. This prevents sulcal closure by using a higher isovalue\n\
    in sulci, and prevent cutting gyri by using a lower isovalue in gyri."},

  {"-strength-gyrimask", ARGV_FLOAT, (char *) TRUE, (char *) &strength_gyri_mask,
    "Define the strength of isovalue correction using the smooth.\n\
     gyri/sulci mask from the label map."},
  
  {"-thresh", ARGV_FLOAT, (char *) TRUE, (char *) &min_threshold,
    "Define the volume threshold, also known as the isovalue.\n\
     This value is crucial for initial image thresholding."},
  
  {"-pre-fwhm", ARGV_FLOAT, (char *) TRUE, (char *) &pre_fwhm,
    "Specify the Full Width Half Maximum (FWHM) for the preprocessing\n\
     smoothing filter. This helps in preserving gyri and sulci by\n\
     creating a weighted average between original and smoothed\n\
     images based on the gradient of the input image. Areas with \n\
     topology artefacts are often characterized by large gradients,\n\
     thus smoothing in these areas tries to prevent these artefacts."},

  {"-iter-laplacian", ARGV_INT, (char *) TRUE, (char *) &iter_laplacian,
    "Set number of iterations for Laplacian surface smoothing. This aids\n\
     in removing noise from the mesh"},
  
  {"-dist-morph", ARGV_FLOAT, (char *) TRUE, (char *) &dist_morph,
    "Apply initial morphological opening or closing step. Closing is used\n\
     by a value around 1.0 and opening by negative values around -1.0.\n\
     The default automatically estimates the optimal value"},
  
  {"-median-filter", ARGV_INT, (char *) TRUE, (char *) &n_median_filter,
    "Specify the number of iterations to apply a median filter to areas\n\
     where the gradient of the thresholded image indicates larger clusters.\n\
     These clusters may point to potential topology artifacts and regions\n\
     with high local variations. This process helps to smooth these areas, \n\
     improving the quality of the surface reconstruction in subsequent steps."},
  
  {"-iter", ARGV_INT, (char *) TRUE, (char *) &n_iter,
    "Number of iterations."},
  
  {"-fast", ARGV_CONSTANT, (char *) TRUE, (char *) &fast,
    "Enable fast processing without any preprocessing, smoothing, or topology \n\
    correction. Only the final Laplacian surface smoothing is applied if defined."},

  {"-verbose", ARGV_CONSTANT, (char *) TRUE, (char *) &verbose,
    "Enable verbose mode for detailed output during processing."},

    {NULL, ARGV_END, NULL, NULL, NULL}
};


private void
usage(
    char *executable)
{
    char *usage_str = "\n\
Usage: CAT_VolMarchingCubes [options] input.nii output_surface_file [change_map.nii]\n\
\n\
    This method generates a mesh with an Euler number of 2 (genus 0) from the\n\
    thresholded volume. The process involves:\n\
    \n\
    1. **Preprocessing with Smoothing Filter:**\n\
       - Apply a smoothing filter to the input image to remove outliers.\n\
       - Use a weighted average of the original and smoothed images to\n\
         preserve gyri and sulci.\n\
       - Weighting is based on the gradient of the input image.\n\
       - Weights range from 0 (areas with low gradient) to 1 (areas\n\
         with large gradient), with intermediate values scaled linearly.\n\
       - Weighting effect is enhanced by squaring the value.\n\
    \n\
    2. **Preprocessing with Median Filter:**\n\
       - Apply an iterative median filter to remove noise.\n\
    \n\
    3. **Optionally Use of smooth mask of gyri and sulci:**\n\
        - This prevents sulcal closure by using a higher isovalue in sulci,\n\
          and prevent cutting gyri by using a lower isovalue in gyri\n\
          An optional label maps must be defined to use this approach.\n\
    \n\
    4. **Morphological Opening:**\n\
       - Apply additional morphological opening or closing, defined by\n\
         `-dist-morph`, to minimize changes due to topology correction.\n\
       - Closing is used for positive values (e.g. 1.0) and opening\n\
         for negative values. The default is to automatically estimate\n\
         the optimal value to minimize issues due to topology correction.\n\
    \n\
    5. **Extraction of the Largest Component:**\n\
       - Extract the largest component for further processing.\n\
    \n\
    6. **Mesh Smoothing:**\n\
       - Smooth the extracted mesh with a Laplacian filter.\n\n";

    fprintf(stderr,"%s\n %s\n",usage_str, executable);
}

int main(int argc, char *argv[]) {
    float *input_float, *label;
    char out_diff[1024];
    char *input_filename, *output_filename;
    object_struct *object;

    /* get the arguments from the command line */
    if (ParseArgv(&argc, argv, argTable, 0)) {
        usage(argv[0]);
        fprintf(stderr, "     %s -help\n\n", argv[0]);
        return EXIT_FAILURE;
    }

    initialize_argument_processing(argc, argv);

    if (!get_string_argument(NULL, &input_filename) || !get_string_argument(NULL, &output_filename)) {
        usage(argv[0]);
        fprintf(stderr, "Usage: CAT_VolMarchingCubes input.nii output_surface_file [options] [change_map.nii]\n");
        return EXIT_FAILURE;
    }

    /* Read input volume */
    nifti_image *nii_ptr = read_nifti_float(input_filename, &input_float, 0);
    if (!nii_ptr) {
        fprintf(stderr, "Error reading %s.\n", input_filename);
        return EXIT_FAILURE;
    }

    if (label_filename) {
        nifti_image *nii_ptr2 = read_nifti_float(label_filename, &label, 0);
        if (!nii_ptr2) {
            fprintf(stderr, "Error reading %s.\n", label_filename);
            return EXIT_FAILURE;
        }
        if ((nii_ptr->nx != nii_ptr2->nx) ||
            (nii_ptr->ny != nii_ptr2->ny) ||
            (nii_ptr->nz != nii_ptr2->nz)) {
            fprintf(stderr, "Error: Label image must have the same dimensions as the input image.\n");
            return EXIT_FAILURE;
        }
    } else label = NULL;

    if (fast) {
        object = apply_marching_cubes_fast(input_float, nii_ptr, 
                    min_threshold, iter_laplacian, verbose);
    } else {
        object = apply_marching_cubes(input_float, nii_ptr, label, 
                    min_threshold, pre_fwhm, iter_laplacian, dist_morph, n_median_filter, 
                    n_iter, strength_gyri_mask, verbose);
    }
    if (object) {
        output_graphics_any_format(output_filename, ASCII_FORMAT, 1, &object, NULL);
    } else {
        fprintf(stderr, "Error generating surface.\n");
    }

    if(argc > 3) {
        double voxelsize[N_DIMENSIONS];
        int dims[MAX_DIMENSIONS];
        dims[0] = nii_ptr->nx;
        dims[1] = nii_ptr->ny;
        dims[2] = nii_ptr->nz;
        voxelsize[0] = nii_ptr->dx;
        voxelsize[1] = nii_ptr->dy;
        voxelsize[2] = nii_ptr->dz;
        (void) sprintf(out_diff, "%s", argv[3]); 
        if (!write_nifti_float(out_diff, input_float, DT_FLOAT32, 1.0, dims, voxelsize, nii_ptr)) 
            exit(EXIT_FAILURE);
    }

    free(input_float);
    delete_marching_cubes_table();
    return 0;
}
