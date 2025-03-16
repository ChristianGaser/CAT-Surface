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

/* the argument table */
static ArgvInfo argTable[] = {
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

  {"-post-fwhm", ARGV_FLOAT, (char *) TRUE, (char *) &post_fwhm,
    "Set FWHM for surface smoothing. This aids in correcting the mesh\n\
     in folded areas like gyri and sulci. Note: Do not use smoothing\n\
     size > 3 mm for reliable compensation in these areas."},
  
  {"-dist-morph", ARGV_FLOAT, (char *) TRUE, (char *) &dist_morph,
    "Apply initial morphological open or close step. Close is used\n\
     by a value around 1.0 and open by negative values around -1.0.\n\
     The default automatically estimates the optimal value"},
  
  {"-median-filter", ARGV_INT, (char *) TRUE, (char *) &n_median_filter,
    "Specify the number of iterations to apply a median filter to areas\n\
     where the gradient of the thresholded image indicates larger clusters.\n\
     These clusters may point to potential topology artifacts and regions\n\
     with high local variations. This process helps to smooth these areas, \n\
     improving the quality of the surface reconstruction in subsequent steps."},
  
  {"-iter", ARGV_INT, (char *) TRUE, (char *) &n_iter,
    "Number of iterations."},
  
  {"-local-smoothing", ARGV_FLOAT, (char *) TRUE, (char *) &local_smoothing,
    "Apply local surface smoothing to resulting surface in areas where the distance\n\
     between the surface and a shifted surface is below the expected distance,\n\
     which often happens due to self intersections of the surface."},
  
  {"-verbose", ARGV_CONSTANT, (char *) TRUE, (char *) &verbose,
    "Enable verbose mode for detailed output during processing."},

    {NULL, ARGV_END, NULL, NULL, NULL}
};


private void
usage(
    char *executable)
{
    char *usage_str = "\n\
Usage: CAT_VolMarchingCubes input.nii output_surface_file [change_map.nii thresholded_map.nii]\n\
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
       - Apply an iterative median filter to areas where the gradient of \n\
         the thresholded image indicates larger clusters.\n\
       - Use a weighted average of the original and median filterd images.\n\
       - Weighting is estimated using gradient of the input image and.\n\
         morphological operations to find larger clusters\n\
    \n\
    3. **Morphological Opening:**\n\
       - Apply additional morphological opening, scaled by `scl_open`,\n\
         to prevent gyri fusion and minimize local artifacts.\n\
       - Opening strength is determined by analyzing the impact of\n\
         different `scl_open` values and tracking RMSE changes.\n\
    \n\
    4. **Extraction of the Largest Component:**\n\
       - Extract the largest component for further processing.\n\
    \n\
    5. **Mesh Smoothing:**\n\
       - Smooth the extracted mesh.\n\
    \n\
    6. **Mesh Correction in Folded Areas:**\n\
       - Correct the mesh in areas with folds, particularly in gyri and\n\
         sulci, to counterbalance the averaging effect from smoothing.\n\
       - Use mean curvature average as a folding measure to estimate\n\
         necessary compensation.\n\
       - Compensation degree is auto-calculated based on deviation\n\
         from the defined isovalue.\n\
    7. **Mesh Correction in Areas with Self Intersections:**\n\
       - Apply local surface smoothing to resulting surface in areas where\n\
         the distance between the surface and a shifted surface is below \n\
         the expected distance, which often happens due to self intersections\n\
         of the surface.\n\n";

    fprintf(stderr,"%s\n %s\n",usage_str, executable);
}

int main(int argc, char *argv[]) {
    float *input_float;
    char out_diff[1024];

    initialize_argument_processing(argc, argv);

    /* get the arguments from the command line */
    if (ParseArgv(&argc, argv, argTable, 0)) {
        usage(argv[0]);
        fprintf(stderr, "     %s -help\n\n", argv[0]);
        return EXIT_FAILURE;
    }

    char *input_filename, *output_filename;
    if (!get_string_argument(NULL, &input_filename) || !get_string_argument(NULL, &output_filename)) {
        usage(argv[0]);
        fprintf(stderr, "Usage: CAT_VolMarchingCubes input.nii output_surface_file\n");
        return EXIT_FAILURE;
    }

    /* Read input volume */
    nifti_image *nii_ptr = read_nifti_float(input_filename, &input_float, 0);
    if (!nii_ptr) {
        fprintf(stderr, "Error reading %s.\n", input_filename);
        return EXIT_FAILURE;
    }

    object_struct *object = apply_marching_cubes(input_float, nii_ptr);
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
