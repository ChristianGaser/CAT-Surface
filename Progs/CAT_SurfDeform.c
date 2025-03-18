/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry, University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 *
 * Surface Deformation Algorithm with Diffeomorphic Constraints
 */

#include <ParseArgv.h>
#include "CAT_SurfaceIO.h"
#include "CAT_Deform.h"
#include "CAT_NiftiLib.h"

/* Default parameter values */
double w1 = 0.1;  // Internal smoothness force
double w2 = 0.1;   // Gradient alignment force
double w3 = 1.0;   // Balloon force (expansion/contraction)
double isovalue = 0.5;  // Target isosurface value
double sigma = 0.2;
int iterations = 100;    // Number of iterations
int verbose = 0;         // Verbose mode

/* Command-line argument table */
static ArgvInfo argTable[] = {
    {"-isovalue", ARGV_FLOAT, (char *) TRUE, (char *) &isovalue, "Define isovalue (target intensity)."},
    {"-w1", ARGV_FLOAT, (char *) TRUE, (char *) &w1, "Set internal smoothness weight (w1)."},
    {"-w2", ARGV_FLOAT, (char *) TRUE, (char *) &w2, "Set gradient alignment weight (w2)."},
    {"-w3", ARGV_FLOAT, (char *) TRUE, (char *) &w3, "Set balloon force weight (w3)."},
    {"-sigma", ARGV_FLOAT, (char *) TRUE, (char *) &sigma, "Define sigma for smoothing the displacement field."},
    {"-iter", ARGV_INT, (char *) TRUE, (char *) &iterations, "Set number of deformation iterations."},
    {"-verbose", ARGV_CONSTANT, (char *) TRUE, (char *) &verbose, "Enable verbose output."},
    {NULL, ARGV_END, NULL, NULL, NULL}
};

/**
 * @brief Displays usage instructions for the surface deformation tool.
 *
 * This function provides an overview of the deformation process and guides users
 * on how to correctly execute the program.
 *
 * **Approach Overview:**
 *  - The method deforms a 3D mesh based on image intensity and gradient forces.
 *  - Ensures a balance between smoothness, gradient attraction, and volume expansion.
 *  - Enforces diffeomorphic constraints to prevent self-intersections.
 *
 * **Usage:**
 * ```
 *  CAT_SurfDeform [options] <volume.nii> <input_surface.gii> <output_surface.gii>
 * ```
 *
 * **Arguments:**
 *  - `<volume.nii>`: 3D image file (NIfTI format).
 *  - `<input_surface.gii>`: Initial mesh surface file.
 *  - `<output_surface.gii>`: Output deformed mesh surface file.
 *
 * **Options:**
 *  - `-isovalue <float>` → Target isosurface value (default: `0.5`).
 *  - `-w1 <float>` → Internal smoothing weight (`0.01`).
 *  - `-w2 <float>` → Gradient alignment weight (`0.1`).
 *  - `-w3 <float>` → Balloon force weight (`0.2`).
 *  - `-iter <int>` → Number of iterations (`100`).
 *  - `-verbose` → Enable detailed output.
 */
void usage(char *executable) {
    char *usage_str = "\n"
    "--------------------------------------------------------------------------------\n"
    " Surface Deformation Tool - Mesh Adaptation Based on Image Intensity\n"
    "--------------------------------------------------------------------------------\n"
    "\n"
    "Usage:\n"
    "  %s [options] <volume.nii> <input_surface.gii> <output_surface.gii>\n"
    "\n"
    "Description:\n"
    "  This program deforms a 3D surface mesh using an image-driven approach.\n"
    "  It balances smoothness, edge attraction, and volume preservation while\n"
    "  ensuring a diffeomorphic transformation (preventing self-intersections).\n"
    "  The latter is ensured by smoothing the displacement field at each iteration.\n"
    "\n"
    "Required Arguments:\n"
    "  <volume.nii>         3D image file (NIfTI format) used for deformation.\n"
    "  <input_surface.gii>  Input mesh surface file.\n"
    "  <output_surface.gii> Output deformed mesh surface file.\n"
    "\n"
    "Optional Parameters:\n"
    "  -isovalue <float>    Target isosurface intensity value (default: 0.5).\n"
    "  -w1 <float>          Internal smoothing weight (default: 0.1).\n"
    "  -w2 <float>          Gradient alignment weight (default: 0.1).\n"
    "  -w3 <float>          Balloon force weight (default: 1.0).\n"
    "  -sigma <float>       Smoothing the displacement field at each iteration (default: 0.2).\n"
    "  -iter <int>          Number of deformation iterations (default: 100).\n"
    "  -verbose             Enable verbose mode for detailed logs.\n"
    "\n"
    "Deformation Forces:\n"
    "  - w1: Encourages local smoothness by averaging neighboring vertices.\n"
    "  - w2: Attracts the surface toward edges in the image (gradient force).\n"
    "  - w3: Expands or contracts the surface, assisting in volume growth/shrinkage.\n"
    "\n"
    "Example Usage:\n"
    "  %s -isovalue 0.6 -w1 0.1 -w2 0.1 -w3 1.0 -iter 150 \\\n"
    "     brain.nii.gz white_surface.gii white_surface_deformed.gii\n"
    "\n"
    "Expected Output:\n"
    "  - The output surface mesh is deformed according to the image intensities.\n"
    "  - A valid diffeomorphic transformation is maintained to minimize self-intersections.\n"
    "--------------------------------------------------------------------------------\n";

    fprintf(stderr, usage_str, executable, executable);
}
int
main(int argc, char *argv[])
{
    char *volume_file = NULL;
    char *input_file = NULL, *output_surface_file = NULL;
    float *input;
    int n_objects;
    File_formats file_format;
    object_struct **object_list;
    polygons_struct *polygons;
    nifti_image *nii_ptr;

    initialize_argument_processing(argc, argv);

    /* get the arguments from the command line */
    if (ParseArgv(&argc, argv, argTable, 0)) {
        usage(argv[0]);
        fprintf(stderr, "     %s -help\n\n", argv[0]);
        return EXIT_FAILURE;
    }

    if (!get_string_argument(NULL, &volume_file) || 
        !get_string_argument(NULL, &input_file)  ||
        !get_string_argument(NULL, &output_surface_file)) {
        usage(argv[0]);
        fprintf(stderr, "Usage: CAT_SurfDeform [options] volume_file surface_file output_surface_file\n");
        return EXIT_FAILURE;
    }

    /* read first image to get image parameters */
    nii_ptr = read_nifti_float(volume_file, &input, 0);
    if(nii_ptr == NULL) {
        fprintf(stderr,"Error reading %s.\n", volume_file);
        return(EXIT_FAILURE);
    }

    if (input_graphics_any_format(input_file, &file_format,
              &n_objects, &object_list) == ERROR ||
              n_objects != 1 || object_list[0]->object_type != POLYGONS) {
        fprintf(stderr, "File must contain 1 polygons struct.\n");
        exit(EXIT_FAILURE);
    }

    polygons = get_polygons_ptr(object_list[0]);    

    double weights[3] = {w1, w2, w3};
    surf_deform(polygons, input, nii_ptr, weights, sigma, isovalue, iterations, verbose);
    
    if (output_graphics_any_format(output_surface_file, ASCII_FORMAT,
                     n_objects, object_list, NULL) == ERROR)
        exit(EXIT_FAILURE);

    return(EXIT_SUCCESS);
}
