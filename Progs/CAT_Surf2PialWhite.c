/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Surface reconstruction pipeline: Estimate pial and white matter surfaces
 * from a central surface using cortical thickness and volume information.
 *
 * Copyright Christian Gaser, University of Jena.
 */

#include <bicpl.h>
#include <ParseArgv.h>

#include "CAT_SurfaceIO.h"
#include "CAT_NiftiLib.h"
#include "CAT_SurfPialWhite.h"

/* -----------------------------------------------
 * Default arguments (map onto CAT_PialWhiteOptions)
 * ----------------------------------------------- */
double w1 = 0.05;
double w2 = 0.05;
double w3 = 0.05;
double sigma = 0.2;
int verbose = 0;
int iterations = 100;
int gradient_iterations = 00;
int method = 0;

/* Argument table for command-line parsing */
static ArgvInfo argTable[] = {
    {"-w1", ARGV_FLOAT, (char *)TRUE, (char *)&w1,
     "Set internal smoothness weight (w1)."},
    {"-w2", ARGV_FLOAT, (char *)TRUE, (char *)&w2,
     "Set gradient alignment weight (w2)."},
    {"-w3", ARGV_FLOAT, (char *)TRUE, (char *)&w3,
     "Set balloon force weight (w3)."},
    {"-sigma", ARGV_FLOAT, (char *)TRUE, (char *)&sigma,
     "Define sigma for smoothing the displacement field."},
    {"-iter", ARGV_INT, (char *)TRUE, (char *)&iterations,
     "Set number of deformation iterations."},
    {"-giter", ARGV_INT, (char *)TRUE, (char *)&gradient_iterations,
     "Set number of gradient refinement iterations (0 to disable)."},
    {"-method", ARGV_INT, (char *)TRUE, (char *)&method,
     "Method: 0 = deformation (default), 1 = Laplacian, 2 = ADE."},
    {"-verbose", ARGV_CONSTANT, (char *)TRUE, (char *)&verbose,
     "Enable verbose output."},
    {NULL, ARGV_END, NULL, NULL, NULL}};

/* -----------------------------------------------
 * Print usage/help text
 * ----------------------------------------------- */
static void
usage(const char *executable)
{
    fprintf(stderr,
            "\nUsage: %s surface_file thickness_file label_file "
            "output_pial_file output_white_file\n"
            "\n"
            "Estimate pial and white matter surfaces from a central surface using:\n"
            "- Cortical thickness values.\n"
            "- A label image that encodes tissue classes.\n\n"
            "This tool performs the following steps:\n"
            "1. Estimate preliminary pial and white surfaces using thickness.\n"
            "2. Smooth pial surface with curvature-guided correction.\n"
            "3. Perform joint deformation of both surfaces using the image\n"
            "   intensity and gradient field, while preserving topology and\n"
            "   maintaining surface distance (cortical thickness).\n\n"
            "Key deformation forces:\n"
            "  -w1     Internal smoothness term (e.g. 0.1).\n"
            "  -w2     Gradient alignment force (edges attraction).\n"
            "  -w3     Balloon force, based on isovalue distance.\n"
            "  -sigma  Controls displacement smoothing.\n"
            "  -iter   Number of iterations (e.g. 50).\n\n",
            executable);
}

/* -----------------------------------------------
 * Main entry point
 * ----------------------------------------------- */
int main(int argc, char *argv[])
{
    int n_objects, n_values;
    char *src_file, *pial_file, *white_file, *values_file, *label_file;
    float *labels;
    double *thickness_values;
    File_formats format;
    nifti_image *nii_ptr;
    polygons_struct *polygons;
    object_struct **object_list;

    /* Parse optional flags */
    if (ParseArgv(&argc, argv, argTable, 0))
    {
        usage(argv[0]);
        exit(EXIT_FAILURE);
    }

    initialize_argument_processing(argc, argv);

    /* Positional arguments */
    if (!get_string_argument(NULL, &src_file) ||
        !get_string_argument(NULL, &values_file) ||
        !get_string_argument(NULL, &label_file) ||
        !get_string_argument(NULL, &pial_file) ||
        !get_string_argument(NULL, &white_file))
    {
        usage(argv[0]);
        exit(EXIT_FAILURE);
    }

    /* Load central surface */
    if (input_graphics_any_format(src_file, &format, &n_objects,
                                  &object_list) != OK)
        exit(EXIT_FAILURE);

    if (n_objects > 1)
    {
        fprintf(stderr, "Error: Only one surface object allowed.\n");
        exit(EXIT_FAILURE);
    }
    polygons = get_polygons_ptr(object_list[0]);

    /* Load thickness values */
    if (input_values_any_format(values_file, &n_values,
                                &thickness_values) != OK)
        exit(EXIT_FAILURE);

    if (polygons->n_points != n_values)
    {
        fprintf(stderr,
                "Error: Number of surface vertices does not match "
                "number of thickness values.\n");
        exit(EXIT_FAILURE);
    }

    /* Load NIfTI label volume */
    nii_ptr = read_nifti_float(label_file, &labels, 0);
    if (!nii_ptr)
    {
        fprintf(stderr, "Error reading label volume: %s.\n", label_file);
        return EXIT_FAILURE;
    }

    /* Prepare output surface objects */
    object_struct *object_pial = create_object(POLYGONS);
    object_struct *object_white = create_object(POLYGONS);
    polygons_struct *pial_poly = get_polygons_ptr(object_pial);
    polygons_struct *white_poly = get_polygons_ptr(object_white);

    /* Fill option struct from command-line globals */
    CAT_PialWhiteOptions opts;
    CAT_PialWhiteOptionsInit(&opts);
    opts.w1 = w1;
    opts.w2 = w2;
    opts.w3 = w3;
    opts.sigma = sigma;
    opts.iterations = iterations;
    opts.gradient_iterations = gradient_iterations;
    opts.method = method;
    opts.verbose = verbose;

    /* Run the library estimation */
    if (CAT_SurfEstimatePialWhite(polygons, thickness_values, labels,
                                  nii_ptr, pial_poly, white_poly,
                                  &opts) != 0)
    {
        fprintf(stderr, "Error: Pial/white estimation failed.\n");
        exit(EXIT_FAILURE);
    }

    /* Save output surfaces */
    if (output_graphics_any_format(pial_file, format, 1,
                                   &object_pial, NULL) != OK ||
        output_graphics_any_format(white_file, format, 1,
                                   &object_white, NULL) != OK)
    {
        fprintf(stderr, "Error writing output surfaces.\n");
        exit(EXIT_FAILURE);
    }

    /* Cleanup */
    delete_object_list(n_objects, object_list);
    free(thickness_values);
    delete_object(object_pial);
    delete_object(object_white);

    return EXIT_SUCCESS;
}
