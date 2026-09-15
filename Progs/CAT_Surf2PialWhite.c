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
int gradient_iterations = 0;
int method = 0;
int remove_intersect = 0;
int legacy_pial = 0;
/* Profile placement: negative values keep the library defaults */
double valley_depth = -1.0;
double search_out = -1.0;
int pial_iterations = -1;

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
     "Method: 0 = deformation (default), 1 = ADE, 2 = deformation:pial | ADE:white."},
    {"-legacy-pial", ARGV_CONSTANT, (char *)TRUE, (char *)&legacy_pial,
     "Deform the pial surface with the balloon-force deformation instead of\n\
                 placing it by profile search."},
    {"-valley-depth", ARGV_FLOAT, (char *)TRUE, (char *)&valley_depth,
     "Rise above the running minimum of the label profile that ends a valley\n\
                 (glued sulcus) in profile placement (default 0.05)."},
    {"-pial-search", ARGV_FLOAT, (char *)TRUE, (char *)&search_out,
     "Outward search distance along the normal in mm for profile placement\n\
                 (default 2.0)."},
    {"-pial-iter", ARGV_INT, (char *)TRUE, (char *)&pial_iterations,
     "Number of iterations of profile placement (default 60)."},
    {"-remove_intersect", ARGV_CONSTANT, (char *)TRUE, (char *)&remove_intersect,
     "Remove self-intersections of the resulting pial and white surfaces.\n\
                 The mesh topology is preserved, i.e. both surfaces keep their\n\
                 vertex correspondence with the central surface."},
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
            "3. Deform the white surface using the image intensity and\n"
            "   gradient field.\n"
            "4. Place the pial surface by searching the label profile along\n"
            "   each normal for the CSF/GM boundary (1.5), or for the bottom\n"
            "   of the valley where a glued sulcus never reaches it.  Facing\n"
            "   sulcal walls meet in the middle.  Use -legacy-pial to deform\n"
            "   the pial surface together with the white surface instead.\n\n"
            "Key deformation forces:\n"
            "  -w1     Internal smoothness term (e.g. 0.1).\n"
            "  -w2     Gradient alignment force (edges attraction).\n"
            "  -w3     Balloon force, based on isovalue distance.\n"
            "  -sigma  Controls displacement smoothing.\n"
            "  -method Controls general approach (ADE or deformation).\n"
            "  -iter   Number of iterations (e.g. 50).\n\n"
            "Use -remove_intersect to repair self-intersections of the resulting\n"
            "pial and white surfaces.  This preserves the mesh topology, so the\n"
            "vertex correspondence with the central surface (and thus the\n"
            "per-vertex thickness) stays valid.\n\n",
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
    opts.remove_intersect = remove_intersect;
    opts.pial_profile = !legacy_pial;
    if (valley_depth >= 0.0)
        opts.profile.valley_depth = valley_depth;
    if (search_out >= 0.0)
        opts.profile.search_out = search_out;
    if (pial_iterations >= 0)
        opts.profile.iterations = pial_iterations;
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
