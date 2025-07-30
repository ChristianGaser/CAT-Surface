/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 *
 */

/*
 * Laplacian smooth with Humphrey’s Classes to preserve volume: 
 * https://doi.org/10.1111/1467-8659.00334
 */

#include <bicpl.h>
#include <ParseArgv.h>

#include "CAT_SurfaceIO.h"
#include "CAT_Surf.h"

/* argument defaults */
int iter = 10;
double alpha = 0.1;
double beta  = 0.5;

/* the argument table */
static ArgvInfo argTable[] = {
  {"-iter", ARGV_INT, (char *) TRUE, (char *) &iter,
   "Iterations for Laplacian smoothing."},
  {"-alpha", ARGV_FLOAT, (char *) TRUE, (char *) &alpha,
   "alpha."},
  {"-beta", ARGV_FLOAT, (char *) TRUE, (char *) &beta,
   "beta."},
   {NULL, ARGV_END, NULL, NULL, NULL}
};

void
usage(char *executable)
{
    char *usage_str = "\n\
Usage: %s [options] surface_file surface_output_file\n\n\
   Improved Laplacian Smoothing of Noisy Surface Meshes\n\
   https://doi.org/10.1111/1467-8659.00334\n\n";
    fprintf(stderr, usage_str, executable);
}

int
main(int argc, char *argv[])
{
    char       *input_file, *output_surface_file;
    int        i, p, n_objects;
    File_formats   format;
    object_struct  **object_list;
    polygons_struct  *polygons, *polygons_orig;

    /* get the arguments from the command line */
    if (ParseArgv(&argc, argv, argTable, 0)) {
        usage(argv[0]);
        fprintf(stderr, "   %s -help\n\n", argv[0]);
        exit(EXIT_FAILURE);
    }

    initialize_argument_processing(argc, argv);

    if (!get_string_argument(NULL, &input_file) ||
      !get_string_argument(NULL, &output_surface_file)) {
        usage(argv[0]);
        fprintf(stderr, "   %s -help\n\n", argv[0]);
        exit(EXIT_FAILURE);
    }

    if (input_graphics_any_format(input_file, &format, &n_objects,
                    &object_list) != OK ||
      n_objects != 1 || get_object_type(object_list[0]) != POLYGONS) {
        fprintf(stderr, "Error reading %s.\n", input_file);
        exit(EXIT_FAILURE);
    }

    polygons = get_polygons_ptr(object_list[0]);

    if (smooth_laplacian(polygons, iter, alpha, beta) != OK) {
        fprintf(stderr, "Error smoothing %s.\n", input_file);
        exit(EXIT_FAILURE);
    }

    compute_polygon_normals(polygons);

    if(output_graphics_any_format(output_surface_file, format, 1, 
            object_list, NULL) != OK)
        exit(EXIT_FAILURE);

    return(EXIT_SUCCESS);
}
