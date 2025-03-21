/**
 * CAT_AddValuesToSurf.c
 * This program adds values to a surface file in the gifti format.
 *
 * The program takes three command line arguments: the input surface file, the input values file,
 * and the output surface file. It reads the surface file and the values file, and then adds the values
 * to the surface. The resulting surface with values is saved in the gifti format.
 *
 * Usage: ./CAT_AddValuesToSurf surface_file values_file output_surface_file.gii
 *
 */

/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 *
 */

#include <bicpl.h>
#include <float.h>

#include "CAT_SurfaceIO.h"

void
usage(char *executable)
{
    char *usage_str = "\n\
Usage: %s  surface_file values_file output_surface_file.gii\n\n";

    fprintf(stderr, usage_str, executable);
}
  
int
main(int argc, char *argv[])
{
    double         value, *values;
    Status         status;
    char         *src_file, *dest_file, *values_file;
    int          i, p, n_objects, n_pts, n_values;
    Point        *pts;
    File_formats     format;
    object_struct    **object_list;

    initialize_argument_processing(argc, argv);

    if (!get_string_argument( NULL, &src_file) ||
      !get_string_argument( NULL, &values_file) ||
      !get_string_argument( NULL, &dest_file)) {
        usage(argv[0]);
        exit(EXIT_FAILURE);
    }

    if (filename_extension_matches(dest_file, "gii") != 1) {
        fprintf(stderr,"Only gifti output format allowed (use .gii as extension).\n");
        exit(EXIT_FAILURE);
    }
    
    if (input_graphics_any_format(src_file, &format, &n_objects,
                    &object_list) != OK)
        exit(EXIT_FAILURE);

    if (input_values_any_format(values_file, &n_values, &values) != OK)
        exit(EXIT_FAILURE);
    
    if (n_objects > 1) {
        fprintf(stderr,"Only one object allowed.\n");
        exit(EXIT_FAILURE);
    }

    n_pts = get_object_points(object_list[0], &pts);
    if (n_pts != n_values) {
        fprintf(stderr,"Number of points differs from number of values.\n");
        exit(EXIT_FAILURE);
    }

    status = output_graphics_any_format(dest_file, format,
                      n_objects, object_list, values);

    return(status != OK);
}
