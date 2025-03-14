/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 *
 */

#include <bicpl.h>
#include <ParseArgv.h>

#include "CAT_Surf.h"
#include "CAT_Curvature.h"
#include "CAT_SurfaceIO.h"
  
BOOLEAN use_sqrt = 0;

/* the argument table */
ArgvInfo argTable[] = {
  { "-sqrt", ARGV_CONSTANT, (char *) 1, (char *) &use_sqrt,
  "Obtain sqrt-transformed values to render data more normally distributed." },
  { NULL, ARGV_END, NULL, NULL, NULL }
};

void
usage(char *executable)
{
    char *usage_str = "\n\
Usage: %s  surface_file sphere_file output_values_file [-sqrt]\n\n\
   Calculate sulcus depth based on the euclidean distance between the central surface\n\
   and its convex hull.\n\n";

    fprintf(stderr, usage_str, executable);
}

int
main(int argc, char *argv[])
{
    char         *object_file, *sphere_file, *output_surface_file;
    File_formats     format;
    int          n_objects, i;
    object_struct    **objects;
    polygons_struct    *polygons, *sphere;
    double         *depth_values;
    Point        points[MAX_POINTS_PER_POLYGON];

    initialize_argument_processing(argc, argv);

    if (ParseArgv(&argc, argv, argTable, 0) || argc < 3) {
        usage(argv[0]);
        fprintf(stderr, "   %s -help\n\n", argv[0]);
        exit(EXIT_FAILURE);
    }

    if (!get_string_argument(NULL, &object_file) ||
      !get_string_argument(NULL, &sphere_file) ||
      !get_string_argument(NULL, &output_surface_file)) {
        usage(argv[0]);
        exit(EXIT_FAILURE);
    }

    if (input_graphics_any_format(object_file, &format,
                    &n_objects, &objects) != OK)
        exit(EXIT_FAILURE);

    polygons = get_polygons_ptr(objects[0]);

    if (n_objects != 1 || get_object_type(objects[0]) != POLYGONS) {
        printf("File must contain 1 polygons object.\n");
        exit(EXIT_FAILURE);
    }

    if (input_graphics_any_format(sphere_file, &format,
                    &n_objects, &objects) != OK)
        exit(EXIT_FAILURE);

    sphere = get_polygons_ptr(objects[0]);

    depth_values = (double *) malloc(sizeof(double) * polygons->n_points);
    compute_sulcus_depth(polygons, depth_values);
  
    if (use_sqrt) {
        for(i=0; i<polygons->n_points; i++) 
            depth_values[i] = sqrt(depth_values[i]);
    }

    output_values_any_format(output_surface_file, polygons->n_points,
                     depth_values, TYPE_DOUBLE);
    
    delete_object_list(n_objects, objects);

    return(EXIT_SUCCESS);
}
