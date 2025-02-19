/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * most of the code is modified from
 * caret/caret_brain_set/BrainModelSurface.cxx.
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 *
 */

#include <bicpl.h>
#include <float.h>

#include "CAT_Surf.h"
#include "CAT_SurfaceIO.h"
#include "CAT_ConvexHull.h"

void
usage(char *executable)
{
    char *usage_str = "\n\
Usage: %s surface_file sphere_file output_surface_file \n\n\
   Extracts a convex hull of the surface.\n";

    fprintf(stderr, usage_str, executable);
}

int
main(int argc, char *argv[])
{
    char       *input_file, *sphere_file, *output_surface_file;
    int        n_objects, i;
    File_formats   format;
    object_struct  **object_list;
    polygons_struct  *polygons, *sphere, *convex;
  
    initialize_argument_processing(argc, argv);

    if (!get_string_argument(NULL, &input_file) ||
      !get_string_argument(NULL, &sphere_file) ||
      !get_string_argument(NULL, &output_surface_file)) {
        usage(argv[0]);
        exit(EXIT_FAILURE);
    }

    if (input_graphics_any_format(input_file, &format, &n_objects,
                    &object_list) != OK || n_objects != 1 ||
      get_object_type(object_list[0]) != POLYGONS) {
        fprintf(stderr, "Error reading %s.\n", input_file);
        exit(EXIT_FAILURE);
    }
    polygons = get_polygons_ptr(object_list[0]);
    
    if (input_graphics_any_format(sphere_file, &format, &n_objects,
                    &object_list) != OK || n_objects != 1 ||
      get_object_type(object_list[0]) != POLYGONS) {
        fprintf(stderr, "Error reading %s.\n", input_file);
        exit(EXIT_FAILURE);
    }
    sphere = get_polygons_ptr(object_list[0]);

    /* get convex hull */
    object_list = surface_get_convex_hull(polygons, sphere);
    convex = get_polygons_ptr(object_list[0]);
   
    if(output_graphics_any_format(output_surface_file, format, 1, 
            object_list, NULL) != OK)
        exit(EXIT_FAILURE);
    return(EXIT_SUCCESS);
}
