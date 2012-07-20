/* Christian Gaser - christian.gaser@uni-jena.de
/* Department of Psychiatry
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

void
usage(char *executable)
{
        char *usage_str = "\n\
Usage: %s surface.obj sphere.obj [stop_at] [increase_iterations_by_factor]\n\n\
     Maps a surface to a sphere using the caret inflating approach.\n\n\
     The inflating can be limited using stop_at (default 5), where\n\
       1 - Low smooth\n\
       2 - Inflating\n\
       3 - Very inflating\n\
       4 - High smoothing\n\
       5 - Ellipsoid\n\n\
     Number of iterations can be increased by an additional factor to deal with high resoluted surfaces.";

        fprintf(stderr, usage_str, executable);
}

int
main(int argc, char *argv[])
{
        char             *input_file, *output_file;
        int              n_objects, i, stop_at, increase_iterations_by_factor;
        File_formats     format;
        object_struct    **object_list;
        polygons_struct  *polygons;
        BOOLEAN          enableFingerSmoothing = 1;
        int              fingerSmoothingIters;
        double           surfarea;
    
        initialize_argument_processing(argc, argv);

        if (!get_string_argument(NULL, &input_file) ||
            !get_string_argument(NULL, &output_file)) {
                usage(argv[0]);
                exit(EXIT_FAILURE);
        }

        get_int_argument(5, &stop_at);
        get_int_argument(1, &increase_iterations_by_factor);
    
        if (input_graphics_any_format(input_file, &format, &n_objects,
                                      &object_list) != OK || n_objects != 1 ||
            get_object_type(object_list[0]) != POLYGONS) {
                fprintf(stderr, "Error reading %s.\n", input_file);
                exit(EXIT_FAILURE);
        }

        polygons = get_polygons_ptr(object_list[0]);

        if (euler_characteristic(polygons) != 2) {
                fprintf(stderr, "Euler characteristic of %s must be 2.\n",
                            input_file);
        }
        
        surf_to_sphere(polygons, stop_at, increase_iterations_by_factor);
     
        if(output_graphics_any_format(output_file, format, 1, 
                        object_list) != OK)
                exit(EXIT_FAILURE);
        return(EXIT_SUCCESS);
}
