/* Rachel Yotter - rachel.yotter@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 */

#include "CAT_NiftiIO.h"
#include "CAT_SurfaceIO.h"
#include "CAT_DeformPolygons.h"

void
usage(char *executable)
{
        static char *usage_str = "\n\
Usage: %s surface_file target_surface_file output_surface_file\n\
Roughly deform a surface onto the border of another (similar) surface.\n\n\n";

       fprintf(stderr, usage_str, executable);
}


int
main(int argc, char *argv[])
{
        char                 *surface_file, *target_surface_file, *output_surface_file;
        File_formats         format;
        int                  n_objects;
        polygons_struct      *surface;
        object_struct        **surf_objects, **target_objects;

        initialize_argument_processing(argc, argv);

        if (!get_string_argument(NULL, &surface_file) ||
            !get_string_argument(NULL, &target_surface_file) ||
            !get_string_argument(NULL, &output_surface_file)) {
                usage(argv[0]);
                exit(EXIT_FAILURE);
        }

        if (input_graphics_any_format(surface_file, &format,
                                      &n_objects, &surf_objects) != OK)
                exit(EXIT_FAILURE);

        /* check that the surface file contains a polyhedron */
        if (n_objects != 1 || get_object_type(surf_objects[0]) != POLYGONS) {
                fprintf(stderr,"Surface file must contain 1 polygons object.\n");
                exit(EXIT_FAILURE);
        }
        /* get a pointer to the surface */
        surface = get_polygons_ptr(surf_objects[0]);
    
        if (input_graphics_any_format(target_surface_file, &format,
                                      &n_objects, &target_objects) != OK)
                exit(EXIT_FAILURE);

        /* check that the surface file contains a polyhedron */
        if (n_objects != 1 || get_object_type(target_objects[0]) != POLYGONS) {
                fprintf(stderr,"Surface file must contain 1 polygons object.\n");
                exit(EXIT_FAILURE);
        }

        deform_surf2object(surface, target_objects[0]);

        if (output_graphics_any_format(output_surface_file, ASCII_FORMAT, 1,
                                       surf_objects, NULL) != OK)
                exit(EXIT_FAILURE);

        /* clean up */

        delete_object_list(1, surf_objects);
        delete_object_list(1, target_objects);

        return(EXIT_SUCCESS);    
}
