/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 *
 */

#include <bicpl.h>

#include "CAT_Surf.h"
#include "CAT_SurfaceIO.h"

void
usage(char *executable)
{
        static char *usage_str = "\n\
Usage: %s  surface_file output_values_file [radius]\n\
    Computes surface ratio based on the method of Toro et al. 2008.\n\
    A default radius of 20mm is used.\n\
\n\n";

        fprintf(stderr, usage_str, executable);
}

int
main(int argc, char *argv[])
{
        char                 *object_file, *output_surface_file;
        File_formats         format;
        int                  n_objects;
        object_struct        **objects;
        polygons_struct      *polygons;
        double               *area_values;
        double                 radius;

        initialize_argument_processing(argc, argv);

        if (!get_string_argument(NULL, &object_file) ||
            !get_string_argument(NULL, &output_surface_file)) {
                usage(argv[0]);
                exit(EXIT_FAILURE);
        }

        get_real_argument(20.0, &radius);

        if (input_graphics_any_format(object_file, &format,
                                      &n_objects, &objects) != OK)
                exit(EXIT_FAILURE);

        if (n_objects != 1 || get_object_type(objects[0]) != POLYGONS) {
                printf("File must contain 1 polygons object.\n");
                exit(EXIT_FAILURE);
        }

        polygons = get_polygons_ptr(objects[0]);

        area_values = get_surface_ratio(radius, polygons);
    
        output_values_any_format(output_surface_file, polygons->n_points,
                                 area_values, TYPE_DOUBLE);

        delete_object_list(n_objects, objects);
        free(area_values);

        return(EXIT_SUCCESS);
}
