/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 */

#include <volume_io/internal_volume_io.h>
#include <bicpl.h>

#include "CAT_Surf.h"

void
usage(char *executable)
{
        static char *usage_str = "\n\
Usage: %s  surface.obj output_values.txt [radius]\n\
    Computes surface ratio based on the method of Toro et al. 2008.\n\
    A default radius of 20mm is used.\n\
\n\n";

        fprintf(stderr, usage_str, executable);
}

int
main(int argc, char *argv[])
{
        char                 *object_file, *output_file;
        File_formats         format;
        int                  n_objects, ptidx;
        object_struct        **objects;
        polygons_struct      *polygons;
        double               *area_values;
        Real                 radius;

        initialize_argument_processing(argc, argv);

        if (!get_string_argument(NULL, &object_file) ||
            !get_string_argument(NULL, &output_file)) {
                usage(argv[0]);
                return(1);
        }

        get_real_argument(20.0, &radius);

        if (input_graphics_any_format(object_file, &format,
                                      &n_objects, &objects) != OK)
                return(1);

        if (n_objects != 1 || get_object_type(objects[0]) != POLYGONS) {
                printf("File must contain 1 polygons object.\n");
                return(1);
        }

        polygons = get_polygons_ptr(objects[0]);

        ALLOC(area_values, polygons->n_points);
        area_values = get_surface_ratio(radius, polygons);
    
        output_values_any_format(output_file, polygons->n_points, area_values);

        delete_object_list(n_objects, objects);

        return(0);
}
