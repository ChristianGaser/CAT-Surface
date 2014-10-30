/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id: CAT_DumpSurfArea.c 261 2012-04-10 10:24:11Z gaser $
 *
 */

#include <bicpl.h>

#include "CAT_Surf.h"
#include "CAT_SurfaceIO.h"
    
void
usage(char *executable)
{
        char *usage_str = "\n\
Usage: %s  surface_file sphere_file output_values_file \n\n\
     Calculate sulcus depth based on the euclidian distance between the surface\n\
     and its convex hull.\n\n";

        fprintf(stderr, usage_str, executable);
}

int
main(int argc, char *argv[])
{
        char                 *object_file, *sphere_file, *output_surface_file;
        File_formats         format;
        int                  n_objects;
        object_struct        **objects;
        polygons_struct      *polygons, *sphere;
        double               *depth_values;
        Point                points[MAX_POINTS_PER_POLYGON];

        initialize_argument_processing(argc, argv);

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
        get_sulcus_depth(polygons, sphere, depth_values);
    
        output_values_any_format(output_surface_file, polygons->n_points,
                                         depth_values, TYPE_DOUBLE);
        
        delete_object_list(n_objects, objects);

        return(EXIT_SUCCESS);
}
