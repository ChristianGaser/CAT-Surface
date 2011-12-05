/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 *
 */

#include <volume_io/internal_volume_io.h>
#include <bicpl.h>

#include "CAT_Blur2d.h"
#include "CAT_Surf.h"
#include "CAT_SurfaceIO.h"
#include "CAT_Metrics.h"
    
int
main(int argc, char *argv[])
{
        char                 *object_file, *output_file;
        File_formats         format;
        int                  poly, n_objects, **neighbours, *n_neighbours;
        object_struct        **objects;
        polygons_struct      *polygons;
        double               *convexity;

        initialize_argument_processing(argc, argv);

        if (!get_string_argument(NULL, &object_file)) {
                fprintf(stderr, "Usage: %s  object_file [output_file]\n",
                        argv[0]);
                exit(EXIT_FAILURE);
        }

        if (input_graphics_any_format(object_file, &format,
                                      &n_objects, &objects) != OK)
                exit(EXIT_FAILURE);

        if (get_string_argument(NULL, &output_file) != OK)
                exit(EXIT_FAILURE);

        if (n_objects != 1 || get_object_type(objects[0]) != POLYGONS) {
                printf("File must contain 1 polygons object.\n");
                exit(EXIT_FAILURE);
        }

        polygons = get_polygons_ptr(objects[0]);

        create_polygon_point_neighbours(polygons, TRUE, &n_neighbours,
                                        &neighbours, NULL, NULL);

        ALLOC(convexity, polygons->n_points);
        calc_convexity(polygons, n_neighbours, neighbours, 1.0f, convexity);
    
        output_values_any_format(output_file, polygons->n_points, convexity,
                                 TYPE_DOUBLE);
        
        delete_object_list(n_objects, objects);

        return(EXIT_SUCCESS);
}
