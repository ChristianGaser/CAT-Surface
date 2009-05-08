/* Rachel Yotter - rachel.yotter@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id: CAT_DumpCurv.c 85 2009-01-14 23:29:48Z raytrace $
 */

#include <volume_io/internal_volume_io.h>
#include <bicpl.h>

#include "CAT_Blur2d.h"
#include "CAT_Curvature.h"

void
usage(char *executable)
{
        char *usage_str = "\n\
Usage: %s  object_file output_file\n\n\
     Dump the sharpness values.\n";

        fprintf(stderr, usage_str, executable);
}

int
main(int argc, char *argv[])
{
        char                 *object_file, *output_file;
        FILE                 *fp;
        File_formats         format;
        int                  n_objects;
        int                  *n_neighbours, **neighbours;
        object_struct        **objects;
        polygons_struct      *polygons;
        double               *sharpness;

        initialize_argument_processing(argc, argv);

        if (!get_string_argument(NULL, &object_file) ||
            !get_string_argument(NULL, &output_file)) {
                usage(argv[0]);
                return(1);
        }

        if (input_graphics_any_format(object_file, &format, &n_objects,
                                      &objects) != OK)
                return(1);

        if (n_objects != 1 || get_object_type(objects[0]) != POLYGONS) {
                printf("File must contain 1 polygons object.\n");
                return(1);
        }

        polygons = get_polygons_ptr(objects[0]);
        get_all_polygon_point_neighbours(polygons, &n_neighbours, &neighbours);

        ALLOC(sharpness, polygons->n_points);
        compute_local_sharpness(polygons, n_neighbours, neighbours, sharpness);

        output_values_any_format(output_file, polygons->n_points, sharpness);

        delete_object_list(n_objects, objects);
        FREE(sharpness);
    
        return(0);
}