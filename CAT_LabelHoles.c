/* Rachel Yotter - rachel.yotter@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id: CAT_LabelHoles.c 89 2009-01-27 14:43:59Z raytrace $
 */

#include <volume_io/internal_volume_io.h>
#include <bicpl.h>
#include <bicpl/deform.h>
#include <ParseArgv.h>

#include "CAT_SheetIO.h"
#include "CAT_Map2d.h"
#include "CAT_Blur2d.h"
#include "CAT_Intersect.h"
#include "CAT_SurfaceIO.h"
#include "CAT_Defect.h"


void
usage(char *executable)
{
        static char *usage_str = "\n\
Usage: %s surface.obj sphere.obj volume.mnc output.txt\n\
Find and label the holes (=1), handles (=2), and large errors (=3) of a surface.\n\n\n";

       fprintf(stderr, usage_str, executable);
}

int
main(int argc, char *argv[])
{
        char                 *surface_file, *sphere_file;
        char                 *volume_file, *output_file;
        File_formats         format;
        int                  n_objects;
        Volume               volume;
        polygons_struct      *surface, *sphere;
        object_struct        **objects, **sphere_objects;
        int                  *n_neighbours, **neighbours;
        int                  p, n_defects, *defects, *polydefects, *holes;
        double               t1_threshold = 0;

        /* Call ParseArgv */
        if (argc != 5) {
                usage(argv[0]);
                fprintf(stderr, "       %s -help\n\n", argv[0]);
                exit(EXIT_FAILURE);
        }

        initialize_argument_processing(argc, argv);

        if (!get_string_argument(NULL, &surface_file) ||
            !get_string_argument(NULL, &sphere_file) ||
            !get_string_argument(NULL, &volume_file) ||
            !get_string_argument(NULL, &output_file)) {
                usage(argv[0]);
                exit(EXIT_FAILURE);
        }
     
        if (input_volume_all(volume_file, 3, File_order_dimension_names,
                             NC_UNSPECIFIED, FALSE, 0.0, 0.0, TRUE,
                             &volume, NULL) == ERROR) {
                fprintf(stderr, "Error opening T1 file: %s\n", volume_file);
                exit(EXIT_FAILURE);
        }

        if (input_graphics_any_format(surface_file, &format,
                                      &n_objects, &objects) != OK)
                exit(EXIT_FAILURE);

        /* check that the surface file contains a polyhedron */
        if (n_objects != 1 || get_object_type(objects[0]) != POLYGONS) {
                fprintf(stderr,"Surface file must contain 1 polygons object.\n");
                exit(EXIT_FAILURE);
        }

        surface = get_polygons_ptr(objects[0]);
    
        if (input_graphics_any_format(sphere_file, &format,
                                      &n_objects, &sphere_objects) != OK)
                exit(EXIT_FAILURE);

        /* check that the surface file contains a polyhedron */
        if (n_objects != 1 || get_object_type(sphere_objects[0]) != POLYGONS) {
                fprintf(stderr,"Sphere file must contain 1 polygons object.\n");
                exit(EXIT_FAILURE);
        }

        sphere = get_polygons_ptr(sphere_objects[0]);

        /* check that surface and sphere are same size */
        if (surface->n_items != sphere->n_items) {
                fprintf(stderr,"Surface and sphere must have same size.\n");
                exit(EXIT_FAILURE);
        }
    

        create_polygon_point_neighbours(sphere, TRUE, &n_neighbours,
                                        &neighbours, NULL, NULL);


        /* find defects in original uncorrected surface */
        defects = (int *) malloc(sizeof(int) * sphere->n_points);
        polydefects = (int *) malloc(sizeof(int) * sphere->n_items);
        n_defects = find_topological_defects(surface, sphere, defects,
                                             polydefects, n_neighbours,
                                             neighbours);

        holes = (int *) malloc(sizeof(int) * sphere->n_points);

        if (n_defects > 0) {
                t1_threshold = get_holes_handles(surface, sphere, defects,
                                                 n_defects, holes, volume,
                                                 n_neighbours, neighbours);
                printf("t1 threshold = %f\n", t1_threshold);
        }

        if (output_values_any_format(output_file, surface->n_points,
                                     holes, TYPE_INTEGER) != OK) {
                exit(EXIT_FAILURE);
        }

        delete_polygon_point_neighbours(sphere, n_neighbours,
                                        neighbours, NULL, NULL);
        free(defects);
        free(holes);

        delete_object_list(1, objects);
        delete_object_list(1, sphere_objects);
    
        return(EXIT_SUCCESS);    
}
