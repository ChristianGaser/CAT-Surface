/* Rachel Yotter - rachel.yotter@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 *
 */

#include <bicpl.h>
#include <ParseArgv.h>

#include "CAT_Curvature.h"
#include "CAT_Smooth.h"
#include "CAT_Octree.h"
#include "CAT_SurfaceIO.h"
#include "CAT_Patch.h"
#include "CAT_Intersect.h"

BOOLEAN dump_patch = FALSE; /* dump patches of defects */

/* the argument table */
ArgvInfo argTable[] = {
  { "-patch", ARGV_CONSTANT, (char *) TRUE,
    (char *) &dump_patch,
    "Dump defect patches." },
  { NULL, ARGV_END, NULL, NULL, NULL }
};


void
usage(char *executable)
{
        char *usage_str =
"\nUsage: %s [options] object_file output_file\n\n\
    Find the number of self-intersections and mark the intersecting\n\
    triangles in the output_file.\n\n";

        fprintf(stderr, usage_str, executable);
}


int
main(int argc, char** argv)
{
        char               *surface_file, *out_file;
        object_struct      **objects, **patch_objects;
        polygons_struct    *polygons, *patch;
        int                *n_neighbours, **neighbours;
        int                *defects, *polydefects;
        File_formats       format;
        int                n_objects, n_intersects, i;
        char               str[80];
        progress_struct    progress;
        FILE               *fp;

        if (ParseArgv(&argc, argv, argTable, 0) || argc != 3) {
                usage(argv[0]);
                fprintf(stderr, "       %s -help\n\n", argv[0]);
                exit(EXIT_FAILURE);
        }

        initialize_argument_processing(argc, argv);
        if (!get_string_argument(NULL, &surface_file) ||
            !get_string_argument(NULL, &out_file)) {
                fprintf(stderr, "\nUsage: %s object_file output_file\n",
                                argv[0]);
                exit(EXIT_FAILURE);
        }

        if (input_graphics_any_format(surface_file, &format, &n_objects,
                                      &objects) != OK) {
                printf("Error reading input file %s\n", surface_file);
                exit(EXIT_FAILURE);
        }

        if (n_objects != 1 || get_object_type(objects[0]) != POLYGONS) {
                printf("Input file must contain one polygon object.\n");
                exit(EXIT_FAILURE);
        }

        polygons = get_polygons_ptr(objects[0]);

        defects = (int *) malloc(sizeof(int) * polygons->n_points);
        polydefects = (int *) malloc(sizeof(int) * polygons->n_items);

        create_polygon_point_neighbours(polygons, TRUE, &n_neighbours,
                                        &neighbours, NULL, NULL);

        n_intersects = find_selfintersections(polygons, defects, polydefects);

        printf("All Triangle Intersections:  %d\n", n_intersects);

        /* consolidate intersections */
        n_intersects = join_intersections(polygons, defects, polydefects,
                                          n_neighbours, neighbours);
        printf("Self Intersections:  %d\n", n_intersects);

        output_values_any_format(out_file, polygons->n_points,
                                 defects, TYPE_INTEGER);

        if (dump_patch == TRUE) {
                for (i = 1; i <= n_intersects; i++) {
                        sprintf(str, "patch_%d.obj\n", i);
                        patch_objects = extract_patch_points(polygons,
                                                             defects, i);
                        patch = get_polygons_ptr(objects[0]);
                        if (output_graphics_any_format(str, ASCII_FORMAT, 1,
                                                       patch_objects) != OK)
                                    exit(EXIT_FAILURE);
                        delete_object_list(1, patch_objects);
                }
        }

        free(defects);
        free(polydefects);

        delete_polygon_point_neighbours(polygons, n_neighbours,
                                        neighbours, NULL, NULL);

        delete_object_list(n_objects, objects);

        return(EXIT_SUCCESS);
}