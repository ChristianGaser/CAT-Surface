/* Rachel Yotter - rachel.yotter@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id: CAT_SelfIntersections.c 121 2009-07-27 14:42:19Z raytrace $
 *
 */

#include <volume_io/internal_volume_io.h>
#include <volume_io/geometry.h>
#include <bicpl.h>
#include <ParseArgv.h>

#include "CAT_Curvature.h"
#include "CAT_Blur2d.h"
#include "CAT_Octree.h"
#include "CAT_SurfaceIO.h"

#define PINF  1.7976931348623157e+308 /* for doubles */
#define NINF -1.7976931348623157e+308 /* for doubles */


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
        char               *in_file, *out_file;
        object_struct      **objects;
        polygons_struct    *polygons;
        int                *n_neighbours, **neighbours;
        int                *defects, *polydefects;
        int                n_objects;
        File_formats       format;
        int                n_intersects, i;
        progress_struct    progress;
        FILE               *fp;

        initialize_argument_processing(argc, argv);
        if (!get_string_argument(NULL, &in_file) ||
            !get_string_argument(NULL, &out_file)) {
                fprintf(stderr, "\nUsage: %s object_file output_file\n",
                                argv[0]);
                return(1);
        }

        if (input_graphics_any_format(in_file, &format, &n_objects,
                                      &objects) != OK) {
                printf("Error reading input file %s\n", in_file);
                return(1);
        }

        if (n_objects != 1 || get_object_type(objects[0]) != POLYGONS) {
                printf("Input file must contain one polygon object.\n");
                return(1);
        }

        polygons = get_polygons_ptr(objects[0]);

        defects = (int *) malloc(sizeof(int) * polygons->n_points);
        polydefects = (int *) malloc(sizeof(int) * polygons->n_items);

        n_intersects = find_selfintersections(polygons, defects, polydefects);

        printf("All Triangle Intersections:  %d\n", n_intersects);

        /* consolidate intersections */
        create_polygon_point_neighbours(polygons, TRUE, &n_neighbours,
                                        &neighbours, NULL, NULL);
        n_intersects = join_intersections(polygons, defects, n_neighbours,
                                          neighbours);

        if (open_file(out_file, WRITE_FILE, ASCII_FORMAT, &fp) != OK)
                exit(0);
        for (i = 0; i < polygons->n_points; i++)
                fprintf(fp, " %0.1f\n", defects[i]);
        fclose(fp);

        printf("Self Intersections:  %d\n", n_intersects);

        free(defects);
        free(polydefects);

        delete_object_list(n_objects, objects);

        return(0);
}
