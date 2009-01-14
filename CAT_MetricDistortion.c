/* Rachel Yotter - rachel.yotter@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 *
 */

#include <volume_io/internal_volume_io.h>
#include <volume_io/geometry.h>
#include <bicpl.h>
#include <ParseArgv.h>

#include "CAT_Curvature.h"
#include "CAT_Blur2d.h"
#include "CAT_Surf.h"
#include "CAT_SurfaceIO.h"

#define PINF  1.7976931348623157e+308 /* for doubles */
#define NINF -1.7976931348623157e+308 /* for doubles */

BOOLEAN Energy = 0; /* by default, use the direct metric distortion */

/* the argument table */
ArgvInfo argTable[] = {
  { "-energy", ARGV_CONSTANT, (char *) 1,
    (char *) &Energy,
    "Calculate metric distortion using energy (Fischl)." },
  { NULL, ARGV_END, NULL, NULL, NULL }
};


void
usage(char *executable)
{
        char *usage_str =
"\nUsage: %s [options] object_file sphericalmap_file output_file\n\n\
    Calculate the metric distortion between a brain surface and its\n\
    spherical map.  Results are saved in the text file output_file.\n\n";

        fprintf(stderr, usage_str, executable);
}

int
main(int argc, char** argv)
{
        char               *in_file, *map_file, *out_file;
        object_struct      **objects, **objects2;
        polygons_struct    *polygons, *polygons2;
        int                *n_neighbours, **neighbours;
        int                n_objects;
        File_formats       format;
        int                n, p, nidx;
        double             *metric_dist;
        double             length, length2, radius, total_dist;
        Point              dir;

        if (ParseArgv(&argc, argv, argTable, 0) || (argc < 4)) {
                usage(argv[0]);
                fprintf(stderr, "       %s -help\n\n", argv[0]);
                return(1);
        }

        initialize_argument_processing(argc, argv);
        if (!get_string_argument(NULL, &in_file) ||
            !get_string_argument(NULL, &map_file) ||
            !get_string_argument(NULL, &out_file)) {
                fprintf(stderr, "\nUsage: %s [options] object_file sphericalmap_file output_file\n", argv[0]);
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

        if (input_graphics_any_format(map_file, &format, &n_objects,
                                      &objects2) != OK) {
                printf("Error reading spherical map file %s\n", map_file);
                return(1);
        }

        if (n_objects != 1 || get_object_type(objects2[0]) != POLYGONS) {
                printf("Input file must contain one polygon object.\n");
                return(1);
        }

        polygons = get_polygons_ptr(objects[0]);
        polygons2 = get_polygons_ptr(objects2[0]);

        if (polygons->n_points != polygons2->n_points ||
            polygons->n_items != polygons2->n_items) {
                printf("The input and spherical map meshes do not match.\n");
                return(1);
        }
    
        create_polygon_point_neighbours(polygons, TRUE, &n_neighbours,
                                         &neighbours, NULL, NULL);

        metric_dist = (double *) malloc(sizeof(double) * polygons->n_points);

        /* map the spherical map to a sphere of the same surface area */
        radius = sqrt(get_polygons_surface_area(polygons) / (4.0 * PI));
        for (p = 0; p < polygons->n_points; p++)
                set_vector_length(&polygons2->points[p], radius);

        for (p = 0; p < polygons->n_points; p++) {
                metric_dist[p] = 0;
                if (n_neighbours[p] <= 1)
                        continue; /* skip this point */

                for (n = 0; n < n_neighbours[p]; n++) {
                        nidx = neighbours[p][n];

                        /* length to neighbouring nodes */
                        SUB_POINTS(dir, polygons->points[nidx],
                                        polygons->points[p]);
                        length = MAGNITUDE(dir);
                        SUB_POINTS(dir, polygons2->points[nidx],
                                        polygons2->points[p]);
                        length2 = MAGNITUDE(dir);

                        if (Energy)
                            metric_dist[p] += pow(length2 - length, 2);
                        else
                            metric_dist[p] += fabs(length2 - length) / length;
                }
                if (!Energy)
                        metric_dist[p] /= n_neighbours[p];
        }

        output_values_any_format(out_file, polygons->n_points, metric_dist);

        total_dist = 0;
        for (p = 0; p < polygons->n_points; p++)
                total_dist += metric_dist[p];


        if (Energy)
                total_dist /= 4 * polygons->n_points;
        else
                total_dist /= polygons->n_points;

        printf("Metric distortion = %f\n", total_dist);

        delete_object_list(n_objects, objects);
        delete_object_list(n_objects, objects2);

        free(metric_dist);

        return(0);
}
