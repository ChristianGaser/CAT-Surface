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
#include "CAT_Surf.h"
#include "CAT_SurfaceIO.h"

#define PINF  1.7976931348623157e+308 /* for doubles */
#define NINF -1.7976931348623157e+308 /* for doubles */

BOOLEAN Energy = 0; /* by default, use the direct metric distortion */
BOOLEAN Global = 0; /* by default, calculate the local metric distortion */
int NumPts = 0; /* by default, use all points */

/* the argument table */
ArgvInfo argTable[] = {
  { "-energy", ARGV_CONSTANT, (char *) 1,
    (char *) &Energy,
    "Calculate metric distortion using energy (Fischl)." },
  { "-global", ARGV_CONSTANT, (char *) 1,
    (char *) &Global,
    "Calculate global metric distortion (Dijkstra algorithm, slow)." },
  { "-random", ARGV_INT, (char *) 1,
    (char *) &NumPts,
    "Use a random subset of the points [0=all]." },
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

void
dijkstra_dist(polygons_struct *polygons, int *n_neighbours, int **neighbours,
              int p, unsigned char *visited, double *dist)
{
        int i, n, nidx, cur_p;
        Point dir;
        double mindist = PINF, length;
        int done = 0;
        int *set, setsize, setidx;

        visited[p] = 1;
        dist[p] = 0;

        set = (int *) malloc(sizeof(int) * polygons->n_points);
        setsize = 0;

        /* get the distances to the 1-neighbors */
        for (n = 0; n < n_neighbours[p]; n++) {
                nidx = neighbours[p][n];

                /* add it to the working set */
                set[setsize] = nidx;
                setsize++;

                SUB_POINTS(dir, polygons->points[nidx],
                                polygons->points[p]);
                dist[nidx] = MAGNITUDE(dir);

                if (dist[nidx] < mindist) {
                        mindist = dist[nidx];
                        cur_p = nidx;
                        setidx = setsize - 1;
                }
        }


        /* do the closest 1-neighbor */
        visited[cur_p] = 1;
        set[setidx] = set[setsize - 1]; /* remove it from working set */
        setsize--;
        for (n = 0; n < n_neighbours[cur_p]; n++) {
                nidx = neighbours[cur_p][n];

                if (visited[nidx] == 0) {
                        SUB_POINTS(dir, polygons->points[nidx],
                                        polygons->points[cur_p]);
                        length = MAGNITUDE(dir) + dist[cur_p];
                        if (dist[nidx] == PINF) {
                                /* add it to the working set */
                                set[setsize] = nidx;
                                setsize++;
                        }

                        if (length < dist[nidx])
                                dist[nidx] = length;
                }
        }

        /* churn through the rest of them */
        while (!done) {
                mindist = PINF;
                for (i = 0; i < setsize; i++) {
                        if (visited[set[i]] == 0 && dist[set[i]] < mindist) {
                                mindist = dist[set[i]];
                                cur_p = set[i];
                                setidx = i;
                        }
                }
                if (mindist == PINF) { /* check the rest of them */
                        for (i = 0; i < polygons->n_points; i++) {
                                if (visited[i] == 0 && dist[i] < mindist) {
                                        mindist = dist[i];
                                        cur_p = i;
                                }
                        }
                        if (mindist == PINF) { /* really really done */
                                done = 1;
                                break;
                        }
                } else { /* remove it from working set */
                        set[setidx] = set[setsize - 1];
                        setsize--;
                }

                visited[cur_p] = 1;
                for (n = 0; n < n_neighbours[cur_p]; n++) {
                        nidx = neighbours[cur_p][n];

                        if (visited[nidx] == 0) {
                                SUB_POINTS(dir, polygons->points[nidx],
                                                polygons->points[cur_p]);
                                length = MAGNITUDE(dir) + dist[cur_p];

                                if (dist[nidx] == PINF) {
                                        /* add it to the working set */
                                        set[setsize] = nidx;
                                        setsize++;
                                        dist[nidx] = length;
                                } else if (length < dist[nidx]) {
                                        dist[nidx] = length;
                                }
                        }
                }
        }
        free(set);
}

int
main(int argc, char** argv)
{
        char               *in_file, *map_file, *out_file;
        object_struct      **objects, **objects2;
        polygons_struct    *polygons, *polygons2;
        int                *n_neighbours, **neighbours;
        unsigned char      *visited, *pointset;
        double             *geo_dist, *geo_dist2;
        int                n_objects;
        File_formats       format;
        int                n, p, nidx, count;
        double             *metric_dist;
        double             length, length2, radius, total_dist;
        Point              dir;
        progress_struct    progress;


        if (ParseArgv(&argc, argv, argTable, 0) || (argc < 4)) {
                usage(argv[0]);
                fprintf(stderr, "       %s -help\n\n", argv[0]);
                exit(EXIT_FAILURE);
        }

        initialize_argument_processing(argc, argv);
        if (!get_string_argument(NULL, &in_file) ||
            !get_string_argument(NULL, &map_file) ||
            !get_string_argument(NULL, &out_file)) {
                fprintf(stderr, "\nUsage: %s [options] object_file sphericalmap_file output_file\n", argv[0]);
                exit(EXIT_FAILURE);
        }

        if (input_graphics_any_format(in_file, &format, &n_objects,
                                      &objects) != OK) {
                printf("Error reading input file %s\n", in_file);
                exit(EXIT_FAILURE);
        }

        if (n_objects != 1 || get_object_type(objects[0]) != POLYGONS) {
                printf("Input file must contain one polygon object.\n");
                exit(EXIT_FAILURE);
        }

        if (input_graphics_any_format(map_file, &format, &n_objects,
                                      &objects2) != OK) {
                printf("Error reading spherical map file %s\n", map_file);
                exit(EXIT_FAILURE);
        }

        if (n_objects != 1 || get_object_type(objects2[0]) != POLYGONS) {
                printf("Input file must contain one polygon object.\n");
                exit(EXIT_FAILURE);
        }

        polygons = get_polygons_ptr(objects[0]);
        polygons2 = get_polygons_ptr(objects2[0]);

        if (polygons->n_points != polygons2->n_points ||
            polygons->n_items != polygons2->n_items) {
                printf("The input and spherical map meshes do not match.\n");
                exit(EXIT_FAILURE);
        }
    
        create_polygon_point_neighbours(polygons, TRUE, &n_neighbours,
                                         &neighbours, NULL, NULL);

        metric_dist = (double *) malloc(sizeof(double) * polygons->n_points);
        memset(metric_dist, 0, sizeof(double) * polygons->n_points);

        if (Global) {
                visited = (unsigned char *) malloc(sizeof(unsigned char) *
                                                   polygons->n_points);
                geo_dist = (double *) malloc(sizeof(double) *
                                             polygons->n_points);
                geo_dist2 = (double *) malloc(sizeof(double) *
                                             polygons->n_points);
        }

        /* map the spherical map to a sphere of the same surface area */
        radius = sqrt(get_polygons_surface_area(polygons) / (4.0 * PI));
        for (p = 0; p < polygons->n_points; p++)
                set_vector_length(&polygons2->points[p], radius);

        /* create a subset of points to visit if -random flag is used */
        pointset = (unsigned char *) malloc(sizeof(unsigned char) *
                                            polygons->n_points);
        if (NumPts > 0) {
                srand(time(NULL));
                memset(pointset, 0, sizeof(unsigned char) * polygons->n_points);
                for (p = 0; p < NumPts; p++) {
                        n = round(rand() % (polygons->n_points - 1));
                        while (pointset[n] == 1)
                                n = round(rand() % (polygons->n_points - 1));
                        pointset[n] = 1;
                }
        } else {
                NumPts = polygons->n_points;
                memset(pointset, 1, sizeof(char) * polygons->n_points);
        }

        if (Global)
                initialize_progress_report(&progress, FALSE, polygons->n_points,
                                           "Dijkstra Distance");

        for (p = 0; p < polygons->n_points; p++) {
                metric_dist[p] = 0;
                if (pointset[p] == 0 || n_neighbours[p] <= 1)
                        continue; /* skip this point */

                if (Global) {
                        for (n = 0; n < polygons->n_points; n++) {
                                geo_dist[n] = geo_dist2[n] = PINF;
                        }

                        geo_dist[p] = geo_dist2[p] = 0;

                        memset(visited, 0,
                               sizeof(unsigned char) * polygons->n_points);
                        dijkstra_dist(polygons, n_neighbours, neighbours,
                                      p, visited, geo_dist);
                        memset(visited, 0,
                               sizeof(unsigned char) * polygons->n_points);

                        dijkstra_dist(polygons2, n_neighbours, neighbours,
                                      p, visited, geo_dist2);

                        count = 0;
                        for (n = 0; n < polygons->n_points; n++) {
                                if (geo_dist[n] == 0) continue;
                                if (Energy) {
                                        metric_dist[p] += pow(geo_dist2[n] -
                                                              geo_dist[n], 2);
                                } else {
                                        metric_dist[p] += fabs(geo_dist2[n] -
                                                               geo_dist[n]) /
                                                          geo_dist[n];
                                        count++;
                                }
                        }
                        if (!Energy && count > 0)
                                metric_dist[p] /= count;
                        update_progress_report(&progress, p);
                } else {
                        for (n = 0; n < n_neighbours[p]; n++) {
                                nidx = neighbours[p][n];

                                /* length to neighbouring nodes */
                                SUB_POINTS(dir, polygons->points[nidx],
                                                polygons->points[p]);
                                length = MAGNITUDE(dir);
                                SUB_POINTS(dir, polygons2->points[nidx],
                                                polygons2->points[p]);
                                length2 = MAGNITUDE(dir);
                                 if (Energy) {
                                         metric_dist[p] += pow(length2 -
                                                               length, 2);
                                 } else {
                                         metric_dist[p] += fabs(length2 -
                                                                length) /
                                                           length;
                                 }
                        }
                        if (!Energy) {
                                metric_dist[p] /= n_neighbours[p];
                        }
                }

        }

        if (Global)
                terminate_progress_report(&progress);

        output_values_any_format(out_file, polygons->n_points,
                                 metric_dist, TYPE_DOUBLE);

        total_dist = 0;
        for (p = 0; p < polygons->n_points; p++)
                total_dist += metric_dist[p];


        if (Energy)
                total_dist /= 4 * NumPts;
        else
                total_dist /= NumPts;

        printf("Metric distortion = %f\n", total_dist);

        delete_object_list(n_objects, objects);
        delete_object_list(n_objects, objects2);

        free(metric_dist);
        free(pointset);

        if (Global) {
                free(visited);
                free(geo_dist);
                free(geo_dist2);
        }

        return(EXIT_SUCCESS);
}
