/* Rachel Yotter - rachel.yotter@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Rachel Yotter, University of Jena.
 */

#include <volume_io/internal_volume_io.h>
#include <bicpl.h>
#include <ParseArgv.h>

#include "CAT_Surf.h"

#define PINF  1.7976931348623157e+308 /* for doubles */
#define NINF -1.7976931348623157e+308 /* for doubles */

BOOLEAN exact = 0; /* 0 - find the closest point, 1 - match point-for-point */

/* the argument table */
ArgvInfo argTable[] = {
  { "-exact", ARGV_CONSTANT, (char *) 1,
    (char *) &exact,
    "Calculate the Hausdorff distance on a point-by-point basis.  Requires that both meshes are of the same brain with the same number of points." },
  { NULL, ARGV_END, NULL, NULL, NULL }
};

/*
 * Go through the neighboring triangles around the closest point and find the
 * closest distance between each triangle and the point p.
 */
double
get_closest_dist(Point *p, Point *p2, int n_neighbours, int *neighbours,
                 polygons_struct *polygons2)
{
        int n;
        Point pts[3], closest;
        double dist, min_dist = PINF;

        pts[0] = *p2;
        for (n = 0; n < n_neighbours - 1; n++) {
                pts[1] = polygons2->points[ neighbours[n] ];
                pts[2] = polygons2->points[ neighbours[n+1] ];
                dist = find_point_polygon_distance_sq(p, 3, pts, &closest);
                if (dist < min_dist)
                        min_dist = dist;
        }
        pts[1] = polygons2->points[ neighbours[n] ];
        pts[2] = polygons2->points[ neighbours[0] ];
        dist = find_point_polygon_distance_sq(p, 3, pts, &closest);
        if (dist < min_dist)
                min_dist = dist;

        return(sqrt(min_dist));
}

/*
 * Calculate the exact Hausdorff distance.  This assumes that the two
 * input meshes are the same size and of the same brain.
 */
double
calc_exact_hausdorff(polygons_struct *p, polygons_struct *p2, double *hd)
{
        int i;
        double max_hd = 0.0, avg_hd = 0.0;

        /* walk through the points */
        for (i = 0; i < p->n_points; i++) {
                hd[i] = distance_between_points(&p->points[i], &p2->points[i]);
                hd[i] = sqrt(hd[i]);

                if (hd[i] > max_hd)
                        max_hd = hd[i];
                avg_hd += hd[i];
        }

        avg_hd /= p->n_points;
        printf("Hausdorff distance: %f\n", max_hd);
        printf("Mean Hausdorff distance: %f\n", avg_hd);

        return(max_hd);
}

/*
 * Calculate the general Hausdorff distance.  The two input meshes do
 * not need to be the same size.
 */
double
calc_hausdorff(polygons_struct *p, polygons_struct *p2, double *hd)
{
        int i, j, n, minj;
        int *n_neighbours, **neighbours;
        int *n_neighbours2, **neighbours2;
        double max_hd, max_revhd;
        double avg_hd, avg_revhd;
        double dist, min_dist;
        int *revpts;
        double *revhd;
        Point closest, pts[3];
        progress_struct progress;

        initialize_progress_report(&progress, FALSE, p->n_points,
                                   "CalcHausdorff");

        create_polygon_point_neighbours(p, TRUE, &n_neighbours,
                                        &neighbours, NULL, NULL);
        create_polygon_point_neighbours(p2, TRUE, &n_neighbours2,
                                        &neighbours2, NULL, NULL);

        revpts = (int *) malloc(sizeof(int) * p2->n_points);
        revhd = (double *) malloc(sizeof(double) * p2->n_points);

        for (j = 0; j < p2->n_points; j++)
                revhd[j] = PINF; /* mark as not found yet */

        /* walk through the points */
        max_hd = 0; avg_hd = 0;
        for (i = 0; i < p->n_points; i++) {
                min_dist = PINF;

                for (j = 0; j < p2->n_points; j++) {
                        dist = sq_distance_between_points(&p->points[i],
                                                          &p2->points[j]);
                        if (dist < min_dist) {
                                min_dist = dist;
                                minj = j;
                        }
                        if (dist < revhd[j]) {
                                revhd[j] = dist; /* save for reverse calc */
                                revpts[j] = i;
                        }
                }
                min_dist = sqrt(min_dist);

                dist = get_closest_dist(&p->points[i], &p2->points[minj],
                                        n_neighbours2[minj], neighbours2[minj],
                                        p2);

                hd[i] = min_dist < dist ? min_dist : dist;
                if (hd[i] > max_hd)
                        max_hd = hd[i];
                avg_hd += hd[i];

                update_progress_report(&progress, i);
        }

        /* calculate the reverse Hausdorff distance */
        max_revhd = 0;
        for (j = 0; j < p2->n_points; j++) {
                revhd[j] = sqrt(revhd[j]);
                dist = get_closest_dist(&p2->points[j], &p->points[revpts[j]],
                                        n_neighbours[revpts[j]],
                                        neighbours[revpts[j]], p);
                if (dist < revhd[j])
                        revhd[j] = dist;
                avg_revhd += revhd[j];
                if (revhd[j] > max_revhd)
                        max_revhd = revhd[j];
        }

        avg_hd /= p->n_points;
        avg_revhd /= p2->n_points;

        printf("Hausdorff distance: %f\n", max_hd);
        printf("Reverse Hausdorff distance: %f\n", max_revhd);
        printf("Mean Hausdorff distance: %f\n", avg_hd);
        printf("Mean reverse Hausdorff distance: %f\n", avg_revhd);

        free(revpts);
        free(revhd);

        return(max_hd);
}

int
main(int argc, char *argv[])
{
        char                 *object_file, *object2_file, *output_file;
        FILE                 *fp;
        File_formats         format;
        int                  n_objects;
        int                  i;
        object_struct        **objects, **objects2;
        polygons_struct      *polygons, *polygons2;
        double               max_hd = 0, *hd;

        /* Call ParseArgv */
        if (ParseArgv(&argc, argv, argTable, 0) || (argc < 3)) {
                fprintf(stderr,"\nUsage: %s [options] object_file object_file2 output_file\n", argv[0]);
                fprintf( stderr,"\nCalculate Hausdorff distance between two surfaces.\n");
                fprintf(stderr, "       %s -help\n\n", argv[0]);
                return(1);
        }

        initialize_argument_processing(argc, argv);

        if (!get_string_argument(NULL, &object_file) ||
            !get_string_argument(NULL, &object2_file) ||
            !get_string_argument(NULL, &output_file)) {
                fprintf(stderr,
                      "Usage: %s  object_file object_file2 output_file\n",
                      argv[0]);
                return(1);
        }

        if (input_graphics_any_format(object_file, &format,
                                      &n_objects, &objects) != OK) {
                return(1);
        }

        if (n_objects != 1 || get_object_type(objects[0]) != POLYGONS) {
                printf("File must contain 1 polygons object.\n");
                return(1);
        }

        if (input_graphics_any_format(object2_file, &format,
                                      &n_objects, &objects2) != OK) {
                return(1);
        }

        if (n_objects != 1 || get_object_type(objects2[0]) != POLYGONS) {
                printf("File must contain 1 polygons object.\n");
                return(1);
        }

        if (open_file(output_file, WRITE_FILE, ASCII_FORMAT, &fp) != OK) {
                return(1);
        }

        polygons = get_polygons_ptr(objects[0]);
        polygons2 = get_polygons_ptr(objects2[0]);

        if (exact && (polygons->n_items != polygons2->n_items ||
                      polygons->n_points != polygons2->n_points)) {
                fprintf(stderr, "Input polygons don't match. Exiting.\n");
                return(1);
        }

        ALLOC(hd, polygons->n_points);

        if (exact) { /* O(n) time */
                max_hd = calc_exact_hausdorff(polygons, polygons2, hd);
        } else { /* O(n*m) time */
                max_hd = calc_hausdorff(polygons, polygons2, hd);
        }

        for (i = 0; i < polygons->n_points; i++) {
                if (output_double(fp, hd[i]) != OK || output_newline(fp) != OK)
                        return(1);
        }

        close_file(fp);

        delete_object_list(n_objects, objects);
        delete_object_list(n_objects, objects2);

        FREE(hd);
        return(0);
}
