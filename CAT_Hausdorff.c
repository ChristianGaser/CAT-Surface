/* Rachel Yotter - rachel.yotter@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Rachel Yotter, University of Jena.
 * $Id$
 *
 */

#include <volume_io/internal_volume_io.h>
#include <bicpl.h>
#include <ParseArgv.h>

#include "CAT_Surf.h"
#include "CAT_Octree.h"

BOOLEAN exact = 0; /* 0 - find the closest point, 1 - match point-for-point */
BOOLEAN point2point = 0; /* 0 - closest surface, 1 - closest mesh point */

/* the argument table */
ArgvInfo argTable[] = {
  { "-exact", ARGV_CONSTANT, (char *) 1, (char *) &exact,
    "Calculate the Hausdorff distance on a point-by-point basis.  Requires that both meshes are of the same brain with the same number of points." },
  { "-point", ARGV_CONSTANT, (char *) 1, (char *) &point2point,
    "Calculate the Hausdorff distance using only the mesh points.  Faster with a slight overestimation." },
  { NULL, ARGV_END, NULL, NULL, NULL }
};

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
 * Calculate the Hausdorff distance using mesh points only.
 */
double
calc_point_hausdorff(polygons_struct *p, polygons_struct *p2, double *hd)
{
        int i, poly;
        double *revhd;
        double max_hd = 0.0, avg_hd = 0.0;
        double max_revhd = 0.0, avg_revhd = 0.0;
        Point closest;

        create_polygons_bintree(p, ROUND((double) p->n_items * 0.5));
        create_polygons_bintree(p2, ROUND((double) p2->n_items * 0.5));

        /* walk through the points */
        for (i = 0; i < p->n_points; i++) {
                poly = find_closest_polygon_point(&p->points[i], p2, &closest);
                hd[i] = distance_between_points(&p->points[i], &closest);

                if (hd[i] > max_hd)
                        max_hd = hd[i];
                avg_hd += hd[i];
        }

        avg_hd /= p->n_points;
        printf("Hausdorff distance: %f\n", max_hd);
        printf("Mean Hausdorff distance: %f\n", avg_hd);

        /* Calculate the reverse Hausdorff */

        revhd = (double *) malloc(sizeof(double) * p2->n_points);

        for (i = 0; i < p2->n_points; i++) {
                poly = find_closest_polygon_point(&p2->points[i], p, &closest);
                revhd[i] = distance_between_points(&p2->points[i], &closest);

                if (revhd[i] > max_revhd)
                        max_revhd = revhd[i];
                avg_hd += revhd[i];
        }

        avg_revhd /= p2->n_points;
        printf("Reverse Hausdorff distance: %f\n", max_revhd);
        printf("Mean Hausdorff Reverse distance: %f\n", avg_revhd);

        return(max_hd);
}

/*
 * Calculate the general Hausdorff distance.  The two input meshes do
 * not need to be the same size.
 */
double
calc_hausdorff(polygons_struct *p, polygons_struct *p2, double *hd)
{
        int i;
        double max_hd, max_revhd;
        double avg_hd, avg_revhd;
        double *revhd;
        double val;
        struct octree *tree;
        progress_struct progress;

        /* find the forward hausdorff distances */
        max_hd = 0; avg_hd = 0;
        tree = build_octree(p2);
        create_polygons_bintree(p2, ROUND((double) p2->n_items * 0.5));

        initialize_progress_report(&progress, FALSE, p->n_points,
                                   "ForwardHausdorff");

        for (i = 0; i < p->n_points; i++) {
                hausdorff_distance(p->points[i], p2, tree, &hd[i]);

                if (hd[i] > max_hd) {
                        max_hd = hd[i];
                }
                avg_hd += hd[i];

                update_progress_report(&progress, i);
        }
        terminate_progress_report(&progress);
        delete_octree(tree);

        avg_hd /= p->n_points;
        printf("Hausdorff distance: %f\n", max_hd);
        printf("Mean Hausdorff distance: %f\n\n", avg_hd);

        /* calculate the reverse Hausdorff distance */
        revhd = (double *) malloc(sizeof(double) * p2->n_points);
        max_revhd = 0; avg_revhd = 0;

        tree = build_octree(p);
        create_polygons_bintree(p, ROUND((double) p->n_items * 0.5));

        initialize_progress_report(&progress, FALSE, p->n_points,
                                   "ReverseHausdorff");

        for (i = 0; i < p2->n_points; i++) {
                hausdorff_distance(p2->points[i], p, tree, &revhd[i]);

                if (revhd[i] > max_revhd)
                        max_revhd = revhd[i];
                avg_revhd += revhd[i];

                update_progress_report(&progress, i);
        }
        terminate_progress_report(&progress);

        avg_revhd /= p2->n_points;
        printf("Reverse Hausdorff distance: %f\n", max_revhd);
        printf("Mean reverse Hausdorff distance: %f\n\n", avg_revhd);

        printf("Mean(Fwd + Rev) Hausdorff distance: %f\n",
               (max_revhd + max_hd)/2);
        printf("Mean(Fwd + Rev) Mean Hausdorff distance: %f\n",
               (avg_revhd + avg_hd)/2);

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
        } else if (point2point) { /* O(n*log(m)) time */
                max_hd = calc_point_hausdorff(polygons, polygons2, hd);
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
