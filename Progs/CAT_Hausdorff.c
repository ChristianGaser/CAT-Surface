/* Rachel Yotter - rachel.yotter@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Rachel Yotter, University of Jena.
 * $Id$
 *
 */

#include <bicpl.h>
#include <ParseArgv.h>

#include "CAT_Surf.h"
#include "CAT_SurfaceIO.h"

BOOLEAN exact = 0; /* 0 - find the closest point, 1 - match point-for-point */

/* the argument table */
ArgvInfo argTable[] = {
  { "-exact", ARGV_CONSTANT, (char *) 1, (char *) &exact,
    "Calculate the Hausdorff distance on a point-by-point basis.  Requires that both meshes are of the same brain with the same number of points." },
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
        printf("Mean distance error: %f\n", avg_hd);

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
        printf("Mean distance error: %f\n", avg_hd);

        /* Calculate the reverse Hausdorff */

        revhd = (double *) malloc(sizeof(double) * p2->n_points);

        for (i = 0; i < p2->n_points; i++) {
                poly = find_closest_polygon_point(&p2->points[i], p, &closest);
                revhd[i] = distance_between_points(&p2->points[i], &closest);

                if (revhd[i] > max_revhd)
                        max_revhd = revhd[i];
                avg_revhd += revhd[i];
        }

        avg_revhd /= p2->n_points;
        printf("Reverse Hausdorff distance: %f\n", max_revhd);
        printf("Mean reverse distance error: %f\n", avg_revhd);

        return(max_hd);
}

int
main(int argc, char *argv[])
{
        char                 *object_file, *object2_file, *output_surface_file;
        FILE                 *fp;
        File_formats         format;
        int                  n_objects;
        int                  i;
        object_struct        **objects, **objects2;
        polygons_struct      *polygons, *polygons2;
        double               max_hd = 0, *hd;

        /* Call ParseArgv */
        if (ParseArgv(&argc, argv, argTable, 0) || (argc < 3)) {
                fprintf(stderr,"\nUsage: %s [options] surface_file surface_file2 output_values_file\n", argv[0]);
                fprintf( stderr,"\nCalculate Hausdorff distance between two surfaces on a point-by-point basis.\n");
                fprintf(stderr, "       %s -help\n\n", argv[0]);
                exit(EXIT_FAILURE);
        }

        initialize_argument_processing(argc, argv);

        if (!get_string_argument(NULL, &object_file) ||
            !get_string_argument(NULL, &object2_file) ||
            !get_string_argument(NULL, &output_surface_file)) {
                fprintf(stderr,
                      "Usage: %s  object_file object_file2 output_surface_file\n",
                      argv[0]);
                exit(EXIT_FAILURE);
        }

        if (input_graphics_any_format(object_file, &format,
                                      &n_objects, &objects) != OK) {
                exit(EXIT_FAILURE);
        }

        if (n_objects != 1 || get_object_type(objects[0]) != POLYGONS) {
                printf("File must contain 1 polygons object.\n");
                exit(EXIT_FAILURE);
        }

        if (input_graphics_any_format(object2_file, &format,
                                      &n_objects, &objects2) != OK) {
                exit(EXIT_FAILURE);
        }

        if (n_objects != 1 || get_object_type(objects2[0]) != POLYGONS) {
                printf("File must contain 1 polygons object.\n");
                exit(EXIT_FAILURE);
        }

        polygons = get_polygons_ptr(objects[0]);
        polygons2 = get_polygons_ptr(objects2[0]);

        if (exact && (polygons->n_items != polygons2->n_items ||
                      polygons->n_points != polygons2->n_points)) {
                fprintf(stderr, "Input polygons don't match. Exiting.\n");
                exit(EXIT_FAILURE);
        }

        ALLOC(hd, polygons->n_points);

        if (exact) { /* exact method */
                max_hd = calc_exact_hausdorff(polygons, polygons2, hd);
        } else { /* point-by-point method */
                max_hd = calc_point_hausdorff(polygons, polygons2, hd);
        }

        if (output_values_any_format(output_surface_file, polygons->n_points,
                                     hd, TYPE_DOUBLE) != OK) {
                exit(EXIT_FAILURE);
        }

        delete_object_list(n_objects, objects);
        delete_object_list(n_objects, objects2);

        FREE(hd);
        return(EXIT_SUCCESS);
}
