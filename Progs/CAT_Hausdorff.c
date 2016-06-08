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

        hd = (double *) malloc(sizeof(double) * polygons->n_points);

        if (exact) { /* exact method */
                max_hd = compute_exact_hausdorff(polygons, polygons2, hd);
        } else { /* point-by-point method */
                max_hd = compute_point_hausdorff(polygons, polygons2, hd, 1);
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
