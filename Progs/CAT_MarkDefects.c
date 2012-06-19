/* Rachel Yotter - rachel.yotter@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 */

#include <bicpl.h>
#include <ParseArgv.h>

#include "CAT_Map.h"
#include "CAT_Smooth.h"
#include "CAT_SPH.h"
#include "CAT_Octree.h"
#include "CAT_SurfaceIO.h"
#include "CAT_Patch.h"

BOOLEAN dump_patch = FALSE; /* dump patches of defects */

/* minimize the size of the defect */
int
test_minimize(polygons_struct *surface, int *defects, int *polydefects, int *n_neighbours, int **neighbours)
{
        int p, n, euler, n_euler, done, *n_defects, defect = 1;

        for (n = 0; n < surface->n_items; n++) polydefects[n] = 1;
        for (n = 0; n < surface->n_points; n++) defects[n] = 1;

        euler = defect_euler(surface, defects, polydefects, 1,
                             n_neighbours, neighbours);

        /* delete triangles that don't change the euler number */
        done = 0;
        while (!done) {
                done = 1;
                for (p = 0; p < surface->n_points; p++) {
                        if (defects[p] != defect)
                                continue;

                        //if (isedge(surface, defects, polydefects,
                        //           n_neighbours, neighbours, p) == 0) {
                        //        continue;
                        //}
                        //printf("poly %d is an edge\n", p);

                        defects[p] = 0;
                        update_polydefects(surface, defects, polydefects);
                        n_euler = defect_euler(surface, defects,
                                               polydefects, defect,
                                               n_neighbours, neighbours);
                        if (euler != n_euler) {
                                printf("poly %d: %d != %d\n", p, n_euler, euler);
                                defects[p] = defect;
                                update_polydefects(surface, defects, polydefects);
                        } else {
                                printf("deleted point %d\n", p);
                                done = 0;
                        }
                }
        }
        update_polydefects(surface, defects, polydefects);
}

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
"\nUsage: %s [options] surface_file sphere_file output_file\n\n\
    Locate and mark topological errors using a spherical mapping.  Output is a text file.\n\n";

        fprintf(stderr, usage_str, executable);
}


int
main(int argc, char *argv[])
{
        int                  *n_neighbours, **neighbours;
        char                 *surface_file, *sphere_file, *out_file;
        File_formats         format;
        int                  p, n_objects, n_defects, *defects, *polydefects;
        polygons_struct      *surface, *sphere, *patch;
        object_struct        **objects, **sphere_objects, **patch_objects;
        char                 str[80];

        if (ParseArgv(&argc, argv, argTable, 0) || argc != 4) {
                usage(argv[0]);
                fprintf(stderr, "       %s -help\n\n", argv[0]);
                exit(EXIT_FAILURE);
        }

        initialize_argument_processing(argc, argv);

        if (!get_string_argument(NULL, &surface_file) ||
            !get_string_argument(NULL, &sphere_file) ||
            !get_string_argument(NULL, &out_file)) {
                usage(argv[0]);
                exit(EXIT_FAILURE);
        }
     
        if (input_graphics_any_format(surface_file, &format,
                                      &n_objects, &objects) != OK)
                exit(EXIT_FAILURE);

        if (n_objects != 1 || get_object_type(objects[0]) != POLYGONS) {
                printf("Surface file must contain 1 polygons object.\n");
                exit(EXIT_FAILURE);
        }
        surface = get_polygons_ptr(objects[0]);

        if (input_graphics_any_format(sphere_file, &format,
                                      &n_objects, &sphere_objects) != OK)
                exit(EXIT_FAILURE);

        if (n_objects != 1 || get_object_type(sphere_objects[0]) != POLYGONS) {
                printf("Surface file must contain 1 polygons object.\n");
                exit(EXIT_FAILURE);
        }

        sphere = get_polygons_ptr(sphere_objects[0]);

        if (surface->n_items != sphere->n_items) {
                fprintf(stderr,"Surface and sphere must have same size.\n");
                exit(EXIT_FAILURE);
        }

        create_polygon_point_neighbours(surface, TRUE, &n_neighbours,
                                        &neighbours, NULL, NULL);

        defects = (int *) malloc(sizeof(int) * sphere->n_points);
        n_defects = find_topological_defects(surface, sphere, defects,
                                             n_neighbours, neighbours);

        //test_minimize(surface, defects, polydefects, n_neighbours, neighbours);

        printf("%d errors found\n", n_defects);

        output_values_any_format(out_file, sphere->n_points,
                                 defects, TYPE_INTEGER);

        if (dump_patch == TRUE) {
                for (p = 1; p <= n_defects; p++) {
                        sprintf(str, "patch_%d.obj\n", p);
                        patch_objects = extract_patch_points(surface,
                                                             defects, p);
                        patch = get_polygons_ptr(objects[0]);
                        if (output_graphics_any_format(str, ASCII_FORMAT, 1,
                                                      patch_objects) != OK)
                                    exit(EXIT_FAILURE);
                        delete_object_list(1, patch_objects);
                }
        }
        return(EXIT_SUCCESS);

        /* clean up */
        free(defects);

        delete_polygon_point_neighbours(sphere, n_neighbours,
                                        neighbours, NULL, NULL);

        delete_object_list(1, sphere_objects);
    
        return(EXIT_SUCCESS);
}