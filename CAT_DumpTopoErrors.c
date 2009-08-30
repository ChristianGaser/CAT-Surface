/* Rachel Yotter - rachel.yotter@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id: CAT_DumpTopoErrors.c 89 2009-01-27 14:43:59Z raytrace $
 */

#include <volume_io/internal_volume_io.h>
#include <bicpl.h>
#include <ParseArgv.h>

#include "CAT_SheetIO.h"
#include "CAT_Map2d.h"
#include "CAT_Blur2d.h"
#include "CAT_SPH.h"
#include "CAT_Octree.h"


void
usage(char *executable)
{
        char *usage_str =
"\nUsage: %s [options] sphere_file output_file\n\n\
    Locate and mark topological errors using a spherical mapping.  Output is a text file.\n\n";

        fprintf(stderr, usage_str, executable);
}

int
main(int argc, char *argv[])
{
        int                  *n_neighbours, **neighbours;
        char                 *sphere_file, *out_file;
        File_formats         format;
        int                  p, n_defects, *defects, n_objects;
        polygons_struct      *sphere;
        object_struct        **sphere_objects;
        FILE                 *fp;

        initialize_argument_processing(argc, argv);

        if (!get_string_argument(NULL, &sphere_file) ||
            !get_string_argument(NULL, &out_file)) {
                usage(argv[0]);
                exit(EXIT_FAILURE);
        }
     
        if (input_graphics_any_format(sphere_file, &format,
                                      &n_objects, &sphere_objects) != OK)
                exit(EXIT_FAILURE);

        if (n_objects != 1 || get_object_type(sphere_objects[0]) != POLYGONS) {
                printf("Surface file must contain 1 polygons object.\n");
                exit(EXIT_FAILURE);
        }
        /* get a pointer to the surface */
        sphere = get_polygons_ptr(sphere_objects[0]);

        create_polygon_point_neighbours(sphere, TRUE, &n_neighbours,
                                        &neighbours, NULL, NULL);
    
        defects = (int *) malloc(sizeof(int) * sphere->n_points);
        n_defects = find_topological_defects(sphere, defects,
                                             n_neighbours, neighbours);
        expand_defects(sphere, defects, 0, 2, n_neighbours, neighbours);

        printf("%d errors found\n", n_defects);

        if (open_file(out_file, WRITE_FILE, ASCII_FORMAT, &fp) != OK)
                exit(0);

        for (p = 0; p < sphere->n_points; p++)
                fprintf(fp, " %d.0\n", defects[p]);
        fclose(fp);

        /* clean up */
        free(defects);
        delete_polygon_point_neighbours(sphere, n_neighbours,
                                        neighbours, NULL, NULL);
        delete_object_list(1, sphere_objects);
    
        return(EXIT_SUCCESS);
}
