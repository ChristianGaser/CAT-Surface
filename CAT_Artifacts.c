/* Rachel Yotter - rachel.yotter@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 */

#include <volume_io/internal_volume_io.h>
#include <bicpl.h>
#include <ParseArgv.h>

#include "CAT_Map2d.h"
#include "CAT_Blur2d.h"
#include "CAT_SPH.h"
#include "CAT_Octree.h"
#include "CAT_SurfaceIO.h"
#include "CAT_Patch.h"

double dist = 5.0f; /* mm */

/* the argument table */
ArgvInfo argTable[] = {
  {"-dist", ARGV_FLOAT, (char *) 0, (char *) &dist,
    "Minimum distance to define as an artifact"},
  { NULL, ARGV_END, NULL, NULL, NULL }
};


void
usage(char *executable)
{
        char *usage_str =
"\nUsage: %s [options] surf_with_artifacts surf_wo_artifacts output_file\n\n\
    Locate and mark artifacts, using (usually) a smoothed surface for the second argument.  Output is a text file.\n\n";

        fprintf(stderr, usage_str, executable);
}


int
main(int argc, char *argv[])
{
        int                  *n_neighbours, **neighbours;
        char                 *asurf_file, *surf_file, *out_file;
        File_formats         format;
        int                  p, n_objects, n_artifacts, *artifacts;
        polygons_struct      *asurf, *smsurf;
        object_struct        **aobjects, **objects;
        char                 str[80];

        if (ParseArgv(&argc, argv, argTable, 0) || argc != 4) {
                usage(argv[0]);
                fprintf(stderr, "       %s -help\n\n", argv[0]);
                exit(EXIT_FAILURE);
        }

        initialize_argument_processing(argc, argv);

        if (!get_string_argument(NULL, &asurf_file) ||
            !get_string_argument(NULL, &surf_file) ||
            !get_string_argument(NULL, &out_file)) {
                usage(argv[0]);
                exit(EXIT_FAILURE);
        }

        if (input_graphics_any_format(asurf_file, &format,
                                      &n_objects, &aobjects) != OK)
                exit(EXIT_FAILURE);

        if (n_objects != 1 || get_object_type(aobjects[0]) != POLYGONS) {
                printf("Surface file must contain 1 polygons object.\n");
                exit(EXIT_FAILURE);
        }
        asurf = get_polygons_ptr(aobjects[0]);

        if (input_graphics_any_format(surf_file, &format,
                                      &n_objects, &objects) != OK)
                exit(EXIT_FAILURE);

        if (n_objects != 1 || get_object_type(objects[0]) != POLYGONS) {
                printf("Surface file must contain 1 polygons object.\n");
                exit(EXIT_FAILURE);
        }

        smsurf = get_polygons_ptr(objects[0]);

        create_polygon_point_neighbours(asurf, TRUE, &n_neighbours,
                                        &neighbours, NULL, NULL);

        artifacts = (int *) malloc(sizeof(int) * asurf->n_points);
        n_artifacts = find_artifacts(asurf, smsurf, artifacts,
                                     n_neighbours, neighbours, dist);

        printf("%d artifacts found\n", n_artifacts);

        output_values_any_format(out_file, asurf->n_points,
                                 artifacts, TYPE_INTEGER);

        /* clean up */
        free(artifacts);

        delete_polygon_point_neighbours(asurf, n_neighbours,
                                        neighbours, NULL, NULL);

        delete_object_list(1, objects);
        delete_object_list(1, aobjects);
    
        return(EXIT_SUCCESS);
}
