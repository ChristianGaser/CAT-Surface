/* Rachel Yotter - rachel.yotter@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 *
 */

/* Program to make a conformal map more isometric. */

#include "CAT_Isometrics.h"

#define PRE_AND_POST 0
#define PRE_ONLY     1
#define POST_ONLY    2

BOOLEAN prepostflag = PRE_AND_POST; /* default: both pre + post processing */
BOOLEAN quiet = 0; /* turn progress reports on and off */

static ArgvInfo argTable[] = {
  {"-preonly", ARGV_CONSTANT, (char *) PRE_ONLY, (char *) &prepostflag,
    "execute only the preprocessing pipeline" },
  {"-postonly", ARGV_CONSTANT, (char *) POST_ONLY, (char *) &prepostflag,
    "execute only the postprocessing pipeline" },
  {"-quiet", ARGV_CONSTANT, (char *) TRUE, (char *) &quiet,
    "turn off progress reports" },
   {NULL, ARGV_END, NULL, NULL, NULL}
};

void
usage(char *executable)
{
        char *usage_str =
"\nUsage: %s [options] infile.obj conformalmap.obj outfile.obj\n\n\
    Adjust a spherical map to be more isometric (area preserving).\n\
    The original mesh is infile.obj and the conformal map of the mesh is\n\
    conformalmap.obj.  Results are saved in outfile.obj.\n\n";

        fprintf(stderr, usage_str, executable);
}

int
main(int argc, char** argv)
{
        char               *input_file, *cmap_file, *output_file;
        object_struct      **objects;
        polygons_struct    *polygons, *map;
        struct metricdata  *brain;
        File_formats       format;
        int                n_objects;

        if (ParseArgv(&argc, argv, argTable, 0) || (argc < 3)) {
                usage(argv[0]);
                fprintf(stderr, "       %s -help\n\n", argv[0]);
                return(1);
        }

        initialize_argument_processing(argc, argv);
        if (!get_string_argument(NULL, &input_file) ||
            !get_string_argument(NULL, &cmap_file) ||
            !get_string_argument(NULL, &output_file)) {
                fprintf(stderr, "\nUsage: %s [options] infile.obj conformalmap.obj outfile.obj\n", argv[0]);
                return(1);
        }

        if (input_graphics_any_format(input_file, &format, &n_objects,
                                      &objects) != OK) {
                printf("Error reading input file\n");
                return(1);
        }

        if (n_objects != 1 || get_object_type(objects[0]) != POLYGONS) {
                printf("Input file must contain one polygon object.\n");
                return(1);
        }

        polygons = get_polygons_ptr(objects[0]);
        compute_polygon_normals(polygons);

        if (input_graphics_any_format(cmap_file, &format, &n_objects,
                                      &objects) != OK) {
                printf("Error reading conformal map file\n");
                return(1);
        }

        if (n_objects != 1 || get_object_type(objects[0]) != POLYGONS) {
                printf("Conformal map file must contain one polygon object.\n");
                return(1);
        }

        map = get_polygons_ptr(objects[0]);
        compute_polygon_normals(map);

        if (polygons->n_points != map->n_points
            || polygons->n_items != map->n_items) {
                printf("Input mesh and conformal map mesh do not match.\n");
                return(1);
        }
    
        brain = getmetricdata(polygons);

        if (prepostflag != POST_ONLY)
                distortcorrect(brain, map, 15000, SELECT_OFF, quiet);
        if (prepostflag != PRE_ONLY) {
                smooth(brain, map, 100, SELECT_ON, quiet);
                stretch(brain, map, 100, SELECT_OFF, quiet, LARGE_ONLY);
                distortcorrect(brain, map, 100, SELECT_ON, quiet);
                stretch(brain, map, 100, SELECT_ON, quiet, ~LARGE_ONLY);
                distortcorrect(brain, map, 100, SELECT_ON, quiet);
                stretch(brain, map, 100, SELECT_ON, quiet, ~LARGE_ONLY);
                distortcorrect(brain, map, 100, SELECT_ON, quiet);
        }

        compute_polygon_normals(map);
        output_graphics_any_format(output_file, format, 1, objects);

        delete_object_list(n_objects, objects);

        return(0);
}
