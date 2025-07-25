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
#define PRE_ONLY   1
#define POST_ONLY  2

#define TOLERANCE  0.001

BOOLEAN prepostflag = PRE_AND_POST; /* default: both pre + post processing */

static ArgvInfo argTable[] = {
  {"-preonly", ARGV_CONSTANT, (char *) PRE_ONLY, (char *) &prepostflag,
  "execute only the preprocessing pipeline" },
  {"-postonly", ARGV_CONSTANT, (char *) POST_ONLY, (char *) &prepostflag,
  "execute only the postprocessing pipeline" },
   {NULL, ARGV_END, NULL, NULL, NULL}
};

void
usage(char *executable)
{
    char *usage_str =
"\nUsage: %s [options] surface_file conformalmap_file output_surface_file\n\n\
  Adjust a spherical map to be more isometric (area preserving).\n\
  The original mesh is surfcae_file and the conformal map of the mesh is\n\
  conformalmap_file.  Results are saved in output_surface_file.\n\n";

    fprintf(stderr, usage_str, executable);
}

int
main(int argc, char** argv)
{
    char         *input_file, *cmap_file, *output_surface_file;
    object_struct    **objects;
    polygons_struct    *polygons, *map;
    struct metricdata  *brain;
    File_formats     format;
    int          n_objects;
    int          p, iters, count;
    double         radius;

    if (ParseArgv(&argc, argv, argTable, 0) || (argc < 3)) {
        usage(argv[0]);
        fprintf(stderr, "   %s -help\n\n", argv[0]);
        exit(EXIT_FAILURE);
    }

    initialize_argument_processing(argc, argv);
    if (!get_string_argument(NULL, &input_file) ||
      !get_string_argument(NULL, &cmap_file) ||
      !get_string_argument(NULL, &output_surface_file)) {
        fprintf(stderr, "\nUsage: %s [options] surface_file conformalmap_file output_surface_file\n", argv[0]);
        exit(EXIT_FAILURE);
    }

    if (input_graphics_any_format(input_file, &format, &n_objects,
                    &objects) != OK) {
        printf("Error reading input file\n");
        exit(EXIT_FAILURE);
    }

    if (n_objects != 1 || get_object_type(objects[0]) != POLYGONS) {
        printf("Input file must contain one polygon object.\n");
        exit(EXIT_FAILURE);
    }

    polygons = get_polygons_ptr(objects[0]);
    compute_polygon_normals(polygons);

    if (input_graphics_any_format(cmap_file, &format, &n_objects,
                    &objects) != OK) {
        printf("Error reading conformal map file\n");
        exit(EXIT_FAILURE);
    }

    if (n_objects != 1 || get_object_type(objects[0]) != POLYGONS) {
        printf("Conformal map file must contain one polygon object.\n");
        exit(EXIT_FAILURE);
    }

    map = get_polygons_ptr(objects[0]);
    compute_polygon_normals(map);

    if (polygons->n_points != map->n_points
      || polygons->n_items != map->n_items) {
        printf("Input mesh and conformal map mesh do not match.\n");
        exit(EXIT_FAILURE);
    }
  
    radius = sqrt(get_polygons_surface_area(map) / (4.0 * PI));
    for (p = 0; p < map->n_points; p++)
        set_vector_length(&map->points[p], radius);

    brain = getmetricdata(polygons);

    if (prepostflag != POST_ONLY) {
        distortcorrect(brain, map, 100000, SELECT_OFF, TOLERANCE);
    }

    if (prepostflag != PRE_ONLY) {
        count = 0;
        do {
            iters = 0;
            count++;
            iters += smooth(brain, map, 1000, SELECT_ON, TOLERANCE);
            stretch(brain, map, 1000, SELECT_ON, ~LARGE_ONLY,
                TOLERANCE);
        } while (iters > 0 && count < 5);

        do {
            iters = 0;
            count++;
            iters += stretch(brain, map, 1000, SELECT_OFF,
                     LARGE_ONLY, TOLERANCE);
            stretch(brain, map, 1000, SELECT_ON, ~LARGE_ONLY,
                TOLERANCE);
        } while (iters > 0 && count < 10);

        do {
            iters = 0;
            count++;
            iters += distortcorrect(brain, map, 1000, SELECT_ON,
                        TOLERANCE);
            stretch(brain, map, 1000, SELECT_ON, ~LARGE_ONLY,
                TOLERANCE);
        } while (iters > 0 && count < 15);
    }

    compute_polygon_normals(map);
    if(output_graphics_any_format(output_surface_file, format, 1, 
             objects, NULL) != OK)
        exit(EXIT_FAILURE);

    delete_object_list(n_objects, objects);

    return(EXIT_SUCCESS);
}
