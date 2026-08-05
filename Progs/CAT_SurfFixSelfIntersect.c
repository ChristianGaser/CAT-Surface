/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 *
 */

#include <bicpl.h>
#include <ParseArgv.h>

#include "CAT_SurfaceIO.h"
#include "CAT_Octree.h"
#include "CAT_Intersect.h"

/* argument defaults */
static int max_passes = 10;
static int maxiter = 50;
static double near_threshold = 0.0;
static int verbose = 0;

static ArgvInfo argTable[] = {
    {"-passes", ARGV_INT, (char *)1, (char *)&max_passes,
     "Maximum number of detect/smooth passes."},
    {"-iter", ARGV_INT, (char *)1, (char *)&maxiter,
     "Maximum smoothing iterations per pass."},
    {"-near", ARGV_FLOAT, (char *)1, (char *)&near_threshold,
     "Additionally push apart non-adjacent vertices closer than this multiple\n\
                 of the average edge length (0 = off)."},
    {"-verbose", ARGV_CONSTANT, (char *)1, (char *)&verbose,
     "Print progress of the correction."},
    {NULL, ARGV_END, NULL, NULL, NULL}};

static void
usage(char *executable)
{
    char *usage_str =
        "\nUsage: %s [options] surface_file output_surface_file\n\n\
    Remove self-intersections from a triangle mesh by locally smoothing the\n\
    intersecting regions.  Intersecting triangles are detected, grouped into\n\
    defect regions and then relaxed until they no longer intersect; regions\n\
    that do not improve are expanded before the next attempt.\n\n\
    In contrast to CAT_SurfRemoveIntersections (MeshFix) the mesh topology is\n\
    preserved: the number of vertices and faces and their connectivity are\n\
    unchanged, so per-vertex data such as thickness stays valid.\n\n\
    Use CAT_SurfSelfIntersect to only detect and mark self-intersections.\n\n\
    The corrected surface is always written; the exit code is 1 if any\n\
    self-intersections remain, which allows to check the result in scripts.\n\n";

    fprintf(stderr, usage_str, executable);
}

int
main(int argc, char *argv[])
{
    char *surface_file, *output_file;
    object_struct **objects;
    polygons_struct *polygons;
    int *defects, *polydefects;
    int *n_neighbours, **neighbours;
    int n_objects, n_before, n_after, n_remaining;
    File_formats format;

    if (ParseArgv(&argc, argv, argTable, 0) || argc != 3)
    {
        usage(argv[0]);
        fprintf(stderr, "     %s -help\n\n", argv[0]);
        exit(EXIT_FAILURE);
    }

    initialize_argument_processing(argc, argv);

    if (!get_string_argument(NULL, &surface_file) ||
        !get_string_argument(NULL, &output_file))
    {
        usage(argv[0]);
        exit(EXIT_FAILURE);
    }

    if (max_passes < 1 || maxiter < 1)
    {
        fprintf(stderr, "Number of passes and iterations must be positive.\n");
        exit(EXIT_FAILURE);
    }

    if (input_graphics_any_format(surface_file, &format, &n_objects,
                                  &objects) != OK)
    {
        fprintf(stderr, "Error reading input file %s\n", surface_file);
        exit(EXIT_FAILURE);
    }

    if (n_objects != 1 || get_object_type(objects[0]) != POLYGONS)
    {
        fprintf(stderr, "Input file must contain one polygon object.\n");
        exit(EXIT_FAILURE);
    }

    polygons = get_polygons_ptr(objects[0]);

    /* count self-intersections before the correction */
    defects = (int *)malloc(sizeof(int) * polygons->n_points);
    polydefects = (int *)malloc(sizeof(int) * polygons->n_items);
    create_polygon_point_neighbours(polygons, TRUE, &n_neighbours,
                                    &neighbours, NULL, NULL);

    find_selfintersections(polygons, defects, polydefects, 1);
    n_before = join_intersections(polygons, defects, polydefects,
                                  n_neighbours, neighbours);

    printf("Self intersections before:  %d\n", n_before);

    if (n_before > 0)
    {
        n_remaining = remove_intersections_iter(polygons, max_passes, maxiter,
                                                verbose);
        if (n_remaining > 0)
            fprintf(stderr, "Warning: %d self intersection(s) could not be "
                            "corrected within %d passes.\n",
                    n_remaining, max_passes);
    }

    if (near_threshold > 0.0)
        remove_near_intersections(polygons, near_threshold, verbose);

    /* count self-intersections after the correction */
    find_selfintersections(polygons, defects, polydefects, 1);
    n_after = join_intersections(polygons, defects, polydefects,
                                 n_neighbours, neighbours);

    printf("Self intersections after:   %d\n", n_after);

    if (output_graphics_any_format(output_file, format, n_objects,
                                   objects, NULL) != OK)
    {
        fprintf(stderr, "Error writing output file %s\n", output_file);
        exit(EXIT_FAILURE);
    }

    free(defects);
    free(polydefects);

    delete_polygon_point_neighbours(polygons, n_neighbours,
                                    neighbours, NULL, NULL);

    delete_object_list(n_objects, objects);

    return (n_after > 0);
}
