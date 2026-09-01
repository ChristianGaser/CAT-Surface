/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 *
 */

/*
 * CAT_SurfInfo.c
 *
 * Command-line tool that prints geometric and topological information about a
 * surface: number of vertices, faces and edges, surface area, Euler number and
 * the number of self-intersections.  Calls compute_surf_info() from
 * CAT_SurfInfo.c.
 */

#include <stdio.h>
#include <stdlib.h>
#include <bicpl.h>
#include <ParseArgv.h>

#include "CAT_SurfInfo.h"
#include "CAT_SurfaceIO.h"

/* argument defaults */
static int no_intersect = 0;
static int tabular = 0;
static int verbose = 0;

/* the argument table */
static ArgvInfo argTable[] = {
    {"-no-selfintersect", ARGV_CONSTANT, (char *) TRUE, (char *) &no_intersect,
        "Skip the self-intersection test, which is the expensive part."},
    {"-tab", ARGV_CONSTANT, (char *) TRUE, (char *) &tabular,
        "Print tab-separated key/value pairs instead of a report."},
    {"-verbose", ARGV_CONSTANT, (char *) TRUE, (char *) &verbose,
        "Print progress information."},
    {NULL, ARGV_END, NULL, NULL, NULL}
};

/**
 * \brief Print usage for CAT_SurfInfo.
 *
 * \param executable (in) name the tool was invoked with
 * \return void (no return value).
 */
static void
usage(const char *executable)
{
    static char *usage_str = "\n\
Usage: %s [options] surface_file\n\n\
Print information about a surface: number of vertices, faces and edges,\n\
surface area, enclosed volume, bounding box, Euler number, genus, connected\n\
components and the number of self-intersections.  For a GIFTI file the\n\
DataArrays it carries besides the mesh are listed as well, so embedded\n\
per-vertex data (thickness, curvature, labels) becomes visible.\n\n\
The self-intersection test dominates the runtime on large meshes, use\n\
-no-selfintersect to skip it.  With -tab the same numbers are written as\n\
tab-separated key/value pairs, one per line, for use in scripts.\n\n\
Example:\n\
  %s -no-selfintersect lh.central.gii\n\n";

    fprintf(stderr, usage_str, executable, executable);
}

int
main(int argc, char *argv[])
{
    char *surface_file;
    File_formats format;
    int n_objects;
    object_struct **objects;
    polygons_struct *polygons;
    surf_info_struct info;
    gifti_darray_info *darrays = NULL;
    int n_darrays = 0;

    initialize_argument_processing(argc, argv);

    if (ParseArgv(&argc, argv, argTable, 0)) {
        usage(argv[0]);
        fprintf(stderr, "   %s -help\n\n", argv[0]);
        exit(EXIT_FAILURE);
    }

    if (!get_string_argument(NULL, &surface_file)) {
        usage(argv[0]);
        fprintf(stderr, "   %s -help\n\n", argv[0]);
        exit(EXIT_FAILURE);
    }

    if (input_graphics_any_format(surface_file, &format, &n_objects,
                                  &objects) != OK) {
        fprintf(stderr, "Error reading input file %s\n", surface_file);
        exit(EXIT_FAILURE);
    }

    if (n_objects != 1 || get_object_type(objects[0]) != POLYGONS) {
        fprintf(stderr, "File must contain 1 polygons object.\n");
        exit(EXIT_FAILURE);
    }

    polygons = get_polygons_ptr(objects[0]);

    compute_surf_info(polygons, &info, !no_intersect, verbose);

    if (filename_extension_matches(surface_file, "gii")) {
        if (input_gifti_darrays(surface_file, &n_darrays, &darrays) != OK)
            fprintf(stderr, "Warning: cannot list the GIFTI data arrays of %s\n",
                    surface_file);
    }

    if (tabular) {
        print_surf_info_tabular(&info, surface_file, stdout);
        if (darrays != NULL)
            print_gifti_darrays_tabular(darrays, n_darrays, stdout);
    } else {
        print_surf_info(&info, surface_file, stdout);
        if (darrays != NULL)
            print_gifti_darrays(darrays, n_darrays, stdout);
        printf("\n");
    }

    free(darrays);

    delete_object_list(n_objects, objects);

    return EXIT_SUCCESS;
}
