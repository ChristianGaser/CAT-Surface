/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 *
 */

#include <bicpl.h>

#include "CAT_MeshClean.h"
#include "CAT_SurfaceIO.h"

static void
usage(char *executable)
{
    char *usage_str = "\n\
Usage: %s surface_file output_file [options]\n\n\
    Remove self-intersections and degeneracies from a triangle mesh\n\
    using the MeshFix algorithm.\n\n\
    Options:\n\
        -iter N      Maximum iterations for cleaning (default: 10)\n\
        -inner N     Inner loop iterations (default: 3)\n\
        -nofill      Do not fill boundary holes\n\
        -count       Only count intersections, do not fix (no output)\n\
        -v           Verbose output\n\n";

    fprintf(stderr, usage_str, executable);
}

int
main(int argc, char *argv[])
{
    char *input_file, *output_file;
    File_formats format;
    int n_objects;
    object_struct **objects;
    polygons_struct *polygons;
    CAT_MeshCleanOptions opts;
    int count_only = 0;
    int i;

    /* Initialize options */
    CAT_MeshCleanOptionsInit(&opts);

    /* Parse arguments */
    if (argc < 2) {
        usage(argv[0]);
        return 1;
    }

    input_file = argv[1];
    output_file = NULL;

    /* Check for -count first (no output file needed) */
    for (i = 2; i < argc; i++) {
        if (strcmp(argv[i], "-count") == 0) {
            count_only = 1;
        }
    }

    if (!count_only) {
        if (argc < 3) {
            usage(argv[0]);
            return 1;
        }
        output_file = argv[2];
    }

    /* Parse remaining options */
    for (i = (count_only ? 2 : 3); i < argc; i++) {
        if (strcmp(argv[i], "-iter") == 0 && i + 1 < argc) {
            opts.max_iters = atoi(argv[++i]);
        } else if (strcmp(argv[i], "-inner") == 0 && i + 1 < argc) {
            opts.inner_loops = atoi(argv[++i]);
        } else if (strcmp(argv[i], "-nofill") == 0) {
            opts.fill_holes = 0;
        } else if (strcmp(argv[i], "-count") == 0) {
            count_only = 1;
        } else if (strcmp(argv[i], "-v") == 0) {
            opts.verbose = 1;
        } else if (argv[i][0] == '-') {
            fprintf(stderr, "Unknown option: %s\n", argv[i]);
            usage(argv[0]);
            return 1;
        }
    }

    /* Load input surface */
    if (input_graphics_any_format(input_file, &format, &n_objects, &objects) != OK ||
        n_objects != 1 || get_object_type(objects[0]) != POLYGONS) {
        fprintf(stderr, "Error reading surface file: %s\n", input_file);
        return 1;
    }

    polygons = get_polygons_ptr(objects[0]);

    if (opts.verbose) {
        printf("Input: %d vertices, %d polygons\n",
               polygons->n_points, polygons->n_items);
    }

    if (count_only) {
        /* Just count intersections */
        int n_intersect = CAT_SurfCountIntersections(polygons);
        if (n_intersect < 0) {
            fprintf(stderr, "Error counting intersections\n");
            delete_object_list(n_objects, objects);
            return 1;
        }
        printf("%d\n", n_intersect);
        delete_object_list(n_objects, objects);
        return 0;
    }

    /* Run mesh cleaning */
    int result = CAT_SurfMeshClean(polygons, &opts);
    
    if (result < 0) {
        fprintf(stderr, "Error during mesh cleaning\n");
        delete_object_list(n_objects, objects);
        return 1;
    }

    if (result > 0 && opts.verbose) {
        fprintf(stderr, "Warning: Some issues could not be fixed\n");
    }

    if (opts.verbose) {
        printf("Output: %d vertices, %d polygons\n",
               polygons->n_points, polygons->n_items);
    }

    /* Recompute normals after modification */
    compute_polygon_normals(polygons);

    /* Save output */
    if (output_graphics_any_format(output_file, format, n_objects, objects, NULL) != OK) {
        fprintf(stderr, "Error writing output file: %s\n", output_file);
        delete_object_list(n_objects, objects);
        return 1;
    }

    delete_object_list(n_objects, objects);
    
    return result;
}
