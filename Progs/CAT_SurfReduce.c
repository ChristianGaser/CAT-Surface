/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 *
 */

/*
 * CAT_SurfReduce.c
 * 
 * Command-line tool to reduce mesh complexity using Quadric Error Metrics (QEM).
 * Calls reduce_mesh_quadrics() from CAT_Surf.c.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ParseArgv.h>
#include "CAT_Surf.h"
#include "CAT_SurfaceIO.h"

/* Defaults */
static double ratio = 0.50;          /* target faces as ratio of input triangles */
static double aggressiveness = 5.0;  /* larger => stronger decimation */
static int preserve_sharp = 1;       /* prevent sharp-edge collapses */
static int verbose = 0;

/* Command-line argument table (CAT style) */
static ArgvInfo argTable[] = {
    {"-ratio", ARGV_FLOAT, (char *) TRUE, (char *) &ratio,
        "Target triangle ratio in (0,1], e.g., 0.5 keeps ~50% of faces (default 0.5)."},
    {"-aggr", ARGV_FLOAT, (char *) TRUE, (char *) &aggressiveness,
        "Aggressiveness of simplification (default 7.0)."},
    {"-preserve_sharp", ARGV_CONSTANT, (char *) TRUE, (char *) &preserve_sharp,
        "Preserve sharp features (default on)."},
    {"-no_preserve_sharp", ARGV_CONSTANT, (char *) FALSE, (char *) &preserve_sharp,
        "Allow collapses across sharp features."},
    {"-verbose", ARGV_CONSTANT, (char *) TRUE, (char *) &verbose,
        "Enable verbose output."},
    {NULL, ARGV_END, NULL, NULL, NULL}
};

/**
 * @brief Print usage for CAT_SurfReduce.
 *
 * Usage:
 *   CAT_SurfReduce [options] <input_surface> <output_surface>
 *
 * Options:
 *   -ratio <float>           Target triangle ratio in (0,1] (default 0.5).
 *   -aggr <float>            Aggressiveness of simplification (default 7.0).
 *   -preserve_sharp          Prevent aggressive collapses across sharp edges (default).
 *   -no_preserve_sharp       Allow collapses across sharp features.
 *   -verbose                 Verbose logging.
 *
 * Notes:
 *   - The absolute target face count is computed as round(ratio * #input_triangles).
 *   - Non-triangle polygons are fan-triangulated internally by reduce_mesh_quadrics().
 *   - Input/output formats are auto-detected by CAT_SurfaceIO.
 */
static void usage(const char *prog)
{
    fprintf(stderr,
        "\n"
        "-----------------------------------------------------------------------------\n"
        " CAT_SurfReduce - Quadric Error Metrics (QEM) Mesh Reduction\n"
        "-----------------------------------------------------------------------------\n"
        "\n"
        "Usage:\n"
        "  %s [options] <input_surface> <output_surface>\n"
        "\n"
        "Options:\n"
        "  -ratio <float>         Target triangle ratio in (0,1] (default 0.5)\n"
        "  -aggr <float>          Aggressiveness of simplification (default 7.0)\n"
        "  -preserve_sharp        Preserve sharp features (default ON)\n"
        "  -no_preserve_sharp     Disable sharp-feature preservation\n"
        "  -verbose               Verbose output\n"
        "\n"
        "Example:\n"
        "  %s -ratio 0.25 -aggr 8.0 -preserve_sharp input.gii output.gii\n"
        "\n"
        "Notes:\n"
        "  * Target faces = round(ratio * input_triangles). Limits: [1, input_triangles].\n"
        "  * Non-triangle faces are triangulated internally.\n"
        "  * Supported formats depend on CAT_SurfaceIO (e.g., GIfTI, MNI OBJ, etc.).\n"
        "-----------------------------------------------------------------------------\n",
        prog, prog);
}

int
main(int argc, char *argv[])
{
    char *in_surface = NULL;
    char *out_surface = NULL;
    int n_objects = 0;
    int target_faces;
    File_formats fmt;
    object_struct **out_object, **obj_list = NULL;
    polygons_struct *polygons_out, *polygons = NULL;

    initialize_argument_processing(argc, argv);

    if (ParseArgv(&argc, argv, argTable, 0)) {
        usage(argv[0]);
        fprintf(stderr, "       %s -help\n\n", argv[0]);
        return EXIT_FAILURE;
    }

    if (!get_string_argument(NULL, &in_surface) ||
        !get_string_argument(NULL, &out_surface)) {
        usage(argv[0]);
        fprintf(stderr, "Usage: %s [options] input_surface output_surface\n", argv[0]);
        return EXIT_FAILURE;
    }

    if (ratio <= 0.0) {
        fprintf(stderr, "Error: -ratio must be a value >0, got %.6f\n", ratio);
        return EXIT_FAILURE;
    }

    /* Read input surface (must be a single POLYGONS object) */
    if (input_graphics_any_format(in_surface, &fmt, &n_objects, &obj_list) == ERROR ||
        n_objects != 1 || obj_list[0]->object_type != POLYGONS) {
        fprintf(stderr, "Error: file must contain exactly one polygons object: %s\n", in_surface);
        return EXIT_FAILURE;
    }

    polygons = get_polygons_ptr(obj_list[0]);

    /* Compute absolute target from ratio */
    int nf0 =  polygons->n_items;
    if (ratio < 1.0)
        target_faces = (int)round(ratio * nf0);
    else
        target_faces = (int)ratio;
        
    //if (target_faces > nf0) target_faces = nf0;

    if (verbose) {
        fprintf(stderr, "[Reduce] input faces (triangulated): %d\n", nf0);
        fprintf(stderr, "[Reduce] ratio=%.6f -> target=%d faces\n", ratio, target_faces);
        fprintf(stderr, "[Reduce] aggressiveness=%.3f, preserve_sharp=%d\n", aggressiveness, preserve_sharp);
    }

    int en, en0 = euler_characteristic(polygons);
    en = en0;
    out_object  = (object_struct **) malloc(sizeof(object_struct *));
    *out_object = create_object(POLYGONS);
    polygons_out = get_polygons_ptr(*out_object);

    /* Perform simplification */
    while ((en == en0) && (polygons->n_items > target_faces)) {
        if (reduce_mesh_quadrics( polygons, ROUND(0.9*polygons->n_items), aggressiveness, preserve_sharp, verbose) != 0) {
            fprintf(stderr, "Error: reduce_mesh_quadrics failed.\n");
            return EXIT_FAILURE;
        }
    
        (void)remove_duplicate_edges(polygons, verbose);
        en = euler_characteristic(polygons);
        if (en == en0) copy_polygons(polygons, polygons_out);
    }
    compute_polygon_normals(polygons_out);
        
    /* Write output surface */
    if (output_graphics_any_format(out_surface, ASCII_FORMAT,
                                   n_objects, obj_list, NULL) == ERROR) {
        fprintf(stderr, "Error: writing output surface: %s\n", out_surface);
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
