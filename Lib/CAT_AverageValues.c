/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 *
 */

/**
 * @file CAT_AverageValues.c
 * @brief Implementation of surface value/texture averaging.
 *
 * Averages per-vertex surface values (and optionally the underlying
 * mesh geometry) from multiple GIFTI mesh files.  Also provides
 * standard-deviation and z-score computation.
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include "CAT_AverageValues.h"
#include "CAT_SurfaceIO.h"

/* ================================================================== */
/* CAT_AverageValuesCompute                                            */
/* ================================================================== */
int
CAT_AverageValuesCompute(
    char **infiles,
    int nfiles,
    int verbose,
    double **avg_values_out,
    int *n_avg_values_out,
    object_struct **out_object_out)
{
    int i, j;
    int n_pts, n_avg_pts = 0, n_values, n_avg_values = 0;
    int n_objects;
    double *avg_values = NULL, *input_values = NULL;
    File_formats format;
    object_struct **object_list = NULL;
    object_struct *out_object = NULL;
    Point *pts, *avg_pts;

    if (!infiles || nfiles < 1)
        return -1;

    for (i = 0; i < nfiles; i++) {
        if (verbose)
            fprintf(stdout, "%5d/%d: %s\n", i + 1, nfiles, infiles[i]);

        /* check GIFTI extension */
        if (!filename_extension_matches(infiles[i], "gii")) {
            fprintf(stderr, "CAT_AverageValuesCompute: "
                    "%s is not a GIFTI (.gii) file.\n", infiles[i]);
            free(avg_values);
            return -2;
        }

        /* read mesh + texture */
        if (input_gifti_mesh_and_texture(infiles[i], &format,
                                         &n_objects, &object_list,
                                         &n_values, &input_values) != OK) {
            fprintf(stderr, "CAT_AverageValuesCompute: "
                    "could not read mesh and texture from %s\n", infiles[i]);
            free(avg_values);
            return -3;
        }

        n_pts = get_object_points(object_list[0], &pts);

        if (i == 0) {
            /* first file: initialise accumulators */
            out_object   = object_list[0];
            n_avg_pts    = n_pts;
            n_avg_values = n_values;
            get_object_points(out_object, &avg_pts);

            avg_values = (double *) calloc(n_avg_values, sizeof(double));
            if (!avg_values) {
                fprintf(stderr, "CAT_AverageValuesCompute: "
                        "memory allocation failed.\n");
                return -4;
            }

            for (j = 0; j < n_avg_values; j++)
                avg_values[j] = input_values[j] / (double) nfiles;

            free(input_values);
        } else {
            /* subsequent files: check dimensions, accumulate */
            if (n_pts != n_avg_pts) {
                fprintf(stderr, "CAT_AverageValuesCompute: "
                        "mesh point count mismatch in %s (%d vs %d).\n",
                        infiles[i], n_pts, n_avg_pts);
                free(avg_values);
                return -5;
            }
            if (n_values != n_avg_values) {
                fprintf(stderr, "CAT_AverageValuesCompute: "
                        "texture value count mismatch in %s (%d vs %d).\n",
                        infiles[i], n_values, n_avg_values);
                free(avg_values);
                return -6;
            }

            /* accumulate mesh coordinates */
            get_object_points(out_object, &avg_pts);
            for (j = 0; j < n_pts; j++)
                ADD_POINTS(avg_pts[j], avg_pts[j], pts[j]);

            /* accumulate texture values */
            for (j = 0; j < n_avg_values; j++)
                avg_values[j] += input_values[j] / (double) nfiles;

            free(input_values);
            delete_object_list(n_objects, object_list);
        }
    }

    /* finalise mesh average */
    get_object_points(out_object, &avg_pts);
    for (j = 0; j < n_avg_pts; j++)
        SCALE_POINT(avg_pts[j], avg_pts[j], (1.0 / nfiles));

    if (get_object_type(out_object) == POLYGONS)
        compute_polygon_normals(get_polygons_ptr(out_object));

    *avg_values_out   = avg_values;
    *n_avg_values_out = n_avg_values;
    *out_object_out   = out_object;

    return 0;
}

/* ================================================================== */
/* CAT_AverageValuesStd                                                */
/* ================================================================== */
int
CAT_AverageValuesStd(
    char **infiles,
    int nfiles,
    const double *avg_values,
    int n_avg_values,
    double **std_values_out)
{
    int i, j;
    int n_objects, n_values;
    double *sum_squares, *input_values;
    File_formats format;
    object_struct **object_list;

    if (nfiles < 2)
        return -1;

    sum_squares = (double *) calloc(n_avg_values, sizeof(double));
    if (!sum_squares)
        return -2;

    for (i = 0; i < nfiles; i++) {
        if (input_gifti_mesh_and_texture(infiles[i], &format,
                                         &n_objects, &object_list,
                                         &n_values, &input_values) != OK) {
            free(sum_squares);
            return -3;
        }

        for (j = 0; j < n_avg_values; j++)
            sum_squares[j] += input_values[j] * input_values[j];

        free(input_values);
        delete_object_list(n_objects, object_list);
    }

    /* std = sqrt( (sum_x2 - N * mean^2) / (N-1) ) */
    for (j = 0; j < n_avg_values; j++)
        sum_squares[j] = sqrt(
            1.0 / ((double) nfiles - 1.0)
            * (sum_squares[j] - (double) nfiles
               * avg_values[j] * avg_values[j]));

    *std_values_out = sum_squares;
    return 0;
}

/* ================================================================== */
/* CAT_AverageValuesZscore                                             */
/* ================================================================== */
int
CAT_AverageValuesZscore(
    char **infiles,
    int nfiles,
    const double *avg_values,
    const double *std_values,
    int n_avg_values,
    int verbose,
    double *zscores)
{
    int i, j;
    int n_objects, n_values, n_valid;
    double *input_values, diff;
    File_formats format;
    object_struct **object_list;

    for (i = 0; i < nfiles; i++) {
        if (verbose)
            fprintf(stdout, "zscore %5d/%d: %s\n",
                    i + 1, nfiles, infiles[i]);

        if (input_gifti_mesh_and_texture(infiles[i], &format,
                                         &n_objects, &object_list,
                                         &n_values, &input_values) != OK) {
            return -1;
        }

        zscores[i] = 0.0;
        n_valid = 0;
        for (j = 0; j < n_avg_values; j++) {
            if (std_values[j] == 0.0 || isnan(std_values[j]))
                continue;
            diff = input_values[j] - avg_values[j];
            zscores[i] += diff * diff / std_values[j];
            n_valid++;
        }
        if (n_valid > 0)
            zscores[i] /= (double) n_valid;

        free(input_values);
        delete_object_list(n_objects, object_list);
    }

    return 0;
}
