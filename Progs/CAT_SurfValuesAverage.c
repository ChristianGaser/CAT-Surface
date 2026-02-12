/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 *
 */

/* Average surface values/textures (and underlying meshes) from resampled
 * GIFTI mesh files.  This is a thin CLI wrapper around the library functions
 * in CAT_AverageValues.
 *
 * Input  : one or more GIFTI mesh files that each contain a mesh (POINTSET +
 *          TRIANGLE) and at least one SHAPE data array (texture/value).
 *          Files are typically resampled surfaces whose names contain "mesh".
 *          Only GIFTI files (.gii) with embedded texture data are accepted.
 * Output : a single GIFTI file with the averaged mesh and the averaged
 *          texture values.
 */

#include <bicpl.h>
#include <ParseArgv.h>

#include "CAT_SurfaceIO.h"
#include "CAT_AverageValues.h"
#include "CAT_SafeAlloc.h"

#include <stdlib.h>
#include <string.h>

/* ------------------------------------------------------------------ */
/* argument handling                                                   */
/* ------------------------------------------------------------------ */
int verbose = 0;

char *std_filename        = NULL;
char *zscore_filename     = NULL;
char *zscore_txt_filename = NULL;

static ArgvInfo argTable[] = {
    {"-std", ARGV_STRING, (char *) 1, (char *) &std_filename,
          "Write standard deviation of values to file."},
    {"-zscore", ARGV_STRING, (char *) 1, (char *) &zscore_filename,
          "Write filenames and z-scores in csv-file."},
    {"-zscore-txt", ARGV_STRING, (char *) 1, (char *) &zscore_txt_filename,
          "Write z-scores in txt-file."},
    {"-v", ARGV_CONSTANT, (char *) 1, (char *) &verbose,
          "Be verbose."},
    {NULL, ARGV_END, NULL, NULL, NULL}
};

/* ================================================================== */
/* main                                                                */
/* ================================================================== */
int
main(int argc, char *argv[])
{
    char **infiles, *outfile;
    int i, nfiles, n_avg_values;
    double *avg_values, *std_values, *zscore;
    object_struct *out_object;
    FILE *fid;

    /* ---- parse arguments ---- */
    if (ParseArgv(&argc, argv, argTable, 0) || (argc < 3)) {
        fprintf(stderr,
            "\nUsage: %s [options] in1.gii ... inN.gii out.gii\n\n"
            "Average values/textures (and meshes) from resampled GIFTI mesh\n"
            "files.  Only GIFTI files that contain both a surface mesh and\n"
            "texture (SHAPE) data arrays are accepted.\n\n"
            "Options:\n"
            "  -std <file>         Write standard deviation of values.\n"
            "  -zscore <file.csv>  Write filenames and z-scores as CSV.\n"
            "  -zscore-txt <file>  Write z-scores as plain text.\n"
            "  -v                  Verbose output.\n\n"
            "  %s -help\n\n",
            argv[0], argv[0]);
        exit(EXIT_FAILURE);
    }

    nfiles  = argc - 2;
    infiles = &argv[1];
    outfile = argv[argc - 1];

    if (nfiles == 0) {
        fprintf(stderr, "Error: No input files specified.\n");
        exit(EXIT_FAILURE);
    }

    /* ---- compute average values and mesh ---- */
    if (CAT_AverageValuesCompute(infiles, nfiles, verbose,
                                 &avg_values, &n_avg_values,
                                 &out_object) != 0) {
        fprintf(stderr, "Error: Failed computing average.\n");
        exit(EXIT_FAILURE);
    }

    /* ---- standard deviation ---- */
    if (std_filename || zscore_filename || zscore_txt_filename) {
        if (CAT_AverageValuesStd(infiles, nfiles,
                                 avg_values, n_avg_values,
                                 &std_values) != 0) {
            fprintf(stderr, "Error: Failed computing standard deviation.\n");
            exit(EXIT_FAILURE);
        }

        if (std_filename) {
            output_values_any_format(std_filename, n_avg_values,
                                     std_values, TYPE_DOUBLE);
        }

        /* ---- z-scores ---- */
        if (zscore_filename || zscore_txt_filename) {
            zscore = (double *) malloc(sizeof(double) * nfiles);

            if (CAT_AverageValuesZscore(infiles, nfiles,
                                        avg_values, std_values,
                                        n_avg_values, verbose,
                                        zscore) != 0) {
                fprintf(stderr, "Error: Failed computing z-scores.\n");
                exit(EXIT_FAILURE);
            }

            if (zscore_filename) {
                fid = SAFE_FOPEN(zscore_filename, "w");
                fprintf(fid, "filename,z-score\n");
                for (i = 0; i < nfiles; i++)
                    fprintf(fid, "%s,%g\n", infiles[i], zscore[i]);
                fclose(fid);
            }
            if (zscore_txt_filename) {
                fid = SAFE_FOPEN(zscore_txt_filename, "w");
                for (i = 0; i < nfiles; i++)
                    fprintf(fid, "%g\n", zscore[i]);
                fclose(fid);
            }

            free(zscore);
        }

        free(std_values);
    }

    /* ---- write output: averaged mesh + averaged values ---- */
    if (output_graphics_any_format(outfile, ASCII_FORMAT,
                                   1, &out_object, avg_values) != OK) {
        fprintf(stderr, "Error: Could not write output to %s\n", outfile);
        free(avg_values);
        exit(EXIT_FAILURE);
    }

    free(avg_values);

    return EXIT_SUCCESS;
}
