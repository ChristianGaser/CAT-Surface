/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 *
 */

#ifndef _CAT_AVERAGE_VALUES_H_
#define _CAT_AVERAGE_VALUES_H_

/**
 * @file CAT_AverageValues.h
 * @brief Surface value/texture averaging library.
 *
 * Provides functions for computing the mean, standard deviation, and
 * per-file z-scores of vertex-wise surface values across multiple
 * GIFTI mesh files (i.e. files that carry both a mesh and SHAPE data).
 *
 * The underlying meshes are averaged simultaneously.
 */

#include <bicpl.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Compute the average of vertex-wise surface values (and their mesh)
 * from multiple GIFTI files.
 *
 * Each input file must be a GIFTI (.gii) file containing a mesh
 * (POINTSET + TRIANGLE) and at least one SHAPE data array.
 *
 * @param infiles       Array of input GIFTI filenames.
 * @param nfiles        Number of input files.
 * @param verbose       If non-zero, print progress to stdout.
 * @param avg_values    Output: averaged vertex values (caller must free).
 * @param n_avg_values  Output: number of averaged values.
 * @param out_object    Output: averaged mesh object (caller must free).
 *
 * @return 0 on success, non-zero on error.
 */
int CAT_AverageValuesCompute(
    char **infiles,
    int nfiles,
    int verbose,
    double **avg_values,
    int *n_avg_values,
    object_struct **out_object
);

/**
 * Compute per-vertex standard deviation from pre-computed average.
 *
 * Re-reads the input files to accumulate sum-of-squares, then computes
 * the sample standard deviation  sqrt( (sum_x2 - N*mean^2) / (N-1) ).
 *
 * @param infiles       Array of input GIFTI filenames.
 * @param nfiles        Number of input files.
 * @param avg_values    Average values from CAT_AverageValuesCompute().
 * @param n_avg_values  Number of average values.
 * @param std_values    Output: per-vertex standard deviation (caller must free).
 *
 * @return 0 on success, non-zero on error.
 */
int CAT_AverageValuesStd(
    char **infiles,
    int nfiles,
    const double *avg_values,
    int n_avg_values,
    double **std_values
);

/**
 * Compute per-file z-score from pre-computed average and std.
 *
 * For each file, z-score = mean_vertex( (x - avg)^2 / std ) over
 * vertices where std > 0.
 *
 * @param infiles       Array of input GIFTI filenames.
 * @param nfiles        Number of input files.
 * @param avg_values    Average values from CAT_AverageValuesCompute().
 * @param std_values    Standard deviation from CAT_AverageValuesStd().
 * @param n_avg_values  Number of values per file.
 * @param verbose       If non-zero, print progress to stdout.
 * @param zscores       Output: per-file z-scores (caller must allocate nfiles doubles).
 *
 * @return 0 on success, non-zero on error.
 */
int CAT_AverageValuesZscore(
    char **infiles,
    int nfiles,
    const double *avg_values,
    const double *std_values,
    int n_avg_values,
    int verbose,
    double *zscores
);

#ifdef __cplusplus
}
#endif

#endif /* _CAT_AVERAGE_VALUES_H_ */
