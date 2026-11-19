/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 *
 */

#ifndef _CAT_VOLUMESTATS_H_
#define _CAT_VOLUMESTATS_H_

/**
 * \file CAT_VolumeStats.h
 * \brief Voxel-wise mean, standard deviation and z-scores of NIfTI volumes.
 *
 * A sample is one 3D volume. A 3D file contributes one sample and a 4D file
 * one sample per volume, so the statistics of a 4D series (for instance the
 * z-maps of a normative model, one volume per subject) are taken over its
 * volumes. All files must have the same matrix size; they may be mixed.
 *
 * Non-finite values of float images are read as 0 by the NIfTI library, so a
 * NaN-masked background shows up as 0 in the mean and the standard deviation.
 *
 * Typical use: CAT_VolAverageScan(), optionally CAT_VolAverageSelect() to use
 * a part of the samples only, CAT_VolAverageCompute(), optionally
 * CAT_VolAverageZscore(), write res.mean / res.std with write_nifti_double()
 * using res.header as template, then CAT_VolAverageFree().
 * CAT_VolAverageRead() gives the samples themselves, one pointer per volume,
 * for callers that compute something else across them (CAT_VolCalc).
 */

#include <stddef.h>

#include "CAT_NiftiLib.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * \brief Geometry of the inputs and the voxel-wise statistics.
 *
 * Filled by CAT_VolAverageScan() (everything but mean, std and selected),
 * CAT_VolAverageSelect() (selected and n_selected) and
 * CAT_VolAverageCompute() (mean and std). Release with CAT_VolAverageFree().
 */
typedef struct
{
    int          dims[3];      /**< nx, ny, nz of every volume */
    double       voxelsize[3]; /**< dx, dy, dz of the first input */
    size_t       n_voxels;     /**< voxels of one volume, nx*ny*nz */
    int          n_samples;    /**< volumes over all input files */
    int         *n_volumes;    /**< volumes of each input file (nfiles entries) */
    nifti_image *header;       /**< header of the first input, template for writing results */
    double      *mean;         /**< voxel-wise mean (n_voxels values) */
    double      *std;          /**< voxel-wise standard deviation (n_voxels values) or NULL */
    unsigned char *selected;   /**< one flag per sample, 1 for the samples that are used;
                                    NULL (no selection) uses all */
    int          n_selected;   /**< samples that are used; n_samples without a selection */
} CAT_VolAverageResult;

/**
 * \brief Read the headers of all inputs, check them and count the samples.
 *
 * Cheap: no voxel data is read. Fails if a file cannot be read, has more than
 * four dimensions or a different matrix size than the first file. A voxel size
 * that differs from the first file is only reported (with verbose).
 *
 * \param infiles (in)  input file names
 * \param nfiles  (in)  number of input files (at least 1)
 * \param verbose (in)  non-zero to report differing voxel sizes
 * \param res     (out) dims, voxelsize, n_voxels, n_samples, n_volumes and header
 *                      are set; release with CAT_VolAverageFree(). Zeroed on error.
 * \return 0 on success, -1 on error (a message is written to stderr)
 */
int CAT_VolAverageScan(char **infiles, int nfiles, int verbose,
                       CAT_VolAverageResult *res);

/**
 * \brief Use only a part of the samples, chosen by a label per sample.
 *
 * The selection is given as "FILE" or "FILE=VALUE[,VALUE...]". FILE is a text
 * file with one label per sample, one per line, in the order of the samples:
 * by input file and, within a 4D file, by volume. Empty lines are skipped and
 * white space around a label is ignored. A sample is used if its label is one
 * of the values; without values, if its label is a number other than 0. Labels
 * and values are compared as numbers if all of them are numbers (so 2 matches
 * 2.0), otherwise as text. If a file of the full name exists, the whole
 * argument is the file name, also when it contains '='.
 *
 * Afterwards CAT_VolAverageCompute(), CAT_VolAverageZscore() and
 * CAT_VolAverageRead() work on the selected samples only, as if the others had
 * not been given. A new selection replaces the previous one and discards the
 * mean and the standard deviation.
 *
 * \param selection (in)     "FILE" or "FILE=VALUE[,VALUE...]"
 * \param res       (in/out) result of CAT_VolAverageScan(); selected and
 *                           n_selected are set
 * \return 0 on success; -1 if the file cannot be read, does not have one label
 *         per sample, has labels that are not numbers although no values are
 *         given, or if no sample is selected (a message is written to stderr
 *         and the previous selection, if any, is kept)
 */
int CAT_VolAverageSelect(const char *selection, CAT_VolAverageResult *res);

/**
 * \brief Compute the voxel-wise mean and, optionally, standard deviation.
 *
 * Every volume is one sample, the standard deviation uses n - 1. It is clamped
 * at 0 where the one-pass formula comes out slightly negative from rounding, as
 * it does for voxels that are constant across the samples. After
 * CAT_VolAverageSelect() only the selected samples are used.
 *
 * \param infiles  (in)     input file names, as given to CAT_VolAverageScan()
 * \param nfiles   (in)     number of input files
 * \param want_std (in)     non-zero to compute the standard deviation as well;
 *                          needs at least two (selected) samples
 * \param verbose  (in)     non-zero to print progress to stdout
 * \param res      (in/out) result of CAT_VolAverageScan(); mean (and std) are set
 * \return 0 on success, -1 on error (a message is written to stderr)
 */
int CAT_VolAverageCompute(char **infiles, int nfiles, int want_std, int verbose,
                          CAT_VolAverageResult *res);

/**
 * \brief Compute one z-score per sample.
 *
 * The z-score of a sample is the fourth root of the mean of z^4, with
 * z = (x - mean) / std, over the voxels where std is positive and the mean
 * exceeds a quarter of the average of the mean image over the voxels brighter
 * than an eighth of its non-zero average (the background is left out).
 * Samples are ordered by file and, within a 4D file, by volume. After
 * CAT_VolAverageSelect() the samples that are not selected get NaN.
 *
 * \param infiles (in)  input file names, as given to CAT_VolAverageScan()
 * \param nfiles  (in)  number of input files
 * \param res     (in)  result of CAT_VolAverageCompute() with the standard deviation
 * \param verbose (in)  non-zero to print progress to stdout
 * \param zscores (out) res->n_samples z-scores, allocated by the caller
 * \return 0 on success, -1 on error (a message is written to stderr)
 */
int CAT_VolAverageZscore(char **infiles, int nfiles, const CAT_VolAverageResult *res,
                         int verbose, double *zscores);

/**
 * \brief Read the samples into memory, one pointer per volume.
 *
 * Every input file is read once as doubles; a 4D file gives one sample per
 * volume. After CAT_VolAverageSelect() only the selected samples are returned,
 * and files without a selected volume are not read. As everywhere in this
 * module, non-finite values arrive as 0.
 *
 * \param infiles (in)  input file names, as given to CAT_VolAverageScan()
 * \param nfiles  (in)  number of input files
 * \param res     (in)  result of CAT_VolAverageScan()
 * \param samples (out) res->n_selected pointers to res->n_voxels values each,
 *                      ordered by file and volume; they point into buffers
 * \param buffers (out) nfiles buffers with the data of the files (NULL for a
 *                      file that was not read)
 * \return 0 on success, -1 on error (a message is written to stderr and both
 *         outputs are NULL). Release with CAT_VolAverageReadFree().
 */
int CAT_VolAverageRead(char **infiles, int nfiles, const CAT_VolAverageResult *res,
                       double ***samples, double ***buffers);

/**
 * \brief Release the samples of CAT_VolAverageRead().
 *
 * \param samples (in) pointers returned by CAT_VolAverageRead(); may be NULL
 * \param buffers (in) buffers returned by CAT_VolAverageRead(); may be NULL
 * \param nfiles  (in) number of input files
 */
void CAT_VolAverageReadFree(double **samples, double **buffers, int nfiles);

/**
 * \brief Release everything held by a result and zero it.
 *
 * \param res (in/out) result of CAT_VolAverageScan(); may be zeroed already
 */
void CAT_VolAverageFree(CAT_VolAverageResult *res);

#ifdef __cplusplus
}
#endif

#endif /* _CAT_VOLUMESTATS_H_ */
