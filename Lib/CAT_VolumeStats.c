/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 *
 */

/**
 * \file CAT_VolumeStats.c
 * \brief Voxel-wise mean, standard deviation and z-scores of NIfTI volumes.
 *
 * Every volume is one sample: a 3D file contributes one and a 4D file one per
 * volume. See CAT_VolumeStats.h.
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "CAT_SafeAlloc.h"
#include "CAT_VolumeStats.h"

/**
 * \brief Number of volumes of a NIfTI header (1 for a 3D image).
 *
 * \param hdr (in) NIfTI header
 * \return number of volumes, at least 1
 */
static int
volumes_of(const nifti_image *hdr)
{
    return (hdr->nt > 1) ? hdr->nt : 1;
}

/**
 * \brief Read one input file as doubles and check it against the scan.
 *
 * \param infile    (in) file name
 * \param n_voxels  (in) voxels of one volume found by CAT_VolAverageScan()
 * \param n_volumes (in) volumes found by CAT_VolAverageScan()
 * \return n_volumes * n_voxels values, volume after volume (free with free());
 *         NULL on error
 */
static double *
read_volumes(const char *infile, size_t n_voxels, int n_volumes)
{
    nifti_image *nii;
    double *data = NULL;

    nii = read_nifti_double(infile, &data, 0);
    if (!nii)
    {
        fprintf(stderr, "Error reading %s.\n", infile);
        return NULL;
    }

    if (nii->nvox != n_voxels || volumes_of(nii) != n_volumes)
    {
        fprintf(stderr, "Error: %s changed since it was scanned.\n", infile);
        free(data);
        nifti_image_free(nii);
        return NULL;
    }

    nifti_image_free(nii);
    return data;
}

/**
 * \brief Print the progress line of a pass over the input files.
 *
 * \param pass      (in) "run1" or "run2"
 * \param f         (in) index of the file
 * \param nfiles    (in) number of files
 * \param infile    (in) name of the file
 * \param n_volumes (in) volumes of the file
 */
static void
print_progress(const char *pass, int f, int nfiles, const char *infile,
               int n_volumes)
{
    if (n_volumes > 1)
        fprintf(stdout, "%s %5d/%d:\t%s (%d volumes)\n", pass, f + 1, nfiles,
                infile, n_volumes);
    else
        fprintf(stdout, "%s %5d/%d:\t%s\n", pass, f + 1, nfiles, infile);
}

/**
 * \brief Turn the sum of squares into the standard deviation, in place.
 *
 * sqrt((sum_x2 - n*mean^2) / (n-1)). The radicand is clamped at 0: for a voxel
 * that is constant across the samples the two terms are equal up to rounding,
 * and when the difference comes out negative sqrt() would return NaN.
 *
 * \param sum_squares (in/out) sum of squares of the samples; the std on return
 * \param mean        (in)     mean of the samples
 * \param n_voxels    (in)     number of voxels
 * \param n           (in)     number of samples, at least 2
 */
static void
std_from_sums(double *sum_squares, const double *mean, size_t n_voxels, double n)
{
    double var;
    size_t j;

    for (j = 0; j < n_voxels; j++)
    {
        var = sum_squares[j] - n * mean[j] * mean[j];
        if (var < 0.0)
            var = 0.0;
        sum_squares[j] = sqrt(1.0 / (n - 1.0) * var);
    }
}

/**
 * \brief Mean value above which a voxel is taken into the z-score.
 *
 * A quarter of the average of the mean image over the voxels brighter than an
 * eighth of its non-zero average, which leaves out the background.
 *
 * \param mean     (in) mean image
 * \param n_voxels (in) number of voxels
 * \return threshold for the mean image
 */
static double
zscore_threshold(const double *mean, size_t n_voxels)
{
    double mean_nz = 0.0, mean_ab = 0.0, thresh8;
    size_t j, cnt_nz = 0, cnt_ab = 0;

    for (j = 0; j < n_voxels; j++)
    {
        if (mean[j] != 0.0)
        {
            mean_nz += mean[j];
            cnt_nz++;
        }
    }
    if (cnt_nz > 0)
        mean_nz /= (double)cnt_nz;
    thresh8 = mean_nz / 8.0;

    for (j = 0; j < n_voxels; j++)
    {
        if (mean[j] > thresh8)
        {
            mean_ab += mean[j];
            cnt_ab++;
        }
    }
    if (cnt_ab > 0)
        mean_ab /= (double)cnt_ab;

    return 0.25 * mean_ab;
}

/**
 * \brief Z-score of one volume: the fourth root of the mean of z^4.
 *
 * \param volume   (in) voxel values of the volume
 * \param mean     (in) mean image
 * \param std      (in) standard deviation image
 * \param n_voxels (in) number of voxels
 * \param thresh   (in) voxels with a mean up to this value are left out
 * \return z-score, 0 if no voxel qualifies
 */
static double
volume_zscore(const double *volume, const double *mean, const double *std,
              size_t n_voxels, double thresh)
{
    double sum = 0.0, z;
    size_t j, n = 0;

    for (j = 0; j < n_voxels; j++)
    {
        if (std[j] == 0.0 || isnan(std[j]) || mean[j] <= thresh)
            continue;
        z = (volume[j] - mean[j]) / std[j];
        sum += z * z * z * z;
        n++;
    }

    return (n > 0) ? pow(sum / (double)n, 0.25) : 0.0;
}

/**
 * \brief Release everything held by a result and zero it.
 *
 * \param res (in/out) result of CAT_VolAverageScan(); may be zeroed already
 */
void
CAT_VolAverageFree(CAT_VolAverageResult *res)
{
    if (!res)
        return;

    free(res->n_volumes);
    free(res->mean);
    free(res->std);
    if (res->header)
        nifti_image_free(res->header);

    memset(res, 0, sizeof(*res));
}

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
int
CAT_VolAverageScan(char **infiles, int nfiles, int verbose,
                   CAT_VolAverageResult *res)
{
    nifti_image *hdr;
    int f;

    memset(res, 0, sizeof(*res));

    if (!infiles || nfiles < 1)
    {
        fprintf(stderr, "CAT_VolAverageScan: No input files specified.\n");
        return -1;
    }

    res->n_volumes = SAFE_CALLOC(int, nfiles);

    for (f = 0; f < nfiles; f++)
    {
        hdr = nifti_image_read(infiles[f], 0);
        if (!hdr)
        {
            fprintf(stderr, "Error reading %s.\n", infiles[f]);
            CAT_VolAverageFree(res);
            return -1;
        }

        if (hdr->nu > 1 || hdr->nv > 1 || hdr->nw > 1)
        {
            fprintf(stderr, "Error: %s has more than four dimensions.\n",
                    infiles[f]);
            nifti_image_free(hdr);
            CAT_VolAverageFree(res);
            return -1;
        }

        res->n_volumes[f] = volumes_of(hdr);
        res->n_samples += res->n_volumes[f];

        if (f == 0)
        {
            res->header = hdr;
            res->dims[0] = hdr->nx;
            res->dims[1] = hdr->ny;
            res->dims[2] = hdr->nz;
            res->voxelsize[0] = hdr->dx;
            res->voxelsize[1] = hdr->dy;
            res->voxelsize[2] = hdr->dz;
            res->n_voxels = (size_t)hdr->nx * (size_t)hdr->ny * (size_t)hdr->nz;
        }
        else
        {
            if (hdr->nx != res->dims[0] || hdr->ny != res->dims[1] ||
                hdr->nz != res->dims[2])
            {
                fprintf(stderr,
                        "Error: %s has %dx%dx%d voxels, but %s has %dx%dx%d.\n",
                        infiles[f], hdr->nx, hdr->ny, hdr->nz, infiles[0],
                        res->dims[0], res->dims[1], res->dims[2]);
                nifti_image_free(hdr);
                CAT_VolAverageFree(res);
                return -1;
            }
            if (verbose)
                equal_image_dimensions(res->header, hdr);
            nifti_image_free(hdr);
        }
    }

    return 0;
}

/**
 * \brief Compute the voxel-wise mean and, optionally, standard deviation.
 *
 * Every volume is one sample, the standard deviation uses n - 1. It is clamped
 * at 0 where the one-pass formula comes out slightly negative from rounding, as
 * it does for voxels that are constant across the samples.
 *
 * \param infiles  (in)     input file names, as given to CAT_VolAverageScan()
 * \param nfiles   (in)     number of input files
 * \param want_std (in)     non-zero to compute the standard deviation as well;
 *                          needs at least two samples
 * \param verbose  (in)     non-zero to print progress to stdout
 * \param res      (in/out) result of CAT_VolAverageScan(); mean (and std) are set
 * \return 0 on success, -1 on error (a message is written to stderr)
 */
int
CAT_VolAverageCompute(char **infiles, int nfiles, int want_std, int verbose,
                      CAT_VolAverageResult *res)
{
    double *data, *sum_squares = NULL, n;
    const double *volume;
    size_t j;
    int f, v;

    if (!infiles || !res || !res->n_volumes || nfiles < 1)
    {
        fprintf(stderr, "CAT_VolAverageCompute: Run CAT_VolAverageScan() first.\n");
        return -1;
    }

    if (want_std && res->n_samples < 2)
    {
        fprintf(stderr, "Error: The standard deviation needs at least two "
                        "volumes, but only %d found.\n", res->n_samples);
        return -1;
    }

    free(res->mean);
    free(res->std);
    res->std = NULL;
    res->mean = SAFE_CALLOC(double, res->n_voxels);
    if (want_std)
        sum_squares = SAFE_CALLOC(double, res->n_voxels);

    n = (double)res->n_samples;

    for (f = 0; f < nfiles; f++)
    {
        if (verbose)
            print_progress("run1", f, nfiles, infiles[f], res->n_volumes[f]);

        data = read_volumes(infiles[f], res->n_voxels, res->n_volumes[f]);
        if (!data)
        {
            free(sum_squares);
            free(res->mean);
            res->mean = NULL;
            return -1;
        }

        for (v = 0; v < res->n_volumes[f]; v++)
        {
            volume = data + (size_t)v * res->n_voxels;
            for (j = 0; j < res->n_voxels; j++)
                res->mean[j] += volume[j] / n;
            if (sum_squares)
            {
                for (j = 0; j < res->n_voxels; j++)
                    sum_squares[j] += volume[j] * volume[j];
            }
        }
        free(data);
    }

    if (sum_squares)
    {
        std_from_sums(sum_squares, res->mean, res->n_voxels, n);
        res->std = sum_squares;
    }

    return 0;
}

/**
 * \brief Compute one z-score per sample.
 *
 * The z-score of a sample is the fourth root of the mean of z^4, with
 * z = (x - mean) / std, over the voxels where std is positive and the mean
 * exceeds a quarter of the average of the mean image over the voxels brighter
 * than an eighth of its non-zero average (the background is left out).
 * Samples are ordered by file and, within a 4D file, by volume.
 *
 * \param infiles (in)  input file names, as given to CAT_VolAverageScan()
 * \param nfiles  (in)  number of input files
 * \param res     (in)  result of CAT_VolAverageCompute() with the standard deviation
 * \param verbose (in)  non-zero to print progress to stdout
 * \param zscores (out) res->n_samples z-scores, allocated by the caller
 * \return 0 on success, -1 on error (a message is written to stderr)
 */
int
CAT_VolAverageZscore(char **infiles, int nfiles, const CAT_VolAverageResult *res,
                     int verbose, double *zscores)
{
    double *data, thresh;
    int f, v, k = 0;

    if (!infiles || !res || !res->mean || !res->std || !zscores || nfiles < 1)
    {
        fprintf(stderr, "CAT_VolAverageZscore: Run CAT_VolAverageCompute() with "
                        "the standard deviation first.\n");
        return -1;
    }

    thresh = zscore_threshold(res->mean, res->n_voxels);

    for (f = 0; f < nfiles; f++)
    {
        if (verbose)
            print_progress("run2", f, nfiles, infiles[f], res->n_volumes[f]);

        data = read_volumes(infiles[f], res->n_voxels, res->n_volumes[f]);
        if (!data)
            return -1;

        for (v = 0; v < res->n_volumes[f]; v++)
            zscores[k++] = volume_zscore(data + (size_t)v * res->n_voxels,
                                         res->mean, res->std, res->n_voxels,
                                         thresh);
        free(data);
    }

    return 0;
}
