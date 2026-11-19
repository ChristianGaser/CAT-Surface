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

#include <ctype.h>
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
    free(res->selected);
    if (res->header)
        nifti_image_free(res->header);

    memset(res, 0, sizeof(*res));
}

/**
 * \brief Whether a sample is used: all are without a selection.
 *
 * \param res (in) result of CAT_VolAverageScan()
 * \param k   (in) index of the sample
 * \return non-zero if the sample is used
 */
static int
sample_used(const CAT_VolAverageResult *res, int k)
{
    return !res->selected || res->selected[k];
}

/**
 * \brief Whether a file has a sample that is used.
 *
 * \param res       (in) result of CAT_VolAverageScan()
 * \param first     (in) index of the first sample of the file
 * \param n_volumes (in) volumes of the file
 * \return non-zero if at least one of its volumes is used
 */
static int
file_used(const CAT_VolAverageResult *res, int first, int n_volumes)
{
    int v;

    for (v = 0; v < n_volumes; v++)
        if (sample_used(res, first + v))
            return 1;
    return 0;
}

/**
 * \brief Strip the white space at both ends of a string, in place.
 *
 * \param s (in/out) string; a terminator is written after its last character
 * \return pointer to its first character that is not white space
 */
static char *
strip(char *s)
{
    char *end;

    while (*s && isspace((unsigned char)*s))
        s++;
    end = s + strlen(s);
    while (end > s && isspace((unsigned char)end[-1]))
        end--;
    *end = '\0';
    return s;
}

/**
 * \brief Value of a text that is a number as a whole.
 *
 * \param text  (in)  text without surrounding white space
 * \param value (out) its value, if it is a number
 * \return 1 if the whole text is a number, 0 otherwise
 */
static int
as_number(const char *text, double *value)
{
    char *end;

    if (*text == '\0')
        return 0;
    *value = strtod(text, &end);
    return *end == '\0';
}

/**
 * \brief Larger copy of a block, for a list that grows.
 *
 * \param old      (in) block to copy and free; may be NULL
 * \param old_size (in) its size in bytes
 * \param new_size (in) size of the new block in bytes
 * \return the new block (never NULL)
 */
static void *
grown(void *old, size_t old_size, size_t new_size)
{
    char *block = SAFE_MALLOC(char, new_size);

    if (old)
    {
        memcpy(block, old, old_size);
        free(old);
    }
    return block;
}

/**
 * \brief Release a list of labels.
 *
 * \param labels (in) labels of read_labels(); may be NULL
 * \param n      (in) number of labels
 */
static void
free_labels(char **labels, int n)
{
    int i;

    if (!labels)
        return;
    for (i = 0; i < n; i++)
        free(labels[i]);
    free(labels);
}

/**
 * \brief Read the labels of a text file, one per line.
 *
 * Lines of any length and with any line ending; empty lines are skipped and
 * the white space around a label is dropped.
 *
 * \param file   (in)  file name
 * \param labels (out) the labels; release with free_labels()
 * \return number of labels, -1 if the file cannot be opened
 */
static int
read_labels(const char *file, char ***labels)
{
    FILE *fp = fopen(file, "r");
    char **list = NULL, *line = NULL, *label;
    size_t len = 0, cap = 0;
    int n = 0, n_cap = 0, c;

    *labels = NULL;
    if (!fp)
        return -1;

    do
    {
        c = fgetc(fp);
        if (c != EOF && c != '\n')
        {
            if (len + 2 > cap)
            {
                line = grown(line, len, cap ? 2 * cap : 64);
                cap = cap ? 2 * cap : 64;
            }
            line[len++] = (char)c;
            continue;
        }
        if (len > 0)
        {
            line[len] = '\0';
            label = strip(line);
            if (*label)
            {
                if (n == n_cap)
                {
                    list = grown(list, sizeof(char *) * (size_t)n,
                                 sizeof(char *) * (size_t)(n_cap ? 2 * n_cap : 256));
                    n_cap = n_cap ? 2 * n_cap : 256;
                }
                list[n] = SAFE_MALLOC(char, strlen(label) + 1);
                strcpy(list[n++], label);
            }
            len = 0;
        }
    } while (c != EOF);

    free(line);
    fclose(fp);
    *labels = list;
    return n;
}

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
int
CAT_VolAverageSelect(const char *selection, CAT_VolAverageResult *res)
{
    char *spec, *file, *values = NULL, *sep, **labels = NULL, **wanted = NULL;
    unsigned char *selected = NULL;
    double *label_value = NULL, *wanted_value = NULL;
    int n_labels = 0, n_wanted = 0, numeric = 1, n_selected = 0, shown = 0;
    int i, j, rc = -1;
    FILE *fp;

    if (!selection || !res || !res->n_volumes || res->n_samples < 1)
    {
        fprintf(stderr, "CAT_VolAverageSelect: Run CAT_VolAverageScan() first.\n");
        return -1;
    }

    /* "FILE=VALUES", unless a file of the whole name exists */
    spec = SAFE_MALLOC(char, strlen(selection) + 1);
    strcpy(spec, selection);
    file = spec;
    fp = fopen(spec, "r");
    if (fp)
        fclose(fp);
    else if ((sep = strrchr(spec, '=')) != NULL)
    {
        *sep = '\0';
        values = sep + 1;
    }

    n_labels = read_labels(file, &labels);
    if (n_labels < 0)
    {
        fprintf(stderr, "Error: Cannot read the labels of %s.\n", file);
        n_labels = 0;
        goto done;
    }
    if (n_labels != res->n_samples)
    {
        fprintf(stderr, "Error: %s has %d label%s, but there %s %d volume%s.\n",
                file, n_labels, (n_labels == 1) ? "" : "s",
                (res->n_samples == 1) ? "is" : "are", res->n_samples,
                (res->n_samples == 1) ? "" : "s");
        goto done;
    }

    label_value = SAFE_CALLOC(double, n_labels);
    for (i = 0; i < n_labels; i++)
        numeric = as_number(labels[i], &label_value[i]) && numeric;

    if (values)
    {
        /* the values, split at the commas in place */
        n_wanted = 1;
        for (sep = values; *sep; sep++)
            n_wanted += (*sep == ',');
        wanted = SAFE_CALLOC(char *, n_wanted);
        wanted_value = SAFE_CALLOC(double, n_wanted);
        for (j = 0, sep = values; j < n_wanted; j++)
        {
            wanted[j] = sep;
            sep = strchr(sep, ',');
            if (sep)
                *sep++ = '\0';
            wanted[j] = strip(wanted[j]);
            if (*wanted[j] == '\0')
            {
                fprintf(stderr, "Error: Empty value in the selection %s.\n", selection);
                goto done;
            }
            numeric = as_number(wanted[j], &wanted_value[j]) && numeric;
        }
    }
    else if (!numeric)
    {
        fprintf(stderr, "Error: The labels of %s are not numbers: give the values "
                        "to select as %s=VALUE[,VALUE...].\n", file, file);
        goto done;
    }

    selected = SAFE_CALLOC(unsigned char, n_labels);
    for (i = 0; i < n_labels; i++)
    {
        if (!values)
            selected[i] = (label_value[i] != 0.0);
        for (j = 0; j < n_wanted && !selected[i]; j++)
            selected[i] = numeric ? (label_value[i] == wanted_value[j])
                                  : (strcmp(labels[i], wanted[j]) == 0);
        n_selected += selected[i];
    }

    if (n_selected == 0)
    {
        fprintf(stderr, "Error: The selection %s matches none of the %d volumes "
                        "(labels:", selection, n_labels);
        for (i = 0; i < n_labels && shown < 8; i++)
        {
            for (j = 0; j < i && strcmp(labels[i], labels[j]) != 0; j++)
                ;
            if (j == i)
                fprintf(stderr, "%s %s", shown++ ? "," : "", labels[i]);
        }
        fprintf(stderr, "%s).\n", (i < n_labels) ? ", ..." : "");
        goto done;
    }

    /* a mean of other samples does not belong to this selection */
    free(res->selected);
    free(res->mean);
    free(res->std);
    res->mean = res->std = NULL;
    res->selected = selected;
    res->n_selected = n_selected;
    selected = NULL;
    rc = 0;

done:
    free(selected);
    free(wanted);
    free(wanted_value);
    free(label_value);
    free_labels(labels, n_labels);
    free(spec);
    return rc;
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

    res->n_selected = res->n_samples;

    return 0;
}

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
int
CAT_VolAverageCompute(char **infiles, int nfiles, int want_std, int verbose,
                      CAT_VolAverageResult *res)
{
    double *data, *sum_squares = NULL, n;
    const double *volume;
    size_t j;
    int f, v, k;

    if (!infiles || !res || !res->n_volumes || nfiles < 1)
    {
        fprintf(stderr, "CAT_VolAverageCompute: Run CAT_VolAverageScan() first.\n");
        return -1;
    }

    if (want_std && res->n_selected < 2)
    {
        fprintf(stderr, "Error: The standard deviation needs at least two "
                        "volumes, but only %d %s.\n", res->n_selected,
                res->selected ? "selected" : "found");
        return -1;
    }

    free(res->mean);
    free(res->std);
    res->std = NULL;
    res->mean = SAFE_CALLOC(double, res->n_voxels);
    if (want_std)
        sum_squares = SAFE_CALLOC(double, res->n_voxels);

    n = (double)res->n_selected;

    /* k is the first sample of file f */
    for (f = 0, k = 0; f < nfiles; k += res->n_volumes[f], f++)
    {
        if (!file_used(res, k, res->n_volumes[f]))
            continue;

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
            if (!sample_used(res, k + v))
                continue;
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
        if (!file_used(res, k, res->n_volumes[f]))
        {
            for (v = 0; v < res->n_volumes[f]; v++)
                zscores[k++] = NAN;
            continue;
        }

        if (verbose)
            print_progress("run2", f, nfiles, infiles[f], res->n_volumes[f]);

        data = read_volumes(infiles[f], res->n_voxels, res->n_volumes[f]);
        if (!data)
            return -1;

        for (v = 0; v < res->n_volumes[f]; v++, k++)
            zscores[k] = sample_used(res, k)
                             ? volume_zscore(data + (size_t)v * res->n_voxels,
                                             res->mean, res->std, res->n_voxels,
                                             thresh)
                             : NAN;
        free(data);
    }

    return 0;
}

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
int
CAT_VolAverageRead(char **infiles, int nfiles, const CAT_VolAverageResult *res,
                   double ***samples, double ***buffers)
{
    double **list, **data;
    int f, v, k, n = 0;

    if (samples)
        *samples = NULL;
    if (buffers)
        *buffers = NULL;

    if (!infiles || !res || !res->n_volumes || !samples || !buffers ||
        nfiles < 1 || res->n_selected < 1)
    {
        fprintf(stderr, "CAT_VolAverageRead: Run CAT_VolAverageScan() first.\n");
        return -1;
    }

    list = SAFE_CALLOC(double *, res->n_selected);
    data = SAFE_CALLOC(double *, nfiles);

    /* k is the first sample of file f */
    for (f = 0, k = 0; f < nfiles; k += res->n_volumes[f], f++)
    {
        if (!file_used(res, k, res->n_volumes[f]))
            continue;

        data[f] = read_volumes(infiles[f], res->n_voxels, res->n_volumes[f]);
        if (!data[f])
        {
            CAT_VolAverageReadFree(list, data, nfiles);
            return -1;
        }

        for (v = 0; v < res->n_volumes[f]; v++)
            if (sample_used(res, k + v))
                list[n++] = data[f] + (size_t)v * res->n_voxels;
    }

    *samples = list;
    *buffers = data;
    return 0;
}

/**
 * \brief Release the samples of CAT_VolAverageRead().
 *
 * \param samples (in) pointers returned by CAT_VolAverageRead(); may be NULL
 * \param buffers (in) buffers returned by CAT_VolAverageRead(); may be NULL
 * \param nfiles  (in) number of input files
 */
void
CAT_VolAverageReadFree(double **samples, double **buffers, int nfiles)
{
    int f;

    if (buffers)
    {
        for (f = 0; f < nfiles; f++)
            free(buffers[f]);
        free(buffers);
    }
    free(samples);
}
