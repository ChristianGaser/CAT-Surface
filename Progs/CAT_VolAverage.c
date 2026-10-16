/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 *
 */

/* Mean, standard deviation and z-scores of volumes.  This is a thin CLI
 * wrapper around the library functions in CAT_VolAverage.
 *
 * Input  : one or more NIfTI files, all with the same matrix size.  A 3D file
 *          is one sample, a 4D file one sample per volume.
 * Output : whatever is asked for with -mean, -std, -zscore and -zscore-txt.
 */

#include "ParseArgv.h"
#include "CAT_NiftiLib.h"

#include <stdio.h>
#include <stdlib.h>
#include "CAT_SafeAlloc.h"
#include "CAT_ArgvCheck.h"
#include "CAT_VolAverage.h"

int verbose = 0;

char *mean_filename = NULL;
char *std_filename = NULL;
char *zscore_filename = NULL;
char *zscore_txt_filename = NULL;

static ArgvInfo argTable[] = {
    {"-mean", ARGV_STRING, (char *) 1, (char *) &mean_filename,
          "Write mean image."},
    {"-avg", ARGV_STRING, (char *) 1, (char *) &mean_filename,
          "Same as -mean."},
    {"-std", ARGV_STRING, (char *) 1, (char *) &std_filename,
          "Write standard deviation (needs at least 2 volumes)."},
    {"-zscore", ARGV_STRING, (char *) 1, (char *) &zscore_filename,
          "Write filenames and z-scores in csv-file (one row per volume)."},
    {"-zscore-txt", ARGV_STRING, (char *) 1, (char *) &zscore_txt_filename,
          "Write z-scores in txt-file (one line per volume)."},
    {"-v", ARGV_CONSTANT, (char *) 1, (char *) &verbose,
          "Be verbose."},
     {NULL, ARGV_END, NULL, NULL, NULL}
};

static void usage(const char *executable)
{
    (void) fprintf(stderr,
    "\nUsage: %s [-mean mean.nii] [-std std.nii] [-zscore zscore.csv]\n"
    "          [-zscore-txt zscore.txt] in1.nii ... inx.nii\n\n"
    "All files are inputs; the mean image is only written with -mean.\n"
    "A 4D file counts as one sample per volume, so\n"
    "     %s -mean mean.nii -std std.nii zmaps4d.nii\n"
    "gives the mean and standard deviation over its volumes.\n"
    "     %s -help\n\n", executable, executable, executable);
}

/* Main program */

int main(int argc, char *argv[])
{
    CAT_VolAverageResult res;
    char **infiles;
    int f, v, k, nfiles, need_std;
    double *zscore = NULL;
    FILE *fid;

    /* Get arguments */
    if (ParseArgv(&argc, argv, argTable, 0) || (argc < 2)) {
        usage(argv[0]);
        exit(EXIT_FAILURE);
    }

    /* a leftover that looks like an option is a typo or a removed flag;
       ParseArgv() keeps it, so it would be read as a file name below */
    if (cat_check_unknown_options(argc, argv, argv[0])) {
        usage(argv[0]);
        exit(EXIT_FAILURE);
    }

    /* every positional argument is an input, nothing is overwritten unless
       it is named after an output option */
    nfiles = argc - 1;
    infiles = &argv[1];

    need_std = (std_filename != NULL) || (zscore_filename != NULL) ||
               (zscore_txt_filename != NULL);

    if ((mean_filename == NULL) && !need_std) {
        (void) fprintf(stderr,
        "Error: No output requested. Write the mean image with -mean (the last\n"
        "file of the command line is an input and is not overwritten).\n");
        usage(argv[0]);
        exit(EXIT_FAILURE);
    }

    if (CAT_VolAverageScan(infiles, nfiles, verbose, &res))
        exit(EXIT_FAILURE);

    fprintf(stdout, "%d volume%s in %d file%s\n", res.n_samples,
            (res.n_samples == 1) ? "" : "s", nfiles, (nfiles == 1) ? "" : "s");

    if (CAT_VolAverageCompute(infiles, nfiles, need_std, verbose, &res))
        exit(EXIT_FAILURE);

    if (std_filename != NULL) {
        if (!write_nifti_double( std_filename, res.std, DT_FLOAT32, 1.0, res.dims, res.voxelsize, res.header))
            exit(EXIT_FAILURE);
    }

    if ((zscore_filename != NULL) || (zscore_txt_filename != NULL)) {
        zscore = SAFE_CALLOC(double, res.n_samples);
        if (CAT_VolAverageZscore(infiles, nfiles, &res, verbose, zscore))
            exit(EXIT_FAILURE);

        if (zscore_filename != NULL) {
            fid = SAFE_FOPEN(zscore_filename, "w");
            fprintf(fid,"filename,z-score\n");
            for (f=0, k=0; f<nfiles; f++) {
                /* the volumes of a 4D file are numbered from 1 */
                for (v=0; v<res.n_volumes[f]; v++, k++) {
                    if (res.n_volumes[f] > 1)
                        fprintf(fid,"%s[%d],%g\n",infiles[f],v+1,zscore[k]);
                    else
                        fprintf(fid,"%s,%g\n",infiles[f],zscore[k]);
                }
            }
            fclose(fid);
        }
        if (zscore_txt_filename != NULL) {
            fid = SAFE_FOPEN(zscore_txt_filename, "w");
            for (k=0; k<res.n_samples; k++)
                fprintf(fid,"%g\n",zscore[k]);
            fclose(fid);
        }

        free(zscore);
    }

    if (mean_filename != NULL) {
        if (!write_nifti_double( mean_filename, res.mean, DT_FLOAT32, 1.0, res.dims, res.voxelsize, res.header))
            exit(EXIT_FAILURE);
    }

    CAT_VolAverageFree(&res);

    return(EXIT_SUCCESS);

}
