/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 */

/* Head-motion correction of a 4D series with mcflirt's command line and
 * outputs, for the options fMRIPrep's head-motion workflow uses.  This is a
 * thin CLI wrapper around the library functions in CAT_Realign.
 *
 * Input  : -in, a 4D NIfTI series; optionally -reffile, a 3D reference.
 * Output : the realigned series (always), and with -mats, -plots, -rmsabs
 *          and -rmsrel mcflirt's matrices, parameters and RMS files.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "ParseArgv.h"
#include "CAT_NiftiLib.h"
#include "CAT_ArgvCheck.h"
#include "CAT_Realign.h"

static CAT_VolRealignOpts opts; /* defaults set in main() before parsing */

static char *in_file = NULL;
static char *out_base = NULL;
static char *ref_file = NULL;
static int save_mats = 0;
static int save_plots = 0;
static int save_rms_abs = 0;
static int save_rms_rel = 0;
static int no_iscale = 0;

static ArgvInfo argTable[] = {
    {"-in", ARGV_STRING, (char *)1, (char *)&in_file,
     "4D series to motion-correct (required)."},
    {"-out", ARGV_STRING, (char *)1, (char *)&out_base,
     "Output base name; default <in>_mcf.  The other outputs append to it."},
    {"-o", ARGV_STRING, (char *)1, (char *)&out_base,
     "Same as -out."},
    {"-reffile", ARGV_STRING, (char *)1, (char *)&ref_file,
     "Separate 3D reference volume; the output is resampled into its grid."},
    {"-r", ARGV_STRING, (char *)1, (char *)&ref_file,
     "Same as -reffile."},
    {"-refvol", ARGV_INT, (char *)1, (char *)&opts.refvol,
     "Reference volume of the series (0-based) without -reffile; default nt/2."},
    {"-mats", ARGV_CONSTANT, (char *)1, (char *)&save_mats,
     "Write the transforms to <out>.mat/MAT_%04d (FSL convention)."},
    {"-plots", ARGV_CONSTANT, (char *)1, (char *)&save_plots,
     "Write the motion parameters to <out>.par."},
    {"-rmsabs", ARGV_CONSTANT, (char *)1, (char *)&save_rms_abs,
     "Write the RMS displacement from the reference to <out>_abs.rms."},
    {"-rmsrel", ARGV_CONSTANT, (char *)1, (char *)&save_rms_rel,
     "Write the RMS displacement between volumes to <out>_rel.rms."},
    {"-levels", ARGV_INT, (char *)1, (char *)&opts.reg.n_levels,
     "Gaussian pyramid levels of the registration."},
    {"-iter", ARGV_INT, (char *)1, (char *)&opts.reg.max_iter,
     "IRLS iterations per pyramid level."},
    {"-sat", ARGV_FLOAT, (char *)1, (char *)&opts.reg.sat_k,
     "Tukey saturation multiplier."},
    {"-no-iscale", ARGV_CONSTANT, (char *)1, (char *)&no_iscale,
     "Do not estimate a global intensity scale per volume."},
    {"-verbose", ARGV_CONSTANT, (char *)1, (char *)&opts.verbose,
     "Print the estimate of every volume."},
    {NULL, ARGV_END, NULL, NULL, NULL}};

static void
usage(const char *executable)
{
    CAT_VolRealignOpts d;

    CAT_VolRealignOptionsInit(&d);
    fprintf(stderr,
            "\nUsage: %s -in <series> [-out <base>] [-reffile <ref> | -refvol <n>]\n"
            "          [-mats] [-plots] [-rmsabs] [-rmsrel] [options]\n\n"
            "Head-motion correction of a 4D series with mcflirt's command line and\n"
            "outputs, for the options fMRIPrep's head-motion workflow uses:\n\n"
            "    %s -in bold.nii.gz -out bold_mcf.nii.gz -reffile boldref.nii.gz \\\n"
            "        -mats -plots -rmsabs -rmsrel\n\n"
            "Every volume is registered rigidly to the reference (robust Gauss-Newton\n"
            "with a global intensity scale), outward from the middle volume.  Outputs,\n"
            "named by appending to the -out base as mcflirt does:\n"
            "  <out>                realigned series (trilinear; .nii.gz appended to a\n"
            "                       base without .nii/.nii.gz)\n"
            "  <out>.mat/MAT_%%04d   volume -> reference, FSL convention   (-mats)\n"
            "  <out>.par            rx ry rz (rad) tx ty tz (mm)           (-plots)\n"
            "  <out>_abs.rms        RMS displacement from the reference    (-rmsabs)\n"
            "  <out>_rel.rms        RMS displacement between volumes       (-rmsrel)\n"
            "  <out>_abs_mean.rms, <out>_rel_mean.rms  their means\n\n"
            "Options beyond mcflirt's (defaults in brackets):\n"
            "  -levels <n>   pyramid levels [%d]   -iter <n>  iterations per level [%d]\n"
            "  -sat <k>      Tukey multiplier [%g]  -no-iscale  keep the intensity scale at 1\n"
            "  -verbose      print the estimate of every volume\n\n",
            executable, executable, d.reg.n_levels, d.reg.max_iter, d.reg.sat_k);
}

/* <in> without .nii/.nii.gz, plus _mcf: mcflirt's default output base */
static char *
default_base(const char *in)
{
    size_t n = strlen(in);
    char *base = (char *)malloc(n + 5);

    strcpy(base, in);
    if (n > 7 && strcmp(base + n - 7, ".nii.gz") == 0)
        base[n - 7] = '\0';
    else if (n > 4 && strcmp(base + n - 4, ".nii") == 0)
        base[n - 4] = '\0';
    strcat(base, "_mcf");
    return base;
}

int
main(int argc, char *argv[])
{
    CAT_VolRealignResult res;
    nifti_image *nii, *ref_nii = NULL, *out_nii;
    float *series = NULL, *ref = NULL, *out = NULL;
    int dims[3], ref_dims[3], nt;
    double vox[3];
    char *base, image[4096];

    CAT_VolRealignOptionsInit(&opts);
    if (ParseArgv(&argc, argv, argTable, 0) || argc != 1 || !in_file)
    {
        /* a leftover that looks like an option is a typo or a removed flag */
        cat_check_unknown_options(argc, argv, argv[0]);
        usage(argv[0]);
        exit(EXIT_FAILURE);
    }
    opts.reg.iscale = !no_iscale;
    base = out_base ? out_base : default_base(in_file);

    nii = read_nifti_float(in_file, &series, 0);
    if (!nii)
    {
        fprintf(stderr, "Error reading %s.\n", in_file);
        exit(EXIT_FAILURE);
    }
    dims[0] = nii->nx;
    dims[1] = nii->ny;
    dims[2] = nii->nz;
    nt = (nii->nt > 1) ? nii->nt : 1;

    if (ref_file)
    {
        ref_nii = read_nifti_float(ref_file, &ref, 0);
        if (!ref_nii)
        {
            fprintf(stderr, "Error reading %s.\n", ref_file);
            exit(EXIT_FAILURE);
        }
        if (ref_nii->nt > 1)
            fprintf(stderr, "Warning: %s has %d volumes; the first is the reference.\n",
                    ref_file, ref_nii->nt);
        ref_dims[0] = ref_nii->nx;
        ref_dims[1] = ref_nii->ny;
        ref_dims[2] = ref_nii->nz;
    }

    if (CAT_VolRealignEstimate(series, nii, dims, nt, ref, ref_nii, ref_dims,
                               &opts, &res) != 0)
        exit(EXIT_FAILURE);

    /* the realigned series, in the grid of the reference */
    out_nii = ref_nii ? ref_nii : nii;
    if (!ref_nii)
    {
        ref_dims[0] = dims[0];
        ref_dims[1] = dims[1];
        ref_dims[2] = dims[2];
    }
    out = (float *)malloc(sizeof(float) * (size_t)ref_dims[0] * ref_dims[1] *
                          ref_dims[2] * nt);
    if (!out || CAT_VolRealignResample(series, nii, dims, nt, out_nii, ref_dims,
                                       &res, out) != 0)
    {
        fprintf(stderr, "Error: out of memory for the realigned series.\n");
        exit(EXIT_FAILURE);
    }
    out_nii->dt = nii->dt; /* the repetition time of the series */
    vox[0] = out_nii->dx;
    vox[1] = out_nii->dy;
    vox[2] = out_nii->dz;
    CAT_VolRealignImageName(base, image, sizeof(image));
    if (!write_nifti_float_4d(image, out, DT_FLOAT32, 0.0, ref_dims, nt, vox, out_nii))
        exit(EXIT_FAILURE);

    if (CAT_VolRealignWrite(base, &res, save_mats, save_plots, save_rms_abs,
                            save_rms_rel) != 0)
        exit(EXIT_FAILURE);

    if (opts.verbose)
        printf("mean RMS displacement: %.4f mm from the reference, %.4f mm"
               " between volumes\n", res.rms_abs_mean, res.rms_rel_mean);

    CAT_VolRealignFree(&res);
    free(out);
    free(series);
    free(ref);
    nifti_image_free(nii);
    if (ref_nii)
        nifti_image_free(ref_nii);
    if (!out_base)
        free(base);
    return EXIT_SUCCESS;
}
