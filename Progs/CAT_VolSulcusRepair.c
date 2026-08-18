/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 *
 */

#include <float.h>
#include <stdlib.h>
#include <string.h>

#if !defined(_WIN32) && !defined(_WIN64)
#include <libgen.h>
#endif

#include <bicpl.h>
#include <ParseArgv.h>
#include "CAT_NiftiLib.h"
#include "CAT_Vol.h"
#include "CAT_Math.h"
#include "CAT_Sheetness.h"
#include "CAT_SulcusRepair.h"

int recover_csf = 1;
int reconnect_gyri = 1;
int refine_pve = 0;
int verbose = 0;

double sheet_sigma_min = 0.3;
double sheet_sigma_max = 1.0;
int sheet_n_scales = 3;

double csf_min_dist = 1.5;
double csf_min_wmdist = 0.75;
double csf_thresh = 0.3;
double csf_strength = 0.8;

double wm_thresh = 0.3;
double wm_strength = 0.8;
int wm_max_gap = 2;

double band_min_dist = 1.5;
int band_window = 4;
double band_strength = 0.7;

char *out_sheetness = NULL;

static ArgvInfo argTable[] = {
    {"-verbose", ARGV_CONSTANT, (char *)1, (char *)&verbose,
     "Report what each step changed."},

    {"-no-recover-csf", ARGV_CONSTANT, (char *)0, (char *)&recover_csf,
     "Skip the sulcal CSF recovery."},

    {"-no-reconnect-gyri", ARGV_CONSTANT, (char *)0, (char *)&reconnect_gyri,
     "Skip the gyral WM reconnection."},

    {"-refine-pve", ARGV_CONSTANT, (char *)1, (char *)&refine_pve,
     "Also re-estimate partial volume in the narrow band where the label map\n\
    reports no CSF but the intensity image still shows a dip across the sheet.\n\
    Off by default because it is the most aggressive of the three steps."},

    {"-sheet-sigma-min", ARGV_FLOAT, (char *)1, (char *)&sheet_sigma_min,
     "Smallest Gaussian scale of the sheetness filter, in mm (default 0.3)."},

    {"-sheet-sigma-max", ARGV_FLOAT, (char *)1, (char *)&sheet_sigma_max,
     "Largest Gaussian scale of the sheetness filter, in mm (default 1.0).\n\
    The range should bracket the structures being looked for: a sulcal CSF\n\
    sheet or a gyral WM blade that survives at all is one to three voxels\n\
    thick at 0.5 mm. Larger scales start responding to the cortical ribbon\n\
    itself, which is not what this filter is for."},

    {"-sheet-scales", ARGV_INT, (char *)1, (char *)&sheet_n_scales,
     "Number of log-spaced scales between the two sigmas (default 3)."},

    {"-csf-min-dist", ARGV_FLOAT, (char *)1, (char *)&csf_min_dist,
     "Act on sulcal CSF only where the label map finds no CSF within this many\n\
    mm (default 1.5). This is what restricts the correction to places where\n\
    the classifier really did commit to solid gray matter."},

    {"-csf-min-wmdist", ARGV_FLOAT, (char *)1, (char *)&csf_min_wmdist,
     "Never carve closer than this many mm to the WM boundary (default 0.75).\n\
    This is the T > 1 guard of ACE and is what keeps thin cortex from being\n\
    eaten from the inside."},

    {"-csf-thresh", ARGV_FLOAT, (char *)1, (char *)&csf_thresh,
     "Dark-sheet response below this is ignored (default 0.3)."},

    {"-csf-strength", ARGV_FLOAT, (char *)1, (char *)&csf_strength,
     "How far a label may be blended towards the intensity-implied value,\n\
    0..1 (default 0.8). The correction is one-sided: it can only lower a\n\
    label, and never below what the intensity itself supports."},

    {"-wm-thresh", ARGV_FLOAT, (char *)1, (char *)&wm_thresh,
     "Bright-sheet response below this is ignored (default 0.3)."},

    {"-wm-strength", ARGV_FLOAT, (char *)1, (char *)&wm_strength,
     "How far a label may be blended towards WM, 0..1 (default 0.8)."},

    {"-wm-max-gap", ARGV_INT, (char *)1, (char *)&wm_max_gap,
     "Largest gap half-width, in voxels, that the two-sided test will bridge\n\
    (default 2). A candidate is only accepted when WM is present on two\n\
    opposite sides within this distance, which is why the step cannot fire\n\
    inside a sulcus."},

    {"-band-min-dist", ARGV_FLOAT, (char *)1, (char *)&band_min_dist,
     "Narrow-band refit acts only outside this distance to detected CSF, in mm\n\
    (default 1.5)."},

    {"-band-window", ARGV_INT, (char *)1, (char *)&band_window,
     "Half-width in voxels of the window used to re-estimate local class\n\
    levels (default 4). Local levels make the refit immune to residual bias."},

    {"-band-strength", ARGV_FLOAT, (char *)1, (char *)&band_strength,
     "Blend weight towards the refitted value, 0..1 (default 0.7)."},

    {"-sheetness", ARGV_STRING, (char *)1, (char *)&out_sheetness,
     "Also write the sheetness map of the last executed step to this file.\n\
    Useful for checking the scale range and thresholds on a new protocol, and\n\
    it can be fed to CAT_VolAmap -mrf-aniso as a precomputed field."},

    {NULL, ARGV_END, NULL, NULL, NULL}};

static void usage(char *executable)
{
    char *usage_str = "\n\
Usage: %s [options] t1.nii label.nii [output_label.nii]\n\
\n\
    Anatomy-aware repair of a PVE label map, to be run before PBT.\n\
\n\
    Three failure modes of a tissue classifier break central-surface\n\
    extraction, and all three are failures of evidence rather than of\n\
    smoothness, which is why no local filter fixes them:\n\
\n\
      1. Glued sulci. Two banks of a tight sulcus end up as one thick GM band\n\
         because no CSF was detected between them. Typical in the occipital\n\
         midline, where cortex is thin and contrast is poorest.\n\
      2. Broken gyri. A thin gyral WM blade is interrupted by a small\n\
         missegmentation, which corrupts the WM distance map and with it the\n\
         thickness and the central surface along the whole blade.\n\
      3. Residual partial-volume error in the band where (1) happens: the\n\
         label map says there is no CSF anywhere near, while the T1 still\n\
         shows an intensity dip across the sulcus.\n\
\n\
    Regularization cannot repair any of these, because it cannot create\n\
    evidence -- it can only redistribute what the classifier already committed\n\
    to. Each step here goes back to the intensity image instead and recovers\n\
    evidence the classifier discarded, using a Hessian sheetness filter as the\n\
    shape prior. A sheetness filter keeps thin sheets and ignores blobs, and it\n\
    shrinks nothing, because unlike a median or an MRF it makes no statement\n\
    about boundary length.\n\
\n\
    Every operation is one-sided and gated on several independent pieces of\n\
    evidence, so none of them can produce the failure mode it is meant to fix:\n\
    the CSF recovery can only lower a label and never within one voxel of the\n\
    WM boundary; the gyral reconnection only fires where WM is present on two\n\
    opposite sides; the narrow-band refit only touches voxels with a genuine\n\
    local minimum across the sheet.\n\
\n\
    The label map is a PVE image in [0..3] with CSF = 1, GM = 2, WM = 3, the\n\
    same convention CAT_VolThicknessPbt expects. The intensity image should be\n\
    bias corrected but may be in arbitrary units; it is rescaled internally.\n\
\n\
    Options:\n\
    -no-recover-csf            Skip step 1.\n\
    -no-reconnect-gyri         Skip step 2.\n\
    -refine-pve                Also run step 3 (off by default).\n\
    -sheetness <file>          Write the sheetness map for inspection.\n\
    -verbose                   Report what each step changed.\n\
\n\
    Example:\n\
      %s -verbose t1_corr.nii label.nii label_repaired.nii\n\
      CAT_VolThicknessPbt label_repaired.nii gmt.nii ppm.nii\n\
\n\
    References:\n\
      Han et al., Proc SPIE Med Imag 4322:194-203, 2001 (ACE).\n\
      Han et al., NeuroImage 23(3):997-1012, 2004 (CRUISE).\n\
      Kim et al., NeuroImage 27(1):210-221, 2005 (CLASP, CSF skeleton).\n\
      Descoteaux et al., Med Image Anal 10(4):638-651, 2006 (sheetness).\n\n";

    fprintf(stderr, usage_str, executable, executable);
}

int main(int argc, char *argv[])
{
    char out_label[1024];
    int dims[3], nvox, rc;
    float *t1 = NULL, *label = NULL, *sheetness = NULL;
    double voxelsize[3], slope;
    CAT_SulcusRepairOpts opts;

    initialize_argument_processing(argc, argv);

    if (ParseArgv(&argc, argv, argTable, 0) || (argc < 3))
    {
        usage(argv[0]);
        fprintf(stderr, "     %s -help\n\n", argv[0]);
        exit(EXIT_FAILURE);
    }

    char *t1file = argv[1];
    char *labelfile = argv[2];

    if (argc >= 4)
        (void)sprintf(out_label, "%s", argv[3]);
    else
    {
#if !defined(_WIN32) && !defined(_WIN64)
        (void)sprintf(out_label, "%s/repaired_%s", dirname(labelfile), basename(labelfile));
#else
        fprintf(stderr, "\nUsage: %s t1.nii label.nii output_label.nii\n\n", argv[0]);
        return (1);
#endif
    }

    nifti_image *t1_ptr = read_nifti_float(t1file, &t1, 0);
    if (!t1_ptr)
    {
        fprintf(stderr, "Error reading %s.\n", t1file);
        return (EXIT_FAILURE);
    }

    nifti_image *label_ptr = read_nifti_float(labelfile, &label, 0);
    if (!label_ptr)
    {
        fprintf(stderr, "Error reading %s.\n", labelfile);
        return (EXIT_FAILURE);
    }

    if (t1_ptr->nvox != label_ptr->nvox ||
        t1_ptr->nx != label_ptr->nx ||
        t1_ptr->ny != label_ptr->ny ||
        t1_ptr->nz != label_ptr->nz)
    {
        fprintf(stderr, "Intensity and label image must have the same dimensions.\n");
        return (EXIT_FAILURE);
    }

    nvox = label_ptr->nvox;
    dims[0] = label_ptr->nx;
    dims[1] = label_ptr->ny;
    dims[2] = label_ptr->nz;
    voxelsize[0] = label_ptr->dx;
    voxelsize[1] = label_ptr->dy;
    voxelsize[2] = label_ptr->dz;

    CAT_SulcusRepairOptionsInit(&opts);
    opts.sheet_sigma_min = sheet_sigma_min;
    opts.sheet_sigma_max = sheet_sigma_max;
    opts.sheet_n_scales = sheet_n_scales;
    opts.csf_min_dist = csf_min_dist;
    opts.csf_min_wmdist = csf_min_wmdist;
    opts.csf_thresh = csf_thresh;
    opts.csf_strength = csf_strength;
    opts.wm_thresh = wm_thresh;
    opts.wm_strength = wm_strength;
    opts.wm_max_gap = wm_max_gap;
    opts.band_min_dist = band_min_dist;
    opts.band_window = band_window;
    opts.band_strength = band_strength;
    opts.verbose = verbose;

    if (out_sheetness)
    {
        sheetness = (float *)malloc(sizeof(float) * nvox);
        if (!sheetness)
        {
            fprintf(stderr, "Memory allocation error\n");
            exit(EXIT_FAILURE);
        }
    }

    /* Order matters: the CSF recovery changes the label map that the gyral
       reconnection reads, and both change the CSF distance the narrow-band
       refit is gated on. */
    if (recover_csf)
    {
        rc = CAT_VolRecoverSulcalCSF(t1, label, sheetness, dims, voxelsize, &opts);
        if (rc != 0)
        {
            fprintf(stderr, "Sulcal CSF recovery failed (%d).\n", rc);
            exit(EXIT_FAILURE);
        }
    }

    if (reconnect_gyri)
    {
        rc = CAT_VolReconnectGyri(t1, label, sheetness, dims, voxelsize, &opts);
        if (rc != 0)
        {
            fprintf(stderr, "Gyral reconnection failed (%d).\n", rc);
            exit(EXIT_FAILURE);
        }
    }

    if (refine_pve)
    {
        rc = CAT_VolRefinePveNarrowBand(t1, label, dims, voxelsize, &opts);
        if (rc != 0)
        {
            fprintf(stderr, "Narrow-band PVE refit failed (%d).\n", rc);
            exit(EXIT_FAILURE);
        }
    }

    slope = 1.0;
    if (!write_nifti_float(out_label, label, DT_FLOAT32, slope, dims, voxelsize,
                           label_ptr))
        exit(EXIT_FAILURE);

    if (out_sheetness && sheetness)
    {
        if (!write_nifti_float(out_sheetness, sheetness, DT_FLOAT32, slope, dims,
                               voxelsize, label_ptr))
            exit(EXIT_FAILURE);
    }

    free(t1);
    free(label);
    if (sheetness)
        free(sheetness);

    return (EXIT_SUCCESS);
}
