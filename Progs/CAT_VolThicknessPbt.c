/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 *
 */

#include <float.h>
#include <math.h>
#include <stdlib.h>

#if !defined(_WIN32) && !defined(_WIN64)
#include <libgen.h>
#endif

#include <bicpl.h>
#include <ParseArgv.h>
#include "CAT_NiftiLib.h"
#include "CAT_Vol.h"
#include "CAT_VolPbt.h"
#include "CAT_Math.h"

int fast = 0;
int verbose = 0;
int n_avgs = -1;
int n_median_filter = -1;
int median_subsample = -1;
int blood_vessel_correction = 1;
double range = -1.0;
double downsample = 0.0;
double fill_thresh = -1.0;
double correct_thickness = NAN;
double sulcal_width = -1.0;
int pve_distance = 0;
int oriented_filter = 0;
int sulcal_barrier = 0;
double barrier_q = -1.0;
double barrier_dmin = -1.0;
double barrier_gmtmax = -1.0;
double barrier_gmtfactor = -1.0;
double barrier_gmtpct = -1.0;
double barrier_ramp = -1.0;
double barrier_local = -1.0;
double barrier_tmin = -1.0;
double barrier_halfwidth = -1.0;
double oriented_strength = -1.0;
double oriented_cutoff = -1.0;
char *thickness_offset_file = NULL;

static ArgvInfo argTable[] = {
    {"-verbose", ARGV_CONSTANT, (char *)1, (char *)&verbose,
     "Enable verbose mode. Provides detailed output during processing for debugging\n\
    and monitoring."},

    {"-fast", ARGV_CONSTANT, (char *)1, (char *)&fast,
     "Enable fast mode in order to get a very quick and rougher estimate of thickness only."},

    {"-no-blood-vessel-correction", ARGV_CONSTANT, (char *)0, (char *)&blood_vessel_correction,
     "Enable blood-vessel correction before thickness estimation (0 disables, >0 enables)."},

    {"-no-bvc", ARGV_CONSTANT, (char *)0, (char *)&blood_vessel_correction,
     "Enable blood-vessel correction before thickness estimation (0 disables, >0 enables)."},

    {"-n-avgs", ARGV_INT, (char *)1, (char *)&n_avgs,
     "Number of averages for distance estimation (library default 5). Used for averaging\n\
    the distances in White Matter (WM) and Cerebrospinal Fluid (CSF) to obtain a\n\
    less noisy measure. A higher number results in smoother but potentially less\n\
    accurate measures."},

    {"-fill-holes", ARGV_FLOAT, (char *)1, (char *)&fill_thresh,
     "Fill remaining holes in the PPM image using the defined threshold.\n\
    To maximize the filling-effect, this threshold should be the same as used for the\n\
    subsequent Marching Cubes approach (e.g. 0.5). Set to '0' to disable filling."},

    {"-range", ARGV_FLOAT, (char *)1, (char *)&range,
     "Extend range for masking of euclidean distance estimation. A slight increase\n\
    of range (i.e 0.3) helps in obtaining a more stable distance estimation."},

    {"-downsample", ARGV_FLOAT, (char *)1, (char *)&downsample,
     "Downsample PPM and GMT image to defined resolution since we do not need that 0.5mm\n\
    spatial resolution for the subsequent steps. Set to '0' to disable downsampling."},

    {"-median-filter", ARGV_INT, (char *)TRUE, (char *)&n_median_filter,
     "Specify the number of iterations for weighted local median filtering of the\n\
        final PPM image. The filter is not applied uniformly: CAT first estimates a\n\
        topology-artifact weight map from the positive residual PPM - smooth(PPM),\n\
        keeps only high-residual voxels in sufficiently thick cortex, regularizes this\n\
        mask by close/open/dilate and smoothing, and then blends original PPM with the\n\
        locally median-filtered PPM. Higher weights mean stronger median-filter\n\
        influence. Set to 0 to disable this cleanup."},

    {"-median-subsample", ARGV_INT, (char *)TRUE, (char *)&median_subsample,
     "Specify the size of subsampling for the median filter to smooth local\n\
     thickness values"},

    {"-sulcal-barrier", ARGV_CONSTANT, (char *)1, (char *)&sulcal_barrier,
     "Stop the CSF distance at the sulcal medial surface.\n\
     Where the classifier lost the CSF in a sulcus there is no boundary for the\n\
     CSF distance to stop at, so the front from one bank runs through the fused\n\
     grey matter into the other; the thickness follows it, the PPM never drops\n\
     below the isovalue, and the buried sulcus is created there -- in the\n\
     distance map, long before marching cubes is asked to draw it.\n\
     The midline the front should have stopped at is recovered from geometry\n\
     rather than intensity: it is where the fronts from the two banks collide.\n\
     No intensity image and no per-subject threshold are involved.\n\
     The bound is applied as a minimum, so where CSF was segmented properly it\n\
     is nearer than the midline and the result is identical to leaving this off;\n\
     the distance can only shrink, so an overestimated thickness can be fixed\n\
     but a correct one can never be inflated."},

    {"-barrier-local", ARGV_FLOAT, (char *)1, (char *)&barrier_local,
     "FWHM in mm over which the gate follows regional thickness (library default\n\
     0, i.e. one global gate; the regional form measured neutral on three subjects). Cortex is not one thickness -- occipital runs\n\
     near 2 mm while insular and temporal reach 3.5 mm -- so a single gate is at\n\
     once too high for the thin regions, where two glued 2 mm banks imply only\n\
     4 mm and never reach it, and too low for the thick ones, where it catches\n\
     ordinary tissue and leaves a visible seam. Gating against the local value\n\
     fixes both, and a smoothly varying gate leaves no boundary to see."},

    {"-barrier-ramp", ARGV_FLOAT, (char *)1, (char *)&barrier_ramp,
     "Width over which the gate fades in, as a fraction of the gate itself\n\
     (library default 0.5; 0 restores a hard threshold). A hard threshold\n\
     switches the correction on between one voxel and its neighbour, so on a\n\
     thick cortex -- where much of the band sits near the gate -- capped and\n\
     uncapped tissue end up side by side and the seam shows in the thickness\n\
     map as a step. Fading the cap in over a range removes the seam without\n\
     changing what is corrected well above the gate or left alone well below."},

    {"-barrier-gmtpct", ARGV_FLOAT, (char *)1, (char *)&barrier_gmtpct,
     "Percentile below which the thickness proxy is averaged (library default\n\
     90; 100 gives a plain mean). The proxy dist_WM + dist_CSF runs high\n\
     because the glued sulci the gate exists to find sit in its upper tail, and\n\
     a median only limits their influence -- it still sits inside a distribution\n\
     they have skewed. Cutting the tail off and averaging what is left tracks\n\
     the cortex more closely: measured against the reported GMT on four\n\
     hemispheres from two datasets, the ratio spans 0.087 for this estimator\n\
     against 0.102 for a median."},

    {"-barrier-gmtfactor", ARGV_FLOAT, (char *)1, (char *)&barrier_gmtfactor,
     "Multiple of the median thickness at which the gate sits (library default 1.75; 0 or\n\
     less disables the gate). This is the criterion in its natural form: a glued\n\
     sulcus is two cortices back to back, so the threshold belongs at twice the\n\
     typical thickness of the brain being processed, not at a fixed millimetre\n\
     value that is only right for the cortex it was tuned on. The median is\n\
     taken over dist_WM + dist_CSF inside the GM band -- for a band of locally\n\
     constant thickness those are complementary and sum to it exactly -- and a\n\
     median is unmoved by the glued minority the gate exists to catch. Use\n\
     -barrier-gmtmax to override it with an absolute value."},

    {"-barrier-gmtmax", ARGV_FLOAT, (char *)1, (char *)&barrier_gmtmax,
     "Absolute override for the gate, in mm (default 0 = derive it from\n\
     -barrier-gmtfactor and the data). This is the gate that matters. A glued sulcus\n\
     does not merely look thick, it looks like two cortices back to back -- 5-6\n\
     mm where 2-3 mm is normal -- so the implied thickness at the voxel,\n\
     dist_WM + dist_CSF, separates the two populations cleanly. Gating on the\n\
     CSF distance alone does not: plenty of ordinary voxels sit more than 2 mm\n\
     from CSF simply by being near the white matter."},

    {"-barrier-dmin", ARGV_FLOAT, (char *)1, (char *)&barrier_dmin,
     "Only bound a CSF distance that is already implausible for real cortex, in\n\
     mm (default 2.0; 0 disables the gate). Bounding a distance always shrinks\n\
     the thickness, so without this the barrier lowers GMT wherever it happens\n\
     to fire and the mean thickness of the whole brain becomes a readout of how\n\
     often that was -- which makes the mean move with -barrier-q, exactly the\n\
     parameter dependence the barrier was meant to remove. A voxel in the middle\n\
     of a 2.5-4 mm band sits 1.2-2 mm from CSF; only a front that ran across a\n\
     glued sulcus comes back with more, so this separates the two and keeps\n\
     normally-thick cortex untouched."},

    {"-barrier-q", ARGV_FLOAT, (char *)1, (char *)&barrier_q,
     "Shock threshold of the medial set (library default 0.7).\n\
     Lower is stricter and gives a thinner midline."},

    {"-barrier-tmin", ARGV_FLOAT, (char *)1, (char *)&barrier_tmin,
     "Ignore collisions closer than this to the WM boundary, in mm (default\n\
     0.5). This is what stops the barrier carving into thin cortex from the\n\
     inside."},

    {"-barrier-halfwidth", ARGV_FLOAT, (char *)1, (char *)&barrier_halfwidth,
     "Half the width of the CSF sheet the barrier stands for, in mm (default 0,\n\
     which selects half a voxel). The distance transform measures to the medial\n\
     voxel centre while the grey matter ends at that sheet's surface, so the raw\n\
     distance is short by half its width."},

    {"-oriented-filter", ARGV_CONSTANT, (char *)1, (char *)&oriented_filter,
     "Replace the isotropic 3x3x3 median filters by sheetness-oriented ones.\n\
     An isotropic median penalizes boundary area, so it removes thin structures\n\
     regardless of which side of the label boundary they lie on: the same filter\n\
     that opens a glued sulcus closes a cerebellar fissure, and tuning trades one\n\
     against the other. The oriented variant estimates a Hessian sheetness field\n\
     once and then admits only those neighbours that lie in the plane of the local\n\
     sheet, so it averages along a thin structure and never across it. Where no\n\
     sheet is detected every neighbour is admitted and the filter is identical to\n\
     the isotropic one, so nothing changes away from thin structures."},

    {"-oriented-strength", ARGV_FLOAT, (char *)1, (char *)&oriented_strength,
     "Overall gain on the sheetness before the oriented filters use it (default\n\
     1.0). 0 reproduces the isotropic filters exactly. Values above 1 amplify a\n\
     response too weak to matter: the oriented median admits every neighbour\n\
     unless the sheetness exceeds 0.5, so a map peaking at 0.5 leaves it\n\
     bit-identical to the isotropic median. Inspect the map with CAT_VolSheetness\n\
     before raising this, since the noise floor is amplified along with it."},

    {"-oriented-cutoff", ARGV_FLOAT, (char *)1, (char *)&oriented_cutoff,
     "Admission cutoff of the oriented medians (default 0.10). A neighbour at\n\
     offset d is admitted when sheetness*(dhat.n)^2 < cutoff, so the 6 face\n\
     neighbours drop out at sheetness = cutoff, the 12 edge neighbours at\n\
     2*cutoff and the 8 corners at 3*cutoff, while the 9 offsets in the sheet\n\
     plane are always admitted. A one-voxel-thick sheet is preserved from\n\
     2*cutoff upwards. Set it to about half the sheetness your data actually\n\
     reaches; inspect that with CAT_VolSheetness."},

    {"-correct-thickness", ARGV_FLOAT, (char *)1, (char *)&correct_thickness,
     "Additive thickness correction in mm, compensating the systematic shift of the\n\
     GM/WM and GM/CSF borders produced by the segmentation."},

    {"-pve-distance", ARGV_CONSTANT, (char *)1, (char *)&pve_distance,
     "EXPERIMENTAL. Correct the WM/CSF distance maps for partial volume. The distance\n\
     transform measures to the nearest source voxel centre, but the tissue boundary\n\
     lies a sub-voxel distance beyond that centre; how far is given by the partial\n\
     volume of the source voxel and the width of the PVE ramp, which is estimated\n\
     from the image. Mirrors the correction in CAT12's cat_vol_pbtsimpleCS4.\n\
     This shifts the thickness calibration, so '-correct-thickness' has to be\n\
     re-derived when this is enabled."},

    {"-sulcal-width", ARGV_FLOAT, (char *)1, (char *)&sulcal_width,
     "Maximum distance from CSF boundary (mm) for sulcal PPM correction. Voxels in\n\
     the GM/CSF partial-volume zone within this distance are attenuated to prevent\n\
     sulcal gluing in the central surface. Set to 0 to disable."},

    {"-thickness-offset", ARGV_STRING, (char *)1, (char *)&thickness_offset_file,
     "Per-voxel additive thickness correction in mm, as a volume on the input\n\
     grid. A spatially varying counterpart to -correct-thickness, for the case\n\
     where the border shift is not one number for the brain -- notably\n\
     myelinated cortex, where the classifier puts the GM/WM boundary too far\n\
     out and the ribbon comes back too thin. Produce one with\n\
     CAT_VolBoundaryOffset. Applied after the PPM, so it moves the thickness\n\
     and leaves the surfaces alone."},
    {NULL, ARGV_END, NULL, NULL, NULL}};

private void usage(char *executable)
{
    char *usage_str = "\n\
Usage: %s [options] <input.nii> <output_GMT.nii> <output_PPM.nii>\n\
                    [<output_WMD.nii> <output_CSD.nii>]\n\
\n\
    Projection-based cortical thickness and percentage position mapping (PPM)\n\
    from a PVE label image (CSF=1, GM=2, WM=3), after:\n\
\n\
        Dahnke R, Yotter RA, Gaser C.\n\
        Cortical thickness and central surface estimation.\n\
        Neuroimage. 2013 Jan 15;65:336-48.\n\
\n\
    The input is a PVE label image; GMT (thickness, in mm) and PPM (relative\n\
    position in the cortical ribbon, 0 at the pial boundary and 1 in the white\n\
    matter) are always written. The distance maps are written as well when the\n\
    two optional output names are given.\n\
\n\
    The steps are:\n\
\n\
    1. Distance estimation. Distances to the WM and CSF boundaries are measured\n\
       by shifting the GM/WM and GM/CSF thresholds and averaging over -n-avgs\n\
       levels, which makes the measure less noisy.\n\
\n\
    2. Thickness estimation. A sulcal and a gyral estimate are projected through\n\
       the cortical band and the smaller of the two is kept.\n\
\n\
    3. PPM construction, followed by hole filling and an optional weighted local\n\
       median cleanup where a topology-artifact map is high.\n\
\n\
    Two optional corrections address the case where the classifier lost the CSF\n\
    inside a sulcus, which otherwise inflates the thickness and buries the\n\
    sulcus in the surface:\n\
\n\
    -sulcal-barrier bounds the CSF distance by the sulcal midline, located from\n\
    geometry alone by front collision, so no intensity image and no per-subject\n\
    threshold are involved. It is applied as a minimum, so where CSF was found\n\
    the result is unchanged, and it can only ever reduce an overestimated\n\
    thickness. -barrier-gmtfactor keeps it to voxels whose implied thickness is\n\
    impossible for cortex; run with -verbose to see the fraction of the GM band\n\
    it touched, which should be a few percent.\n\
\n\
    -oriented-filter replaces the isotropic median filters with sheetness-\n\
    oriented ones, which cannot close a thin structure. Where no sheet is\n\
    detected it is identical to the isotropic filter.\n\
\n\
    Every option is listed under 'Command-specific options' above. The values\n\
    shown there are the library defaults from CAT_PbtOptionsInit(), which is\n\
    the single source of truth -- an option left unset keeps whatever the\n\
    library specifies rather than a number duplicated in this tool.\n\
\n\
Examples:\n\
    %s input.nii gmt.nii ppm.nii\n\
    %s -sulcal-barrier -verbose input.nii gmt.nii ppm.nii\n\
    %s -sulcal-barrier -oriented-filter input.nii gmt.nii ppm.nii wmd.nii csd.nii\n\n";

    fprintf(stderr, usage_str, executable, executable, executable, executable);
}

int main(int argc, char *argv[])
{
    char out_GMT[1024], out_PPM[1024], out_CSD[1024], out_WMD[1024];
    int i, j, dims[3], dims_reduced[3];
    float *src;
    float *offset_data = NULL;
    double voxelsize[3], voxelsize_reduced[3], samp[3], s[3], slope;

    initialize_argument_processing(argc, argv);

    if (ParseArgv(&argc, argv, argTable, 0) || (argc < 2))
    {
        usage(argv[0]);
        fprintf(stderr, "     %s -help\n\n", argv[0]);
        exit(EXIT_FAILURE);
    }

    char *infile = argv[1];

    /* Determine output filenames based on input filename or command-line arguments */
    if (argc >= 6)
    {
        (void)sprintf(out_WMD, "%s", argv[4]);
        (void)sprintf(out_CSD, "%s", argv[5]);
    }
    if (argc >= 4)
    {
        (void)sprintf(out_GMT, "%s", argv[2]);
        (void)sprintf(out_PPM, "%s", argv[3]);
    }
    else
    {
#if !defined(_WIN32) && !defined(_WIN64)
        (void)sprintf(out_GMT, "%s/gmt_%s", dirname(infile), basename(infile));
        (void)sprintf(out_PPM, "%s/ppm_%s", dirname(infile), basename(infile));
#else
        fprintf(stderr, "\nUsage: %s input.nii GMT.nii PPM.nii\n\n", argv[0]);
        return (1);
#endif
    }

    /* read source image */
    nifti_image *src_ptr = read_nifti_float(infile, &src, 0);
    if (!src_ptr)
    {
        fprintf(stderr, "Error reading %s.\n", infile);
        return (EXIT_FAILURE);
    }

    /* Optional per-voxel thickness correction, on the input grid */
    if (thickness_offset_file)
    {
        nifti_image *off_ptr = read_nifti_float(thickness_offset_file,
                                                &offset_data, 0);
        if (!off_ptr)
        {
            fprintf(stderr, "Error reading %s.\n", thickness_offset_file);
            return (EXIT_FAILURE);
        }
        if (off_ptr->nx != src_ptr->nx || off_ptr->ny != src_ptr->ny ||
            off_ptr->nz != src_ptr->nz)
        {
            fprintf(stderr, "Error: %s does not match the input dimensions.\n",
                    thickness_offset_file);
            return (EXIT_FAILURE);
        }
    }

    /* Number of voxels */
    int nvox = src_ptr->nvox;

    /* Prepare output NIfTI images */
    nifti_image *out_ptr = nifti_copy_nim_info(src_ptr);

    /* Retrieve dimensions and voxel size from source image */
    voxelsize[0] = src_ptr->dx;
    voxelsize[1] = src_ptr->dy;
    voxelsize[2] = src_ptr->dz;
    dims[0] = src_ptr->nx;
    dims[1] = src_ptr->ny;
    dims[2] = src_ptr->nz;
    float *dist_CSF = (float *)malloc(sizeof(float) * nvox);
    float *dist_WM = (float *)malloc(sizeof(float) * nvox);
    float *GMT = (float *)malloc(sizeof(float) * nvox);
    float *PPM = (float *)malloc(sizeof(float) * nvox);

    /* check for memory faults */
    if (!dist_CSF || !dist_WM || !GMT || !PPM)
    {
        fprintf(stderr, "Memory allocation error\n");
        exit(EXIT_FAILURE);
    }

    /* Change some defaults for fast option */
    if (fast)
    {
        n_avgs /= 2;
        n_median_filter = 0;
        fill_thresh = 0.0;
        downsample = 0.0;
    }

    /* Ensure that n_avgs is at least 1 */
    n_avgs = (n_avgs < 1) ? 1 : n_avgs;

    /* Optional blood-vessel correction before any other operation */
    if (blood_vessel_correction > 0.0)
    {
        if (verbose)
            fprintf(stderr, "Apply blood-vessel correction on input PVE map.\n");
        blood_vessel_correction_pve_float(src, dims, voxelsize);
    }

    CAT_PbtOptions opts;
    CAT_PbtOptionsInit(&opts);
    /* CAT_PbtOptionsInit() is the single source of truth for the defaults; a
       sentinel here means the user did not ask for anything and the library
       default stands.  Duplicating the numbers in the front-ends is how the CLI
       and the Python binding drifted apart before (median_subsample was 2 in one
       and 4 in the other).  correct_thickness uses NaN rather than a negative
       sentinel because a negative correction is a legitimate value. */
    if (n_avgs            >= 0)   opts.n_avgs = n_avgs;
    if (n_median_filter   >= 0)   opts.n_median_filter = n_median_filter;
    if (median_subsample  >= 0)   opts.median_subsample = median_subsample;
    if (range             >= 0.0) opts.range = range;
    if (fill_thresh       >= 0.0) opts.fill_thresh = fill_thresh;
    if (!isnan(correct_thickness)) opts.correct_thickness = correct_thickness;
    if (sulcal_width      >= 0.0) opts.sulcal_width = sulcal_width;
    opts.pve_distance = pve_distance;
    opts.fast = fast;
    opts.verbose = verbose;
    opts.oriented_filter = oriented_filter;
    opts.sulcal_barrier = sulcal_barrier;
    if (barrier_q         >= 0.0) opts.barrier_q = barrier_q;
    if (barrier_dmin      >= 0.0) opts.barrier_dmin = barrier_dmin;
    if (barrier_gmtmax    >= 0.0) opts.barrier_gmtmax = barrier_gmtmax;
    if (barrier_gmtfactor >= 0.0) opts.barrier_gmtfactor = barrier_gmtfactor;
    if (barrier_gmtpct    >= 0.0) opts.barrier_gmtpct = barrier_gmtpct;
    if (barrier_ramp      >= 0.0) opts.barrier_ramp = barrier_ramp;
    if (barrier_local     >= 0.0) opts.barrier_local = barrier_local;
    if (barrier_tmin      >= 0.0) opts.barrier_tmin = barrier_tmin;
    if (barrier_halfwidth >= 0.0) opts.barrier_halfwidth = barrier_halfwidth;
    if (oriented_strength >= 0.0) opts.oriented_strength = oriented_strength;
    if (oriented_cutoff   >= 0.0) opts.oriented_cutoff = oriented_cutoff;

    if (CAT_VolComputePbt(src, GMT, PPM, dist_CSF, dist_WM, offset_data,
                          dims, voxelsize, &opts) != 0)
    {
        fprintf(stderr, "Error computing projection-based thickness.\n");
        exit(EXIT_FAILURE);
    }

    /* Downsample images */
    slope = 1.0;
    if (downsample > 0.0)
    {

        nifti_image *out_ptr_reduced = nifti_copy_nim_info(src_ptr);

        for (i = 0; i < 3; i++)
        {
            s[i] = 1.2;
            voxelsize_reduced[i] = downsample;
            samp[i] = voxelsize_reduced[i] / voxelsize[i];
        }

        /*  Define grid dimensions */
        for (i = 0; i < 3; i++)
            dims_reduced[i] = (int)ceil((dims[i] - 1) / ((double)samp[i])) + 1;

        out_ptr_reduced->nvox = dims_reduced[0] * dims_reduced[1] * dims_reduced[2];

        /* Correct affine matrix */
        for (i = 0; i < 3; i++)
        {
            for (j = 0; j < 3; j++)
            {
                out_ptr_reduced->sto_xyz.m[i][j] = out_ptr->sto_xyz.m[i][j] * samp[i];
            }
        }
        out_ptr_reduced->sto_ijk = nifti_mat44_inverse(out_ptr_reduced->sto_xyz);

        smooth3(GMT, dims, voxelsize, s, 0, DT_FLOAT32);
        smooth3(PPM, dims, voxelsize, s, 0, DT_FLOAT32);

        /* Save GMT and PPM image */
        float *vol_reduced = (float *)malloc(sizeof(float) * out_ptr_reduced->nvox);

        subsample3(GMT, vol_reduced, dims, dims_reduced, DT_FLOAT32);
        if (!write_nifti_float(out_GMT, vol_reduced, DT_FLOAT32, slope, dims_reduced, voxelsize_reduced, out_ptr_reduced))
            exit(EXIT_FAILURE);

        subsample3(PPM, vol_reduced, dims, dims_reduced, DT_FLOAT32);
        if (!write_nifti_float(out_PPM, vol_reduced, DT_FLOAT32, slope, dims_reduced, voxelsize_reduced, out_ptr_reduced))
            exit(EXIT_FAILURE);

        free(vol_reduced);
    }
    else
    {
        if (!write_nifti_float(out_GMT, GMT, DT_FLOAT32, slope, dims, voxelsize, out_ptr))
            exit(EXIT_FAILURE);

        if (!write_nifti_float(out_PPM, PPM, DT_FLOAT32, slope, dims, voxelsize, out_ptr))
            exit(EXIT_FAILURE);
    }

    if (argc >= 6)
    {
        if (!write_nifti_float(out_CSD, dist_CSF, DT_FLOAT32, slope, dims, voxelsize, out_ptr))
            exit(EXIT_FAILURE);
        if (!write_nifti_float(out_WMD, dist_WM, DT_FLOAT32, slope, dims, voxelsize, out_ptr))
            exit(EXIT_FAILURE);
    }

    free(dist_CSF);
    free(dist_WM);
    free(GMT);
    free(PPM);
    free(src);

    return (EXIT_SUCCESS);
}
