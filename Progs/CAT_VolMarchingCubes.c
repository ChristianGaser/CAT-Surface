/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 *
 */

#include <ParseArgv.h>
#include "CAT_MarchingCubes.h"
#include "CAT_Sheetness.h"
#include "CAT_SurfaceIO.h"

/* argument defaults */
char *label_filename = NULL;
double min_threshold = 0.5;
double pre_fwhm = 2.0;
double dist_morph = FLT_MAX;
double strength_gyri_mask = 0.1;
double strength_sulci = -1.0;
double sulci_cutoff = -1.0;
double sulci_sheet_strength = -1.0;
double sulci_thresh = -1.0;
double sulci_band = -1.0;
double sulci_normalize = -1.0;
int sulci_skeleton = -1;
double sulci_sigma_factor = -1.0;
double sulci_sigma_min = -1.0;
double sulci_sigma_max = -1.0;
int sulci_scales = -1;
double sulci_offset = -1.0;
double sulci_offset_gyri = -1.0;
int iter_laplacian = 50;
int n_median_filter = 2;
int verbose = 0;
int n_iter = 5;
int fast = 0;

/* the argument table */
static ArgvInfo argTable[] = {
  {"-label", ARGV_STRING, (char *) NULL, (char *) &label_filename, 
    "File containing segmentation labels for creating smooth mask of gyral\n\
    and sulcal areas. This prevents sulcal closure by using a higher isovalue\n\
    in sulci, and prevent cutting gyri by using a lower isovalue in gyri."},

  {"-strength-gyrimask", ARGV_FLOAT, (char *) TRUE, (char *) &strength_gyri_mask,
    "Define the strength of isovalue correction using the smooth.\n\
     gyri/sulci mask from the label map."},

  {"-strength-sulci", ARGV_FLOAT, (char *) TRUE, (char *) &strength_sulci,
    "Strength of the buried-sulcus correction on the PPM (default 0 = off).\n\
     Raising it enables the correction; the remaining -sulci-* options are the\n\
     values tuned on real data and take effect only once this is set.\n\
     A buried sulcus is a valley in the PPM whose floor never drops below the\n\
     isovalue, so the two banks fuse when the isosurface is extracted. No\n\
     intensity image is needed here: crossing a sulcus the PPM runs 1 -> 0.5 ->\n\
     ~0 -> 0.5 -> 1 and crossing a gyral blade it runs 0 -> 0.5 -> ~1 -> 0.5 ->\n\
     0, so a sulcus is a valley and a blade a ridge, and a Hessian sheetness\n\
     filter separates them by polarity alone. The detected valley field is used\n\
     three times: floors sitting just above the isovalue are pushed below it;\n\
     the gyral boost of -strength-gyrimask is damped there, because\n\
     strengthening a thin white-matter finger otherwise lifts the neighbouring\n\
     sulcal floor back over the isovalue; and the median filter is oriented\n\
     along the sheet so it cannot re-close what was opened. The correction only\n\
     ever lowers a value, acts only on values within 0.25 of the isovalue, and\n\
     requires a local minimum along the sheet normal -- which a gyral crown can\n\
     never satisfy, so it cannot cut a gyrus."},

  {"-sulci-sheet-strength", ARGV_FLOAT, (char *) TRUE, (char *) &sulci_sheet_strength,
    "Gain on the sheetness response used by -strength-sulci (library default 30).\n\
     The automatic noise scale of the sheetness filter is half the largest\n\
     Hessian norm in the volume, which on real data is set by the cortical\n\
     ribbon itself. A thin sulcal valley is far weaker than that, so the raw\n\
     response typically sits an order of magnitude below -sulci-thresh and\n\
     -strength-sulci does nothing at all at a gain of 1. This is the same\n\
     reason the sheetness gain exists at all, and the value\n\
     does NOT carry over: that one is measured on the intensity image, this\n\
     one on the PPM. Run with -verbose -- it reports the p99 and maximum of\n\
     the response next to the threshold, and warns when the gain is too low."},

  {"-sulci-sigma-factor", ARGV_FLOAT, (char *) TRUE, (char *) &sulci_sigma_factor,
    "Largest sheetness scale as a multiple of the median cortical thickness\n\
     (default 1.25; 0 or less keeps -sulci-sigma-max as given). The structure\n\
     the filter has to find is a valley whose width is set by how far apart the\n\
     two banks stand, so the scale belongs at a multiple of this brain's\n\
     thickness rather than at a fixed millimetre value. The thickness is read\n\
     out of the PPM itself -- the map climbs by 1 across the ribbon, so its\n\
     gradient is 1/thickness -- so no thickness map, label map or intensity\n\
     image is needed, which matters because none of them exist at this point.\n\
     Run with -verbose to see the thickness measured and the scale derived."},

  {"-sulci-sigma-min", ARGV_FLOAT, (char *) TRUE, (char *) &sulci_sigma_min,
    "Smallest sheetness scale in mm (library default 0.3)."},

  {"-sulci-sigma-max", ARGV_FLOAT, (char *) TRUE, (char *) &sulci_sigma_max,
    "Largest sheetness scale in mm. Overrides -sulci-sigma-factor when set.\n\
     lever on what gets found: a sulcus wider than the largest scale is not seen\n\
     as a sheet at all. Raise it if glued sulci persist, and use -sulci-skeleton\n\
     so the wider response still lands on the midline."},

  {"-sulci-scales", ARGV_INT, (char *) TRUE, (char *) &sulci_scales,
    "Number of log-spaced sheetness scales (library default 3)."},

  {"-sulci-thresh", ARGV_FLOAT, (char *) TRUE, (char *) &sulci_thresh,
    "Sheetness below this is ignored (default 0.1). Aim for a gain that puts\n\
     the p99 of the response near twice this value, matching the convention\n\
     the other sheetness-gated steps use."},

  {"-sulci-band", ARGV_FLOAT, (char *) TRUE, (char *) &sulci_band,
    "Only values in (isovalue, isovalue + band) are opened (default 0.25).\n\
     Raising it reaches more deeply buried sulci at the cost of touching\n\
     tissue further from the surface."},

  {"-sheet-offset", ARGV_FLOAT, (char *) TRUE, (char *) &sulci_offset,
    "Signed sheetness offset applied to the PPM before the surface is extracted,\n\
     in map units (library default 0.6; 0 disables it). A sulcus is a valley in\n\
     the PPM and a gyral\n\
     blade a ridge, so a signed sheetness is negative on one and positive on the\n\
     other; adding it lowers the map along sulci and raises it along blades in a\n\
     single pass. This is what lowering -thresh globally cannot do -- that moves\n\
     every voxel the same way, so it severs thin white-matter fingers at the same\n\
     rate as it opens glued sulci. The offset is scaled by the response, so it is\n\
     strongest where the structure is clearest and vanishes where nothing was\n\
     found. Requires -strength-sulci to be set (it shares the same sheetness\n\
     scales, gain, normalization and skeleton settings). Start around 0.05-0.1:\n\
     the isovalue is 0.5, so 0.1 is a fifth of the way to either extreme."},

  {"-sheet-offset-gyri", ARGV_FLOAT, (char *) TRUE, (char *) &sulci_offset_gyri,
    "Offset for the raising half of the signed map -- the gyral blades -- in map\n\
     units. Unset (default) means the same value as -sheet-offset, i.e. the signed\n\
     map is applied whole.\n\
     Use it when sulci stay glued at an offset that is already large. The two\n\
     halves do different jobs and only one of them can re-glue banks: lowering a\n\
     valley opens a sulcus, while raising a ridge protects a thin blade but also\n\
     lifts the sulcal floor beside it. On a real PPM the raising half is the\n\
     larger of the two -- at offset 0.6 on an OASIS subject it lifted 134214\n\
     voxels over the isovalue against 85450 pushed under, so the balanced offset\n\
     added tissue on net. Dropping this to 0 on that subject opened exactly as\n\
     many glued voxels (the raising half contributes none) and removed the\n\
     inflation. Lower it before raising -sheet-offset, which past saturation does\n\
     nothing. -verbose reports both crossing counts."},

  {"-no-sulci-skeleton", ARGV_CONSTANT, (char *) FALSE, (char *) &sulci_skeleton,
    "Do not thin the sulcal valley field to its medial sheet. Thinning is off by\n\
     default because the offset and the barrier should act at the structure\n\
     rather than across a band as wide as the Gaussian that found it; this\n\
     switch is here so the two can be compared. Note the thinned field is a\n\
     weaker correction at the same -sheet-offset -- roughly a third of the\n\
     isovalue crossings -- so raise the offset by about 3x when comparing."},

  {"-sulci-skeleton", ARGV_CONSTANT, (char *) TRUE, (char *) &sulci_skeleton,
    "Thin the sulcal valley field to its medial sheet before any threshold is\n\
     applied. Off by default.\n\
     applied. The plate response is as wide as the Gaussian that produced it, so\n\
     a large -sulci-sigma-max locates the valleys well but answers several voxels\n\
     into the banks on either side; the per-voxel gates then reach into that\n\
     tissue too. Suppressing everything that is not a maximum along its own\n\
     normal collapses the band onto its ridge line -- one voxel at any scale --\n\
     and leaves the value on the ridge unchanged."},

  {"-sulci-normalize", ARGV_FLOAT, (char *) TRUE, (char *) &sulci_normalize,
    "Value the p99.9 of the valley response is scaled to (default 1.0); pass 0\n\
     to keep the raw response. This is what makes -sulci-thresh mean the same\n\
     thing on every dataset."},

  {"-sulci-cutoff", ARGV_FLOAT, (char *) TRUE, (char *) &sulci_cutoff,
    "Admission cutoff of the oriented median used with -strength-sulci\n\
     (library default 0.4; 0 selects CAT_ORIENTED_MEDIAN_CUTOFF). It is half the sheetness at\n\
     which a one-voxel-thick sheet starts being preserved, so it has to match\n\
     the response the PPM actually produces -- measure it with CAT_VolSheetness\n\
     -polarity -1 on the PPM and halve what the sulci reach."},
  
  {"-thresh", ARGV_FLOAT, (char *) TRUE, (char *) &min_threshold,
    "Define the volume threshold, also known as the isovalue.\n\
     This value is crucial for initial image thresholding."},
  
  {"-pre-fwhm", ARGV_FLOAT, (char *) TRUE, (char *) &pre_fwhm,
    "Specify the Full Width Half Maximum (FWHM) for the preprocessing\n\
     smoothing filter. This helps in preserving gyri and sulci by\n\
     creating a weighted average between original and smoothed\n\
     images based on the gradient of the input image. Areas with \n\
     topology artefacts are often characterized by large gradients,\n\
     thus smoothing in these areas tries to prevent these artefacts."},

  {"-iter-laplacian", ARGV_INT, (char *) TRUE, (char *) &iter_laplacian,
    "Set number of iterations for Laplacian surface smoothing. This aids\n\
     in removing noise from the mesh"},
  
  {"-dist-morph", ARGV_FLOAT, (char *) TRUE, (char *) &dist_morph,
    "Apply initial morphological opening or closing step. Closing is used\n\
     by a value around 1.0 and opening by negative values around -1.0.\n\
     The default automatically estimates the optimal value"},
  
  {"-median-filter", ARGV_INT, (char *) TRUE, (char *) &n_median_filter,
    "Specify the number of iterations to apply a median filter to areas\n\
     where the gradient of the thresholded image indicates larger clusters.\n\
     These clusters may point to potential topology artifacts and regions\n\
     with high local variations. This process helps to smooth these areas, \n\
     improving the quality of the surface reconstruction in subsequent steps."},
  
  {"-iter", ARGV_INT, (char *) TRUE, (char *) &n_iter,
    "Maximum number of iterations."},
  
  {"-fast", ARGV_CONSTANT, (char *) TRUE, (char *) &fast,
    "Enable fast processing without any preprocessing, smoothing, or topology \n\
    correction. Only the final Laplacian surface smoothing is applied if defined."},

  {"-verbose", ARGV_CONSTANT, (char *) TRUE, (char *) &verbose,
    "Enable verbose mode for detailed output during processing."},

    {NULL, ARGV_END, NULL, NULL, NULL}
};


private void
usage(
    char *executable)
{
    char *usage_str = "\n\
Usage: %s [options] <input.nii> <output_surface_file> [<change_map.nii>]\n\
\n\
    This method generates a mesh with an Euler number of 2 (genus 0) from the\n\
    thresholded volume. The process involves:\n\
    \n\
    1. **Preprocessing with Smoothing Filter:**\n\
       - Apply a smoothing filter to the input image to remove outliers.\n\
       - Use a weighted average of the original and smoothed images to\n\
         preserve gyri and sulci.\n\
       - Weighting is based on the gradient of the input image.\n\
       - Weights range from 0 (areas with low gradient) to 1 (areas\n\
         with large gradient), with intermediate values scaled linearly.\n\
       - Weighting effect is enhanced by squaring the value.\n\
    \n\
    2. **Preprocessing with Median Filter:**\n\
       - Apply an iterative median filter to remove noise.\n\
    \n\
    3. **Optionally Use of smooth mask of gyri and sulci:**\n\
        - This prevents sulcal closure by using a higher isovalue in sulci,\n\
          and prevent cutting gyri by using a lower isovalue in gyri\n\
          An optional label maps must be defined to use this approach.\n\
    \n\
    4. **Morphological Opening:**\n\
       - Apply additional morphological opening or closing, defined by\n\
         `-dist-morph`, to minimize changes due to topology correction.\n\
       - Closing is used for positive values (e.g. 1.0) and opening\n\
         for negative values. The default is to automatically estimate\n\
         the optimal value to minimize issues due to topology correction.\n\
    \n\
    5. **Extraction of the Largest Component:**\n\
       - Extract the largest component for further processing.\n\
    \n\
    6. **Mesh Smoothing:**\n\
       - Smooth the extracted mesh with a Laplacian filter.\n\n";

    fprintf(stderr, usage_str, executable);
}

int main(int argc, char *argv[]) {
    float *input_float, *label;
    char out_diff[1024];
    char *input_filename, *output_filename;
    object_struct *object;

    /* get the arguments from the command line */
    if (ParseArgv(&argc, argv, argTable, 0)) {
        usage(argv[0]);
        fprintf(stderr, "     %s -help\n\n", argv[0]);
        return EXIT_FAILURE;
    }

    initialize_argument_processing(argc, argv);

    if (!get_string_argument(NULL, &input_filename) || !get_string_argument(NULL, &output_filename)) {
        usage(argv[0]);
        fprintf(stderr, "Usage: CAT_VolMarchingCubes input.nii output_surface_file [options] [change_map.nii]\n");
        return EXIT_FAILURE;
    }

    /* Read input volume */
    nifti_image *nii_ptr = read_nifti_float(input_filename, &input_float, 0);
    if (!nii_ptr) {
        fprintf(stderr, "Error reading %s.\n", input_filename);
        return EXIT_FAILURE;
    }

    if (label_filename) {
        nifti_image *nii_ptr2 = read_nifti_float(label_filename, &label, 0);
        if (!nii_ptr2) {
            fprintf(stderr, "Error reading %s.\n", label_filename);
            return EXIT_FAILURE;
        }
        if ((nii_ptr->nx != nii_ptr2->nx) ||
            (nii_ptr->ny != nii_ptr2->ny) ||
            (nii_ptr->nz != nii_ptr2->nz)) {
            fprintf(stderr, "Error: Label image must have the same dimensions as the input image.\n");
            return EXIT_FAILURE;
        }
    } else label = NULL;

    if (fast) {
        object = apply_marching_cubes_fast(input_float, nii_ptr, 
                    min_threshold, iter_laplacian, verbose);
    } else {
        CAT_PpmSulciOpts sulci_opts;
        CAT_PpmSulciOpts *sulci_ptr = NULL;
        if (strength_sulci != 0.0) {
            /* CAT_PpmSulciOptionsInit() is the single source of truth for the
               defaults; a negative value here means the user did not ask for
               anything and the library default stands.  Duplicating the numbers
               in the front-ends is how the CLI and the Python binding drifted
               apart in the first place. */
            CAT_PpmSulciOptionsInit(&sulci_opts);
            if (sulci_sheet_strength >= 0.0) sulci_opts.sheet_strength = sulci_sheet_strength;
            if (sulci_normalize      >= 0.0) sulci_opts.sheet_normalize = sulci_normalize;
            if (sulci_skeleton       >= 0)   sulci_opts.sheet_skeleton = sulci_skeleton;
            if (sulci_thresh         >= 0.0) sulci_opts.thresh = sulci_thresh;
            if (sulci_band           >= 0.0) sulci_opts.band = sulci_band;
            if (sulci_sigma_factor   >= 0.0) sulci_opts.sigma_factor = sulci_sigma_factor;
            if (sulci_sigma_min      >= 0.0) sulci_opts.sigma_min = sulci_sigma_min;
            if (sulci_sigma_max      >= 0.0)
            {
                sulci_opts.sigma_max = sulci_sigma_max;
                /* An explicit sigma_max means an explicit sigma_max.  The
                   derivation from the PPM's own thickness runs whenever
                   sigma_factor is positive and overwrites it, so asking for a
                   scale on the command line has to switch the derivation off --
                   otherwise the option is silently ignored.  Passing both keeps
                   the factor, since that is the only way to ask for the derived
                   value with a floor. */
                if (sulci_sigma_factor < 0.0)
                    sulci_opts.sigma_factor = 0.0;
            }
            if (sulci_scales         >= 1)   sulci_opts.n_scales = sulci_scales;
            if (sulci_offset  >= 0.0) sulci_opts.offset = sulci_offset;
            if (sulci_offset_gyri >= 0.0)
                sulci_opts.offset_gyri = sulci_offset_gyri;
            if (strength_sulci >= 0.0) sulci_opts.strength = strength_sulci;
            if (sulci_cutoff   >= 0.0) sulci_opts.cutoff = sulci_cutoff;
            sulci_opts.verbose = verbose;
            sulci_ptr = &sulci_opts;
        }

        object = apply_marching_cubes(input_float, nii_ptr, label, 
                    min_threshold, pre_fwhm, iter_laplacian, dist_morph, n_median_filter, 
                    n_iter, strength_gyri_mask, sulci_ptr, verbose);
    }
    if (object) {
        output_graphics_any_format(output_filename, ASCII_FORMAT, 1, &object, NULL);
    } else {
        fprintf(stderr, "Error generating surface.\n");
    }

    if(argc > 3) {
        double voxelsize[N_DIMENSIONS];
        int dims[MAX_DIMENSIONS];
        dims[0] = nii_ptr->nx;
        dims[1] = nii_ptr->ny;
        dims[2] = nii_ptr->nz;
        voxelsize[0] = nii_ptr->dx;
        voxelsize[1] = nii_ptr->dy;
        voxelsize[2] = nii_ptr->dz;
        (void) sprintf(out_diff, "%s", argv[3]); 
        if (!write_nifti_float(out_diff, input_float, DT_FLOAT32, 1.0, dims, voxelsize, nii_ptr)) 
            exit(EXIT_FAILURE);
    }

    free(input_float);
    delete_marching_cubes_table();
    return 0;
}
