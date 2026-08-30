/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 */

/**
 * \file CAT_VolBoundaryOffset.c
 * \brief CLI tool writing the myelination-induced GM/WM boundary displacement.
 *
 * Usage:
 *   CAT_VolBoundaryOffset [options] label.nii t1w.nii offset.nii
 *
 * This is a diagnostic before it is a correction: run it on a subject and
 * look at whether the map picks out the central sulcus and V1.  Only a map
 * that does is worth feeding to PBT.
 */

#include <stdlib.h>
#include <stdio.h>

#include <bicpl.h>
#include <ParseArgv.h>

#include "CAT_NiftiLib.h"
#include "CAT_VolBoundaryOffset.h"

double ref_fwhm = 10.0;
double erosion_mm = 3.0;
double t_lo = 0.25;
double t_hi = 0.75;
double search_mm = 4.0;
double step_mm = 0.25;
double width_pct = 25.0;
double gain = 0.5;
double max_offset_mm = 1.5;
double smooth_fwhm = 8.0;
char *width_file = NULL;
int verbose = 0;

static ArgvInfo argTable[] = {
    {"-ref-fwhm", ARGV_FLOAT, (char *)TRUE, (char *)&ref_fwhm,
     "FWHM (mm) of the local GM/WM intensity references (default: 10.0)."},
    {"-erosion", ARGV_FLOAT, (char *)TRUE, (char *)&erosion_mm,
     "Erosion (mm) of the WM intensity control set (default: 3.0)."},
    {"-t-lo", ARGV_FLOAT, (char *)TRUE, (char *)&t_lo,
     "Lower level crossing of the normalized profile (default: 0.25)."},
    {"-t-hi", ARGV_FLOAT, (char *)TRUE, (char *)&t_hi,
     "Upper level crossing (default: 0.75). The distance between the two\n"
     "is the transition width, which is the actual observable."},
    {"-search", ARGV_FLOAT, (char *)TRUE, (char *)&search_mm,
     "Half-length (mm) of the profile searched along the normal (default: 4.0)."},
    {"-step", ARGV_FLOAT, (char *)TRUE, (char *)&step_mm,
     "Sampling step (mm) along the profile (default: 0.25)."},
    {"-width-pct", ARGV_FLOAT, (char *)TRUE, (char *)&width_pct,
     "Percentile of the measured widths used as this brain's healthy\n"
     "reference (default: 25.0)."},
    {"-gain", ARGV_FLOAT, (char *)TRUE, (char *)&gain,
     "Displacement per unit excess width (default: 0.5)."},
    {"-max-offset", ARGV_FLOAT, (char *)TRUE, (char *)&max_offset_mm,
     "Clamp on the displacement in mm (default: 1.5)."},
    {"-smooth", ARGV_FLOAT, (char *)TRUE, (char *)&smooth_fwhm,
     "FWHM (mm) of smoothing within the boundary sheet (default: 8.0).\n"
     "0 disables it and shows the raw per-voxel estimate."},
    {"-width-out", ARGV_STRING, (char *)1, (char *)&width_file,
     "Also write the raw transition width in mm to this volume. This is the\n"
     "map to look at first: the displacement is only a rescaling of it."},
    {"-verbose", ARGV_CONSTANT, (char *)TRUE, (char *)&verbose,
     "Enable verbose output."},
    {NULL, ARGV_END, NULL, NULL, NULL}};

static void
usage(const char *executable)
{
    fprintf(stderr,
            "\nUsage: %s [options] label.nii t1w.nii offset.nii\n\n"
            "Estimate how far myelination has displaced the GM/WM boundary,\n"
            "in millimetres, for use as an additive correction to PBT's WM\n"
            "distance map.\n\n"
            "In the motor strip and along the line of Gennari in V1 the deep\n"
            "cortical layers approach WM intensity, so the boundary is placed\n"
            "too far out and the cortex comes back too thin. The label map is\n"
            "derived from the intensity, so the two boundaries agree by\n"
            "construction and their difference measures nothing; what does\n"
            "distinguish myelinated cortex is the *width* of the intensity\n"
            "transition across the boundary. Healthy cortex steps from grey to\n"
            "white over about the partial-volume width, a myelinated deep layer\n"
            "turns that step into a ramp one to three millimetres wide, and the\n"
            "excess over this brain's own healthy majority is the signal.\n\n",
            executable);
}

int main(int argc, char *argv[])
{
    char *label_file, *t1w_file, *out_file;
    float *label_data, *t1w_data, *offset_data, *width_data = NULL;
    nifti_image *label_nii, *t1w_nii;
    int dims[3];
    double voxelsize[3];
    long rc;

    if (ParseArgv(&argc, argv, argTable, 0))
    {
        usage(argv[0]);
        exit(EXIT_FAILURE);
    }

    initialize_argument_processing(argc, argv);

    if (!get_string_argument(NULL, &label_file) ||
        !get_string_argument(NULL, &t1w_file) ||
        !get_string_argument(NULL, &out_file))
    {
        usage(argv[0]);
        exit(EXIT_FAILURE);
    }

    label_nii = read_nifti_float(label_file, &label_data, 0);
    if (!label_nii)
    {
        fprintf(stderr, "Error reading label volume: %s\n", label_file);
        return EXIT_FAILURE;
    }

    t1w_nii = read_nifti_float(t1w_file, &t1w_data, 0);
    if (!t1w_nii)
    {
        fprintf(stderr, "Error reading T1w volume: %s\n", t1w_file);
        return EXIT_FAILURE;
    }

    if (label_nii->nx != t1w_nii->nx ||
        label_nii->ny != t1w_nii->ny ||
        label_nii->nz != t1w_nii->nz)
    {
        fprintf(stderr, "Error: label and T1w volumes have different "
                        "dimensions.\n");
        return EXIT_FAILURE;
    }

    dims[0] = label_nii->nx;
    dims[1] = label_nii->ny;
    dims[2] = label_nii->nz;
    voxelsize[0] = label_nii->dx;
    voxelsize[1] = label_nii->dy;
    voxelsize[2] = label_nii->dz;

    offset_data = (float *)calloc((size_t)dims[0] * dims[1] * dims[2],
                                  sizeof(float));
    if (width_file)
        width_data = (float *)calloc((size_t)dims[0] * dims[1] * dims[2],
                                     sizeof(float));
    if (!offset_data || (width_file && !width_data))
    {
        fprintf(stderr, "Memory allocation error.\n");
        return EXIT_FAILURE;
    }

    CAT_BoundaryOffsetOpts opts;
    CAT_BoundaryOffsetOptionsInit(&opts);
    opts.ref_fwhm = ref_fwhm;
    opts.erosion_mm = erosion_mm;
    opts.t_lo = t_lo;
    opts.t_hi = t_hi;
    opts.search_mm = search_mm;
    opts.step_mm = step_mm;
    opts.width_pct = width_pct;
    opts.gain = gain;
    opts.max_offset_mm = max_offset_mm;
    opts.smooth_fwhm = smooth_fwhm;
    opts.verbose = verbose;

    rc = CAT_VolBoundaryOffset(label_data, t1w_data, offset_data, width_data,
                               dims, voxelsize, &opts);
    if (rc < 0)
    {
        fprintf(stderr, "Error during boundary-offset estimation.\n");
        return EXIT_FAILURE;
    }

    nifti_image *out_nii = nifti_copy_nim_info(label_nii);
    if (!write_nifti_float(out_file, offset_data, DT_FLOAT32, 1.0,
                           dims, voxelsize, out_nii))
    {
        fprintf(stderr, "Error writing output: %s\n", out_file);
        return EXIT_FAILURE;
    }

    if (width_data)
    {
        nifti_image *w_nii = nifti_copy_nim_info(label_nii);
        if (!write_nifti_float(width_file, width_data, DT_FLOAT32, 1.0,
                               dims, voxelsize, w_nii))
        {
            fprintf(stderr, "Error writing width map: %s\n", width_file);
            return EXIT_FAILURE;
        }
    }

    if (verbose)
        fprintf(stderr, "Boundary offset written to %s\n", out_file);

    free(label_data);
    free(t1w_data);
    free(offset_data);
    free(width_data);

    return EXIT_SUCCESS;
}
