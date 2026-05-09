/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * CAT_SurfBBReg — Boundary-Based Registration of a volume to a surface.
 *
 * Registers a moving NIfTI volume to one or two WM surfaces (left and/or right
 * hemisphere) using the BBR cost function (Greve & Fischl 2009, NeuroImage
 * 48:63-72).  Outputs a plain 4×4 RAS-to-RAS transform matrix compatible
 * with FSL (flirt -init) and ANTs.
 *
 * Copyright Christian Gaser, University of Jena.
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ParseArgv.h>
#include <bicpl.h>

#include "CAT_NiftiLib.h"
#include "CAT_SurfaceIO.h"
#include "CAT_BBReg.h"
#include "CAT_VolumeReg.h"
#include "CAT_Vol.h"

/* -----------------------------------------------------------------------
 * Argument defaults
 * ----------------------------------------------------------------------- */
static double wm_dist = 2.0;         /* mm inside WM                        */
static double gm_dist = 0.5;         /* mm inside GM (absolute, fallback)   */
static double gm_proj_frac = 0.0;    /* fraction of thickness for GM offset */
static double slope = 0.5;           /* BBR cost slope                      */
static int invert_contrast = -1;     /* -1=auto-detect, 0=T1/FLAIR, 1=T2/BOLD */
static double grid_range_mm = 4.0;   /* Stage-1 translation range (mm)      */
static double grid_range_rad = 0.07; /* Stage-1 rotation range (rad, ~4°)   */
static int grid_steps = 2;           /* Stage-1 steps per DOF               */
static int max_iter = 200;           /* Powell iterations                   */
static double tol = 1e-5;            /* Powell convergence tolerance        */
static int verbose = 0;
static double fwhm = 0.0;               /* Gaussian pre-smoothing FWHM (mm, 0=off) */
static char *output_volume_file = NULL; /* optional coregistered volume output */

/* Surface / label / thickness file names (NULL = not provided) */
static char *lh_surf_file = NULL;
static char *rh_surf_file = NULL;
static char *lh_label_file = NULL;
static char *rh_label_file = NULL;
static char *lh_thick_file = NULL;
static char *rh_thick_file = NULL;

/* Initial parameters (identity transform) */
static double init_tx = 0.0, init_ty = 0.0, init_tz = 0.0;
static double init_rx = 0.0, init_ry = 0.0, init_rz = 0.0;

/* Volume registration (robust init) options */
static int    vol_reg_levels  = 4;    /* Gaussian pyramid levels          */
static int    vol_reg_bins    = 64;   /* NMI histogram bins per axis      */
static double vol_reg_sat     = 4.685;/* Tukey saturation (ROB mode only) */
static int    vol_reg_iter    = 30;   /* Powell/IRLS iterations per level */
static int    no_vol_reg      = 0;    /* skip volume registration init    */
static int    vol_reg_rob     = 0;    /* use IRLS/ROB instead of NMI      */
static char  *ref_vol_file    = NULL; /* separate reference volume (T1w)  */

/* -----------------------------------------------------------------------
 * Argument table
 * ----------------------------------------------------------------------- */
static ArgvInfo argTable[] =
    {
        {"-lh", ARGV_STRING, (char *)1, (char *)&lh_surf_file,
         "Left-hemisphere white-matter surface file."},

        {"-rh", ARGV_STRING, (char *)1, (char *)&rh_surf_file,
         "Right-hemisphere white-matter surface file."},

        {"-lh-label", ARGV_STRING, (char *)1, (char *)&lh_label_file,
         "Left-hemisphere cortex label (scalar .txt/.mnc/.func.mgh). Vertices"
         " with value <= 0.5 are excluded from the cost. Default: all vertices."},

        {"-rh-label", ARGV_STRING, (char *)1, (char *)&rh_label_file,
         "Right-hemisphere cortex label. Default: all vertices."},

        {"-lh-thickness", ARGV_STRING, (char *)1, (char *)&lh_thick_file,
         "Left-hemisphere cortical thickness scalar file (mm).  Used with"
         " -gm-proj-frac."},

        {"-rh-thickness", ARGV_STRING, (char *)1, (char *)&rh_thick_file,
         "Right-hemisphere cortical thickness scalar file (mm).  Used with"
         " -gm-proj-frac."},

        {"-gm-proj-frac", ARGV_FLOAT, (char *)1, (char *)&gm_proj_frac,
         "Fraction of cortical thickness to use as GM sampling offset."
         " Requires -lh-thickness / -rh-thickness.  0 = use -gm-dist. Default: 0"},

        {"-wm-dist", ARGV_FLOAT, (char *)1, (char *)&wm_dist,
         "Sampling distance inside white matter along inward surface normal (mm). Default: 2.0"},

        {"-gm-dist", ARGV_FLOAT, (char *)1, (char *)&gm_dist,
         "Absolute GM sampling offset (mm); used when -gm-proj-frac is 0 or"
         " no thickness file is given. Default: 0.5"},

        {"-slope", ARGV_FLOAT, (char *)1, (char *)&slope,
         "Saturation slope of the BBR cost function. Default: 0.5"},

        {"-t2", ARGV_CONSTANT, (char *)1, (char *)&invert_contrast,
         "Force T2/BOLD contrast (WM darker than GM). Skips auto-detection."
         " Default: auto-detect from surface WM/GM intensities."},

        {"-t1", ARGV_CONSTANT, (char *)0, (char *)&invert_contrast,
         "Force T1/FLAIR contrast (WM brighter than GM). Skips auto-detection."
         " Default: auto-detect from surface WM/GM intensities."},

        {"-grid-range-mm", ARGV_FLOAT, (char *)1, (char *)&grid_range_mm,
         "Half-width of Stage-1 brute-force translation search (mm). Default: 4.0"},

        {"-grid-range-rad", ARGV_FLOAT, (char *)1, (char *)&grid_range_rad,
         "Half-width of Stage-1 brute-force rotation search (radians). Default: 0.07"},

        {"-grid-steps", ARGV_INT, (char *)1, (char *)&grid_steps,
         "Number of grid steps per translation DOF in Stage-1 search. Default: 2"},

        {"-max-iter", ARGV_INT, (char *)1, (char *)&max_iter,
         "Maximum Powell iterations in Stage-2 optimisation. Default: 200"},

        {"-tol", ARGV_FLOAT, (char *)1, (char *)&tol,
         "Fractional convergence tolerance for Powell optimizer. Default: 1e-5"},

        {"-init-tx", ARGV_FLOAT, (char *)1, (char *)&init_tx,
         "Initial x-translation (mm). Default: 0"},
        {"-init-ty", ARGV_FLOAT, (char *)1, (char *)&init_ty,
         "Initial y-translation (mm). Default: 0"},
        {"-init-tz", ARGV_FLOAT, (char *)1, (char *)&init_tz,
         "Initial z-translation (mm). Default: 0"},
        {"-init-rx", ARGV_FLOAT, (char *)1, (char *)&init_rx,
         "Initial x-rotation (radians). Default: 0"},
        {"-init-ry", ARGV_FLOAT, (char *)1, (char *)&init_ry,
         "Initial y-rotation (radians). Default: 0"},
        {"-init-rz", ARGV_FLOAT, (char *)1, (char *)&init_rz,
         "Initial z-rotation (radians). Default: 0"},

        {"-fwhm", ARGV_FLOAT, (char *)1, (char *)&fwhm,
         "Gaussian smoothing FWHM applied to the moving volume before computing the"
         " BBR cost (mm).  Equivalent to FreeSurfer mri_segreg --fwhm.  Each"
         " timepoint of a 4D series is smoothed independently.  Default: 0 (off)."},

        {"-verbose", ARGV_CONSTANT, (char *)1, (char *)&verbose,
         "Print optimisation progress. Default: off"},

        {"-output-volume", ARGV_STRING, (char *)1, (char *)&output_volume_file,
         "Output filename for the coregistered volume.  Only the NIfTI sform/qform"
         " header is updated to reflect the registration transform;"
         " voxel data is not resampled."},

        /* ----- robust volume registration (init) options ----- */
        {"-ref", ARGV_STRING, (char *)1, (char *)&ref_vol_file,
         "Reference (fixed) volume for NMI-based volume registration initialisation"
         " (e.g. T1w anatomical). When surfaces are provided this initialises BBR;"
         " without surfaces only volume registration is performed."},

        {"-vol-reg-levels", ARGV_INT, (char *)1, (char *)&vol_reg_levels,
         "Number of Gaussian pyramid levels for volume registration (default: 4)."},

        {"-vol-reg-bins", ARGV_INT, (char *)1, (char *)&vol_reg_bins,
         "Number of histogram bins per axis for NMI registration (default: 64)."},

        {"-vol-reg-iter", ARGV_INT, (char *)1, (char *)&vol_reg_iter,
         "Maximum Powell/IRLS iterations per pyramid level (default: 30)."},

        {"-vol-reg-sat", ARGV_FLOAT, (char *)1, (char *)&vol_reg_sat,
         "Tukey saturation multiplier used only with -vol-reg-rob (default: 4.685)."},

        {"-vol-reg-rob", ARGV_CONSTANT, (char *)1, (char *)&vol_reg_rob,
         "Use same-modality IRLS/ROB cost instead of NMI for volume registration."
         " Use this for T1w-to-T1w or EPI-to-EPI alignment."},

        {"-no-vol-reg", ARGV_CONSTANT, (char *)1, (char *)&no_vol_reg,
         "Skip volume registration initialisation entirely; start BBR from"
         " -init-tx/ty/tz/rx/ry/rz values (or identity)."},

        {NULL, ARGV_END, NULL, NULL, NULL}};

/* -----------------------------------------------------------------------
 * Usage
 * ----------------------------------------------------------------------- */
static void
usage(const char *exe)
{
    fprintf(stdout,
            "\nUsage: %s [options] volume_file output_matrix_file\n\n"
            "  Boundary-Based Registration (BBR) of a NIfTI volume to one or both\n"
            "  white-matter surfaces.  Surfaces must be in the same RAS coordinate\n"
            "  space as the volume.  The output is a plain-text 4x4 RAS-to-RAS rigid\n"
            "  transform matrix compatible with FSL (flirt -init) and ANTs.\n\n"
            "  At least one surface must be supplied via -lh and/or -rh.\n\n"
            "  Cortex label files (-lh-label / -rh-label) should contain one float\n"
            "  value per vertex; vertices with value <= 0.5 are excluded from the\n"
            "  cost (equivalent to FreeSurfer ?h.cortex.label masking).\n\n"
            "  Thickness files (-lh-thickness / -rh-thickness) enable fractional GM\n"
            "  projection via -gm-proj-frac: d_gm = frac * thickness[i] per vertex.\n\n"
            "  Contrast type (T1/FLAIR vs T2/BOLD) is auto-detected from the WM/GM\n"
            "  intensity ratio at the initial position.  Override with -t1 or -t2.\n\n"
            "Options:\n",
            exe);
    fprintf(stdout,
            "  -lh <file>              Left-hemisphere WM surface\n"
            "  -rh <file>              Right-hemisphere WM surface\n"
            "  -lh-label <file>        Left cortex label (per-vertex float, 0/1)\n"
            "  -rh-label <file>        Right cortex label (per-vertex float, 0/1)\n"
            "  -lh-thickness <file>    Left cortical thickness (mm, per-vertex)\n"
            "  -rh-thickness <file>    Right cortical thickness (mm, per-vertex)\n"
            "  -gm-proj-frac <frac>    GM offset = frac * thickness (default: 0)\n"
            "  -wm-dist <mm>           WM sampling offset (default: 2.0)\n"
            "  -gm-dist <mm>           Absolute GM sampling offset (default: 0.5)\n"
            "  -slope <val>            BBR cost saturation slope (default: 0.5)\n"
            "  -t1                     Force T1/FLAIR contrast (WM > GM); skip auto-detect\n"
            "  -t2                     Force T2/BOLD  contrast (WM < GM); skip auto-detect\n"
            "  -grid-range-mm <mm>     Stage-1 translation search range (default: 4.0)\n"
            "  -grid-range-rad <rad>   Stage-1 rotation search range (default: 0.07)\n"
            "  -grid-steps <n>         Stage-1 steps per DOF (default: 2)\n"
            "  -max-iter <n>           Powell max iterations (default: 200)\n"
            "  -tol <val>              Powell convergence tolerance (default: 1e-5)\n"
            "  -init-tx/ty/tz <mm>     Initial translation (default: 0)\n"
            "  -init-rx/ry/rz <rad>    Initial rotation (default: 0)\n"
            "  -fwhm <mm>              Gaussian pre-smoothing FWHM (default: 0 = off)\n"
            "  -verbose                Print optimisation progress\n"
            "  -output-volume <file>   Save coregistered volume (header update only)\n"
            "\n"
            "NMI volume registration (initialisation for BBR):\n"
            "  -ref <file>             T1w reference volume (fixed). Without surfaces,\n"
            "                          only volume registration is performed (no BBR).\n"
            "  -vol-reg-levels <n>     Gaussian pyramid levels (default: 4)\n"
            "  -vol-reg-bins <n>       NMI histogram bins per axis (default: 64)\n"
            "  -vol-reg-iter <n>       Powell iterations per level (default: 30)\n"
            "  -vol-reg-rob            Use IRLS/ROB cost (same-modality) instead of NMI\n"
            "  -vol-reg-sat <val>      Tukey multiplier for -vol-reg-rob (default: 4.685)\n"
            "  -no-vol-reg             Skip volume-registration initialisation\n\n");
}

/* -----------------------------------------------------------------------
 * Helper: load a surface from file
 * ----------------------------------------------------------------------- */
static polygons_struct *
load_surface(const char *filename)
{
    int n_objects;
    File_formats format;
    object_struct **objects;

    if (input_graphics_any_format((char *)filename, &format, &n_objects, &objects) != OK || n_objects < 1 || get_object_type(objects[0]) != POLYGONS)
    {
        fprintf(stderr, "Error reading surface: %s\n", filename);
        return NULL;
    }
    polygons_struct *poly = get_polygons_ptr(objects[0]);
    compute_polygon_normals(poly);
    return poly;
}

/* -----------------------------------------------------------------------
 * Helper: load per-vertex scalar float array
 * ----------------------------------------------------------------------- */
static float *
load_scalars(const char *filename, int expected_n, const char *desc)
{
    int n_vals = 0;
    double *dvals = NULL;
    float *vals;
    int i;

    if (input_values_any_format((char *)filename, &n_vals, &dvals) != OK)
    {
        fprintf(stderr, "Error reading %s file: %s\n", desc, filename);
        return NULL;
    }
    if (n_vals != expected_n)
    {
        fprintf(stderr,
                "Warning: %s file has %d values but surface has %d vertices; "
                "ignoring.\n",
                desc, n_vals, expected_n);
        free(dvals);
        return NULL;
    }
    vals = (float *)malloc(n_vals * sizeof(float));
    if (!vals)
    {
        fprintf(stderr, "Out of memory loading %s\n", desc);
        free(dvals);
        return NULL;
    }
    for (i = 0; i < n_vals; i++)
        vals[i] = (float)dvals[i];
    free(dvals);
    return vals;
}

/* -----------------------------------------------------------------------
 * main
 * ----------------------------------------------------------------------- */
int main(int argc, char *argv[])
{
    char *volume_file, *output_file;
    nifti_image *nii_ptr;
    float *vol = NULL;
    int dims[3];
    CAT_RigidParams p_init, p_best;
    double m_best[16]; /* surface(T1)→EPI; internal optimizer direction */
    double m_inv[16];  /* EPI→surface(T1); standard output / header direction */
    double final_cost;

    CAT_SurfData surfs[2];
    int n_surfs = 0;

    polygons_struct *lh_poly = NULL;
    polygons_struct *rh_poly = NULL;
    float *lh_label = NULL;
    float *rh_label = NULL;
    float *lh_thick = NULL;
    float *rh_thick = NULL;

    if (ParseArgv(&argc, argv, argTable, 0) || argc != 3)
    {
        usage(argv[0]);
        exit(EXIT_FAILURE);
    }

    volume_file = argv[1];
    output_file = argv[2];

    /* --- Require either at least one surface OR a reference volume --- */
    if (!lh_surf_file && !rh_surf_file && !ref_vol_file)
    {
        fprintf(stderr,
                "Error: provide at least one WM surface (-lh/-rh) for BBR, or"
                " a reference volume (-ref) for standalone volume registration.\n");
        usage(argv[0]);
        exit(EXIT_FAILURE);
    }

    /* --- Load surfaces and optional per-vertex data --- */
    if (lh_surf_file)
    {
        lh_poly = load_surface(lh_surf_file);
        if (!lh_poly)
            exit(EXIT_FAILURE);

        if (lh_label_file)
            lh_label = load_scalars(lh_label_file, lh_poly->n_points, "lh-label");

        if (lh_thick_file)
            lh_thick = load_scalars(lh_thick_file, lh_poly->n_points, "lh-thickness");

        surfs[n_surfs].surface = lh_poly;
        surfs[n_surfs].cortex_mask = lh_label;
        surfs[n_surfs].thickness = lh_thick;
        surfs[n_surfs].gm_proj_frac = (lh_thick && gm_proj_frac > 0.0)
                                          ? gm_proj_frac
                                          : 0.0;
        n_surfs++;
    }

    if (rh_surf_file)
    {
        rh_poly = load_surface(rh_surf_file);
        if (!rh_poly)
            exit(EXIT_FAILURE);

        if (rh_label_file)
            rh_label = load_scalars(rh_label_file, rh_poly->n_points, "rh-label");

        if (rh_thick_file)
            rh_thick = load_scalars(rh_thick_file, rh_poly->n_points, "rh-thickness");

        surfs[n_surfs].surface = rh_poly;
        surfs[n_surfs].cortex_mask = rh_label;
        surfs[n_surfs].thickness = rh_thick;
        surfs[n_surfs].gm_proj_frac = (rh_thick && gm_proj_frac > 0.0)
                                          ? gm_proj_frac
                                          : 0.0;
        n_surfs++;
    }

    /* --- Load volume --- */
    nii_ptr = read_nifti_float(volume_file, &vol, 0);
    if (!nii_ptr || !vol)
    {
        fprintf(stderr, "Error reading volume: %s\n", volume_file);
        exit(EXIT_FAILURE);
    }
    dims[0] = nii_ptr->nx;
    dims[1] = nii_ptr->ny;
    dims[2] = nii_ptr->nz;

    /* --- Optional Gaussian pre-smoothing (applied per timepoint) --- */
    if (fwhm > 0.0)
    {
        double voxelsize[3];
        double fwhm_vox[3];
        int t, nvox3 = dims[0] * dims[1] * dims[2];
        int nt = (nii_ptr->nt > 1) ? nii_ptr->nt : 1;

        voxelsize[0] = nii_ptr->dx > 0.0 ? nii_ptr->dx : 1.0;
        voxelsize[1] = nii_ptr->dy > 0.0 ? nii_ptr->dy : 1.0;
        voxelsize[2] = nii_ptr->dz > 0.0 ? nii_ptr->dz : 1.0;

        /* smooth3 takes FWHM in voxel units */
        fwhm_vox[0] = fwhm / voxelsize[0];
        fwhm_vox[1] = fwhm / voxelsize[1];
        fwhm_vox[2] = fwhm / voxelsize[2];

        if (verbose)
            printf("Smoothing: FWHM %.2f mm (%.2f x %.2f x %.2f vox), %d frame(s)\n",
                   fwhm, fwhm_vox[0], fwhm_vox[1], fwhm_vox[2], nt);

        for (t = 0; t < nt; t++)
            smooth3(vol + t * nvox3, dims, voxelsize, fwhm_vox, 0, DT_FLOAT32);
    }

    /* --- Robust volume registration initialisation --- */
    p_init.tx = init_tx;
    p_init.ty = init_ty;
    p_init.tz = init_tz;
    p_init.rx = init_rx;
    p_init.ry = init_ry;
    p_init.rz = init_rz;

    if (!no_vol_reg && ref_vol_file)
    {
        /* Load the reference (T1w) volume */
        float *ref_vol = NULL;
        nifti_image *ref_nii = read_nifti_float(ref_vol_file, &ref_vol, 0);
        if (!ref_nii || !ref_vol)
        {
            fprintf(stderr, "Error reading reference volume: %s\n", ref_vol_file);
            free(vol);
            nii_ptr->data = NULL;
            nifti_image_free(nii_ptr);
            exit(EXIT_FAILURE);
        }
        int ref_dims[3];
        ref_dims[0] = ref_nii->nx;
        ref_dims[1] = ref_nii->ny;
        ref_dims[2] = ref_nii->nz;

        if (verbose)
            printf("VolumeReg init (%s): fixed=%s (%dx%dx%d)  moving=%s"
                   " (%dx%dx%d)\n",
                   vol_reg_rob ? "ROB/IRLS" : "NMI",
                   ref_vol_file, ref_dims[0], ref_dims[1], ref_dims[2],
                   volume_file, dims[0], dims[1], dims[2]);

        double vol_res;
        if (vol_reg_rob)
        {
            /* Same-modality: linearised intensity matching + Tukey biweight */
            vol_res = CAT_VolumeReg_register(
                ref_vol,  ref_nii,  ref_dims,
                vol,      nii_ptr,  dims,
                &p_init,
                vol_reg_levels, vol_reg_sat, vol_reg_iter, verbose);
            if (verbose)
                printf("VolumeReg done (ROB): residual=%.4f  "
                       "tx=%.2f ty=%.2f tz=%.2f  rx=%.4f ry=%.4f rz=%.4f\n",
                       vol_res,
                       p_init.tx, p_init.ty, p_init.tz,
                       p_init.rx, p_init.ry, p_init.rz);
        }
        else
        {
            /* Cross-modal (default): Normalised Mutual Information */
            vol_res = CAT_VolumeReg_register_NMI(
                ref_vol,  ref_nii,  ref_dims,
                vol,      nii_ptr,  dims,
                &p_init,
                vol_reg_levels, vol_reg_bins, vol_reg_iter, verbose);
            if (verbose)
                printf("VolumeReg done (NMI=%.4f):  "
                       "tx=%.2f ty=%.2f tz=%.2f  rx=%.4f ry=%.4f rz=%.4f\n",
                       vol_res,
                       p_init.tx, p_init.ty, p_init.tz,
                       p_init.rx, p_init.ry, p_init.rz);
        }

        free(ref_vol);
        ref_nii->data = NULL;
        nifti_image_free(ref_nii);

        /* If no surfaces are provided, write the volume-registration result
         * directly and exit (standalone mode). */
        if (n_surfs == 0)
        {
            double m_vol[16], m_vol_inv[16];
            CAT_BBReg_params_to_matrix(&p_init, m_vol);
            CAT_BBReg_invert_matrix(m_vol, m_vol_inv);

            if (CAT_BBReg_write_matrix(output_file, m_vol_inv) != 0)
            {
                fprintf(stderr, "Error writing output matrix: %s\n", output_file);
                free(vol);
                nii_ptr->data = NULL;
                nifti_image_free(nii_ptr);
                exit(EXIT_FAILURE);
            }
            if (vol_reg_rob)
                printf("VolumeReg done (standalone ROB).  Residual: %.6f\n",
                       vol_res);
            else
                printf("VolumeReg done (standalone NMI).  NMI: %.6f\n",
                       vol_res);
            printf("Transform written to: %s\n", output_file);
            printf("Parameters (ref->moving, internal): "
                   "tx=%.3f ty=%.3f tz=%.3f  rx=%.4f ry=%.4f rz=%.4f\n",
                   p_init.tx, p_init.ty, p_init.tz,
                   p_init.rx, p_init.ry, p_init.rz);
            free(vol);
            nii_ptr->data = NULL;
            nifti_image_free(nii_ptr);
            return EXIT_SUCCESS;
        }
    }
    else if (!no_vol_reg && n_surfs == 0)
    {
        fprintf(stderr,
                "Error: no surfaces supplied (-lh/-rh) and no reference volume"
                " (-ref) given.\n"
                "Provide at least one surface for BBR, or supply -ref for"
                " standalone volume registration.\n");
        exit(EXIT_FAILURE);
    }

    if (verbose)
    {
        int total_verts = 0;
        int s;
        for (s = 0; s < n_surfs; s++)
            total_verts += surfs[s].surface->n_points;

        printf("Surfaces: %d hemisphere(s), %d total vertices\n",
               n_surfs, total_verts);
        if (lh_surf_file)
            printf("  LH: %s (%d verts%s%s)\n", lh_surf_file,
                   lh_poly->n_points,
                   lh_label ? ", label masked" : "",
                   lh_thick ? ", thickness proj" : "");
        if (rh_surf_file)
            printf("  RH: %s (%d verts%s%s)\n", rh_surf_file,
                   rh_poly->n_points,
                   rh_label ? ", label masked" : "",
                   rh_thick ? ", thickness proj" : "");
        printf("Volume:  %s (%dx%dx%d)\n",
               volume_file, dims[0], dims[1], dims[2]);
        printf("WM dist: %.2f mm  ", wm_dist);
        if (gm_proj_frac > 0.0)
            printf("GM frac: %.2f of thickness", gm_proj_frac);
        else
            printf("GM dist: %.2f mm", gm_dist);
        printf("  slope: %.3f  %s\n", slope,
               invert_contrast == 1 ? "T2/BOLD (forced)"  :
               invert_contrast == 0 ? "T1/FLAIR (forced)" : "auto-detect");
        printf("Stage-1 grid: +/-%.1f mm, +/-%.3f rad, %d steps per DOF\n",
               grid_range_mm, grid_range_rad, grid_steps);
    }

    /* --- Auto-detect contrast type if not forced by -t1 / -t2 --- */
    if (invert_contrast < 0)
    {
        int detected = CAT_BBReg_detect_contrast(&p_init,
                                                 surfs, n_surfs,
                                                 vol, nii_ptr, dims,
                                                 wm_dist, gm_dist,
                                                 verbose);
        invert_contrast = (detected >= 0) ? detected : 0;
        if (detected < 0 && verbose)
            printf("Contrast auto-detection inconclusive — defaulting to T1.\n");
    }

    /* --- Optimise --- */
    final_cost = CAT_BBReg_optimise(&p_init, &p_best,
                                    surfs, n_surfs,
                                    vol, nii_ptr, dims,
                                    wm_dist, gm_dist, slope, invert_contrast,
                                    grid_range_mm, grid_range_rad, grid_steps,
                                    max_iter, tol, verbose);

    /* --- Convert best parameters to 4x4 matrix and write ---
     * m_best maps fixed-surface(T1) RAS → moving-volume(EPI) RAS; that is the
     * direction the optimizer uses to transform surface sample points into the
     * EPI before trilinear interpolation.  The standard output convention
     * (matching bbregister / bbreg) is the *inverse*: EPI RAS → T1 RAS.
     * Applying the inverse to the original sform (voxel→EPI RAS) then gives
     * the updated sform (voxel→T1 RAS). */
    CAT_BBReg_params_to_matrix(&p_best, m_best);
    CAT_BBReg_invert_matrix(m_best, m_inv);

    if (CAT_BBReg_write_matrix(output_file, m_inv) != 0)
    {
        fprintf(stderr, "Error writing output matrix: %s\n", output_file);
        free(vol);
        nii_ptr->data = NULL;
        nifti_image_free(nii_ptr);
        exit(EXIT_FAILURE);
    }

    printf("BBReg done.  Final cost: %.6f\n", final_cost);
    printf("Transform written to: %s\n", output_file);
    printf("Parameters (T1->EPI, internal): tx=%.3f ty=%.3f tz=%.3f  rx=%.4f ry=%.4f rz=%.4f\n",
           p_best.tx, p_best.ty, p_best.tz,
           p_best.rx, p_best.ry, p_best.rz);

    /* --- Optionally write coregistered volume (header update only) --- */
    if (output_volume_file)
    {
        int i, j, k;
        mat44 new_sto;
        float qb, qc, qd, qx, qy, qz, dx, dy, dz, qfac;

        /* new_sto_xyz = m_inv * old_sto_xyz
         * m_inv maps EPI RAS → T1/surface RAS (standard output convention).
         * Composing it with the original sform (voxel → EPI RAS) gives the
         * updated sform (voxel → T1/surface RAS). */
        for (i = 0; i < 4; i++)
            for (j = 0; j < 4; j++)
            {
                double sum = 0.0;
                for (k = 0; k < 4; k++)
                    sum += m_inv[i * 4 + k] * (double)nii_ptr->sto_xyz.m[k][j];
                new_sto.m[i][j] = (float)sum;
            }
        nii_ptr->sto_xyz = new_sto;

        /* Update qform with the same transform if it was set */
        if (nii_ptr->qform_code > 0)
        {
            mat44 new_qto;
            for (i = 0; i < 4; i++)
                for (j = 0; j < 4; j++)
                {
                    double sum = 0.0;
                    for (k = 0; k < 4; k++)
                        sum += m_inv[i * 4 + k] * (double)nii_ptr->qto_xyz.m[k][j];
                    new_qto.m[i][j] = (float)sum;
                }
            nii_ptr->qto_xyz = new_qto;

            /* Recompute quaternion scalar fields from the updated matrix */
            nifti_mat44_to_quatern(new_qto,
                                   &qb, &qc, &qd,
                                   &qx, &qy, &qz,
                                   &dx, &dy, &dz, &qfac);
            nii_ptr->quatern_b = qb;
            nii_ptr->quatern_c = qc;
            nii_ptr->quatern_d = qd;
            nii_ptr->qoffset_x = qx;
            nii_ptr->qoffset_y = qy;
            nii_ptr->qoffset_z = qz;
            nii_ptr->qfac = qfac;
        }

        /* Write original (float-converted) voxel data with updated header.
         * read_nifti_float() divided nvox by nt to give the 3D slice count, but
         * vol[] holds all nt timepoints.  Restore the full 4D voxel count so
         * nifti_image_write() writes every timepoint, not just the first. */
        if (nii_ptr->nt > 1)
            nii_ptr->nvox = (size_t)dims[0] * (size_t)dims[1] *
                            (size_t)dims[2] * (size_t)nii_ptr->nt;

        nii_ptr->datatype = DT_FLOAT32;
        nii_ptr->nbyper = sizeof(float);
        nii_ptr->scl_slope = 0.0;
        nii_ptr->scl_inter = 0.0;
        nii_ptr->data = (void *)vol;

        if (nifti_set_filenames(nii_ptr, output_volume_file, 0, 1) != 0)
            fprintf(stderr, "Error setting output volume filename: %s\n",
                    output_volume_file);
        else
        {
            nifti_image_write(nii_ptr);
            printf("Coregistered volume written to: %s\n", output_volume_file);
        }

        nii_ptr->data = NULL; /* prevent double-free; vol is freed below */
    }

    free(vol);
    nii_ptr->data = NULL;
    nifti_image_free(nii_ptr);
    if (lh_label)
        free(lh_label);
    if (rh_label)
        free(rh_label);
    if (lh_thick)
        free(lh_thick);
    if (rh_thick)
        free(rh_thick);
    return EXIT_SUCCESS;
}
