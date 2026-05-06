/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * CAT_SurfBBReg — Boundary-Based Registration of a volume to a surface.
 *
 * Registers a moving NIfTI volume to a WM surface using the BBR cost function
 * (Greve & Fischl 2009, NeuroImage 48:63-72).  Outputs a plain 4×4 RAS-to-RAS
 * transform matrix compatible with FSL (flirt -init) and ANTs.
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

/* -----------------------------------------------------------------------
 * Argument defaults
 * ----------------------------------------------------------------------- */
static double wm_dist = 2.0;         /* mm inside WM */
static double gm_dist = 0.5;         /* mm inside GM */
static double slope = 0.5;           /* BBR cost slope */
static int invert_contrast = 0;      /* 0=T1, 1=T2/BOLD */
static double grid_range_mm = 4.0;   /* Stage-1 translation range (mm) */
static double grid_range_rad = 0.07; /* Stage-1 rotation range (rad, ~4°) */
static int grid_steps = 2;           /* Stage-1 steps per DOF */
static int max_iter = 200;           /* Powell iterations */
static double tol = 1e-5;            /* Powell convergence tolerance */
static int verbose = 0;

/* Initial parameters (identity transform) */
static double init_tx = 0.0, init_ty = 0.0, init_tz = 0.0;
static double init_rx = 0.0, init_ry = 0.0, init_rz = 0.0;

/* -----------------------------------------------------------------------
 * Argument table
 * ----------------------------------------------------------------------- */
static ArgvInfo argTable[] =
    {
        {"-wm-dist", ARGV_FLOAT, (char *)1, (char *)&wm_dist,
         "Sampling distance inside white matter along inward surface normal (mm). Default: 2.0"},

        {"-gm-dist", ARGV_FLOAT, (char *)1, (char *)&gm_dist,
         "Sampling distance inside grey matter along outward surface normal (mm). Default: 0.5"},

        {"-slope", ARGV_FLOAT, (char *)1, (char *)&slope,
         "Saturation slope of the BBR cost function. Default: 0.5"},

        {"-t2", ARGV_CONSTANT, (char *)1, (char *)&invert_contrast,
         "Invert contrast (use for T2-weighted or BOLD EPI input). Default: off (T1)"},

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

        {"-verbose", ARGV_CONSTANT, (char *)1, (char *)&verbose,
         "Print optimisation progress. Default: off"},

        {NULL, ARGV_END, NULL, NULL, NULL}};

static void
usage(const char *exe)
{
    fprintf(stdout,
            "\nUsage: %s [options] surface_file volume_file output_matrix_file\n\n"
            "  Boundary-Based Registration (BBR) of a NIfTI volume to a white-matter\n"
            "  surface.  Both inputs must be in the same RAS coordinate space.\n"
            "  The output is a plain-text 4x4 RAS-to-RAS rigid transform matrix\n"
            "  compatible with FSL (flirt -init) and ANTs.\n\n"
            "  Positive wm-dist moves the WM sample point inward along the inward\n"
            "  normal (into white matter); positive gm-dist moves the GM sample\n"
            "  point outward (into grey matter / CSF).\n\n"
            "  For T1 input WM > GM intensity; for T2/BOLD use -t2 to invert.\n\n"
            "Options:\n",
            exe);
    fprintf(stdout,
            "  -wm-dist <mm>         WM sampling offset (default: 2.0)\n"
            "  -gm-dist <mm>         GM sampling offset (default: 0.5)\n"
            "  -slope <val>          BBR cost saturation slope (default: 0.5)\n"
            "  -t2                   Invert contrast (T2/BOLD input)\n"
            "  -grid-range-mm <mm>   Stage-1 translation search range (default: 4.0)\n"
            "  -grid-range-rad <rad> Stage-1 rotation search range (default: 0.07)\n"
            "  -grid-steps <n>       Stage-1 steps per DOF (default: 2)\n"
            "  -max-iter <n>         Powell max iterations (default: 200)\n"
            "  -tol <val>            Powell convergence tolerance (default: 1e-5)\n"
            "  -init-tx/ty/tz <mm>   Initial translation (default: 0)\n"
            "  -init-rx/ry/rz <rad>  Initial rotation (default: 0)\n"
            "  -verbose              Print optimisation progress\n\n");
}

int main(int argc, char *argv[])
{
    char *surface_file, *volume_file, *output_file;
    int n_objects;
    File_formats format;
    object_struct **objects;
    polygons_struct *polygons;
    nifti_image *nii_ptr;
    float *vol = NULL;
    int dims[3];
    CAT_RigidParams p_init, p_best;
    double m_best[16];
    double final_cost;

    if (ParseArgv(&argc, argv, argTable, 0) || argc != 4)
    {
        usage(argv[0]);
        exit(EXIT_FAILURE);
    }

    surface_file = argv[1];
    volume_file = argv[2];
    output_file = argv[3];

    /* --- Load surface --- */
    if (input_graphics_any_format(surface_file, &format,
                                  &n_objects, &objects) != OK ||
        n_objects < 1 ||
        get_object_type(objects[0]) != POLYGONS)
    {
        fprintf(stderr, "Error reading surface: %s\n", surface_file);
        exit(EXIT_FAILURE);
    }
    polygons = get_polygons_ptr(objects[0]);
    compute_polygon_normals(polygons);

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

    /* --- Initial parameters --- */
    p_init.tx = init_tx;
    p_init.ty = init_ty;
    p_init.tz = init_tz;
    p_init.rx = init_rx;
    p_init.ry = init_ry;
    p_init.rz = init_rz;

    if (verbose)
    {
        printf("Surface: %s (%d vertices)\n", surface_file, polygons->n_points);
        printf("Volume:  %s (%dx%dx%d)\n",
               volume_file, dims[0], dims[1], dims[2]);
        printf("WM dist: %.2f mm  GM dist: %.2f mm  slope: %.3f  %s\n",
               wm_dist, gm_dist, slope,
               invert_contrast ? "T2/BOLD" : "T1");
        printf("Stage-1 grid: ±%.1f mm, ±%.3f rad, %d steps per DOF\n",
               grid_range_mm, grid_range_rad, grid_steps);
    }

    /* --- Optimise --- */
    final_cost = CAT_BBReg_optimise(&p_init, &p_best,
                                    polygons, vol, nii_ptr, dims,
                                    wm_dist, gm_dist, slope, invert_contrast,
                                    grid_range_mm, grid_range_rad, grid_steps,
                                    max_iter, tol, verbose);

    /* --- Convert best parameters to 4×4 matrix and write --- */
    CAT_BBReg_params_to_matrix(&p_best, m_best);

    if (CAT_BBReg_write_matrix(output_file, m_best) != 0)
    {
        fprintf(stderr, "Error writing output matrix: %s\n", output_file);
        free(vol);
        nifti_image_free(nii_ptr);
        exit(EXIT_FAILURE);
    }

    printf("BBReg done.  Final cost: %.6f\n", final_cost);
    printf("Transform written to: %s\n", output_file);
    printf("Parameters: tx=%.3f ty=%.3f tz=%.3f  rx=%.4f ry=%.4f rz=%.4f\n",
           p_best.tx, p_best.ty, p_best.tz,
           p_best.rx, p_best.ry, p_best.rz);

    free(vol);
    nifti_image_free(nii_ptr);
    return EXIT_SUCCESS;
}
