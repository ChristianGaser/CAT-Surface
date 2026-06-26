/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 *
 * Thin CLI over the CAT_WarpDemons library module. Parses arguments, loads the
 * source/template surfaces and their spheres, and runs demon-based spherical
 * registration (default: 2-stage Spherical Demons, sulcal depth -> mean
 * curvature). All registration logic lives in Lib/CAT_WarpDemons.c.
 */

#include <bicpl.h>
#include <ParseArgv.h>

#include "CAT_WarpDemons.h"
#include "CAT_SurfaceIO.h"

/* defaults (kept in sync with CAT_WarpDemonsDefaults) */
char  *src_file            = NULL;
char  *src_sphere_file     = NULL;
char  *trg_file            = NULL;
char  *trg_sphere_file     = NULL;
char  *output_surface_file = NULL;
char  *output_sphere_file  = NULL;

int    n_points   = 20480;  /* finest pyramid level; coarser levels are 1/4 each */
int    rotate     = 1;
int    curvtype0  = 5;   /* level 0: sulcal-depth-like (coarsest) */
int    curvtype1  = 5;   /* level 1: sulcal-depth-like */
int    curvtype2  = 5;   /* level 2: sulcal-depth-like (only with -steps 3+) */
int    curvtype3  = 0;   /* level 3: mean curvature (only with -steps 4) */
int    n_steps    = 2;
int    debug      = 0;
int    iters      = 100;
int    verbose    = 0;
double rate       = 0.97;
double fwhm_flow  = 30.0;
double fwhm_curv  = 6.0;
double fwhm_disp  = 10.0;
double max_step_deg = 10.0;
double sigma_x_default = 2.0;  /* SD max_step = 2 */
int    smooth_velocity = 0;    /* SD default: velocity smoothing off */
int    smooth_displacement = 1;
int    use_hessian = 1;
int    use_line_search = 1;
int    use_expmap  = 1;
double step_factor = 1.0;

static ArgvInfo argTable[] = {
  {"-i", ARGV_STRING, (char *) 1, (char *) &src_file,
   "Input file."},
  {"-is", ARGV_STRING, (char *) 1, (char *) &src_sphere_file,
   "Input sphere file."},
  {"-t", ARGV_STRING, (char *) 1, (char *) &trg_file,
   "Template file."},
  {"-ts", ARGV_STRING, (char *) 1, (char *) &trg_sphere_file,
   "Template sphere file."},
  {"-w", ARGV_STRING, (char *) 1, (char *) &output_surface_file,
   "Warped brain."},
  {"-ws", ARGV_STRING, (char *) 1, (char *) &output_sphere_file,
   "Warped input sphere."},
  {"-npoints", ARGV_INT, (char *) 1, (char *) &n_points,
   "Finest pyramid resolution (e.g. 81920); coarser levels use 1/4 the points each."},
  {"-fwhm-flow", ARGV_FLOAT, (char *) 1, (char *) &fwhm_flow,
   "Filter size for velocity update in FWHM."},
  {"-fwhm", ARGV_FLOAT, (char *) 1, (char *) &fwhm_curv,
   "Filter size for curvature map in FWHM."},
  {"-rate", ARGV_FLOAT, (char *) 1, (char *) &rate,
   "Change of fwhm-flow for each iteration."},
  {"-max-step-deg", ARGV_FLOAT, (char *) 1, (char *) &max_step_deg,
   "Clamp per-iteration |dtheta,dphi| to this many degrees (<=0 disables)."},
  {"-fwhm-disp", ARGV_FLOAT, (char *) 1, (char *) &fwhm_disp,
   "Filter size for displacement field smoothing (elastic prior) in FWHM."},
  {"-smooth-velocity", ARGV_CONSTANT, (char *) TRUE, (char *) &smooth_velocity,
   "Enable velocity update smoothing (fluid prior, default OFF as in original SD)."},
  {"-no-smooth-displacement", ARGV_CONSTANT, (char *) FALSE, (char *) &smooth_displacement,
   "Disable displacement field smoothing (elastic prior, default ON)."},
  {"-sigma-x", ARGV_FLOAT, (char *) 1, (char *) &sigma_x_default,
   "Regularization weight sigma_x for Spherical Demons (default 1.0)."},
  {"-step-factor", ARGV_FLOAT, (char *) 1, (char *) &step_factor,
   "Global step size factor (default 1.0)."},
  {"-no-hessian", ARGV_CONSTANT, (char *) FALSE, (char *) &use_hessian,
   "Disable per-vertex Hessian-based update (use gradient only)."},
  {"-no-line-search", ARGV_CONSTANT, (char *) FALSE, (char *) &use_line_search,
   "Disable line search during integration."},
  {"-no-expmap", ARGV_CONSTANT, (char *) FALSE, (char *) &use_expmap,
   "Disable diffeomorphic scaling-and-squaring exponential map; use additive flow."},
  {"-maxiters", ARGV_INT, (char *) 1, (char *) &iters,
   "Maximum number of iterations per stage."},
  {"-steps", ARGV_INT, (char *) 1, (char *) &n_steps,
   "Number of multi-resolution pyramid levels (1-3, coarse to fine)."},
  {"-norot", ARGV_CONSTANT, (char *) FALSE, (char *) &rotate,
   "Don't rotate input surface before warping."},
  {"-type0", ARGV_INT, (char *) 1, (char *) &curvtype0,
   "Curvature type for level 1 (coarsest)\n\t0 - mean curvature (averaged over 3mm, in degrees)\n\t1 - gaussian curvature\n\t2 - curvedness\n\t3 - shape index\n\t4 - mean curvature (in radians)\n\t5 - sulcal depth like estimator."},
  {"-type1", ARGV_INT, (char *) 1, (char *) &curvtype1,
   "Curvature type for level 2 (see -type0 for values)."},
  {"-type2", ARGV_INT, (char *) 1, (char *) &curvtype2,
   "Curvature type for level 3 (see -type0 for values)."},
  {"-type3", ARGV_INT, (char *) 1, (char *) &curvtype3,
   "Curvature type for level 4 (finest; see -type0 for values)."},
  {"-verbose", ARGV_CONSTANT, (char *) TRUE, (char *) &verbose,
   "Be verbose."},
  {"-debug", ARGV_CONSTANT, (char *) TRUE, (char *) &debug,
   "Save debug files."},
   {NULL, ARGV_END, NULL, NULL, NULL}
};

int
main(int argc, char *argv[])
{
    File_formats     format;
    polygons_struct  *src, *trg, *src_sphere, *trg_sphere, *warped_src_sphere;
    int              n_objects;
    object_struct    **objects;
    object_struct    **src_objects;
    CAT_WarpDemonsOptions opt;

    if (ParseArgv(&argc, argv, argTable, 0) ||
        src_file == NULL || trg_file == NULL ||
        src_sphere_file == NULL || trg_sphere_file == NULL ||
        (output_surface_file == NULL && output_sphere_file == NULL)) {
        fprintf(stderr, "\nUsage: %s [options]\n", argv[0]);
        fprintf(stderr, "       %s -help\n\n", argv[0]);
        return EXIT_FAILURE;
    }

    if (input_graphics_any_format(trg_file, &format, &n_objects, &objects) != OK)
        return EXIT_FAILURE;
    if (n_objects != 1 || get_object_type(objects[0]) != POLYGONS) {
        fprintf(stderr, "Template file must contain 1 polygons object.\n");
        return EXIT_FAILURE;
    }
    trg = get_polygons_ptr(objects[0]);

    if (input_graphics_any_format(src_file, &format, &n_objects, &src_objects) != OK)
        return EXIT_FAILURE;
    if (n_objects != 1 || get_object_type(src_objects[0]) != POLYGONS) {
        fprintf(stderr, "Surface file must contain 1 polygons object.\n");
        return EXIT_FAILURE;
    }
    src = get_polygons_ptr(src_objects[0]);

    if (input_graphics_any_format(src_sphere_file, &format, &n_objects, &objects) != OK)
        return EXIT_FAILURE;
    src_sphere = get_polygons_ptr(objects[0]);

    if (input_graphics_any_format(trg_sphere_file, &format, &n_objects, &objects) != OK)
        return EXIT_FAILURE;
    trg_sphere = get_polygons_ptr(objects[0]);

    /* fill options from defaults, then apply parsed overrides */
    CAT_WarpDemonsDefaults(&opt);
    opt.n_points            = n_points;
    opt.n_steps             = n_steps;
    opt.curvtype[0]         = curvtype0;
    opt.curvtype[1]         = curvtype1;
    opt.curvtype[2]         = curvtype2;
    opt.curvtype[3]         = curvtype3;
    opt.iters               = iters;
    opt.rotate              = rotate;
    opt.smooth_velocity     = smooth_velocity;
    opt.smooth_displacement = smooth_displacement;
    opt.use_hessian         = use_hessian;
    opt.use_line_search     = use_line_search;
    opt.use_expmap          = use_expmap;
    opt.fwhm_flow           = fwhm_flow;
    opt.fwhm_curv           = fwhm_curv;
    opt.fwhm_disp           = fwhm_disp;
    opt.rate                = rate;
    opt.max_step_deg        = max_step_deg;
    opt.sigma_x             = sigma_x_default;
    opt.step_factor         = step_factor;
    opt.verbose             = verbose;
    opt.debug               = debug;

    if (opt.n_steps < 1) opt.n_steps = 1;
    if (opt.n_steps > CAT_WARP_DEMONS_MAX_STEPS) opt.n_steps = CAT_WARP_DEMONS_MAX_STEPS;

    /* Build the coarse-to-fine pyramid: finest level = n_points, each coarser
     * level uses 1/4 of the points (next tetrahedral subdivision down). */
    {
        int lvl, np = n_points;
        for (lvl = opt.n_steps - 1; lvl >= 0; lvl--) {
            opt.level_points[lvl] = np;
            np /= 4;
        }
    }

    warped_src_sphere = (polygons_struct *) malloc(sizeof(polygons_struct));

    if (CAT_WarpDemonsRegister(src, src_sphere, trg, trg_sphere,
                               warped_src_sphere, &opt) != OK) {
        fprintf(stderr, "Registration failed.\n");
        return EXIT_FAILURE;
    }

    if (output_sphere_file != NULL) {
        object_struct *out = create_object(POLYGONS);
        *get_polygons_ptr(out) = *warped_src_sphere;
        if (output_graphics_any_format(output_sphere_file, format, 1, &out, NULL) != OK)
            return EXIT_FAILURE;
    }

    if (output_surface_file != NULL) {
        /* The deformed sphere defines the new spherical parameterization; write
         * the input surface geometry alongside it for downstream resampling. */
        if (output_graphics_any_format(output_surface_file, format, 1, src_objects, NULL) != OK)
            return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
