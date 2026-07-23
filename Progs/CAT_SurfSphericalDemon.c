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
char *src_file = NULL;
char *src_sphere_file = NULL;
char *trg_file = NULL;
char *trg_sphere_file = NULL;
char *output_surface_file = NULL;
char *output_sphere_file = NULL;
char *mask_file = NULL;

int n_points = 20480; /* finest pyramid level; coarser levels are 1/4 each */
int rotate = 1;
int curvtype0 = 1000; /* level 0: depth-potential (coarsest) */
int curvtype1 = 250;  /* level 1: depth-potential */
int curvtype2 = 125;  /* level 2: depth-potential (only with -steps 3+) */
int curvtype3 = 15;   /* level 3: depth-potential (only with -steps 4) */
int n_steps = 4;
int debug = 0;
int iters = 150;
int verbose = 0;
double fwhm_flow = 16.0;
double fwhm_curv = 16.0;
double fwhm_disp = 6.0;
double max_step_deg = 50.0;
double sigma_x_default = 20.0; /* SD max_step = 2 */
double l_dist = 0.6;           /* metric-distortion regularizer weight (prototype) */
double coarse_stiffness = 1.5; /* extra coarse-level regularization (1 = off) */
double rot_max_degrees = 64.0; /* rotation search: initial half-span */
double rot_min_degrees = 1.0;  /* rotation search: finest span */
int rot_nangles = 4;           /* rotation search: samples per axis per pass */
int use_unfold = 0;            /* relax folded triangles in the final warp */

static ArgvInfo argTable[] = {
    {"-i", ARGV_STRING, (char *)1, (char *)&src_file,
     "Input file."},
    {"-is", ARGV_STRING, (char *)1, (char *)&src_sphere_file,
     "Input sphere file."},
    {"-t", ARGV_STRING, (char *)1, (char *)&trg_file,
     "Template file."},
    {"-ts", ARGV_STRING, (char *)1, (char *)&trg_sphere_file,
     "Template sphere file."},
    {"-mask", ARGV_STRING, (char *)1, (char *)&mask_file,
     "Per-vertex cortex mask on the TEMPLATE mesh (0 = exclude, e.g. medial wall;\n\t>0 = include). Excludes non-cortex from the data term (FreeSurfer-style).\n\tMust match the template vertex count; resampled to each pyramid level."},
    {"-dist", ARGV_FLOAT, (char *)1, (char *)&l_dist,
     "Metric-distortion regularizer weight (FreeSurfer-style distance term): each\n\titeration takes a gradient step pulling warped neighbour distances back toward\n\tthe original sphere metric, resisting local stretch/fold. 0 = off (default).\n\tTry small values (e.g. 0.05-0.2); independent of -fwhm-disp smoothing."},
    {"-unfold", ARGV_CONSTANT, (char *)TRUE, (char *)&use_unfold,
     "Post-step: relax folded (negative-area) triangles in the final warp until\n\torientations are restored. Removes folds introduced when up-sampling the\n\twarp onto an irregular full-resolution mesh. Default off."},
    {"-w", ARGV_STRING, (char *)1, (char *)&output_surface_file,
     "Warped brain."},
    {"-ws", ARGV_STRING, (char *)1, (char *)&output_sphere_file,
     "Warped input sphere."},
    {"-npoints", ARGV_INT, (char *)1, (char *)&n_points,
     "Finest pyramid resolution (e.g. 81920); coarser levels use 1/4 the points each."},
    {"-fwhm-flow", ARGV_FLOAT, (char *)1, (char *)&fwhm_flow,
     "Filter size for velocity update in FWHM."},
    {"-fwhm", ARGV_FLOAT, (char *)1, (char *)&fwhm_curv,
     "Filter size for curvature map in FWHM."},
    {"-max-step-deg", ARGV_FLOAT, (char *)1, (char *)&max_step_deg,
     "Clamp per-iteration |dtheta,dphi| to this many degrees (<=0 disables)."},
    {"-fwhm-disp", ARGV_FLOAT, (char *)1, (char *)&fwhm_disp,
     "Filter size for displacement field smoothing (elastic prior) in FWHM."},
    {"-sigma-x", ARGV_FLOAT, (char *)1, (char *)&sigma_x_default,
     "Regularization weight sigma_x for Spherical Demons (default 1.0)."},
    {"-stiffness", ARGV_FLOAT, (char *)1, (char *)&coarse_stiffness,
     "Extra Dartel-like stiffness on the coarser pyramid levels: the flow and\n\tdisplacement smoothing FWHM are multiplied by a factor that is this value at\n\tthe coarsest level and decays to 1.0 at the finest. A stiffer coarse warp\n\tmoves whole folds together and resists a sulcus slipping one wavelength into\n\tits neighbour. 1.0 = off (default); try 1.5-2.5."},
    {"-maxiters", ARGV_INT, (char *)1, (char *)&iters,
     "Maximum number of iterations per stage."},
    {"-steps", ARGV_INT, (char *)1, (char *)&n_steps,
     "Number of multi-resolution pyramid levels (1-4, coarse to fine)."},
    {"-norot", ARGV_CONSTANT, (char *)FALSE, (char *)&rotate,
     "Don't rotate input surface before warping."},
    {"-rot-max-deg", ARGV_FLOAT, (char *)1, (char *)&rot_max_degrees,
     "Initial rotation search: half-width of the initial angular span in degrees\n\t(default 64, matching FreeSurfer). This is the capture range."},
    {"-rot-min-deg", ARGV_FLOAT, (char *)1, (char *)&rot_min_degrees,
     "Initial rotation search: stop once the angular span falls below this many\n\tdegrees (default 1)."},
    {"-rot-nangles", ARGV_INT, (char *)1, (char *)&rot_nangles,
     "Initial rotation search: grid samples per axis per pass (default 4). Cost\n\tgrows as (nangles+1)^3 per pass; raise for a denser search, lower for speed."},
    {"-type0", ARGV_INT, (char *)1, (char *)&curvtype0,
     "Curvature type for level 1 (coarsest)\n\t0 - mean curvature (averaged over 3mm, in degrees)\n\t1 - gaussian curvature\n\t2 - curvedness\n\t3 - shape index\n\t4 - mean curvature (in radians)\n\t5 - sulcal depth like estimator\n\t>5 - depth potential with parameter alpha = 1/curvtype."},
    {"-type1", ARGV_INT, (char *)1, (char *)&curvtype1,
     "Curvature type for level 2 (see -type0 for values)."},
    {"-type2", ARGV_INT, (char *)1, (char *)&curvtype2,
     "Curvature type for level 3 (see -type0 for values)."},
    {"-type3", ARGV_INT, (char *)1, (char *)&curvtype3,
     "Curvature type for level 4 (finest; see -type0 for values)."},
    {"-verbose", ARGV_CONSTANT, (char *)TRUE, (char *)&verbose,
     "Be verbose."},
    {"-debug", ARGV_CONSTANT, (char *)TRUE, (char *)&debug,
     "Save debug files."},
    {NULL, ARGV_END, NULL, NULL, NULL}};

int main(int argc, char *argv[])
{
    File_formats format;
    polygons_struct *src, *trg, *src_sphere, *trg_sphere, *warped_src_sphere;
    int n_objects;
    object_struct **objects;
    object_struct **src_objects;
    CAT_WarpDemonsOptions opt;

    if (ParseArgv(&argc, argv, argTable, 0) ||
        src_file == NULL || trg_file == NULL ||
        src_sphere_file == NULL || trg_sphere_file == NULL ||
        (output_surface_file == NULL && output_sphere_file == NULL))
    {
        fprintf(stderr, "\nUsage: %s [options]\n", argv[0]);
        fprintf(stderr, "       %s -help\n\n", argv[0]);
        return EXIT_FAILURE;
    }

    if (input_graphics_any_format(trg_file, &format, &n_objects, &objects) != OK)
        return EXIT_FAILURE;
    if (n_objects != 1 || get_object_type(objects[0]) != POLYGONS)
    {
        fprintf(stderr, "Template file must contain 1 polygons object.\n");
        return EXIT_FAILURE;
    }
    trg = get_polygons_ptr(objects[0]);

    if (input_graphics_any_format(src_file, &format, &n_objects, &src_objects) != OK)
        return EXIT_FAILURE;
    if (n_objects != 1 || get_object_type(src_objects[0]) != POLYGONS)
    {
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
    opt.n_points = n_points;
    opt.n_steps = n_steps;
    opt.curvtype[0] = curvtype0;
    opt.curvtype[1] = curvtype1;
    opt.curvtype[2] = curvtype2;
    opt.curvtype[3] = curvtype3;
    opt.iters = iters;
    opt.rotate = rotate;
    opt.smooth_velocity = 1;
    opt.smooth_displacement = 1;
    opt.use_hessian = 1;
    opt.use_line_search = 0;
    opt.use_expmap = 1;
    opt.use_tangent = 1;
    opt.l_dist = l_dist;
    opt.coarse_stiffness = coarse_stiffness;
    opt.rot_max_degrees = rot_max_degrees;
    opt.rot_min_degrees = rot_min_degrees;
    opt.rot_nangles = rot_nangles;
    opt.geodesic = 1;
    opt.unfold = use_unfold;
    opt.fwhm_flow = fwhm_flow;
    opt.fwhm_curv = fwhm_curv;
    opt.fwhm_disp = fwhm_disp;
    opt.rate = 1.0;
    opt.max_step_deg = max_step_deg;
    opt.sigma_x = sigma_x_default;
    opt.step_factor = 1.0;
    opt.verbose = verbose;
    opt.debug = debug;

    /* Optional template cortex mask (independent of the std map). */
    if (mask_file != NULL)
    {
        int n_mask;
        double *mask_values;
        if (input_values_any_format(mask_file, &n_mask, &mask_values) != OK)
            return EXIT_FAILURE;
        if (n_mask != trg->n_points)
        {
            fprintf(stderr, "Cortex mask has %d values but template has %d points.\n",
                    n_mask, trg->n_points);
            return EXIT_FAILURE;
        }
        opt.cortex_mask = mask_values;
    }

    if (opt.n_steps < 1)
        opt.n_steps = 1;
    if (opt.n_steps > CAT_WARP_DEMONS_MAX_STEPS)
        opt.n_steps = CAT_WARP_DEMONS_MAX_STEPS;

    /* Build the coarse-to-fine pyramid: finest level = n_points, each coarser
     * level uses 1/4 of the points (next tetrahedral subdivision down). */
    {
        int lvl, np = n_points;
        for (lvl = opt.n_steps - 1; lvl >= 0; lvl--)
        {
            opt.level_points[lvl] = np;
            np /= 4;
        }
    }

    warped_src_sphere = (polygons_struct *)malloc(sizeof(polygons_struct));

    if (CAT_WarpDemonsRegister(src, src_sphere, trg, trg_sphere,
                               warped_src_sphere, &opt) != OK)
    {
        fprintf(stderr, "Registration failed.\n");
        return EXIT_FAILURE;
    }

    if (output_sphere_file != NULL)
    {
        object_struct *out = create_object(POLYGONS);
        *get_polygons_ptr(out) = *warped_src_sphere;
        if (output_graphics_any_format(output_sphere_file, format, 1, &out, NULL) != OK)
            return EXIT_FAILURE;
    }

    if (output_surface_file != NULL)
    {
        /* The deformed sphere defines the new spherical parameterization; write
         * the input surface geometry alongside it for downstream resampling. */
        if (output_graphics_any_format(output_surface_file, format, 1, src_objects, NULL) != OK)
            return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
