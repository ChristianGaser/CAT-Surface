/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 *
 */

#include <bicpl.h>
#include <float.h>
#include <ParseArgv.h>

#include "CAT_Warp.h"
#include "CAT_Map.h"
#include "CAT_Surf.h"
#include "CAT_Curvature.h"
#include "CAT_SurfaceIO.h"
#include "CAT_Smooth.h"
#include "CAT_Resample.h"
#include "dartel.h"
#include "CAT_SafeAlloc.h"

/* argument defaults */
char *param_file     = NULL;
char *source_file    = NULL;
char *source_sphere_file = NULL;
char *target_file    = NULL;
char *target_sphere_file = NULL;
char *jacdet_file    = NULL;
char *output_sphere_file = NULL;
char *pgm_file       = NULL;

int rotate     = 1;
int avg      = 0;
int code     = 1;
int loop     = 6;
int verbose    = 0;
int rtype    = 1;
int curvtype0  = 5;
int curvtype1  = 5;
int curvtype2  = 2;
int muchange   = 4;
int distortion_correction = 0;
int multires_levels = 0;  /* 0 = disabled, 1-3 = number of coarse levels */
int sz_map[2]  = {512, 256};
int n_triangles  = 81920;
int n_steps    = 2;
int n_runs     = 2;
int debug    = 0;
int debug_dc_weights = 0;
double murate  = 1.25;
double lambda  = 0;
double mu    = 0.125;
double lmreg   = 1e-3;
double fwhm    = 5.0;
double fwhm_surf = 0.0;

static ArgvInfo argTable[] = {
  {"-i", ARGV_STRING, (char *) 1, (char *) &source_file, 
   "Input surface file."},
  {"-is", ARGV_STRING, (char *) 1, (char *) &source_sphere_file, 
   "Input sphere file."},
  {"-t", ARGV_STRING, (char *) 1, (char *) &target_file, 
   "Template surface file."},
  {"-ts", ARGV_STRING, (char *) 1, (char *) &target_sphere_file, 
   "Template sphere file."},
  {"-ws", ARGV_STRING, (char *) 1, (char *) &output_sphere_file, 
   "Warped input sphere."},
  {"-o", ARGV_STRING, (char *) 1, (char *) &pgm_file, 
   "Warped map as pgm file."},
  {"-j", ARGV_STRING, (char *) 1, (char *) &jacdet_file, 
   "Save Jacobian determinant values (subtract 1 to ease the use of relative volume changes) of the surface."},
  {"-p", ARGV_STRING, (char *) 1, (char *) &param_file, 
   "Parameter file."},
  {"-code", ARGV_INT, (char *) 1, (char *) &code,
   "Objective function (code): 0 - sum of squares; 1 - symmetric sum of squares; 2 - multinomial."},
  {"-rtype", ARGV_INT, (char *) 1, (char *) &rtype,
   "Regularization type: 0 - linear elastic energy; 1 - membrane energy; 2 - bending energy."},
  {"-mu", ARGV_FLOAT, (char *) 1, (char *) &mu,
   "Regularization parameter mu."},
  {"-muchange", ARGV_INT, (char *) TRUE, (char *) &muchange,
   "Decrease mu after muchange loops."},
  {"-murate", ARGV_FLOAT, (char *) TRUE, (char *) &murate,
   "Divide mu after muchange loops with murate."},
  {"-lambda", ARGV_FLOAT, (char *) 1, (char *) &lambda,
   "Regularization parameter lambda."},
  {"-lmreg", ARGV_FLOAT, (char *) 1, (char *) &lmreg,
   "LM regularization."},
  {"-fwhm", ARGV_FLOAT, (char *) 1, (char *) &fwhm,
   "Filter size for curvature map in FWHM. This filter size is decreased by factor 3 with each step."},
  {"-fwhm-surf", ARGV_FLOAT, (char *) 1, (char *) &fwhm_surf,
   "Filter size for smoothing surface in FWHM. This filter size is decreased by factor 3 with each step."},
  {"-loop", ARGV_INT, (char *) 1, (char *) &loop,
   "Number of outer Dartel loops for default parameters (max. 6)."},
  {"-steps", ARGV_INT, (char *) 1, (char *) &n_steps,
   "Number of Dartel steps (max. 3):\n\t1 - Inflated surface\n\t2 - High smoothed surface\n\t3 - Low smoothed surface."},
  {"-runs", ARGV_INT, (char *) 1, (char *) &n_runs,
   "Number of runs (repetitions) for whole Dartel approach."},
  {"-size", ARGV_INT, (char *) 2, (char *) sz_map,
   "Size of curvature map for warping."},
  {"-norot", ARGV_CONSTANT, (char *) FALSE, (char *) &rotate,
   "Don't rotate input surface before warping."},
  {"-avg", ARGV_CONSTANT, (char *) TRUE, (char *) &avg,
   "Average together two weighted DARTEL solutions into final mesh."},
  {"-type0", ARGV_INT, (char *) 1, (char *) &curvtype0,
   "Curvature type for 1st step\n\t0 - mean curvature (averaged over 3mm, in degrees)\n\t1 - gaussian curvature\n\t2 - curvedness\n\t3 - shape index\n\t4 - mean curvature (in radians)\n\t5 - sulcal depth like estimator\n\t>5 - depth potential with parameter alpha = 1/curvtype."},
  {"-type1", ARGV_INT, (char *) 1, (char *) &curvtype1,
   "Curvature type for the 2nd step\n\t0 - mean curvature (averaged over 3mm, in degrees)\n\t1 - gaussian curvature\n\t2 - curvedness\n\t3 - shape index\n\t4 - mean curvature (in radians)\n\t5 - sulcal depth like estimator\n\t>5 - depth potential with parameter alpha = 1/curvtype."},
  {"-type2", ARGV_INT, (char *) 1, (char *) &curvtype2,
   "Curvature type for the 3rd step\n\t0 - mean curvature (averaged over 3mm, in degrees)\n\t1 - gaussian curvature\n\t2 - curvedness\n\t3 - shape index\n\t4 - mean curvature (in radians)\n\t5 - sulcal depth like estimator\n\t>5 - depth potential with parameter alpha = 1/curvtype."},
  {"-distortion-correction", ARGV_CONSTANT, (char *) TRUE, (char *) &distortion_correction,
   "Apply distortion correction weighting (sin(theta)) to 2D registration to reduce bias at poles. References: Fischl et al. 2008 (https://pmc.ncbi.nlm.nih.gov/articles/PMC7784120/), Schuh et al. 2024 (https://www.sciencedirect.com/science/article/pii/S1361841524002172)."},
  {"-multires", ARGV_INT, (char *) 1, (char *) &multires_levels,
   "Multi-resolution levels (0-3). Start registration at coarser resolution and progressively refine.\n\t0 - disabled (default)\n\t1 - 2 levels (256x128 -> 512x256)\n\t2 - 3 levels (128x64 -> 256x128 -> 512x256)\n\t3 - 4 levels (64x32 -> 128x64 -> 256x128 -> 512x256)."},
  {"-verbose", ARGV_CONSTANT, (char *) TRUE, (char *) &verbose,
   "Be verbose."},
  {"-debug", ARGV_CONSTANT, (char *) TRUE, (char *) &debug,
   "Save debug files."},
    {"-debug-dc-weights", ARGV_CONSTANT, (char *) TRUE, (char *) &debug_dc_weights,
     "Save dc_weights.pgm (requires -distortion-correction)."},
   {NULL, ARGV_END, NULL, NULL, NULL}
};

void
solve_dartel_flow(polygons_struct *src, polygons_struct *src_sphere,
          polygons_struct *trg, polygons_struct *trg_sphere,
          struct dartel_prm *prm, int dm[3], int n_steps,
          double rot[3], double *flow, int n_loops, int distortion_correction)
{
    int        step, i, x, y, it, it0, it1, xy_size, it_scratch, curvtype;
    int        res_level, n_res_levels, loops_per_level;
    int        cur_dm[3], prev_dm[3];
    polygons_struct  *sm_src, *sm_trg, *sm_src_sphere, *sm_trg_sphere;
    double       rotation_matrix[9];
    double       *flow1, *inflow, *map_src, *map_trg;
    double       *scratch, *jd, *jd1, *values;
    double       *dc_weights, *prev_flow, *cur_flow;
    double       *map_src_full, *map_trg_full;
    double       ll[3];
    double       cur_fwhm, cur_fwhm_surf;

    xy_size = dm[0] * dm[1];

    /* Determine number of resolution levels */
    n_res_levels = multires_levels + 1;
    if (n_res_levels > 4) n_res_levels = 4;
    if (n_res_levels < 1) n_res_levels = 1;
    
    /* Distribute loops across resolution levels */
    loops_per_level = (n_loops + n_res_levels - 1) / n_res_levels;

    flow1 = (double *) malloc(sizeof(double) * xy_size * 2);
    inflow = (double *) malloc(sizeof(double) * xy_size * 2);
    map_src_full = (double *) malloc(sizeof(double) * xy_size);
    map_trg_full = (double *) malloc(sizeof(double) * xy_size);
    sm_src = (polygons_struct *) malloc(sizeof(polygons_struct));
    sm_src_sphere = (polygons_struct *) malloc(sizeof(polygons_struct));
    sm_trg = (polygons_struct *) malloc(sizeof(polygons_struct));
    sm_trg_sphere = (polygons_struct *) malloc(sizeof(polygons_struct));

    /* Initialize flow to zero */
    for (i = 0; i < xy_size * 2; i++) inflow[i] = 0.0;
    
    for (step = 0; step < n_steps; step++) {
        cur_fwhm = fwhm;
        cur_fwhm_surf = fwhm_surf;
        
        /* resample source and target surface */
        resample_spherical_surface(src, src_sphere, sm_src, NULL, NULL,
                       n_triangles);
        resample_spherical_surface(trg, trg_sphere, sm_trg, NULL, NULL,
                       n_triangles);

        /* initialization */
        if (step == 0) {
             curvtype = curvtype0;
             if (cur_fwhm_surf > 0) {
                 smooth_heatkernel(sm_src, NULL, cur_fwhm_surf);
                 smooth_heatkernel(sm_trg, NULL, cur_fwhm_surf);
             }
             resample_spherical_surface(src_sphere, src_sphere,
                           sm_src_sphere, NULL, NULL,
                           n_triangles);
             resample_spherical_surface(trg_sphere, trg_sphere,
                           sm_trg_sphere, NULL, NULL,
                           n_triangles);

            /* initial rotation if n_loops < 0 */
            if (n_loops < 0) {
                rotate_polygons_to_atlas(sm_src, sm_src_sphere,
                             sm_trg, sm_trg_sphere,
                             cur_fwhm, curvtype0, rot, verbose, distortion_correction);
                rotation_to_matrix(rotation_matrix,
                           rot[0], rot[1], rot[2]);
        
                /* rotate source sphere */
                rotate_polygons(src_sphere,
                        NULL, rotation_matrix);

                resample_spherical_surface(src_sphere,
                               src_sphere,
                               sm_src_sphere, NULL,
                               NULL, n_triangles);
            }

            for (i = 0; i < xy_size*2; i++)  inflow[i] = 0.0;

        } else if (step == 1) {
             curvtype = curvtype1;
             if (cur_fwhm_surf > 0) {
                 smooth_heatkernel(sm_src, NULL, cur_fwhm_surf);
                 smooth_heatkernel(sm_trg, NULL, cur_fwhm_surf);
             }
             
        } else if (step == 2) {
             curvtype = curvtype2;
            
             if (cur_fwhm_surf > 0) {
                 smooth_heatkernel(sm_src, NULL, cur_fwhm_surf);
                 smooth_heatkernel(sm_trg, NULL, cur_fwhm_surf);
             }
        }

        /* get curvatures at full resolution */
        map_sphere_values_to_sheet(sm_trg, sm_trg_sphere, (double *)0,
                         map_trg_full, cur_fwhm, dm, curvtype);
        map_sphere_values_to_sheet(sm_src, sm_src_sphere, (double *)0, 
                         map_src_full, cur_fwhm, dm, curvtype);
        
        if (debug) {
            if (write_pgm("source.pgm", map_src_full, dm[0], dm[1]) != 0)
                exit(EXIT_FAILURE);
            if (write_pgm("target.pgm", map_trg_full, dm[0], dm[1]) != 0)
                exit(EXIT_FAILURE);
        }

        /* Multi-resolution loop: from coarse to fine */
        for (res_level = n_res_levels - 1; res_level >= 0; res_level--) {
            int cur_xy_size, level_loops;
            
            /* Calculate current resolution dimensions */
            cur_dm[0] = dm[0] >> res_level;  /* Divide by 2^res_level */
            cur_dm[1] = dm[1] >> res_level;
            cur_dm[2] = 1;
            cur_xy_size = cur_dm[0] * cur_dm[1];
            
            /* Allocate arrays for current resolution */
            map_src = (double *) malloc(sizeof(double) * cur_xy_size);
            map_trg = (double *) malloc(sizeof(double) * cur_xy_size);
            cur_flow = (double *) malloc(sizeof(double) * cur_xy_size * 2);
            
            /* Distortion correction weights for current resolution */
            dc_weights = (double *) 0;
            if (distortion_correction) {
                dc_weights = (double *) malloc(sizeof(double) * cur_xy_size);
                if (!dc_weights) {
                    fprintf(stderr, "Memory allocation error for distortion-correction weights.\n");
                    exit(EXIT_FAILURE);
                }
                for (y = 0; y < cur_dm[1]; y++) {
                    double v = ((double) y + 0.5) / (double) cur_dm[1];
                    double theta = v * PI;
                    double w = sin(theta);
                    if (w < 1e-10) w = 1e-10;
                    for (x = 0; x < cur_dm[0]; x++) {
                        dc_weights[x + cur_dm[0] * y] = w;
                    }
                }
            }
            
            /* Downsample curvature maps if at coarser level */
            if (res_level > 0) {
                int tmp_dm[3];
                double *tmp_src, *tmp_trg;
                int level;
                
                /* Start with full resolution */
                tmp_dm[0] = dm[0];
                tmp_dm[1] = dm[1];
                tmp_dm[2] = 1;
                tmp_src = map_src_full;
                tmp_trg = map_trg_full;
                
                /* Progressively downsample to current resolution */
                for (level = 0; level < res_level; level++) {
                    int next_dm[3];
                    double *next_src, *next_trg;
                    
                    next_dm[0] = tmp_dm[0] / 2;
                    next_dm[1] = tmp_dm[1] / 2;
                    next_dm[2] = 1;
                    
                    if (level == res_level - 1) {
                        /* Last step: downsample directly to map_src/map_trg */
                        downsample_image(tmp_src, tmp_dm, map_src, next_dm);
                        downsample_image(tmp_trg, tmp_dm, map_trg, next_dm);
                    } else {
                        /* Intermediate step: allocate temp arrays */
                        next_src = (double *) malloc(sizeof(double) * next_dm[0] * next_dm[1]);
                        next_trg = (double *) malloc(sizeof(double) * next_dm[0] * next_dm[1]);
                        downsample_image(tmp_src, tmp_dm, next_src, next_dm);
                        downsample_image(tmp_trg, tmp_dm, next_trg, next_dm);
                        
                        if (tmp_src != map_src_full) free(tmp_src);
                        if (tmp_trg != map_trg_full) free(tmp_trg);
                        
                        tmp_src = next_src;
                        tmp_trg = next_trg;
                    }
                    
                    tmp_dm[0] = next_dm[0];
                    tmp_dm[1] = next_dm[1];
                }
                
                /* Free intermediate buffers if any remain */
                if (tmp_src != map_src_full && tmp_src != map_src) free(tmp_src);
                if (tmp_trg != map_trg_full && tmp_trg != map_trg) free(tmp_trg);
            } else {
                /* At full resolution, just copy */
                for (i = 0; i < cur_xy_size; i++) {
                    map_src[i] = map_src_full[i];
                    map_trg[i] = map_trg_full[i];
                }
            }
            
            /* Handle flow field initialization/upsampling */
            if (res_level == n_res_levels - 1) {
                /* At coarsest level: reallocate inflow to current size and initialize to zero */
                free(inflow);
                inflow = (double *) malloc(sizeof(double) * cur_xy_size * 2);
                for (i = 0; i < cur_xy_size * 2; i++) {
                    cur_flow[i] = 0.0;
                    inflow[i] = 0.0;
                }
            } else {
                /* Upsample flow from previous coarser level */
                double *new_inflow = (double *) malloc(sizeof(double) * cur_xy_size * 2);
                upsample_flow_field(inflow, prev_dm, new_inflow, cur_dm);
                free(inflow);
                inflow = new_inflow;
                /* Copy to cur_flow as starting point */
                for (i = 0; i < cur_xy_size * 2; i++) {
                    cur_flow[i] = inflow[i];
                }
            }
            
            /* Determine loops for this level */
            level_loops = loops_per_level;
            if (res_level == 0) {
                /* Final level gets remaining loops */
                level_loops = n_loops - loops_per_level * (n_res_levels - 1);
                if (level_loops < 1) level_loops = 1;
            }
            
            if (verbose && n_res_levels > 1) {
                fprintf(stdout, "Resolution level %d: %dx%d (%d loops)\n", 
                        n_res_levels - res_level, cur_dm[0], cur_dm[1], level_loops);
            }

            /* go through dartel steps at current resolution */
            for (it = 0, it0 = 0; it0 < level_loops; it0++) {
                it_scratch = dartel_scratchsize((int *)cur_dm,
                                 prm[it0].code);

                scratch = (double *) malloc(sizeof(double)*it_scratch);
                for (it1 = 0; it1 < prm[it0].its; it1++) {
                    it++;
                    /* map target onto source */
                    if (INVERSE_WARPING) {
                        dartel(prm[it0], cur_dm, inflow, map_src,
                            map_trg, dc_weights, cur_flow, ll, scratch);
                    } else {
                        dartel(prm[it0], cur_dm, inflow, map_trg,
                            map_src, dc_weights, cur_flow, ll, scratch);
                    }
                    if (verbose) {
                        if (n_res_levels > 1)
                            fprintf(stdout,"\r  %02d-%02d-%02d: %8.2f", step+1, n_res_levels - res_level, it, ll[0]);
                        else
                            fprintf(stdout,"\r%02d-%02d: %8.2f", step+1, it, ll[0]);
                        fflush(stdout);
                    }

                    for (i = 0; i < cur_xy_size*2; i++)
                        inflow[i] = cur_flow[i];
                }
                free(scratch);
            }
            
            /* Save current dimensions for next level's upsampling */
            prev_dm[0] = cur_dm[0];
            prev_dm[1] = cur_dm[1];
            prev_dm[2] = cur_dm[2];
            
            /* Free current level's temporary arrays */
            free(map_src);
            free(map_trg);
            if (dc_weights) free(dc_weights);
            
            /* At final level, copy result to output flow */
            if (res_level == 0) {
                for (i = 0; i < xy_size * 2; i++) {
                    flow[i] = cur_flow[i];
                }
            }
            free(cur_flow);
        }

        /* use smaller FWHM for next steps */
        fwhm /= 3.0;
        fwhm_surf /= 3.0;
    }
    if (verbose) fprintf(stdout,"\n");

    free(sm_src);
    free(sm_trg);
    free(sm_src_sphere);
    free(sm_trg_sphere);
    free(map_src_full);
    free(map_trg_full);
    
    /* get deformations and jacobian det. from flow field */
    if (jacdet_file != NULL) {
        fprintf(stdout,"Warning: Saving jacobians not working\n");
        if (rotate)
             fprintf(stdout,"Warning: Rotation not yet considered\n");
        jd  = (double *) malloc(sizeof(double) * xy_size);
        jd1 = (double *) malloc(sizeof(double) * xy_size);

        /* Use the flow field stored during the final level */
        for (i = 0; i < xy_size * 2; i++) inflow[i] = flow[i];
        
        expdefdet(dm, 10, inflow, flow, flow1, jd, jd1);
        
        /* subtract 1 to get values around 0 instead of 1 and invert */
        for (i = 0; i < xy_size; i++) {
            jd1[i] -= 1;
            jd1[i] *= -1;
        }

        values = (double *) malloc(sizeof(double) *
                       src_sphere->n_points);

        map_sheet2d_to_sphere(jd1, values, src_sphere, 1, dm);

        output_values_any_format(jacdet_file, src_sphere->n_points,
                     values, TYPE_DOUBLE);

        free(values);
        free(jd1);
        free(jd);
        free(inflow);
    } else {
        for (i = 0; i < xy_size * 2; i++) inflow[i] = flow[i];
        expdef(dm, 10, inflow, flow, flow1,
             (double *) 0, (double *) 0);
        free(inflow);
    }

    free(flow1);
}

int
main(int argc, char *argv[])
{
    File_formats   format;
    FILE       *fp;
    char       line[1024];
    polygons_struct  *src, *trg, *src_sphere, *trg_sphere;
    polygons_struct  *rsrc, *rs_sph, *rtrg, *rt_sph;
    polygons_struct  *as_sph;
    int        i, j, run;
    int        n_objects, xy_size;
    double       *flow, *flow2, *data;
    object_struct  **objects;
    static double  param[2] = {1.0, 1.0};
    int        dm[3];
    double       rotation_matrix[9];
    struct       dartel_prm* prm;
    double       rot[3];


    /* get the arguments from the command line */
    if (ParseArgv(&argc, argv, argTable, 0) ||
      source_file == NULL || target_file == NULL || 
      source_sphere_file == NULL || target_sphere_file == NULL ||
      (jacdet_file == NULL && pgm_file == NULL && output_sphere_file == NULL)) {
        fprintf(stderr, "\nUsage: %s [options]\n", argv[0]);
        fprintf(stderr, "   %s -help\n\n", argv[0]);
        exit(EXIT_FAILURE);
    }

    /* size of curvature map for warping */
    dm[0] = sz_map[0];
    dm[1] = sz_map[1];
    dm[2] = 1;

    if (input_graphics_any_format(target_file, &format,
                    &n_objects, &objects) != OK)
        exit(EXIT_FAILURE);

    /* get a pointer to the surface */
    trg = get_polygons_ptr(objects[0]);

    /* check that the surface file contains a polyhedron */
    if (n_objects != 1 || get_object_type(objects[0]) != POLYGONS) {
        fprintf(stdout,"Surface file must contain 1 polygons object.\n");
        exit(EXIT_FAILURE);
    }

    if (input_graphics_any_format(source_file, &format,
                    &n_objects, &objects) != OK)
        exit(EXIT_FAILURE);

    /* get a pointer to the surface */
    src = get_polygons_ptr(objects[0]);

    /* check that the surface file contains a polyhedron */
    if (n_objects != 1 || get_object_type(objects[0]) != POLYGONS) {
        fprintf(stdout,"Surface file must contain 1 polygons object.\n");
        exit(EXIT_FAILURE);
    }
  
    /* read sphere for input surface */
    if (input_graphics_any_format(source_sphere_file, &format,
                &n_objects, &objects) != OK)
        exit(EXIT_FAILURE);
    src_sphere = get_polygons_ptr(objects[0]);

    /* read sphere for template surface */
    if (input_graphics_any_format(target_sphere_file, &format,
                &n_objects, &objects) != OK)
        exit(EXIT_FAILURE);
    trg_sphere = get_polygons_ptr(objects[0]);

    translate_to_center_of_mass(src_sphere);
    for (i = 0; i < src_sphere->n_points; i++)
        set_vector_length(&src_sphere->points[i], 1.0);
    translate_to_center_of_mass(trg_sphere);
    for (i = 0; i < trg_sphere->n_points; i++)
        set_vector_length(&trg_sphere->points[i], 1.0);

    prm = (struct dartel_prm*) malloc(sizeof(struct dartel_prm) * 100);
    
    /* first two entries of param are equal */
    for (j = 0; j < loop; j++) {
        for (i = 0; i < 2; i++)
            prm[j].rparam[i] = param[i];
    }

    /* read values from parameter file */
    if (param_file != NULL) {
        fp = SAFE_FOPEN(param_file, "r");
  
        j = 0;

        fprintf(stdout,"Read parameters from %s\n", param_file);
        while (fgets(line, sizeof(line), fp)) {
            /* check for 9 values in each line */
            if (sscanf(line, "%d %lf %lf %lf %lf %d %d %d %d",
                   &prm[j].rtype, &prm[j].rparam[2],
                   &prm[j].rparam[3], &prm[j].rparam[4],
                   &prm[j].lmreg, &prm[j].cycles,
                   &prm[j].its, &prm[j].k,
                   &prm[j].code) != 9)
                continue;
            j++;
        }
        fclose(fp);

        if (j == 0) {
            fprintf(stderr, "Could not read parameter file %s. Check that each line contains 9 values\n", param_file);
            exit(EXIT_FAILURE);
        }

    } else { /* use predefined values */
        for (j = 0; j < loop; j++) {
            /* some entries are equal */
            prm[j].rtype = rtype;
            prm[j].cycles = 3;
            prm[j].its = 3;
            prm[j].code = code;
            prm[j].lmreg = lmreg;
            prm[j].rparam[2] = mu;
            prm[j].rparam[3] = lambda;
            prm[j].rparam[4] = lambda/2.0;
            prm[j].k = j;
            if ((j+1) % muchange == 0) mu /= murate;
            lambda /= 5.0;
        }
    }

    if (verbose) {
        fprintf(stdout,"___________________________________");
        fprintf(stdout,"________________________________________\n");
        fprintf(stdout,"Parameters\n");
        fprintf(stdout,"___________________________________");
        fprintf(stdout,"________________________________________\n");
        fprintf(stdout,"Regularization (0 - elastic; 1 - membrane; ");
        fprintf(stdout,"2 - bending):\t\t%d\n", prm[0].rtype);
        fprintf(stdout,"Number of cycles for full multi grid (FMG):");
        fprintf(stdout,"\t\t\t\t%d\n", prm[0].cycles);
        fprintf(stdout,"Number of relaxation iterations in each ");
        fprintf(stdout,"Multigrid cycle:\t\t%d\n", prm[0].its);
        fprintf(stdout,"Objective function (0 - sum of squares; ");
        fprintf(stdout,"1 - sym. sum of squares):\t%d\n", prm[0].code);
        fprintf(stdout,"Levenberg-Marquardt regularization:");
        fprintf(stdout,"\t\t\t\t\t%g\n", prm[0].lmreg);
        fprintf(stdout,"Curvature types:\t\t\t\t\t\t\t%d", curvtype0);
        if (n_steps > 1) fprintf(stdout,"/%d",curvtype1);
        if (n_steps > 2) fprintf(stdout,"/%d",curvtype2);
        fprintf(stdout,"\n%d Iterative loops\n", loop);
        fprintf(stdout,"\nRegularization parameter mu:\t\t");
        for (i = 0; i < loop; i++)
            fprintf(stdout,"%8g\t", prm[i].rparam[2]);
        fprintf(stdout,"\n");
        fprintf(stdout,"Regularization parameter lambda:\t");
        for (i = 0; i < loop; i++)
            fprintf(stdout,"%8g\t", prm[i].rparam[3]);
        fprintf(stdout,"\n");
        fprintf(stdout,"Regularization parameter id:\t\t");
        for (i = 0; i < loop; i++)
            fprintf(stdout,"%8g\t", prm[i].rparam[4]);
        fprintf(stdout,"\n");
        fprintf(stdout,"Time steps for solving the PDE:\t\t");
        for (i = 0; i < loop; i++)
            fprintf(stdout,"%8d\t",prm[i].k);
        fprintf(stdout,"\n\n");
    }
  
    xy_size = dm[0] * dm[1];
    flow  = (double *) malloc(sizeof(double) * xy_size * 2);
    
    /* estimate rotation only */
    if (rotate) {
        solve_dartel_flow(src, src_sphere, trg, trg_sphere, prm, dm,
                  n_steps, rot, flow, -1, distortion_correction);
    }
    
    if (debug && rotate) {
        data = (double *) malloc(sizeof(double) * dm[0] * dm[1]);

        map_sphere_values_to_sheet(src, src_sphere, (double *)0,
                         data, 0.0, dm, curvtype0);

        if (write_pgm("source_rotated.pgm", data, dm[0], dm[1]) != 0)
            exit(EXIT_FAILURE);

        free(data);
    }

    if (avg) {
        flow2 = (double *) malloc(sizeof(double) * xy_size * 2);
        rsrc = (polygons_struct *) malloc(sizeof(polygons_struct));
        rs_sph  = (polygons_struct *) malloc(sizeof(polygons_struct));
        rtrg = (polygons_struct *) malloc(sizeof(polygons_struct));
        rt_sph  = (polygons_struct *) malloc(sizeof(polygons_struct));
        as_sph  = (polygons_struct *) malloc(sizeof(polygons_struct));

        rotation_to_matrix(rotation_matrix, 0.0, PI/2.0, 0.0);
        rotate_polygons(trg, rtrg, rotation_matrix);
        rotate_polygons(trg_sphere, rt_sph, rotation_matrix);
    }
    
    /* run dartel */
    for (run = 0; run < n_runs; run++) {
        solve_dartel_flow(src, src_sphere, trg, trg_sphere, prm, dm, n_steps,
              rot, flow, loop, distortion_correction);

        /* solve again, but rotated to change pole location */
        if (avg && (run==(n_runs-1))) {
            rotation_to_matrix(rotation_matrix, 0.0, PI/2.0, 0.0);
            rotate_polygons(src, rsrc, rotation_matrix);
            rotate_polygons(src_sphere, rs_sph, rotation_matrix);

            solve_dartel_flow(rsrc, rs_sph, rtrg, rt_sph, prm,
                  dm, n_steps, rot, flow2, loop, distortion_correction);

            apply_warp(src_sphere, src_sphere, flow, dm, !INVERSE_WARPING);
            apply_warp(rs_sph, rs_sph, flow2, dm, !INVERSE_WARPING); 
            
            rotation_to_matrix(rotation_matrix, 0.0, -PI/2.0, 0.0);
            rotate_polygons(rs_sph, NULL, rotation_matrix);

            average_xz_surf(rs_sph, src_sphere, as_sph);
            copy_polygons(as_sph, src_sphere);
        } else {
            apply_warp(src_sphere, src_sphere, flow, dm,
                   !INVERSE_WARPING);
        }
    }

    if (avg) {
        free(flow2);
        free(rsrc);
        free(rs_sph);
        free(rtrg);
        free(rt_sph);
        free(as_sph);
    }

    if (output_sphere_file != NULL) {
        compute_polygon_normals(src_sphere);
        *get_polygons_ptr(objects[0]) = *src_sphere;
        if (output_graphics_any_format(output_sphere_file, format,
                         n_objects, objects, NULL) != OK)
            exit(EXIT_FAILURE);
    }

    if (pgm_file != NULL) {
        data = (double *) malloc(sizeof(double) * dm[0] * dm[1]);

        map_sphere_values_to_sheet(src, src_sphere, (double *)0,
                         data, 0.0, dm, curvtype0);

        if (write_pgm(pgm_file, data, dm[0], dm[1]) != 0)
            exit(EXIT_FAILURE);

        free(data);
    }


    delete_object_list(n_objects, objects);
    free(flow);
    free(prm);

    return(EXIT_SUCCESS);
}
