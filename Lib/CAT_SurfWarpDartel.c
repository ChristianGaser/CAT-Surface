/*
 * CAT_SurfWarpDartel.c
 */

#include <bicpl.h>
#include <float.h>
#include <stdio.h>
#include <stdlib.h>

#include "CAT_SurfWarpDartel.h"

#include "CAT_Warp.h"
#include "CAT_Map.h"
#include "CAT_Smooth.h"
#include "CAT_Resample.h"
#include "CAT_SafeAlloc.h"

Status CAT_SurfWarpSolveDartelFlow(
    polygons_struct *src,
    polygons_struct *src_sphere,
    polygons_struct *trg,
    polygons_struct *trg_sphere,
    struct dartel_prm *prm,
    int dm[3],
    int n_steps,
    double rot[3],
    double *flow,
    int n_loops,
    const CAT_SurfWarpDartelOptions *opt)
{
    int step, i, it, it0, it1, xy_size, it_scratch, curvtype;
    int res_level, n_res_levels, loops_per_level;
    int cur_dm[3], prev_dm[3];
    polygons_struct *sm_src, *sm_trg, *sm_src_sphere, *sm_trg_sphere;
    double rotation_matrix[9];
    double *flow1, *inflow, *map_src, *map_trg;
    double *scratch, *jd, *jd1, *values;
    double *dc_weights, *cur_flow;
    double *map_src_full, *map_trg_full;
    double ll[3];
    double cur_fwhm, cur_fwhm_surf;

    if (!src || !src_sphere || !trg || !trg_sphere || !prm || !dm || !rot || !flow || !opt) {
        fprintf(stderr, "CAT_SurfWarpSolveDartelFlow: invalid arguments.\n");
        return ERROR;
    }
    if (!opt->fwhm || !opt->fwhm_surf) {
        fprintf(stderr, "CAT_SurfWarpSolveDartelFlow: fwhm pointers must not be NULL.\n");
        return ERROR;
    }

    xy_size = dm[0] * dm[1];

    /* Determine number of resolution levels */
    n_res_levels = opt->multires_levels + 1;
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

    if (!flow1 || !inflow || !map_src_full || !map_trg_full || !sm_src || !sm_src_sphere || !sm_trg || !sm_trg_sphere) {
        fprintf(stderr, "CAT_SurfWarpSolveDartelFlow: out of memory.\n");
        return ERROR;
    }

    /* Initialize flow to zero */
    for (i = 0; i < xy_size * 2; i++) inflow[i] = 0.0;

    for (step = 0; step < n_steps; step++) {
        cur_fwhm = *opt->fwhm;
        cur_fwhm_surf = *opt->fwhm_surf;

        /* resample source and target surface */
        resample_spherical_surface(src, src_sphere, sm_src, NULL, NULL,
                                   opt->n_triangles);
        resample_spherical_surface(trg, trg_sphere, sm_trg, NULL, NULL,
                                   opt->n_triangles);

        /* initialization */
        if (step == 0) {
            curvtype = opt->curvtype0;
            if (cur_fwhm_surf > 0) {
                smooth_heatkernel(sm_src, NULL, cur_fwhm_surf);
                smooth_heatkernel(sm_trg, NULL, cur_fwhm_surf);
            }
            resample_spherical_surface(src_sphere, src_sphere,
                                       sm_src_sphere, NULL, NULL,
                                       opt->n_triangles);
            resample_spherical_surface(trg_sphere, trg_sphere,
                                       sm_trg_sphere, NULL, NULL,
                                       opt->n_triangles);

            /* initial rotation if n_loops < 0 */
            if (n_loops < 0) {
                rotate_polygons_to_atlas(sm_src, sm_src_sphere,
                                         sm_trg, sm_trg_sphere,
                                         cur_fwhm, opt->curvtype0, rot, opt->verbose);
                rotation_to_matrix(rotation_matrix,
                                   rot[0], rot[1], rot[2]);

                /* rotate source sphere */
                rotate_polygons(src_sphere,
                                NULL, rotation_matrix);

                resample_spherical_surface(src_sphere,
                                           src_sphere,
                                           sm_src_sphere, NULL,
                                           NULL, opt->n_triangles);
            }

            for (i = 0; i < xy_size * 2; i++) inflow[i] = 0.0;

        } else if (step == 1) {
            curvtype = opt->curvtype1;
            if (cur_fwhm_surf > 0) {
                smooth_heatkernel(sm_src, NULL, cur_fwhm_surf);
                smooth_heatkernel(sm_trg, NULL, cur_fwhm_surf);
            }

        } else {
            curvtype = opt->curvtype2;

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

        if (opt->debug) {
            if (write_pgm("source.pgm", map_src_full, dm[0], dm[1]) != 0)
                exit(EXIT_FAILURE);
            if (write_pgm("target.pgm", map_trg_full, dm[0], dm[1]) != 0)
                exit(EXIT_FAILURE);
        }

        /* Multi-resolution loop: from coarse to fine */
        for (res_level = n_res_levels - 1; res_level >= 0; res_level--) {
            int cur_xy_size, level_loops;

            /* Calculate current resolution dimensions */
            cur_dm[0] = dm[0] >> res_level; /* Divide by 2^res_level */
            cur_dm[1] = dm[1] >> res_level;
            cur_dm[2] = 1;
            cur_xy_size = cur_dm[0] * cur_dm[1];

            /* Allocate arrays for current resolution */
            map_src = (double *) malloc(sizeof(double) * cur_xy_size);
            map_trg = (double *) malloc(sizeof(double) * cur_xy_size);
            cur_flow = (double *) malloc(sizeof(double) * cur_xy_size * 2);

            /* Distortion correction weights for current resolution */
            dc_weights = (double *)0;

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

            if (opt->verbose && n_res_levels > 1) {
                fprintf(stdout, "  Resolution level %d: %dx%d (%d loops)\n",
                        n_res_levels - res_level, cur_dm[0], cur_dm[1], level_loops);
            }

            /* go through dartel steps at current resolution */
            for (it = 0, it0 = 0; it0 < level_loops; it0++) {
                it_scratch = dartel_scratchsize((int *)cur_dm,
                                               prm[it0].code);

                scratch = (double *) malloc(sizeof(double) * it_scratch);
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
                    if (opt->verbose) {
                        if (n_res_levels > 1)
                            fprintf(stdout, "\r  %02d-%02d-%02d: %8.2f", step + 1, n_res_levels - res_level, it, ll[0]);
                        else
                            fprintf(stdout, "\r%02d-%02d: %8.2f", step + 1, it, ll[0]);
                        fflush(stdout);
                    }

                    for (i = 0; i < cur_xy_size * 2; i++)
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
        *opt->fwhm /= 3.0;
        *opt->fwhm_surf /= 3.0;
    }
    if (opt->verbose) fprintf(stdout, "\n");

    free(sm_src);
    free(sm_trg);
    free(sm_src_sphere);
    free(sm_trg_sphere);
    free(map_src_full);
    free(map_trg_full);

    /* get deformations and jacobian det. from flow field */
    if (opt->jacdet_file != NULL) {
        fprintf(stdout, "Warning: Saving jacobians not working\n");
        if (opt->rotate)
            fprintf(stdout, "Warning: Rotation not yet considered\n");
        jd = (double *) malloc(sizeof(double) * xy_size);
        jd1 = (double *) malloc(sizeof(double) * xy_size);

        /* Use the flow field stored during the final level */
        for (i = 0; i < xy_size * 2; i++) inflow[i] = flow[i];

        expdefdet(dm, 10, inflow, flow, flow1, jd, jd1);

        /* subtract 1 to get values around 0 instead of 1 and invert */
        for (i = 0; i < xy_size; i++) {
            jd1[i] -= 1;
            jd1[i] *= -1;
        }

        values = (double *) malloc(sizeof(double) * src_sphere->n_points);

        map_sheet2d_to_sphere(jd1, values, src_sphere, 1, dm);

        output_values_any_format(opt->jacdet_file, src_sphere->n_points,
                                 values, TYPE_DOUBLE);

        free(values);
        free(jd1);
        free(jd);
        free(inflow);
    } else {
        for (i = 0; i < xy_size * 2; i++) inflow[i] = flow[i];
        expdef(dm, 10, inflow, flow, flow1,
               (double *)0, (double *)0);
        free(inflow);
    }

    free(flow1);

    return OK;
}
