/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 *
 */

#include <volume_io/internal_volume_io.h>
#include <bicpl.h>
#include <float.h>
#include <ParseArgv.h>

#include "CAT_SheetIO.h"
#include "CAT_Map2d.h"
#include "CAT_Surf.h"
#include "CAT_Curvature.h"
#include "CAT_SurfaceIO.h"
  
struct dartel_prm {
  int rtype;         /* regularization type: 0 - linear elastic energy; */
                     /* 1 - membrane energy; 2 - bending energy */
  double rparam[5];  /* regularization parameters: ?? ?? mu lambda id */
  double lmreg;      /* LM regularization */
  int cycles;        /* # of cycles for full multi grid (FMG) */
  int its;           /* # of relaxation iterations in each multigrid cycle */
  int k;             /* time steps for solving the PDE */
  int code;          /* objective function: 0 - sum of squares; */
                     /* 1 - symmetric sum of squares */
};

/* defaults */
char *param_file         = NULL;
char *source_file        = NULL;
char *source_sphere_file = NULL;
char *target_file        = NULL;
char *target_sphere_file = NULL;
char *weight_file        = NULL;
char *jacdet_file        = NULL;
char *deform_file        = NULL;
char *output_file        = NULL;
char *pgm_file           = NULL;
char *outflow_file       = NULL;
char *inflow_file        = NULL;

int translate = 0;
int code      = 1;
int loop      = 6;
int verbose   = 0;
int rtype     = 1;
int curvtype  = 0;
int muchange  = 4;
int sz_map[2] = {512, 256};
double murate = 1.25;
double lambda = 0.0;
double mu     = 0.25;
double lmreg  = 0.0;
double fwhm   = 10.0;

static ArgvInfo argTable[] = {
  {"-i", ARGV_STRING, (char *) 1, (char *) &source_file, 
     "Input file."},
  {"-is", ARGV_STRING, (char *) 1, (char *) &source_sphere_file, 
     "Input sphere file."},
  {"-t", ARGV_STRING, (char *) 1, (char *) &target_file, 
     "Template file."},
  {"-ts", ARGV_STRING, (char *) 1, (char *) &target_sphere_file, 
     "Template sphere file."},
  {"-w", ARGV_STRING, (char *) 1, (char *) &output_file, 
     "Warped brain."},
  {"-o", ARGV_STRING, (char *) 1, (char *) &pgm_file, 
     "Warped map as pgm file."},
  {"-j", ARGV_STRING, (char *) 1, (char *) &jacdet_file, 
     "Save Jacobian determinant values (-1 to ease the use of relative volume changes) of the surface."},
  {"-s", ARGV_STRING, (char *) 1, (char *) &weight_file, 
     "Weight warping with the inverse from values in this file (e.g. std)."},
  {"-d", ARGV_STRING, (char *) 1, (char *) &deform_file, 
     "Save amplitude of deformations (displacements) of the surface."},
  {"-if", ARGV_STRING, (char *) 1, (char *) &inflow_file, 
     "Input flow field."},
  {"-of", ARGV_STRING, (char *) 1, (char *) &outflow_file, 
     "Output flow field."},
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
     "Filter size for curvature map in FWHM."},
  {"-loop", ARGV_INT, (char *) 1, (char *) &loop,
     "Number of outer loops for default parameters (max. 6)."},
  {"-size", ARGV_INT, (char *) 2, (char *) &sz_map,
     "Size of curvature map for warping."},
  {"-shift", ARGV_CONSTANT, (char *) TRUE, (char *) &translate,
     "Shift map before warping."},
  {"-type", ARGV_INT, (char *) 1, (char *) &curvtype,
     "Curvature type\n\t0 - mean curvature (averaged over 3mm, in degrees)\n\t1 - gaussian curvature\n\t2 - curvedness\n\t3 - shape index\n\t4 - mean curvature (in radians)."},
  {"-v", ARGV_CONSTANT, (char *) TRUE, (char *) &verbose,
     "Be verbose."},
   {NULL, ARGV_END, NULL, NULL, NULL}
};

void
translate_to_template(double *source, double *target, int *size_curv,
                      int *shift)
{
        int     mn = -10, mx = 10;
        int     i, sx, sy, x, y;
        double  *shifted, sdiff, diff, sdiff0 = FLT_MAX;
  
        shifted = (double *) malloc(sizeof(double) * size_curv[0] * size_curv[1]);

        /* restricted to shifts on x-scale only */
        for (sx = mn; sx < mx; sx++) {
                for (sy = 0; sy < 1; sy++) {
                        sdiff = 0.0;
                        for (x = 0; x < size_curv[0]; x++) {
                                for (y = 0; y < size_curv[1]; y++) {
                                        i = x + (size_curv[0] * y);
                                        diff = target[i] -
                                               source[bound(x + sx, y + sy,
                                                            size_curv)];
                                        sdiff += diff*diff;
                                }
                        }
                        if (sdiff < sdiff0) {
                                shift[0] = sx;
                                shift[1] = sy;
                                sdiff0 = sdiff;
                        }
                }
        }
  
        for (x = 0; x < size_curv[0]; x++) {
                for (y = 0; y < size_curv[1]; y++) {
                        i = x + (size_curv[0] * y);
                        shifted[i] = source[bound(x + shift[0], y + shift[1],
                                                  size_curv)];
                }
        }

        for (i = 0; i < size_curv[0]*size_curv[1]; i++)
                source[i] = shifted[i];

        free(shifted);
}


int
main(int argc, char *argv[])
{
        File_formats     format;
        FILE             *fp, *fp_flow;
        char             line[1024];
        polygons_struct  *source, *target, *source_sphere, *target_sphere;
        int              x, y, i, j, it, it0, it1;
        int              n_objects, it_scratch, xy_size, n_weights;
        double           *weights;
        double           *map_source, *map_target;
        double           *map_warp, *map_weights, *map_source0, *deform;
        double           *flow, *flow1, *inflow, *scratch, *jd, *jd1, *values;
        object_struct    **objects;
        double           ll[3];
        static double    param[2] = {1.0, 1.0};
        int              size_curv[3], shift[2];
        double           xp, yp, xm, ym;
        double           H00, H01, H10, H11;
        struct dartel_prm* prm;

        /* get the arguments from the command line */

        if (ParseArgv(&argc, argv, argTable, 0) ||
             source_file == NULL || target_file == NULL ||
            (jacdet_file == NULL && output_file == NULL &&
             pgm_file == NULL   && outflow_file == NULL &&
             deform_file == NULL)) {
                fprintf(stderr, "\nUsage: %s [options]\n", argv[0]);
                fprintf(stderr, "     %s -help\n\n", argv[0]);
                exit(EXIT_FAILURE);
        }

        if (input_graphics_any_format(target_file, &format,
                                      &n_objects, &objects) != OK)
                exit(EXIT_FAILURE);

        /* check that the surface file contains a polyhedron */
        if (n_objects != 1 || get_object_type(objects[0]) != POLYGONS) {
                printf("Surface file must contain 1 polygons object.\n");
                exit(EXIT_FAILURE);
        }
        /* get a pointer to the surface */
        target = get_polygons_ptr(objects[0]);

        if (input_graphics_any_format(source_file, &format,
                                      &n_objects, &objects) != OK)
                exit(EXIT_FAILURE);

        /* get a pointer to the surface */
        source = get_polygons_ptr(objects[0]);

        /* check that the surface file contains a polyhedron */
        if (n_objects != 1 || get_object_type(objects[0]) != POLYGONS) {
                printf("Surface file must contain 1 polygons object.\n");
                exit(EXIT_FAILURE);
        }
  
        /* read sphere for input surface */
        if (source_sphere_file != NULL) {
                if (input_graphics_any_format(source_sphere_file, &format,
                                      &n_objects, &objects) != OK)
                        exit(EXIT_FAILURE);
                source_sphere = get_polygons_ptr(objects[0]);
        } else 
                source_sphere = (polygons_struct *) 0;

        /* read sphere for template surface */
        if (target_sphere_file != NULL) {
                if (input_graphics_any_format(target_sphere_file, &format,
                                      &n_objects, &objects) != OK)
                        exit(EXIT_FAILURE);
                target_sphere = get_polygons_ptr(objects[0]);
        } else 
                target_sphere = (polygons_struct *) 0;

        /* read weights */
        if (weight_file != NULL) {
                if (input_values_any_format(weight_file, &n_weights,
                                         &weights) != OK)
                        exit(EXIT_FAILURE);
        }

        prm = (struct dartel_prm*) malloc(sizeof(struct dartel_prm) * 100);
        
        /* size of curvature map for warping */
        size_curv[0] = sz_map[0];
        size_curv[1] = sz_map[1];
        size_curv[2] = 1;

        /* first two entries of param are equal */
        for (j = 0; j < loop; j++) {
                for (i = 0; i < 2; i++)
                        prm[j].rparam[i] = param[i];
        }

        /* read values from parameter file */
        if (param_file != NULL) {
                if ((fp = fopen(param_file, "r")) == 0) {
                        fprintf(stderr, "Couldn't open parameter file %s.\n",
                                param_file);
                        exit(EXIT_FAILURE);
                }
    
                loop = 0;

                fprintf(stderr, "Read parameters from %s\n", param_file);
                while (fgets(line, sizeof(line), fp)) {
                        /* check for 9 values in each line */
                        if (sscanf(line, "%d %lf %lf %lf %lf %d %d %d %d",
                                   &prm[loop].rtype, &prm[loop].rparam[2],
                                   &prm[loop].rparam[3], &prm[loop].rparam[4],
                                   &prm[loop].lmreg, &prm[loop].cycles,
                                   &prm[loop].its, &prm[loop].k,
                                   &prm[loop].code) != 9)
                                continue;
                        loop++;
                }
                fclose(fp);

                if (loop == 0) {
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
                }
                for (i = 0; i < 24; i++) {
                        prm[i].rparam[2] = mu;
                        prm[i].rparam[3] = lambda;
                        prm[i].rparam[4] = lambda/2.0;
                        prm[i].k = i;
                        if ((i+1) % muchange == 0) mu /= murate;
                        lambda /= 5.0;
                }
        }

        if (verbose) {
                fprintf(stderr, "___________________________________");
                fprintf(stderr, "________________________________________\n");
                fprintf(stderr, "Parameters\n");
                fprintf(stderr, "___________________________________");
                fprintf(stderr, "________________________________________\n");
                fprintf(stderr, "Regularization (0 - elastic; 1 - membrane; ");
                fprintf(stderr, "2 - bending):\t\t%d\n", prm[0].rtype);
                fprintf(stderr, "Number of cycles for full multi grid (FMG):");
                fprintf(stderr, "\t\t\t\t%d\n", prm[0].cycles);
                fprintf(stderr, "Number of relaxation iterations in each ");
                fprintf(stderr, "multigrid cycle:\t\t%d\n", prm[0].its);
                fprintf(stderr, "Objective function (0 - sum of squares; ");
                fprintf(stderr, "1 - sym. sum of squares):\t%d\n", prm[0].code);
                fprintf(stderr, "Levenberg-Marquardt regularization:");
                fprintf(stderr, "\t\t\t\t\t%g\n", prm[0].lmreg);
                fprintf(stderr, "\n%d Iterative loops\n", loop);
                fprintf(stderr, "\nRegularization parameter mu:\t\t");
                for (i = 0; i < loop; i++)
                        fprintf(stderr, "%8g\t", prm[i].rparam[2]);
                fprintf(stderr,"\n");
                fprintf(stderr, "Regularization parameter lambda:\t");
                for (i = 0; i < loop; i++)
                        fprintf(stderr, "%8g\t", prm[i].rparam[3]);
                fprintf(stderr,"\n");
                fprintf(stderr, "Regularization parameter id:\t\t");
                for (i = 0; i < loop; i++)
                        fprintf(stderr, "%8g\t", prm[i].rparam[4]);
                fprintf(stderr,"\n");
                fprintf(stderr, "Time steps for solving the PDE:\t\t");
                for (i = 0; i < loop; i++)
                        fprintf(stderr,"%8d\t",prm[i].k);
                fprintf(stderr,"\n\n");
        }
    
        xy_size = size_curv[0] * size_curv[1];

        flow        = (double *) malloc(sizeof(double) * xy_size * 2);
        flow1       = (double *) malloc(sizeof(double) * xy_size * 2);
        inflow      = (double *) malloc(sizeof(double) * xy_size * 2);
        map_source  = (double *) malloc(sizeof(double) * xy_size);
        map_target  = (double *) malloc(sizeof(double) * xy_size);
        map_warp    = (double *) malloc(sizeof(double) * xy_size);

        map_smoothed_curvature_to_sphere(source, source_sphere, (double *)0, map_source, fwhm,
                                         size_curv, curvtype);
        map_smoothed_curvature_to_sphere(target, target_sphere, (double *)0, map_target, fwhm,
                                         size_curv, curvtype);

        if (verbose) {
                if (write_pgm("source.pgm", map_source,
                              size_curv[0], size_curv[1]) != 0)
                        exit(EXIT_FAILURE);

                if (write_pgm("target.pgm", map_target,
                              size_curv[0], size_curv[1]) != 0)
                        exit(EXIT_FAILURE);
        }
        
        /* weight maps of target and source by the inverse std */
        if (weight_file != NULL) {
                map_weights = (double *) malloc(sizeof(double) * xy_size);
                map_source0 = (double *) malloc(sizeof(double) * xy_size);
                map_smoothed_curvature_to_sphere(source, source_sphere, weights, map_weights,
                                                 fwhm, size_curv, curvtype);
                for (i = 0; i < xy_size; i++) {
                        map_source0[i] = map_source[i];
                        map_weights[i] += 1.0;
                        map_source[i] /= map_weights[i];
                        map_target[i] /= map_weights[i];
                }
                free(map_weights);
        }
  
        if (translate) {
                translate_to_template(map_source, map_target, size_curv, shift);
                fprintf(stderr, "%d %d\n", shift[0], shift[1]);
        } else {
                shift[0] = 0; shift[1] = 0;
        }
  
        for (i = 0; i < xy_size; i++)
                map_warp[i]  = 0.0;

        if (inflow_file != NULL) {  
                if ((fp_flow = fopen(inflow_file, "rb")) == NULL) {
                        fprintf(stderr, "Error: Couldn't read file %s.\n",
                                inflow_file);
                        exit(EXIT_FAILURE);
                }

                fprintf(stderr, "Shift is not considered!!!\n");
                fread(&size_curv, 2, sizeof(int), fp_flow);
                fread(&shift, 2, sizeof(int), fp_flow);
                fread(inflow, xy_size*2, sizeof(double), fp_flow);
                fclose(fp_flow);
                size_curv[2] = 1;
        } else {
                for (i = 0; i < xy_size*2; i++) {
                        inflow[i] = 0.0;
                }
        }

        for (it = 0, it0 = 0; it0 < loop; it0++) {
                it_scratch = dartel_scratchsize((int *)size_curv,
                                                prm[it0].code);
                scratch = (double *) malloc(sizeof(double) * it_scratch);

                for (it1 = 0; it1 < prm[it0].its; it1++) {
                        it++;
                        dartel(size_curv, prm[it0].k, inflow, map_target,
                               map_source, (double *)0, prm[it0].rtype, 
                               prm[it0].rparam, prm[it0].lmreg,
                               prm[it0].cycles, prm[it0].its, prm[it0].code,
                               flow, ll, scratch);
                        fprintf(stderr, "%02d:\t%8.2f\t%5.2f\t%8.2f\t%5.2f\n",
                                it, ll[0], ll[1], ll[0]+ll[1], ll[2]);
                        for (i = 0; i < xy_size*2; i++)
                                inflow[i] = flow[i];
                }
                free(scratch);
        }

        if (outflow_file != NULL) {  
                if ((fp_flow = fopen(outflow_file, "wb")) == NULL) {
                        fprintf(stderr, "Error: Couldn't write file %s.\n",
                                outflow_file);
                        exit(EXIT_FAILURE);
                }
                fwrite(&size_curv, 2, sizeof(int), fp_flow);
                fwrite(&shift, 2, sizeof(int), fp_flow);
                fwrite(inflow, xy_size*2, sizeof(double), fp_flow);
                fclose(fp_flow);
        }

        /* get deformations and jacobian det. from flow field */
        if (jacdet_file != NULL) {
                jd  = (double *) malloc(sizeof(double) * xy_size);
                jd1 = (double *) malloc(sizeof(double) * xy_size);

                expdefdet(size_curv, 10, inflow, flow, flow1, jd, jd1);
                
                /* subtract 1 to get values around 0 instead of 1 */
                for (i = 0; i < xy_size; i++) {
                        jd1[i] -= 1;
                }

                values = (double *) malloc(sizeof(double) *
                                           source->n_points);

                map_sheet2d_to_sphere(jd1, values, source,
                                      1, size_curv);

                output_values_any_format(jacdet_file, source->n_points,
                                         values, TYPE_DOUBLE);

                free(values);
                free(jd);
                free(jd1);
        } else {
                expdef(size_curv, 10, inflow, flow, flow1,
                       (double *) 0, (double *) 0);
        }

        /* get amplitude of deformations from flow field */
        if (deform_file != NULL) {
                deform  = (double *) malloc(sizeof(double) * xy_size);

                for (i = 0; i < xy_size; i++) {
                        x = (int) flow[i] - 1.0;
                        y = (int) flow[i + xy_size] - 1.0;
                        xp = flow[i] - 1.0 - x;
                        yp = flow[i + xy_size] - 1.0 - y;
                        xm = 1.0 - xp;
                        ym = 1.0 - yp;
                        deform[i] = sqrt(xp*xp + yp*yp);
                }
                
                values = (double *) malloc(sizeof(double) *
                                           source->n_points);

                map_sheet2d_to_sphere(deform, values, source, 1, size_curv);

                output_values_any_format(deform_file, source->n_points,
                                         values, TYPE_DOUBLE);

                free(values);
                free(deform);
        }

        free(flow1);
        free(inflow);

        if (pgm_file != NULL) {
                if (weight_file != NULL) { /* rescue old unweighted map data */
                        for (i = 0; i < xy_size; i++)
                                map_source[i] = map_source0[i];
                        free(map_source0);
                }
                for (i = 0; i < xy_size; i++) {
                        x = (int) flow[i] - 1.0;
                        y = (int) flow[i + xy_size] - 1.0;
                        xp = flow[i] - 1.0 - x;
                        yp = flow[i + xy_size] - 1.0 - y;
                        xm = 1.0 - xp;
                        ym = 1.0 - yp;
                        H00 = map_source[bound(x,  y,  size_curv)];
                        H01 = map_source[bound(x,  y+1,size_curv)];
                        H10 = map_source[bound(x+1,y,  size_curv)];
                        H11 = map_source[bound(x+1,y+1,size_curv)];

                        map_warp[i] = ym * (xm * H00 + xp * H10) +
		                      yp * (xm * H01 + xp * H11);
                }
                if (write_pgm(pgm_file, map_warp, size_curv[0],
                              size_curv[1]) != 0)
                        exit(EXIT_FAILURE);
        }

        if (output_file != NULL) {
                apply_warp(source, source_sphere, flow, size_curv, shift);  
  
                if (output_graphics_any_format(output_file, format, n_objects,
                                               objects) != OK)
                        exit(EXIT_FAILURE);
        }

        delete_object_list(n_objects, objects);
        free(map_source);
        free(map_warp);
        free(map_target);
        free(flow);
        free(prm);

        return(EXIT_SUCCESS);
}
