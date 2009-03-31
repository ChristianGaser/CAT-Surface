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
  
struct dartel_prm {
  int rtype;         /* regularization type: 0 - linear elastic energy; */
                     /* 1 - membrane energy; 2 - bending energy */
  double rparam[5];  /* regularization parameters */
  double lmreg;      /* LM regularization */
  int cycles;        /* # of cycles for full multi grid (FMG) */
  int its;           /* # of relaxation iterations in each multigrid cycle */
  int k;             /* time steps for solving the PDE */
  int code;          /* objective function: 0 - sum of squares; */
                     /* 1 - symmetric sum of squares */
};

/* defaults */
char *param_file   = NULL;
char *source_file  = NULL;
char *target_file  = NULL;
char *weight_file  = NULL;
char *jacdet_file  = NULL;
char *output_file  = NULL;
char *pgm_file     = NULL;
char *outflow_file = NULL;
char *inflow_file  = NULL;
int translate = 0;
int code = 1;
int loop = 6;
int verbose = 0;
int rtype = 1;
int curvtype = 0;
double reg   = 0.0001;
double lmreg = 0.0001;
double fwhm  = 10.0;

static ArgvInfo argTable[] = {
  {"-i", ARGV_STRING, (char *) 1, (char *) &source_file, 
     "Input file."},
  {"-t", ARGV_STRING, (char *) 1, (char *) &target_file, 
     "Template file."},
  {"-j", ARGV_STRING, (char *) 1, (char *) &jacdet_file, 
     "Jacobian determinant values on the surface."},
  {"-w", ARGV_STRING, (char *) 1, (char *) &output_file, 
     "Warped input."},
  {"-o", ARGV_STRING, (char *) 1, (char *) &pgm_file, 
     "Warped map as pgm file."},
  {"-s", ARGV_STRING, (char *) 1, (char *) &weight_file, 
     "Weight warping with the inverse from values in this file (e.g. std)."},
  {"-if", ARGV_STRING, (char *) 1, (char *) &inflow_file, 
     "Warped input."},
  {"-of", ARGV_STRING, (char *) 1, (char *) &outflow_file, 
     "Warped input."},
  {"-p", ARGV_STRING, (char *) 1, (char *) &param_file, 
     "Parameter file."},
  {"-code", ARGV_INT, (char *) 1, (char *) &code,
     "Objective function (code): 0 - sum of squares; 1 - symmetric sum of squares."},
  {"-rtype", ARGV_INT, (char *) 1, (char *) &rtype,
     "Regularization type: 0 - linear elastic energy; 1 - membrane energy; 2 - bending energy."},
  {"-reg", ARGV_FLOAT, (char *) 1, (char *) &reg,
     "Regularization parameter."},
  {"-lmreg", ARGV_FLOAT, (char *) 1, (char *) &lmreg,
     "LM regularization."},
  {"-fwhm", ARGV_FLOAT, (char *) 1, (char *) &fwhm,
     "Filter size for curvature map in FWHM."},
  {"-loop", ARGV_INT, (char *) 1, (char *) &loop,
     "Number of outer loops for default parameters (max. 6)."},
  {"-shift", ARGV_CONSTANT, (char *) TRUE, (char *) &translate,
     "Shift map before warping."},
  {"-type", ARGV_INT, (char *) 1, (char *) &curvtype,
     "Objective function (code): 0 - sum of squares; 1 - symmetric sum of squares."},
  {"-v", ARGV_CONSTANT, (char *) TRUE, (char *) &verbose,
     "Be verbose."},
   {NULL, ARGV_END, NULL, NULL, NULL}
};

void
translate_to_template(double *source, double *target, int *size_map,
                      int *shift)
{
        int     mn = -10, mx = 10;
        int     i, sx, sy, x, y;
        double  *shifted, sdiff, diff, sdiff0 = FLT_MAX;
  
        shifted = (double *) malloc(sizeof(double) * size_map[0] * size_map[1]);

        /* restricted to shifts on x-scale only */
        for (sx = mn; sx < mx; sx++) {
                for (sy = 0; sy < 1; sy++) {
                        sdiff = 0.0;
                        for (x = 0; x < size_map[0]; x++) {
                                for (y = 0; y < size_map[1]; y++) {
                                        i = x + (size_map[0] * y);
                                        diff = target[i] -
                                               source[bound(x + sx, y + sy,
                                                            size_map)];
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
  
        for (x = 0; x < size_map[0]; x++) {
                for (y = 0; y < size_map[1]; y++) {
                        i = x + (size_map[0] * y);
                        shifted[i] = source[bound(x + shift[0], y + shift[1],
                                                  size_map)];
                }
        }

        for (i = 0; i < size_map[0]*size_map[1]; i++)
                source[i] = shifted[i];

        free(shifted);
}


int
main(int argc, char *argv[])
{
        File_formats     format;
        FILE             *fp, *fp_flow;
        char             line[1024];
        polygons_struct  *polygons_source, *polygons_target;
        int              x, y, i, j, it, it0, it1;
        int              n_objects, it_scratch, xy_size, n_weights;
        double           *weights;
        double           *map_source, *map_target;
        double           *map_warp, *map_weights, *map_source0;
        double           *flow, *flow1, *inflow, *scratch, *jd, *jd1, *values;
        object_struct    **objects, *object;
        double           ll[3];
        static double    param[3] = {1.0, 1.0, 0.25};
        int              size_map[3] = {512, 256, 1}, shift[2];
        double           xp, yp, xm, ym;
        double           H00, H01, H10, H11;
        struct dartel_prm* prm;

        /* get the arguments from the command line */

        if (ParseArgv(&argc, argv, argTable, 0) ||
            source_file == NULL || target_file == NULL ||
            (jacdet_file == NULL && output_file == NULL &&
             pgm_file == NULL  && outflow_file == NULL)) {
                fprintf(stderr, "\nUsage: %s [options]\n", argv[0]);
                fprintf(stderr, "     %s -help\n\n", argv[0]);
                exit(EXIT_FAILURE);
        }

        if (input_graphics_any_format(target_file, &format,
                                      &n_objects, &objects) != OK)
                return(1);

        /* check that the surface file contains a polyhedron */
        if (n_objects != 1 || get_object_type(objects[0]) != POLYGONS) {
                printf("Surface file must contain 1 polygons object.\n");
                return(1);
        }
        /* get a pointer to the surface */
        polygons_target = get_polygons_ptr(objects[0]);

        if (input_graphics_any_format(source_file, &format,
                                      &n_objects, &objects) != OK)
                return(1);

        /* check that the surface file contains a polyhedron */
        if (n_objects != 1 || get_object_type(objects[0]) != POLYGONS) {
                printf("Surface file must contain 1 polygons object.\n");
                return(1);
        }
  
        /* read weights */
        if (weight_file != NULL) {
                if (input_values_any_format(weight_file, &n_weights,
                                         &weights) != OK)
                        return(1);
        }

        prm = (struct dartel_prm*) malloc(sizeof(struct dartel_prm) * 100);

        /* first three entrys of param are equal */
        for (j = 0; j < loop; j++) {
                for (i = 0; i < 3; i++)
                        prm[j].rparam[i] = param[i];
        }

        /* read values from parameter file */
        if (param_file != NULL) {
                if ((fp = fopen(param_file, "r")) == 0) {
                        fprintf(stderr, "Couldn't open parameter file %s.\n",
                                param_file);
                        return(0);
                }
    
                loop = 0;

                fprintf(stderr, "Read parameters from %s\n", param_file);
                while (fgets(line, sizeof(line), fp)) {
                        /* check for 9 values in each line */
                        if (sscanf(line, "%d %lf %lf %lf %d %d %d %d",
                                   &prm[loop].rtype, &prm[loop].rparam[3],
                                   &prm[loop].rparam[4], &prm[loop].lmreg,
                                   &prm[loop].cycles, &prm[loop].its,
                                   &prm[loop].k, &prm[loop].code) != 8)
                                continue;
                        loop++;
                }
                fclose(fp);

                if (loop == 0) {
                        fprintf(stderr, "Could not read parameter file %s. Check that each line contains 9 values\n", param_file);
                        return(0);
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
                        prm[i].rparam[3] = reg;
                        prm[i].rparam[4] = reg/2;
                        prm[i].k = i;
                        reg /= 5;
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
                fprintf(stderr, "\nRegularization parameters 1:\t");
                for (i = 0; i < loop; i++)
                        fprintf(stderr, "%8g\t", prm[i].rparam[3]);
                fprintf(stderr,"\n");
                fprintf(stderr, "Regularization parameters 2:\t");
                for (i = 0; i < loop; i++)
                        fprintf(stderr, "%8g\t", prm[i].rparam[4]);
                fprintf(stderr,"\n");
                fprintf(stderr, "Time steps for solving the PDE:\t");
                for (i = 0; i < loop; i++)
                        fprintf(stderr,"%8d\t",prm[i].k);
                fprintf(stderr,"\n\n");
        }
  
        /* get a pointer to the surface */
        polygons_source = get_polygons_ptr(objects[0]);
  
        xy_size = size_map[0] * size_map[1];

        flow        = (double *) malloc(sizeof(double) * xy_size * 2);
        flow1       = (double *) malloc(sizeof(double) * xy_size * 2);
        inflow      = (double *) malloc(sizeof(double) * xy_size * 2);
        map_source  = (double *) malloc(sizeof(double) * xy_size);
        map_target  = (double *) malloc(sizeof(double) * xy_size);
        map_warp    = (double *) malloc(sizeof(double) * xy_size);

        map_smoothed_curvature_to_sphere(polygons_source, (double *)0,
                                         map_source, fwhm, size_map, curvtype);
        map_smoothed_curvature_to_sphere(polygons_target, (double *)0,
                                         map_target, fwhm, size_map, curvtype);

        if (write_pgm("source.pgm", map_source, size_map[0], size_map[1]) != 0)
                return(1);

        if (write_pgm("target.pgm", map_target, size_map[0], size_map[1]) != 0)
                return(1);

        /* weight maps of target and source by the inverse std */
        if (weight_file != NULL) {
                map_weights = (double *) malloc(sizeof(double) * xy_size);
                map_source0 = (double *) malloc(sizeof(double) * xy_size);
                map_smoothed_curvature_to_sphere(polygons_source, weights,
                                                 map_weights, fwhm, size_map, curvtype);
                for (i = 0; i < xy_size; i++) {
                        map_source0[i] = map_source[i];
                        map_weights[i] += 1.0;
                        map_source[i] /= map_weights[i];
                        map_target[i] /= map_weights[i];
                }
                free(map_weights);
        }
  
        if (translate) {
                translate_to_template(map_source, map_target, size_map, shift);
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
                        return(1);
                }
                fprintf(stderr,"Shift is not considered!!!\n");
                fread(&size_map, 2, sizeof(int), fp_flow);
                fread(&shift, 2, sizeof(int), fp_flow);
                fread(inflow, xy_size*2, sizeof(double), fp_flow);
                fclose(fp_flow);
                size_map[2] = 1;
        } else {
                for (i = 0; i < xy_size*2; i++) {
                        inflow[i] = 0.0;
                }
        }

        for (it = 0, it0 = 0; it0 < loop; it0++) {
                it_scratch = dartel_scratchsize((int *)size_map, prm[it0].code);
                scratch = (double *) malloc(sizeof(double) * it_scratch);

                for (it1 = 0; it1 < prm[it0].its; it1++) {
                        it++;
                        dartel(size_map, prm[it0].k, inflow, map_target,
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
                        return(1);
                }
                fwrite(&size_map, 2, sizeof(int), fp_flow);
                fwrite(&shift, 2, sizeof(int), fp_flow);
                fwrite(inflow, xy_size*2, sizeof(double), fp_flow);
                fclose(fp_flow);
        }

        /* get deformations and jacobian det. from flow field */
        if (jacdet_file != NULL) {
                jd  = (double *) malloc(sizeof(double) * xy_size);
                jd1 = (double *) malloc(sizeof(double) * xy_size);

                expdefdet(size_map, 10, inflow, flow, flow1, jd, jd1);

                values = (double *) malloc(sizeof(double) *
                                           polygons_source->n_points);

                map_sheet2d_to_sphere(jd1, values, polygons_source,
                                      1, size_map);

                output_values_any_format(jacdet_file, polygons_source->n_points, values);

                free(values);
                free(jd);
                free(jd1);
        } else {
                expdef(size_map, 10, inflow, flow, flow1,
                       (double *) 0, (double *) 0);
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
                        H00 = map_source[bound(x,  y,  size_map)];
                        H01 = map_source[bound(x,  y+1,size_map)];
                        H10 = map_source[bound(x+1,y,  size_map)];
                        H11 = map_source[bound(x+1,y+1,size_map)];

                        map_warp[i] = ym * (xm * H00 + xp * H10) +
		                      yp * (xm * H01 + xp * H11);
                }
                if (write_pgm(pgm_file, map_warp, size_map[0],
                              size_map[1]) != 0)
                        return(1);
        }

        if (output_file != NULL) {
                apply_warp(polygons_source, flow, size_map, shift);  
  
                if (output_graphics_any_format(output_file, format, n_objects,
                                               objects) != OK)
                        return(1);
        }

        delete_object_list(n_objects, objects);
        free(map_source);
        free(map_warp);
        free(map_target);
        free(flow);
        free(prm);

        return(0);
}
