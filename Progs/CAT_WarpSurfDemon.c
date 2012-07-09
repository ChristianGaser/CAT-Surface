/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 *
 */

#include <bicpl.h>
#include <float.h>
#include <ParseArgv.h>

#include "CAT_Map.h"
#include "CAT_Surf.h"
#include "CAT_Curvature.h"
#include "CAT_SurfaceIO.h"
#include "CAT_Interpolate.h"
#include "dartel/dartel.h"

#define INVERSE_WARPING 0
#define RADIANS(deg) ((PI * (double)(deg)) / 180.0)
#define DEGREES(rad) ((180.0 * (double)(rad)) / PI)


/* defaults */
char *param_file         = NULL;
char *source_file        = NULL;
char *source_sphere_file = NULL;
char *target_file        = NULL;
char *target_sphere_file = NULL;
char *output_file        = NULL;
char *output_sphere_file = NULL;

int rotate      = 1;
int code        = 1;
int loop        = 6;
int verbose     = 0;
int rtype       = 1;
int curvtype    = 3;
int muchange    = 4;
int n_triangles = 81920;
int n_steps     = 1;
int debug       = 0;
double murate   = 1.25;
double lambda   = 0.0;
double mu       = 0.25;
double lmreg    = 0.0;
double fwhm     = 10.0;

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
  {"-ws", ARGV_STRING, (char *) 1, (char *) &output_sphere_file, 
     "Warped input sphere."},
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
  {"-steps", ARGV_INT, (char *) 1, (char *) &n_steps,
     "Number of Dartel steps (max. 6)."},
  {"-norot", ARGV_CONSTANT, (char *) FALSE, (char *) &rotate,
     "Don't rotate input surface before warping."},
  {"-type", ARGV_INT, (char *) 1, (char *) &curvtype,
     "Curvature type\n\t0 - mean curvature (averaged over 3mm, in degrees)\n\t1 - gaussian curvature\n\t2 - curvedness\n\t3 - shape index\n\t4 - mean curvature (in radians)."},
  {"-v", ARGV_CONSTANT, (char *) TRUE, (char *) &verbose,
     "Be verbose."},
  {"-debug", ARGV_CONSTANT, (char *) TRUE, (char *) &debug,
     "Save debug files."},
   {NULL, ARGV_END, NULL, NULL, NULL}
};

void
gradient_poly(polygons_struct *polygons, struct dartel_poly *dpoly,
             double f[], double dtheta[], double dphi[])
{
    int i, mm = polygons->n_points;
    double kxm, kxp, kym, kyp;

    for (i = 0; i < mm; i++) {

        kxm   = interp_point_unit_sphere(polygons, f, dpoly->ntheta[i]);
        kxp   = interp_point_unit_sphere(polygons, f, dpoly->ptheta[i]);
        kym   = interp_point_unit_sphere(polygons, f, dpoly->nphi[i]);
        kyp   = interp_point_unit_sphere(polygons, f, dpoly->pphi[i]);
        dtheta[i] = (kxp - kxm)/2.0;
        dphi[i]   = (kyp - kym)/2.0;

    }

    if (polygons->bintree != NULL) delete_the_bintree(&polygons->bintree);
}

void
rotate_polygons(polygons_struct *polygons, polygons_struct *rotated_polygons,
                double *rotation_matrix)
{
        int i;
        double x, y, z;
        
        if (rotated_polygons != NULL)
                copy_polygons(polygons, rotated_polygons);
        
        for (i = 0; i < polygons->n_points; i++) {
                x = Point_x(polygons->points[i])*rotation_matrix[0] 
                  + Point_y(polygons->points[i])*rotation_matrix[1]
                  + Point_z(polygons->points[i])*rotation_matrix[2];
                y = Point_x(polygons->points[i])*rotation_matrix[3] 
                  + Point_y(polygons->points[i])*rotation_matrix[4]
                  + Point_z(polygons->points[i])*rotation_matrix[5];
                z = Point_x(polygons->points[i])*rotation_matrix[6] 
                  + Point_y(polygons->points[i])*rotation_matrix[7]
                  + Point_z(polygons->points[i])*rotation_matrix[8];
                if (rotated_polygons != NULL) {
                        fill_Point(rotated_polygons->points[i], x, y, z);
                } else  fill_Point(polygons->points[i], x, y, z);
        }
}

void
rotation_to_matrix(double *rotation_matrix, double alpha, double beta,
                   double gamma)
{
        int i, j, k;        
        double sum, rot[9];
        
        /* rotation matrices */
        double rot_x[9] = {1.0, 0.0,         0.0,
                           0.0, cos(alpha),  sin(alpha),
                           0.0, -sin(alpha), cos(alpha)}; 
        double rot_y[9] = {cos(beta),  0.0, sin(beta),
                           0.0,        1.0, 0.0,
                           -sin(beta), 0.0, cos(beta)}; 
        double rot_z[9] = {cos(gamma),  sin(gamma), 0.0,
                           -sin(gamma), cos(gamma), 0.0,
                           0.0,         0.0,        1.0}; 

        /* combine x and y r tation */
        for (i = 0; i < 3; i++) {
                for (j = 0; j < 3; j++) {
                        sum = 0.0;
                        for (k = 0; k < 3; k++)
                                sum += rot_y[i + 3*k] * rot_x[k + 3*j];
                        rot[i + 3*j] = sum;
                }
        }

        /* combine with z rotation */
        for (i = 0; i < 3; i++) {
                for (j = 0; j < 3; j++) {
                        sum = 0.0;
                        for (k = 0; k < 3; k++) 
                                sum += rot_z[i + 3*k] * rot[k + 3*j];
                        rotation_matrix[i + 3*j] = sum;
                }
        }
}

int
main(int argc, char *argv[])
{
        File_formats     format;
        FILE             *fp;
        char             line[1024], buffer[1024];
        polygons_struct  *source, *target, *src_sphere, *trg_sphere;
        polygons_struct  *sm_source, *sm_target;
        int              x, y, i, j, it, it0, it1, step;
        int              n_objects, it_scratch, xy_size;
        double           *values;
        object_struct    **objects;
        double           ll[3];
        static double    param[2] = {UTHETA, VPHI};
        int              size_curv[3], shift[2] = {0, 0};
        double           xp, yp, xm, ym, idiff, denom;
        double           H00, H01, H10, H11, rotation_matrix[9];
        struct           dartel_prm* prm;
        double           rot[3], *curv_tmp, *curv_target, *curv_source, *dtheta, *dphi, *u, *v;
        struct           dartel_poly *dpoly;
        int              *n_neighbours, **neighbours;

        /* get the arguments from the command line */

        if (ParseArgv(&argc, argv, argTable, 0) ||
            source_file == NULL || target_file == NULL || 
            source_sphere_file == NULL || target_sphere_file == NULL ||
            (output_file == NULL && output_sphere_file == NULL)) {
                fprintf(stderr, "\nUsage: %s [options]\n", argv[0]);
                fprintf(stderr, "     %s -help\n\n", argv[0]);
                exit(EXIT_FAILURE);
        }

        if (input_graphics_any_format(target_file, &format,
                                      &n_objects, &objects) != OK)
                exit(EXIT_FAILURE);

        /* get a pointer to the surface */
        target = get_polygons_ptr(objects[0]);

        /* check that the surface file contains a polyhedron */
        if (n_objects != 1 || get_object_type(objects[0]) != POLYGONS) {
                printf("Surface file must contain 1 polygons object.\n");
                exit(EXIT_FAILURE);
        }

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
                        prm[j].cycles = 1; //3;
                        prm[j].its = 1; //3;
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
    
        xy_size = source->n_points;

        dpoly = (struct dartel_poly *) malloc(sizeof(struct dartel_poly));
        init_dartel_poly(src_sphere, dpoly);

        sm_source   = (polygons_struct *) malloc(sizeof(polygons_struct));
        sm_target   = (polygons_struct *) malloc(sizeof(polygons_struct));

        curv_tmp    = (double *) malloc(sizeof(double) * target->n_points);
        curv_source = (double *) malloc(sizeof(double) * source->n_points);
        curv_target = (double *) malloc(sizeof(double) * source->n_points);
        dtheta      = (double *) malloc(sizeof(double) * source->n_points);
        dphi        = (double *) malloc(sizeof(double) * source->n_points);
        u           = (double *) malloc(sizeof(double) * source->n_points);
        v           = (double *) malloc(sizeof(double) * source->n_points);
        
        get_all_polygon_point_neighbours(source, &n_neighbours, &neighbours);
        get_polygon_vertex_curvatures_cg(source, n_neighbours, neighbours,
                                         0.0, curvtype, curv_source);

        get_all_polygon_point_neighbours(target, &n_neighbours, &neighbours);
        get_polygon_vertex_curvatures_cg(target, n_neighbours, neighbours,
                                         0.0, curvtype, curv_tmp);

        /* resample curvature of target to source space */
        resample_noscale(trg_sphere, src_sphere, curv_tmp, curv_target);
                 
        init_dartel_poly(trg_sphere, dpoly);

        for (step = 0; step < n_steps; step++) {
                /* resample source and target surface */
                copy_polygons(source, sm_source);
                copy_polygons(target, sm_target);

                /* initialization */
                if (step == 0) {
                        /* inflate surfaces */
                        inflate_surface_and_smooth_fingers(sm_source,
                                                           1, 0.2, 50, 1.0,
                                                           3.0, 1.0, 0);
                        inflate_surface_and_smooth_fingers(sm_target,
                                                           1, 0.2, 50, 1.0,
                                                           3.0, 1.0, 0);
                        inflate_surface_and_smooth_fingers(sm_source,
                                                           2, 1.0, 30, 1.4,
                                                           3.0, 1.0, 0);
                        inflate_surface_and_smooth_fingers(sm_target,
                                                           2, 1.0, 30, 1.4,
                                                           3.0, 1.0, 0);
                        
                } else if (step == 1) {
                        /* smooth surfaces */
                        inflate_surface_and_smooth_fingers(sm_source,
                                                           1, 0.2, 30, 1.0,
                                                           3.0, 1.0, 0);
                        inflate_surface_and_smooth_fingers(sm_target,
                                                           1, 0.2, 30, 1.0,
                                                           3.0, 1.0, 0);
                        inflate_surface_and_smooth_fingers(sm_source,
                                                           2, 1.0, 10, 1.4,
                                                           3.0, 1.0, 0);
                        inflate_surface_and_smooth_fingers(sm_target,
                                                           2, 1.0, 10, 1.4,
                                                           3.0, 1.0, 0);
                } else if (step == 2) {
                        /* smooth surfaces */
                        inflate_surface_and_smooth_fingers(sm_source,
                                                           1, 0.2, 20, 1.0,
                                                           3.0, 1.0, 0);
                        inflate_surface_and_smooth_fingers(sm_target,
                                                           1, 0.2, 20, 1.0,
                                                           3.0, 1.0, 0);
                }
                
                gradient_poly(trg_sphere, dpoly, curv_target, dtheta, dphi);
                for (i = 0; i < source->n_points; i++) {
                        idiff = curv_source[i] - curv_target[i];
                        denom = -((dtheta[i]*dtheta[i] + dphi[i]*dphi[i]) + idiff*idiff);
                        if (denom == 0.0) {
                                u[i] = 0.0;
                                v[i] = 0.0;
                        } else {
                                u[i] = idiff*dtheta[i]/denom;
                                v[i] = idiff*dphi[i]/denom;                        
                        }
                }
                
                /* lowpass filter displacements */
                smooth_heatkernel(source, u, 5);
                smooth_heatkernel(source, v, 5);

                /* use smaller FWHM for next steps */
                fwhm /= 2.0;
        }
        fprintf(stderr,"\n");
        
        output_values_any_format("u.txt", source->n_points, u, TYPE_DOUBLE);
        output_values_any_format("v.txt", source->n_points, v, TYPE_DOUBLE);
        output_values_any_format("curv.txt", source->n_points, curv_target, TYPE_DOUBLE);

        if (output_sphere_file != NULL) {
  
                /* get a pointer to the surface */
                *get_polygons_ptr(objects[0]) = *src_sphere;
                if (output_graphics_any_format(output_sphere_file, format,
                                               n_objects, objects) != OK)
                        exit(EXIT_FAILURE);
        }

        delete_object_list(n_objects, objects);
        free(prm);

        return(EXIT_SUCCESS);
}
