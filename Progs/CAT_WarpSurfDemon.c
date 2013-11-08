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
#include "CAT_Resample.h"
#include "dartel/dartel.h"

/* defaults */
char *src_file        = NULL;
char *sphere_src_file = NULL;
char *trg_file        = NULL;
char *sphere_trg_file = NULL;
char *output_file        = NULL;
char *output_sphere_file = NULL;

int rotate      = 1;
int curvtype    = 3;
int n_steps     = 1;
int debug       = 0;
double fwhm     = 20.0;

static ArgvInfo argTable[] = {
  {"-i", ARGV_STRING, (char *) 1, (char *) &src_file, 
     "Input file."},
  {"-is", ARGV_STRING, (char *) 1, (char *) &sphere_src_file, 
     "Input sphere file."},
  {"-t", ARGV_STRING, (char *) 1, (char *) &trg_file, 
     "Template file."},
  {"-ts", ARGV_STRING, (char *) 1, (char *) &sphere_trg_file, 
     "Template sphere file."},
  {"-w", ARGV_STRING, (char *) 1, (char *) &output_file, 
     "Warped brain."},
  {"-ws", ARGV_STRING, (char *) 1, (char *) &output_sphere_file, 
     "Warped input sphere."},
  {"-fwhm", ARGV_FLOAT, (char *) 1, (char *) &fwhm,
     "Filter size for curvature map in FWHM."},
  {"-steps", ARGV_INT, (char *) 1, (char *) &n_steps,
     "Number of multigrid steps."},
  {"-norot", ARGV_CONSTANT, (char *) FALSE, (char *) &rotate,
     "Don't rotate input surface before warping."},
  {"-type", ARGV_INT, (char *) 1, (char *) &curvtype,
     "Curvature type\n\t0 - mean curvature (averaged over 3mm, in degrees)\n\t1 - gaussian curvature\n\t2 - curvedness\n\t3 - shape index\n\t4 - mean curvature (in radians)."},
  {"-debug", ARGV_CONSTANT, (char *) TRUE, (char *) &debug,
     "Save debug files."},
   {NULL, ARGV_END, NULL, NULL, NULL}
};

/* not used */
void
scale_values(double values[], int n_values)
{
    int i;
    double mn = FLT_MAX, mx = -FLT_MAX;

    for (i = 0; i < n_values; i++) {
        mn = MIN(values[i], mn);
        mx = MAX(values[i], mx);
    }

    for (i = 0; i < n_values; i++) 
        values[i] = (values[i] - mn)/(mx - mn);
        
}

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

/* not used */
void
laplace_poly(polygons_struct *polygons, struct dartel_poly *dpoly,
             double f[], double laplace[])
{
    int i, mm = polygons->n_points;
    double kxm, kxp, kym, kyp;

    for (i = 0; i < mm; i++) {

        kxm   = interp_point_unit_sphere(polygons, f, dpoly->ntheta[i]);
        kxp   = interp_point_unit_sphere(polygons, f, dpoly->ptheta[i]);
        kym   = interp_point_unit_sphere(polygons, f, dpoly->nphi[i]);
        kyp   = interp_point_unit_sphere(polygons, f, dpoly->pphi[i]);
        laplace[i] = kxp + kxm + kyp + kym - 4*f[i];

    }

    if (polygons->bintree != NULL) delete_the_bintree(&polygons->bintree);
}

/* not used */
void
gradient_vector_flow(polygons_struct *polygons, struct dartel_poly *dpoly, double f[], double dtheta[], double dphi[], double u[], double v[], double mu, int iterations)
{
    int i, j, mm = polygons->n_points;
    double fmin = FLT_MAX, fmax = -FLT_MAX;
    double *b, *c1, *c2, *Lu, *Lv;

    b  = (double *) malloc(sizeof(double) * polygons->n_points);
    c1 = (double *) malloc(sizeof(double) * polygons->n_points);
    c2 = (double *) malloc(sizeof(double) * polygons->n_points);
    Lu = (double *) malloc(sizeof(double) * polygons->n_points);
    Lv = (double *) malloc(sizeof(double) * polygons->n_points);

    for (i = 0; i < mm; i++) {
        fmin = MIN(f[i], fmin);
        fmax = MAX(f[i], fmax);
    }

    /* scale f to a range of 0..1 and initialize parameters */
    for (i = 0; i < mm; i++) {
        f[i]  = (f[i] - fmin)/(fmax - fmin);
        u[i]  = dtheta[i];
        v[i]  = dphi[i];
        b[i]  = dtheta[i]*dtheta[i] + dphi[i]*dphi[i];
        c1[i] = b[i]*dtheta[i];
        c2[i] = b[i]*dphi[i];
    }
      
    for (j = 0; j < iterations; j++) {
        laplace_poly(polygons, dpoly, u, Lu);
        laplace_poly(polygons, dpoly, v, Lv);

        for (i = 0; i < mm; i++) {
            u[i] = (1 - b[i])*u[i] + mu*Lu[i] + c1[i];
            v[i] = (1 - b[i])*v[i] + mu*Lv[i] + c2[i];
        }
        
        smooth_heatkernel(polygons, u, fwhm);
        smooth_heatkernel(polygons, v, fwhm);

        output_values_any_format("u.txt", polygons->n_points, u, TYPE_DOUBLE);
        output_values_any_format("v.txt", polygons->n_points, v, TYPE_DOUBLE);
        fprintf(stderr, "GVF iteration %d\n", j);
    }
    
    free(b);
    free(c1);
    free(c2);
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
        polygons_struct  *src, *trg, *sphere_src, *sphere_trg;
        polygons_struct  *sm_src, *sm_trg, *warped_sphere_src;
        int              i, step;
        int              n_objects, xy_size;
        object_struct    **objects;
        int              iter_demon, count_break;
        double           idiff, denom_trg, denom_src, squared_diff, old_squared_diff;
        double           *curv_trg_tmp, *curv_trg, *curv_src;
        double           *curv_src_tmp, *dtheta_trg, *dphi_trg, *dtheta_src, *dphi_src, *u, *v, *Utheta, *Uphi;
        struct           dartel_poly *dpoly_src, *dpoly_trg;
        int              *n_neighbours, **neighbours;

        /* get the arguments from the command line */

        if (ParseArgv(&argc, argv, argTable, 0) ||
            src_file == NULL || trg_file == NULL || 
            sphere_src_file == NULL || sphere_trg_file == NULL ||
            (output_file == NULL && output_sphere_file == NULL)) {
                fprintf(stderr, "\nUsage: %s [options]\n", argv[0]);
                fprintf(stderr, "     %s -help\n\n", argv[0]);
                exit(EXIT_FAILURE);
        }

        if (input_graphics_any_format(trg_file, &format,
                                      &n_objects, &objects) != OK)
                exit(EXIT_FAILURE);

        /* get a pointer to the surface */
        trg = get_polygons_ptr(objects[0]);

        /* check that the surface file contains a polyhedron */
        if (n_objects != 1 || get_object_type(objects[0]) != POLYGONS) {
                printf("Surface file must contain 1 polygons object.\n");
                exit(EXIT_FAILURE);
        }

        if (input_graphics_any_format(src_file, &format,
                                      &n_objects, &objects) != OK)
                exit(EXIT_FAILURE);

        /* get a pointer to the surface */
        src = get_polygons_ptr(objects[0]);

        /* check that the surface file contains a polyhedron */
        if (n_objects != 1 || get_object_type(objects[0]) != POLYGONS) {
                printf("Surface file must contain 1 polygons object.\n");
                exit(EXIT_FAILURE);
        }
  
        /* read sphere for input surface */
        if (input_graphics_any_format(sphere_src_file, &format,
                              &n_objects, &objects) != OK)
                exit(EXIT_FAILURE);
        sphere_src = get_polygons_ptr(objects[0]);

        /* read sphere for template surface */
        if (input_graphics_any_format(sphere_trg_file, &format,
                              &n_objects, &objects) != OK)
                exit(EXIT_FAILURE);
        sphere_trg = get_polygons_ptr(objects[0]);

        translate_to_center_of_mass(sphere_src);
        for (i = 0; i < sphere_src->n_points; i++)
                set_vector_length(&sphere_src->points[i], 1.0);
        translate_to_center_of_mass(sphere_trg);
        for (i = 0; i < sphere_trg->n_points; i++)
                set_vector_length(&sphere_trg->points[i], 1.0);

        xy_size = src->n_points;

        dpoly_src         = (struct dartel_poly *) malloc(sizeof(struct dartel_poly));
        dpoly_trg         = (struct dartel_poly *) malloc(sizeof(struct dartel_poly));

        sm_src            = (polygons_struct *) malloc(sizeof(polygons_struct));
        sm_trg            = (polygons_struct *) malloc(sizeof(polygons_struct));
        warped_sphere_src = (polygons_struct *) malloc(sizeof(polygons_struct));

        curv_trg_tmp      = (double *) malloc(sizeof(double) * trg->n_points);
        curv_src_tmp      = (double *) malloc(sizeof(double) * src->n_points);
        curv_src          = (double *) malloc(sizeof(double) * src->n_points);
        curv_trg          = (double *) malloc(sizeof(double) * src->n_points);
        dtheta_trg        = (double *) malloc(sizeof(double) * src->n_points);
        dphi_trg          = (double *) malloc(sizeof(double) * src->n_points);
        dtheta_src       = (double *) malloc(sizeof(double) * src->n_points);
        dphi_src          = (double *) malloc(sizeof(double) * src->n_points);
        u                 = (double *) malloc(sizeof(double) * src->n_points);
        v                 = (double *) malloc(sizeof(double) * src->n_points);
        Utheta            = (double *) malloc(sizeof(double) * src->n_points);
        Uphi              = (double *) malloc(sizeof(double) * src->n_points);
        
        init_dartel_poly(sphere_src, dpoly_src);
        init_dartel_poly(sphere_trg, dpoly_trg);

        copy_polygons(sphere_src, warped_sphere_src);

        for (i = 0; i < src->n_points; i++) {
                u[i] = 0.0;
                v[i] = 0.0;
        }
        
        for (step = 0; step < n_steps; step++) {
                /* resample src and trg surface */
                copy_polygons(src, sm_src);
                copy_polygons(trg, sm_trg);

                /* initialization */
                if (step == 0) {
                        /* inflate surfaces */
                        inflate_surface_and_smooth_fingers(sm_src,
                                                           1, 0.2, 50, 1.0,
                                                           3.0, 1.0, 0);
                        inflate_surface_and_smooth_fingers(sm_trg,
                                                           1, 0.2, 50, 1.0,
                                                           3.0, 1.0, 0);
                        inflate_surface_and_smooth_fingers(sm_src,
                                                           2, 1.0, 30, 1.4,
                                                           3.0, 1.0, 0);
                        inflate_surface_and_smooth_fingers(sm_trg,
                                                           2, 1.0, 30, 1.4,
                                                           3.0, 1.0, 0);
                        
                } else if (step == 1) {
                        /* inflate and smooth surfaces */
                        inflate_surface_and_smooth_fingers(sm_src,
                                                           1, 0.2, 30, 1.0,
                                                           3.0, 1.0, 0);
                        inflate_surface_and_smooth_fingers(sm_trg,
                                                           1, 0.2, 30, 1.0,
                                                           3.0, 1.0, 0);
                        inflate_surface_and_smooth_fingers(sm_src,
                                                           2, 1.0, 10, 1.4,
                                                           3.0, 1.0, 0);
                        inflate_surface_and_smooth_fingers(sm_trg,
                                                           2, 1.0, 10, 1.4,
                                                           3.0, 1.0, 0);
                } else if (step == 2) {
                        /* inflate and smooth surfaces */
                        inflate_surface_and_smooth_fingers(sm_src,
                                                           1, 0.2, 20, 1.0,
                                                           3.0, 1.0, 0);
                        inflate_surface_and_smooth_fingers(sm_trg,
                                                           1, 0.2, 20, 1.0,
                                                           3.0, 1.0, 0);
                }
                
                get_all_polygon_point_neighbours(sm_trg, &n_neighbours, &neighbours);
                get_polygon_vertex_curvatures_cg(sm_trg, n_neighbours, neighbours,
                                         0.0, curvtype, curv_trg_tmp);

                get_all_polygon_point_neighbours(sm_src, &n_neighbours, &neighbours);
                get_polygon_vertex_curvatures_cg(sm_src, n_neighbours, neighbours,
                                         0.0, curvtype, curv_src);

                /* resample curvature of trg to src space */
                resample_noscale(sphere_trg, sphere_src, curv_trg_tmp, curv_trg);

                //smooth_heatkernel(src, curv_src, 10);
                //smooth_heatkernel(src, curv_trg, 10);
                
                /* scale values to a range of 0..1 */
//                scale_values(curv_src, src->n_points);
//                scale_values(curv_trg, src->n_points);

                output_values_any_format("curv.txt", sm_src->n_points, curv_src, TYPE_DOUBLE);
                output_values_any_format("curv_trg.txt", sm_trg->n_points, curv_trg, TYPE_DOUBLE);
                
                /* gradient of static image */
                gradient_poly(sphere_trg, dpoly_trg, curv_trg, dtheta_trg, dphi_trg);
                
                /* scale with grid size of coordinate system */
                double scale_theta = THETA, scale_phi = PHI;
                old_squared_diff = 1e15;

                double alpha = 5.0;
                int method = 3;
                /* method 1: Thirion J P 
                   Image matching as a diffusion process: an analogy with Maxwell’s demons
                   Med. Image Anal. 2 243–60, 1998 */
                /* method 2: Wang H, Dong L, O’Daniel J, Mohan R, Garden A S, Ang K K, Kuban D A, Bonnen M, Chang J Y and Cheung R
                   Validation of an accelerated ‘demons’ algorithm for deformable image registration in radiation therapy
                   Phys. Med. Biol. 50 2887–905, 2005 */
                /* method 3: Yang D, Li H, Low D A, Deasy J O and El Naqa I 
                   A fast inverse consistent deformable image registration method based on symmetric optical flow computation 
                   Phys. Med. Biol. 53 6143–65, 2008 */
                
                count_break = 0;

                for (iter_demon = 0; iter_demon < 1; iter_demon++) {

                        /* gradient of moving image */
                        if (method > 1)
                                gradient_poly(sphere_src, dpoly_src, curv_src, dtheta_src, dphi_src);

                        squared_diff = 0.0;
                        for (i = 0; i < src->n_points; i++) {
                                idiff = curv_src[i] - curv_trg[i];
                                double idiff2 = idiff*idiff;
                                squared_diff += idiff2;
                                
                                if (method == 3) {
                                        denom_trg  = (dtheta_trg[i]*dtheta_trg[i] + dphi_trg[i]*dphi_trg[i] + dtheta_src[i]*dtheta_src[i] + dphi_src[i]*dphi_src[i]) + alpha*alpha*idiff2;
                                        if (denom_trg == 0.0) {
                                                Utheta[i] = 0.0;
                                                Uphi[i] = 0.0;
                                        } else {
                                                /* passive force */
                                                Utheta[i] = idiff*dtheta_trg[i]*dtheta_src[i]/denom_trg*scale_theta;
                                                Uphi[i]   = idiff*dphi_trg[i]*dphi_src[i]/denom_trg*scale_phi;   
                                        }
                                } else if (method < 3) {
                                        /* denom for passive force */
                                        denom_trg  = ((dtheta_trg[i]*dtheta_trg[i] + dphi_trg[i]*dphi_trg[i]) + alpha*alpha*idiff2);
                                        /* denom for active force */
                                        if (method == 2)
                                                denom_src = ((dtheta_src[i]*dtheta_src[i] + dphi_src[i]*dphi_src[i]) + alpha*alpha*idiff2);
                                        else denom_src = 1.0;
                                
                                        if ((denom_trg == 0.0) || (denom_src == 0.0)) {
                                                Utheta[i] = 0.0;
                                                Uphi[i] = 0.0;
                                        } else {
                                                /* passive force */
                                                Utheta[i] = idiff*dtheta_trg[i]/denom_trg*scale_theta;
                                                Uphi[i]   = idiff*dphi_trg[i]/denom_trg*scale_phi;  
                                                 
                                                /* active force */
                                                if (method == 2) {
                                                        Utheta[i] += idiff*dtheta_src[i]/denom_src*scale_theta;
                                                        Uphi[i]   += idiff*dphi_src[i]/denom_src*scale_phi;   
                                                }
                                        }
                                } 
                        }
                        
                        fprintf(stderr,"squared diff: %g\n",squared_diff);
                
                        /* lowpass filtering of displacements */
                        smooth_heatkernel(src, Utheta, fwhm);
                        smooth_heatkernel(src, Uphi, fwhm);
                        
                        /* sum up deformations */
                        for (i = 0; i < src->n_points; i++) {
                                u[i] += Utheta[i];
                                v[i] += Uphi[i];
                        }
                        
                        apply_uv_warp(warped_sphere_src, sphere_src, u, v, 0);

                        resample_values(sphere_src, warped_sphere_src, curv_src, curv_src_tmp);
                        for (i = 0; i < src->n_points; i++) curv_src[i] = curv_src_tmp[i];
                        
                        if ((iter_demon > 4) && (old_squared_diff < squared_diff)) {
                                count_break++;
                                if (count_break > 0) break;

                        }
                        old_squared_diff = squared_diff;   
                        
                        /* use larger constrain with each step to limit max. deformations to 1/(2*alpha) */
                        alpha += 5.0;        

                }

                /* use smaller FWHM for next steps */
                //fwhm *= 0.5;
        }
        fprintf(stderr,"\n");
        
        if (debug) {
                output_values_any_format("u.txt", src->n_points, u, TYPE_DOUBLE);
                output_values_any_format("v.txt", src->n_points, v, TYPE_DOUBLE);
                output_values_any_format("warped_curv.txt", src->n_points, curv_src, TYPE_DOUBLE);
        }

        if (output_sphere_file != NULL) {
  
                /* get a pointer to the surface */
                *get_polygons_ptr(objects[0]) = *warped_sphere_src;
                if (output_graphics_any_format(output_sphere_file, format,
                                               n_objects, objects, NULL) != OK)
                        exit(EXIT_FAILURE);
        }

        delete_object_list(n_objects, objects);

        return(EXIT_SUCCESS);
}
