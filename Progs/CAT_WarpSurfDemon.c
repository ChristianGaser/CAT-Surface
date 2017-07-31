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

#define RADIANS(deg) ((PI * (double)(deg)) / 180.0)
#define DEGREES(rad) ((180.0 * (double)(rad)) / PI)

/* defaults */
char *src_file           = NULL;
char *sphere_src_file    = NULL;
char *trg_file           = NULL;
char *sphere_trg_file    = NULL;
char *output_surface_file        = NULL;
char *output_sphere_file = NULL;

int rotate       = 1;
int curvtype0    = 0;
int curvtype     = 0;
int n_steps      = 1;
int debug        = 0;
int iters        = 100;
int method       = 3;
double rate      = 1.05;
double fwhm_flow = 25.0;
double fwhm_curv = 6.0;
double alpha     = 0.7;

static ArgvInfo argTable[] = {
  {"-i", ARGV_STRING, (char *) 1, (char *) &src_file, 
     "Input file."},
  {"-is", ARGV_STRING, (char *) 1, (char *) &sphere_src_file, 
     "Input sphere file."},
  {"-t", ARGV_STRING, (char *) 1, (char *) &trg_file, 
     "Template file."},
  {"-ts", ARGV_STRING, (char *) 1, (char *) &sphere_trg_file, 
     "Template sphere file."},
  {"-w", ARGV_STRING, (char *) 1, (char *) &output_surface_file, 
     "Warped brain."},
  {"-ws", ARGV_STRING, (char *) 1, (char *) &output_sphere_file, 
     "Warped input sphere."},
  {"-fwhm-flow", ARGV_FLOAT, (char *) 1, (char *) &fwhm_flow,
     "Filter size for displacement map in FWHM."},
  {"-fwhm", ARGV_FLOAT, (char *) 1, (char *) &fwhm_curv,
     "Filter size for curvature map in FWHM."},
  {"-rate", ARGV_FLOAT, (char *) 1, (char *) &rate,
     "Change of fwhm and alpha for each iteration."},
  {"-alpha", ARGV_FLOAT, (char *) 1, (char *) &alpha,
     "ALPHA."},
  {"-maxiters", ARGV_INT, (char *) 1, (char *) &iters,
     "Maximum number of iterations."},
  {"-steps", ARGV_INT, (char *) 1, (char *) &n_steps,
     "Number of multigrid steps."},
  {"-method", ARGV_INT, (char *) 1, (char *) &method,
     "Demon method \n\t1 - Thirions approach using passive force only (Thirion JP Med. Image Anal. 2 243–60, 1998)\n\t2 - Accelerated demon using active and passive forces (Wang et al. Phys. Med. Biol. 50 2887–905, 2005)\n\t3 - Fast inverse consistent demon (Yang et al. Phys. Med. Biol. 53 6143–65, 2008)"},
  {"-norot", ARGV_CONSTANT, (char *) FALSE, (char *) &rotate,
     "Don't rotate input surface before warping."},
  {"-type0", ARGV_INT, (char *) 1, (char *) &curvtype0,
     "Curvature type\n\t0 - mean curvature (averaged over 3mm, in degrees)\n\t1 - gaussian curvature\n\t2 - curvedness\n\t3 - shape index\n\t4 - mean curvature (in radians)\n\t5 - sulcal depth like estimator."},
  {"-type", ARGV_INT, (char *) 1, (char *) &curvtype,
     "Curvature type\n\t0 - mean curvature (averaged over 3mm, in degrees)\n\t1 - gaussian curvature\n\t2 - curvedness\n\t3 - shape index\n\t4 - mean curvature (in radians)\n\t5 - sulcal depth like estimator."},
  {"-debug", ARGV_CONSTANT, (char *) TRUE, (char *) &debug,
     "Save debug files."},
   {NULL, ARGV_END, NULL, NULL, NULL}
};


/* Compose two vector fields */
void
compose_field(double ax[], double ay[], double bx[], double by[], int n_values)
{
/* This is just a simplified version (simple addition) for test purposes !!!!x*/                
        int i;
        
        for (i = 0; i < n_values; i++) {
                ax[i] += bx[i];
                ay[i] += by[i];
        }
}

/* Exponentiate vector field */
void
expfield(double vx[], double vy[], int n_values)
{
        int i, j;
        double n, mx_normv2 = -FLT_MAX;
        
        for (i = 0; i < n_values; i++)
                mx_normv2 = MAX(mx_normv2, (vx[i]*vx[i] + vy[i]*vy[i]));
                
        mx_normv2 = sqrt(mx_normv2);
        if (mx_normv2 == 0.0)
                n = 0;
        else {
                n = ceil(log2(2.0*mx_normv2)); /* n big enough so max(v * 2^-n) < 0.5 pixel) */
                n = MAX(n,0); /* avoid null values */
        }
                
printf("%g %g %g\n",n,mx_normv2,log2(2.0*mx_normv2));
        if (n != 0) {         
                for (i = 0; i < n_values; i++) {
                        vx[i] *= pow(2,-n);
                        vy[i] *= pow(2,-n);
                }
        }
    
        for (j = 0; j < n; j++)
                compose_field(vx, vy, vx, vy, n_values);

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

        /* combine x and y rotation */
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

/* This function uses the approach of the matlab function
 * SD_rotateAtlas2Sphere.m from the Spherical Demon software of Thomas Yeo
 * and Mert Sabuncu */
void
rotate_polygons_to_atlas(polygons_struct *src, polygons_struct *src_sphere,
                         polygons_struct *trg, polygons_struct *trg_sphere,
                         double fwhm, int curvtype, double *rot)
{
        int             i;
        int             n_angles;
        double          alpha, beta, gamma, sum_sq, min_sum_sq;
        double          degrees, delta, d;
        double          curr_alpha = 0.0, curr_beta = 0.0, curr_gamma = 0.0;
        double          best_alpha, best_beta, best_gamma;
        polygons_struct rot_src_sphere;
        double          rotation_tmp[9], min_degrees, max_degrees;
        double          *orig_trg, *map_trg, *map_src;
                
        min_degrees = RADIANS(1.0);
        max_degrees = RADIANS(32.0);
        degrees = max_degrees;
        n_angles  = 4;
        
        min_sum_sq = 1e15;

        orig_trg = (double *) malloc(sizeof(double) * trg->n_points);
        map_trg  = (double *) malloc(sizeof(double) * src->n_points);
        map_src  = (double *) malloc(sizeof(double) * src->n_points);

        get_smoothed_curvatures(trg, orig_trg,
                                fwhm, curvtype);
        get_smoothed_curvatures(src, map_src,
                                fwhm, curvtype);

        for (degrees = max_degrees; degrees >= min_degrees; degrees /= 2.0f) {
                delta = 2.0*degrees / (double) n_angles;
                for (alpha = curr_alpha - degrees;
                     alpha < curr_alpha + degrees; alpha += delta) {
                        for (beta = curr_beta - degrees;
                             beta < curr_beta + degrees; beta += delta) {
                                for (gamma = curr_gamma - degrees;
                                     gamma < curr_gamma + degrees;
                                     gamma += delta) {
                                
                                        /* rotate source sphere */
                                        rotation_to_matrix(rotation_tmp, alpha,
                                                           beta, gamma);
                                        rotate_polygons(src_sphere,
                                                        &rot_src_sphere,
                                                        rotation_tmp);
                                        resample_values_sphere(trg_sphere,
                                                         &rot_src_sphere,
                                                         orig_trg, map_trg, 1);

                                        /* estimate squared difference between
                                         * rotated source map and target map */
                                        sum_sq = 0.0;
                                        for (i = 0; i < src->n_points; i++) {
                                                d = map_src[i] - map_trg[i];
                                                sum_sq += d*d;
                                        }

                                        if(sum_sq < min_sum_sq) {
                                                min_sum_sq = sum_sq;
                                                best_alpha = alpha;
                                                best_beta = beta;
                                                best_gamma = gamma;
                                                rot[0] = best_alpha;
                                                rot[1] = best_beta;
                                                rot[2] = best_gamma;
                                                printf("alpha: %5.3f\tbeta: %5.3f\tgamma: %5.3f\tsquared difference: %5.3f",
                                                        DEGREES(alpha),
                                                        DEGREES(beta),
                                                        DEGREES(gamma), sum_sq);
                                                printf( "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b");
                                                printf( "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b");
                                                printf( "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b");
                                                printf( "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b");

                                        }
                                }
                        }
                }
                
                /* save best estimates */
                curr_alpha = best_alpha;
                curr_beta  = best_beta;
                curr_gamma = best_gamma;
        }
        printf("\n");
        free(orig_trg);
        free(map_trg);
        free(map_src);
}

/* qicksort */
void swap(double *a, double *b)
{
  double t=*a; *a=*b; *b=t;
}

void sort(double arr[], int beg, int end)
{
  if (end > beg + 1)
  {
    double piv = arr[beg];
    int l = beg + 1, r = end;
    while (l < r)
    {
      if (arr[l] <= piv)
        l++;
      else
        swap(&arr[l], &arr[--r]);
    }
    swap(&arr[--l], &arr[beg]);
    sort(arr, beg, l);
    sort(arr, r, end);
  }
}

double
median(double *data, int length)
{
        double *data_sort;
        int i;
    
        data_sort = (double *) malloc(sizeof(double) * length);

        for (i = 0; i < length; i++) 
               data_sort[i] = data[i];
 
        sort(data_sort, 0, length);
        return data_sort[(int)(length/2)];
}

double
stdev(double *data, int length)
{
        double avg, sum2, data0;
        int i;
    
        avg = sum2 = 0.0;
        for (i = 0; i < length; i++) 
               avg  += data[i];
               
        avg /= (double)length;
                
        for (i = 0; i < length; i++) {
               data0 = data[i] - avg; 
               sum2  += data0*data0;
        }

        return sqrt(sum2);
}

void
normalizeVector(double *data, int length)
{
        double median_data, stdev_data;
        int i;
    
        median_data = median(data, length);

        for (i = 0; i < length; i++) 
                data[i] -= median_data;
                
        stdev_data = stdev(data, length);

        for (i = 0; i < length; i++)
                data[i] /= stdev_data;

        for (i = 0; i < length; i++) {
                if (data[i] < -3.0)
                        data[i] = -3.0 - (1.0 - exp(3.0 - fabs(data[i])));
                if (data[i] >  3.0)
                        data[i] =  3.0 + (1.0 - exp(3.0 - fabs(data[i])));
        }
}

double
Correlation(double *x, double *y, int N)
{
    double EX, EY, EXY, EX2, EY2;
    int i;
    
    EX = EY = EXY = EX2 = EY2 = 0.0;
    for (i = 0; i < N; i++)
    {
        EX += x[i];
        EY += y[i];
        EXY += x[i]*y[i];
        EX2 += x[i]*x[i];
        EY2 += y[i]*y[i];
    }

    return (N*EXY - EX*EY) / sqrt((N*EX2 - EX*EX) * (N*EY2 - EY*EY));
}

void
WarpDemon(polygons_struct *src, polygons_struct *sphere_src, polygons_struct *trg, 
                         polygons_struct *sphere_trg, polygons_struct *warped_sphere_src, 
                         struct dartel_poly *dpoly_src, struct dartel_poly *dpoly_trg, int type)
{
        int              *n_neighbours, **neighbours;
        int              i,it, count_break;
        double           idiff, denom_trg, denom_src, sum_diff2;
        double           *curv_trg0, *curv_trg, *curv_src0, *curv_src, cc, old_cc, distance;
        double           *dtheta_trg, *dphi_trg, *dtheta_src, *dphi_src, *u, *v, *Utheta, *Uphi;
        polygons_struct  *warped_sphere_trg;
        
        if (src->n_points != trg->n_points) {
                fprintf(stderr,"Source and target have different size!");
                return; 
        }
        
        curv_src0         = (double *) malloc(sizeof(double) * src->n_points);
        curv_src          = (double *) malloc(sizeof(double) * src->n_points);
        curv_trg0         = (double *) malloc(sizeof(double) * src->n_points);
        curv_trg          = (double *) malloc(sizeof(double) * src->n_points);
        dtheta_trg        = (double *) malloc(sizeof(double) * src->n_points);
        dphi_trg          = (double *) malloc(sizeof(double) * src->n_points);
        dtheta_src        = (double *) malloc(sizeof(double) * src->n_points);
        dphi_src          = (double *) malloc(sizeof(double) * src->n_points);
        u                 = (double *) malloc(sizeof(double) * src->n_points);
        v                 = (double *) malloc(sizeof(double) * src->n_points);
        Utheta            = (double *) malloc(sizeof(double) * src->n_points);
        Uphi              = (double *) malloc(sizeof(double) * src->n_points);

        warped_sphere_trg = (polygons_struct *) malloc(sizeof(polygons_struct));

        if (type == 0)
                distance = 3.0;
        else distance = 0.0;

        get_all_polygon_point_neighbours(trg, &n_neighbours, &neighbours);
        get_polygon_vertex_curvatures_cg(trg, n_neighbours, neighbours,
                                         distance, type, curv_trg0);

        get_all_polygon_point_neighbours(src, &n_neighbours, &neighbours);
        get_polygon_vertex_curvatures_cg(src, n_neighbours, neighbours,
                                         distance, type, curv_src0);

        normalizeVector(curv_trg0, src->n_points);
        normalizeVector(curv_src0, src->n_points);
                        
        for (i = 0; i < src->n_points; i++)  {
                curv_src[i] = curv_src0[i];
                curv_trg[i] = curv_trg0[i];
                u[i] = 0.0;
                v[i] = 0.0;
        }

        if (debug) {
                output_values_any_format("curv.txt", src->n_points, curv_src, TYPE_DOUBLE);
                output_values_any_format("curv_trg.txt", src->n_points, curv_trg, TYPE_DOUBLE);
        }
        
        /* gradient of static image */
        gradient_poly(sphere_trg, dpoly_trg, curv_trg, dtheta_trg, dphi_trg);
                
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
        old_cc = -FLT_MAX;
                
        copy_polygons(sphere_src, warped_sphere_src);
        copy_polygons(sphere_trg, warped_sphere_trg);

        for (it = 0; it < iters; it++) {

                /* gradient of moving image */
                if (method > 1)
                        gradient_poly(sphere_src, dpoly_src, curv_src, dtheta_src, dphi_src);

                if (method == 3)
                        gradient_poly(sphere_trg, dpoly_trg, curv_trg, dtheta_trg, dphi_trg);
                        
                sum_diff2 = 0.0;
                for (i = 0; i < src->n_points; i++) {
                        idiff = curv_src[i] - curv_trg[i];
                        double idiff2 = idiff*idiff;      
                        sum_diff2 += idiff2;                          
                                
                        if (method == 3) {
                                denom_trg  = (dtheta_trg[i]*dtheta_trg[i] + dphi_trg[i]*dphi_trg[i] + 
                                              dtheta_src[i]*dtheta_src[i] + dphi_src[i]*dphi_src[i]) + alpha*alpha*idiff2;
                                if (denom_trg == 0.0) {
                                        Utheta[i] = 0.0;
                                        Uphi[i] = 0.0;
                                } else {
                                        /* inverse consistent force */
                                        Utheta[i] = idiff*(dtheta_trg[i]+dtheta_src[i])/denom_trg;
                                        Uphi[i]   = idiff*(dphi_trg[i]+dphi_src[i])/denom_trg;   
                                }
                        } else {
                                /* denom for passive force */
                                denom_trg = ((dtheta_trg[i]*dtheta_trg[i] + dphi_trg[i]*dphi_trg[i]) + alpha*alpha*idiff2);
                                /* denom for active force */
                                if (method == 2)
                                        denom_src  = ((dtheta_src[i]*dtheta_src[i] + dphi_src[i]*dphi_src[i]) + alpha*alpha*idiff2);
                                else denom_src = 1.0;
                                
                                if ((denom_trg == 0.0) || (denom_src == 0.0)) {
                                        Utheta[i] = 0.0;
                                        Uphi[i] = 0.0;
                                } else {
                                        /* passive force */
                                        Utheta[i] = idiff*dtheta_trg[i]/denom_trg;
                                        Uphi[i]   = idiff*dphi_trg[i]/denom_trg;  
                                                 
                                        /* active force */
                                        if (method == 2) {
                                                Utheta[i] += idiff*dtheta_src[i]/denom_src;
                                                Uphi[i]   += idiff*dphi_src[i]/denom_src;   
                                        }
                                }
                        } 
                }
                            
                /* lowpass filtering of displacements */
                smooth_heatkernel(src, Utheta, fwhm_flow);
                smooth_heatkernel(src, Uphi, fwhm_flow);
                        
                /* sum up deformations */
                for (i = 0; i < src->n_points; i++) {
                        Utheta[i] *= THETA;
                        Uphi[i]   *= PHI;
                        u[i]      += Utheta[i];
                        v[i]      += Uphi[i];
                }
                        
                apply_uv_warp(warped_sphere_src, warped_sphere_src, Utheta, Uphi, 1);
                resample_values_sphere(sphere_src, warped_sphere_src, curv_src0, curv_src, 1);                        
                        
                normalizeVector(curv_src, src->n_points);

                if (method == 3) {
                        apply_uv_warp(warped_sphere_trg, warped_sphere_trg, Utheta, Uphi, 0);
                        resample_values_sphere(sphere_trg, warped_sphere_trg, curv_trg0, curv_trg, 1);                                   
                        normalizeVector(curv_trg, src->n_points);
                }

                cc = Correlation(curv_src, curv_trg, src->n_points);
//                fprintf(stderr, "%02d: CC=%5.4f diff=%g\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b", it+1, cc, sum_diff2);
                printf( "%02d: CC=%5.4f diff=%g alpha=%g fwhm-flow=%g\n", it+1, cc, sum_diff2,alpha,fwhm_flow);

                /* use at least 5 iterations */
                if ((it >= 5) && (cc/old_cc < 1.0025)) {
                        count_break++;
                        if (count_break > 0) break;

                }
                old_cc = cc;   
                alpha *= rate;
                fwhm_flow *= rate;

        }
        
        if (debug) {
                output_values_any_format("u.txt", src->n_points, u, TYPE_DOUBLE);
                output_values_any_format("v.txt", src->n_points, v, TYPE_DOUBLE);
                output_values_any_format("warped_curv.txt", src->n_points, curv_src, TYPE_DOUBLE);
        }

        /* invert deformation because we need inverse transformation */
        apply_uv_warp(sphere_src, warped_sphere_src, u, v, 0);

        free(curv_src0);
        free(curv_src);
        free(curv_trg0);
        free(curv_trg);
        free(dtheta_trg);
        free(dphi_trg);
        free(dtheta_src);
        free(dphi_src);
        free(u);
        free(v);
        free(Utheta);
        free(Uphi);
}

int
main(int argc, char *argv[])
{
        File_formats     format;
        polygons_struct  *src, *trg, *sphere_src, *sphere_trg;
        polygons_struct  *sm_src, *sm_trg, *sm_sphere_src, *sm_sphere_trg, *warped_sphere_src;
        int              i, step;
        int              n_objects;
        object_struct    **objects;
        double           rotation_matrix[9], rot[3];
        struct           dartel_poly *dpoly_src, *dpoly_trg;
        int              *n_neighbours, **neighbours;

        /* get the arguments from the command line */

        if (ParseArgv(&argc, argv, argTable, 0) ||
            src_file == NULL || trg_file == NULL || 
            sphere_src_file == NULL || sphere_trg_file == NULL ||
            (output_surface_file == NULL && output_sphere_file == NULL)) {
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

        dpoly_src         = (struct dartel_poly *) malloc(sizeof(struct dartel_poly));
        dpoly_trg         = (struct dartel_poly *) malloc(sizeof(struct dartel_poly));

        sm_src            = (polygons_struct *) malloc(sizeof(polygons_struct));
        sm_trg            = (polygons_struct *) malloc(sizeof(polygons_struct));
        sm_sphere_src     = (polygons_struct *) malloc(sizeof(polygons_struct));
        sm_sphere_trg     = (polygons_struct *) malloc(sizeof(polygons_struct));
        warped_sphere_src = (polygons_struct *) malloc(sizeof(polygons_struct));
        
        int n_points = 81920;
//        n_points = 20480;

        if (n_steps > 1)
                printf("Warning: Multistep approach is not yet working!\n");

        for (step = 0; step < n_steps; step++) {

                objects = resample_surface( src,  sphere_src, n_points,  NULL,  NULL);
                sm_src = get_polygons_ptr(objects[0]);
                objects = resample_surface( sphere_src,  sphere_src, n_points,  NULL,  NULL);
                sm_sphere_src = get_polygons_ptr(objects[0]);
                objects = resample_surface( trg,  sphere_trg, n_points,  NULL,  NULL);
                sm_trg = get_polygons_ptr(objects[0]);
                objects = resample_surface( sphere_trg,  sphere_trg, n_points,  NULL,  NULL);
                sm_sphere_trg = get_polygons_ptr(objects[0]);
                
                init_dartel_poly(sm_sphere_src, dpoly_src);
                init_dartel_poly(sm_sphere_trg, dpoly_trg);
                
                printf("Resample surfaces to %d points\n",sm_src->n_points);

                /* initialization */
                if (step == 0) {
                        smooth_heatkernel(sm_src, NULL, fwhm_curv);
                        smooth_heatkernel(sm_trg, NULL, fwhm_curv);
                        
                        /* initial rotation */
                        if (rotate) {
                                rotate_polygons_to_atlas(sm_src, sm_sphere_src,
                                                         sm_trg, sm_sphere_trg,
                                                         10.0, 2, rot);

                                rotation_to_matrix(rotation_matrix,
                                                   rot[0], rot[1], rot[2]);
                
                                /* rotate source sphere */
                                rotate_polygons(sm_sphere_src,
                                                NULL, rotation_matrix);
                        }



                } else if (step == 1) {
                        smooth_heatkernel(sm_src, NULL, fwhm_curv*rate);
                        smooth_heatkernel(sm_trg, NULL, fwhm_curv*rate);

                        curvtype0 = curvtype;
                } else if (step == 2) {
                        smooth_heatkernel(sm_src, NULL, fwhm_curv*rate*rate);
                        smooth_heatkernel(sm_trg, NULL, fwhm_curv*rate*rate);

                        /* use default at final step */
                        curvtype0 = curvtype;
                }
                                
                WarpDemon(sm_src, sm_sphere_src, sm_trg, 
                         sm_sphere_trg, warped_sphere_src, dpoly_src, dpoly_trg, curvtype0);

                objects = resample_surface_to_target_sphere(sm_sphere_src, warped_sphere_src, sphere_src, NULL, NULL);
                sphere_src = get_polygons_ptr(objects[0]);

                /* use smaller FWHM for next steps */
                fwhm_flow *= rate;
        }
        printf("\n");
        
        if (output_sphere_file != NULL) {
  
                if (output_graphics_any_format(output_sphere_file, format,
                                               n_objects, objects, NULL) != OK)
                        exit(EXIT_FAILURE);
        }

        return(EXIT_SUCCESS);
}
