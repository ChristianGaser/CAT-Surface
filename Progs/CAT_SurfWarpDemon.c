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

#include "CAT_Warp.h"
#include "CAT_Map.h"
#include "CAT_Vol.h"
#include "CAT_Surf.h"
#include "CAT_Curvature.h"
#include "CAT_SurfaceIO.h"
#include "CAT_Interpolate.h"
#include "CAT_Resample.h"
#include "dartel.h"

#define RADIANS(deg) ((PI * (double)(deg)) / 180.0)
#define DEGREES(rad) ((180.0 * (double)(rad)) / PI)

/* defaults */
char *src_file       = NULL;
char *src_sphere_file  = NULL;
char *trg_file       = NULL;
char *trg_sphere_file  = NULL;
char *output_surface_file    = NULL;
char *output_sphere_file = NULL;

int n_points = 20480;
int rotate     = 1;
int curvtype0  = 5;
int curvtype   = 5;
int n_steps    = 1;
int debug    = 0;
int iters    = 100;
int method     = 3;
int verbose    = 0;
double rate    = 0.97;  /* rate < 1: fwhm_flow DECREASES (coarse-to-fine like multi-resolution) */
double fwhm_flow = 30.0;  /* Start with broader smoothing, then refine */
double fwhm_curv = 6.0;
double fwhm_disp = 5.0;  /* FWHM for displacement field smoothing (elastic prior) */
double alpha0  = 0.5;  /* Reduced alpha for smoother registration */
double max_step_deg = 10.0;  /* Reduced max step for stability */
double sigma_x_default = 1.0;  /* regularization weight for Spherical Demons */
int smooth_velocity = 1;  /* smooth velocity UPDATE (fluid prior) - NOW ON by default */
int smooth_displacement = 1;  /* smooth accumulated DISPLACEMENT (elastic prior) */
int use_hessian = 1;  /* use per-vertex Hessian-based update (Newton-Raphson) */
int use_line_search = 1;  /* use line search during integration */
double step_factor = 1.0;  /* global step size factor */

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
   "Number of points for resampling (e.g. 81920 or 20480)."},
  {"-fwhm-flow", ARGV_FLOAT, (char *) 1, (char *) &fwhm_flow,
   "Filter size for displacement map in FWHM."},
  {"-fwhm", ARGV_FLOAT, (char *) 1, (char *) &fwhm_curv,
   "Filter size for curvature map in FWHM."},
  {"-rate", ARGV_FLOAT, (char *) 1, (char *) &rate,
   "Change of fwhm and alpha for each iteration."},
  {"-alpha", ARGV_FLOAT, (char *) 1, (char *) &alpha0,
   "ALPHA."},
  {"-max-step-deg", ARGV_FLOAT, (char *) 1, (char *) &max_step_deg,
   "Clamp per-iteration |dtheta,dphi| to this many degrees (<=0 disables)."},
  {"-fwhm-disp", ARGV_FLOAT, (char *) 1, (char *) &fwhm_disp,
   "Filter size for displacement field smoothing (elastic prior) in FWHM."},
  {"-no-smooth-velocity", ARGV_CONSTANT, (char *) FALSE, (char *) &smooth_velocity,
   "Disable velocity UPDATE smoothing (fluid prior, default ON)."},
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
  {"-maxiters", ARGV_INT, (char *) 1, (char *) &iters,
   "Maximum number of iterations."},
  {"-steps", ARGV_INT, (char *) 1, (char *) &n_steps,
   "Number of multigrid steps."},
  {"-method", ARGV_INT, (char *) 1, (char *) &method,
   "Demon method \n\t1 - Thirions approach using passive force only (Thirion JP Med. Image Anal. 2 243–60, 1998)\n\t2 - Accelerated demon using active and passive forces (Wang et al. Phys. Med. Biol. 50 2887–905, 2005)\n\t3 - Fast inverse consistent demon (Yang et al. Phys. Med. Biol. 53 6143–65, 2008)\n\t4 - Spherical Demons (Yeo et al. IEEE TMI 29 469-486, 2010)"},
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

/**
 * Project a 3D tangent vector onto the tangent plane of the sphere at a given vertex.
 * This ensures that the update vector lies in the tangent plane of the sphere.
 */
void
project_onto_tangent_plane(polygons_struct *sphere, double *vec_x, double *vec_y, double *vec_z)
{
    int i;
    double nx, ny, nz, dot, len;
    
    for (i = 0; i < sphere->n_points; i++) {
        /* Get vertex position (which is also the normal for a unit sphere) */
        nx = Point_x(sphere->points[i]);
        ny = Point_y(sphere->points[i]);
        nz = Point_z(sphere->points[i]);
        
        /* Normalize the normal vector */
        len = sqrt(nx*nx + ny*ny + nz*nz);
        if (len > 1e-10) {
            nx /= len;
            ny /= len;
            nz /= len;
        }
        
        /* Project: v_tangent = v - (v · n) * n */
        dot = vec_x[i]*nx + vec_y[i]*ny + vec_z[i]*nz;
        vec_x[i] -= dot * nx;
        vec_y[i] -= dot * ny;
        vec_z[i] -= dot * nz;
    }
}

/**
 * Normalize sphere vertices to have a constant radius.
 * This ensures vertices stay on the sphere after warping.
 */
void
normalize_sphere_radius(polygons_struct *sphere, double radius)
{
    int i;
    double len;
    
    for (i = 0; i < sphere->n_points; i++) {
        len = sqrt(Point_x(sphere->points[i]) * Point_x(sphere->points[i]) +
                   Point_y(sphere->points[i]) * Point_y(sphere->points[i]) +
                   Point_z(sphere->points[i]) * Point_z(sphere->points[i]));
        if (len > 1e-10) {
            double scale = radius / len;
            Point_x(sphere->points[i]) *= scale;
            Point_y(sphere->points[i]) *= scale;
            Point_z(sphere->points[i]) *= scale;
        }
    }
}

/**
 * Compute orthonormal tangent plane basis vectors (e1, e2) for each vertex.
 * These are used to project 3D vectors onto the tangent plane.
 * e1 is in the theta direction, e2 is in the phi direction.
 */
void
compute_tangent_basis(polygons_struct *sphere, double **e1x, double **e1y, double **e1z,
                      double **e2x, double **e2y, double **e2z)
{
    int i;
    int n = sphere->n_points;
    
    *e1x = (double *)malloc(sizeof(double) * n);
    *e1y = (double *)malloc(sizeof(double) * n);
    *e1z = (double *)malloc(sizeof(double) * n);
    *e2x = (double *)malloc(sizeof(double) * n);
    *e2y = (double *)malloc(sizeof(double) * n);
    *e2z = (double *)malloc(sizeof(double) * n);
    
    for (i = 0; i < n; i++) {
        double x = Point_x(sphere->points[i]);
        double y = Point_y(sphere->points[i]);
        double z = Point_z(sphere->points[i]);
        double r = sqrt(x*x + y*y + z*z);
        double rxy = sqrt(x*x + y*y);
        
        if (r < 1e-10) r = 1e-10;
        
        /* e1 is in the theta direction (longitude) */
        /* For spherical coords: d/dtheta = (-sin(theta)*sin(phi), cos(theta)*sin(phi), 0) */
        /* Simplified: tangent along longitude circles */
        if (rxy > 1e-10) {
            (*e1x)[i] = -y / rxy;
            (*e1y)[i] = x / rxy;
            (*e1z)[i] = 0.0;
        } else {
            /* At poles, pick arbitrary tangent direction */
            (*e1x)[i] = 1.0;
            (*e1y)[i] = 0.0;
            (*e1z)[i] = 0.0;
        }
        
        /* e2 is in the phi direction (latitude), orthogonal to e1 and normal */
        /* e2 = n x e1 where n = (x,y,z)/r */
        double nx = x / r, ny = y / r, nz = z / r;
        (*e2x)[i] = ny * (*e1z)[i] - nz * (*e1y)[i];
        (*e2y)[i] = nz * (*e1x)[i] - nx * (*e1z)[i];
        (*e2z)[i] = nx * (*e1y)[i] - ny * (*e1x)[i];
        
        /* Normalize e2 */
        double e2len = sqrt((*e2x)[i]*(*e2x)[i] + (*e2y)[i]*(*e2y)[i] + (*e2z)[i]*(*e2z)[i]);
        if (e2len > 1e-10) {
            (*e2x)[i] /= e2len;
            (*e2y)[i] /= e2len;
            (*e2z)[i] /= e2len;
        }
    }
}

/**
 * Compute 3D Cartesian gradient of a scalar field on the mesh using neighbor differences.
 * Returns gradient projected onto tangent plane as (grad_e1, grad_e2).
 */
void
compute_mesh_gradient_3d(polygons_struct *polygons, double *f,
                         int **neighbours, int *n_neighbours,
                         double *e1x, double *e1y, double *e1z,
                         double *e2x, double *e2y, double *e2z,
                         double *grad_e1, double *grad_e2)
{
    int i, j, n;
    double gx, gy, gz, total_weight;
    
    for (i = 0; i < polygons->n_points; i++) {
        gx = gy = gz = 0.0;
        total_weight = 0.0;
        
        double xi = Point_x(polygons->points[i]);
        double yi = Point_y(polygons->points[i]);
        double zi = Point_z(polygons->points[i]);
        
        /* Compute weighted gradient from neighbors */
        n = n_neighbours[i];
        for (j = 0; j < n; j++) {
            int nb = neighbours[i][j];
            double xj = Point_x(polygons->points[nb]);
            double yj = Point_y(polygons->points[nb]);
            double zj = Point_z(polygons->points[nb]);
            
            double dx = xj - xi;
            double dy = yj - yi;
            double dz = zj - zi;
            double dist = sqrt(dx*dx + dy*dy + dz*dz);
            
            if (dist > 1e-10) {
                double df = f[nb] - f[i];
                double weight = 1.0 / dist;
                gx += weight * df * dx / dist;
                gy += weight * df * dy / dist;
                gz += weight * df * dz / dist;
                total_weight += weight;
            }
        }
        
        if (total_weight > 1e-10) {
            gx /= total_weight;
            gy /= total_weight;
            gz /= total_weight;
        }
        
        /* Project onto tangent plane basis */
        grad_e1[i] = gx * e1x[i] + gy * e1y[i] + gz * e1z[i];
        grad_e2[i] = gx * e2x[i] + gy * e2y[i] + gz * e2z[i];
    }
}

/**
 * Compute per-vertex 2x2 Hessian inverse and apply to gradient.
 * This implements the Newton-Raphson update from Spherical Demons.
 * H = G*G^T + sigma_x^-2 * I, update = H^-1 * gradient
 */
void
apply_hessian_update(double *g1, double *g2, double *diff, int n,
                     double sigma_x, double min_step_sq,
                     double *u1, double *u2)
{
    int i;
    double reg = 1.0 / (sigma_x * sigma_x * min_step_sq + 1e-10);
    
    for (i = 0; i < n; i++) {
        /* Form 2x2 Hessian: H = [g1*g1 + reg, g1*g2; g1*g2, g2*g2 + reg] */
        double h11 = g1[i] * g1[i] + reg;
        double h12 = g1[i] * g2[i];
        double h22 = g2[i] * g2[i] + reg;
        
        /* Compute determinant */
        double det = h11 * h22 - h12 * h12;
        
        /* Residual: diff * gradient */
        double r1 = diff[i] * g1[i];
        double r2 = diff[i] * g2[i];
        
        if (fabs(det) > 1e-10) {
            /* H^-1 * residual */
            u1[i] = (h22 * r1 - h12 * r2) / det;
            u2[i] = (-h12 * r1 + h11 * r2) / det;
        } else {
            /* Fall back to gradient descent */
            u1[i] = r1 / (reg + 1e-10);
            u2[i] = r2 / (reg + 1e-10);
        }
    }
}

/**
 * Spherical exponential map using scaling and squaring.
 * This implements diffeomorphic integration of the velocity field.
 * The velocity field is scaled down, then composed with itself repeatedly.
 */
void
spherical_exp_map(polygons_struct *sphere, polygons_struct *src_sphere,
                  double *u, double *v, int n_points, double min_dist, int num_compositions)
{
    int i, j, n_comp;
    double max_vec_norm = 0.0;
    double *u_scaled, *v_scaled;
    
    /* Find maximum vector magnitude */
    for (i = 0; i < n_points; i++) {
        double norm = sqrt(u[i]*u[i] + v[i]*v[i]);
        if (norm > max_vec_norm) max_vec_norm = norm;
    }
    
    /* Determine number of compositions needed for diffeomorphic integration */
    /* n big enough so max(v * 2^-n) < 0.5 * min_dist */
    if (num_compositions <= 0) {
        if (max_vec_norm < 1e-10 || min_dist < 1e-10) {
            n_comp = 0;
        } else {
            n_comp = (int)ceil(log2(2.0 * max_vec_norm / min_dist)) + 3;
            n_comp = MAX(n_comp, 0);
        }
    } else {
        n_comp = num_compositions;
    }
    
    if (n_comp == 0) {
        /* Direct application */
        apply_uv_warp(sphere, src_sphere, u, v, 1);
        return;
    }
    
    /* Allocate scaled velocity field */
    u_scaled = (double *)malloc(sizeof(double) * n_points);
    v_scaled = (double *)malloc(sizeof(double) * n_points);
    
    /* Scale down the velocity field */
    double scale = pow(2.0, -n_comp);
    for (i = 0; i < n_points; i++) {
        u_scaled[i] = u[i] * scale;
        v_scaled[i] = v[i] * scale;
    }
    
    /* Apply first small step */
    copy_polygons(src_sphere, sphere);
    apply_uv_warp(sphere, sphere, u_scaled, v_scaled, 1);
    normalize_sphere_radius(sphere, 100.0);
    
    /* Compose with itself n_comp times (scaling and squaring) */
    for (j = 0; j < n_comp; j++) {
        /* Resample the current deformation and compose */
        apply_uv_warp(sphere, sphere, u_scaled, v_scaled, 1);
        normalize_sphere_radius(sphere, 100.0);
    }
    
    free(u_scaled);
    free(v_scaled);
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
normalizeVector(double *data, int length)
{
    double median_data, stdev_data;
    int i;
  
    median_data = get_median_double(data, length, 0);

    for (i = 0; i < length; i++) 
        data[i] -= median_data;
        
    stdev_data = get_std_double(data, length, 0);

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
WarpDemon(polygons_struct *src, polygons_struct *src_sphere, polygons_struct *orig_sphere,
             polygons_struct *trg, polygons_struct *trg_sphere, polygons_struct *warped_src_sphere, 
             struct dartel_poly *dpoly_src, struct dartel_poly *dpoly_trg, int type)
{
    int        *n_neighbours, **neighbours;
    int        i,it, count_break;
    double       idiff, denom_trg, denom_src, sum_diff2;
    double       *curv_trg0, *curv_trg, *curv_src0, *curv_src, cc, old_cc, distance;
    double       *dtheta_trg, *dphi_trg, *dtheta_src, *dphi_src, *u, *v, *Utheta, *Uphi;
    polygons_struct  *warped_trg_sphere;
    
    if (src->n_points != trg->n_points) {
        fprintf(stderr,"Source and target have different size!");
        return; 
    }
    
    curv_src0     = (double *) malloc(sizeof(double) * src->n_points);
    curv_src      = (double *) malloc(sizeof(double) * src->n_points);
    curv_trg0     = (double *) malloc(sizeof(double) * src->n_points);
    curv_trg      = (double *) malloc(sizeof(double) * src->n_points);
    dtheta_trg      = (double *) malloc(sizeof(double) * src->n_points);
    dphi_trg      = (double *) malloc(sizeof(double) * src->n_points);
    dtheta_src      = (double *) malloc(sizeof(double) * src->n_points);
    dphi_src      = (double *) malloc(sizeof(double) * src->n_points);
    u         = (double *) malloc(sizeof(double) * src->n_points);
    v         = (double *) malloc(sizeof(double) * src->n_points);
    Utheta        = (double *) malloc(sizeof(double) * src->n_points);
    Uphi        = (double *) malloc(sizeof(double) * src->n_points);

    warped_trg_sphere = (polygons_struct *) malloc(sizeof(polygons_struct));

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
        curv_trg[i] = curv_trg0[i];
        u[i] = 0.0;
        v[i] = 0.0;
    }

    /* If src_sphere is different from orig_sphere (i.e., already warped from previous step),
     * we need to resample the curvatures through the existing warp first */
    if (src_sphere != orig_sphere) {
        /* src_sphere already contains the warp from previous step */
        resample_values_sphere(orig_sphere, src_sphere, curv_src0, curv_src, 0, 0);
        normalizeVector(curv_src, src->n_points);
    } else {
        /* First step: just copy curvatures */
        for (i = 0; i < src->n_points; i++) {
            curv_src[i] = curv_src0[i];
        }
    }

    if (debug) {
        output_values_any_format("curv.txt", src->n_points, curv_src, TYPE_DOUBLE);
        output_values_any_format("curv_trg.txt", src->n_points, curv_trg, TYPE_DOUBLE);
    }
    
    /* gradient of static image */
    gradient_poly(trg_sphere, dpoly_trg, curv_trg, dtheta_trg, dphi_trg);
        
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
        
    copy_polygons(src_sphere, warped_src_sphere);
    copy_polygons(trg_sphere, warped_trg_sphere);

    for (it = 0; it < iters; it++) {

        /* gradient of moving image */
        if (method > 1) {
             gradient_poly(src_sphere, dpoly_src, curv_src, dtheta_src, dphi_src);
        }

        if (method == 3)
            gradient_poly(trg_sphere, dpoly_trg, curv_trg, dtheta_trg, dphi_trg);
            
        sum_diff2 = 0.0;
        for (i = 0; i < src->n_points; i++) {
            idiff = curv_src[i] - curv_trg[i];
            double idiff2 = idiff*idiff;    
            sum_diff2 += idiff2;              
                
            if (method == 3) {
                denom_trg  = (dtheta_trg[i]*dtheta_trg[i] + dphi_trg[i]*dphi_trg[i] + 
                        dtheta_src[i]*dtheta_src[i] + dphi_src[i]*dphi_src[i]) + alpha0*alpha0*idiff2;
                if (denom_trg == 0.0) {
                    Utheta[i] = 0.0;
                    Uphi[i] = 0.0;
                } else {
                    /* inverse consistent force */
                    Utheta[i] = idiff*(dtheta_trg[i]+dtheta_src[i])/denom_trg;
                    Uphi[i]   = idiff*(dphi_trg[i]+dphi_src[i])/denom_trg;   
                }
            } else if (method == 4) {
                 /* Spherical Demons (Yeo et al. IEEE TMI 2010) */
                 /* Per-vertex Newton-Raphson update with Hessian approximation */
                 double g1 = dtheta_trg[i] + dtheta_src[i];  /* symmetric gradient e1 component */
                 double g2 = dphi_trg[i] + dphi_src[i];      /* symmetric gradient e2 component */
                 double G2 = g1*g1 + g2*g2;
                 
                 if (use_hessian && G2 > 1e-9) {
                     /* Form per-vertex 2x2 Hessian: H = G*G^T + reg*I */
                     /* Regularization controls step size - smaller reg = larger steps */
                     double sigma_x_sq = sigma_x_default * sigma_x_default;
                     double reg = idiff2 / sigma_x_sq + 0.1;  /* Similar form to other methods */
                     
                     double h11 = g1 * g1 + reg;
                     double h12 = g1 * g2;
                     double h22 = g2 * g2 + reg;
                     double det = h11 * h22 - h12 * h12;
                     
                     /* Residual: idiff * G */
                     double r1 = idiff * g1;
                     double r2 = idiff * g2;
                     
                     if (fabs(det) > 1e-12) {
                         /* u = H^-1 * residual */
                         Utheta[i] = (h22 * r1 - h12 * r2) / det;
                         Uphi[i]   = (-h12 * r1 + h11 * r2) / det;
                     } else {
                         /* Fallback to simple gradient */
                         double denom = G2 + reg;
                         Utheta[i] = idiff * g1 / denom;
                         Uphi[i]   = idiff * g2 / denom;
                     }
                 } else {
                     /* Simple gradient-based update */
                     double sigma_x_sq = sigma_x_default * sigma_x_default;
                     double denom = G2 + idiff2/sigma_x_sq + 0.01;
                     
                     Utheta[i] = idiff * g1 / denom;
                     Uphi[i]   = idiff * g2 / denom;
                 }
            } else {
                /* denom for passive force */
                denom_trg = ((dtheta_trg[i]*dtheta_trg[i] + dphi_trg[i]*dphi_trg[i]) + alpha0*alpha0*idiff2);
                /* denom for active force */
                if (method == 2)
                    denom_src  = ((dtheta_src[i]*dtheta_src[i] + dphi_src[i]*dphi_src[i]) + alpha0*alpha0*idiff2);
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
              
        /* lowpass filtering of velocity UPDATE (fluid prior) */
        if (smooth_velocity) {
            smooth_heatkernel(src, Utheta, fwhm_flow);
            smooth_heatkernel(src, Uphi, fwhm_flow);
        }

        /* clip per-vertex step to avoid overshoot/folding */
        if (max_step_deg > 0.0) {
            double max_step = max_step_deg * (PI / 180.0);
            for (i = 0; i < src->n_points; i++) {
                double step = sqrt(Utheta[i]*Utheta[i] + Uphi[i]*Uphi[i]);
                if (step > max_step && step > 0.0) {
                    double scale = max_step / step;
                    Utheta[i] *= scale;
                    Uphi[i]   *= scale;
                }
            }
        }
            
        /* sum up deformations with step factor */
        for (i = 0; i < src->n_points; i++) {
            Utheta[i] *= THETA * step_factor;
            Uphi[i]   *= PHI * step_factor;
            u[i]    += Utheta[i];
            v[i]    += Uphi[i];
        }
        
        /* Smooth the accumulated DISPLACEMENT field (elastic prior) */
        /* This is the key difference from the original SD - smooth AFTER integration */
        if (smooth_displacement && fwhm_disp > 0.0) {
            smooth_heatkernel(src, u, fwhm_disp);
            smooth_heatkernel(src, v, fwhm_disp);
        }
        
        /* Apply total accumulated warp using diffeomorphic exponential map */
        copy_polygons(src_sphere, warped_src_sphere);
        apply_uv_warp(warped_src_sphere, warped_src_sphere, u, v, 1);
        
        /* Normalize radius to ensure vertices stay on sphere */
        normalize_sphere_radius(warped_src_sphere, 100.0);
        
        /* Resample curvatures through the accumulated warp:
         * Use orig_sphere as source to handle multi-step properly */
        resample_values_sphere(orig_sphere, warped_src_sphere, curv_src0, curv_src, 0, 0);              
            
        normalizeVector(curv_src, src->n_points);

        if (method == 3) {
            copy_polygons(trg_sphere, warped_trg_sphere);
            apply_uv_warp(warped_trg_sphere, warped_trg_sphere, u, v, 0);
            normalize_sphere_radius(warped_trg_sphere, 100.0);
            resample_values_sphere(trg_sphere, warped_trg_sphere, curv_trg0, curv_trg, 0, 0);                   
            normalizeVector(curv_trg, src->n_points);
        }

        cc = Correlation(curv_src, curv_trg, src->n_points);
        printf("%02d: CC=%5.4f diff=%g sigma_x=%g fwhm-flow=%g step=%g\n", 
               it+1, cc, sum_diff2, sigma_x_default, fwhm_flow, step_factor);

        /* Adaptive step size: reduce if no improvement */
        if (it >= 3) {
            if (cc < old_cc) {
                /* Step was too large, backtrack */
                count_break++;
                if (use_line_search && step_factor > 0.25) {
                    step_factor *= 0.5;
                    printf("  -> Reducing step factor to %g\n", step_factor);
                }
            } else if (cc/old_cc < 1.00025) {
                count_break++;
            } else {
                count_break = 0;
                /* Gradually increase step if making progress */
                if (step_factor < 1.0) step_factor *= 1.1;
            }
            
            if (count_break > 2) break;
        }
        
        old_cc = cc;
        /* Coarse-to-fine strategy: rate < 1 decreases fwhm_flow each iteration
         * (mimics multi-resolution approach of original Spherical Demons) */
        if (rate != 0.0 && rate != 1.0) {
            fwhm_flow *= rate;
        }

    }
    
    if (debug) {
        output_values_any_format("u.txt", src->n_points, u, TYPE_DOUBLE);
        output_values_any_format("v.txt", src->n_points, v, TYPE_DOUBLE);
        output_values_any_format("warped_curv.txt", src->n_points, curv_src, TYPE_DOUBLE);
    }

    /* invert deformation because we need inverse transformation */
    if (method != 4)
        apply_uv_warp(src_sphere, warped_src_sphere, u, v, 0);

    //delete_polygon_point_neighbours(src, n_neighbours,
    //                neighbours, NULL, NULL);
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
    File_formats   format;
    polygons_struct  *src, *trg, *src_sphere, *trg_sphere;
    polygons_struct  *sm_src, *sm_trg, *sm_src_sphere, *sm_trg_sphere, *warped_src_sphere;
    int        i, step;
    int        n_objects;
    object_struct  **objects;
    double       rotation_matrix[9], rot[3];
    struct       dartel_poly *dpoly_src, *dpoly_trg;
    int        *n_neighbours, **neighbours;
    polygons_struct *orig_sm_src_sphere;

    /* get the arguments from the command line */

    if (ParseArgv(&argc, argv, argTable, 0) ||
      src_file == NULL || trg_file == NULL || 
      src_sphere_file == NULL || trg_sphere_file == NULL ||
      (output_surface_file == NULL && output_sphere_file == NULL)) {
        fprintf(stderr, "\nUsage: %s [options]\n", argv[0]);
        fprintf(stderr, "   %s -help\n\n", argv[0]);
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
    if (input_graphics_any_format(src_sphere_file, &format,
                &n_objects, &objects) != OK)
        exit(EXIT_FAILURE);
    src_sphere = get_polygons_ptr(objects[0]);

    /* read sphere for template surface */
    if (input_graphics_any_format(trg_sphere_file, &format,
                &n_objects, &objects) != OK)
        exit(EXIT_FAILURE);
    trg_sphere = get_polygons_ptr(objects[0]);

    translate_to_center_of_mass(src_sphere);
    for (i = 0; i < src_sphere->n_points; i++)
        set_vector_length(&src_sphere->points[i], 100.0);
    translate_to_center_of_mass(trg_sphere);
    for (i = 0; i < trg_sphere->n_points; i++)
        set_vector_length(&trg_sphere->points[i], 100.0);

    dpoly_src     = (struct dartel_poly *) malloc(sizeof(struct dartel_poly));
    dpoly_trg     = (struct dartel_poly *) malloc(sizeof(struct dartel_poly));

    sm_src        = (polygons_struct *) malloc(sizeof(polygons_struct));
    sm_trg        = (polygons_struct *) malloc(sizeof(polygons_struct));
    sm_src_sphere   = (polygons_struct *) malloc(sizeof(polygons_struct));
    sm_trg_sphere   = (polygons_struct *) malloc(sizeof(polygons_struct));
    warped_src_sphere = (polygons_struct *) malloc(sizeof(polygons_struct));
    orig_sm_src_sphere = (polygons_struct *) malloc(sizeof(polygons_struct));
    
    resample_spherical_surface(trg, trg_sphere, sm_trg, NULL, NULL, n_points);
    resample_spherical_surface(src_sphere, src_sphere, sm_src_sphere, NULL, NULL, n_points);
    resample_spherical_surface(trg_sphere, trg_sphere, sm_trg_sphere, NULL, NULL, n_points);
    resample_spherical_surface(src, src_sphere, sm_src, NULL, NULL, n_points);
    init_dartel_poly(sm_trg_sphere, dpoly_trg);
    init_dartel_poly(sm_src_sphere, dpoly_src);

    printf("Resample surfaces to %d points\n",sm_src->n_points);
    /* Keep copy of original sphere for final resampling */
    copy_polygons(sm_src_sphere, orig_sm_src_sphere);

    int fwhm_flow0 = fwhm_flow;
    
    for (step = 0; step < n_steps; step++) {       

        /* initialization */
        if (step == 0) {
            smooth_heatkernel(sm_src, NULL, fwhm_curv);
            
            /* initial rotation */
            if (rotate) {
                rotate_polygons_to_atlas(sm_src, sm_src_sphere,
                             sm_trg, sm_trg_sphere,
                             10.0, 5, rot, verbose);

                rotation_to_matrix(rotation_matrix,
                           rot[0], rot[1], rot[2]);
        
                /* rotate source sphere */
                rotate_polygons(src_sphere,
                        NULL, rotation_matrix);

                resample_spherical_surface(src_sphere,
                               src_sphere,
                               sm_src_sphere, NULL,
                               NULL, n_points);
            }


        } else if (step > 0)
            curvtype0 = curvtype;
                
        WarpDemon(sm_src, sm_src_sphere, orig_sm_src_sphere, sm_trg, 
             sm_trg_sphere, warped_src_sphere, dpoly_src, dpoly_trg, curvtype0);

        /* update sphere for next step */
        if (n_steps > 0) copy_polygons(warped_src_sphere, sm_src_sphere);
        fwhm_flow = fwhm_flow0;

    }
    printf("\n");
    
    objects = resample_surface_to_target_sphere(orig_sm_src_sphere, warped_src_sphere, src_sphere, NULL, NULL, 0, 0);
    src_sphere = get_polygons_ptr(objects[0]);
    if (output_sphere_file != NULL) {
  
        if (output_graphics_any_format(output_sphere_file, format,
                         n_objects, objects, NULL) != OK)
            exit(EXIT_FAILURE);
    }

    return(EXIT_SUCCESS);
}
