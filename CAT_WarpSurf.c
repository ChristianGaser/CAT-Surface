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
char *jacdet_file        = NULL;
char *output_file        = NULL;
char *output_sphere_file = NULL;
char *pgm_file           = NULL;

int rotate      = 1;
int avg         = 0;
int code        = 1;
int loop        = 6;
int verbose     = 0;
int rtype       = 1;
int curvtype    = 3;
int muchange    = 4;
int sz_map[2]   = {500, 300};
int n_triangles = 81920;
int n_steps     = 3;
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
     "Filter size for curvature map in FWHM."},
  {"-loop", ARGV_INT, (char *) 1, (char *) &loop,
     "Number of outer loops for default parameters (max. 6)."},
  {"-steps", ARGV_INT, (char *) 1, (char *) &n_steps,
     "Number of Dartel steps (max. 6)."},
  {"-size", ARGV_INT, (char *) 2, (char *) &sz_map,
     "Size of curvature map for warping."},
  {"-norot", ARGV_CONSTANT, (char *) FALSE, (char *) &rotate,
     "Don't rotate input surface before warping."},
  {"-avg", ARGV_CONSTANT, (char *) TRUE, (char *) &avg,
     "Average together two weighted DARTEL solutions into final mesh."},
  {"-type", ARGV_INT, (char *) 1, (char *) &curvtype,
     "Curvature type\n\t0 - mean curvature (averaged over 3mm, in degrees)\n\t1 - gaussian curvature\n\t2 - curvedness\n\t3 - shape index\n\t4 - mean curvature (in radians)."},
  {"-v", ARGV_CONSTANT, (char *) TRUE, (char *) &verbose,
     "Be verbose."},
  {"-debug", ARGV_CONSTANT, (char *) TRUE, (char *) &debug,
     "Save debug files."},
   {NULL, ARGV_END, NULL, NULL, NULL}
};

void
resample_spherical_surface(polygons_struct *polygons,
                           polygons_struct *poly_src_sphere,
                           polygons_struct *resampled_source,
                           double *input_values, double *output_values,
                           int n_triangles)
{
        int    i, k, poly, n_points;
        int    *n_neighbours, **neighbours;
        Point  centre, point_on_src_sphere, scaled_point;
        Point  poly_points[MAX_POINTS_PER_POLYGON];
        Point  poly_points_src[MAX_POINTS_PER_POLYGON];
        Point  *new_points;
        Real   weights[MAX_POINTS_PER_POLYGON];
        double sphereRadius, r, bounds[6];

        /*
         * Determine radius for the output sphere.  The sphere is not always
         * perfectly spherical, thus use average radius
         */
        sphereRadius = 0.0;
        for (i = 0; i < poly_src_sphere->n_points; i++) {
                r = 0.0;
                for (k = 0; k < 3; k++) 
                        r += Point_coord(poly_src_sphere->points[i], k) *
                             Point_coord(poly_src_sphere->points[i], k);
                sphereRadius += sqrt(r);
        }
        sphereRadius /= poly_src_sphere->n_points;

        /* Calc. sphere center based on bounds of input (correct for shifts) */
        get_bounds(poly_src_sphere, bounds);
        fill_Point(centre, bounds[0]+bounds[1],
                           bounds[2]+bounds[3], bounds[4]+bounds[5]);
    
        /*
         * Make radius slightly smaller to get sure that the
         * inner side of handles will be found as nearest point on the surface
         */
        sphereRadius *= 0.975;
        create_tetrahedral_sphere(&centre, sphereRadius, sphereRadius,
                                  sphereRadius, n_triangles, resampled_source);

        create_polygons_bintree(poly_src_sphere,
                                ROUND((Real) poly_src_sphere->n_items * 0.5));

        ALLOC(new_points, resampled_source->n_points);
        if (input_values != NULL)
                ALLOC(output_values, resampled_source->n_points);

        for (i = 0; i < resampled_source->n_points; i++) {
                poly = find_closest_polygon_point(&resampled_source->points[i],
                                                  poly_src_sphere,
                                                  &point_on_src_sphere);
		
                n_points = get_polygon_points(poly_src_sphere, poly,
                                              poly_points_src);
                get_polygon_interpolation_weights(&point_on_src_sphere,
                                                  n_points, poly_points_src,
                                                  weights);

                if (get_polygon_points(polygons, poly, poly_points) != n_points)
                        handle_internal_error("map_point_between_polygons");

                fill_Point(new_points[i], 0.0, 0.0, 0.0);
                if (input_values != NULL)
                        output_values[i] = 0.0;

                for (k = 0; k < n_points; k++) {
                        SCALE_POINT(scaled_point, poly_points[k], weights[k]);
                        ADD_POINTS(new_points[i], new_points[i], scaled_point);
                        if (input_values != NULL)
                                output_values[i] += weights[k] *
                                    input_values[polygons->indices[
                                    POINT_INDEX(polygons->end_indices,poly,k)]];
                }
       }

        create_polygon_point_neighbours(resampled_source, TRUE, &n_neighbours,
                                        &neighbours, NULL, NULL);

        for (i = 0; i < resampled_source->n_points; i++) {
                resampled_source->points[i] = new_points[i];
        }
		
        compute_polygon_normals(resampled_source);
        free(new_points);
}

/* input 2 surfaces w/ weighting along x- or z-axes, output weighted average */
void
average_xz_surf(polygons_struct *xsurf, polygons_struct *zsurf,
                polygons_struct *surface)
{
        double xx, xy, xz, zx, zy, zz, phi, wx, wz, wtot, dot;
        int p, *yflag, *zflag;

        copy_polygons(zsurf, surface);
        compute_polygon_normals(xsurf);
        compute_polygon_normals(zsurf);

        for (p = 0; p < surface->n_points; p++) {
                xx = Point_x(xsurf->points[p]);
                xy = Point_y(xsurf->points[p]);
                xz = Point_z(xsurf->points[p]);
                zx = Point_x(zsurf->points[p]);
                zy = Point_y(zsurf->points[p]);
                zz = Point_z(zsurf->points[p]);

                phi = acos(xx) / PI;
                wx = 0.9 - 16.0 * pow(phi - 0.5, 4.0);
                if (wx < 0.0) wx = 0.0;
                phi = acos(zz) / PI;
                wz = 0.9 - 16.0 * pow(phi - 0.5, 4.0);
                if (wz < 0.0) wz = 0.0;

                wtot = wx + wz;
                if (wtot > 0) {
                        wx /= wtot;
                        wz /= wtot;
                } else { /* just average */
                        wx = 0.5;
                        wz = 0.5;
                }

                fill_Point(surface->points[p], wx*xx + wz*zx,
                           wx*xy + wz*zy, wx*xz + wz*zz);
                set_vector_length(&surface->points[p], 1.0);
        }
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
        map_trg = (double *) malloc(sizeof(double) * src->n_points);
        map_src = (double *) malloc(sizeof(double) * src->n_points);

        get_smoothed_curvatures(trg, trg_sphere, orig_trg,
                                fwhm, curvtype);
        get_smoothed_curvatures(src, src_sphere, map_src,
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
                                        resample_noscale(trg_sphere,
                                                         &rot_src_sphere,
                                                         orig_trg, map_trg);

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
                                                fprintf(stderr,"alpha: %5.3f\tbeta: %5.3f\tgamma: %5.3f\tsquared difference: %5.3f",
                                                        DEGREES(alpha),
                                                        DEGREES(beta),
                                                        DEGREES(gamma), sum_sq);
                                                fprintf(stderr, "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b");
                                                fprintf(stderr, "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b");
                                                fprintf(stderr, "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b");
                                                fprintf(stderr, "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b");
                                        }
                                }
                        }
                }
                
                /* save best estimates */
                curr_alpha = best_alpha;
                curr_beta  = best_beta;
                curr_gamma = best_gamma;
        }
        fprintf(stderr,"\n");
        free(orig_trg);
        free(map_trg);
        free(map_src);
}

void
solve_dartel_flow(polygons_struct *src, polygons_struct *src_sphere,
                  polygons_struct *trg, polygons_struct *trg_sphere,
                  struct dartel_prm *prm, int dm[3], int n_steps,
                  double rot[3], double *flow)
{
        int step, i, x, y, it, it0, it1, xy_size, it_scratch;
        polygons_struct *sm_src, *sm_trg, *sm_src_sphere, *sm_trg_sphere;
        double           H00, H01, H10, H11, rotation_matrix[9];
        double *flow1, *inflow, *map_src, *map_trg, *map_warp;
        polygons_struct *rot_sphere;
        double  *scratch, *jd, *jd1, *values;
        char buffer[1024];
        double           ll[3];
        double           xp, yp, xm, ym;

        xy_size = dm[0] * dm[1];

        flow1       = (double *) malloc(sizeof(double) * xy_size * 2);
        inflow      = (double *) malloc(sizeof(double) * xy_size * 2);
        map_src  = (double *) malloc(sizeof(double) * xy_size);
        map_trg  = (double *) malloc(sizeof(double) * xy_size);
        map_warp    = (double *) malloc(sizeof(double) * xy_size);
        sm_src     = (polygons_struct *) malloc(sizeof(polygons_struct));
        sm_src_sphere = (polygons_struct *) malloc(sizeof(polygons_struct));
        sm_trg     = (polygons_struct *) malloc(sizeof(polygons_struct));
        sm_trg_sphere = (polygons_struct *) malloc(sizeof(polygons_struct));
        rot_sphere    = (polygons_struct *) malloc(sizeof(polygons_struct));

        copy_polygons(src_sphere, rot_sphere);

        for (step = 0; step < n_steps; step++) {
                /* resample source and target surface */
                resample_spherical_surface(src, rot_sphere, sm_src, NULL, NULL,
                                           n_triangles);
                resample_spherical_surface(trg, trg_sphere,
                                           sm_trg, NULL, NULL,
                                           n_triangles);

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
                        
                        /* save pgm for debugging */
                        if (debug) {
                                get_smoothed_curvatures(sm_src, rot_sphere,
                                                        map_src, fwhm,
                                                        curvtype);
                                sprintf(buffer,"source%d.txt",curvtype);
                                if (output_values_any_format(buffer,
                                                   src->n_points, map_src,
                                                   TYPE_DOUBLE) != 0)
                                        exit(EXIT_FAILURE);
                        }

                        resample_spherical_surface(rot_sphere, rot_sphere,
                                                   sm_src_sphere, NULL, NULL,
                                                   n_triangles);
                        resample_spherical_surface(trg_sphere, trg_sphere,
                                                   sm_trg_sphere, NULL, NULL,
                                                   n_triangles);

                        /* initial rotation */
                        if (rotate) {
                                rotate_polygons_to_atlas(sm_src, sm_src_sphere,
                                                         sm_trg, sm_trg_sphere,
                                                         fwhm, curvtype, rot);
                                rotation_to_matrix(rotation_matrix,
                                                   rot[0], rot[1], rot[2]);
                
                                /* rotate source sphere */
                                rotate_polygons(rot_sphere,
                                                NULL, rotation_matrix);

                                resample_spherical_surface(rot_sphere,
                                                           rot_sphere,
                                                           sm_src_sphere, NULL,
                                                           NULL, n_triangles);
                                //*get_polygons_ptr(objects[0]) = *rot_sphere;
                                //output_graphics_any_format("rotated.obj",
                                                         //format, 1, objects);
                        }

                        /* init warps and flows with zeros */
                        for (i = 0; i < xy_size;   i++)  map_warp[i]  = 0.0;
                        for (i = 0; i < xy_size*2; i++)  inflow[i] = 0.0;

                } else if (step == 1) {
                        /* smooth surfaces */
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
                        /* smooth surfaces */
                        inflate_surface_and_smooth_fingers(sm_src,
                                                           1, 0.2, 20, 1.0,
                                                           3.0, 1.0, 0);
                        inflate_surface_and_smooth_fingers(sm_trg,
                                                           1, 0.2, 20, 1.0,
                                                           3.0, 1.0, 0);
                }

                /* get curvatures */
                map_smoothed_curvature_to_sphere(sm_trg, NULL, (double *)0,
                                                 map_trg, fwhm, dm, curvtype);
                map_smoothed_curvature_to_sphere(sm_src, sm_src_sphere,
                                                 (double *)0, map_src, fwhm,
                                                 dm, curvtype);

                /* save pgm images for debugging */
                if (debug) {
                        values = (double *) malloc(sizeof(double) *
                                                   rot_sphere->n_points);

                        sprintf(buffer,"source%d%d.pgm",curvtype,step+1);
                        if (write_pgm(buffer, map_src, dm[0], dm[1]) != 0)
                                exit(EXIT_FAILURE);
                        sprintf(buffer,"source%d%d.txt",curvtype,step+1);
                        map_sheet2d_to_unit_sphere(map_src, values,
                                                   rot_sphere, 1, dm);
                        if (output_values_any_format(buffer,
                                                     rot_sphere->n_points,
                                                     values, TYPE_DOUBLE) != 0)
                                exit(EXIT_FAILURE);

                        sprintf(buffer,"target%d%d.pgm",curvtype,step+1);
                        if (write_pgm(buffer, map_trg, dm[0], dm[1]) != 0)
                                exit(EXIT_FAILURE);

                        sprintf(buffer,"target%d%d.txt",curvtype,step+1);
                        map_sheet2d_to_unit_sphere(map_trg, values,
                                                   rot_sphere, 1, dm);
                        if (output_values_any_format(buffer,
                                                     rot_sphere->n_points,
                                                     values, TYPE_DOUBLE) != 0)
                                exit(EXIT_FAILURE);
                        free(values);
                }

                /* go through dartel steps */
                for (it = 0, it0 = 0; it0 < loop; it0++) {
                        it_scratch = dartel_scratchsize((int *)dm,
                                                         prm[it0].code);
                        scratch = (double *) malloc(sizeof(double)*it_scratch);
                        for (it1 = 0; it1 < prm[it0].its; it1++) {
                                it++;
                                /* map target onto source */
                                if (INVERSE_WARPING) {
                                        dartel(prm[it0], dm, inflow, map_src,
                                              map_trg, NULL, flow, ll, scratch);
                                } else {
                                        dartel(prm[it0], dm, inflow, map_trg,
                                              map_src, NULL, flow, ll, scratch);
                                }
                                fprintf(stderr, "%02d-%02d: %8.2f\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b", step+1, it, ll[0]);
                                for (i = 0; i < xy_size*2; i++)
                                        inflow[i] = flow[i];
                        }
                        free(scratch);
                }

                /* use smaller FWHM for next steps */
                //if (fwhm > 4) fwhm -= 2.0;
                fwhm /= 2.0;
        }
        fprintf(stderr,"\n");

        /* get deformations and jacobian det. from flow field */
        if (jacdet_file != NULL) {
                fprintf(stderr,"Warning: Saving jacobians not working\n");
                if (rotate)
                       fprintf(stderr,"Warning: Rotation not yet considered\n");
                jd  = (double *) malloc(sizeof(double) * xy_size);
                jd1 = (double *) malloc(sizeof(double) * xy_size);

                expdefdet(dm, 10, inflow, flow, flow1, jd, jd1);
                
                /* subtract 1 to get values around 0 instead of 1 and invert */
                for (i = 0; i < xy_size; i++) {
                        jd1[i] -= 1;
                        jd1[i] *= -1;
                }

                values = (double *) malloc(sizeof(double) *
                                           src_sphere->n_points);

                map_sheet2d_to_sphere(jd1, values, rot_sphere, 1, dm);

                output_values_any_format(jacdet_file, rot_sphere->n_points,
                                         values, TYPE_DOUBLE);

                free(values);
                free(jd1);
                free(jd);
        } else {
                expdef(dm, 10, inflow, flow, flow1,
                       (double *) 0, (double *) 0);
        }

        free(flow1);
        free(inflow);

        if (pgm_file != NULL) {
                fprintf(stderr,"Warning: Inversion of deformation field is not prepared.\n");
                for (i = 0; i < xy_size; i++) {
                        x = (int) flow[i] - 1.0;
                        y = (int) flow[i + xy_size] - 1.0;
                        xp = flow[i] - x - 1.0;
                        yp = flow[i + xy_size] - y - 1.0;
                        xm = 1.0 - xp;
                        ym = 1.0 - yp;
                        H00 = map_src[bound(x,  y,  dm)];
                        H01 = map_src[bound(x,  y+1,dm)];
                        H10 = map_src[bound(x+1,y,  dm)];
                        H11 = map_src[bound(x+1,y+1,dm)];

                        map_warp[i] = ym * (xm * H00 + xp * H10) +
		                      yp * (xm * H01 + xp * H11);
                }
                if (write_pgm(pgm_file, map_warp, dm[0], dm[1]) != 0)
                        exit(EXIT_FAILURE);
        }

        free(map_src);
        free(map_warp);
        free(map_trg);
}

int
main(int argc, char *argv[])
{
        File_formats     format;
        FILE             *fp;
        char             line[1024], buffer[1024];
        polygons_struct  *src, *trg, *src_sphere, *trg_sphere;
        polygons_struct  *rsrc, *rtrg, *rs_sph, *rt_sph;
        polygons_struct  *asrc, *as_sph;
        int              x, y, i, j, p;
        int              n_objects, xy_size;
        double           *flow, *flow2;
        object_struct    **objects;
        static double    param[2] = {1.0, 1.0};
        int              dm[3];
        double           xp, yp, xm, ym;
        double           H00, H01, H10, H11, rotation_matrix[9];
        struct           dartel_prm* prm;
        double           rot[3];


        /* get the arguments from the command line */
        if (ParseArgv(&argc, argv, argTable, 0) ||
            source_file == NULL || target_file == NULL || 
            source_sphere_file == NULL || target_sphere_file == NULL ||
            (jacdet_file == NULL && output_file == NULL &&
             pgm_file == NULL && output_sphere_file == NULL)) {
                fprintf(stderr, "\nUsage: %s [options]\n", argv[0]);
                fprintf(stderr, "     %s -help\n\n", argv[0]);
                exit(EXIT_FAILURE);
        }

        if (input_graphics_any_format(target_file, &format,
                                      &n_objects, &objects) != OK)
                exit(EXIT_FAILURE);

        /* get a pointer to the surface */
        trg = get_polygons_ptr(objects[0]);

        /* check that the surface file contains a polyhedron */
        if (n_objects != 1 || get_object_type(objects[0]) != POLYGONS) {
                printf("Surface file must contain 1 polygons object.\n");
                exit(EXIT_FAILURE);
        }

        if (input_graphics_any_format(source_file, &format,
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
        
        /* size of curvature map for warping */
        dm[0] = sz_map[0];
        dm[1] = sz_map[1];
        dm[2] = 1;

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
    
        xy_size = dm[0] * dm[1];
        flow    = (double *) malloc(sizeof(double) * xy_size * 2);

        solve_dartel_flow(src, src_sphere, trg, trg_sphere, prm, dm, n_steps,
                          rot, flow);
        rotation_to_matrix(rotation_matrix, rot[0], rot[1], rot[2]);
        rotate_polygons(src_sphere, NULL, rotation_matrix);
        rotate = 0;

        /* solve again, but rotated to change pole location */
        if (avg) {
                flow2 = (double *) malloc(sizeof(double) * xy_size * 2);
                rsrc = (polygons_struct *) malloc(sizeof(polygons_struct));
                rs_sph  = (polygons_struct *) malloc(sizeof(polygons_struct));

                rotation_to_matrix(rotation_matrix, 0.0, PI/2.0, 0.0);

                rotate_polygons(src, rsrc, rotation_matrix);
                rotate_polygons(src_sphere, rs_sph, rotation_matrix);
                rotate_polygons(trg, NULL, rotation_matrix);
                rotate_polygons(trg_sphere, NULL, rotation_matrix);

                solve_dartel_flow(rsrc, rs_sph, trg, trg_sphere, prm,
                                  dm, n_steps, rot, flow2);

                rotation_to_matrix(rotation_matrix, 0.0, -PI/2.0, 0.0);

                if (output_file != NULL) {
                        apply_warp(rsrc, rs_sph, flow2, dm, INVERSE_WARPING); 
  
                        asrc = (polygons_struct *)
                                                malloc(sizeof(polygons_struct));
                        rotate_polygons(rsrc, NULL, rotation_matrix);
                        average_xz_surf(rsrc, src, asrc);
                        copy_polygons(asrc, src);
                        free(asrc);
                }

                if (output_sphere_file != NULL) {
                        apply_warp(rs_sph, rs_sph, flow2, dm, INVERSE_WARPING); 

                        as_sph  = (polygons_struct *)
                                                malloc(sizeof(polygons_struct));
                        rotate_polygons(rs_sph, NULL, rotation_matrix);
                        average_xz_surf(rs_sph, src_sphere, as_sph);
                        copy_polygons(as_sph, src_sphere);
                        free(as_sph);
                }
                free(flow2);
                free(rsrc);
                free(rs_sph);
        } else {
                if (output_file != NULL) {
                        /* apply inverse deformations */
                        apply_warp(src, src_sphere, flow, dm,
                                   INVERSE_WARPING); 
                }

                if (output_sphere_file != NULL) {
                        apply_warp(src_sphere, src_sphere, flow, dm,
                                   INVERSE_WARPING);
                }
        }

        if (output_file != NULL) {
                compute_polygon_normals(src);
                *get_polygons_ptr(objects[0]) = *src;
                if (output_graphics_any_format(output_file, format, n_objects,
                                               objects) != OK)
                        exit(EXIT_FAILURE);
        }

        if (output_sphere_file != NULL) {
                compute_polygon_normals(src_sphere);
                *get_polygons_ptr(objects[0]) = *src_sphere;
                if (output_graphics_any_format(output_sphere_file, format,
                                               n_objects, objects) != OK)
                        exit(EXIT_FAILURE);
        }

        delete_object_list(n_objects, objects);
        free(flow);
        free(prm);

        return(EXIT_SUCCESS);
}
