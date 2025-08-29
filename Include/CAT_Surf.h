/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 *
 */

#ifndef _CAT_SURF_H_
#define _CAT_SURF_H_

#include <stdbool.h>
#include <float.h>
#include "bicpl.h"
#include "CAT_Vol.h"
#include "CAT_Math.h"
#include "CAT_Map.h"
#include "CAT_Smooth.h"
#include "CAT_Resample.h"
#include "CAT_Intersect.h"
#include "CAT_Curvature.h"
#include "CAT_Defect.h"
#include "CAT_Deform.h"
#include "CAT_SurfaceIO.h"

#define BINTREE_FACTOR 0.5
#define BW 1024
#define FWHM 30.0
#define DATAFORMAT 1 /* 1 = double data, 0 = complex data */
#define  MAX_NEIGHBOURS  2000
#define _PI 3.14159265358979323846264338327510

/* edge incidence */
typedef struct { int a,b,c, idx; } keytri_t;
typedef struct {
    int a,b;      /* a<b undirected edge key */
    int tri;      /* triangle index */
    double area;  /* triangle area for ranking */
} edge_occ;

int bound(int, int, int *);
void indices(int *, int *, int *);
void get_bounds(polygons_struct *, double [6]);
void apply_warp(polygons_struct *, polygons_struct *, double *, int *, int);
void apply_uv_warp(polygons_struct *, polygons_struct *, double *,
              double *, int );
double get_area_of_points_normalized_to_sphere(polygons_struct *, polygons_struct *, double *);
double * get_surface_ratio(double, polygons_struct *, int);
double get_area_of_points(polygons_struct *, double *);
double get_area_of_polygons(polygons_struct *, double *);
void get_radius_of_points(polygons_struct *, double *);
void localstat_surface_double(polygons_struct *polygons,
                              double *input,
                              unsigned char *mask,
                              int stat_func,
                              int iters);
void translate_to_center_of_mass(polygons_struct *);
void correct_bounds_to_target(polygons_struct *, polygons_struct *);
void correct_bounds_to_target_with_scaling(polygons_struct *, polygons_struct *);
double get_sphere_radius(polygons_struct *);
void set_vector_length(Point *, double);
double compute_point_hausdorff(polygons_struct *, polygons_struct *, double *, int);
double compute_exact_hausdorff(polygons_struct *, polygons_struct *p, double *, int);
double compute_point_distance(polygons_struct *, polygons_struct *, double *, int);
double compute_point_distance_mean(polygons_struct *, polygons_struct *, double *, int);
int euler_characteristic(polygons_struct *, int verbose);
void convert_ellipsoid_to_sphere_with_surface_area(polygons_struct *, double);
void linear_smoothing(polygons_struct *, double, int, int, int *, int);
void areal_smoothing(polygons_struct *, double, int, int, int *, int);
void distance_smoothing(polygons_struct *, double, int, int, int *, int);
void inflate_surface_and_smooth_fingers(polygons_struct *, const int,
                                        const double, const int, const double,
                                        const double, const double, const int);                                        
void surf_to_sphere(polygons_struct *, int, int);
object_struct ** central_to_new_pial(polygons_struct *, double *, double *, 
                                    int, double, int, int);
void central_to_pial(polygons_struct *, double *, double *, int, double, int, int);
double get_area_of_points_central_to_pial(polygons_struct *, double *, double *, double);
int correct_mesh_folding(polygons_struct *, polygons_struct *, float *, 
                         nifti_image *, double);
int reduce_mesh_quadrics(polygons_struct *, int, double, int, int);

#endif
