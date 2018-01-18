/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 *
 */

#include <bicpl.h>

#define BINTREE_FACTOR 0.5
#define BW 1024
#define FWHM 30.0
#define DATAFORMAT 1 /* 1 = double data, 0 = complex data */

int bound(int, int, int *);
void indices(int *, int *, int *);
void get_bounds(polygons_struct *, double [6]);
void apply_warp(polygons_struct *, polygons_struct *, double *, int *, int);
void apply_uv_warp(polygons_struct *, polygons_struct *, double *,
              double *, int );
double get_area_of_points_normalized_to_sphere(polygons_struct *, polygons_struct *, double *);
double * get_surface_ratio(double, polygons_struct *);
double get_area_of_points(polygons_struct *, double *);
double get_area_of_polygons(polygons_struct *, double *);
void get_radius_of_points(polygons_struct *, double *);
void translate_to_center_of_mass(polygons_struct *);
void correct_bounds_to_target(polygons_struct *, polygons_struct *);
double get_sphere_radius(polygons_struct *);
void set_vector_length(Point *, double);
int count_edges(polygons_struct *, int [], int *[]);
double compute_point_hausdorff(polygons_struct *, polygons_struct *, double *, int);
double compute_exact_hausdorff(polygons_struct *, polygons_struct *p, double *);
void compute_point_distance(polygons_struct *, polygons_struct *, double *);
int euler_characteristic(polygons_struct *);
void convert_ellipsoid_to_sphere_with_surface_area(polygons_struct *, double);
void linear_smoothing(polygons_struct *, double, int, int, int *, int);
void areal_smoothing(polygons_struct *, double, int, int, int *, int);
void distance_smoothing(polygons_struct *, double, int, int, int *, int);
void inflate_surface_and_smooth_fingers(polygons_struct *, const int,
                                        const double, const int, const double,
                                        const double, const double, const int);                                        
void surf_to_sphere(polygons_struct *, int);
object_struct ** central_to_new_pial(polygons_struct *, double *, double *);
void central_to_pial(polygons_struct *, double *, double *);
