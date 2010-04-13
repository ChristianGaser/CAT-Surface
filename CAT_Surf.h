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

#define BINTREE_FACTOR 0.5

int bound(int, int, int *);
void get_bounds(polygons_struct *, double [6]);

double * get_surface_ratio(double, polygons_struct *);
double get_area_of_points(polygons_struct *, double *);
double get_area_of_polygons(polygons_struct *, double *);
void get_radius_of_points(polygons_struct *, double *);

void translate_to_center_of_mass(polygons_struct *);
double get_largest_dist(polygons_struct *);

void set_vector_length(Point *, double);
int count_edges(polygons_struct *, int [], int *[]);

int euler_characteristic(polygons_struct *);
void convert_ellipsoid_to_sphere_with_surface_area(polygons_struct *, double);

void linear_smoothing(polygons_struct *, double, int, int, int *, int);
void areal_smoothing(polygons_struct *, double, int, int, int *, int);
void distance_smoothing(polygons_struct *, double, int, int, int *, int);

void inflate_surface_and_smooth_fingers(polygons_struct *, const int,
                                        const double, const int, const double,
                                        const double, const double, const int);
