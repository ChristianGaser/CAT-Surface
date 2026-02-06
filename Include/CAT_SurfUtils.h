/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 *
 */

#ifndef _CAT_SURFUTILS_H_
#define _CAT_SURFUTILS_H_

#include <bicpl.h>

/*
 * Lightweight surface-utility declarations.
 *
 * This header intentionally contains only small, broadly used helpers that are
 * implemented in Lib/CAT_Surf.c, but are needed by multiple other library
 * modules. It exists to avoid pulling the full umbrella header CAT_Surf.h
 * (and its large include/dependency tree) into modules that only need one or
 * two simple utilities.
 */

#ifndef BINTREE_FACTOR
#define BINTREE_FACTOR 0.5
#endif

#ifndef MAX_NEIGHBOURS
#define MAX_NEIGHBOURS 2000
#endif

int bound(int i, int j, int dm[]);
void get_bounds(polygons_struct *polygons, double bounds[6]);
void set_vector_length(Point *p, double newLength);

void translate_to_center_of_mass(polygons_struct *polygons);
double get_area_of_polygons(polygons_struct *polygons, double *areas);
double get_area_of_points_normalized_to_sphere(polygons_struct *polygons,
											   polygons_struct *sphere,
											   double *areas);
void surf_to_sphere(polygons_struct *polygons, int n_triangles, int verbose);

#endif
