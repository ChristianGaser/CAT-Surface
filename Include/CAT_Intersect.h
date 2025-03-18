/* Rachel Yotter - rachel.yotter@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 *
 */

#ifndef _CAT_INTERSECT_H_
#define _CAT_INTERSECT_H_

#include <bicpl.h>

#define PINF  1.7976931348623157e+308 /* for doubles */
#define NINF -1.7976931348623157e+308 /* for doubles */

#define BINTREE_FACTOR 0.5

int intersect_poly_poly(int, int, polygons_struct *);
int intersect_triangle_triangle(int [3], int [3], polygons_struct *);
int intersect_segment_triangle(Point, Point, int [3], polygons_struct *);
int find_selfintersections(polygons_struct *, int *, int *, int);
int find_selfintersections_masked(polygons_struct *, int *, int *);
int correct_simple_selfintersections(polygons_struct *, int *, int *, int *, int **);
int join_intersections(polygons_struct *, int *, int *, int *, int **);
int find_remaining_intersections(polygons_struct *, int *, int *, int *,
                                 int **);
int patch_selfintersections(polygons_struct *, polygons_struct *, int *, int *,
                            int, int *, int **);
int smooth_selfintersections(polygons_struct *, int *, int *, int, int *,
                             int **, int);
int has_selfintersections(polygons_struct *, int *, int);
void remove_intersections(polygons_struct *, int);

#endif