/* Rachel Yotter - rachel.yotter@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id: CAT_Intersect.h 121 2009-07-27 14:42:19Z raytrace $
 *
 */

#include <volume_io/internal_volume_io.h>
#include <volume_io/geometry.h>
#include <bicpl.h>

#define PINF  1.7976931348623157e+308 /* for doubles */
#define NINF -1.7976931348623157e+308 /* for doubles */

#define BINTREE_FACTOR 0.5

int intersect_triangle_triangle(int [3], int [3], polygons_struct *);
int intersect_segment_triangle(Point, Point, int [3], polygons_struct *);
int find_selfintersections(polygons_struct *, int *);
int join_intersections(polygons_struct *, int *, int *, int **);
int join_intersections_poly(polygons_struct *, int *, int *, int *, int **);
int find_topological_defects(polygons_struct *, int *, int *, int **);
void expand_defects(polygons_struct *, int *, int, int, int *, int **);
void expand_defects_poly(polygons_struct *, int *, int *, int, int, int *,
                         int **);
void update_polydefects(polygons_struct *, int *, int *);
void remap_defect_points(polygons_struct *, int *, polygons_struct *, int *);
void repair_selfintersections(polygons_struct *, int *, int **);
int count_defect_intersections(polygons_struct *, int *, int);
