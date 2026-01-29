/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 *
*/

#ifndef _CAT_MAP_H_
#define _CAT_MAP_H_

typedef struct {
    long x;
    long y;
} Header;

typedef struct {
    double x;
    double y;
} Vector2D;

void point_to_uv(Point *, double *, double *);
void uv_to_point(double, double, Point *);
void map_sphere_values_to_sheet(polygons_struct *, polygons_struct *,
                                      double *, double *, double, int *, int);
void map_sheet2d_to_sphere(double *, double *, polygons_struct *, int, int *);
void upsample_flow_field(double *src_flow, int *src_dm, double *dst_flow, int *dst_dm);
void downsample_image(double *src, int *src_dm, double *dst, int *dst_dm);

#endif
