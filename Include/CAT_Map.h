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

/**
 * \brief Public API for point_to_uv.
 *
 * This function is part of the CAT-Surface public library interface and is used by command-line tools.
 *
 * \param param (in/out) Parameter of point_to_uv.
 * \param param (in/out) Parameter of point_to_uv.
 * \param param (in/out) Parameter of point_to_uv.
 * \return void (no return value).
 */
void point_to_uv(Point *, double *, double *);
/**
 * \brief Public API for uv_to_point.
 *
 * This function is part of the CAT-Surface public library interface and is used by command-line tools.
 *
 * \param double (in/out) Parameter of uv_to_point.
 * \param double (in/out) Parameter of uv_to_point.
 * \param param (in/out) Parameter of uv_to_point.
 * \return void (no return value).
 */
void uv_to_point(double, double, Point *);
void map_sphere_values_to_sheet(polygons_struct *, polygons_struct *,
                                      double *, double *, double, int *, int);
/**
 * \brief Public API for map_sheet2d_to_sphere.
 *
 * This function is part of the CAT-Surface public library interface and is used by command-line tools.
 *
 * \param param (in/out) Parameter of map_sheet2d_to_sphere.
 * \param param (in/out) Parameter of map_sheet2d_to_sphere.
 * \param param (in/out) Parameter of map_sheet2d_to_sphere.
 * \param int (in/out) Parameter of map_sheet2d_to_sphere.
 * \param param (in/out) Parameter of map_sheet2d_to_sphere.
 * \return void (no return value).
 */
void map_sheet2d_to_sphere(double *, double *, polygons_struct *, int, int *);
void upsample_flow_field(double *src_flow, int *src_dm, double *dst_flow, int *dst_dm);
void downsample_image(double *src, int *src_dm, double *dst, int *dst_dm);

#endif
