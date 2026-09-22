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

#include <bicpl.h>

typedef struct {
    long x;
    long y;
} Header;

typedef struct {
    double x;
    double y;
} Vector2D;

/**
 * \brief Convert 3D Cartesian point on unit sphere to 2D (u,v) spherical coordinates.
 *
 * \param point (in)  3D point on unit sphere (x² + y² + z² = 1)
 * \param u     (out) horizontal coordinate [0, 1] (azimuth/longitude)
 * \param v     (out) vertical coordinate [0, 1] (elevation/latitude)
 */
void point_to_uv(Point *point, double *u, double *v);
/**
 * \brief Convert 2D (u,v) spherical coordinates to 3D point on unit sphere.
 *
 * \param u     (in)  horizontal coordinate [0, 1] (azimuth/longitude)
 * \param v     (in)  vertical coordinate [0, 1] (elevation/latitude)
 * \param point (out) 3D point on unit sphere with x² + y² + z² = 1
 */
void uv_to_point(double u, double v, Point *point);
/**
 * \brief Map per-vertex values of a surface onto a 2D (u,v) sheet via its sphere.
 *
 * \param polygons      (in)     surface the values belong to
 * \param sphere        (in)     its spherical mapping; NULL creates a
 *                                   tetrahedral unit sphere with the same number of
 *                                   triangles
 * \param sphere_values (in/out) per-vertex values; NULL computes smoothed
 *                                   curvatures of type curvtype, otherwise they are
 *                                   smoothed in-place when fwhm > 0
 * \param mapped_data   (out)    dm[0]*dm[1] sheet values in [0, 1]
 * \param fwhm          (in)     heat-kernel FWHM in mm applied first (0: none)
 * \param dm            (in)     sheet dimensions {nu, nv}
 * \param curvtype      (in)     curvature type used when sphere_values is NULL
 */
void map_sphere_values_to_sheet(polygons_struct *polygons,
                                polygons_struct *sphere, double *sphere_values,
                                double *mapped_data, double fwhm, int *dm,
                                int curvtype);
/**
 * \brief Sample a 2D (u,v) sheet at the vertices of a sphere.
 *
 * \param sheet2d     (in)  dm[0]*dm[1] sheet values
 * \param values      (out) polygons->n_points values, allocated by the caller
 * \param polygons    (in)  mesh whose triangle count defines the sphere
 * \param interpolate (in)  non-zero for bilinear interpolation, zero for the
 *                              nearest pixel
 * \param dm          (in)  sheet dimensions {nu, nv}
 */
void map_sheet2d_to_sphere(double *sheet2d, double *values,
                           polygons_struct *polygons, int interpolate, int *dm);
/**
 * \brief Upsample a 2D flow field by factor of 2 using bilinear interpolation.
 *
 * \param src_flow (in)  source flow field (size src_dm[0]*src_dm[1])
 * \param src_dm   (in)  source dimensions [width, height]
 * \param dst_flow (out) destination flow field (size dst_dm[0]*dst_dm[1])
 * \param dst_dm   (in)  destination dimensions [width, height]
 */
void upsample_flow_field(double *src_flow, int *src_dm, double *dst_flow,
                         int *dst_dm);
/**
 * \brief Downsample a 2D image by factor of 2 using area averaging.
 *
 * \param src (in)  source image (size src_dm[0]*src_dm[1])
 * \param src_dm (in) source dimensions [width, height]
 * \param dst (out) destination image (size dst_dm[0]*dst_dm[1])
 * \param dst_dm (in) destination dimensions [width, height]
 */
void downsample_image(double *src, int *src_dm, double *dst, int *dst_dm);

#endif
