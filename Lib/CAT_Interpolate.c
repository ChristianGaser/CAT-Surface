/* Rachel Yotter - rachel.yotter@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 *
 */

#include "CAT_Interpolate.h"
#include "CAT_Map.h"

/**
 * \brief Interpolate scalar value at a 3D point on a unit sphere mesh.
 *
 * Uses barycentric interpolation within the closest polygon to compute
 * a scalar value at an arbitrary 3D coordinate. The sphere is assumed
 * to have unit radius. Builds a spatial index (bintree) on first call.
 *
 * \param sphere    (in)  polygon mesh representing the unit sphere
 * \param values    (in)  double[sphere->n_points]; scalar values at each vertex
 * \param pt        (in)  3D point at which to interpolate (should be near sphere)
 * \return               Interpolated scalar value using barycentric weighting
 */
double
interp_point_unit_sphere(polygons_struct *sphere, double *values, Point pt)
{
    Point on_sphere_pt, poly_pts[1000];
    double value;
    double weights[1000];
    int i, poly, size, ind;

    if (sphere->bintree == NULL)
        create_polygons_bintree(sphere, ROUND((double)sphere->n_items *
                                              BINTREE_FACTOR));

    poly = find_closest_polygon_point(&pt, sphere, &on_sphere_pt);
    size = get_polygon_points(sphere, poly, poly_pts);
    get_polygon_interpolation_weights(&on_sphere_pt, size, poly_pts,
                                      weights);
    value = 0.0;
    for (i = 0; i < size; i++)
    {
        ind = sphere->indices[POINT_INDEX(sphere->end_indices,
                                          poly, i)];
        value += (double)weights[i] * values[ind];
    }

    return value;
}

/**
 * \brief Interpolate scalar value at a 3D point on a sphere mesh of arbitrary radius.
 *
 * Similar to interp_point_unit_sphere() but works with spheres of any radius.
 * Internally normalizes the sphere to unit radius, performs interpolation,
 * then normalizes back. Useful for working with scaled surface meshes.
 *
 * \param sphere    (in)  polygon mesh representing a sphere (any radius)
 * \param values    (in)  double[sphere->n_points]; scalar values at each vertex
 * \param pt        (in)  3D point at which to interpolate
 * \return               Interpolated scalar value
 */
double
interp_point_sphere(polygons_struct *sphere, double *values, Point pt)
{
    polygons_struct *unit_sphere;
    double value;
    int i;

    unit_sphere = (polygons_struct *)malloc(sizeof(polygons_struct));
    copy_polygons(sphere, unit_sphere);

    /* set radius to 1 */
    for (i = 0; i < unit_sphere->n_points; i++)
        set_vector_length(&unit_sphere->points[i], 1.0);

    value = interp_point_unit_sphere(unit_sphere, values, pt);
    free(unit_sphere);
    return (value);
}

/**
 * \brief Interpolate scalar value at (u,v) latitude/longitude coordinates on unit sphere.
 *
 * Converts (u,v) spherical coordinates to 3D Cartesian coordinates,
 * then interpolates using barycentric weighting on the unit sphere mesh.
 *
 * \param sphere    (in)  polygon mesh representing the unit sphere
 * \param values    (in)  double[sphere->n_points]; scalar values at each vertex
 * \param u         (in)  first spherical coordinate (0..1 or other range)
 * \param v         (in)  second spherical coordinate (0..1 or other range)
 * \return               Interpolated scalar value at (u,v)
 */
double
interp_uv_unit_sphere(polygons_struct *sphere, double *values,
                      double u, double v)
{
    Point pt;

    uv_to_point(u, v, &pt);
    return (interp_point_unit_sphere(sphere, values, pt));
}

/**
 * \brief Interpolate scalar value at (u,v) latitude/longitude coordinates on arbitrary sphere.
 *
 * Wrapper around spherical interpolation that works with meshes of any radius.
 * Converts (u,v) coordinates to 3D Cartesian, then interpolates on the sphere.
 *
 * \param sphere    (in)  polygon mesh representing a sphere (any radius)
 * \param values    (in)  double[sphere->n_points]; scalar values at each vertex
 * \param u         (in)  first spherical coordinate
 * \param v         (in)  second spherical coordinate
 * \return               Interpolated scalar value at (u,v)
 */
double
interp_uv_sphere(polygons_struct *sphere, double *values, double u, double v)
{
    Point pt;

    uv_to_point(u, v, &pt);
    return (interp_point_sphere(sphere, values, pt));
}

/**
 * \brief Interpolate scalar value at Cartesian coordinates (x,y,z) on arbitrary sphere.
 *
 * Performs barycentric interpolation at given 3D Cartesian coordinates
 * on a sphere mesh of any radius. Wraps interp_point_sphere().
 *
 * \param sphere    (in)  polygon mesh representing a sphere (any radius)
 * \param values    (in)  double[sphere->n_points]; scalar values at each vertex
 * \param x         (in)  x-coordinate of interpolation point
 * \param y         (in)  y-coordinate of interpolation point
 * \param z         (in)  z-coordinate of interpolation point
 * \return               Interpolated scalar value at (x,y,z)
 */
double
interp_xyz_sphere(polygons_struct *sphere, double *values, double x, double y,
                  double z)
{
    Point pt;

    fill_Point(pt, x, y, z);
    return (interp_point_sphere(sphere, values, pt));
}

/**
 * \brief Interpolate scalar value at Cartesian coordinates (x,y,z) on unit sphere.
 *
 * Performs barycentric interpolation at given 3D Cartesian coordinates
 * on a unit radius sphere mesh. Wraps interp_point_unit_sphere().
 *
 * \param sphere    (in)  polygon mesh representing the unit sphere
 * \param values    (in)  double[sphere->n_points]; scalar values at each vertex
 * \param x         (in)  x-coordinate of interpolation point
 * \param y         (in)  y-coordinate of interpolation point
 * \param z         (in)  z-coordinate of interpolation point
 * \return               Interpolated scalar value at (x,y,z)
 */
double
interp_xyz_unit_sphere(polygons_struct *sphere, double *values,
                       double x, double y, double z)
{
    Point pt;

    fill_Point(pt, x, y, z);
    return (interp_point_unit_sphere(sphere, values, pt));
}
