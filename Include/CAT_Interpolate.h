/* Rachel Yotter - rachel.yotter@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 *
 */

#ifndef _CAT_INTERPOLATE_H_
#define _CAT_INTERPOLATE_H_

#include <bicpl.h>
#include <float.h>

#include "CAT_Surf.h"
#include "CAT_Smooth.h"
#include "CAT_Curvature.h"

#define  BINTREE_FACTOR   0.5
#define  NEW_COORDINATE_SYSTEM   1

/**
 * \brief Interpolate scalar value at a 3D point on a unit sphere mesh.
 *
 * \param sphere    (in)  polygon mesh representing the unit sphere
 * \param values    (in)  double[sphere->n_points]; scalar values at each vertex
 * \param pt        (in)  3D point at which to interpolate (should be near sphere)
 * \return               Interpolated scalar value using barycentric weighting
 */
double interp_point_unit_sphere(polygons_struct *sphere, double *values,
                                Point pt);
/**
 * \brief Interpolate scalar value at a 3D point on a sphere mesh of arbitrary radius.
 *
 * \param sphere    (in)  polygon mesh representing a sphere (any radius)
 * \param values    (in)  double[sphere->n_points]; scalar values at each vertex
 * \param pt        (in)  3D point at which to interpolate
 * \return               Interpolated scalar value
 */
double interp_point_sphere(polygons_struct *sphere, double *values, Point pt);
/**
 * \brief Interpolate scalar value at (u,v) latitude/longitude coordinates on unit sphere.
 *
 * \param sphere    (in)  polygon mesh representing the unit sphere
 * \param values    (in)  double[sphere->n_points]; scalar values at each vertex
 * \param u         (in)  first spherical coordinate (0..1 or other range)
 * \param v         (in)  second spherical coordinate (0..1 or other range)
 * \return               Interpolated scalar value at (u,v)
 */
double interp_uv_unit_sphere(polygons_struct *sphere, double *values, double u,
                             double v);
/**
 * \brief Interpolate scalar value at (u,v) latitude/longitude coordinates on arbitrary sphere.
 *
 * \param sphere    (in)  polygon mesh representing a sphere (any radius)
 * \param values    (in)  double[sphere->n_points]; scalar values at each vertex
 * \param u         (in)  first spherical coordinate
 * \param v         (in)  second spherical coordinate
 * \return               Interpolated scalar value at (u,v)
 */
double interp_uv_sphere(polygons_struct *sphere, double *values, double u,
                        double v);
/**
 * \brief Interpolate scalar value at Cartesian coordinates (x,y,z) on arbitrary sphere.
 *
 * \param sphere    (in)  polygon mesh representing a sphere (any radius)
 * \param values    (in)  double[sphere->n_points]; scalar values at each vertex
 * \param x         (in)  x-coordinate of interpolation point
 * \param y         (in)  y-coordinate of interpolation point
 * \param z         (in)  z-coordinate of interpolation point
 * \return               Interpolated scalar value at (x,y,z)
 */
double interp_xyz_sphere(polygons_struct *sphere, double *values, double x,
                         double y, double z);
/**
 * \brief Interpolate scalar value at Cartesian coordinates (x,y,z) on unit sphere.
 *
 * \param sphere    (in)  polygon mesh representing the unit sphere
 * \param values    (in)  double[sphere->n_points]; scalar values at each vertex
 * \param x         (in)  x-coordinate of interpolation point
 * \param y         (in)  y-coordinate of interpolation point
 * \param z         (in)  z-coordinate of interpolation point
 * \return               Interpolated scalar value at (x,y,z)
 */
double interp_xyz_unit_sphere(polygons_struct *sphere, double *values, double x,
                              double y, double z);

#endif
