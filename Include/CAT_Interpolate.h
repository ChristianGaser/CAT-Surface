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
 * \brief Public API for interp_point_unit_sphere.
 *
 * This function is part of the CAT-Surface public library interface and is used by command-line tools.
 *
 * \param param (in/out) Parameter of interp_point_unit_sphere.
 * \param param (in/out) Parameter of interp_point_unit_sphere.
 * \param Point (in/out) Parameter of interp_point_unit_sphere.
 * \return Return value of interp_point_unit_sphere.
 */
double interp_point_unit_sphere(polygons_struct *, double *, Point);
double interp_point_sphere(polygons_struct *, double *, Point);
double interp_uv_unit_sphere(polygons_struct *, double *, double, double);
double interp_uv_sphere(polygons_struct *, double *, double, double);
double interp_xyz_sphere(polygons_struct *, double *, double, double, double);
double interp_xyz_unit_sphere(polygons_struct *, double *, double, double, double z);

#endif
