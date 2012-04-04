/* Rachel Yotter - rachel.yotter@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 *
 */

#include <bicpl.h>
#include <float.h>

#include "CAT_Surf.h"
#include "CAT_Smooth.h"
#include "CAT_Curvature.h"

#define  BINTREE_FACTOR   0.5
#define  NEW_COORDINATE_SYSTEM   1

double interp_point_unit_sphere(polygons_struct *, double *, Point);
double interp_point_sphere(polygons_struct *, double *, Point);
double interp_uv_unit_sphere(polygons_struct *, double *, double, double);
double interp_uv_sphere(polygons_struct *, double *, double, double);
double interp_xyz_sphere(polygons_struct *, double *, double, double, double);
double interp_xyz_unit_sphere(polygons_struct *, double *, double, double,
                              double z);
