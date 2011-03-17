/* Rachel Yotter - rachel.yotter@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id: CAT_Interpolate.h 175 2010-05-29 18:49:24Z raytrace $
 *
 */

#include <volume_io/internal_volume_io.h>
#include <bicpl.h>
#include <float.h>

#include "CAT_Surf.h"
#include "CAT_Blur2d.h"
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
