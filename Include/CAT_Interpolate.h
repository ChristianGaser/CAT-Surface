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


#endif
