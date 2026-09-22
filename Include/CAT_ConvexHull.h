/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 *
 */

#ifndef _CAT_CONVEXHULL_H_
#define _CAT_CONVEXHULL_H_

extern int dbg;
extern int dbg2;

#include  <bicpl.h>

#include "CAT_SurfaceIO.h"
#include "CAT_Refine.h"
#include "CAT_Resample.h"

#define  TOLERANCE_2D   1.0e-3
#define  TOLERANCE_DISTANCE   1.0e-6

#define  POINT_USED_IN_CONVEX_HULL  1
#define  POINT_DISCARDED            2

/**
 * \brief Compute convex hull surface and optionally resample to a target sphere.
 *
 * \param polygons        (in)  input surface mesh
 * \param polygons_sphere (in)  target spherical mesh (NULL to skip resampling)
 * \return Newly allocated object list containing the convex hull surface
 */
object_struct ** surface_get_convex_hull(polygons_struct *polygons,
                                         polygons_struct *polygons_sphere);

#endif
