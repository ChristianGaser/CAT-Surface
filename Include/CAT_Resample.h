/* Rachel Yotter - rachel.yotter@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 *
 */

#include <bicpl.h>

#include "CAT_Surf.h"
#include "CAT_SurfaceIO.h"
#include "CAT_Complexity.h"

object_struct ** resample_tetrahedron(polygons_struct *, polygons_struct *,
                                      int, double *, double *);
object_struct ** resample_surface(polygons_struct *, polygons_struct *,
                                  int, double *, double *);
object_struct ** resample_surface_sphere(polygons_struct *, polygons_struct *);
void resample_noscale(polygons_struct *, polygons_struct *, double *, double *);
void resample_values(polygons_struct *, polygons_struct *, double *, double *);

