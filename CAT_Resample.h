/* Rachel Yotter - rachel.yotter@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id: CAT_Resample.h 139 2009-10-02 15:39:28Z raytrace $
 *
 */

#include <volume_io/internal_volume_io.h>
#include <bicpl.h>

#include "Cat_Surf.h"
#include "CAT_SurfaceIO.h"
#include "CAT_Complexity.h"

object_struct ** resample_tetrahedron(polygons_struct *, polygons_struct *,
                                      int, double *, double *);
void resample_noscale(polygons_struct *, polygons_struct *, double *, double *);
void resample_values(polygons_struct *, polygons_struct *, double *, double *);
