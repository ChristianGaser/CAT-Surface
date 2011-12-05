/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 *
 */

#include <volume_io/internal_volume_io.h>
#include <volume_io/geometry.h>
#include <bicpl.h>
#include <float.h>

void calc_convexity(polygons_struct *, int [], int *[], double, double *);

