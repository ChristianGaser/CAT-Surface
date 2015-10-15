/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id: CAT_SurfParameters.h 333 2015-01-27 10:46:13Z gaser $
 *
 */


#include <bicpl.h>

void compute_sulcus_depth(polygons_struct *, polygons_struct *, double *);
void compute_local_sharpness(polygons_struct *, int [], int * [], double *);
void compute_convexity(polygons_struct *, int [], int *[], double, double *);
