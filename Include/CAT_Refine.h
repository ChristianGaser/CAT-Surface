/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 *
 */

#include <bicpl.h>

#ifndef MAX
#define MAX(x,y) (((x) > (y)) ? (x) : (y)) 
#endif

int  refine_mesh(Point **, polygons_struct *, double, polygons_struct *, double);
