/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 *
 */

#ifndef _CAT_REFINE_H_
#define _CAT_REFINE_H_

#include <bicpl.h>

#ifndef MAX
#define MAX(x,y) (((x) > (y)) ? (x) : (y)) 
#endif

int  refine_mesh(Point **, polygons_struct *, double, polygons_struct *, double);

#endif