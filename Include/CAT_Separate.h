/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 *
 */

#ifndef _CAT_SEPARATE_H_
#define _CAT_SEPARATE_H_

#include <bicpl.h>

int separate_polygons(polygons_struct *, int, object_struct ***);
void  triangulate_polygons(polygons_struct *, polygons_struct *);

#endif