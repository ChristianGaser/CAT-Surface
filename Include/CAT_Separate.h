/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id: CAT_Separate.h 247 2012-03-16 11:13:59Z gaser $
 *
 */

#include <bicpl.h>

int   separate_polygons(
    polygons_struct    *polygons,
    int                desired_index,
    object_struct      **out[] );

void  triangulate_polygons(
    polygons_struct  *polygons,
    polygons_struct  *triangles );
