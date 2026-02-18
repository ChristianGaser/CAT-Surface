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

/**
 * \brief Public API for separate_polygons.
 *
 * This function is part of the CAT-Surface public library interface and is used by command-line tools.
 *
 * \param param (in/out) Parameter of separate_polygons.
 * \param int (in/out) Parameter of separate_polygons.
 * \param param (in/out) Parameter of separate_polygons.
 * \return Return value of separate_polygons.
 */
int separate_polygons(polygons_struct *, int, object_struct ***);
/**
 * \brief Public API for triangulate_polygons.
 *
 * This function is part of the CAT-Surface public library interface and is used by command-line tools.
 *
 * \param param (in/out) Parameter of triangulate_polygons.
 * \param param (in/out) Parameter of triangulate_polygons.
 * \return void (no return value).
 */
void  triangulate_polygons(polygons_struct *, polygons_struct *);

#endif
