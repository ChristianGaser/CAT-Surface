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

/**
 * \brief Public API for refine_mesh.
 *
 * This function is part of the CAT-Surface public library interface and is used by command-line tools.
 *
 * \param param (in/out) Parameter of refine_mesh.
 * \param param (in/out) Parameter of refine_mesh.
 * \param double (in/out) Parameter of refine_mesh.
 * \param param (in/out) Parameter of refine_mesh.
 * \param double (in/out) Parameter of refine_mesh.
 * \return Return value of refine_mesh.
 */
int  refine_mesh(Point **, polygons_struct *, double, polygons_struct *, double);

#endif
