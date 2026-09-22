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
 * \brief Refine a mesh by subdividing long edges.
 *
 * \param length_points   (in/out) points used for length measurement
 * \param polygons        (in)  input mesh
 * \param max_length      (in)  maximum allowed edge length
 * \param new_polygons    (out) refined mesh
 * \param weight_curvature (in) curvature weighting (0 disables)
 * \return Number of new polygons added
 */
int refine_mesh(Point *length_points[], polygons_struct *polygons,
                double max_length, polygons_struct *new_polygons,
                double weight_curvature);

#endif
