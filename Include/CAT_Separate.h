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
 * \brief Separate a mesh into connected components.
 *
 * \param polygons      (in)  input mesh
 * \param desired_index (in)  component index to extract, or -1 for all
 * \param out           (out) array of output objects (allocated)
 * \return Number of output objects
 */
int separate_polygons(polygons_struct *polygons, int desired_index,
                      object_struct **out[]);
/**
 * \brief Triangulate all polygons in a mesh.
 *
 * \param polygons  (in)  input mesh with polygons
 * \param triangles (out) triangulated mesh
 */
void triangulate_polygons(polygons_struct *polygons,
                          polygons_struct *triangles);

#endif
