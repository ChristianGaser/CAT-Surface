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
/**
 * \brief Label the clusters of positive values on a mesh.
 *
 * A cluster is a set of vertices with a value > 0 that are connected along
 * mesh edges.  The clusters are numbered by size, 1 being the largest;
 * clusters of equal size keep the order of their lowest vertex index.
 *
 * \param polygons (in)  mesh
 * \param values   (in)  double[n_points]; value per vertex
 * \param labels   (out) int[n_points]; cluster number per vertex, 0 where the
 *                       value is <= 0 or NaN
 * \return number of clusters
 */
int label_value_clusters(polygons_struct *polygons, const double values[],
                         int labels[]);
/**
 * \brief Extract the polygons of one label as a new mesh.
 *
 * \param polygons (in)  mesh
 * \param labels   (in)  int[n_points]; label per vertex
 * \param label    (in)  label to extract
 * \param out      (out) the polygons whose vertices all carry the label, with
 *                       the vertices renumbered; free with delete_polygons()
 * \return number of polygons in out
 */
int extract_label_polygons(polygons_struct *polygons, const int labels[],
                           int label, polygons_struct *out);

#endif
