/* Rachel Yotter - rachel.yotter@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 *
 */

#ifndef _CAT_DEFECT_H_
#define _CAT_DEFECT_H_

#include <bicpl.h>

#define BINTREE_FACTOR 0.5

#define HOLE 1
#define HANDLE 2

/**
 * \brief Calculate average normal direction vector of a topological defect patch.
 *
 * \param surface (in) input mesh
 * \param defects (in) per-vertex defect labels (0=no defect, >0=defect ID)
 * \param defect  (in) defect ID to analyze
 * \return Vector pointing in average direction of defect normals
 */
Vector defect_direction(polygons_struct *surface, int *defects, int defect);
/**
 * \brief Compute Euler characteristic of a defect patch within mesh topology.
 *
 * \param surface     (in) input mesh
 * \param defects     (in) per-vertex defect labels
 * \param polydefects (in/out) per-polygon defect labels (allocated if NULL)
 * \param defect      (in) defect ID to analyze
 * \param n_neighbours (in) per-vertex neighbor counts
 * \param neighbours  (in) per-vertex neighbor lists
 * \return Euler characteristic (typically 0, 2, -2, etc.)
 */
int defect_euler(polygons_struct *surface, int *defects, int *polydefects,
                 int defect, int *n_neighbours, int **neighbours);
/**
 * \brief Check if a triangle polygon lies on the boundary of a defect patch.
 *
 * \param surface      (in) input mesh
 * \param defects      (in) per-vertex defect labels
 * \param polydefects  (in) per-polygon defect labels
 * \param n_neighbours (in) per-vertex neighbor counts
 * \param neighbours   (in) per-vertex neighbor lists
 * \param polygon      (in) polygon index to test
 * \return 1 if polygon is on defect boundary; 0 otherwise
 */
int isedge(polygons_struct *surface, int *defects, int *polydefects,
           int *n_neighbours, int **neighbours, int polygon);
/**
 * \brief Detect topological defects (holes and handles) in a cortical surface mesh.
 *
 * \param surface      (in)  input mesh to analyze
 * \param sphere       (in)  reference sphere for validation
 * \param defects      (out) allocated per-vertex defect labels
 * \param n_neighbours (in)  per-vertex neighbor counts
 * \param neighbours   (in)  per-vertex neighbor lists
 * \return number of defects found
 */
int find_topological_defects(polygons_struct *surface, polygons_struct *sphere,
                             int *defects, int *n_neighbours, int **neighbours);
/**
 * \brief Find surface artifacts as regions far from a reference surface.
 *
 * \param surface      (in)  surface to examine
 * \param sph          (in)  reference surface, e.g. a smoothed copy of surface
 * \param artifacts    (out) per-vertex artifact labels (0 = none, 1..n)
 * \param n_neighbours (in)  per-vertex neighbour counts
 * \param neighbours   (in)  per-vertex neighbour lists
 * \param dist         (in)  distance above which a vertex is an artifact candidate
 * \return number of artifacts found
 */
int find_artifacts(polygons_struct *surface, polygons_struct *sph,
                   int *artifacts, int *n_neighbours, int **neighbours,
                   double dist);
/**
 * \brief Expand defect region by propagating labels to neighboring vertices.
 *
 * \param surface      (in) input mesh
 * \param defects      (in/out) per-vertex defect labels (modified in place)
 * \param polydefects  (in/out) per-polygon defect labels (updated)
 * \param defect       (in) defect ID to expand (0 for all)
 * \param level        (in) number of dilation iterations
 * \param n_neighbours (in) per-vertex neighbor counts
 * \param neighbours   (in) per-vertex neighbor lists
 */
void expand_defects(polygons_struct *surface, int *defects, int *polydefects,
                    int defect, int level, int *n_neighbours, int **neighbours);
/**
 * \brief Update per-vertex defect labels from per-polygon (triangle) defect labels.
 *
 * \param surface    (in)  input mesh
 * \param polydefects (in) per-polygon defect labels
 * \param defects    (out) per-vertex defect labels (overwritten)
 */
void update_defects(polygons_struct *surface, int *polydefects, int *defects);
/**
 * \brief Update per-polygon (triangle) defect labels from per-vertex defect labels.
 *
 * \param surface     (in) input mesh
 * \param defects     (in) per-vertex defect labels
 * \param polydefects (out) per-polygon defect labels (updated)
 */
void update_polydefects(polygons_struct *surface, int *defects,
                        int *polydefects);
/**
 * \brief Calculate the centroid (center point) of a defect patch.
 *
 * \param surface (in) input mesh
 * \param defects (in) per-vertex defect labels
 * \param defect  (in) defect ID to analyze
 * \return 3D point at center of mass of defect vertices
 */
Point get_defect_center(polygons_struct *surface, int *defects, int defect);
/**
 * \brief Calculate the relative size of each defect patch as fraction of total surface.
 *
 * \param surface       (in) input mesh
 * \param defects       (in) per-vertex defect labels
 * \param n_defects     (in) total number of defects
 * \param defect_size   (out) per-vertex size values (normalized 0-1)
 */
void get_defect_size(polygons_struct *surface, int *defects, int n_defects,
                     double *defect_size);
/**
 * \brief Bisect topological defects by sulcal depth to separate holes and handles.
 *
 * \param surface      (in)  original mesh
 * \param sphere       (in)  spherical reference mesh
 * \param defects      (in)  per-vertex defect labels and IDs
 * \param n_defects    (in)  number of distinct defects
 * \param holes        (in/out) array identifying hole vs handle classification
 * \param bisected     (out) per-vertex marking indicating bisection result
 * \param detect_euler (in)  flag to compute Euler characteristic for classification
 */
void bisect_defects(polygons_struct *surface, polygons_struct *sphere,
                    int *defects, int n_defects, int *holes, int *bisected,
                    int detect_euler);
/**
 * \brief Remap topological defects from one spherical surface to another reference map.
 *
 * \param sphere           (in)  source spherical surface with defects
 * \param defects          (in)  per-vertex defect labels on source
 * \param polydefects      (in)  per-polygon defect labels on source
 * \param remap            (in)  target spherical reference surface (bintree built if missing)
 * \param remap_defects    (out) per-vertex defect labels on target
 * \param remap_polydefects (out) per-polygon defect labels on target
 */
void remap_defect(polygons_struct *sphere, int *defects, int *polydefects,
                  polygons_struct *remap, int *remap_defects,
                  int *remap_polydefects);
/**
 * \brief Inflate surface while preserving topology defects for visualization and analysis.
 *
 * \param polygons (in/out) surface mesh modified in place by inflation and smoothing operations
 */
void inflate_surface_with_topology_defects(polygons_struct *polygons);

#endif
