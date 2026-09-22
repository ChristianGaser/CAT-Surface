/* Rachel Yotter - rachel.yotter@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 *
 */

#ifndef _CAT_INTERSECT_H_
#define _CAT_INTERSECT_H_

#include <bicpl.h>

#define PINF  1.7976931348623157e+308 /* for doubles */
#define NINF -1.7976931348623157e+308 /* for doubles */

#define BINTREE_FACTOR 0.5

/**
 * \brief Test geometric intersection between two triangles in 3D space.
 *
 * \param pidx0    (in) int[3]; vertex indices of first triangle
 * \param pidx1    (in) int[3]; vertex indices of second triangle
 * \param surface (in) mesh containing vertex coordinates
 * \return 1 if triangles intersect, 0 if disjoint
 */
int intersect_triangle_triangle(int pidx0[3], int pidx1[3],
                                polygons_struct *surface);
/**
 * \brief Test if a line segment intersects a triangle in 3D space.
 *
 * \param p0 (in) first endpoint of line segment
 * \param p1 (in) second endpoint of line segment
 * \param tpidx    (in) int[3]; vertex indices of triangle
 * \param surface (in) mesh with vertex coordinates
 * \return -1 degenerate triangle, 0 no intersection, 1 unique intersection, 2 coplanar
 */
int intersect_segment_triangle(Point p0, Point p1, int tpidx[3],
                               polygons_struct *surface);
/**
 * \brief Find self-intersecting triangles of a mesh.
 *
 * \param polygons    (in)     triangle mesh
 * \param defects     (out)    per-vertex labels (0 = no intersection)
 * \param polydefects (in/out) per-triangle labels; with init == 0, triangles
 *                                 with a negative label are skipped
 * \param init        (in)     non-zero to clear polydefects first
 * \return number of intersecting triangle pairs
 */
int find_selfintersections(polygons_struct *polygons, int *defects,
                           int *polydefects, int init);
/**
 * \brief Consolidate spatially-connected self-intersection regions into components.
 *
 * \param surface (in) mesh structure with neighborhood info
 * \param defects (in/out) per-vertex defect labels (remapped and consolidated)
 * \param polydefects (in/out) per-polygon defect labels (updated from consolidated vertex labels)
 * \param n_neighbours (in) neighbor counts per vertex
 * \param neighbours (in) neighbor lists per vertex
 * \return number of consolidated defect components
 */
int join_intersections(polygons_struct *surface, int *defects, int *polydefects,
                       int *n_neighbours, int **neighbours);
/**
 * \brief Re-test marked self-intersections to identify which ones persist after correction.
 *
 * \param surface (in) mesh after partial correction attempt
 * \param defects (in/out) per-vertex labels (recomputed from remaining intersections)
 * \param polydefects (in/out) per-polygon labels (cleared for resolved, updated for persistent)
 * \param n_neighbours (in) vertex connectivity
 * \param neighbours (in) neighbor lists
 * \return count of remaining unresolved intersection groups
 */
int find_remaining_intersections(polygons_struct *surface, int *defects,
                                 int *polydefects, int *n_neighbours,
                                 int **neighbours);
/**
 * \brief Replace intersection regions with corrected coordinates from reference patch surface.
 *
 * \param surface (in/out) mesh to repair
 * \param patch (in) reference surface with corrected geometry
 * \param defects (in/out) per-vertex defect labels
 * \param polydefects (in/out) per-polygon defect labels
 * \param n_defects (in) number of distinct defects (for context)
 * \param n_neighbours (in) vertex connectivity
 * \param neighbours (in) neighbor lists
 * \return remaining self-intersections after patching
 */
int patch_selfintersections(polygons_struct *surface, polygons_struct *patch,
                            int *defects, int *polydefects, int n_defects,
                            int *n_neighbours, int **neighbours);
/**
 * \brief Iteratively smooth defect regions to resolve self-intersections via Laplacian relaxation.
 *
 * \param surface (in/out) mesh modified by iterative smoothing
 * \param defects (in/out) per-vertex defect labels
 * \param polydefects (in/out) per-polygon defect labels
 * \param n_defects (in) number of distinct defects to track
 * \param n_neighbours (in) vertex neighbor counts
 * \param neighbours (in/out) neighbor lists (may be updated)
 * \param maxiter (in) maximum smoothing iterations (typical 200-500)
 * \return number of defects successfully repaired
 */
int smooth_selfintersections(polygons_struct *surface, int *defects,
                             int *polydefects, int n_defects, int *n_neighbours,
                             int **neighbours, int maxiter);
/**
 * \brief Test all defect regions for remaining self-intersections in a single pass.
 *
 * \param polygons (in) mesh to check
 * \param polydefects (in) per-polygon defect labels
 * \param n_defects (in) largest defect ID in use
 * \param siflags (out) array of n_defects+1 entries; siflags[d] is set to 1 if defect
 *                          d still self-intersects, 0 otherwise (index 0 is unused)
 * \return number of defects that still self-intersect
 */
int find_intersecting_defects(polygons_struct *polygons, int *polydefects,
                              int n_defects, int *siflags);
/**
 * \brief Remove all self-intersections from a mesh with explicit iteration limits.
 *
 * \param polygons (in/out) mesh to repair
 * \param max_passes (in) maximum number of detect/smooth passes (default 10)
 * \param maxiter (in) maximum smoothing iterations per pass (default 50)
 * \param verbose (in) 1 for progress output; 0 for silent
 * \return number of self-intersecting defect regions that remain (0 = fully repaired)
 */
int remove_intersections_iter(polygons_struct *polygons, int max_passes,
                              int maxiter, int verbose);

/** Share of the way back to the reference per retreat step. */
#define CAT_RETREAT_FRACTION 0.25
/** Rings of neighbours moved together with a remaining defect. */
#define CAT_RETREAT_RINGS 2
/** Maximum number of retreat steps. */
#define CAT_RETREAT_STEPS 16

/**
 * \brief Remove self-intersections, retreating stubborn defects towards a reference.
 *
 * Runs remove_intersections_iter() and, where defects survive it, moves their
 * vertices and CAT_RETREAT_RINGS rings of neighbours by CAT_RETREAT_FRACTION of
 * the way back to the reference positions before repairing again, for at most
 * CAT_RETREAT_STEPS steps, stopping early when two steps bring no progress.
 * Local smoothing cannot separate two sheets that were driven through each
 * other, e.g. the two sides of a thin gyral blade; the surface a deformation
 * started from gives the way back.
 *
 * \param polygons   (in/out) mesh to repair
 * \param reference  (in)     reference positions, one per vertex of polygons
 *                            (same topology, e.g. the start of the
 *                            deformation); NULL makes this identical to
 *                            remove_intersections_iter()
 * \param max_passes (in)     detect/smooth passes of each repair (default 10)
 * \param maxiter    (in)     smoothing iterations per pass (default 50)
 * \param verbose    (in)     1 for progress output; 0 for silent
 * \return number of self-intersecting defect regions that remain (0 = fully repaired)
 */
int remove_intersections_ref(polygons_struct *polygons, const Point *reference,
                             int max_passes, int maxiter, int verbose);
/**
 * \brief Find vertices closer to a non-adjacent vertex than a distance threshold.
 *
 * \param polygons         (in)  source 3D polygonal mesh (normals are not used)
 * \param threshold_factor (in)  multiplier for average edge length to define search radius
 * \param n_hits_out       (out) number of flagged vertices; may be NULL
 * \return Allocated array of flags (length = n_points, 1 = near hit), caller must free
 */
int *find_near_self_intersections(polygons_struct *polygons,
                                  double threshold_factor, int *n_hits_out);

/**
 * \brief Find vertices close to a facing sheet of the same mesh.
 *
 * Distance test of find_near_self_intersections() restricted to vertices with
 * opposing normals (n_i . n_j < -min_opposition), as across a sulcus or a thin
 * blade, so that 2-ring neighbours of the same sheet on irregular meshes are not
 * reported.
 *
 * \param polygons         (in)  mesh with current normals
 * \param threshold_factor (in)  search radius as multiple of the mean edge length
 * \param min_opposition   (in)  required opposition of the normals (e.g. 0.3)
 * \param n_hits_out       (out) number of flagged vertices; may be NULL
 * \return allocated flag array (length n_points), caller must free
 */
int *find_near_facing_intersections(polygons_struct *polygons, double threshold_factor,
                                    double min_opposition, int *n_hits_out);
/**
 * \brief Remove near-intersecting vertices by iterative vertex repositioning.
 *
 * \param polygons   (in)  source 3D polygonal mesh to be modified in-place
 * \param threshold  (in)  distance threshold for near-intersection detection (typically 0.05-0.20 times edge length)
 * \param verbose    (in)  1 to print progress messages to stdout, 0 for silent operation
 */
void remove_near_intersections(polygons_struct *polygons, double threshold,
                               int verbose);

#endif
