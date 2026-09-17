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

int intersect_poly_poly(int, int, polygons_struct *);
int intersect_triangle_triangle(int [3], int [3], polygons_struct *);
int intersect_segment_triangle(Point, Point, int [3], polygons_struct *);
/**
 * \brief Public API for find_selfintersections.
 *
 * This function is part of the CAT-Surface public library interface and is used by command-line tools.
 *
 * \param param (in/out) Parameter of find_selfintersections.
 * \param param (in/out) Parameter of find_selfintersections.
 * \param param (in/out) Parameter of find_selfintersections.
 * \param int (in/out) Parameter of find_selfintersections.
 * \return Return value of find_selfintersections.
 */
int find_selfintersections(polygons_struct *, int *, int *, int);
/**
 * \brief Public API for join_intersections.
 *
 * This function is part of the CAT-Surface public library interface and is used by command-line tools.
 *
 * \param param (in/out) Parameter of join_intersections.
 * \param param (in/out) Parameter of join_intersections.
 * \param param (in/out) Parameter of join_intersections.
 * \param param (in/out) Parameter of join_intersections.
 * \param param (in/out) Parameter of join_intersections.
 * \return Return value of join_intersections.
 */
int join_intersections(polygons_struct *, int *, int *, int *, int **);
int find_remaining_intersections(polygons_struct *, int *, int *, int *,
                                 int **);
int patch_selfintersections(polygons_struct *, polygons_struct *, int *, int *,
                            int, int *, int **);
int smooth_selfintersections(polygons_struct *, int *, int *, int, int *,
                             int **, int);
int has_selfintersections(polygons_struct *, int *, int);
int find_intersecting_defects(polygons_struct *, int *, int, int *);
void remove_intersections(polygons_struct *, int);
int remove_intersections_iter(polygons_struct *, int, int, int);

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
int *find_near_self_intersections(polygons_struct *polygons, double threshold_factor, 
                            int *n_hits_out);

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
void remove_near_intersections(polygons_struct *polygons, double threshold, int verbose);

#endif
