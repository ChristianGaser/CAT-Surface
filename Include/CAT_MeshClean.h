/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 *
 */

#ifndef _CAT_MESH_CLEAN_H_
#define _CAT_MESH_CLEAN_H_

/**
 * @file CAT_MeshClean.h
 * @brief Surface mesh cleaning and intersection removal.
 *
 * Provides functions for removing self-intersections and degeneracies
 * from triangle meshes using the MeshFix library.
 */

#include <bicpl.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Options for mesh cleaning operations.
 */
typedef struct {
    int max_iters;      /**< Maximum iterations for meshclean (default: 10) */
    int inner_loops;    /**< Inner loop iterations (default: 3) */
    int fill_holes;     /**< Whether to fill boundary holes (default: 1) */
    int verbose;        /**< Verbosity level (default: 0) */
} CAT_MeshCleanOptions;

/**
 * \brief Public API for CAT_MeshCleanOptionsInit.
 *
 * This function is part of the CAT-Surface public library interface and is used by command-line tools.
 *
 * \param opts (in/out) Parameter of CAT_MeshCleanOptionsInit.
 * \return void (no return value).
 */
void CAT_MeshCleanOptionsInit(CAT_MeshCleanOptions *opts);

/**
 * Check and fix mesh intersections and degeneracies.
 *
 * Uses MeshFix's meshclean algorithm to iteratively remove
 * self-intersections and degenerate triangles from the mesh.
 * The input surface is modified in place.
 *
 * @param polygons Input/output surface mesh. Modified in place.
 * @param opts     Options controlling the cleaning process (NULL for defaults).
 *
 * @return 0 on success (all issues fixed),
 *         1 if some issues remain,
 *        -1 on error.
 */
int CAT_SurfMeshClean(
    polygons_struct *polygons,
    const CAT_MeshCleanOptions *opts
);

/**
 * \brief Public API for CAT_SurfCountIntersections.
 *
 * This function is part of the CAT-Surface public library interface and is used by command-line tools.
 *
 * \param polygons (in/out) Parameter of CAT_SurfCountIntersections.
 * \return Return value of CAT_SurfCountIntersections.
 */
int CAT_SurfCountIntersections(polygons_struct *polygons);

/**
 * Fast self-intersection mask using MeshFix.
 *
 * Returns an int* mask of length n_items (triangles), with 1 for
 * intersecting triangles, 0 otherwise.  Caller must free() the result.
 *
 * @param polygons  Input surface mesh.
 * @param n_hits_out Output: number of intersecting triangles (can be NULL).
 * @return Allocated mask array, or NULL on error.
 */
int *find_self_intersections_meshfix(polygons_struct *polygons, int *n_hits_out);

#ifdef __cplusplus
}
#endif

#endif /* _CAT_MESH_CLEAN_H_ */
