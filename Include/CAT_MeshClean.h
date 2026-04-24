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
 * @brief Surface mesh cleaning and intersection removal (non-Windows only).
 *
 * Provides functions for removing self-intersections and degeneracies
 * from triangle meshes using the MeshFix library.
 *
 * MeshFix requires C++ and is not supported on Windows.  This header
 * provides no declarations when compiled on Windows; callers must guard
 * any inclusion and all call sites with
 * #if !defined(_WIN32) && !defined(_WIN64).
 */

#if !defined(_WIN32) && !defined(_WIN64)

#include <bicpl.h>

#ifdef __cplusplus
extern "C"
{
#endif

    /**
     * Options for mesh cleaning operations.
     */
    typedef struct
    {
        int max_iters;   /**< Maximum iterations for meshclean (default: 10) */
        int inner_loops; /**< Inner loop iterations (default: 3) */
        int fill_holes;  /**< Whether to fill boundary holes (default: 1) */
        int verbose;     /**< Verbosity level (default: 0) */
    } CAT_MeshCleanOptions;

    /**
     * \brief Initialise a CAT_MeshCleanOptions struct with default values.
     *
     * \param opts (in/out) Options struct to initialise.
     */
    void CAT_MeshCleanOptionsInit(CAT_MeshCleanOptions *opts);

    /**
     * \brief Clean a surface mesh using MeshFix.
     *
     * Uses MeshFix's meshclean algorithm to iteratively remove
     * self-intersections and degenerate triangles from the mesh.
     * The input surface is modified in place.
     *
     * \param polygons (in/out) Surface mesh to clean.
     * \param opts     (in)     Cleaning options, or NULL for defaults.
     * \return 0 on success (all issues fixed), 1 if some issues remain, -1 on error.
     */
    int CAT_SurfMeshClean(
        polygons_struct *polygons,
        const CAT_MeshCleanOptions *opts);

    /**
     * \brief Count self-intersecting triangle pairs using MeshFix.
     *
     * \param polygons (in) Surface mesh to inspect.
     * \return Number of intersecting triangle pairs, or -1 on error.
     */
    int CAT_SurfCountIntersections(polygons_struct *polygons);

    /**
     * \brief Build a per-triangle self-intersection mask using MeshFix.
     *
     * Returns an int* mask of length polygons->n_items, with 1 for
     * intersecting triangles, 0 otherwise.  Caller must free() the result.
     *
     * \param polygons    (in)  Surface mesh to inspect.
     * \param n_hits_out  (out) Number of intersecting triangles (may be NULL).
     * \return Allocated mask array, or NULL on error.
     */
    int *find_self_intersections_meshfix(polygons_struct *polygons, int *n_hits_out);

#ifdef __cplusplus
}
#endif

#endif /* !_WIN32 */

#endif /* _CAT_MESH_CLEAN_H_ */
