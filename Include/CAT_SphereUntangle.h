/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Removal of the folds of a spherical map.
 *
 * Copyright Christian Gaser, University of Jena.
 */

#ifndef _CAT_SPHEREUNTANGLE_H_
#define _CAT_SPHEREUNTANGLE_H_

#include <bicpl.h>

#ifdef __cplusplus
extern "C"
{
#endif

    /**
     * \brief Remove the folded (negative-area) triangles of a spherical map.
     *
     * A triangle counts as bad when its signed height -- twice its area over
     * its longest edge -- is below 1e-6 of the radius: folded, or so thin that
     * writing the sphere in single precision could flip it.  Bad triangles
     * that share a vertex form a cluster.  Every cluster is re-embedded on a
     * disk of the mesh around it, grown by 2, 4, 8, 16 or 32 rings until it
     * works, with the vertices just outside the disk fixed: in the gnomonic
     * plane of the disk the barrier energy of Garanzha et al. (2021,
     * "Foldover-free maps in 50 lines of code") is minimised by L-BFGS while
     * its regularisation eps is lowered until no triangle is folded.  The
     * energy compares every triangle with its shape on \p surface (degenerate
     * ones regularised), so the disk is untangled without the slivers that a
     * Tutte or Laplacian re-embedding leaves where the map is compressed --
     * which is where the folds of an inflated sphere sit.  A disk is accepted
     * only if all its triangles end up good, and no vertex outside the
     * accepted disks moves.
     *
     * A map without bad triangles is left untouched, bit for bit.
     *
     * \param surface (in)     the surface the sphere maps (e.g. the central
     *                         surface), same mesh as \p sphere
     * \param sphere  (in/out) spherical map of \p surface
     * \param verbose (in)     print the clusters and what was done with them
     * \return the number of bad triangles left (0: none), or -1 if the meshes
     *         differ, are not triangulated, or memory runs out
     */
    int CAT_SphereRemoveFolds(polygons_struct *surface, polygons_struct *sphere,
                              int verbose);

#ifdef __cplusplus
}
#endif

#endif
