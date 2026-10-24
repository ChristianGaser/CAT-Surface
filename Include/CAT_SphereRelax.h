/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Metric relaxation of a spherical map, as FreeSurfer's mris_sphere does it.
 *
 * Copyright Christian Gaser, University of Jena.
 */

#ifndef _CAT_SPHERERELAX_H_
#define _CAT_SPHERERELAX_H_

#include <bicpl.h>

#ifdef __cplusplus
extern "C"
{
#endif

    /**
     * \brief Options of the metric relaxation of a spherical map.
     *
     * Set the defaults with CAT_SphereRelaxOptionsInit() and change what you
     * need.
     */
    typedef struct
    {
        int    n_rings;  /**< neighbourhood of the distance term in rings of
                              the mesh (default 2) */
        int    max_iter; /**< iterations per averaging scale (default 30) */
        double l_nlarea; /**< weight of FreeSurfer's nonlinear area term, which
                              keeps triangles from collapsing and pushes folded
                              ones back (default 1000) */
        double tol;      /**< relative decrease of the energy below which a
                              scale ends (default 1e-4) */
        int    verbose;  /**< print the energy of every scale */
    } CAT_SphereRelaxOpts;

    /**
     * \brief Set the default options of the metric relaxation.
     *
     * \param opts (out) options to initialise
     */
    void CAT_SphereRelaxOptionsInit(CAT_SphereRelaxOpts *opts);

    /**
     * \brief Reduce the metric distortion of a spherical map of a surface.
     *
     * Moves the vertices of \p sphere on the sphere so that the distances
     * between every vertex and its neighbours within opts->n_rings rings match
     * the distances on \p surface, scaled to the area of the sphere -- the
     * distance term of FreeSurfer's mris_sphere.  The gradient is averaged over
     * 1024, 256, 64, 16, 4, 1 and 0 neighbour iterations in turn, as in
     * FreeSurfer's MRISunfold: averaged over many neighbours it moves area
     * between whole lobes, which a local step cannot, and the inflation leaves
     * the poles of the brain compressed to a third of their area and the
     * medial wall expanded.  FreeSurfer's nonlinear area term keeps triangles
     * from collapsing and pushes folded ones back; folds left after each scale
     * are removed by a local Laplacian.
     *
     * \param surface (in)     the surface the sphere maps (e.g. the central
     *                         surface), same mesh as \p sphere
     * \param sphere  (in/out) spherical map of \p surface, centred on the
     *                         origin on return
     * \param opts    (in)     options; NULL for CAT_SphereRelaxOptionsInit()
     * \return 0 on success, -1 if the meshes differ, are not triangulated, or
     *         memory runs out
     */
    int CAT_SphereRelaxMetric(polygons_struct *surface, polygons_struct *sphere,
                              const CAT_SphereRelaxOpts *opts);

#ifdef __cplusplus
}
#endif

#endif
