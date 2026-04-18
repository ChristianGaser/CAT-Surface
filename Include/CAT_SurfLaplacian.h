/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 *
 */

#ifndef _CAT_SURF_LAPLACIAN_H_
#define _CAT_SURF_LAPLACIAN_H_

/**
 * @file CAT_SurfLaplacian.h
 * @brief Pial/white surface placement via ADE streamline tracing.
 *
 * Solves the Adaptive Diffusion Equation (div(f*grad(phi)) = 0) in the
 * cortical ribbon between tissue boundary isovalues, then traces
 * streamlines from each central surface vertex toward both boundaries.
 * Because streamlines of a harmonic function cannot cross, the resulting
 * surfaces are guaranteed to be non-self-intersecting (up to
 * discretization).
 *
 * Reference:
 *   Joshi AA et al. (2025).  Robust Cortical Thickness Estimation in the
 *   Presence of Partial Volumes using Adaptive Diffusion Equation.
 *   J Neurosci Methods 423:110552.
 */

#include <bicpl.h>
#include "CAT_NiftiLib.h"

#ifdef __cplusplus
extern "C"
{
#endif

    /**
     * \brief Compute pial and white surfaces from a central surface using
     *        Adaptive Diffusion Equation (ADE) streamline tracing.
     *
     * Solves the ADE  div(f(x) grad(phi)) = 0  where f(x) is the
     * gray-matter fraction derived from the tissue label volume.  This
     * makes the streamlines sensitive to partial-volume content, producing
     * more accurate placement in regions with thin cortex or strong
     * partial-volume effects.
     *
     * The ADE is solved on a wide ribbon from CSF to WM so the phi field
     * has enough room for smooth gradients.  The target isovalues
     * lim_pial and lim_white are mapped to phi stop values, so changing
     * them moves the resulting surfaces.
     *
     * \param central          (in)  central surface mesh (must have normals)
     * \param labels           (in)  tissue label / intensity volume
     * \param nii_ptr          (in)  NIfTI header for coordinate transforms
     * \param lim_pial         (in)  pial boundary isovalue (e.g. CGM=1.5)
     * \param lim_white        (in)  white boundary isovalue (e.g. GWM=2.5)
     * \param thickness_values (in)  per-vertex cortical thickness (for max
     *                               travel distance and fallback placement)
     * \param pial_out         (out) output pial surface (copy of central,
     *                               vertices displaced)
     * \param white_out        (out) output white surface (copy of central,
     *                               vertices displaced)
     * \param verbose          (in)  enable progress output
     * \return 0 on success, non-zero on error
     */
    int surf_ade_pial_white(
        polygons_struct *central,
        float *labels,
        nifti_image *nii_ptr,
        float lim_pial,
        float lim_white,
        const double *thickness_values,
        polygons_struct *pial_out,
        polygons_struct *white_out,
        int verbose);

#ifdef __cplusplus
}
#endif

#endif /* _CAT_SURF_LAPLACIAN_H_ */
