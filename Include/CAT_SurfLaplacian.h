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
 * @brief Pial/white surface placement via Laplacian or ADE streamline tracing.
 *
 * Solves either the Laplace equation (nabla^2 phi = 0) or the Adaptive
 * Diffusion Equation (div(f*grad(phi)) = 0) in the cortical ribbon
 * between tissue boundary isovalues, then traces streamlines from each
 * central surface vertex toward both boundaries.  Because streamlines
 * of a harmonic function cannot cross, the resulting surfaces are
 * guaranteed to be non-self-intersecting (up to discretization).
 *
 * References:
 *   Jones SE, Buchbinder BR, Aharon I (2000).  Three-dimensional mapping
 *   of cortical thickness using Laplace's equation.  Human Brain Mapping
 *   11:12-32.
 *
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
     *        Laplacian streamline tracing through the cortical ribbon.
     *
     * Solves the Laplace equation nabla^2 phi = 0 in the cortical ribbon
     * (voxels with tissue label between lim_pial and lim_white), with
     * boundary conditions phi = 1 at the pial (CSF) side and phi = 0 at
     * the white (WM) side.  For each vertex of the central surface,
     * streamlines are traced along +grad(phi) (toward pial) and
     * -grad(phi) (toward white).  The endpoints define the output
     * surface positions.
     *
     * \param central          (in)  central surface mesh (must have normals)
     * \param labels           (in)  tissue label / intensity volume
     * \param nii_ptr          (in)  NIfTI header for coordinate transforms
     * \param lim_pial         (in)  intensity threshold for the pial boundary
     * \param lim_white        (in)  intensity threshold for the white boundary
     * \param thickness_values (in)  per-vertex cortical thickness (for max
     *                               travel distance and fallback placement)
     * \param pial_out         (out) output pial surface (copy of central,
     *                               vertices displaced)
     * \param white_out        (out) output white surface (copy of central,
     *                               vertices displaced)
     * \param verbose          (in)  enable progress output
     * \return 0 on success, non-zero on error
     */
    int surf_laplacian_pial_white(
        polygons_struct *central,
        float *labels,
        nifti_image *nii_ptr,
        float lim_pial,
        float lim_white,
        const double *thickness_values,
        polygons_struct *pial_out,
        polygons_struct *white_out,
        int verbose);

    /**
     * \brief Compute pial and white surfaces from a central surface using
     *        Adaptive Diffusion Equation (ADE) streamline tracing.
     *
     * Like \ref surf_laplacian_pial_white, but solves the ADE
     * div(f(x) grad(phi)) = 0 where f(x) is the gray-matter fraction
     * derived from the tissue label volume.  This makes the streamlines
     * sensitive to partial-volume content, producing more accurate
     * placement in regions with thin cortex or strong partial-volume
     * effects.
     *
     * Parameters are identical to \ref surf_laplacian_pial_white.
     *
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
