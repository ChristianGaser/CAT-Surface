/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 */

#ifndef CAT_VOLUMEREG_H
#define CAT_VOLUMEREG_H

#include <bicpl.h>
#include "CAT_NiftiLib.h"
#include "CAT_BBReg.h"

#ifdef __cplusplus
extern "C"
{
#endif

    /**
     * \brief Robust multi-resolution volume-to-volume rigid registration.
     *
     * Implements a linearised intensity-matching criterion with Tukey biweight
     * M-estimation, analogous to FreeSurfer's mri_robust_register --cost ROB.
     * A Gaussian image pyramid is used for coarse-to-fine initialisation.
     *
     * Reference: Reuter M, Rosas HD, Fischl B (2010), NeuroImage 53:1181-1196.
     *
     * Both volumes must already be loaded as float arrays.  The transform
     * convention follows CAT_BBReg: p_out maps fixed (T1w / reference) RAS
     * coordinates to moving (EPI / BOLD) RAS coordinates, i.e. the same
     * direction used internally by CAT_BBReg_cost().
     *
     * \param fixed_vol   (in)  reference volume data (float, row-major x-fastest)
     * \param fixed_nii   (in)  NIfTI header of the reference volume
     * \param fixed_dims  (in)  [nx, ny, nz] of the reference volume
     * \param moving_vol  (in)  moving volume data (float)
     * \param moving_nii  (in)  NIfTI header of the moving volume
     * \param moving_dims (in)  [nx, ny, nz] of the moving volume
     * \param p_out       (out) estimated 6-DOF rigid parameters (T1w→EPI)
     * \param n_levels    (in)  Gaussian pyramid levels (3–4 recommended)
     * \param sat_k       (in)  Tukey saturation multiplier (4.685 = Gaussian 95%)
     * \param max_iter    (in)  IRLS iterations per level (10–20)
     * \param verbose     (in)  non-zero to print progress
     * \return            final normalised residual (lower = better alignment)
     */
    double CAT_VolumeReg_register(float *fixed_vol,
                                  nifti_image *fixed_nii,
                                  int fixed_dims[3],
                                  float *moving_vol,
                                  nifti_image *moving_nii,
                                  int moving_dims[3],
                                  CAT_RigidParams *p_out,
                                  int n_levels,
                                  double sat_k,
                                  int max_iter,
                                  int verbose);

#ifdef __cplusplus
}
#endif

#endif /* CAT_VOLUMEREG_H */
