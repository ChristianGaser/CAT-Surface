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
     * \brief Gaussian image pyramid of one volume, the input of the registration.
     *
     * Built once per image by CAT_VolumeReg_pyramid(), so a reference that many
     * volumes are registered to (motion correction of a time series) is blurred
     * and subsampled once rather than once per volume.  Opaque; release with
     * CAT_VolumeReg_pyramid_free().
     */
    typedef struct CAT_VolPyramid CAT_VolPyramid;

    /**
     * \brief Build the Gaussian image pyramid of one volume.
     *
     * Level 0 is a copy of the volume; every further level blurs the previous
     * one by 0.5 voxels and keeps every second voxel along each axis.  The
     * volume is not referenced after the call.
     *
     * \param vol      (in) volume data (float, x fastest)
     * \param nii      (in) header whose sto_xyz maps voxels to world coordinates
     *                      (as read_nifti_float() returns it)
     * \param dims     (in) [nx, ny, nz]
     * \param n_levels (in) number of levels, clamped to 1..6
     * \return the pyramid, or NULL if out of memory
     */
    CAT_VolPyramid *CAT_VolumeReg_pyramid(const float *vol,
                                          const nifti_image *nii,
                                          const int dims[3],
                                          int n_levels);

    /**
     * \brief Release a pyramid from CAT_VolumeReg_pyramid().
     *
     * \param pyr (in/out) pyramid to free; NULL is ignored
     */
    void CAT_VolumeReg_pyramid_free(CAT_VolPyramid *pyr);

    /**
     * \brief Options of the robust rigid registration.
     *
     * Set the defaults with CAT_VolumeRegOptionsInit() and change what you need.
     */
    typedef struct
    {
        int    n_levels; /**< Gaussian pyramid levels (default 4); used by
                              CAT_VolumeReg_register(), prepared pyramids
                              bring their own */
        double sat_k;    /**< Tukey saturation multiplier (default 4.685, the
                              95% efficiency point for Gaussian residuals) */
        int    max_iter; /**< IRLS iterations per pyramid level (default 20) */
        int    iscale;   /**< non-zero: estimate a global intensity scale
                              (default 1) */
        int    verbose;  /**< 0 quiet (default), 1 level summaries, >1 every
                              iteration */
    } CAT_VolumeRegOpts;

    /**
     * \brief Set the default options of the robust rigid registration.
     *
     * The defaults are the source of truth for every front-end: 4 pyramid
     * levels, a Tukey multiplier of 4.685, 20 iterations per level, and a
     * global intensity scale estimated along with the motion.
     *
     * \param opts (out) options to initialise
     */
    void CAT_VolumeRegOptionsInit(CAT_VolumeRegOpts *opts);

    /**
     * \brief Robust rigid registration of two prepared pyramids.
     *
     * The registration of CAT_VolumeReg_register(), run coarse to fine over the
     * levels both pyramids have.
     *
     * \param fixed  (in)     pyramid of the reference volume
     * \param moving (in)     pyramid of the moving volume
     * \param p      (in/out) initial and estimated rigid parameters (fixed→moving)
     * \param scale  (in/out) intensity scale s of the model s * moving = fixed;
     *                        in: a positive start value, or <= 0 to fit one by
     *                        least squares at the coarsest level; out: the
     *                        estimate, 1 without opts->iscale.  May be NULL.
     * \param opts   (in)     options; NULL for the defaults
     * \return final normalised residual (lower = better alignment)
     */
    double CAT_VolumeReg_register_pyramids(const CAT_VolPyramid *fixed,
                                           const CAT_VolPyramid *moving,
                                           CAT_RigidParams *p,
                                           double *scale,
                                           const CAT_VolumeRegOpts *opts);

    /**
     * \brief Robust multi-resolution volume-to-volume rigid registration.
     *
     * Implements a linearised intensity-matching criterion with Tukey biweight
     * M-estimation, analogous to FreeSurfer's mri_robust_register --cost ROB.
     * A Gaussian image pyramid is used for coarse-to-fine initialisation.
     *
     * By default a global intensity scale s is estimated along with the motion
     * (the residual is s * moving - fixed, like mri_robust_register --iscale):
     * a reference brighter than the series, such as the median of its
     * non-steady-state volumes or a single-band reference, otherwise turns
     * most of the brain into Tukey outliers.  Without it the registration of
     * such a series was off by 1.4 mm on average, with it by 0.04 mm.
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
     * \param p_out       (in/out) initial and estimated 6-DOF rigid parameters (T1w→EPI)
     * \param scale       (in/out) intensity scale, see
     *                             CAT_VolumeReg_register_pyramids(); may be NULL
     * \param opts        (in)  options; NULL for the defaults
     * \return            final normalised residual (lower = better alignment)
     */
    double CAT_VolumeReg_register(float *fixed_vol,
                                  nifti_image *fixed_nii,
                                  int fixed_dims[3],
                                  float *moving_vol,
                                  nifti_image *moving_nii,
                                  int moving_dims[3],
                                  CAT_RigidParams *p_out,
                                  double *scale,
                                  const CAT_VolumeRegOpts *opts);

    /**
     * \brief Cross-modal rigid registration using Normalised Mutual Information.
     *
     * Equivalent to FreeSurfer's mri_coreg, the initialisation step used by
     * fMRIPrep before bbregister.  A Gaussian image pyramid is used for
     * coarse-to-fine registration; at each level NMI is computed from a joint
     * intensity histogram (bilinear partial-volume interpolation) and maximised
     * with Powell's direction-set method + Brent line search.
     *
     * Appropriate for cross-modal pairs such as BOLD EPI ↔ T1w.  For
     * same-modality registration use CAT_VolumeReg_register() instead.
     *
     * \param fixed_vol   (in)  reference volume data (float, x-fastest) — T1w
     * \param fixed_nii   (in)  NIfTI header of the reference volume
     * \param fixed_dims  (in)  [nx, ny, nz] of the reference volume
     * \param moving_vol  (in)  moving volume data (float) — mean EPI/BOLD
     * \param moving_nii  (in)  NIfTI header of the moving volume
     * \param moving_dims (in)  [nx, ny, nz] of the moving volume
     * \param p_out       (out) estimated 6-DOF rigid parameters (fixed→moving)
     * \param n_levels    (in)  Gaussian pyramid levels (3–4 recommended)
     * \param n_bins      (in)  histogram bins per axis (64 recommended)
     * \param max_iter    (in)  Powell iterations per level (20–50)
     * \param verbose     (in)  ≥1 level summaries; ≥2 per-iteration output
     * \return            final NMI value (higher = better alignment)
     */
    double CAT_VolumeReg_register_NMI(float *fixed_vol,
                                      nifti_image *fixed_nii,
                                      int fixed_dims[3],
                                      float *moving_vol,
                                      nifti_image *moving_nii,
                                      int moving_dims[3],
                                      CAT_RigidParams *p_out,
                                      int n_levels,
                                      int n_bins,
                                      int max_iter,
                                      int verbose);

#ifdef __cplusplus
}
#endif

#endif /* CAT_VOLUMEREG_H */
