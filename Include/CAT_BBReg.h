/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 */

#ifndef CAT_BBREG_H
#define CAT_BBREG_H

#include <bicpl.h>
#include "CAT_NiftiLib.h"

#ifdef __cplusplus
extern "C"
{
#endif

    /* -----------------------------------------------------------------------
     * Rigid-body transform: 6 DOF (tx, ty, tz, rx, ry, rz in radians).
     * The convention follows FreeSurfer bbregister:
     *   T_full = T_trans * R_z * R_y * R_x
     * ----------------------------------------------------------------------- */
    typedef struct
    {
        double tx, ty, tz; /* translation (mm) */
        double rx, ry, rz; /* rotation angles (radians) */
    } CAT_RigidParams;

    /**
     * \brief Convert 6-DOF rigid parameters to a 4×4 homogeneous matrix (row-major).
     *
     * The resulting matrix maps points expressed in the *moving* volume's RAS space
     * to the *fixed* surface's RAS space:
     *
     *   M = T_trans * Rz * Ry * Rx
     *
     * where R* are elementary rotation matrices and T_trans is a pure translation.
     * The 16 elements are stored row-major: m[row*4 + col].
     *
     * \param p  (in)  rigid parameters (tx, ty, tz mm; rx, ry, rz radians)
     * \param m  (out) 4×4 matrix, row-major, 16 doubles
     */
    void CAT_BBReg_params_to_matrix(const CAT_RigidParams *p, double m[16]);

    /**
     * \brief Apply a 4×4 RAS transform to every vertex of a surface (in-place).
     *
     * Transforms each vertex position of \p surface by the homogeneous matrix \p m
     * (row-major, 16 doubles).  Normals are NOT updated; call
     * compute_polygon_normals() afterwards if needed.
     *
     * \param surface (in/out) surface mesh whose vertices are modified
     * \param m       (in)     4×4 transform, row-major
     */
    void CAT_BBReg_apply_matrix(polygons_struct *surface, const double m[16]);

    /**
     * \brief Invert a 4×4 rigid homogeneous matrix.
     *
     * Because the matrix is rigid (R|t), the inverse is computed analytically as
     *   R^T | -R^T t
     * without LU decomposition.
     *
     * \param m    (in)  4×4 row-major rigid matrix
     * \param inv  (out) 4×4 row-major inverse
     */
    void CAT_BBReg_invert_matrix(const double m[16], double inv[16]);

    /**
     * \brief Compute the BBR cost for a given set of rigid parameters.
     *
     * For each surface vertex the volume is sampled at two points along the
     * outward surface normal:
     *   - WM side : vertex - wm_dist  * normal  (inside white matter)
     *   - GM side : vertex + gm_dist  * normal  (inside grey matter)
     *
     * The contrast at vertex i is defined as
     *   c_i = (I_WM - I_GM) / (0.5 * (I_WM + I_GM) + eps)
     *
     * and the cost is the mean of a robust saturation function applied to c_i:
     *   cost = mean_i { q(c_i) }
     *
     * where q(c) = 1 - tanh^2(slope * c).  For T1 images WM > GM so c_i should
     * be positive at correctly registered vertices; cost is therefore minimised
     * when the surface aligns with the WM/GM boundary.  Pass \p invert_contrast
     * non-zero for T2/BOLD data (WM < GM).
     *
     * \param p               (in)  current 6-DOF rigid parameters
     * \param surface         (in)  WM surface in surface-native RAS space
     * \param vol             (in)  floating-point volume data (moving image)
     * \param nii_ptr         (in)  NIfTI header of the moving volume
     * \param dims            (in)  volume dimensions [nx, ny, nz]
     * \param wm_dist         (in)  sampling distance inside WM along inward normal (mm, positive)
     * \param gm_dist         (in)  sampling distance inside GM along outward normal (mm, positive)
     * \param slope           (in)  saturation slope of the robust cost function
     * \param invert_contrast (in)  non-zero for T2/BOLD (negates contrast sign)
     * \return mean BBR cost (lower = better alignment)
     */
    double CAT_BBReg_cost(const CAT_RigidParams *p,
                          const polygons_struct *surface,
                          float *vol,
                          nifti_image *nii_ptr,
                          int dims[3],
                          double wm_dist,
                          double gm_dist,
                          double slope,
                          int invert_contrast);

    /**
     * \brief Optimise the BBR cost function using a two-stage strategy.
     *
     * Stage 1 — brute-force grid search over a coarse ±range neighbourhood of
     *            \p p_init, evaluating the cost at each grid point.  The best
     *            parameter set becomes the starting point for Stage 2.
     *
     * Stage 2 — Powell's direction-set method refines the result from Stage 1
     *            until the fractional cost improvement falls below \p tol or
     *            \p max_iter iterations are reached.
     *
     * \param p_init          (in)     initial 6-DOF parameters (e.g. identity)
     * \param p_best          (out)    optimised 6-DOF parameters
     * \param surface         (in)     WM surface in surface-native RAS space
     * \param vol             (in)     floating-point volume data (moving image)
     * \param nii_ptr         (in)     NIfTI header of the moving volume
     * \param dims            (in)     volume dimensions [nx, ny, nz]
     * \param wm_dist         (in)     WM sampling offset (mm)
     * \param gm_dist         (in)     GM sampling offset (mm)
     * \param slope           (in)     BBR cost saturation slope
     * \param invert_contrast (in)     non-zero for T2/BOLD
     * \param grid_range_mm   (in)     half-width of brute-force translation grid (mm)
     * \param grid_range_rad  (in)     half-width of brute-force rotation grid (rad)
     * \param grid_steps      (in)     number of steps per dimension in grid search
     *                                 (total evaluations = (2*grid_steps+1)^6)
     * \param max_iter        (in)     maximum Powell iterations in Stage 2
     * \param tol             (in)     fractional convergence tolerance for Powell
     * \param verbose         (in)     non-zero to print progress
     * \return final BBR cost after optimisation
     */
    double CAT_BBReg_optimise(const CAT_RigidParams *p_init,
                              CAT_RigidParams *p_best,
                              const polygons_struct *surface,
                              float *vol,
                              nifti_image *nii_ptr,
                              int dims[3],
                              double wm_dist,
                              double gm_dist,
                              double slope,
                              int invert_contrast,
                              double grid_range_mm,
                              double grid_range_rad,
                              int grid_steps,
                              int max_iter,
                              double tol,
                              int verbose);

    /**
     * \brief Write a 4×4 RAS-to-RAS transform to a plain-text file.
     *
     * Format (compatible with FSL flirt -init / ANTs):
     *
     *   # CAT_BBReg RAS-to-RAS transform
     *   m00 m01 m02 m03
     *   m10 m11 m12 m13
     *   m20 m21 m22 m23
     *   0   0   0   1
     *
     * \param filename (in) output file path
     * \param m        (in) 4×4 row-major matrix
     * \return 0 on success, -1 on I/O error
     */
    int CAT_BBReg_write_matrix(const char *filename, const double m[16]);

#ifdef __cplusplus
}
#endif

#endif /* CAT_BBREG_H */
