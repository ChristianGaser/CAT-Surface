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

    /* -----------------------------------------------------------------------
     * One surface with optional per-vertex data.
     *
     * cortex_mask  : array of n_points floats; vertex included if > 0.5.
     *                Pass NULL to include all vertices.
     * thickness    : array of n_points floats with cortical thickness (mm).
     *                Used for fractional GM projection when gm_proj_frac > 0.
     *                Pass NULL to use absolute gm_dist instead.
     * gm_proj_frac : fraction of local thickness for GM sampling offset
     *                (0 = use gm_dist_abs, 1 = full thickness).
     *                Ignored when thickness is NULL.
     * ----------------------------------------------------------------------- */
    typedef struct
    {
        polygons_struct *surface;  /* white-matter surface mesh */
        float           *cortex_mask;  /* per-vertex cortex label (or NULL) */
        float           *thickness;    /* per-vertex thickness in mm (or NULL) */
        double           gm_proj_frac; /* fraction of thickness for GM offset */
    } CAT_SurfData;

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
     * \brief Compute the BBR cost over one or more surfaces.
     *
     * For each vertex (not masked out by cortex_mask) the volume is sampled at
     * two points along the outward surface normal:
     *   - WM side : vertex − wm_dist * normal
     *   - GM side : vertex + d_gm   * normal
     *
     * where d_gm = (sd->gm_proj_frac * thickness[i]) when thickness data is
     * available, otherwise d_gm = gm_dist.
     *
     * The normalised contrast at vertex i is:
     *   c_i = (I_WM − I_GM) / (0.5*(I_WM + I_GM) + eps)
     *
     * Cost = mean_i { 1 − tanh^2(slope * c_i) }  (lower = better alignment).
     * For T2/BOLD pass invert_contrast != 0 to negate c_i.
     *
     * \param p               (in)  current 6-DOF rigid parameters
     * \param surfs           (in)  array of CAT_SurfData structures
     * \param n_surfs         (in)  number of elements in surfs[]
     * \param vol             (in)  floating-point volume data (moving image)
     * \param nii_ptr         (in)  NIfTI header of the moving volume
     * \param dims            (in)  volume dimensions [nx, ny, nz]
     * \param wm_dist         (in)  WM sampling offset (mm, positive)
     * \param gm_dist         (in)  absolute GM sampling offset (mm); overridden
     *                              per-vertex when thickness + gm_proj_frac are set
     * \param slope           (in)  saturation slope of the BBR cost function
     * \param invert_contrast (in)  non-zero for T2/BOLD
     * \return mean BBR cost (lower = better alignment)
     */
    double CAT_BBReg_cost(const CAT_RigidParams *p,
                          const CAT_SurfData *surfs,
                          int n_surfs,
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
     * Stage 1 — brute-force translation grid search over ±grid_range_mm.
     * Stage 2 — Powell's direction-set method for full 6-DOF refinement.
     *
     * \param p_init          (in)     initial 6-DOF parameters
     * \param p_best          (out)    optimised 6-DOF parameters
     * \param surfs           (in)     array of CAT_SurfData (lh + rh, or just one)
     * \param n_surfs         (in)     number of elements in surfs[]
     * \param vol             (in)     floating-point volume data (moving image)
     * \param nii_ptr         (in)     NIfTI header of the moving volume
     * \param dims            (in)     volume dimensions [nx, ny, nz]
     * \param wm_dist         (in)     WM sampling offset (mm)
     * \param gm_dist         (in)     GM absolute sampling offset (mm)
     * \param slope           (in)     BBR cost saturation slope
     * \param invert_contrast (in)     non-zero for T2/BOLD
     * \param grid_range_mm   (in)     half-width of brute-force translation grid (mm)
     * \param grid_range_rad  (in)     half-width of brute-force rotation grid (rad)
     * \param grid_steps      (in)     steps per DOF in grid search
     * \param max_iter        (in)     maximum Powell iterations
     * \param tol             (in)     convergence tolerance for Powell
     * \param verbose         (in)     non-zero to print progress
     * \return final BBR cost after optimisation
     */
    double CAT_BBReg_optimise(const CAT_RigidParams *p_init,
                              CAT_RigidParams *p_best,
                              const CAT_SurfData *surfs,
                              int n_surfs,
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
