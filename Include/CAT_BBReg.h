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

    /**
     * \brief Rigid-body transform, 6 DOF.
     *
     * The convention follows FreeSurfer bbregister:
     * T_full = T_trans * R_z * R_y * R_x
     */
    typedef struct
    {
        double tx; /**< translation along x (mm) */
        double ty; /**< translation along y (mm) */
        double tz; /**< translation along z (mm) */
        double rx; /**< rotation about x (radians) */
        double ry; /**< rotation about y (radians) */
        double rz; /**< rotation about z (radians) */
    } CAT_RigidParams;

    /**
     * \brief One surface with optional per-vertex data.
     *
     * cortex_mask  : array of n_points floats; vertex included if > 0.5.
     *                Pass NULL to include all vertices.
     * thickness    : array of n_points floats with cortical thickness (mm).
     *                Used for fractional GM projection when gm_proj_frac > 0.
     *                Pass NULL to use absolute gm_dist instead.
     * gm_proj_frac : fraction of local thickness for GM sampling offset
     *                (0 = use gm_dist_abs, 1 = full thickness).
     *                Ignored when thickness is NULL.
     */
    typedef struct
    {
        polygons_struct *surface; /**< white-matter surface mesh */
        float *cortex_mask;       /**< per-vertex cortex label (or NULL) */
        float *thickness;         /**< per-vertex thickness in mm (or NULL) */
        double gm_proj_frac;      /**< fraction of thickness for GM offset */
    } CAT_SurfData;

    /**
     * \brief Affine transform, 12 DOF: a rigid core with scales and shears.
     *
     *   M = T_trans * Rz * Ry * Rx * T(centre) * K * T(-centre)
     *
     *       | sx  hxy hxz |
     *   K = | 0   sy  hyz |
     *       | 0   0   sz  |
     *
     * M maps the fixed surface (T1) RAS space to the moving volume (EPI) RAS
     * space, as the matrix of CAT_RigidParams does.  The scales and shears act
     * along the axes of the fixed space about centre, which they leave where
     * the rigid core puts it; with unit scales and no shear M is the rigid
     * matrix of \p rigid exactly.  9 DOF are the rigid core and the scales,
     * 12 DOF add the shears.
     */
    typedef struct
    {
        CAT_RigidParams rigid; /**< rotation (about the origin) and translation */
        double sx, sy, sz;     /**< scales along the fixed x, y, z axes (1 = none) */
        double hxy, hxz, hyz;  /**< shears (0 = none) */
        double centre[3];      /**< fixed RAS point the scales and shears keep */
    } CAT_AffineParams;

    /**
     * \brief Settings of a boundary-based registration as CAT_SurfBBReg runs it.
     *
     * The sampling, cost, search, smoothing and volume-initialisation settings
     * of the whole pipeline, shared by CAT_SurfBBReg and cat_surf.bbreg.  Set
     * the defaults with CAT_BBRegOptionsInit(), the single source of truth for
     * both, and change what you need.
     */
    typedef struct
    {
        double wm_dist;         /**< WM sampling distance inside the surface (mm) */
        double gm_dist;         /**< GM sampling distance outside the surface (mm),
                                     used without thickness or gm_proj_frac */
        double gm_proj_frac;    /**< GM distance as a fraction of the thickness
                                     (0 = gm_dist) */
        double slope;           /**< saturation slope of the cost */
        int    invert_contrast; /**< -1 auto-detect, 0 T1/FLAIR, 1 T2/BOLD */
        double grid_range_mm;   /**< half-width of the translation grid search (mm) */
        double grid_range_rad;  /**< half-width of a rotation grid (rad); the
                                     grid searches translations only */
        int    grid_steps;      /**< grid steps on either side of the start */
        int    max_iter;        /**< Powell iterations */
        double tol;             /**< Powell convergence tolerance */
        int    dof;             /**< 6, 9 or 12, see CAT_BBReg_optimise_affine() */
        double fwhm;            /**< smoothing of both volumes for the volume
                                     initialisation (mm, 0 = off) */
        double fwhm_bbr;        /**< smoothing of the moving volume for BBR
                                     (mm, 0 = off) */
        int    vol_reg_levels;  /**< pyramid levels of the volume initialisation */
        int    vol_reg_bins;    /**< joint-histogram bins per axis */
        int    vol_reg_iter;    /**< Powell iterations per pyramid level */
    } CAT_BBRegOpts;

    /**
     * \brief Set the default settings of a boundary-based registration.
     *
     * WM and GM are sampled 0.5 mm inside and outside the surface, with a
     * cost slope of 0.5 and the contrast detected from the image; a grid of
     * translations of up to 4 mm in 1 mm steps precedes Powell's method (200
     * iterations, tolerance 1e-6), rigidly.  The volume initialisation (NMI, 4
     * levels, 64 bins, 30 iterations) sees both volumes smoothed by 5 mm, BBR
     * the moving one by 2 mm.
     *
     * \param opts (out) settings to initialise
     */
    void CAT_BBRegOptionsInit(CAT_BBRegOpts *opts);

    /**
     * \brief Convert 6-DOF rigid parameters to a 4×4 homogeneous matrix (row-major).
     *
     * The resulting matrix maps points expressed in the *fixed* surface's (T1) RAS
     * space to the *moving* volume's (EPI) RAS space.  This is the direction used
     * internally by the BBR cost function: surface sample points are transformed
     * by M before trilinear sampling of the EPI.
     *
     *   M = T_trans * Rz * Ry * Rx
     *
     * To obtain the standard output convention (EPI → T1, matching bbregister /
     * bbreg), invert this matrix with CAT_BBReg_invert_matrix().
     *
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
     * \brief Invert a 4×4 affine homogeneous matrix (last row 0 0 0 1).
     *
     * Unlike CAT_BBReg_invert_matrix() also for matrices with scales and
     * shears, such as those of CAT_BBReg_affine_to_matrix().
     *
     * \param m    (in)  4×4 row-major affine matrix
     * \param inv  (out) 4×4 row-major inverse
     */
    void CAT_BBReg_invert_affine(const double m[16], double inv[16]);

    /**
     * \brief Convert affine parameters to a 4×4 homogeneous matrix (row-major).
     *
     * M = T_trans * Rz * Ry * Rx * T(centre) * K * T(-centre), see
     * CAT_AffineParams; the direction is that of CAT_BBReg_params_to_matrix().
     *
     * \param p  (in)  affine parameters
     * \param m  (out) 4×4 matrix, row-major, 16 doubles
     */
    void CAT_BBReg_affine_to_matrix(const CAT_AffineParams *p, double m[16]);

    /**
     * \brief Affine parameters of a rigid transform, centred on the surfaces.
     *
     * Unit scales, no shear, and as centre the mean position of the vertices
     * the cost uses (cortex_mask > 0.5), so that scaling does not move the
     * brain.  The start of CAT_BBReg_optimise_affine() from a rigid result.
     *
     * \param rigid   (in)  rigid core, e.g. from CAT_BBReg_optimise()
     * \param surfs   (in)  array of CAT_SurfData structures
     * \param n_surfs (in)  number of elements in surfs[]
     * \param p       (out) affine parameters
     */
    void CAT_BBReg_affine_init(const CAT_RigidParams *rigid,
                               const CAT_SurfData *surfs, int n_surfs,
                               CAT_AffineParams *p);

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
     * The percent contrast at vertex i is (×100 following neuroreg/bbregister):
     *   c_i = 100 * (I_WM − I_GM) / (0.5*(I_WM + I_GM) + eps)
     *
     * Cost = mean_i { 1 − tanh^2(slope * c_i) }  (lower = better alignment).
     * With slope=0.5 (default) the cost saturates (≈0) once |c_i| ≳ 6 (6% contrast).
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
     * \brief Compute the BBR cost at any 4×4 transform.
     *
     * CAT_BBReg_cost() for a matrix instead of rigid parameters, e.g. one of
     * CAT_BBReg_affine_to_matrix().  The sampling positions are taken along
     * the normals of the fixed surface and then transformed, so an affine
     * transform also scales the sampling distances.
     *
     * \param m               (in)  4×4 row-major transform, fixed (T1) RAS →
     *                              moving (EPI) RAS
     * \param surfs           (in)  array of CAT_SurfData structures
     * \param n_surfs         (in)  number of elements in surfs[]
     * \param vol             (in)  floating-point volume data (moving image)
     * \param nii_ptr         (in)  NIfTI header of the moving volume
     * \param dims            (in)  volume dimensions [nx, ny, nz]
     * \param wm_dist         (in)  WM sampling offset (mm, positive)
     * \param gm_dist         (in)  absolute GM offset (mm); used when no thickness
     * \param slope           (in)  saturation slope of the cost function
     * \param invert_contrast (in)  non-zero for T2/BOLD
     * \return mean BBR cost (lower = better alignment)
     */
    double CAT_BBReg_cost_matrix(const double m[16],
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
     * \brief Refine a BBR estimate with scales (9 DOF) or scales and shears (12 DOF).
     *
     * Powell's direction-set method over the first dof parameters
     * (translations, rotations, scales, shears; the others keep their values
     * from p_init).  Without a grid search, so p_init has to be close:
     * normally CAT_BBReg_affine_init() of the result of CAT_BBReg_optimise().
     * Scales and shears start with Powell steps of 0.005.
     *
     * \param p_init          (in)  initial parameters (and the centre)
     * \param p_best          (out) optimised parameters
     * \param dof             (in)  6, 9 or 12
     * \param surfs           (in)  array of CAT_SurfData (lh + rh, or just one)
     * \param n_surfs         (in)  number of elements in surfs[]
     * \param vol             (in)  floating-point volume data (moving image)
     * \param nii_ptr         (in)  NIfTI header of the moving volume
     * \param dims            (in)  volume dimensions [nx, ny, nz]
     * \param wm_dist         (in)  WM sampling offset (mm)
     * \param gm_dist         (in)  GM absolute sampling offset (mm)
     * \param slope           (in)  BBR cost saturation slope
     * \param invert_contrast (in)  non-zero for T2/BOLD
     * \param max_iter        (in)  maximum Powell iterations
     * \param tol             (in)  convergence tolerance for Powell
     * \param verbose         (in)  non-zero to print progress
     * \return final BBR cost, or -1 for an unsupported dof (p_best = p_init)
     */
    double CAT_BBReg_optimise_affine(const CAT_AffineParams *p_init,
                                     CAT_AffineParams *p_best,
                                     int dof,
                                     const CAT_SurfData *surfs,
                                     int n_surfs,
                                     float *vol,
                                     nifti_image *nii_ptr,
                                     int dims[3],
                                     double wm_dist,
                                     double gm_dist,
                                     double slope,
                                     int invert_contrast,
                                     int max_iter,
                                     double tol,
                                     int verbose);

    /**
     * \brief Automatically detect image contrast type from WM/GM intensity.
     *
     * Samples volume intensities at WM-side and GM-side positions along the
     * surface normals using the supplied transform \p p.  The sign of the
     * resulting WM/GM contrast determines the image type:
     *
     *   - T1 / FLAIR : WM brighter than GM  → returns 0
     *   - T2 / BOLD  : WM darker  than GM   → returns 1
     *
     * The result can be passed directly as the \p invert_contrast argument to
     * CAT_BBReg_cost() and CAT_BBReg_optimise().
     *
     * A call at the initial (pre-registration) transform is sufficient for
     * data that are roughly aligned.  At least 100 valid vertices must fall
     * inside the volume FOV and the absolute WM/GM contrast must exceed 1 %
     * for a reliable detection; otherwise -1 is returned and the caller
     * should fall back to a default (usually 0 = T1).
     *
     * \param p        (in)  rigid parameters used to transform surface points
     * \param surfs    (in)  array of CAT_SurfData structures
     * \param n_surfs  (in)  number of elements in surfs[]
     * \param vol      (in)  floating-point volume data (moving image)
     * \param nii_ptr  (in)  NIfTI header of the moving volume
     * \param dims     (in)  volume dimensions [nx, ny, nz]
     * \param wm_dist  (in)  WM sampling offset (mm, positive, inward)
     * \param gm_dist  (in)  absolute GM sampling offset (mm)
     * \param verbose  (in)  non-zero to print diagnostic line
     * \return 0 (T1/FLAIR), 1 (T2/BOLD), -1 (undetermined / too few samples)
     */
    int CAT_BBReg_detect_contrast(const CAT_RigidParams *p,
                                  const CAT_SurfData *surfs,
                                  int n_surfs,
                                  float *vol,
                                  nifti_image *nii_ptr,
                                  int dims[3],
                                  double wm_dist,
                                  double gm_dist,
                                  int verbose);

    /**
     * \brief Write a 4×4 RAS-to-RAS transform to a plain-text file.
     *
     * Format:
     *
     *   # CAT_BBReg RAS-to-RAS transform
     *   m00 m01 m02 m03
     *   m10 m11 m12 m13
     *   m20 m21 m22 m23
     *   0   0   0   1
     *
     * This is a world-space (scanner RAS) matrix, as FreeSurfer's RAS2RAS
     * LTAs and nibabel affines are; it is NOT what FSL or ANTs read.  flirt
     * -init/-applyxfm expect FSL's scaled-voxel convention (convert with
     * ras_to_fsl_matrix() and write with write_fsl_matrix(), or use
     * CAT_SurfBBReg -fsl-matrix), and ANTs reads ITK transforms in LPS.
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
