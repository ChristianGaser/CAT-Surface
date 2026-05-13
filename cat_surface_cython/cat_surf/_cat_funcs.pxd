# cython: language_level=3
"""
Cython declarations for libCAT public functions.
"""
from cat_surf._bic_types cimport (
    polygons_struct, object_struct, Status, File_formats, Point, Vector,
    Colour, Object_types, Smallest_int
)
from cat_surf._nifti_types cimport nifti_image


# ---------------------------------------------------------------------------
# CAT_SurfaceIO.h — Surface and value I/O
# ---------------------------------------------------------------------------
cdef extern from "CAT_SurfaceIO.h":
    Status input_graphics_any_format(char *filename,
                                     File_formats *fmt,
                                     int *n_objects,
                                     object_struct ***objects)

    Status output_graphics_any_format(char *filename,
                                      File_formats fmt,
                                      int n_objects,
                                      object_struct **objects,
                                      double *opacity)

    Status input_values_any_format(char *filename,
                                   int *n_values,
                                   double **values)

    Status output_values_any_format(const char *filename,
                                    int n_values,
                                    void *values,
                                    int datatype)

    int input_freesurfer(char *, File_formats *, int *, object_struct ***)
    int output_freesurfer(char *, File_formats, int, object_struct **)
    int input_freesurfer_curv(char *, int *, double **)
    int output_freesurfer_curv(char *, int, double *)

    # ATABLE is the FreeSurfer .annot colour-table entry struct.  We
    # only ever pass it through opaquely between read_/write_annotation_table.
    ctypedef struct ATABLE:
        pass
    int read_annotation_table(char *filename, int *n_values,
                              int **annot_out, int *n_table,
                              ATABLE **atable_out)
    int write_annotation_table(char *filename, int n_values, int *annot,
                               int n_table, ATABLE *atable)

    int input_gifti(char *, File_formats *, int *, object_struct ***,
                    int *, double **)
    int output_gifti(char *, File_formats, int, object_struct **, double *)
    int input_gifti_curv(char *, int *, double **)
    int output_gifti_curv(char *, int, double *)


# ---------------------------------------------------------------------------
# CAT_Surf.h — Core surface operations
# ---------------------------------------------------------------------------
cdef extern from "CAT_Surf.h":
    double get_area_of_points(polygons_struct *poly, double *area_values)
    double get_vertex_areas(polygons_struct *poly, double *areas)
    void   get_radius_of_points(polygons_struct *poly, double *radii)
    double get_sphere_radius(polygons_struct *poly)

    int    euler_characteristic(polygons_struct *poly, int verbose)

    double compute_point_distance(polygons_struct *s1, polygons_struct *s2,
                                  double *distances, int symmetric)
    double compute_point_distance_mean(polygons_struct *s1,
                                       polygons_struct *s2,
                                       double *distances, int symmetric)
    double compute_point_hausdorff(polygons_struct *s1, polygons_struct *s2,
                                   double *distances, int symmetric)
    double compute_exact_hausdorff(polygons_struct *s1, polygons_struct *s2,
                                   double *distances, int symmetric)

    void correct_bounds_to_target(polygons_struct *source,
                                  polygons_struct *target)
    void correct_bounds_to_target_with_scaling(polygons_struct *source,
                                               polygons_struct *target)

    void apply_warp(polygons_struct *source, polygons_struct *target,
                    double *warp_array, int *warp_dims, int n_steps)

    void areal_smoothing(polygons_struct *poly, double fwhm,
                         int n_iter, int itertype, int *mask, int verbose)
    void linear_smoothing(polygons_struct *poly, double fwhm,
                          int n_iter, int itertype, int *mask, int verbose)
    void distance_smoothing(polygons_struct *poly, double fwhm,
                            int n_iter, int itertype, int *mask, int verbose)

    int  reduce_mesh_quadrics(polygons_struct *poly, int target_faces,
                              double aggressiveness, int preserve_sharp,
                              int verbose)

    void central_to_pial(polygons_struct *central, double *thickness,
                         double *extents, int check_intersect,
                         double sigma, int iterations, int verbose)

    object_struct ** central_to_new_pial(polygons_struct *polygons,
                                         double *thickness,
                                         double *extents,
                                         int check_intersect,
                                         double sigma, int iterations,
                                         int verbose)

    void inflate_surface_and_smooth_fingers(
        polygons_struct *poly, const int iterations,
        const double max_thickness, const int laplace_iter,
        const double laplace_lambda, const double laplace_mu,
        const double laplace_dt, const int verbose)


# ---------------------------------------------------------------------------
# CAT_Smooth.h — Smoothing
# ---------------------------------------------------------------------------
cdef extern from "CAT_Smooth.h":
    void get_all_polygon_point_neighbours(polygons_struct *poly,
                                          int *n_neighbours_ptr[],
                                          int **neighbours_ptr[])
    void smooth_heatkernel(polygons_struct *poly, double *values,
                           double kernel_size)
    int  smooth_laplacian(polygons_struct *poly, int iters,
                          double alpha, double beta)


# ---------------------------------------------------------------------------
# CAT_Curvature.h — Curvature measures
# ---------------------------------------------------------------------------
cdef extern from "CAT_Curvature.h":
    void get_smoothed_curvatures(polygons_struct *poly, double *curvatures,
                                 double fwhm, int n_iter)
    void compute_sulcus_depth(polygons_struct *poly, double *depth)
    void compute_convexity(polygons_struct *poly, int n_neighbours[],
                           int *neighbours[], double *convexity)
    void compute_local_sharpness(polygons_struct *poly, int n_neighbours[],
                                 int *neighbours[], double *sharpness)
    # Mirror of CAT_SurfCurvature's main routine.  curvtype values are
    # the same as the CLI's; distance is 3.0 for curvtype 0, else 0.0.
    void get_polygon_vertex_curvatures_cg(
        polygons_struct *poly, int n_neighbours[], int *neighbours[],
        double distance, int curvtype, double *curvatures)


# ---------------------------------------------------------------------------
# CAT_SurfUtils.h — Sphere mapping, area normalization
# ---------------------------------------------------------------------------
cdef extern from "CAT_SurfUtils.h":
    void surf_to_sphere(polygons_struct *polygons, int n_triangles,
                        int verbose)
    double get_area_of_points_normalized_to_sphere(
        polygons_struct *polygons, polygons_struct *sphere,
        double *areas)
    void translate_to_center_of_mass(polygons_struct *polygons)
    void set_vector_length(Point *p, double newLength)


# ---------------------------------------------------------------------------
# CAT_CorrectThicknessFolding.h
# ---------------------------------------------------------------------------
cdef extern from "CAT_CorrectThicknessFolding.h":
    Status CAT_CorrectThicknessFoldingWeighted(
        polygons_struct *polygons, int n_vals,
        double *thickness, double slope)
    Status CAT_CorrectThicknessFolding(
        polygons_struct *polygons, int n_vals,
        double *thickness)


# ---------------------------------------------------------------------------
# CAT_Deform.h — Surface deformation
# ---------------------------------------------------------------------------
cdef extern from "CAT_Deform.h":
    void surf_deform(polygons_struct *polygons, float *input,
                     nifti_image *nii_ptr, double w[3], double sigma,
                     float lim, int it, int remove_intersections,
                     int verbose)


# ---------------------------------------------------------------------------
# CAT_MeshClean.h — Self-intersection removal (MeshFix)
# ---------------------------------------------------------------------------
cdef extern from "CAT_MeshClean.h":
    ctypedef struct CAT_MeshCleanOptions:
        int max_iters
        int inner_loops
        int fill_holes
        int verbose

    void CAT_MeshCleanOptionsInit(CAT_MeshCleanOptions *opts)
    int  CAT_SurfMeshClean(polygons_struct *polygons,
                           const CAT_MeshCleanOptions *opts)
    int  CAT_SurfCountIntersections(polygons_struct *polygons)


# ---------------------------------------------------------------------------
# CAT_Resample.h — Spherical resampling
# ---------------------------------------------------------------------------
cdef extern from "CAT_Resample.h":
    object_struct ** resample_surface_to_target_sphere(
        polygons_struct *polygons,
        polygons_struct *polygons_sphere,
        polygons_struct *target_sphere,
        double *input_values,
        double *output_values,
        int label_interpolation,
        int areal_interpolation)


# ---------------------------------------------------------------------------
# CAT_SurfPialWhite.h — Pial / white estimation
# ---------------------------------------------------------------------------
cdef extern from "CAT_SurfPialWhite.h":
    ctypedef struct CAT_PialWhiteOptions:
        double w1
        double w2
        double w3
        double sigma
        int iterations
        int gradient_iterations
        int method
        int verbose

    void CAT_PialWhiteOptionsInit(CAT_PialWhiteOptions *opts)
    int  CAT_SurfEstimatePialWhite(
             polygons_struct *central,
             const double *thickness_values,
             float *labels,
             nifti_image *nii_ptr,
             polygons_struct *pial_out,
             polygons_struct *white_out,
             const CAT_PialWhiteOptions *opts)


# ---------------------------------------------------------------------------
# CAT_SurfWarpDartel.h — declared in _dartel.pxd
# ---------------------------------------------------------------------------
# dartel.h has no include guards AND defines an unhygienic `#define S 1.0`
# macro that breaks any file including stdint-based unions. To avoid
# polluting every cimport, the DARTEL declarations live in cat_surf/_dartel.pxd
# and are only cimported by cat_surf/_surf_warp.pyx.


# ---------------------------------------------------------------------------
# CAT_Warp.h — Warp/rotation primitives used by SurfWarp
# ---------------------------------------------------------------------------
cdef extern from "CAT_Warp.h":
    int INVERSE_WARPING
    void rotate_polygons(polygons_struct *src, polygons_struct *dst,
                         double *rotation_matrix)
    void rotation_to_matrix(double *matrix, double alpha, double beta,
                            double gamma)
    void average_xz_surf(polygons_struct *a, polygons_struct *b,
                         polygons_struct *out)


# ---------------------------------------------------------------------------
# CAT_Vol2SurfUtils.h — helpers for CAT_Vol2Surf
# ---------------------------------------------------------------------------
cdef extern from "CAT_Vol2SurfUtils.h":
    double CAT_Vol2SurfEvaluateFunction(const double *val_array, int n_val,
                                        int map_func, const double *kernel,
                                        int *index_out)
    void   CAT_Vol2SurfBuildExpKernel(const double *length_array, int n,
                                      double exp_half, double *kernel_out)
    void   CAT_Vol2SurfBuildGaussianKernel50(int grid_steps, int grid_steps1,
                                             double *kernel_out)


# ---------------------------------------------------------------------------
# CAT_Math.h — F_* mapping function enum + clip_data
# ---------------------------------------------------------------------------
cdef extern from "CAT_Math.h":
    int F_MEAN
    int F_MIN
    int F_MAX
    int F_STD
    int F_SUM
    int F_MAXABS
    int F_EXP
    int F_MEDIAN
    int F_RANGE
    int F_WAVERAGE
    int F_MULTI


# ---------------------------------------------------------------------------
# Extras: isoval (used by Vol2Surf), surface normals & areal helpers
# ---------------------------------------------------------------------------
cdef extern from "CAT_Vol.h":
    float isoval(float vol[], float x, float y, float z, int s[],
                 nifti_image *nii_ptr)

cdef extern from "CAT_Surf.h":
    double get_area_of_points_central_to_pial(
        polygons_struct *polygons, double *area_pial_points,
        double *thickness_values, double extent)


# ---------------------------------------------------------------------------
# CAT_MarchingCubes.h — Isosurface extraction
# ---------------------------------------------------------------------------
cdef extern from "CAT_MarchingCubes.h":
    object_struct *apply_marching_cubes(
        float *input_float,
        nifti_image *nii_ptr,
        float *label,
        double min_threshold,
        double pre_fwhm,
        int iter_laplacian,
        double dist_morph,
        int n_median_filter,
        int n_iter,
        double strength_gyri_mask,
        int verbose)

    object_struct *apply_marching_cubes_fast(
        float *input_float,
        nifti_image *nii_ptr,
        double min_threshold,
        int iter_laplacian,
        int verbose)


# ---------------------------------------------------------------------------
# CAT_Nlm.h — Non-local means denoising
# ---------------------------------------------------------------------------
cdef extern from "CAT_Nlm.h":
    void sanlm(float *ima, int v, int f, int is_rician,
               double strength, const int *dims)
    void ornlm(float *ima, int v, int f, float h, float sigma,
               const int *dims)


# ---------------------------------------------------------------------------
# CAT_VolPbt.h — Projection-based thickness
# ---------------------------------------------------------------------------
cdef extern from "CAT_VolPbt.h":
    ctypedef struct CAT_PbtOptions:
        int n_avgs
        int n_median_filter
        int median_subsample
        double range
        double fill_thresh
        double correct_voxelsize
        double sulcal_width
        int fast
        int verbose

    void CAT_PbtOptionsInit(CAT_PbtOptions *opts)
    int  CAT_VolComputePbt(
             const float *src,
             float *GMT_out,
             float *PPM_out,
             float *dist_CSF_out,
             float *dist_WM_out,
             int dims[3],
             double voxelsize[3],
             const CAT_PbtOptions *opts)

    void blood_vessel_correction_pve_float(float *Yp0, int dims[3],
                                           double vx_vol[3])


# ---------------------------------------------------------------------------
# CAT_Amap.h / CAT_Bmap.h — Brain tissue segmentation
# ---------------------------------------------------------------------------
cdef extern from "CAT_Amap.h":
    void Amap(float *src, unsigned char *label, unsigned char *prob,
              double *mean, int nc, int niters, int sub, int *dims,
              int pve, double weight_MRF, double *voxelsize,
              int niters_ICM, int verbose, int use_median,
              const double *mrf_class_weights, int use_multistep)
    void Pve5(float *src, unsigned char *prob, unsigned char *label,
              double *mean, int *dims)

cdef extern from "CAT_Bmap.h":
    void Bmap(float *src, unsigned char *label, unsigned char *prob,
              double *mean, int n_classes, int BG, int niters,
              int a, int b, int c, float *bias, int *dims,
              int pve, int verbose)


# ---------------------------------------------------------------------------
# CAT_NiftiLib.h — NIfTI float loader
# ---------------------------------------------------------------------------
cdef extern from "CAT_NiftiLib.h":
    nifti_image *read_nifti_float(const char *filename,
                                  float **data, int read_data)


# ---------------------------------------------------------------------------
# CAT_Vol.h — Volume smoothing (partial; see full header for other funcs)
# ---------------------------------------------------------------------------
cdef extern from "CAT_Vol.h":
    void smooth3(void *vol, int dims[3], double voxelsize[3],
                 double s[3], int use_mask, int datatype)


# ---------------------------------------------------------------------------
# CAT_BBReg.h — Boundary-Based Registration
# ---------------------------------------------------------------------------
cdef extern from "CAT_BBReg.h":
    ctypedef struct CAT_RigidParams:
        double tx, ty, tz
        double rx, ry, rz

    ctypedef struct CAT_SurfData:
        polygons_struct *surface
        float           *cortex_mask
        float           *thickness
        double           gm_proj_frac

    void   CAT_BBReg_params_to_matrix(const CAT_RigidParams *p,
                                      double m[16])
    void   CAT_BBReg_apply_matrix(polygons_struct *surface,
                                  const double m[16])
    void   CAT_BBReg_invert_matrix(const double m[16], double inv[16])

    double CAT_BBReg_cost(const CAT_RigidParams *p,
                          const CAT_SurfData *surfs, int n_surfs,
                          float *vol, nifti_image *nii_ptr, int dims[3],
                          double wm_dist, double gm_dist,
                          double slope, int invert_contrast)

    double CAT_BBReg_optimise(const CAT_RigidParams *p_init,
                              CAT_RigidParams *p_best,
                              const CAT_SurfData *surfs, int n_surfs,
                              float *vol, nifti_image *nii_ptr, int dims[3],
                              double wm_dist, double gm_dist,
                              double slope, int invert_contrast,
                              double grid_range_mm, double grid_range_rad,
                              int grid_steps, int max_iter, double tol,
                              int verbose)

    int    CAT_BBReg_detect_contrast(const CAT_RigidParams *p,
                                     const CAT_SurfData *surfs, int n_surfs,
                                     float *vol, nifti_image *nii_ptr,
                                     int dims[3],
                                     double wm_dist, double gm_dist,
                                     int verbose)

    int    CAT_BBReg_write_matrix(const char *filename, const double m[16])


# ---------------------------------------------------------------------------
# CAT_VolumeReg.h — Volume-to-volume rigid registration
# ---------------------------------------------------------------------------
cdef extern from "CAT_VolumeReg.h":
    double CAT_VolumeReg_register(float *fixed_vol,
                                  nifti_image *fixed_nii,
                                  int fixed_dims[3],
                                  float *moving_vol,
                                  nifti_image *moving_nii,
                                  int moving_dims[3],
                                  CAT_RigidParams *p_out,
                                  int n_levels, double sat_k,
                                  int max_iter, int verbose)

    double CAT_VolumeReg_register_NMI(float *fixed_vol,
                                      nifti_image *fixed_nii,
                                      int fixed_dims[3],
                                      float *moving_vol,
                                      nifti_image *moving_nii,
                                      int moving_dims[3],
                                      CAT_RigidParams *p_out,
                                      int n_levels, int n_bins,
                                      int max_iter, int verbose)
