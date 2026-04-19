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
                              double quality, int iterations, int verbose)

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


# ---------------------------------------------------------------------------
# CAT_SurfUtils.h — Sphere mapping, area normalization
# ---------------------------------------------------------------------------
cdef extern from "CAT_SurfUtils.h":
    void surf_to_sphere(polygons_struct *polygons, int n_triangles,
                        int verbose)
    double get_area_of_points_normalized_to_sphere(
        polygons_struct *polygons, polygons_struct *sphere,
        double *areas)


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
# CAT_SurfWarpDartel.h — DARTEL-based spherical registration (DEFERRED)
# ---------------------------------------------------------------------------
# dartel.h has no include guards; declaring dartel_prm here causes
# redefinition errors when CAT_SurfWarpDartel.h is processed.
# SurfWarp wrapping is deferred due to DARTEL complexity.


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
