# cython: language_level=3
"""
Cython declarations for libCAT public functions.
"""
from cat_surf._bic_types cimport (
    polygons_struct, object_struct, Status, File_formats, Point, Vector,
    Colour, Object_types, Smallest_int
)


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
                         double *label, int n_labels, double intensity,
                         int check_intersect, int verbose)

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
