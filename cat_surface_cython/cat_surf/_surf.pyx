# cython: language_level=3, boundscheck=False, wraparound=False
"""
Python wrappers for the core libCAT surface operations.

Every function accepts numpy arrays (vertices, faces, values) and
returns numpy arrays.  The conversion to/from the C polygons_struct
is handled transparently.
"""

import numpy as np
cimport numpy as cnp
from libc.stdlib cimport malloc, free
from libc.string cimport memcpy, memset

from cat_surf._bic_types cimport (
    polygons_struct, object_struct, Object_types, Status,
    create_object, get_polygons_ptr, delete_object, get_object_type,
    initialize_polygons, compute_polygon_normals, POLYGONS, OK,
    File_formats,
)
from cat_surf._nifti_types cimport nifti_image, nifti_image_read, nifti_image_free
from cat_surf._convert cimport PolygonsMesh, _wrap_object
from cat_surf._convert import arrays_to_polygons, polygons_to_arrays
from cat_surf._volume cimport VolumeHandle, open_volume

cimport cat_surf._cat_funcs as C

cnp.import_array()


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
cdef PolygonsMesh _ensure_mesh(vertices, faces):
    """Coerce (vertices, faces) into a PolygonsMesh."""
    verts = np.ascontiguousarray(vertices, dtype=np.float64)
    if verts.ndim != 2 or verts.shape[1] != 3:
        raise ValueError("vertices must have shape (V, 3)")
    tris = np.ascontiguousarray(faces, dtype=np.int32)
    if tris.ndim != 2 or tris.shape[1] != 3:
        raise ValueError("faces must have shape (F, 3)")
    return arrays_to_polygons(verts, tris)


# ===================================================================
# Surface area
# ===================================================================
def get_area(vertices, faces):
    """
    Compute per-vertex area of a triangular mesh.

    Parameters
    ----------
    vertices : array_like, shape (V, 3)
    faces    : array_like, shape (F, 3)

    Returns
    -------
    area_values : ndarray, shape (V,), float64
        Per-vertex area contribution.
    total_area  : float
        Sum of all face areas.
    """
    cdef PolygonsMesh mesh = _ensure_mesh(vertices, faces)
    cdef polygons_struct *poly = mesh.ptr()
    cdef int n = poly.n_points

    cdef cnp.ndarray[cnp.float64_t, ndim=1] area = np.zeros(n, dtype=np.float64)
    cdef double total = C.get_area_of_points(poly, <double *>area.data)
    return area, total


# ===================================================================
# Euler characteristic
# ===================================================================
def euler_characteristic(vertices, faces, bint verbose=False):
    """
    Compute the Euler characteristic of a triangular mesh.

    Parameters
    ----------
    vertices : array_like, shape (V, 3)
    faces    : array_like, shape (F, 3)
    verbose  : bool

    Returns
    -------
    int
        Euler characteristic (2 for a closed sphere).
    """
    cdef PolygonsMesh mesh = _ensure_mesh(vertices, faces)
    return C.euler_characteristic(mesh.ptr(), 1 if verbose else 0)


# ===================================================================
# Point distance  (mirrors CAT_SurfDistance)
# ===================================================================
def point_distance(vertices1, faces1, vertices2, faces2,
                   bint symmetric=False, double max_dist=0.0):
    """
    Compute per-vertex linked (exact) distance between two meshes.

    Mirrors the ``-link`` mode of the ``CAT_SurfDistance`` CLI tool.

    Parameters
    ----------
    vertices1, faces1 : mesh 1 arrays
    vertices2, faces2 : mesh 2 arrays
    symmetric : bool
        If True, compute symmetric (bidirectional) distance.
        Default ``False`` to match the CLI.
    max_dist : float
        If > 0, clip distances at this value (mirrors ``-max``).

    Returns
    -------
    distances : ndarray, shape (V1,), float64
    max_distance : float
        Maximum distance returned by the C routine (unclipped).
    """
    cdef PolygonsMesh m1 = _ensure_mesh(vertices1, faces1)
    cdef PolygonsMesh m2 = _ensure_mesh(vertices2, faces2)
    cdef polygons_struct *p1 = m1.ptr()
    cdef int n = p1.n_points

    cdef cnp.ndarray[cnp.float64_t, ndim=1] dist = np.zeros(n, dtype=np.float64)
    cdef double maxd = C.compute_point_distance(
        p1, m2.ptr(), <double *>dist.data, 1 if symmetric else 0,
    )
    if max_dist > 0.0:
        np.minimum(dist, max_dist, out=dist)
    return dist, maxd


def point_distance_mean(vertices1, faces1, vertices2, faces2,
                        bint symmetric=False, double max_dist=0.0):
    """
    Compute mean (Freesurfer Tfs) per-vertex distances between two meshes.

    Mirrors the default ``-mean`` mode of ``CAT_SurfDistance``.

    Parameters
    ----------
    symmetric : bool
        Default ``False`` to match the CLI.
    max_dist : float
        If > 0, clip distances at this value (mirrors ``-max``).

    Returns
    -------
    distances : ndarray, shape (V1,), float64
    mean_distance : float
    """
    cdef PolygonsMesh m1 = _ensure_mesh(vertices1, faces1)
    cdef PolygonsMesh m2 = _ensure_mesh(vertices2, faces2)
    cdef polygons_struct *p1 = m1.ptr()
    cdef int n = p1.n_points

    cdef cnp.ndarray[cnp.float64_t, ndim=1] dist = np.zeros(n, dtype=np.float64)
    cdef double meand = C.compute_point_distance_mean(
        p1, m2.ptr(), <double *>dist.data, 1 if symmetric else 0,
    )
    if max_dist > 0.0:
        np.minimum(dist, max_dist, out=dist)
    return dist, meand


def hausdorff_distance(vertices1, faces1, vertices2, faces2,
                       bint symmetric=True, bint exact=False):
    """
    Compute Hausdorff distance between two meshes.

    Returns
    -------
    distances : ndarray, shape (V1,), float64
    hausdorff : float
    """
    cdef PolygonsMesh m1 = _ensure_mesh(vertices1, faces1)
    cdef PolygonsMesh m2 = _ensure_mesh(vertices2, faces2)
    cdef polygons_struct *p1 = m1.ptr()
    cdef int n = p1.n_points

    cdef cnp.ndarray[cnp.float64_t, ndim=1] dist = np.zeros(n, dtype=np.float64)
    cdef double hd
    if exact:
        hd = C.compute_exact_hausdorff(
            p1, m2.ptr(), <double *>dist.data, 1 if symmetric else 0)
    else:
        hd = C.compute_point_hausdorff(
            p1, m2.ptr(), <double *>dist.data, 1 if symmetric else 0)
    return dist, hd


# ===================================================================
# Heat-kernel smoothing of per-vertex values
# ===================================================================
def smooth_heatkernel(vertices, faces, values, double fwhm):
    """
    Smooth per-vertex scalar data using a heat-kernel approach.

    Parameters
    ----------
    vertices : array_like, shape (V, 3)
    faces    : array_like, shape (F, 3)
    values   : array_like, shape (V,), float64
    fwhm     : float
        Full-width at half-maximum of the smoothing kernel in mm.

    Returns
    -------
    smoothed : ndarray, shape (V,), float64
    """
    cdef PolygonsMesh mesh = _ensure_mesh(vertices, faces)
    cdef polygons_struct *poly = mesh.ptr()
    cdef int n = poly.n_points

    vals = np.ascontiguousarray(values, dtype=np.float64)
    if vals.shape[0] != n:
        raise ValueError(
            f"values length ({vals.shape[0]}) != number of vertices ({n})")

    # smooth_heatkernel modifies the array in-place
    cdef cnp.ndarray[cnp.float64_t, ndim=1] out = vals.copy()
    C.smooth_heatkernel(poly, <double *>out.data, fwhm)
    return out


# ===================================================================
# Curvature
# ===================================================================
def smoothed_curvatures(vertices, faces, double fwhm=0.0, int n_iter=0):
    """
    Compute smoothed mean curvature for each vertex.

    Parameters
    ----------
    vertices : array_like, shape (V, 3)
    faces    : array_like, shape (F, 3)
    fwhm     : float
        Smoothing FWHM (0 = no smoothing).
    n_iter   : int
        Number of smoothing iterations.

    Returns
    -------
    curvatures : ndarray, shape (V,), float64
    """
    cdef PolygonsMesh mesh = _ensure_mesh(vertices, faces)
    cdef polygons_struct *poly = mesh.ptr()
    cdef int n = poly.n_points

    cdef cnp.ndarray[cnp.float64_t, ndim=1] curv = np.zeros(n, dtype=np.float64)
    C.get_smoothed_curvatures(poly, <double *>curv.data, fwhm, n_iter)
    return curv


# ===================================================================
# Per-vertex curvature  (mirrors CAT_SurfCurvature)
# ===================================================================
def surf_curvature(vertices, faces, int curvtype=0, double fwhm=0.0,
                   bint use_abs_values=False, bint invert_values=False):
    """
    Compute per-vertex curvature.

    Mirrors ``CAT_SurfCurvature``.  ``curvtype`` is the same set of
    measures the CLI exposes:

    ===========  ==========================================================
    ``curvtype`` Description
    ===========  ==========================================================
    0            Mean curvature averaged over 3 mm, in degrees: (k1+k2)/2
    1            Gaussian curvature  k1*k2
    2            Curvedness          sqrt(0.5 * (k1^2 + k2^2))
    3            Shape index         atan((k1+k2) / (k2-k1))
    4            Mean curvature in radians: (k1+k2)/2
    5            Sulcal-depth-like estimator
    6            Bending energy      k1^2 + k2^2
    7            Sharpness           (k1 - k2)^2
    8            Folding index       |k1| * (|k1| - |k2|)
    9            Minimum curvature   k2
    10           Maximum curvature   k1
    11           Inflation-based sulcal depth (FreeSurfer sulc-style)
    > 11         Depth potential with alpha = 1 / ``curvtype``
                 (e.g. 650 = recommended)
    ===========  ==========================================================

    Parameters
    ----------
    vertices : array_like, shape (V, 3)
    faces    : array_like, shape (F, 3)
    curvtype : int
        Curvature measure to compute (default 0).
    fwhm : float
        Heat-kernel smoothing FWHM in mm (default 0 = no smoothing).
    use_abs_values : bool
        Take absolute value of every output (default False).
    invert_values : bool
        Negate every output (default False).

    Returns
    -------
    curvatures : ndarray, shape (V,), float64
    """
    cdef PolygonsMesh mesh = _ensure_mesh(vertices, faces)
    cdef polygons_struct *poly = mesh.ptr()
    cdef int n = poly.n_points

    cdef int *n_neighbours = NULL
    cdef int **neighbours  = NULL
    C.get_all_polygon_point_neighbours(poly, &n_neighbours, &neighbours)

    cdef double distance = 3.0 if curvtype == 0 else 0.0

    cdef cnp.ndarray[cnp.float64_t, ndim=1] curv = np.zeros(n, dtype=np.float64)
    C.get_polygon_vertex_curvatures_cg(
        poly, n_neighbours, neighbours,
        distance, curvtype, <double *>curv.data)

    # Mirror the CLI's range clip for the unit-normalised measures.
    if (0 < curvtype < 4) or curvtype == 10:
        np.clip(curv, -1.0, 1.0, out=curv)

    if use_abs_values:
        np.abs(curv, out=curv)
    if invert_values:
        np.negative(curv, out=curv)

    if fwhm > 0.0:
        C.smooth_heatkernel(poly, <double *>curv.data, fwhm)

    return curv


# ===================================================================
# Sulcus depth
# ===================================================================
def sulcus_depth(vertices, faces):
    """
    Compute sulcus depth for each vertex.

    Parameters
    ----------
    vertices : array_like, shape (V, 3)
    faces    : array_like, shape (F, 3)

    Returns
    -------
    depth : ndarray, shape (V,), float64
    """
    cdef PolygonsMesh mesh = _ensure_mesh(vertices, faces)
    cdef polygons_struct *poly = mesh.ptr()
    cdef int n = poly.n_points

    cdef cnp.ndarray[cnp.float64_t, ndim=1] depth = np.zeros(n, dtype=np.float64)
    C.compute_sulcus_depth(poly, <double *>depth.data)
    return depth


# ===================================================================
# Mesh reduction (quadric decimation)  (mirrors CAT_SurfReduce)
# ===================================================================
def reduce_mesh(vertices, faces, int target_faces,
                double aggressiveness=7.0, bint preserve_sharp=True,
                bint verbose=False):
    """
    Reduce the number of faces in a triangular mesh using
    Quadric Error Metrics (QEM) decimation.

    Parameters
    ----------
    vertices : array_like, shape (V, 3)
    faces    : array_like, shape (F, 3)
    target_faces : int
        Desired number of triangles after decimation.
    aggressiveness : float
        Larger values produce stronger decimation (default 7.0).
    preserve_sharp : bool
        Prevent collapses across sharp features (default True).
    verbose  : bool

    Returns
    -------
    new_vertices : ndarray, shape (V', 3), float64
    new_faces    : ndarray, shape (F', 3), int32
    """
    cdef PolygonsMesh mesh = _ensure_mesh(vertices, faces)
    cdef int rc = C.reduce_mesh_quadrics(
        mesh.ptr(), target_faces, aggressiveness,
        1 if preserve_sharp else 0,
        1 if verbose else 0,
    )
    if rc != 0:
        raise RuntimeError(f"reduce_mesh_quadrics returned error code {rc}")

    from cat_surf._convert import polygons_to_arrays
    return polygons_to_arrays(mesh)


# ===================================================================
# Sphere radius
# ===================================================================
def sphere_radius(vertices, faces):
    """
    Compute the average radius of a spherical mesh.

    Parameters
    ----------
    vertices : array_like, shape (V, 3)
    faces    : array_like, shape (F, 3)

    Returns
    -------
    float
        Mean radius.
    """
    cdef PolygonsMesh mesh = _ensure_mesh(vertices, faces)
    return C.get_sphere_radius(mesh.ptr())


# ===================================================================
# Laplacian smoothing of the mesh itself
# ===================================================================
def smooth_mesh(vertices, faces, int iterations=10,
                double alpha=0.5, double beta=-0.53):
    """
    Taubin-style Laplacian smoothing of the mesh geometry.

    Parameters
    ----------
    vertices : array_like, shape (V, 3)
    faces    : array_like, shape (F, 3)
    iterations : int
    alpha, beta : float
        Shrink/inflate parameters.

    Returns
    -------
    smoothed_vertices : ndarray, shape (V, 3), float64
    faces             : ndarray, shape (F, 3), int32
    """
    cdef PolygonsMesh mesh = _ensure_mesh(vertices, faces)
    C.smooth_laplacian(mesh.ptr(), iterations, alpha, beta)
    from cat_surf._convert import polygons_to_arrays
    return polygons_to_arrays(mesh)


# ===================================================================
# Sphere mapping  (mirrors CAT_Surf2Sphere)
# ===================================================================
def surf_to_sphere(vertices, faces, int stop_at=5, bint verbose=False):
    """
    Map a surface to a sphere using iterative inflation.

    Parameters
    ----------
    vertices : array_like, shape (V, 3)
    faces    : array_like, shape (F, 3)
    stop_at  : int
        Stopping stage (1-5+).  Higher → finer resolution.
    verbose  : bool

    Returns
    -------
    sphere_vertices : ndarray, shape (V, 3), float64
    sphere_faces    : ndarray, shape (F, 3), int32
    """
    cdef PolygonsMesh mesh = _ensure_mesh(vertices, faces)
    C.surf_to_sphere(mesh.ptr(), stop_at, 1 if verbose else 0)
    return polygons_to_arrays(mesh)


# ===================================================================
# Area normalized to sphere  (mirrors CAT_SurfArea -sphere)
# ===================================================================
def get_area_normalized(vertices, faces, sphere_vertices, sphere_faces):
    """
    Compute per-vertex area normalized against a reference sphere.

    Parameters
    ----------
    vertices, faces : mesh arrays
    sphere_vertices, sphere_faces : sphere mesh arrays (same topology)

    Returns
    -------
    area_values : ndarray, shape (V,), float64
    total_area  : float
    """
    cdef PolygonsMesh mesh = _ensure_mesh(vertices, faces)
    cdef PolygonsMesh sph  = _ensure_mesh(sphere_vertices, sphere_faces)
    cdef polygons_struct *poly = mesh.ptr()
    cdef int n = poly.n_points

    cdef cnp.ndarray[cnp.float64_t, ndim=1] area = np.zeros(n, dtype=np.float64)
    cdef double total = C.get_area_of_points_normalized_to_sphere(
        poly, sph.ptr(), <double *>area.data)
    return area, total


# ===================================================================
# Surface average  (mirrors CAT_SurfAverage)
# ===================================================================
def surf_average(list surfaces, bint return_rms=False):
    """
    Compute the vertex-wise average of multiple surfaces.

    All surfaces must share the same topology (same faces, same
    vertex count).  Mirrors ``CAT_SurfAverage -avg`` (and ``-rms``).

    Parameters
    ----------
    surfaces : list of (vertices, faces) tuples
        Each element is a pair of arrays ``(V_array, F_array)``.
    return_rms : bool
        If True, also return per-vertex root-mean-square deviation
        (sample standard deviation of the 3-D vertex positions).

    Returns
    -------
    avg_vertices : ndarray, shape (V, 3), float64
    faces        : ndarray, shape (F, 3), int32
    rms          : ndarray, shape (V,), float64
        Only returned if ``return_rms`` is True.
    """
    if len(surfaces) == 0:
        raise ValueError("surfaces list must not be empty")

    cdef int n_sets = len(surfaces)
    v0, f0 = surfaces[0]
    avg = np.ascontiguousarray(v0, dtype=np.float64).copy()
    faces_out = np.ascontiguousarray(f0, dtype=np.int32)
    cdef int n_pts = avg.shape[0]

    sqr = None
    if return_rms:
        sqr = avg * avg

    cdef int i
    for i in range(1, n_sets):
        vi, _ = surfaces[i]
        vi = np.ascontiguousarray(vi, dtype=np.float64)
        if vi.shape[0] != n_pts:
            raise ValueError(
                f"Surface {i} has {vi.shape[0]} vertices, expected {n_pts}")
        avg += vi
        if return_rms:
            sqr += vi * vi

    avg /= n_sets

    if return_rms:
        if n_sets < 2:
            raise ValueError(
                "At least 2 surfaces are required to compute RMS deviation")
        # Mirrors CAT_SurfAverage RMS:
        # rms_i = sqrt( (sqr_i / (N-1)) - (N/(N-1)) * avg_i**2 )  per axis,
        # then sum over axes inside sqrt.
        var = sqr / (n_sets - 1.0) - (n_sets / (n_sets - 1.0)) * avg * avg
        rms = np.sqrt(np.maximum(var.sum(axis=1), 0.0))
        return avg, faces_out, rms

    return avg, faces_out


# ===================================================================
# Correct thickness for folding  (mirrors CAT_SurfCorrectThicknessFolding)
# ===================================================================
def correct_thickness_folding(vertices, faces, thickness,
                              double slope=0.0, double max_dist=6.0):
    """
    Correct cortical thickness values for folding-related bias.

    Implements the Demirci et al. (2025) correction (doi.org/10.1101/
    2025.05.03.651968).  Mirrors the ``CAT_SurfCorrectThicknessFolding``
    CLI tool.

    Parameters
    ----------
    vertices : array_like, shape (V, 3)
    faces    : array_like, shape (F, 3)
    thickness : array_like, shape (V,), float64
        Uncorrected thickness values.
    slope : float
        Thickness-based weighting strength (default 0.0 = disabled).
        Larger thickness values are corrected more strongly when
        ``slope > 0``.  Applied only for positive mean curvature.
    max_dist : float
        Upper clip for output values (default 6.0 mm).  Values exceeding
        this are clipped.  Set to 0 or +inf to disable.

    Returns
    -------
    corrected : ndarray, shape (V,), float64
    """
    cdef PolygonsMesh mesh = _ensure_mesh(vertices, faces)
    cdef polygons_struct *poly = mesh.ptr()
    cdef int n = poly.n_points

    vals = np.ascontiguousarray(thickness, dtype=np.float64)
    if vals.shape[0] != n:
        raise ValueError(
            f"thickness length ({vals.shape[0]}) != vertices ({n})")

    cdef cnp.ndarray[cnp.float64_t, ndim=1] out = vals.copy()
    cdef Status st = C.CAT_CorrectThicknessFoldingWeighted(
        poly, n, <double *>out.data, slope)
    if st != OK:
        raise RuntimeError("CAT_CorrectThicknessFolding failed")

    if max_dist > 0.0:
        np.clip(out, 1e-10, max_dist, out=out)
    return out


# ===================================================================
# Remove self-intersections  (mirrors CAT_SurfRemoveIntersections)
# ===================================================================
def remove_intersections(vertices, faces, int max_iters=10,
                         int inner_loops=3, bint fill_holes=True,
                         bint verbose=False):
    """
    Remove self-intersections from a triangle mesh (MeshFix).

    Parameters
    ----------
    vertices : array_like, shape (V, 3)
    faces    : array_like, shape (F, 3)
    max_iters : int
    inner_loops : int
    fill_holes : bool
    verbose : bool

    Returns
    -------
    new_vertices : ndarray, shape (V', 3), float64
    new_faces    : ndarray, shape (F', 3), int32
    """
    cdef PolygonsMesh mesh = _ensure_mesh(vertices, faces)

    cdef C.CAT_MeshCleanOptions opts
    C.CAT_MeshCleanOptionsInit(&opts)
    opts.max_iters = max_iters
    opts.inner_loops = inner_loops
    opts.fill_holes = 1 if fill_holes else 0
    opts.verbose = 1 if verbose else 0

    cdef int rc = C.CAT_SurfMeshClean(mesh.ptr(), &opts)
    if rc < 0:
        raise RuntimeError(f"CAT_SurfMeshClean returned error code {rc}")

    return polygons_to_arrays(mesh)


def count_intersections(vertices, faces):
    """
    Count the number of self-intersecting triangles.

    Parameters
    ----------
    vertices : array_like, shape (V, 3)
    faces    : array_like, shape (F, 3)

    Returns
    -------
    int
        Number of self-intersections.
    """
    cdef PolygonsMesh mesh = _ensure_mesh(vertices, faces)
    return C.CAT_SurfCountIntersections(mesh.ptr())


# ===================================================================
# Spherical resampling  (mirrors CAT_SurfResample)
# ===================================================================
def resample_to_sphere(vertices, faces, sphere_vertices, sphere_faces,
                       target_sphere_vertices, target_sphere_faces,
                       values=None, bint label_interpolation=False,
                       bint areal_interpolation=False):
    """
    Resample a surface (and optional per-vertex data) from one
    spherical parameterization to another.

    Parameters
    ----------
    vertices, faces : source mesh arrays
    sphere_vertices, sphere_faces : source sphere arrays
    target_sphere_vertices, target_sphere_faces : target sphere arrays
    values : array_like, shape (V_src,), float64, optional
        Per-vertex data to resample.
    label_interpolation : bool
        Use nearest-neighbor (label) interpolation.
    areal_interpolation : bool
        Use areal-weighted interpolation.

    Returns
    -------
    new_vertices : ndarray, shape (V', 3), float64
    new_faces    : ndarray, shape (F', 3), int32
    new_values   : ndarray or None
        Resampled values (only if ``values`` was provided).
    """
    cdef PolygonsMesh src = _ensure_mesh(vertices, faces)
    cdef PolygonsMesh src_sph = _ensure_mesh(sphere_vertices, sphere_faces)
    cdef PolygonsMesh tgt_sph = _ensure_mesh(
        target_sphere_vertices, target_sphere_faces)

    cdef int n_tgt = tgt_sph.ptr().n_points
    cdef double *in_vals = NULL
    cdef double *out_vals = NULL
    cdef cnp.ndarray[cnp.float64_t, ndim=1] in_arr
    cdef cnp.ndarray[cnp.float64_t, ndim=1] out_arr

    if values is not None:
        in_arr = np.ascontiguousarray(values, dtype=np.float64)
        in_vals = <double *>in_arr.data
        out_arr = np.zeros(n_tgt, dtype=np.float64)
        out_vals = <double *>out_arr.data

    cdef object_struct **result = C.resample_surface_to_target_sphere(
        src.ptr(), src_sph.ptr(), tgt_sph.ptr(),
        in_vals, out_vals,
        1 if label_interpolation else 0,
        1 if areal_interpolation else 0,
    )

    if result == NULL:
        raise RuntimeError("resample_surface_to_target_sphere failed")

    cdef PolygonsMesh out_mesh = _wrap_object(result[0], owns=True)
    free(result)  # free the array, not the object

    new_v, new_f = polygons_to_arrays(out_mesh)
    if values is not None:
        return new_v, new_f, out_arr
    return new_v, new_f, None


# ===================================================================
# Annotation resampling — file-to-file, mirrors:
#   CAT_SurfResample -label <src_surf> <src_sphere> <tgt_sphere> NULL \
#                            <annot_in> <annot_out>
# Used by T1Prep's atlas resampling step.
# ===================================================================
def resample_annot(str source_surface_file,
                   str source_sphere_file,
                   str target_sphere_file,
                   str annot_in_file,
                   str annot_out_file):
    """
    Resample a FreeSurfer ``.annot`` label file from one cortical
    spherical parameterisation to another.

    Mirrors the ``-label <annot_in> <annot_out>`` mode of
    ``CAT_SurfResample``.  Inputs and outputs are file paths so the
    full ATABLE round-trip stays inside libCAT.

    Parameters
    ----------
    source_surface_file : str
        Atlas surface (its topology defines the label vertex count).
    source_sphere_file : str
        Atlas sphere corresponding to ``source_surface_file``.
    target_sphere_file : str
        Subject sphere to resample onto.
    annot_in_file : str
        FreeSurfer ``.annot`` file matching ``source_surface_file``.
    annot_out_file : str
        Output ``.annot`` written with the resampled per-vertex labels.

    Returns
    -------
    n_target_vertices : int
        Number of vertices in the resampled output (= target sphere V).
    """
    cdef bytes b_src      = source_surface_file.encode("utf-8")
    cdef bytes b_src_sph  = source_sphere_file.encode("utf-8")
    cdef bytes b_tgt_sph  = target_sphere_file.encode("utf-8")
    cdef bytes b_annot_in = annot_in_file.encode("utf-8")
    cdef bytes b_annot_o  = annot_out_file.encode("utf-8")

    cdef File_formats fmt
    cdef int n_objects = 0
    cdef object_struct **src_objs       = NULL
    cdef object_struct **src_sph_objs   = NULL
    cdef object_struct **tgt_sph_objs   = NULL
    cdef object_struct **resampled_objs = NULL

    cdef int n_values = 0
    cdef int *annot_in = NULL
    cdef int n_table = 0
    cdef C.ATABLE *atable = NULL

    cdef double *in_vals = NULL
    cdef double *out_vals = NULL
    cdef int *annot_out = NULL

    cdef polygons_struct *src_poly = NULL
    cdef polygons_struct *src_sph_poly = NULL
    cdef polygons_struct *tgt_sph_poly = NULL
    cdef int i = 0
    cdef int n_tgt = 0

    try:
        if C.input_graphics_any_format(b_src,     &fmt, &n_objects, &src_objs)     != OK:
            raise IOError(f"Cannot read source surface: {source_surface_file}")
        if n_objects != 1 or get_object_type(src_objs[0]) != POLYGONS:
            raise ValueError(f"{source_surface_file}: must contain exactly one polygons object")
        if C.input_graphics_any_format(b_src_sph, &fmt, &n_objects, &src_sph_objs) != OK:
            raise IOError(f"Cannot read source sphere: {source_sphere_file}")
        if n_objects != 1 or get_object_type(src_sph_objs[0]) != POLYGONS:
            raise ValueError(f"{source_sphere_file}: must contain exactly one polygons object")
        if C.input_graphics_any_format(b_tgt_sph, &fmt, &n_objects, &tgt_sph_objs) != OK:
            raise IOError(f"Cannot read target sphere: {target_sphere_file}")
        if n_objects != 1 or get_object_type(tgt_sph_objs[0]) != POLYGONS:
            raise ValueError(f"{target_sphere_file}: must contain exactly one polygons object")

        src_poly     = get_polygons_ptr(src_objs[0])
        src_sph_poly = get_polygons_ptr(src_sph_objs[0])
        tgt_sph_poly = get_polygons_ptr(tgt_sph_objs[0])
        n_tgt = tgt_sph_poly.n_points

        if C.read_annotation_table(b_annot_in, &n_values, &annot_in,
                                   &n_table, &atable) != 0:
            raise IOError(f"Cannot read annotation table: {annot_in_file}")

        in_vals  = <double *>malloc(sizeof(double) * <size_t>src_sph_poly.n_points)
        out_vals = <double *>malloc(sizeof(double) * <size_t>n_tgt)
        if in_vals == NULL or out_vals == NULL:
            raise MemoryError()
        for i in range(min(n_values, src_sph_poly.n_points)):
            in_vals[i] = <double>annot_in[i]

        resampled_objs = C.resample_surface_to_target_sphere(
            src_poly, src_sph_poly, tgt_sph_poly,
            in_vals, out_vals,
            1,   # label_interpolation (annot => categorical)
            0,   # areal_interpolation
        )
        if resampled_objs == NULL:
            raise RuntimeError("resample_surface_to_target_sphere failed")

        annot_out = <int *>malloc(sizeof(int) * <size_t>n_tgt)
        if annot_out == NULL:
            raise MemoryError()
        for i in range(n_tgt):
            # Round-to-nearest, matching the CLI's `(int)round(...)`.
            annot_out[i] = <int>(out_vals[i] + (0.5 if out_vals[i] >= 0.0 else -0.5))

        if C.write_annotation_table(b_annot_o, n_tgt, annot_out,
                                    n_table, atable) != 0:
            raise IOError(f"Cannot write annotation table: {annot_out_file}")

        return n_tgt

    finally:
        if annot_in is not NULL:
            free(annot_in)
        if annot_out is not NULL:
            free(annot_out)
        if in_vals is not NULL:
            free(in_vals)
        if out_vals is not NULL:
            free(out_vals)
        if resampled_objs is not NULL:
            delete_object(resampled_objs[0])
            free(resampled_objs)
        if src_objs is not NULL:
            delete_object(src_objs[0])
            free(src_objs)
        if src_sph_objs is not NULL:
            delete_object(src_sph_objs[0])
            free(src_sph_objs)
        if tgt_sph_objs is not NULL:
            delete_object(tgt_sph_objs[0])
            free(tgt_sph_objs)


# ===================================================================
# Surface deformation  (mirrors CAT_SurfDeform)
# ===================================================================
def surf_deform(vertices, faces, volume,
                double w1=0.0, double w2=0.2, double w3=1.2,
                double sigma=0.2, float isovalue=0.5,
                int iterations=75, bint remove_intersect=False,
                bint verbose=False):
    """
    Deform a surface toward a volume isovalue.

    Parameters
    ----------
    vertices, faces : mesh arrays
    volume : str | (ndarray, affine) | nibabel-image-like
        Source volume.  Accepts a NIfTI file path, an ``(array, affine)``
        tuple, or any object with ``.affine`` and ``.get_fdata()`` (e.g.
        ``nibabel.Nifti1Image``).  Datatype is auto-converted to float32.
    w1 : float
        Internal smoothness weight (default 0.0).
    w2 : float
        Gradient alignment weight (default 0.2).
    w3 : float
        Balloon force weight (default 1.2).
    sigma : float
        Displacement smoothing sigma (default 0.2).
    isovalue : float
        Target isovalue in the volume (default 0.5).
    iterations : int
        Number of deformation iterations (default 75).
    remove_intersect : bool
        Remove self-intersections via MeshFix at the end (default False).
    verbose : bool

    Returns
    -------
    new_vertices : ndarray, shape (V, 3), float64
    faces        : ndarray, shape (F, 3), int32
    """
    cdef VolumeHandle vh = open_volume(volume)
    cdef PolygonsMesh mesh = _ensure_mesh(vertices, faces)

    cdef double w[3]
    w[0] = w1; w[1] = w2; w[2] = w3

    try:
        C.surf_deform(mesh.ptr(), vh.data, vh.nii,
                      w, sigma, isovalue, iterations,
                      1 if remove_intersect else 0,
                      1 if verbose else 0)
    finally:
        vh.close()
    return polygons_to_arrays(mesh)


# ===================================================================
# Pial / white surface estimation  (mirrors CAT_Surf2PialWhite)
# ===================================================================
def surf_to_pial_white(vertices, faces, thickness,
                       label,
                       double w1=0.05, double w2=0.05, double w3=0.05,
                       double sigma=0.2, int iterations=100,
                       int gradient_iterations=0, int method=0,
                       bint verbose=False):
    """
    Estimate pial and white matter surfaces from a central surface.

    Mirrors the ``CAT_Surf2PialWhite`` CLI tool.

    Parameters
    ----------
    vertices, faces : central surface mesh arrays
    thickness : array_like, shape (V,), float64
        Cortical thickness per vertex.
    label : str | (ndarray, affine) | nibabel-image-like
        Label volume.  See :func:`cat_surf._volume.open_volume` for
        accepted forms.  Auto-converted to float32.
    w1 : float
        Internal smoothness weight (default 0.05).
    w2 : float
        Gradient alignment weight (default 0.05).
    w3 : float
        Balloon force weight (default 0.05).
    sigma : float
        Displacement smoothing sigma (default 0.2).
    iterations : int
        Number of deformation iterations (default 100).
    gradient_iterations : int
        Number of gradient refinement iterations (default 0 = disabled).
    method : int
        0 = deformation (default), 1 = ADE,
        2 = deformation for pial + ADE for white.
    verbose : bool

    Returns
    -------
    pial_vertices  : ndarray, shape (V, 3), float64
    pial_faces     : ndarray, shape (F, 3), int32
    white_vertices : ndarray, shape (V, 3), float64
    white_faces    : ndarray, shape (F, 3), int32
    """
    cdef VolumeHandle vh = open_volume(label)
    cdef PolygonsMesh mesh = _ensure_mesh(vertices, faces)
    cdef polygons_struct *poly = mesh.ptr()
    cdef int n = poly.n_points

    thick = np.ascontiguousarray(thickness, dtype=np.float64)
    if thick.shape[0] != n:
        vh.close()
        raise ValueError(
            f"thickness length ({thick.shape[0]}) != vertices ({n})")
    cdef cnp.ndarray[cnp.float64_t, ndim=1] thick_arr = thick

    cdef object_struct *pial_obj = create_object(POLYGONS)
    cdef object_struct *white_obj = create_object(POLYGONS)
    cdef polygons_struct *pial_poly = get_polygons_ptr(pial_obj)
    cdef polygons_struct *white_poly = get_polygons_ptr(white_obj)

    cdef C.CAT_PialWhiteOptions opts
    C.CAT_PialWhiteOptionsInit(&opts)
    opts.w1 = w1; opts.w2 = w2; opts.w3 = w3
    opts.sigma = sigma
    opts.iterations = iterations
    opts.gradient_iterations = gradient_iterations
    opts.method = method
    opts.verbose = 1 if verbose else 0

    cdef int rc = C.CAT_SurfEstimatePialWhite(
        poly, <double *>thick_arr.data, vh.data, vh.nii,
        pial_poly, white_poly, &opts)
    vh.close()

    if rc != 0:
        delete_object(pial_obj)
        delete_object(white_obj)
        raise RuntimeError(
            f"CAT_SurfEstimatePialWhite returned error code {rc}")

    cdef PolygonsMesh pial_mesh = _wrap_object(pial_obj, owns=True)
    cdef PolygonsMesh white_mesh = _wrap_object(white_obj, owns=True)

    pv, pf = polygons_to_arrays(pial_mesh)
    wv, wf = polygons_to_arrays(white_mesh)
    return pv, pf, wv, wf


# ===================================================================
# Central → pial distance  (part of CAT_SurfDistance -thickness)
# ===================================================================
def central_to_pial(vertices, faces, thickness, double sigma=1.0,
                    int iterations=300, bint check_intersect=True,
                    bint verbose=False):
    """
    Generate a pial surface from a central surface and thickness values.

    Parameters
    ----------
    vertices, faces : central surface mesh arrays
    thickness : array_like, shape (V,), float64
    sigma : float
        Smoothing sigma.
    iterations : int
        Deformation iterations.
    check_intersect : bool
    verbose : bool

    Returns
    -------
    pial_vertices : ndarray, shape (V, 3), float64
    pial_faces    : ndarray, shape (F, 3), int32
    """
    cdef PolygonsMesh mesh = _ensure_mesh(vertices, faces)
    cdef polygons_struct *poly = mesh.ptr()
    cdef int n = poly.n_points

    thick = np.ascontiguousarray(thickness, dtype=np.float64)
    if thick.shape[0] != n:
        raise ValueError(
            f"thickness length ({thick.shape[0]}) != vertices ({n})")
    cdef cnp.ndarray[cnp.float64_t, ndim=1] thick_arr = thick

    # extents: use thickness as extents (half thickness each side)
    cdef cnp.ndarray[cnp.float64_t, ndim=1] extents = thick_arr * 0.5

    cdef object_struct **result = C.central_to_new_pial(
        poly, <double *>thick_arr.data, <double *>extents.data,
        1 if check_intersect else 0, sigma, iterations,
        1 if verbose else 0)

    if result == NULL:
        raise RuntimeError("central_to_new_pial failed")

    cdef PolygonsMesh out_mesh = _wrap_object(result[0], owns=True)
    free(result)

    return polygons_to_arrays(out_mesh)


# ===================================================================
# Single-unit ROI statistics  (building block for CAT_Surf2ROIMulti)
# ===================================================================
def surf2roi_unit(str src_sphere_file, str trg_sphere_file,
                  str annot_file, str vals_file):
    """
    Resample annotation labels onto a target sphere and compute
    per-ROI means from a values array.

    This is the single-unit building block for ``CAT_Surf2ROIMulti``.
    Call it once per hemisphere / atlas and aggregate the results with
    :func:`cat_surf.surf2roi_multi`.

    Parameters
    ----------
    src_sphere_file : str
        Source sphere whose vertex count matches the annotation.
    trg_sphere_file : str
        Target sphere to resample onto.
    annot_file : str
        FreeSurfer ``.annot`` annotation file on the source sphere.
    vals_file : str
        Per-vertex values file defined on the *target* sphere.

    Returns
    -------
    stats : list of dict
        Each dict has keys ``"id"`` (int), ``"name"`` (str),
        ``"sum"`` (float), and ``"n"`` (int).  The per-ROI mean is
        ``sum / n``.
    """
    cdef bytes b_src_sph  = src_sphere_file.encode("utf-8")
    cdef bytes b_trg_sph  = trg_sphere_file.encode("utf-8")
    cdef bytes b_annot    = annot_file.encode("utf-8")
    cdef bytes b_vals     = vals_file.encode("utf-8")

    cdef File_formats fmt
    cdef int n_objects = 0
    cdef object_struct **src_sph_objs = NULL
    cdef object_struct **trg_sph_objs = NULL
    cdef polygons_struct *src_sph_poly = NULL
    cdef polygons_struct *trg_sph_poly = NULL

    cdef int n_values = 0
    cdef int *annot_in = NULL
    cdef int n_table = 0
    cdef C.ATABLE *atable = NULL
    cdef int *labels_trg = NULL

    cdef int n_vals = 0
    cdef double *vals = NULL

    cdef C.CAT_ROIStat *stats = NULL
    cdef int n_stats = 0
    cdef int i
    cdef Status st

    try:
        if C.input_graphics_any_format(b_src_sph, &fmt, &n_objects,
                                        &src_sph_objs) != OK:
            raise IOError(f"Cannot read source sphere: {src_sphere_file}")
        if n_objects != 1 or get_object_type(src_sph_objs[0]) != POLYGONS:
            raise ValueError(
                f"{src_sphere_file}: must contain exactly one polygons object")
        src_sph_poly = get_polygons_ptr(src_sph_objs[0])

        if C.input_graphics_any_format(b_trg_sph, &fmt, &n_objects,
                                        &trg_sph_objs) != OK:
            raise IOError(f"Cannot read target sphere: {trg_sphere_file}")
        if n_objects != 1 or get_object_type(trg_sph_objs[0]) != POLYGONS:
            raise ValueError(
                f"{trg_sphere_file}: must contain exactly one polygons object")
        trg_sph_poly = get_polygons_ptr(trg_sph_objs[0])

        if C.read_annotation_table(b_annot, &n_values, &annot_in,
                                   &n_table, &atable) != 0:
            raise IOError(f"Cannot read annotation: {annot_file}")

        # Resample integer annotation labels src → trg
        st = C.CAT_ResampleAnnotationLabels(
            src_sph_poly, trg_sph_poly, annot_in, &labels_trg)
        if st != OK:
            raise RuntimeError("CAT_ResampleAnnotationLabels failed")

        # Read per-vertex values on target sphere
        if C.input_values_any_format(b_vals, &n_vals, &vals) != OK:
            raise IOError(f"Cannot read values: {vals_file}")
        if n_vals != trg_sph_poly.n_points:
            raise ValueError(
                f"#values ({n_vals}) != #points in target sphere "
                f"({trg_sph_poly.n_points})")

        # Compute per-label means
        st = C.CAT_ComputeROIMeansFromLabels(
            labels_trg, vals, trg_sph_poly.n_points,
            atable, n_table, &stats, &n_stats)
        if st != OK:
            raise RuntimeError("CAT_ComputeROIMeansFromLabels failed")

        result = []
        for i in range(n_stats):
            nm = stats[i].name.decode("utf-8") if stats[i].name else "unknown"
            result.append({
                "id":   stats[i].id,
                "name": nm,
                "sum":  stats[i].sum,
                "n":    stats[i].n,
            })
        return result

    finally:
        if annot_in  != NULL: free(annot_in)
        if labels_trg != NULL: free(labels_trg)
        if vals      != NULL: free(vals)
        if stats     != NULL: free(stats)
        if atable    != NULL: free(atable)
        if src_sph_objs != NULL:
            delete_object(src_sph_objs[0])
            free(src_sph_objs)
        if trg_sph_objs != NULL:
            delete_object(trg_sph_objs[0])
            free(trg_sph_objs)


# ===================================================================
# Surface ratio  (mirrors CAT_SurfRatio)
# ===================================================================
def surf_ratio(vertices, faces, double radius=20.0, bint normalize=True):
    """
    Compute per-vertex surface ratio based on the method of Toro et al. 2008.

    Mirrors ``CAT_SurfRatio``.  The ratio measures local gyrification by
    comparing the area of a geodesic disc around each vertex to the area
    of the corresponding region on the convex hull.

    Parameters
    ----------
    vertices : array_like, shape (V, 3)
    faces    : array_like, shape (F, 3)
    radius   : float
        Neighbourhood radius in mm (default 20.0).  When negative the
        radius is estimated automatically so that the global surface
        ratio is close to the global gyrification index.
    normalize : bool
        Normalize for individual surface area to make the measure
        scale-invariant (default True).

    Returns
    -------
    ratio_values : ndarray, shape (V,), float64
        Per-vertex surface ratio.
    """
    cdef PolygonsMesh mesh = _ensure_mesh(vertices, faces)
    cdef polygons_struct *poly = mesh.ptr()
    cdef int n = poly.n_points

    cdef double *ratio_ptr = C.get_surface_ratio(
        radius, poly, 1 if normalize else 0)
    if ratio_ptr == NULL:
        raise RuntimeError("get_surface_ratio returned NULL")

    cdef cnp.ndarray[cnp.float64_t, ndim=1] out = np.empty(n, dtype=np.float64)
    memcpy(<void *>out.data, ratio_ptr, sizeof(double) * <size_t>n)
    free(ratio_ptr)
    return out


# ===================================================================
# Fractal dimension  (mirrors CAT_SurfFractalDimension)
# ===================================================================
def surf_fractal_dimension(vertices, faces, sphere_vertices, sphere_faces,
                           int n_triangles=327680, bint use_sph=True,
                           int maxiters=30, bint smooth=True,
                           bint verbose=False):
    """
    Compute local and global fractal dimension of a cortical surface.

    Mirrors ``CAT_SurfFractalDimension``.  Two algorithms are
    available via ``use_sph``:

    * ``use_sph=True`` (default) — SPH bandwidth reduction method
      (Yotter et al.); produces high-quality local FD maps aligned
      to the reparameterized sphere.
    * ``use_sph=False`` — direct area-scaling resampling; faster but
      lower quality.

    Parameters
    ----------
    vertices, faces : array_like
        Central surface mesh (shape ``(V, 3)`` / ``(F, 3)``).
    sphere_vertices, sphere_faces : array_like
        Spherical parameterization of the same surface.
    n_triangles : int
        Number of triangles for the resampled sphere (SPH mode only,
        default 327680).
    use_sph : bool
        Use the spherical-harmonic method (default True).
    maxiters : int
        Number of scale iterations for the non-SPH method (default 30).
    smooth : bool
        Smooth local FD values after estimation (default True).
    verbose : bool
        Print intermediate progress (default False).

    Returns
    -------
    local_fd  : ndarray, float64
        Per-vertex fractal dimension values.  Length equals
        ``sphere.n_points`` (SPH mode) or ``surface.n_points``
        (non-SPH mode).
    global_fd : float
        Global fractal dimension of the surface.
    """
    import tempfile
    import os as _os

    cdef PolygonsMesh surface = _ensure_mesh(vertices, faces)
    cdef PolygonsMesh sphere  = _ensure_mesh(sphere_vertices, sphere_faces)

    # Reparameterization sphere: copy of sphere, centred + unit-radius
    cdef PolygonsMesh reparam = _ensure_mesh(sphere_vertices, sphere_faces)
    cdef polygons_struct *reparam_ptr = reparam.ptr()
    C.translate_to_center_of_mass(reparam_ptr)
    cdef int i
    for i in range(reparam_ptr.n_points):
        C.set_vector_length(&reparam_ptr.points[i], 1.0)

    # Write local FD values to a temporary file, read back as array
    cdef double global_fd = 0.0
    cdef int n_vals = 0
    cdef double *fd_ptr = NULL
    cdef cnp.ndarray[cnp.float64_t, ndim=1] local_fd

    fd_handle, tmp_path = tempfile.mkstemp(suffix=".txt")
    _os.close(fd_handle)
    cdef bytes b_tmp = tmp_path.encode("utf-8")

    try:
        if use_sph:
            global_fd = C.fractal_dimension_sph(
                surface.ptr(), sphere.ptr(), b_tmp,
                n_triangles, reparam_ptr,
                1 if smooth else 0,
                1 if verbose else 0)
        else:
            global_fd = C.fractal_dimension(
                surface.ptr(), sphere.ptr(), maxiters, b_tmp,
                1 if smooth else 0,
                1 if verbose else 0)

        if C.input_values_any_format(b_tmp, &n_vals, &fd_ptr) != OK:
            raise RuntimeError(
                "surf_fractal_dimension: could not read output values")

        local_fd = np.empty(n_vals, dtype=np.float64)
        memcpy(<void *>local_fd.data, fd_ptr,
               sizeof(double) * <size_t>n_vals)
        free(fd_ptr)
        fd_ptr = NULL
    finally:
        if fd_ptr != NULL:
            free(fd_ptr)
        try:
            _os.remove(tmp_path)
        except OSError:
            pass

    return local_fd, global_fd

