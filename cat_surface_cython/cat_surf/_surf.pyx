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
    create_object, get_polygons_ptr, delete_object,
    initialize_polygons, compute_polygon_normals, POLYGONS, OK,
)
from cat_surf._nifti_types cimport nifti_image, nifti_image_read, nifti_image_free
from cat_surf._convert cimport PolygonsMesh, _wrap_object
from cat_surf._convert import arrays_to_polygons, polygons_to_arrays

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
# Point distance
# ===================================================================
def point_distance(vertices1, faces1, vertices2, faces2, bint symmetric=True):
    """
    Compute per-vertex distances between two meshes.

    Parameters
    ----------
    vertices1, faces1 : mesh 1 arrays
    vertices2, faces2 : mesh 2 arrays
    symmetric : bool
        If True, compute symmetric (bidirectional) distance.

    Returns
    -------
    distances : ndarray, shape (V1,), float64
    max_distance : float
    """
    cdef PolygonsMesh m1 = _ensure_mesh(vertices1, faces1)
    cdef PolygonsMesh m2 = _ensure_mesh(vertices2, faces2)
    cdef polygons_struct *p1 = m1.ptr()
    cdef int n = p1.n_points

    cdef cnp.ndarray[cnp.float64_t, ndim=1] dist = np.zeros(n, dtype=np.float64)
    cdef double maxd = C.compute_point_distance(
        p1, m2.ptr(), <double *>dist.data, 1 if symmetric else 0,
    )
    return dist, maxd


def point_distance_mean(vertices1, faces1, vertices2, faces2,
                        bint symmetric=True):
    """
    Compute mean per-vertex distances between two meshes.

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
# Mesh reduction (quadric decimation)
# ===================================================================
def reduce_mesh(vertices, faces, int target_faces, double quality=1.0,
                int iterations=5, bint verbose=False):
    """
    Reduce the number of faces in a triangular mesh using
    quadric error metrics.

    Parameters
    ----------
    vertices : array_like, shape (V, 3)
    faces    : array_like, shape (F, 3)
    target_faces : int
        Desired number of triangles after decimation.
    quality  : float
        Quality parameter (0 = fast, 1 = best; default 1.0).
    iterations : int
        Number of decimation passes.
    verbose  : bool

    Returns
    -------
    new_vertices : ndarray, shape (V', 3), float64
    new_faces    : ndarray, shape (F', 3), int32
    """
    cdef PolygonsMesh mesh = _ensure_mesh(vertices, faces)
    cdef int rc = C.reduce_mesh_quadrics(
        mesh.ptr(), target_faces, quality, iterations,
        1 if verbose else 0,
    )
    if rc != 0:
        raise RuntimeError(f"reduce_mesh_quadrics returned error code {rc}")

    # Extract the modified mesh back to numpy
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
def surf_average(list surfaces):
    """
    Compute the vertex-wise average of multiple surfaces.

    All surfaces must share the same topology (same faces, same
    vertex count).

    Parameters
    ----------
    surfaces : list of (vertices, faces) tuples
        Each element is a pair of arrays ``(V_array, F_array)``.

    Returns
    -------
    avg_vertices : ndarray, shape (V, 3), float64
    faces        : ndarray, shape (F, 3), int32
    """
    if len(surfaces) == 0:
        raise ValueError("surfaces list must not be empty")

    cdef int n_sets = len(surfaces)
    v0, f0 = surfaces[0]
    avg = np.ascontiguousarray(v0, dtype=np.float64).copy()
    faces_out = np.ascontiguousarray(f0, dtype=np.int32)
    cdef int n_pts = avg.shape[0]

    cdef int i
    for i in range(1, n_sets):
        vi, fi = surfaces[i]
        vi = np.ascontiguousarray(vi, dtype=np.float64)
        if vi.shape[0] != n_pts:
            raise ValueError(
                f"Surface {i} has {vi.shape[0]} vertices, expected {n_pts}")
        avg += vi

    avg /= n_sets
    return avg, faces_out


# ===================================================================
# Correct thickness for folding  (mirrors CAT_SurfCorrectThicknessFolding)
# ===================================================================
def correct_thickness_folding(vertices, faces, thickness,
                              double slope=0.0):
    """
    Correct cortical thickness values for folding-related bias.

    Parameters
    ----------
    vertices : array_like, shape (V, 3)
    faces    : array_like, shape (F, 3)
    thickness : array_like, shape (V,), float64
        Uncorrected thickness values.
    slope : float
        Correction slope.  0 = automatic (calls the unweighted variant).

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
    cdef Status st
    if slope == 0.0:
        st = C.CAT_CorrectThicknessFolding(poly, n, <double *>out.data)
    else:
        st = C.CAT_CorrectThicknessFoldingWeighted(
            poly, n, <double *>out.data, slope)
    if st != OK:
        raise RuntimeError("CAT_CorrectThicknessFolding failed")
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
# Surface deformation  (mirrors CAT_SurfDeform)
# ===================================================================
def surf_deform(vertices, faces, str volume_file,
                double w1=0.1, double w2=0.1, double w3=0.1,
                double sigma=0.2, float isovalue=0.5,
                int iterations=100, bint remove_intersect=False,
                bint verbose=False):
    """
    Deform a surface toward a volume isovalue.

    The volume is read from a NIfTI file at the C level.

    Parameters
    ----------
    vertices, faces : mesh arrays
    volume_file : str
        Path to a NIfTI volume.
    w1, w2, w3 : float
        Smoothness, gradient, and balloon force weights.
    sigma : float
        Displacement smoothing sigma.
    isovalue : float
        Target isovalue in the volume.
    iterations : int
    remove_intersect : bool
    verbose : bool

    Returns
    -------
    new_vertices : ndarray, shape (V, 3), float64
    faces        : ndarray, shape (F, 3), int32
    """
    cdef bytes bfname = volume_file.encode("utf-8")
    cdef nifti_image *nii = nifti_image_read(bfname, 1)
    if nii == NULL:
        raise IOError(f"Failed to read NIfTI file '{volume_file}'")

    cdef PolygonsMesh mesh = _ensure_mesh(vertices, faces)

    cdef double w[3]
    w[0] = w1; w[1] = w2; w[2] = w3

    C.surf_deform(mesh.ptr(), <float *>nii.data, nii,
                  w, sigma, isovalue, iterations,
                  1 if remove_intersect else 0,
                  1 if verbose else 0)

    nifti_image_free(nii)
    return polygons_to_arrays(mesh)


# ===================================================================
# Pial / white surface estimation  (mirrors CAT_Surf2PialWhite)
# ===================================================================
def surf_to_pial_white(vertices, faces, thickness,
                       str label_file,
                       double w1=0.05, double w2=0.05, double w3=0.05,
                       double sigma=0.2, int iterations=100,
                       int gradient_iterations=30, int method=0,
                       bint verbose=False):
    """
    Estimate pial and white matter surfaces from a central surface.

    Parameters
    ----------
    vertices, faces : central surface mesh arrays
    thickness : array_like, shape (V,), float64
        Cortical thickness per vertex.
    label_file : str
        Path to the NIfTI label volume.
    w1, w2, w3 : float
        Force weights.
    sigma : float
    iterations, gradient_iterations : int
    method : int
        0 = deformation, 1 = Laplacian, 2 = ADE.
    verbose : bool

    Returns
    -------
    pial_vertices  : ndarray, shape (V, 3), float64
    pial_faces     : ndarray, shape (F, 3), int32
    white_vertices : ndarray, shape (V, 3), float64
    white_faces    : ndarray, shape (F, 3), int32
    """
    cdef bytes bfname = label_file.encode("utf-8")
    cdef nifti_image *nii = nifti_image_read(bfname, 1)
    if nii == NULL:
        raise IOError(f"Failed to read label file '{label_file}'")

    cdef PolygonsMesh mesh = _ensure_mesh(vertices, faces)
    cdef polygons_struct *poly = mesh.ptr()
    cdef int n = poly.n_points

    thick = np.ascontiguousarray(thickness, dtype=np.float64)
    if thick.shape[0] != n:
        nifti_image_free(nii)
        raise ValueError(
            f"thickness length ({thick.shape[0]}) != vertices ({n})")
    cdef cnp.ndarray[cnp.float64_t, ndim=1] thick_arr = thick

    # Allocate output surfaces
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
        poly, <double *>thick_arr.data, <float *>nii.data, nii,
        pial_poly, white_poly, &opts)

    nifti_image_free(nii)

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
