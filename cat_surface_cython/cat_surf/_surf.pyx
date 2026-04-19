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

from cat_surf._bic_types cimport polygons_struct
from cat_surf._convert cimport PolygonsMesh
from cat_surf._convert import arrays_to_polygons

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
