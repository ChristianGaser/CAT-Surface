# cython: language_level=3, boundscheck=False, wraparound=False
"""
Conversion layer between numpy arrays and BIC polygons_struct.

Public helpers
--------------
arrays_to_polygons(vertices, faces)   → opaque PolygonsMesh wrapper
polygons_to_arrays(PolygonsMesh)      → (vertices, faces)

These are the *only* bridge between the Python world (numpy) and the
C world (polygons_struct).  All other .pyx modules work with
PolygonsMesh instances.
"""

import numpy as np
cimport numpy as cnp
from libc.stdlib cimport malloc, free
from libc.string cimport memcpy

from cat_surf._bic_types cimport (
    polygons_struct, object_struct, Object_types, Colour, Surfprop,
    Point, Vector, Smallest_int, bintree_struct_ptr,
    create_object, get_polygons_ptr, delete_object,
    initialize_polygons, compute_polygon_normals,
    ONE_COLOUR, POLYGONS,
)

cnp.import_array()

# Colour constant: WHITE = 0x00FFFFFF (RGBA with alpha=0)
DEF WHITE = 0x00FFFFFF


# ---------------------------------------------------------------------------
# PolygonsMesh — thin wrapper that owns a C polygons_struct
# ---------------------------------------------------------------------------
cdef class PolygonsMesh:
    """
    Opaque wrapper around a BIC ``polygons_struct *``.

    Owns the underlying C memory and frees it on dealloc.
    Not intended for direct Python construction — use
    :func:`arrays_to_polygons` or :func:`read_surface`.
    """

    def __cinit__(self):
        self._obj = NULL
        self._poly = NULL
        self._owns = False

    def __dealloc__(self):
        if self._owns and self._obj != NULL:
            delete_object(self._obj)
            self._obj = NULL
            self._poly = NULL

    @property
    def n_points(self):
        if self._poly == NULL:
            return 0
        return self._poly.n_points

    @property
    def n_items(self):
        """Number of triangular faces."""
        if self._poly == NULL:
            return 0
        return self._poly.n_items

    cdef polygons_struct* ptr(self) noexcept:
        """Return the raw C pointer (for use in other .pyx modules)."""
        return self._poly


# ---------------------------------------------------------------------------
# numpy → polygons_struct
# ---------------------------------------------------------------------------
def arrays_to_polygons(
    cnp.ndarray[cnp.float64_t, ndim=2] vertices not None,
    cnp.ndarray[cnp.int32_t, ndim=2] faces not None,
):
    """
    Build a :class:`PolygonsMesh` from numpy arrays.

    Parameters
    ----------
    vertices : ndarray, shape (V, 3), float64
        Vertex positions.
    faces : ndarray, shape (F, 3), int32
        Triangle indices (0-based).

    Returns
    -------
    PolygonsMesh
        Wrapper owning the C object.
    """
    cdef int n_points = vertices.shape[0]
    cdef int n_items  = faces.shape[0]

    if vertices.shape[1] != 3:
        raise ValueError("vertices must have shape (V, 3)")
    if faces.shape[1] != 3:
        raise ValueError("faces must have shape (F, 3)")

    # --- allocate ---
    cdef object_struct *obj = create_object(POLYGONS)
    if obj == NULL:
        raise MemoryError("create_object(POLYGONS) failed")

    cdef polygons_struct *poly = get_polygons_ptr(obj)
    initialize_polygons(poly, <Colour>WHITE, NULL)

    # --- points ---
    poly.n_points = n_points
    poly.points = <Point *>malloc(n_points * sizeof(Point))
    if poly.points == NULL:
        delete_object(obj)
        raise MemoryError("could not allocate points")

    cdef int i
    cdef double *vrow
    for i in range(n_points):
        vrow = &vertices[i, 0]
        poly.points[i].coords[0] = <float>vrow[0]
        poly.points[i].coords[1] = <float>vrow[1]
        poly.points[i].coords[2] = <float>vrow[2]

    # --- faces (triangles → flat indices + end_indices) ---
    poly.n_items = n_items

    cdef int n_indices = n_items * 3
    poly.end_indices = <int *>malloc(n_items * sizeof(int))
    poly.indices     = <int *>malloc(n_indices * sizeof(int))
    if poly.end_indices == NULL or poly.indices == NULL:
        delete_object(obj)
        raise MemoryError("could not allocate face indices")

    cdef int j, idx = 0
    for i in range(n_items):
        poly.indices[idx]     = <int>faces[i, 0]
        poly.indices[idx + 1] = <int>faces[i, 1]
        poly.indices[idx + 2] = <int>faces[i, 2]
        idx += 3
        poly.end_indices[i] = idx

    # --- normals (must be allocated before compute_polygon_normals) ---
    poly.normals = <Vector *>malloc(n_points * sizeof(Vector))
    if poly.normals == NULL:
        delete_object(obj)
        raise MemoryError("could not allocate normals")
    poly.neighbours = NULL
    poly.visibilities = NULL
    poly.bintree = NULL
    compute_polygon_normals(poly)

    # --- wrap ---
    cdef PolygonsMesh mesh = PolygonsMesh.__new__(PolygonsMesh)
    mesh._obj  = obj
    mesh._poly = poly
    mesh._owns = True
    return mesh


# ---------------------------------------------------------------------------
# polygons_struct → numpy
# ---------------------------------------------------------------------------
def polygons_to_arrays(PolygonsMesh mesh not None):
    """
    Extract numpy arrays from a :class:`PolygonsMesh`.

    Returns
    -------
    vertices : ndarray, shape (V, 3), float64
    faces    : ndarray, shape (F, 3), int32
    """
    cdef polygons_struct *poly = mesh.ptr()
    if poly == NULL:
        raise ValueError("empty PolygonsMesh")

    cdef int n_points = poly.n_points
    cdef int n_items  = poly.n_items

    # --- vertices ---
    cdef cnp.ndarray[cnp.float64_t, ndim=2] verts = np.empty(
        (n_points, 3), dtype=np.float64)
    cdef int i
    for i in range(n_points):
        verts[i, 0] = <double>poly.points[i].coords[0]
        verts[i, 1] = <double>poly.points[i].coords[1]
        verts[i, 2] = <double>poly.points[i].coords[2]

    # --- faces (flat indices → (F, 3)) ---
    cdef cnp.ndarray[cnp.int32_t, ndim=2] tris = np.empty(
        (n_items, 3), dtype=np.int32)
    cdef int start
    for i in range(n_items):
        start = 0 if i == 0 else poly.end_indices[i - 1]
        tris[i, 0] = poly.indices[start]
        tris[i, 1] = poly.indices[start + 1]
        tris[i, 2] = poly.indices[start + 2]

    return verts, tris


# ---------------------------------------------------------------------------
# Helper: build PolygonsMesh from an already-loaded object_struct*
# (used internally by _io.pyx — not exported to Python)
# ---------------------------------------------------------------------------
cdef PolygonsMesh _wrap_object(object_struct *obj, bint owns):
    """Wrap an existing object_struct* into PolygonsMesh."""
    cdef PolygonsMesh mesh = PolygonsMesh.__new__(PolygonsMesh)
    mesh._obj  = obj
    mesh._poly = get_polygons_ptr(obj)
    mesh._owns = owns
    return mesh
