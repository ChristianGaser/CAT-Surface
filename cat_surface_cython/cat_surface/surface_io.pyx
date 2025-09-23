# cython: language_level=3
# distutils: language = c
# cython: boundscheck=False, wraparound=False, initializedcheck=False, cdivision=True

import os
cimport cython
import numpy as np
cimport numpy as cnp

from libc.stdlib cimport malloc, free
from libc.string cimport memcpy
from libc.stdint cimport uintptr_t

from .bicpl_defs cimport (
    Status, File_formats, ASCII_FORMAT, BINARY_FORMAT,
    object_struct, polygons_struct,
    input_graphics_file, get_object_type, get_polygons_ptr,
    delete_object_list,
    POLYGONS
)

ctypedef cnp.float32_t float32
ctypedef cnp.int32_t int32

@cython.cfunc
cdef inline void _copy_points_to_array(polygons_struct* poly, float32* out, Py_ssize_t stride):
    cdef int i
    cdef float x, y, z
    for i in range(poly.n_points):
        x = poly.points[i].coords[0]
        y = poly.points[i].coords[1]
        z = poly.points[i].coords[2]
        out[i*stride + 0] = <float32>x
        out[i*stride + 1] = <float32>y
        out[i*stride + 2] = <float32>z

@cython.cfunc
cdef inline void _copy_faces_to_array(polygons_struct* poly, int32* out, Py_ssize_t stride):
    cdef int p, e, start_idx, end_idx
    cdef int write_idx = 0
    for p in range(poly.n_items):
        end_idx = poly.end_indices[p]
        start_idx = 0 if p == 0 else poly.end_indices[p-1]
        if end_idx - start_idx != 3:
            # Skip non-tri polygons in this minimal example
            continue
        out[write_idx*stride + 0] = <int32>poly.indices[start_idx + 0]
        out[write_idx*stride + 1] = <int32>poly.indices[start_idx + 1]
        out[write_idx*stride + 2] = <int32>poly.indices[start_idx + 2]
        write_idx += 1

@cython.cfunc
cdef inline int _count_triangles(polygons_struct* poly):
    cdef int p, tri = 0
    cdef int start_idx, end_idx
    for p in range(poly.n_items):
        end_idx = poly.end_indices[p]
        start_idx = 0 if p == 0 else poly.end_indices[p-1]
        if end_idx - start_idx == 3:
            tri += 1
    return tri

@cython.boundscheck(False)
@cython.wraparound(False)
cpdef tuple load_bic_surface_arrays(str path):
    """
    Load a BICPL .obj surface and return (vertices, faces).

    - vertices: float32 array, shape (N, 3)
    - faces: int32 array, shape (M, 3) of triangles only (non-tris skipped)
    """
    cdef bytes b = os.fsencode(path)
    cdef char* cpath = b
    cdef File_formats fmt
    cdef int nobj = 0
    cdef object_struct** olist = NULL
    cdef Status st = input_graphics_file(cpath, &fmt, &nobj, &olist)
    if st != 0 or nobj <= 0 or olist == NULL:
        raise RuntimeError(f"input_graphics_file failed or empty: status={st}, nobj={nobj}")

    cdef polygons_struct* polys = NULL
    cdef int i
    for i in range(nobj):
        if get_object_type(olist[i]) == POLYGONS:
            polys = get_polygons_ptr(olist[i])
            break

    if polys == NULL:
        delete_object_list(nobj, olist)
        raise RuntimeError("No POLYGONS object found in file")

    # Create numpy arrays
    cdef Py_ssize_t n_points = polys.n_points
    cdef Py_ssize_t n_items = polys.n_items
    cdef int n_tri = _count_triangles(polys)

    cdef cnp.ndarray[float32, ndim=2] verts = np.empty((n_points, 3), dtype=np.float32)
    cdef cnp.ndarray[int32, ndim=2] faces = np.empty((n_tri, 3), dtype=np.int32)

    _copy_points_to_array(polys, <float32*>verts.data, 3)
    _copy_faces_to_array(polys, <int32*>faces.data, 3)

    delete_object_list(nobj, olist)
    return (verts, faces)

# Optional: small convenience wrapper class for a nibabel-like surface
cdef class Surface:
    cdef cnp.ndarray verts
    cdef cnp.ndarray faces
    cdef dict meta

    def __cinit__(self, cnp.ndarray verts, cnp.ndarray faces, meta=None):
        if verts.ndim != 2 or verts.shape[1] != 3:
            raise ValueError("verts must be (N,3)")
        if faces.ndim != 2 or faces.shape[1] != 3:
            raise ValueError("faces must be (M,3)")
        self.verts = verts.astype(np.float32, copy=False)
        self.faces = faces.astype(np.int32, copy=False)
        self.meta = {} if meta is None else dict(meta)

    @property
    def vertices(self):
        return self.verts

    @property
    def triangles(self):
        return self.faces

    def __repr__(self):
        return f"Surface(n_verts={self.verts.shape[0]}, n_faces={self.faces.shape[0]})"

cpdef Surface load_bic_surface(str path):
    v, f = load_bic_surface_arrays(path)
    return Surface(v, f)
