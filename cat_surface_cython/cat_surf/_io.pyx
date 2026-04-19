# cython: language_level=3, boundscheck=False, wraparound=False
"""
I/O wrappers — read/write surfaces and per-vertex data using
libCAT's multi-format routines (GIFTI, FreeSurfer, BIC obj, …).

These are useful when nibabel is not available or when you want to
use the exact same format-detection logic that the CAT CLI tools use.
"""

import numpy as np
cimport numpy as cnp
from libc.stdlib cimport malloc, free

from cat_surf._bic_types cimport (
    polygons_struct, object_struct, Object_types,
    File_formats, Status, POLYGONS, OK,
    get_polygons_ptr, get_object_type, delete_object,
)
from cat_surf._convert cimport PolygonsMesh, _wrap_object
from cat_surf._convert import polygons_to_arrays

cimport cat_surf._cat_funcs as C

cnp.import_array()


# ===================================================================
# Read surface
# ===================================================================
def read_surface(str filename):
    """
    Read a surface mesh from any format supported by libCAT.

    Supported formats: GIFTI (.gii), FreeSurfer (binary), BIC .obj,
    DFS, DX, OOGL.

    Parameters
    ----------
    filename : str
        Path to the surface file.

    Returns
    -------
    vertices : ndarray, shape (V, 3), float64
    faces    : ndarray, shape (F, 3), int32
    """
    cdef bytes bfname = filename.encode("utf-8")
    cdef char *cfname = bfname

    cdef File_formats fmt
    cdef int n_objects = 0
    cdef object_struct **objects = NULL
    cdef int i, j

    cdef Status st = C.input_graphics_any_format(
        cfname, &fmt, &n_objects, &objects,
    )
    if st != OK or n_objects < 1 or objects == NULL:
        raise IOError(f"Failed to read surface from '{filename}'")

    if get_object_type(objects[0]) != POLYGONS:
        # clean up
        for i in range(n_objects):
            delete_object(objects[i])
        free(objects)
        raise IOError(f"File '{filename}' does not contain a polygon mesh")

    # Wrap first object — we take ownership
    cdef PolygonsMesh mesh = _wrap_object(objects[0], owns=True)

    # Free remaining objects (if any) and the array pointer
    for j in range(1, n_objects):
        delete_object(objects[j])
    free(objects)

    return polygons_to_arrays(mesh)


# ===================================================================
# Write surface
# ===================================================================
def write_surface(str filename, vertices, faces):
    """
    Write a surface mesh to any format supported by libCAT.

    The format is determined by the file extension
    (.gii → GIFTI, no extension → FreeSurfer, .obj → BIC, …).

    Parameters
    ----------
    filename : str
    vertices : array_like, shape (V, 3)
    faces    : array_like, shape (F, 3)
    """
    from cat_surf._convert import arrays_to_polygons

    cdef bytes bfname = filename.encode("utf-8")
    cdef char *cfname = bfname

    cdef PolygonsMesh mesh = arrays_to_polygons(
        np.ascontiguousarray(vertices, dtype=np.float64),
        np.ascontiguousarray(faces, dtype=np.int32),
    )

    cdef object_struct *obj = mesh._obj
    cdef File_formats fmt = <File_formats>0   # BINARY_FORMAT — auto-detected by extension

    cdef Status st = C.output_graphics_any_format(
        cfname, fmt, 1, &obj, NULL,
    )
    if st != OK:
        raise IOError(f"Failed to write surface to '{filename}'")


# ===================================================================
# Read per-vertex values
# ===================================================================
def read_values(str filename):
    """
    Read per-vertex scalar data from any supported format.

    Supports: FreeSurfer curv, GIFTI func.gii, text, etc.

    Parameters
    ----------
    filename : str

    Returns
    -------
    values : ndarray, shape (V,), float64
    """
    cdef bytes bfname = filename.encode("utf-8")
    cdef char *cfname = bfname

    cdef int n_values = 0
    cdef double *vals = NULL

    cdef Status st = C.input_values_any_format(cfname, &n_values, &vals)
    if st != OK or vals == NULL:
        raise IOError(f"Failed to read values from '{filename}'")

    # Copy to numpy and free the C allocation
    cdef cnp.ndarray[cnp.float64_t, ndim=1] out = np.empty(
        n_values, dtype=np.float64)
    cdef int i
    for i in range(n_values):
        out[i] = vals[i]
    free(vals)

    return out


# ===================================================================
# Write per-vertex values
# ===================================================================
def write_values(str filename, values):
    """
    Write per-vertex scalar data.

    Format is inferred from the extension (.gii → GIFTI, etc.).

    Parameters
    ----------
    filename : str
    values   : array_like, shape (V,), float64
    """
    cdef bytes bfname = filename.encode("utf-8")
    cdef char *cfname = bfname

    cdef cnp.ndarray[cnp.float64_t, ndim=1] arr = np.ascontiguousarray(
        values, dtype=np.float64)
    cdef int n = arr.shape[0]

    # TYPE_DOUBLE = 1  (from CAT_SurfaceIO.h)
    cdef Status st = C.output_values_any_format(
        cfname, n, <void *>arr.data, 1,
    )
    if st != OK:
        raise IOError(f"Failed to write values to '{filename}'")
