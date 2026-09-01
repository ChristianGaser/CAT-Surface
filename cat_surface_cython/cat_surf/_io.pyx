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
def write_surface(str filename, vertices, faces, values=None):
    """
    Write a surface mesh to any format supported by libCAT.

    The format is determined by the file extension
    (.gii → GIFTI, no extension → FreeSurfer, .obj → BIC, …).

    Parameters
    ----------
    filename : str
    vertices : array_like, shape (V, 3)
    faces    : array_like, shape (F, 3)
    values   : array_like, shape (V,), optional
        Per-vertex scalar data to embed alongside the mesh.  Mirrors the
        ``values`` argument of ``output_graphics_any_format``: when
        supplied and the extension is ``.gii``/``.dat``, mesh and values
        are written into a single combined GIfTI (or ``.gii``+``.dat``
        pair) so the file can be re-loaded either as a surface (mesh
        only) or as values, matching the legacy CLI behavior.  Ignored
        for formats without a values channel.
    """
    from cat_surf._convert import arrays_to_polygons

    cdef bytes bfname = filename.encode("utf-8")
    cdef char *cfname = bfname

    cdef cnp.ndarray[cnp.float64_t, ndim=2] verts_arr = np.ascontiguousarray(
        vertices, dtype=np.float64)
    cdef PolygonsMesh mesh = arrays_to_polygons(
        verts_arr,
        np.ascontiguousarray(faces, dtype=np.int32),
    )

    cdef object_struct *obj = mesh._obj
    cdef File_formats fmt = <File_formats>0   # BINARY_FORMAT — auto-detected by extension

    cdef cnp.ndarray[cnp.float64_t, ndim=1] vals_arr
    cdef double *vals_ptr = NULL
    if values is not None:
        vals_arr = np.ascontiguousarray(values, dtype=np.float64)
        if vals_arr.shape[0] != verts_arr.shape[0]:
            raise ValueError(
                f"values length ({vals_arr.shape[0]}) does not match "
                f"vertex count ({verts_arr.shape[0]})")
        vals_ptr = <double *>vals_arr.data

    cdef Status st = C.output_graphics_any_format(
        cfname, fmt, 1, &obj, vals_ptr,
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


# ===================================================================
# GIFTI data arrays
# ===================================================================
def read_gifti_darrays(str filename):
    """
    List the DataArrays a GIFTI file contains.

    A ``.gii`` holds the mesh in two DataArrays (POINTSET and TRIANGLE)
    but may embed any number of further arrays next to it -- thickness,
    curvature, labels, time series.  :func:`read_surface` returns only
    the mesh, so this is the way to see what else the file carries.

    Parameters
    ----------
    filename : str
        Path to a GIFTI file.

    Returns
    -------
    list of dict
        One entry per DataArray, in file order, with the keys
        ``intent``, ``intent_name``, ``datatype``, ``encoding``,
        ``external_file``, ``name``, ``dims``, ``n_values``,
        ``n_nonfinite``, ``mean``, ``min`` and ``max``.

        ``mean``, ``min`` and ``max`` are ``None`` for the mesh arrays
        and for any array whose shape or datatype cannot be summarized.
        NaN and infinite entries are counted in ``n_nonfinite`` and left
        out of the statistics, so one masked-out vertex does not turn the
        whole summary into nan.
    """
    cdef bytes bfname = filename.encode("utf-8")
    cdef char *cfname = bfname

    cdef int n_arrays = 0
    cdef C.gifti_darray_info *arrays = NULL
    cdef int i, d

    cdef Status st = C.input_gifti_darrays(cfname, &n_arrays, &arrays)
    if st != OK:
        raise IOError(f"Failed to read GIFTI data arrays from '{filename}'")

    result = []
    try:
        for i in range(n_arrays):
            dims = [arrays[i].dims[d] for d in range(arrays[i].num_dim)]
            entry = {
                "intent": arrays[i].intent,
                "intent_name": arrays[i].intent_name.decode("utf-8"),
                "datatype": arrays[i].datatype_name.decode("utf-8"),
                "encoding": arrays[i].encoding_name.decode("utf-8"),
                "external_file": arrays[i].ext_fname.decode("utf-8") or None,
                "name": arrays[i].name.decode("utf-8") or None,
                "dims": dims,
                "n_values": arrays[i].n_values,
                "n_nonfinite": arrays[i].n_nonfinite,
            }
            if arrays[i].has_range:
                entry["mean"] = arrays[i].mean
                entry["min"] = arrays[i].min
                entry["max"] = arrays[i].max
            else:
                entry["mean"] = None
                entry["min"] = None
                entry["max"] = None
            result.append(entry)
    finally:
        free(arrays)

    return result
