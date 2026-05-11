# cython: language_level=3, boundscheck=False, wraparound=False
"""
Python wrapper for volume-to-surface mapping — mirrors the
``CAT_Vol2Surf`` command-line tool.

Samples a 3-D volume along inward surface normals at a set of grid
positions, optionally relative to per-vertex cortical thickness, then
collapses the per-position samples with a mapping function (mean,
median, (abs-)max, sum, weighted average, exponential decay, multi).

Equi-volume sampling (Bok 1929, Waehnert et al. 2014) is supported and
requires per-vertex thickness.
"""

import numpy as np
cimport numpy as cnp
from libc.stdlib cimport malloc, free
from libc.math cimport sqrt as c_sqrt, isnan as c_isnan

from cat_surf._bic_types cimport (
    polygons_struct, Point, Vector,
    compute_polygon_normals,
)
from cat_surf._nifti_types cimport nifti_image, nifti_image_free
from cat_surf._convert cimport PolygonsMesh
from cat_surf._convert import arrays_to_polygons
from cat_surf._volume cimport VolumeHandle, open_volume

cimport cat_surf._cat_funcs as C

cnp.import_array()


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

# Mapping-function string → libCAT F_* enum value.  Populated lazily
# at module init to avoid Cython-time enum forwarding.
_MAP_FUNCS = None


def _map_funcs():
    global _MAP_FUNCS
    if _MAP_FUNCS is None:
        _MAP_FUNCS = {
            "mean":     C.F_MEAN,
            "avg":      C.F_MEAN,
            "median":   C.F_MEDIAN,
            "max":      C.F_MAX,
            "min":      C.F_MIN,
            "maxabs":   C.F_MAXABS,
            "sum":      C.F_SUM,
            "waverage": C.F_WAVERAGE,
            "exp":      C.F_EXP,
            "multi":    C.F_MULTI,
        }
    return _MAP_FUNCS


def vol2surf(volume,
             vertices,
             faces,
             thickness=None,
             offset=None,
             double offset_value=0.0,
             double grid_start=-0.5,
             double grid_end=0.5,
             int grid_steps=7,
             str map_func="maxabs",
             exp_half=None,
             frange=None,
             bint equivolume=False,
             bint verbose=False):
    """
    Map a volume to a surface by sampling along inward normals.

    Mirrors ``CAT_Vol2Surf``.

    Parameters
    ----------
    volume : str | (ndarray, affine) | nibabel-image-like
        Source volume.  Accepts a NIfTI file path, an ``(array, affine)``
        tuple, or any object with ``.affine`` and ``.get_fdata()``.
        Datatype is auto-converted to float32; 4-D arrays use the
        mid-frame.
    vertices, faces : mesh arrays
        Surface to sample on.  Normals are recomputed before sampling.
    thickness : array_like, shape (V,), float64, optional
        Per-vertex cortical thickness in mm.  When provided, ``grid_start``
        and ``grid_end`` are interpreted as relative positions in the GM
        band (e.g. -0.5 to +0.5 spans the full GM band).
    offset : array_like, shape (V,), float64, optional
        Per-vertex thickness used to offset the sampling origin.  Mutually
        exclusive with ``thickness``.  See ``offset_value``.
    offset_value : float
        Fraction of ``offset`` added to the sampling origin (default 0.0).
        E.g. ``offset_value=0.5`` with a central surface + thickness as
        offset samples the WM surface; ``-0.5`` samples the pial surface.
    grid_start : float
        Start of the sampling grid (default -0.5).  In mm unless
        ``thickness`` is set (then unitless relative to thickness).
    grid_end : float
        End of the sampling grid (default +0.5).
    grid_steps : int
        Number of grid steps (default 7).  Capped at 250.
    map_func : str
        Mapping function applied along normals:
        ``"mean"``/``"avg"``, ``"median"``, ``"max"``, ``"min"``,
        ``"maxabs"`` (default), ``"sum"``, ``"waverage"``
        (Gaussian-weighted average, kernel half-weight at boundaries),
        ``"exp"`` (exponential decay, set ``exp_half`` in mm),
        ``"multi"`` (return per-step values).
    exp_half : float, optional
        Half-decay distance in mm for ``map_func="exp"``.  Required when
        using ``"exp"``; ignored otherwise.
    frange : (float, float), optional
        Clip output values to this range (``-inf, +inf`` if omitted).
    equivolume : bool
        Use Bok (1929) / Waehnert et al. (2014) equi-volume sampling.
        Requires ``thickness``.
    verbose : bool

    Returns
    -------
    values : ndarray, shape (V,), float64
        Mapped scalar value per vertex.  For ``map_func="multi"``
        the shape is ``(V, grid_steps_used)`` where ``grid_steps_used``
        may exceed ``grid_steps`` for ``"waverage"`` (the kernel is
        extended to cover ~95 % of the Gaussian).
    grid_positions : ndarray, shape (grid_steps_used,), float64
        Sampling positions along the normal (mm if no thickness,
        unitless otherwise).
    """
    cdef dict map_table = _map_funcs()
    if map_func not in map_table:
        raise ValueError(
            f"map_func must be one of {sorted(map_table)}, got {map_func!r}")
    cdef int fc = map_table[map_func]

    if equivolume and thickness is None:
        raise ValueError("equivolume requires thickness to be provided")
    if thickness is not None and offset is not None:
        raise ValueError("provide either thickness or offset, not both")
    if map_func == "exp" and exp_half is None:
        raise ValueError("map_func='exp' requires exp_half (mm)")

    # --- Load volume ----------------------------------------------------
    cdef VolumeHandle vh = open_volume(volume)
    cdef int dims[3]
    dims[0] = vh.dims[0]; dims[1] = vh.dims[1]; dims[2] = vh.dims[2]

    # --- Build mesh -----------------------------------------------------
    cdef PolygonsMesh mesh = arrays_to_polygons(
        np.ascontiguousarray(vertices, dtype=np.float64),
        np.ascontiguousarray(faces, dtype=np.int32),
    )
    cdef polygons_struct *poly = mesh.ptr()
    cdef int n_pts = poly.n_points
    compute_polygon_normals(poly)

    # --- Validate thickness / offset arrays -----------------------------
    cdef cnp.ndarray[cnp.float64_t, ndim=1] thick_arr
    cdef double *thick_ptr = NULL
    if thickness is not None:
        thick_arr = np.ascontiguousarray(thickness, dtype=np.float64)
        if thick_arr.shape[0] != n_pts:
            vh.close()
            raise ValueError(
                f"thickness length ({thick_arr.shape[0]}) != "
                f"vertices ({n_pts})")
        thick_ptr = <double *>thick_arr.data
    elif offset is not None:
        thick_arr = np.ascontiguousarray(offset, dtype=np.float64)
        if thick_arr.shape[0] != n_pts:
            vh.close()
            raise ValueError(
                f"offset length ({thick_arr.shape[0]}) != "
                f"vertices ({n_pts})")
        thick_ptr = <double *>thick_arr.data

    # --- Extend grid for waverage to cover the Gaussian tails -----------
    cdef int grid_steps1 = grid_steps
    cdef double grid_start1 = grid_start
    cdef double grid_end1 = grid_end
    cdef double step_size
    cdef int grid_increase

    if fc == C.F_WAVERAGE:
        if grid_steps <= 1:
            raise ValueError("waverage needs grid_steps > 1")
        grid_increase = <int>round(2.0 * (grid_steps - 1.0) / 3.0)
        if grid_increase % 2:
            grid_increase += 1
        step_size = (grid_end - grid_start) / (grid_steps - 1.0)
        grid_steps1 = grid_steps + grid_increase
        grid_start1 = grid_start - grid_increase * step_size / 2.0
        grid_end1 = grid_end + grid_increase * step_size / 2.0

    if grid_steps1 > 250:
        vh.close()
        raise ValueError("effective grid_steps exceeds the 250-sample cap")

    # --- Length array ---------------------------------------------------
    cdef cnp.ndarray[cnp.float64_t, ndim=1] length_arr = np.empty(
        grid_steps1, dtype=np.float64)
    cdef int j
    for j in range(grid_steps1):
        if grid_steps1 > 1:
            length_arr[j] = grid_start1 + (
                <double>j / (grid_steps1 - 1.0)) * (grid_end1 - grid_start1)
        else:
            length_arr[j] = grid_start1

    # --- Build kernel ---------------------------------------------------
    cdef cnp.ndarray[cnp.float64_t, ndim=1] kernel_arr = np.zeros(
        grid_steps1, dtype=np.float64)
    cdef double *kernel_ptr = <double *>kernel_arr.data
    cdef double exp_half_val
    if map_func == "exp":
        exp_half_val = float(exp_half)
        C.CAT_Vol2SurfBuildExpKernel(<double *>length_arr.data,
                                     grid_steps1, exp_half_val, kernel_ptr)
    elif fc == C.F_WAVERAGE:
        C.CAT_Vol2SurfBuildGaussianKernel50(grid_steps, grid_steps1,
                                            kernel_ptr)

    cdef int fc_eval = C.F_EXP if map_func == "exp" else fc

    # --- Equivolume areas (optional) ------------------------------------
    cdef cnp.ndarray[cnp.float64_t, ndim=1] area_outer
    cdef cnp.ndarray[cnp.float64_t, ndim=1] area_inner
    if equivolume:
        area_outer = np.zeros(n_pts, dtype=np.float64)
        area_inner = np.zeros(n_pts, dtype=np.float64)
        C.get_area_of_points_central_to_pial(
            poly, <double *>area_outer.data, thick_ptr, 0.5)
        C.get_area_of_points_central_to_pial(
            poly, <double *>area_inner.data, thick_ptr, -0.5)

    # --- Output value clipping range ------------------------------------
    cdef double frange_lo = -np.inf
    cdef double frange_hi = np.inf
    if frange is not None:
        frange_lo = float(frange[0])
        frange_hi = float(frange[1])
        if frange_lo > frange_hi:
            vh.close()
            raise ValueError("frange[0] > frange[1]")

    # --- Sample ---------------------------------------------------------
    cdef cnp.ndarray[cnp.float64_t, ndim=1] values_1d
    cdef cnp.ndarray[cnp.float64_t, ndim=2] values_2d
    cdef bint is_multi = (fc == C.F_MULTI)

    if is_multi:
        values_2d = np.zeros((n_pts, grid_steps1), dtype=np.float64)
    else:
        values_1d = np.zeros(n_pts, dtype=np.float64)

    cdef double val_array[250]
    cdef int idx_dummy
    cdef int i
    cdef double pos
    cdef double v
    cdef double nx_, ny_, nz_
    cdef double px, py, pz
    cdef double diff_area, pos_rel
    cdef double v_clamped

    for i in range(n_pts):
        # Inward normal (negative outward).
        nx_ = -poly.normals[i].coords[0]
        ny_ = -poly.normals[i].coords[1]
        nz_ = -poly.normals[i].coords[2]
        px = poly.points[i].coords[0]
        py = poly.points[i].coords[1]
        pz = poly.points[i].coords[2]

        for j in range(grid_steps1):
            if thickness is None:
                if offset is not None:
                    pos = length_arr[j] + offset_value * thick_arr[i]
                else:
                    pos = length_arr[j]
            else:
                if equivolume:
                    diff_area = area_outer[i] - area_inner[i]
                    if diff_area != 0.0:
                        pos_rel = length_arr[j] + 0.5
                        pos_rel = (1.0 / diff_area) * (
                            c_sqrt(pos_rel * area_outer[i] * area_outer[i]
                                   + (1.0 - pos_rel) * area_inner[i]
                                   * area_inner[i])
                            - area_inner[i])
                        pos = (pos_rel - 0.5) * thick_arr[i]
                    else:
                        pos = length_arr[j] * thick_arr[i]
                else:
                    pos = length_arr[j] * thick_arr[i]

            v = C.isoval(vh.data,
                         <float>(px + pos * nx_),
                         <float>(py + pos * ny_),
                         <float>(pz + pos * nz_),
                         dims, vh.nii)
            if c_isnan(v):
                v = 0.0
            if is_multi:
                v_clamped = v
                if v_clamped < frange_lo:
                    v_clamped = frange_lo
                elif v_clamped > frange_hi:
                    v_clamped = frange_hi
                values_2d[i, j] = v_clamped
            else:
                val_array[j] = v

        if not is_multi:
            v = C.CAT_Vol2SurfEvaluateFunction(
                val_array, grid_steps1, fc_eval, kernel_ptr, &idx_dummy)
            if v < frange_lo:
                v = frange_lo
            elif v > frange_hi:
                v = frange_hi
            values_1d[i] = v

    vh.close()

    if is_multi:
        return values_2d, length_arr
    return values_1d, length_arr
