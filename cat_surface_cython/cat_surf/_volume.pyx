# cython: language_level=3, boundscheck=False, wraparound=False
"""
Shared NIfTI volume helper.

``open_volume()`` accepts a path, a ``(ndarray, affine)`` tuple, or any
nibabel-image-like object (anything with ``.affine`` and ``.get_fdata()``,
including ``nibabel.Nifti1Image``, ``surfa.Volume``, etc.) and returns a
``VolumeHandle`` that exposes the float32 data pointer and a minimal
in-memory ``nifti_image*`` suitable for the libCAT C functions.

For file paths the data is loaded via ``read_nifti_float`` (datatype
auto-converted, matches the CLI behaviour).  For in-memory inputs the
caller's ndarray is kept alive on the handle (no copy if it is already
float32 / Fortran-order); a fresh ``nifti_image`` is built from the
4x4 affine and the only fields libCAT actually consumes
(``nx, ny, nz``, ``dx, dy, dz``, ``sto_xyz``).

Always use ``with`` semantics or an explicit ``try: ... finally:
vh.close()`` — the handle owns C memory.
"""
import numpy as np
cimport numpy as cnp
from libc.stdlib cimport free
from libc.math cimport sqrt as c_sqrt

from cat_surf._nifti_types cimport (
    nifti_image, nifti_image_free, nifti_simple_init_nim,
    nifti_mat44_inverse, mat44, DT_FLOAT32,
)
cimport cat_surf._cat_funcs as C

cnp.import_array()


cdef class VolumeHandle:
    """Owns/borrows a NIfTI volume usable from C.

    Attributes are private; use :func:`open_volume` to construct.
    Always call :meth:`close` (or use ``try/finally``) when done — the
    handle holds malloc'd memory.
    """

    def __cinit__(self):
        self.nii = NULL
        self.data = NULL
        self._keepalive = None
        self._free_data = False
        self._free_nii = False
        self.dims[0] = 0
        self.dims[1] = 0
        self.dims[2] = 0

    cpdef void close(self) noexcept:
        if self._free_data and self.data != NULL:
            free(self.data)
        self.data = NULL
        if self.nii != NULL:
            # Prevent nifti_image_free from double-freeing borrowed numpy data.
            self.nii.data = NULL
            if self._free_nii:
                nifti_image_free(self.nii)
            self.nii = NULL
        self._keepalive = None

    def __dealloc__(self):
        self.close()


cdef nifti_image *_build_nii(cnp.ndarray arr_f32, cnp.ndarray affine) except NULL:
    """Build a minimal in-memory nifti_image from a 3-D float32 array + 4x4 affine."""
    if arr_f32.ndim != 3:
        raise ValueError("volume array must be 3-D")
    if affine.shape[0] != 4 or affine.shape[1] != 4:
        raise ValueError("affine must be 4x4")

    cdef nifti_image *nii = nifti_simple_init_nim()
    if nii == NULL:
        raise MemoryError("nifti_simple_init_nim failed")

    nii.nx = arr_f32.shape[0]
    nii.ny = arr_f32.shape[1]
    nii.nz = arr_f32.shape[2]
    nii.nt = 1; nii.nu = 1; nii.nv = 1; nii.nw = 1
    nii.ndim = 3
    nii.dim[0] = 3
    nii.dim[1] = nii.nx; nii.dim[2] = nii.ny; nii.dim[3] = nii.nz
    nii.dim[4] = 1; nii.dim[5] = 1; nii.dim[6] = 1; nii.dim[7] = 1
    nii.nvox = <size_t>nii.nx * <size_t>nii.ny * <size_t>nii.nz
    nii.nbyper = 4
    nii.datatype = DT_FLOAT32

    cdef double[:, :] a = affine
    cdef double dx = c_sqrt(a[0, 0] * a[0, 0]
                            + a[1, 0] * a[1, 0]
                            + a[2, 0] * a[2, 0])
    cdef double dy = c_sqrt(a[0, 1] * a[0, 1]
                            + a[1, 1] * a[1, 1]
                            + a[2, 1] * a[2, 1])
    cdef double dz = c_sqrt(a[0, 2] * a[0, 2]
                            + a[1, 2] * a[1, 2]
                            + a[2, 2] * a[2, 2])
    nii.dx = <float>dx
    nii.dy = <float>dy
    nii.dz = <float>dz
    nii.pixdim[0] = 1.0
    nii.pixdim[1] = nii.dx
    nii.pixdim[2] = nii.dy
    nii.pixdim[3] = nii.dz

    # sform = the supplied affine. libCAT reads sto_xyz only.
    nii.sform_code = 1  # NIFTI_XFORM_SCANNER_ANAT
    cdef int i, j
    for i in range(4):
        for j in range(4):
            nii.sto_xyz.m[i][j] = <float>a[i, j]
    nii.sto_ijk = nifti_mat44_inverse(nii.sto_xyz)
    nii.qform_code = 0

    nii.scl_slope = 0.0
    nii.scl_inter = 0.0
    nii.nifti_type = 1  # NIFTI-1 single-file

    # Borrow numpy buffer.
    nii.data = cnp.PyArray_DATA(arr_f32)
    return nii


cpdef VolumeHandle open_volume(object volume):
    """
    Coerce a volume input into a :class:`VolumeHandle`.

    Accepted forms:
    - ``str``: NIfTI file path; loaded with libCAT ``read_nifti_float``
      (auto-converts any datatype to float32).
    - ``(ndarray, affine)``: a 3-D array and a 4x4 affine.  The array is
      cast to float32 / Fortran order; no copy if already so.
    - nibabel-image-like: any object with ``.affine`` and
      ``.get_fdata()`` (nibabel ``Nifti1Image``, ``MGHImage``, surfa
      volumes, etc.).
    """
    cdef VolumeHandle vh = VolumeHandle.__new__(VolumeHandle)
    cdef bytes bpath
    cdef cnp.ndarray arr_f32, affine

    if isinstance(volume, str):
        bpath = volume.encode("utf-8")
        vh.nii = C.read_nifti_float(bpath, &vh.data, 0)
        if vh.nii == NULL or vh.data == NULL:
            raise IOError(f"Failed to read NIfTI '{volume}'")
        vh.dims[0] = vh.nii.nx
        vh.dims[1] = vh.nii.ny
        vh.dims[2] = vh.nii.nz
        vh._free_data = True
        vh._free_nii = True
        return vh

    if hasattr(volume, "affine") and hasattr(volume, "get_fdata"):
        data_py = volume.get_fdata()
        affine_py = volume.affine
    else:
        try:
            data_py, affine_py = volume
        except (TypeError, ValueError) as exc:
            raise TypeError(
                "volume must be a path, an (ndarray, affine) tuple, or a "
                "nibabel-image-like object (with .affine and .get_fdata())"
            ) from exc

    # 4-D series: take the mid-frame (matches read_nifti_float on disk).
    data_np = np.asarray(data_py)
    if data_np.ndim == 4:
        data_np = data_np[..., data_np.shape[3] // 2]

    arr_f32 = np.asfortranarray(data_np, dtype=np.float32)
    affine = np.ascontiguousarray(affine_py, dtype=np.float64)
    if affine.shape[0] != 4 or affine.shape[1] != 4:
        raise ValueError("affine must be a 4x4 matrix")

    vh.nii = _build_nii(arr_f32, affine)
    vh.data = <float *>cnp.PyArray_DATA(arr_f32)
    vh.dims[0] = arr_f32.shape[0]
    vh.dims[1] = arr_f32.shape[1]
    vh.dims[2] = arr_f32.shape[2]
    vh._keepalive = arr_f32   # keep ndarray alive while nii.data borrows it
    vh._free_data = False
    vh._free_nii = True
    return vh
