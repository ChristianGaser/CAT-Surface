# cython: language_level=3
"""
Shared NIfTI volume helper.  Exposes ``VolumeHandle`` and ``open_volume``
so other modules can accept volumes as a path, a (ndarray, affine)
tuple, or a nibabel-image-like object.
"""
cimport numpy as cnp
from cat_surf._nifti_types cimport nifti_image


# Minimal header for a 3-D float32 array and its 4x4 affine (sform); the
# array's buffer is borrowed as data.
cdef nifti_image *_build_nii(cnp.ndarray arr_f32, cnp.ndarray affine) except NULL


cdef class VolumeHandle:
    cdef nifti_image *nii
    cdef float       *data
    cdef int          dims[3]
    cdef object       _keepalive   # numpy array reference, or None
    cdef bint         _free_data
    cdef bint         _free_nii
    cpdef void close(self) noexcept


cpdef VolumeHandle open_volume(object volume)
