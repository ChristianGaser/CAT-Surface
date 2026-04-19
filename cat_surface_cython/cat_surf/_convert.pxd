# cython: language_level=3
"""
Expose _convert internals to other .pyx modules via cimport.
"""
from cat_surf._bic_types cimport polygons_struct, object_struct

cdef class PolygonsMesh:
    cdef object_struct *_obj
    cdef polygons_struct *_poly
    cdef bint _owns
    cdef polygons_struct* ptr(self) noexcept

cdef PolygonsMesh _wrap_object(object_struct *obj, bint owns)
