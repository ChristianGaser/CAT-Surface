# cython: language_level=3

# Low-level C typedefs for volume_io and bicpl needed to read BICPL objects

cdef extern from "volume_io.h":
    ctypedef int BOOLEAN
    ctypedef double Real

    cdef struct VIO_Point:
        float coords[3]
    ctypedef VIO_Point Point

    cdef struct VIO_Vector:
        float coords[3]
    ctypedef VIO_Vector Vector

    ctypedef unsigned int Colour

    cdef struct VIO_Surfprop:
        float a
        float d
        float s
        float se
        float t
    ctypedef VIO_Surfprop Surfprop

    ctypedef int Status

cdef extern from "volume_io/files.h":
    ctypedef enum File_formats:
        ASCII_FORMAT
        BINARY_FORMAT

cdef extern from "bicpl/obj_defs.h":
    cdef enum Colour_flags:
        ONE_COLOUR
        PER_ITEM_COLOURS
        PER_VERTEX_COLOURS

    cdef enum Object_types:
        LINES
        MARKER
        MODEL
        PIXELS
        POLYGONS
        QUADMESH
        TEXT
        N_OBJECT_TYPES

    ctypedef struct polygons_struct:
        Colour_flags colour_flag
        Colour* colours
        Surfprop surfprop
        float line_thickness
        int n_points
        Point* points
        Vector* normals
        int n_items
        int* end_indices
        int* indices
        unsigned char* visibilities
        int* neighbours
        void* bintree

    ctypedef struct object_struct:
        pass

# Ensure BICAPI macro is defined by including top-level bicpl.h
cdef extern from "bicpl.h":
    pass

cdef extern from "bicpl/object_prototypes.h":
    Status input_graphics_file(char* filename,
                               File_formats* format,
                               int* n_objects,
                               object_struct*** object_list)

    Status output_graphics_file(char* filename,
                                File_formats format,
                                int n_objects,
                                object_struct* object_list[])

    Object_types get_object_type(object_struct* object)

    polygons_struct* get_polygons_ptr(object_struct* object)

    void delete_object_list(int n_objects, object_struct* object_list[])
