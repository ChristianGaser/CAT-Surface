# cython: language_level=3
"""
Cython declarations for the BIC volume_io / bicpl-surface types
used by libCAT.

These mirror the C typedefs and structs; they are **not** meant
to be imported from Python — use _convert.pyx for the safe
Python-level wrappers.
"""

# ---------------------------------------------------------------------------
# volume_io  basic types
# ---------------------------------------------------------------------------
cdef extern from "volume_io/basic.h":
    ctypedef int     VIO_BOOL
    ctypedef signed char VIO_SCHAR
    ctypedef double  VIO_Real
    ctypedef char*   VIO_STR

    ctypedef enum VIO_Status:
        OK
        ERROR
        INTERNAL_ERROR
        END_OF_FILE
        QUIT

    # Compat aliases (VIO_PREFIX_NAMES == 0)
    ctypedef VIO_BOOL      BOOLEAN
    ctypedef VIO_SCHAR     Smallest_int
    ctypedef VIO_Status    Status
    ctypedef VIO_Real      Real
    ctypedef VIO_STR       STRING


# ---------------------------------------------------------------------------
# volume_io  geometry structs  — Point, Vector, Colour, Surfprop
# ---------------------------------------------------------------------------
cdef extern from "volume_io/geom_structs.h":
    int VIO_N_DIMENSIONS

    ctypedef float VIO_Point_coord_type

    ctypedef struct VIO_Point:
        VIO_Point_coord_type coords[3]

    ctypedef struct VIO_Vector:
        VIO_Point_coord_type coords[3]

    ctypedef unsigned int VIO_Colour

    ctypedef float VIO_Spr_type

    ctypedef struct VIO_Surfprop:
        VIO_Spr_type a, d, s, se, t

    # Compat aliases
    ctypedef VIO_Point           Point
    ctypedef VIO_Vector          Vector
    ctypedef VIO_Colour          Colour
    ctypedef VIO_Surfprop        Surfprop
    ctypedef VIO_Point_coord_type Point_coord_type


# ---------------------------------------------------------------------------
# volume_io  files — File_formats
# ---------------------------------------------------------------------------
cdef extern from "volume_io/files.h":
    ctypedef enum VIO_File_formats:
        ASCII_FORMAT
        BINARY_FORMAT

    ctypedef VIO_File_formats File_formats


# ---------------------------------------------------------------------------
# bicpl  bintree
# ---------------------------------------------------------------------------
cdef extern from "bicpl/bintree.h":
    ctypedef struct bintree_struct:
        pass
    ctypedef bintree_struct* bintree_struct_ptr


# ---------------------------------------------------------------------------
# bicpl  obj_defs — polygons_struct, object_struct, etc.
# ---------------------------------------------------------------------------
cdef extern from "bicpl/obj_defs.h":
    ctypedef enum Colour_flags:
        ONE_COLOUR
        PER_ITEM_COLOURS
        PER_VERTEX_COLOURS

    ctypedef struct lines_struct:
        Colour_flags   colour_flag
        Colour        *colours
        float          line_thickness
        int            n_points
        Point         *points
        int            n_items
        int           *end_indices
        int           *indices
        bintree_struct_ptr bintree

    ctypedef struct polygons_struct:
        Colour_flags    colour_flag
        Colour         *colours
        Surfprop        surfprop
        float           line_thickness
        int             n_points
        Point          *points
        Vector         *normals
        int             n_items
        int            *end_indices
        int            *indices
        Smallest_int   *visibilities
        int            *neighbours
        bintree_struct_ptr bintree

    ctypedef enum Object_types:
        LINES
        MARKER
        MODEL
        PIXELS
        POLYGONS
        QUADMESH
        TEXT
        N_OBJECT_TYPES

    ctypedef struct object_struct:
        Object_types  object_type
        BOOLEAN       visibility
        # We access the specific union through helper functions,
        # so we don't need to declare the full union here.


# ---------------------------------------------------------------------------
# bicpl  object helpers (from object_prototypes.h)
# ---------------------------------------------------------------------------
cdef extern from "bicpl/object_prototypes.h":
    object_struct  *create_object(Object_types object_type)
    Object_types    get_object_type(object_struct *obj)
    polygons_struct *get_polygons_ptr(object_struct *obj)
    void            delete_object(object_struct *obj)
    void            initialize_polygons(polygons_struct *poly, Colour colour,
                                        Surfprop *surfprop)
    void            compute_polygon_normals(polygons_struct *poly)
