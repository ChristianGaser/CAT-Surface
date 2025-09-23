# cython: language_level=3

from .bicpl_defs cimport (
    BOOLEAN, Real, Point, Vector, Colour, Surfprop, Status,
    File_formats, ASCII_FORMAT, BINARY_FORMAT,
    Colour_flags, ONE_COLOUR, PER_ITEM_COLOURS, PER_VERTEX_COLOURS,
    Object_types, LINES, MARKER, MODEL, PIXELS, POLYGONS, QUADMESH, TEXT,
    polygons_struct, object_struct,
    input_graphics_file, output_graphics_file, get_object_type,
    get_polygons_ptr, delete_object_list
)
