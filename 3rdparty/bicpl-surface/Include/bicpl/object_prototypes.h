#ifndef  DEF_object_prototypes
#define  DEF_object_prototypes

#ifndef  Real
#define  Real float
#endif

#ifdef __cplusplus
extern "C" {
#endif

BICAPI  void  coalesce_object_points(
    int      *n_points,
    Point    *points[],
    int      n_indices,
    int      indices[] );

BICAPI  void  separate_object_points(
    int           *new_n_points,
    Point         *points[],
    int           n_indices,
    int           indices[],
    Colour_flags  colour_flag,
    Colour        *colours[] );

BICAPI  void  get_default_surfprop(
    Surfprop  *surfprop );

BICAPI  BOOLEAN  lookup_colour(
    STRING    colour_name,
    Colour    *col );

BICAPI  BOOLEAN  lookup_colour_name(
    Colour  col,
    STRING  *colour_name );

BICAPI  void  rgb_to_hsl(
    Real    r,
    Real    g,
    Real    b,
    Real    *h,
    Real    *s,
    Real    *l );

BICAPI  void  hsl_to_rgb(
    Real   h,
    Real   sl,
    Real   l,
    Real   *r,
    Real   *g,
    Real   *b );

BICAPI  void  convert_colour_to_hsl(
    Colour   rgb,
    Colour   *hsl );

BICAPI  void  convert_colour_to_rgb(
    Colour   hsl,
    Colour   *rgb );

BICAPI  BOOLEAN  equal_colours(
    Colour  col1,
    Colour  col2 );

BICAPI  STRING  convert_colour_to_string(
    Colour   col );

BICAPI  Colour  convert_string_to_colour(
    STRING     string );

BICAPI  int  get_colour_distance(
    int      r,
    int      g,
    int      b,
    Colour   c2 );

BICAPI  int   find_closest_colour(
    int     r,
    int     g,
    int     b,
    int     n_colours,
    Colour  colours[] );

BICAPI  int  get_Colour_luminance(
    Colour   colour );

BICAPI  Colour  scale_colour(
    Colour   colour,
    Real     factor );

BICAPI  Status   input_graphics_file(
    STRING         filename,
    File_formats   *format,
    int            *n_objects,
    object_struct  ***object_list );

BICAPI  Status   output_graphics_file(
    STRING         filename,
    File_formats   format,
    int            n_objects,
    object_struct  *object_list[] );

BICAPI  Status   input_objects_any_format(
    STRING         filename,
    Colour         marker_colour,
    Real           marker_size,
    Marker_types   marker_type,
    int            *n_objects,
    object_struct  **object_list[] );

BICAPI  STRING  get_default_landmark_file_suffix( void );

BICAPI  Status   input_landmark_file(
    Volume         volume,
    STRING         filename,
    Colour         colour,
    Real           size,
    Marker_types   type,
    int            *n_objects,
    object_struct  **object_list[] );

BICAPI  Status  io_tag_point(
    FILE            *file,
    IO_types        io_direction,
    Volume          volume,
    Real            size,
    marker_struct   *marker );

BICAPI  void  initialize_lines(
    lines_struct    *lines,
    Colour          col );

BICAPI  void  initialize_lines_with_size(
    lines_struct    *lines,
    Colour          col,
    int             size,
    BOOLEAN         closed );

BICAPI  void  delete_lines( lines_struct *lines );

BICAPI  void  start_new_line( lines_struct *lines );

BICAPI  void  add_point_to_line(
    lines_struct   *lines,
    Point          *point );

BICAPI  void  get_line_segment_index(
    lines_struct  *lines,
    int           obj_index,
    int           *line,
    int           *seg );

BICAPI  void  set_bintree_delete_function(
    void  (*func)( bintree_struct_ptr * ) );

BICAPI  void  delete_bintree_if_any(
    bintree_struct_ptr   *bintree );

BICAPI  Real  get_lines_length(
    lines_struct  *lines );

BICAPI  void  get_lines_arc_point(
    lines_struct  *lines,
    Real          arc_length,
    Point         *point );

BICAPI  void   initialize_marker(
    marker_struct     *marker,
    Marker_types      type,
    Colour            colour );

BICAPI  void   delete_marker(
    marker_struct     *marker );

BICAPI  void  initialize_model(
    model_struct  *model );

BICAPI  void  delete_model(
    model_struct   *model );

BICAPI  void  assign_model_extra_ptr(
    model_struct   *model,
    void           *ptr );

BICAPI  void  *get_model_extra_ptr(
    model_struct   *model );

BICAPI  void  add_object_to_list(
    int            *n_objects,
    object_struct  ***object_list,
    object_struct  *new_object );

BICAPI  void  remove_object_from_list(
    int            *n_objects,
    object_struct  **object_list[],
    int            obj_index );

BICAPI  void  delete_object_list(
    int              n_objects,
    object_struct    *object_list[] );

BICAPI  int  find_object_index_in_model(
    model_struct   *model,
    object_struct  *object );

BICAPI  void  add_object_to_model(
    model_struct   *model,
    object_struct  *new_object );

BICAPI  void  insert_object_in_model(
    model_struct   *model,
    object_struct  *new_object,
    int            obj_index );

BICAPI  void  remove_ith_object_from_model(
    model_struct   *model,
    int            obj_index );

BICAPI  void  remove_object_from_model(
    model_struct   *model,
    object_struct  *object );

BICAPI  Status  io_lines(
    FILE                *file,
    IO_types            io_flag,
    File_formats        format,
    lines_struct        *lines );

BICAPI  Status  io_marker(
    FILE                *file,
    IO_types            io_flag,
    File_formats        format,
    marker_struct       *marker );

BICAPI  Status  io_model(
    FILE            *file,
    IO_types        io_flag,
    File_formats    format,
    model_struct    *model );

BICAPI  Status  io_pixels(
    FILE                *file,
    IO_types            io_flag,
    File_formats        format,
    pixels_struct       *pixels );

BICAPI  void  set_use_compressed_polygons_flag(
    BOOLEAN  value );

BICAPI  BOOLEAN  get_use_compressed_polygons_flag( void );

BICAPI  Status  io_polygons(
    FILE                *file,
    IO_types            io_flag,
    File_formats        format,
    polygons_struct     *polygons );

BICAPI  Status  io_quadmesh(
    FILE                *file,
    IO_types            io_flag,
    File_formats        format,
    quadmesh_struct     *quadmesh );

BICAPI  Status  io_text(
    FILE            *file,
    IO_types        io_flag,
    File_formats    format,
    text_struct     *text );

BICAPI  Status  io_point(
    FILE            *file,
    IO_types        io_flag,
    File_formats    format,
    Point           *point );

BICAPI  Status  io_vector(
    FILE            *file,
    IO_types        io_flag,
    File_formats    format,
    Vector          *v );

BICAPI  Status  io_colour(
    FILE            *file,
    IO_types        io_flag,
    File_formats    format,
    Colour          *colour );

BICAPI  Status  io_colours(
    FILE            *file,
    IO_types        io_flag,
    File_formats    format,
    Colour_flags    *colour_flag,
    int             n_items,
    int             n_points,
    Colour          **colours );

BICAPI  Status  io_surfprop(
    FILE            *file,
    IO_types        io_flag,
    File_formats    format,
    Surfprop        *surfprop );

BICAPI  Status  io_object_type(
    FILE           *file,
    IO_types       io_flag,
    File_formats   format,
    Object_types   type );

BICAPI  Status  input_object_type(
    FILE           *file,
    Object_types   *type,
    File_formats   *format,
    BOOLEAN        *eof );

BICAPI  Status  io_pixel_colours(
    FILE            *file,
    IO_types        io_flag,
    File_formats    format,
    int             n,
    Colour          *pixel_colours[] );

BICAPI  Status  io_pixel_colour(
    FILE            *file,
    IO_types        io_flag,
    File_formats    format,
    Colour          *pixel_colour );

BICAPI  Status  input_object(
    STRING         directory,
    FILE           *file,
    File_formats   *format,
    object_struct  **object,
    BOOLEAN        *eof );

BICAPI  Status  output_object(
    FILE           *file,
    File_formats   format,
    object_struct  *object );

BICAPI  object_struct  *create_object(
    Object_types   object_type );

BICAPI  Object_types  get_object_type(
    object_struct   *object );

BICAPI  BOOLEAN  get_object_visibility(
    object_struct  *object );

BICAPI  void  set_object_visibility(
    object_struct  *object,
    BOOLEAN        visibility );

BICAPI  lines_struct  *get_lines_ptr(
    object_struct  *object );

BICAPI  marker_struct  *get_marker_ptr(
    object_struct  *object );

BICAPI  model_struct  *get_model_ptr(
    object_struct  *object );

BICAPI  pixels_struct  *get_pixels_ptr(
    object_struct  *object );

BICAPI  polygons_struct  *get_polygons_ptr(
    object_struct  *object );

BICAPI  quadmesh_struct  *get_quadmesh_ptr(
    object_struct  *object );

BICAPI  text_struct  *get_text_ptr(
    object_struct  *object );

BICAPI  void  delete_object(
    object_struct  *object );

BICAPI  int  get_object_points(
    object_struct   *object,
    Point           *points[] );

BICAPI  int  get_object_normals(
    object_struct   *object,
    Vector          *normals[] );

BICAPI  Colour_flags  *get_object_colours(
    object_struct   *object,
    Colour          *colours[] );

BICAPI  void  set_object_colours(
    object_struct   *object,
    Colour          colours[] );

BICAPI  Surfprop  *get_object_surfprop(
    object_struct   *object );

BICAPI  STRING  get_object_name(
    object_struct  *object );

BICAPI  void  initialize_object_traverse(
    object_traverse_struct  *object_traverse,
    BOOLEAN                 visible_ones_only,
    int                     n_objects,
    object_struct           *object_list[] );

BICAPI  BOOLEAN  get_next_object_traverse(
    object_traverse_struct  *object_traverse,
    object_struct           **object );

BICAPI  void  terminate_object_traverse(
    object_traverse_struct   *object_traverse );

BICAPI  BOOLEAN  get_range_of_object(
    object_struct   *object,
    BOOLEAN         visible_ones_only,
    Point           *min_corner,
    Point           *max_corner );

BICAPI  void  reverse_object_normals(
    object_struct   *object );

BICAPI  BOOLEAN  get_object_colour(
    object_struct  *object,
    Colour         *colour );

BICAPI  void  set_object_colour(
    object_struct  *object,
    Colour         col );

BICAPI  void  set_object_surfprop(
    object_struct  *object,
    Surfprop       *spr );

BICAPI  void  free_colours(
    Colour_flags   colour_flag,
    Colour         colours[],
    int            n_points,
    int            n_items );

BICAPI  BOOLEAN  objects_are_same_topology(
    int    n_points1,
    int    n_items1,
    int    *end_indices1,
    int    *indices1,
    int    n_points2,
    int    n_items2,
    int    *end_indices2,
    int    *indices2 );

BICAPI  int  get_n_colours(
    Colour_flags  colour_flag,
    int           n_points,
    int           n_items );

BICAPI  void  initialize_pixels(
    pixels_struct  *pixels,
    int            x_position,
    int            y_position,
    int            x_size,
    int            y_size,
    Real           x_zoom,
    Real           y_zoom,
    Pixel_types    pixel_type );

BICAPI  void  delete_pixels( pixels_struct *pixels );

BICAPI  void  modify_pixels_size(
    int            *n_pixels_alloced,
    pixels_struct  *pixels,
    int            x_size,
    int            y_size,
    Pixel_types    pixel_type );

BICAPI  void  convert_pixels24_to_pixels8(
    pixels_struct   *pixels_rgb,
    pixels_struct   *pixels_8 );

BICAPI  void  convert_pixels24_to_index8(
    pixels_struct   *pixels_rgb,
    pixels_struct   *pixels_8,
    int             n_colours,
    Colour          colour_table[] );

BICAPI  void  convert_index8_to_pixels24(
    pixels_struct   *pixels_8,
    Colour          colour_table[],
    pixels_struct   *pixels_rgb );

BICAPI  void  convert_pixels24_to_gray_scale(
    pixels_struct   *pixels_rgb,
    pixels_struct   *pixels_8 );

BICAPI  void  convert_pixels24_to_dithered(
    pixels_struct   *pixels_rgb,
    pixels_struct   *pixels_8,
    int             n_colours,
    Colour          colour_table[] );

BICAPI  void  resample_pixels(
    pixels_struct   *pixels,
    Transform_2d    *transform,
    int             n_samples,
    pixels_struct   *new_pixels,
    Colour          background_colour );

BICAPI  void  copy_pixel_region(
    pixels_struct   *pixels,
    int             x_min,
    int             x_max,
    int             y_min,
    int             y_max,
    pixels_struct   *new_pixels );

BICAPI  void  check_polygons_neighbours_computed(
    polygons_struct   *polygons );

BICAPI  void  delete_polygon_point_neighbours(
    polygons_struct  *polygons,
    int              n_point_neighbours[],
    int              *point_neighbours[],
    Smallest_int     interior_flags[],
    int              *point_polygons[] );

BICAPI   void   create_polygon_point_neighbours(
    polygons_struct  *polygons,
    BOOLEAN          across_polygons_flag,
    int              *n_point_neighbours_ptr[],
    int              **point_neighbours_ptr[],
    Smallest_int     *interior_flags_ptr[],
    int              **point_polygons_ptr[] );

BICAPI  void  initialize_polygons(
    polygons_struct   *polygons,
    Colour            col,
    Surfprop          *spr );

BICAPI  void  initialize_polygons_with_size(
    polygons_struct   *polygons,
    Colour            col,
    Surfprop          *spr,
    int               n_points,
    int               n_polygons,
    int               size );

BICAPI  void  free_polygon_neighbours(
    polygons_struct   *polygons );

BICAPI  void  delete_polygons( polygons_struct *polygons );

BICAPI  void  copy_polygons(
    polygons_struct   *src,
    polygons_struct   *dest );

BICAPI  void  set_polygon_per_item_colours(
    polygons_struct   *polygons );

BICAPI  void  create_polygons_visibilities(
    polygons_struct   *polygons );

BICAPI  void  set_polygons_visibilities(
    polygons_struct   *polygons,
    BOOLEAN           state );

BICAPI  void  start_new_polygon(
    polygons_struct *polygons );

BICAPI  void  add_point_to_polygon(
    polygons_struct   *polygons,
    Point             *point,
    Vector            *normal );

BICAPI  int  get_polygon_points(
    polygons_struct   *polygons,
    int               poly,
    Point             points[] );

BICAPI  void  get_polygon_centroid(
    polygons_struct   *polygons,
    int               poly,
    Point             *centroid );

BICAPI  int  find_vertex_index(
    polygons_struct  *polygons,
    int              poly,
    int              point_index );

BICAPI  int  find_edge_index(
    polygons_struct  *polygons,
    int              poly,
    int              point_index1,
    int              point_index2 );

BICAPI  BOOLEAN  find_polygon_with_edge(
    polygons_struct  *polygons,
    int              point_index1,
    int              point_index2,
    int              *poly_containing_edge,
    int              *edge_index );

BICAPI  BOOLEAN  find_polygon_with_vertex(
    polygons_struct   *polygons,
    int               point_index,
    int               *poly_index,
    int               *vertex_index );

BICAPI  BOOLEAN  lookup_polygon_vertex(
    polygons_struct   *polygons,
    Point             *point,
    int               *point_index );

BICAPI  BOOLEAN  find_next_edge_around_point(
    polygons_struct   *polygons,
    int               poly,
    int               index_1,
    int               index_2,
    int               *next_poly,
    int               *next_index_1,
    int               *next_index_2 );

BICAPI  int  get_neighbours_of_point(
    polygons_struct   *polygons,
    int               poly,
    int               vertex_index,
    int               neighbours[],
    int               max_neighbours,
    BOOLEAN           *interior_point );

BICAPI  int  get_polygons_around_vertex(
    polygons_struct   *polygons,
    int               poly,
    int               vertex_index,
    int               poly_indices[],
    int               n_polys_alloced,
    BOOLEAN           *closed_flag );

BICAPI  void  compute_polygon_normal(
    polygons_struct  *polygons,
    int              poly,
    Vector           *normal );

BICAPI  void  compute_polygon_normals(
    polygons_struct  *polygons );

BICAPI  void  average_polygon_normals(
    polygons_struct  *polygons,
    int              n_iters,
    Real             neighbour_weight );

BICAPI  BOOLEAN  get_plane_polygon_intersection(
    Vector           *normal,
    Real             d,
    polygons_struct  *polygons,
    int              poly,
    Point            intersection_points[] );

BICAPI  BOOLEAN  get_plane_segment_intersection(
    Vector           *normal,
    Real             d,
    Point            *p1,
    Point            *p2,
    Point            *intersection_point );

BICAPI  void   reverse_polygons_vertices(
    polygons_struct  *polygons );

BICAPI  void   make_polygons_front_facing(
    polygons_struct  *polygons );

BICAPI  BOOLEAN  polygon_is_back_facing(
    polygons_struct   *polygons,
    int               poly );

BICAPI  void  compute_points_centroid_and_normal(
    polygons_struct  *polygons,
    int              point_index,
    int              n_neighbours,
    int              neighbours[],
    Point            *centroid,
    Vector           *normal,
    Real             *base_length,
    Real             *curvature );

BICAPI  void  compute_polygon_point_centroid(
    polygons_struct  *polygons,
    int              poly,
    int              vertex_index,
    int              point_index,
    Point            *centroid,
    Vector           *normal,
    Real             *base_length,
    Real             *curvature );

BICAPI  Real  compute_polygon_vertex_curvature(
    polygons_struct  *polygons,
    int              point_index );

BICAPI  Real  get_polygon_edge_angle(
    polygons_struct  *polygons,
    int              poly,
    int              edge );

BICAPI  BOOLEAN  polygons_are_same_topology(
    polygons_struct  *p1,
    polygons_struct  *p2 );

BICAPI  void  initialize_quadmesh(
    quadmesh_struct   *quadmesh,
    Colour            col,
    Surfprop          *spr,
    int               m,
    int               n );

BICAPI  void  delete_quadmesh(
    quadmesh_struct *quadmesh );

BICAPI  void  set_quadmesh_point(
    quadmesh_struct  *quadmesh,
    int              i,
    int              j,
    Point            *point,
    Vector           *normal );

BICAPI  BOOLEAN  get_quadmesh_point(
    const quadmesh_struct  *quadmesh,
    int              i,
    int              j,
    Point            *point );

BICAPI  void  get_quadmesh_n_objects(
    const quadmesh_struct  *quadmesh,
    int              *m,
    int              *n );

BICAPI  void  compute_quadmesh_normals(
    quadmesh_struct  *quadmesh );

BICAPI  void  get_quadmesh_patch_indices(
    const quadmesh_struct  *quadmesh,
    int              i,
    int              j,
    int              indices[] );

BICAPI  void  get_quadmesh_patch(
    const quadmesh_struct  *quadmesh,
    int              i,
    int              j,
    Point            points[] );

BICAPI  void  reverse_quadmesh_vertices(
    quadmesh_struct  *quadmesh );

BICAPI  Colour     *get_8bit_rgb_pixel_lookup( void );

BICAPI  int  convert_rgb_pixel_to_8bit_lookup(
    Colour    colour );

BICAPI  Status   input_tag_objects_file(
    STRING         filename,
    Colour         marker_colour,
    Real           default_size,
    Marker_types   default_type,
    int            *n_objects,
    object_struct  **object_list[] );

BICAPI  void   initialize_text(
    text_struct     *text,
    Point           *position,
    Colour          colour,
    Font_types      text_font,
    Real            size );

BICAPI  void   delete_text(
    text_struct     *text );

BICAPI  Status  output_texture_values(
    STRING         filename,
    File_formats   format,
    int            n_values,
    Real           values[] );

BICAPI  Status  input_texture_values(
    STRING         filename,
    int            *n_values,
    Real           *values[] );

#ifdef __cplusplus
}
#endif

#endif
