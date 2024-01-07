#ifndef  DEF_geom_prototypes
#define  DEF_geom_prototypes

#ifdef __cplusplus
extern "C" {
#endif

BICAPI  double  fast_approx_sqrt(
    double  y );

BICAPI  Real  fast_approx_distance_between_points(
    Point  *p1,
    Point  *p2 );

BICAPI  int  clip_polygon_against_box(
    int     n_points,
    Point   points[],
    Real    x_min,
    Real    x_max,
    Real    y_min,
    Real    y_max,
    Real    z_min,
    Real    z_max,
    int     n_output_points,
    Point   output_points[] );

BICAPI  int  clip_polygon_against_plane(
    int     n_points,
    Point   points[],
    Real    plane_constant,
    Vector  *normal,
    Point   output_points[] );

BICAPI  void  split_polygon_with_plane(
    int     n_points,
    Point   points[],
    Real    plane_constant,
    Vector  *normal,
    int     *n_in,
    Point   in_points[],
    int     *n_out,
    Point   out_points[] );

BICAPI  void  get_closest_point_on_line_segment(
    Point  *point,
    Point  *p1,
    Point  *p2,
    Point  *closest_point );

BICAPI  Real  get_distance_to_line_segment(
    Point  *point,
    Point  *p1,
    Point  *p2,
    Real   *alpha );

BICAPI  Real  get_point_object_distance_sq(
    Point                 *point,
    object_struct         *object,
    int                   obj_index,
    Point                 *object_point );

BICAPI  Real  get_point_object_distance(
    Point                 *point,
    object_struct         *object,
    int                   obj_index,
    Point                 *object_point );

BICAPI  Real  get_point_object_vertex_distance(
    Point                 *point,
    object_struct         *object,
    int                   obj_index,
    int                   *object_vertex );

BICAPI  Real  find_closest_point_on_object(
    Point           *point,
    object_struct   *object,
    int             *obj_index,
    Point           *point_on_object );

BICAPI  Real  find_closest_vertex_on_object(
    Point           *point,
    object_struct   *object,
    int             *vertex_on_object );

BICAPI  void  get_polygon_vertex_curvatures(
    polygons_struct   *polygons,
    int               n_neighbours[],
    int               *neighbours[],
    Real              smoothing_distance,
    Real              low_threshold,
    Real              curvatures[] );

BICAPI  void  flatten_around_vertex(
    Point     *vertex,
    int       n_neighbours,
    Point     neighbours[],
    BOOLEAN   closed_flag,
    Real      x_flat[],
    Real      y_flat[] );

BICAPI  void  flatten_around_vertex_to_sphere(
    Real      radius,
    Point     *vertex,
    int       n_neighbours,
    Point     neighbours[],
    Real      x_sphere[],
    Real      y_sphere[],
    Real      z_sphere[] );

BICAPI  int  compute_distances_from_point(
    polygons_struct   *polygons,
    int               n_neighbours[],
    int               *neighbours[],
    Point             *point,
    int               poly,
    Real              max_distance,
    BOOLEAN           distances_initialized,
    float             distances[],
    int               *list[] ) ;

BICAPI  void  find_polygon_normal_no_normalize(
    int      n_points,
    Point    points[],
    Real     *nx,
    Real     *ny,
    Real     *nz );

BICAPI  void  find_polygon_normal(
    int      n_points,
    Point    points[],
    Vector   *normal );

BICAPI  void   get_plane_through_points(
    int      n_points,
    Point    points[],
    Vector   *normal,
    Real     *plane_constant );

BICAPI  Real  distance_from_plane(
    Point    *point,
    Vector   *plane_normal,
    Real     plane_constant );

BICAPI  Real  distance_from_line(
    Point    *point,
    Point    *end_point1,
    Point    *end_point2 );

BICAPI  BOOLEAN   line_segment_intersects_plane(
    Point   *p1,
    Point   *p2,
    Point   *plane_origin,
    Vector  *plane_normal,
    Point   *intersection_point );

BICAPI  BOOLEAN  get_nearest_point_on_lines(
    Point   *origin1,
    Vector  *delta1,
    Point   *origin2,
    Vector  *delta2,
    Point   *nearest_point );

BICAPI  BOOLEAN  clip_line_to_box(
    Point    *origin,
    Vector   *direction,
    Real     x_min,
    Real     x_max,
    Real     y_min,
    Real     y_max,
    Real     z_min,
    Real     z_max,
    Real     *t_min,
    Real     *t_max );

BICAPI  void  create_line_circle(
    Point            *centre,
    int              plane_axis,
    Real             x_radius,
    Real             y_radius,
    int              n_points,
    lines_struct     *lines );

BICAPI  void  get_polygon_interpolation_weights(
    Point       *point,
    int         n_points,
    Point       points[],
    Real        weights[] );

BICAPI  void  map_point_between_polygons(
    polygons_struct  *p1,
    int              poly_index,
    Point            *p1_point,
    polygons_struct  *p2,
    Point            *p2_point );

BICAPI  Real  map_point_to_unit_sphere(
    polygons_struct  *p,
    Point            *point,
    polygons_struct  *unit_sphere,
    Point            *unit_sphere_point );

BICAPI  void  map_unit_sphere_to_point(
    polygons_struct  *unit_sphere,
    Point            *unit_sphere_point,
    polygons_struct  *p,
    Point            *point );

BICAPI  void  polygon_transform_points(
    polygons_struct   *src_polygons,
    polygons_struct   *dest_polygons,
    int               n_points,
    Point             src_points[],
    Point             dest_points[] );

BICAPI  void  map_sphere_to_uv(
    Real    x,
    Real    y,
    Real    z,
    Real    *u,
    Real    *v );

BICAPI  void  map_uv_to_sphere(
    Real    u,
    Real    v,
    Real    *x,
    Real    *y,
    Real    *z );

BICAPI   void   find_path_between_polygons(
    int            polygon1,
    int            polygon2,
    int            n_polygons,
    int            end_indices[],
    Smallest_int   visibilities[],
    int            neighbours[],
    BOOLEAN        *path_exists,
    int            *path_length,
    int            *path[] );

BICAPI  void  create_unit_tetrahedron(
    polygons_struct  *polygons );

BICAPI  void  create_unit_cube(
    polygons_struct  *polygons );

BICAPI  void  create_unit_octohedron(
    polygons_struct  *polygons );

BICAPI  void  create_unit_icosahedron(
    polygons_struct  *polygons );

BICAPI  void   intersect_planes_with_polygons(
    polygons_struct   *polygons,
    Point             *plane_origin,
    Vector            *plane_normal,
    lines_struct      *lines );

BICAPI  void   intersect_planes_with_quadmesh(
    quadmesh_struct   *quadmesh,
    Point             *plane_origin,
    Vector            *plane_normal,
    lines_struct      *lines );

BICAPI  BOOLEAN  null_Point(
    Point   *p );

BICAPI  BOOLEAN  null_Vector(
    Vector   *v );

BICAPI  Real  distance_between_points(
    Point  *p1,
    Point  *p2 );

BICAPI  BOOLEAN  points_within_distance(
    Point  *p1,
    Point  *p2,
    Real   distance );

BICAPI  void  apply_point_to_min_and_max(
    Point   *point,
    Point   *min_point,
    Point   *max_point );

BICAPI  void  expand_min_and_max_points(
    Point   *min_point,
    Point   *max_point,
    Point   *min_to_check,
    Point   *max_to_check );

BICAPI  void  get_range_points(
    int                n_points,
    Point              points[],
    Point              *min_corner,
    Point              *max_corner );

BICAPI  void  get_points_centroid(
    int     n_points,
    Point   points[],
    Point   *centroid );

BICAPI   void     reverse_vectors(
    int       n_vectors,
    Vector    vectors[] );

BICAPI  Real  get_angle_between_points(
    Point  *prev_point,
    Point  *this_point,
    Point  *next_point );

BICAPI  Real  sq_distance_between_points(
    Point  *p1,
    Point  *p2 );

BICAPI  Real  find_point_polygon_distance_sq(
    Point     *point,
    int       n_points,
    Point     poly_points[],
    Point     *closest_point );

BICAPI  Real  find_point_polygon_distance(
    Point     *point,
    int       n_points,
    Point     poly_points[],
    Point     *closest_point );

BICAPI  int  find_closest_polygon_point(
    Point              *point,
    polygons_struct    *polygons,
    Point              *closest_point );

BICAPI  void  create_polygons_sphere(
    Point            *centre,
    Real             x_size,
    Real             y_size,
    Real             z_size,
    int              n_up,
    int              n_around,
    BOOLEAN          subdividing_flag,
    polygons_struct  *polygons );

BICAPI  int  get_sphere_point_index(
    int   up,
    int   around,
    int   n_up,
    int   n_around );

BICAPI  BOOLEAN  is_this_sphere_topology(
    polygons_struct  *polygons );

BICAPI  BOOLEAN  get_tessellation_of_polygons_sphere(
    polygons_struct  *polygons,
    int              *tess );

BICAPI  int  get_tessellation_with_n_points(
    int   n_points );

BICAPI  void  half_sample_sphere_tessellation(
    polygons_struct   *polygons,
    polygons_struct   *half );

BICAPI   void  initialize_intersect_directions( void );

BICAPI  Real  *get_intersect_directions( void );

BICAPI  BOOLEAN  point_within_polygon(
    Point   *pt,
    int     n_points,
    Point   points[],
    Vector  *polygon_normal );

BICAPI  BOOLEAN  line_intersects_ellipsoid(
    Point    *line_origin,
    Vector   *line_direction,
    Point    *sphere_centre,
    Real     x_size,
    Real     y_size,
    Real     z_size,
    Real     *t_min,
    Real     *t_max );

BICAPI  BOOLEAN  ray_intersects_sphere(
    Point       *origin,
    Vector      *direction,
    Point       *centre,
    Real        radius,
    Real        *dist );

BICAPI  void  intersect_ray_object(
    Point                 *origin,
    Vector                *direction,
    object_struct         *object,
    int                   obj_index,
    int                   *closest_obj_index,
    Real                  *closest_dist,
    int                   *n_intersections,
    Real                  *distances[] );

BICAPI  int  intersect_ray_with_object(
    Point           *origin,
    Vector          *direction,
    object_struct   *object,
    int             *obj_index,
    Real            *dist,
    Real            *distances[] );

BICAPI  void   remove_invisible_polygons(
    polygons_struct  *polygons,
    Smallest_int     visibilities[] );

BICAPI  Real  get_smooth_surface_curvature(
    polygons_struct   *polygons,
    int               n_neighbours[],
    int               *neighbours[],
    int               poly,
    int               vertex,
    BOOLEAN           distances_initialized,
    float             distances[],
    Real              smoothing_distance );

BICAPI  void  smooth_lines(
    lines_struct  *lines,
    Real          smooth_length );

BICAPI  void  create_line_spline(
    lines_struct  *lines,
    int           n_curve_segments,
    lines_struct  *new_lines );

BICAPI  void  smooth_polygon(
    polygons_struct  *polygons,
    Real             max_dist_from_original,
    Real             fraction_to_move,
    Real             stop_threshold,
    Real             normal_ratio,
    BOOLEAN          range_flag,
    volume_struct    *volume,
    int              min_value,
    int              max_value );

BICAPI  BOOLEAN  get_interpolation_weights_2d(
    Real   x,
    Real   y,
    int    n_points,
    Real   xs[],
    Real   ys[],
    Real   weights[] );

BICAPI  BOOLEAN  get_prediction_weights_2d(
    Real   x,
    Real   y,
    int    n_points,
    Real   xs[],
    Real   ys[],
    Real   *x_weights[2],
    Real   *x_constant,
    Real   *y_weights[2],
    Real   *y_constant );

BICAPI  BOOLEAN  get_prediction_weights_3d(
    Real   x,
    Real   y,
    Real   z,
    int    n_points,
    Real   xs[],
    Real   ys[],
    Real   zs[],
    Real   *x_weights[3],
    Real   *y_weights[3],
    Real   *z_weights[3] );

BICAPI  BOOLEAN is_single_closed_curve(
    lines_struct   *lines );

BICAPI  void  subdivide_lines(
    lines_struct  *lines );

BICAPI  void  subdivide_polygons(
    polygons_struct  *polygons );

BICAPI  Real  get_polygon_2d_area(
    int      n_points,
    Point    points[] );

BICAPI  Real  get_polygon_surface_area(
    int      n_points,
    Point    points[] );

BICAPI  Real  get_polygons_surface_area(
    polygons_struct  *polygons );

BICAPI  BOOLEAN  is_this_tetrahedral_topology(
    polygons_struct   *polygons );

BICAPI  int  get_tetra_tessellation_with_n_points(
    int   n_points );

BICAPI  void  create_tetrahedral_sphere(
    Point            *centre,
    Real             rx,
    Real             ry,
    Real             rz,
    int              n_triangles,
    polygons_struct  *polygons );

BICAPI  void  half_sample_tetrahedral_tessellation(
    polygons_struct  *polygons,
    polygons_struct  *half );

BICAPI  int  convert_lines_to_tubes(
    lines_struct     *lines,
    int              n_around,
    Real             radius,
    quadmesh_struct  *quadmeshes[] );

BICAPI  void   create_slice_quadmesh(
    Volume           volume,
    int              axis_index,
    Real             voxel_position,
    int              x_tess,
    int              y_tess,
    Real             x_min,
    Real             x_max,
    Real             y_min,
    Real             y_max,
    quadmesh_struct  *quadmesh );

BICAPI  void   create_slice_3d(
    Volume           volume,
    Point            *origin,
    Vector           *normal,
    polygons_struct  *polygons );

#ifdef __cplusplus
}
#endif

#endif
