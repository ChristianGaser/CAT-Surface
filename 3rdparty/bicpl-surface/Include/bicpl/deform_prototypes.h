#ifndef  DEF_deform_prototypes
#define  DEF_deform_prototypes

#ifdef __cplusplus
extern "C" {
#endif

BICAPI  void  deform_lines(
    lines_struct      *lines,
    deform_struct     *deform_parms );

BICAPI  void  deform_lines_one_iteration(
    lines_struct      *lines,
    deform_struct     *deform_parms,
    int               iteration );

BICAPI  void  get_line_equilibrium_point(
    lines_struct                 *lines,
    int                          axis,
    int                          point_index,
    int                          neighbours[],
    Real                         curvature_factors[],
    Real                         max_search_distance,
    int                          degrees_continuity,
    Volume                       volume,
    Volume                       label_volume,
    boundary_definition_struct   *boundary_def,
    deformation_model_struct     *deformation_model,
    Point                        *equilibrium_point,
    Point                        *boundary_point );

BICAPI  int  find_axial_plane(
    lines_struct   *lines );

BICAPI  void  deform_polygons(
    polygons_struct   *polygons,
    deform_struct     *deform_parms );

BICAPI  void  deform_polygons_one_iteration(
    polygons_struct   *polygons,
    deform_struct     *deform_parms,
    int               iteration );

BICAPI  BOOLEAN  find_boundary_in_direction(
    Volume                      volume,
    Volume                      label_volume,
    voxel_coef_struct           *lookup,
    bitlist_3d_struct           *done_bits,
    bitlist_3d_struct           *surface_bits,
    Real                        model_dist,
    Point                       *ray_origin,
    Vector                      *unit_pos_dir,
    Vector                      *unit_neg_dir,
    Real                        max_outwards_search_distance,
    Real                        max_inwards_search_distance,
    int                         degrees_continuity,
    boundary_definition_struct  *boundary_def,
    Real                        *boundary_distance );

BICAPI  int  find_voxel_line_polynomial(
    Real        coefs[],
    int         degrees_continuity,
    int         x,
    int         y,
    int         z,
    Real        line_origin[],
    Real        line_direction[],
    Real        line_poly[] );

BICAPI  int  find_voxel_line_value_intersection(
    Real        coefs[],
    int         degrees_continuity,
    int         x,
    int         y,
    int         z,
    Real        line_origin[],
    Real        line_direction[],
    Real        t_min,
    Real        t_max,
    Real        isovalue,
    Real        distances[3] );

BICAPI  void  initialize_deformation_model(
    deformation_model_struct  *model );

BICAPI  void  print_deformation_model(
    deformation_model_struct   *deformation_model );

BICAPI  Status   add_deformation_model(
    deformation_model_struct   *deformation_model,
    int                        up_to_n_points,
    Real                       model_weight,
    char                       model_filename[],
    Real                       min_curvature_offset,
    Real                       max_curvature_offset );

BICAPI  void  delete_deformation_model(
    deformation_model_struct  *model );

BICAPI  Status  input_original_positions(
    deformation_model_struct  *deform_model,
    char                      position_filename[],
    Real                      max_position_offset,
    int                       n_deforming_points );

BICAPI  BOOLEAN  check_correct_deformation_polygons(
    polygons_struct            *polygons,
    deformation_model_struct   *model );

BICAPI  BOOLEAN  check_correct_deformation_lines(
    lines_struct               *lines,
    deformation_model_struct   *model );

BICAPI  deform_model_struct  *find_relevent_model(
    deformation_model_struct  *model,
    int                       point_index );

BICAPI  void  get_model_shape_point(
    Point    *origin,
    Vector   *pos_model_dir,
    Vector   *neg_model_dir,
    Real     dist,
    Point    *point );

BICAPI  void  compute_equilibrium_point(
    int                       point_index,
    BOOLEAN                   boundary_exists,
    Real                      boundary_dist,
    Real                      base_length,
    Real                      model_dist,
    Vector                    *pos_model_dir,
    Vector                    *neg_model_dir,
    Point                     *centroid,
    deformation_model_struct  *deformation_model,
    Point                     *equilibrium_point );

BICAPI  void  compute_model_dirs(
    Point      *centroid,
    Vector     *normal,
    Real       base_length,
    Point      *model_point,
    Real       *model_dist,
    Point      *search_origin,
    Vector     *pos_model_dir,
    Vector     *neg_model_dir );

BICAPI  void  get_model_point(
    deformation_model_struct  *deformation_model,
    Point                     points[],
    int                       point_index,
    int                       n_neighbours,
    int                       neighbours[],
    Real                      curvatures[],
    Point                     *centroid,
    Vector                    *normal,
    Real                      base_length,
    Point                     *model_point );

BICAPI  void  get_neighbours_of_line_vertex(
    lines_struct    *lines,
    int             vertex_index,
    int             neighbours[2] );

BICAPI  BOOLEAN  deformation_model_includes_average(
    deformation_model_struct   *model );

BICAPI  Real  compute_line_curvature(
    lines_struct    *lines,
    int             axis,
    int             point_index,
    int             prev_point_index,
    int             next_point_index );

BICAPI  Real  deform_point(
    int                        point_index,
    Point                      points[],
    Point                      *equilibrium_point,
    Real                       fractional_step,
    Real                       max_step,
    BOOLEAN                    position_constrained,
    Real                       max_position_offset,
    Point                      original_positions[],
    Point                      *new_point );

BICAPI  void  compute_line_centroid_and_normal(
    lines_struct     *lines,
    int              axis,
    int              prev_point_index,
    int              next_point_index,
    Point            *centroid,
    Vector           *normal,
    Real             *base_length );

BICAPI  int  get_subsampled_neighbours_of_point(
    deformation_model_struct  *deformation_model,
    polygons_struct           *polygons,
    int                       poly,
    int                       vertex_index,
    int                       neighbours[],
    int                       max_neighbours,
    BOOLEAN                   *interior_flag );

BICAPI  BOOLEAN  is_point_inside_surface(
    Volume                      volume,
    Volume                      label_volume,
    int                         continuity,
    Real                        voxel[],
    Vector                      *direction,
    boundary_definition_struct  *boundary_def );

BICAPI  void   get_centre_of_cube(
    Point       *cube,
    int         sizes[3],
    Point       *centre );

BICAPI  BOOLEAN  contains_value(
    Real  values[2][2][2],
    int   sizes[3] );

BICAPI  BOOLEAN  cube_is_small_enough(
    Point     cube[2],
    int       sizes[3],
    Real      min_cube_size );

BICAPI  void  initialize_deform_stats(
    deform_stats  *stats );

BICAPI  void  record_error_in_deform_stats(
    deform_stats   *stats,
    Real           error );

BICAPI  void  print_deform_stats(
    deform_stats   *stats,
    int            n_points );

BICAPI  BOOLEAN   get_max_point_cube_distance(
    Point   cube[2],
    int     sizes[3],
    Point   *point,
    Real    *distance );

BICAPI  void  initialize_deformation_parameters(
    deform_struct  *deform );

BICAPI  void  delete_deformation_parameters(
    deform_struct  *deform );

BICAPI  void  set_boundary_definition(
    boundary_definition_struct  *boundary_def,
    Real                        min_value,
    Real                        max_value,
    Real                        grad_threshold,
    Real                        angle,
    char                        direction,
    Real                        tolerance );

BICAPI  void  initialize_lookup_volume_coeficients(
    voxel_coef_struct  *lookup );

BICAPI  void  lookup_volume_coeficients(
    voxel_coef_struct  *lookup,
    Volume             volume,
    int                degrees_continuity,
    int                x,
    int                y,
    int                z,
    Real               c[] );

BICAPI  void  delete_lookup_volume_coeficients(
    voxel_coef_struct  *lookup );


#ifdef __cplusplus
}
#endif

#endif
