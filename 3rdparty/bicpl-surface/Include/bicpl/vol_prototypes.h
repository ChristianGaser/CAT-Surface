#ifndef  DEF_vol_prototypes
#define  DEF_vol_prototypes

#ifdef __cplusplus
extern "C" {
#endif

BICAPI  Volume  create_box_filtered_volume(
    Volume   volume,
    nc_type  nc_data_type,
    BOOLEAN  sign_flag,
    Real     real_min_value,
    Real     real_max_value,
    Real     x_width,
    Real     y_width,
    Real     z_width );

BICAPI  void  modify_labels_in_range(
    Volume   volume,
    Volume   label_volume,
    int      src_min,
    int      src_max,
    int      dest_label,
    Real     min_threshold,
    Real     max_threshold,
    int      range_changed[2][N_DIMENSIONS] );

BICAPI  void  initialize_colour_coding(
    colour_coding_struct   *colour_coding,
    Colour_coding_types    type,
    Colour                 under_colour,
    Colour                 over_colour,
    Real                   min_value,
    Real                   max_value );

BICAPI  void  delete_colour_coding(
    colour_coding_struct   *colour_coding );

BICAPI  void  set_colour_coding_type(
    colour_coding_struct   *colour_coding,
    Colour_coding_types    type );

BICAPI  Colour_coding_types  get_colour_coding_type(
    colour_coding_struct   *colour_coding );

BICAPI  void  set_colour_coding_min_max(
    colour_coding_struct   *colour_coding,
    Real                   min_value,
    Real                   max_value );

BICAPI  void  get_colour_coding_min_max(
    colour_coding_struct   *colour_coding,
    Real                   *min_value,
    Real                   *max_value );

BICAPI  Colour  get_colour_coding_under_colour(
    colour_coding_struct   *colour_coding );

BICAPI  void  set_colour_coding_under_colour(
    colour_coding_struct   *colour_coding,
    Colour                 under_colour );

BICAPI  Colour  get_colour_coding_over_colour(
    colour_coding_struct   *colour_coding );

BICAPI  void  set_colour_coding_over_colour(
    colour_coding_struct   *colour_coding,
    Colour                 over_colour );

BICAPI  BOOLEAN  define_colour_coding_user_defined(
    colour_coding_struct  *colour_code,
    int                   n_colours,
    Colour                colours[],
    Real                  positions[],
    Colour_spaces         interpolation_space );

BICAPI  Colour  get_colour_code(
    colour_coding_struct  *colour_coding,
    Real                  value );

BICAPI  void  colour_code_object(
    Volume                 volume,
    int                    continuity,
    colour_coding_struct   *colour_coding,
    object_struct          *object );

BICAPI  STRING    get_default_user_def_colour_code_suffix( void );

BICAPI  Status  input_user_defined_colour_coding(
    colour_coding_struct   *colour_coding,
    STRING                 filename );

BICAPI  void  create_volume_slice(
    Volume          volume1,
    Filter_types    filter_type1,
    Real            filter_width1,
    Real            slice_position1[],
    Real            x_axis1[],
    Real            y_axis1[],
    Real            x_translation1,
    Real            y_translation1,
    Real            x_scale1,
    Real            y_scale1,
    Volume          volume2,
    Filter_types    filter_type2,
    Real            filter_width2,
    Real            slice_position2[],
    Real            x_axis2[],
    Real            y_axis2[],
    Real            x_translation2,
    Real            y_translation2,
    Real            x_scale2,
    Real            y_scale2,
    int             x_viewport_size,
    int             y_viewport_size,
    int             x_pixel_start,
    int             x_pixel_end,
    int             y_pixel_start,
    int             y_pixel_end,
    Pixel_types     pixel_type,
    int             degrees_continuity,
    unsigned short  **cmode_colour_map,
    Colour          **rgb_colour_map,
    Colour          empty_colour,
    void            *render_storage,
    BOOLEAN         clip_pixels_flag,
    int             *n_pixels_alloced,
    pixels_struct   *pixels );

BICAPI  void  set_volume_slice_pixel_range(
    Volume          volume1,
    Filter_types    filter_type1,
    Real            filter_width1,
    Real            slice_position1[],
    Real            x_axis1[],
    Real            y_axis1[],
    Real            x_translation1,
    Real            y_translation1,
    Real            x_scale1,
    Real            y_scale1,
    Volume          volume2,
    Filter_types    filter_type2,
    Real            filter_width2,
    Real            slice_position2[],
    Real            x_axis2[],
    Real            y_axis2[],
    Real            x_translation2,
    Real            y_translation2,
    Real            x_scale2,
    Real            y_scale2,
    int             x_viewport_size,
    int             y_viewport_size,
    Pixel_types     pixel_type,
    int             *n_pixels_alloced,
    pixels_struct   *pixels );

BICAPI  BOOLEAN  find_volume_crop_bounds(
    Volume          volume,
    Real            min_crop_threshold,
    Real            max_crop_threshold,
    int             limits[2][MAX_DIMENSIONS] );

BICAPI  Volume  create_cropped_volume(
    Volume          volume,
    int             limits[2][MAX_DIMENSIONS] );

BICAPI  Volume  autocrop_volume(
    Volume    volume );

BICAPI  int  dilate_voxels_3d(
    Volume          volume,
    Volume          label_volume,
    Real            min_inside_label,
    Real            max_inside_label,
    Real            min_inside_value,
    Real            max_inside_value,
    Real            min_outside_label,
    Real            max_outside_label,
    Real            min_outside_value,
    Real            max_outside_value,
    Real            new_label,
    Neighbour_types connectivity,
    int             range_changed[2][N_DIMENSIONS] );

BICAPI  int  get_slice_weights_for_filter(
    Volume         volume,
    Real           voxel_position[],
    Real           voxel_direction[],   /* if filter_type != NEAREST */
    Filter_types   filter_type,
    Real           full_width_half_max,
    Real           ***positions,
    Real           *weights[] );

BICAPI  BOOLEAN  fill_connected_voxels(
    Volume              volume,
    Volume              label_volume,
    Neighbour_types     connectivity,
    int                 voxel[],
    int                 min_label_threshold,
    int                 max_label_threshold,
    int                 desired_label,
    Real                min_threshold,
    Real                max_threshold,
    int                 range_changed[2][N_DIMENSIONS] );

BICAPI  int  get_neighbour_directions(
    Neighbour_types   connectivity,
    int               *dx[],
    int               *dy[] );

BICAPI  int  get_3D_neighbour_directions(
    Neighbour_types   connectivity,
    int               *dx[],
    int               *dy[],
    int               *dz[] );

BICAPI  BOOLEAN  should_change_this_one(
    Volume          volume,
    Volume          label_volume,
    int             voxel[],
    Real            min_threshold,
    Real            max_threshold,
    int             label_min_threshold,
    int             label_max_threshold,
    int             desired_label );

BICAPI  void  interpolate_volume_to_slice(
    Volume          volume1,
    int             n_dims1,
    Real            origin1[],
    Real            x_axis1[],
    Real            y_axis1[],
    Volume          volume2,
    int             n_dims2,
    Real            origin2[],
    Real            x_axis2[],
    Real            y_axis2[],
    int             x_pixel_start,
    int             x_pixel_end,
    int             y_pixel_start,
    int             y_pixel_end,
    int             degrees_continuity,
    unsigned short  **cmode_colour_map,
    Colour          **rgb_colour_map,
    Colour          empty_colour,
    pixels_struct   *pixels );

BICAPI  Status  input_volume_header_only(
    STRING               filename,
    int                  n_dimensions,
    STRING               dim_names[],
    Volume               *volume,
    minc_input_options   *options );

BICAPI  void  set_label_volume_real_range(
    Volume  volume );

BICAPI  Volume  create_label_volume(
    Volume  volume,
    nc_type type );

BICAPI  BOOLEAN  is_label_volume_initialized(
    Volume  volume );

BICAPI  void  set_all_volume_label_data(
    Volume    volume,
    int       value );

BICAPI  void  set_volume_label_data_5d(
    Volume          volume,
    int             v0,
    int             v1,
    int             v2,
    int             v3,
    int             v4,
    int             value );

BICAPI  void  set_volume_label_data(
    Volume          volume,
    int             voxel[],
    int             value );

BICAPI  int  get_volume_label_data(
    Volume          volume,
    int             voxel[] );

BICAPI  int  get_volume_label_data_5d(
    Volume          volume,
    int             v0,
    int             v1,
    int             v2,
    int             v3,
    int             v4 );

BICAPI  int  get_3D_volume_label_data(
    Volume          volume,
    int             x,
    int             y,
    int             z );

BICAPI  BOOLEAN  get_voxel_label_bit(
    Volume          volume,
    int             voxel[],
    int             bit );

BICAPI  void  set_voxel_label_bit(
    Volume          volume,
    int             voxel[],
    int             bit,
    BOOLEAN         value );

BICAPI  void  set_all_volume_label_data_bit(
    Volume         volume,
    int            bit,
    BOOLEAN        value );

BICAPI  BOOLEAN  get_volume_voxel_activity(
    Volume     volume,
    Real       voxel[],
    BOOLEAN    activity_if_mixed );

BICAPI  Status  load_label_volume(
    STRING   filename,
    Volume   label_volume );

BICAPI  Status  save_label_volume(
    STRING   filename,
    STRING   original_filename,
    Volume   label_volume,
    Real     crop_threshold );

BICAPI  Status  input_tags_as_labels(
    FILE    *file,
    Volume  volume,
    Volume  label_volume );

BICAPI  Status  create_label_volume_from_file(
    STRING   filename,
    Volume   volume,
    Volume   *label_volume );

BICAPI  Status  output_labels_as_tags(
    FILE    *file,
    Volume  volume,
    Volume  label_volume,
    int     desired_label,
    Real    size,
    int     patient_id );

BICAPI  Status  input_landmarks_as_labels(
    FILE    *file,
    Volume  volume,
    Volume  label_volume );

BICAPI  int    get_volume_cross_section(
    Volume   volume,
    Real     origin[],
    Real     x_axis[],
    Real     y_axis[],
    Real     clipped_voxels[][MAX_DIMENSIONS] );

BICAPI  void    get_volume_mapping_range(
    Volume   volume,
    Real     origin[],
    Real     x_axis[],
    Real     y_axis[],
    Real     x_trans,
    Real     y_trans,
    Real     x_scale,
    Real     y_scale,
    Real     *x_pixel_start,
    Real     *x_pixel_end,
    Real     *y_pixel_start,
    Real     *y_pixel_end );

BICAPI  void    clip_viewport_to_volume(
    Volume   volume,
    Real     origin[],
    Real     x_axis[],
    Real     y_axis[],
    int      *x_pixel_start,
    int      *x_pixel_end,
    int      *y_pixel_start,
    int      *y_pixel_end );

BICAPI  void  get_mapping(
    Volume          volume,
    Real            origin[],
    Real            x_axis[],
    Real            y_axis[],
    Real            x_translation,
    Real            y_translation,
    Real            x_scale,
    Real            y_scale,
    Real            pix_origin[],
    Real            pix_x_axis[],
    Real            pix_y_axis[] );

BICAPI  void  map_voxel_to_pixel(
    int    n,
    Real   voxel[],
    Real   origin[],
    Real   x_axis[],
    Real   y_axis[],
    Real   *x_pixel,
    Real   *y_pixel );

BICAPI  BOOLEAN  convert_slice_pixel_to_voxel(
    Volume          volume,
    Real            x_pixel,
    Real            y_pixel,
    Real            origin[],
    Real            x_axis[],
    Real            y_axis[],
    Real            x_translation,
    Real            y_translation,
    Real            x_scale,
    Real            y_scale,
    Real            voxel[] );

BICAPI  void  convert_voxel_to_slice_pixel(
    Volume          volume,
    Real            voxel[],
    Real            origin[],
    Real            x_axis[],
    Real            y_axis[],
    Real            x_translation,
    Real            y_translation,
    Real            x_scale,
    Real            y_scale,
    Real            *x_pixel,
    Real            *y_pixel );

BICAPI  void  resize_volume_slice(
    int          old_x_viewport_size,
    int          old_y_viewport_size,
    int          old_used_x_viewport_size,
    int          old_used_y_viewport_size,
    int          new_x_viewport_size,
    int          new_y_viewport_size,
    Real         *x_translation,
    Real         *y_translation,
    Real         *x_scale,
    Real         *y_scale,
    int          *used_x_viewport_size,
    int          *used_y_viewport_size );

BICAPI  void  fit_volume_slice_to_viewport(
    Volume       volume,
    Real         origin[],
    Real         x_axis[],
    Real         y_axis[],
    int          x_viewport_size,
    int          y_viewport_size,
    Real         fraction_oversize,
    Real         *x_translation,
    Real         *y_translation,
    Real         *x_scale,
    Real         *y_scale,
    int          *used_x_viewport_size,
    int          *used_y_viewport_size );

BICAPI  void   scale_slice_about_viewport_centre(
    Real        scale_factor,
    int         x_viewport_size,
    int         y_viewport_size,
    Real        *x_translation,
    Real        *y_translation,
    Real        *x_scale,
    Real        *y_scale );

BICAPI  BOOLEAN  voxel_is_within_volume(
    Volume   volume,
    Real     voxel_position[] );

BICAPI  BOOLEAN  int_voxel_is_within_volume(
    Volume   volume,
    int      indices[] );

BICAPI  void  convert_real_to_int_voxel(
    int      n_dimensions,
    Real     voxel[],
    int      int_voxel[] );

BICAPI  void  convert_int_to_real_voxel(
    int      n_dimensions,
    int      int_voxel[],
    Real     voxel[] );

BICAPI  BOOLEAN  voxel_contains_range(
    Volume   volume,
    int      int_voxel[],
    Real     min_value,
    Real     max_value );

BICAPI  BOOLEAN  volumes_are_same_grid(
    Volume    volume1,
    Volume    volume2 );

BICAPI  Status  output_volume_free_format(
    STRING         prefix,
    Volume         volume,
    int            axis_ordering[] );

BICAPI  void   *initialize_render_storage( void );

BICAPI  void   delete_render_storage(
    void  *ptr );

BICAPI  void  render_volume_to_slice(
    int             n_dims1,
    int             sizes1[],
    void            *volume_data1,
    Data_types      volume1_type,
    int             n_slices1,
    Real            weights1[],
    int             strides1[],
    Real            **origins1,
    Real            x_axis1[],
    Real            y_axis1[],
    int             n_dims2,
    int             sizes2[],
    void            *volume_data2,
    Data_types      volume2_type,
    int             n_slices2,
    Real            weights2[],
    int             strides2[],
    Real            **origins2,
    Real            x_axis2[],
    Real            y_axis2[],
    int             x_pixel_start,
    int             x_pixel_end,
    int             y_pixel_start,
    int             y_pixel_end,
    unsigned short  **cmode_colour_map,
    Colour          **rgb_colour_map,
    Colour          empty_colour,
    void            *render_storage,
    pixels_struct   *pixels );

BICAPI  void  render_one_row (
    void            *volume_data1,
    Data_types      volume1_type,
    int             y,
    int             start_x,
    int             end_x,
    size_t          **y_offsets1,
    size_t          **row_offsets1,
    void            **start_slices1,
    int             n_slices1,
    Real            weights1[],
    void            *volume_data2,
    Data_types      volume2_type,
    size_t          **y_offsets2,
    size_t          **row_offsets2,
    void            **start_slices2,
    int             n_slices2,
    Real            weights2[],
    unsigned short  **cmode_colour_map,
    Colour          **rgb_colour_map,
    pixels_struct   *pixels );

BICAPI  void  initialize_resample_volume(
    resample_struct      *resample,
    Volume               src_volume,
    General_transform    *dest_to_src_transform,
    Volume               dest_volume );

BICAPI  BOOLEAN  do_more_resampling(
    resample_struct  *resample,
    Real             max_seconds,
    Real             *fraction_done );

BICAPI  void  resample_volume(
    Volume                   src_volume,
    General_transform        *dest_to_src_transform,
    Volume                   dest_volume );

BICAPI  void  scan_lines_to_voxels(
    lines_struct     *lines,
    Volume           volume,
    Volume           label_volume,
    int              label );

BICAPI  void  scan_marker_to_voxels(
    marker_struct    *marker,
    Volume           volume,
    Volume           label_volume,
    int              label );

BICAPI  void  scan_object_to_volume(
    object_struct    *object,
    Volume           volume,
    Volume           label_volume,
    int              label,
    Real             max_distance );

BICAPI  void  scan_polygons_to_voxels(
    polygons_struct     *polygons,
    Volume              volume,
    Volume              label_volume,
    int                 label,
    Real                max_distance );

BICAPI  void  scan_a_polygon(
    int                 size,
    Point               vertices[],
    Point               voxels[],
    int                 n_output_vertices,
    Point               output_vertices[],
    Volume              volume,
    Volume              label_volume,
    int                 sizes[],
    int                 label );

BICAPI  void  scan_quadmesh_to_voxels(
    quadmesh_struct     *quadmesh,
    Volume              volume,
    Volume              label_volume,
    int                 label,
    Real                max_distance );

BICAPI  Volume  smooth_resample_volume(
    Volume              volume,
    int                 new_nx,
    int                 new_ny,
    int                 new_nz );

BICAPI  void  convert_voxel_to_talairach(
    Real   x_voxel,
    Real   y_voxel,
    Real   z_voxel,
    int    nx_voxels,
    int    ny_voxels,
    int    nz_voxels,
    Real   *x_tal,
    Real   *y_tal,
    Real   *z_tal );

BICAPI  void  convert_talairach_to_voxel(
    Real   x_tal,
    Real   y_tal,
    Real   z_tal,
    int    nx_voxels,
    int    ny_voxels,
    int    nz_voxels,
    Real   *x_voxel,
    Real   *y_voxel,
    Real   *z_voxel );

BICAPI  void  convert_talairach_to_mm(
    Real   x_tal,
    Real   y_tal,
    Real   z_tal,
    Real   *x_mm,
    Real   *y_mm,
    Real   *z_mm );

BICAPI  void  convert_mm_to_talairach(
    Real   x_mm,
    Real   y_mm,
    Real   z_mm,
    Real   *x_tal,
    Real   *y_tal,
    Real   *z_tal );

#ifdef __cplusplus
}
#endif

#endif
