#ifndef  DEF_trans_prototypes
#define  DEF_trans_prototypes

#ifdef __cplusplus
extern "C" {
#endif

BICAPI  void  get_nonlinear_warp(
   Real     **positions,   /* n_points x n_dims */
   Real     **values,   /* n_points x n_values */
   Real     **INVMLY,   /* n_points+1+n_dims x n_values */
   int      n_points,
   int      n_dims,
   int      n_values );

BICAPI  void  compute_transform_from_tags(
    int                 npoints, 
    Real                **tag_list1, 
    Real                **tag_list2, 
    Trans_type          trans_type,
    General_transform   *transform);

BICAPI  void  build_transformation_matrix(
    Transform   *lt, 
    Real        center[],
    Real        translations[],
    Real        scales[],
    Real        shears[],
    Real        rotations[] );

BICAPI  void  transpose(
    int     rows,
    int     cols,
    Real    **mat,
    Real    **mat_transpose );

BICAPI  void  matrix_multiply(
    int    ldim,
    int    mdim,
    int    ndim, 
    Real   **Amat,
    Real   **Bmat,
    Real   **Cmat );

BICAPI  BOOLEAN  optimize_simplex(
    Real           **pts1,
    Real           **pts2, 
    int            npoints, 
    Trans_type     trans_type,
    Real           center[],
    Real           translations[],
    Real           scales[],
    Real           shears[],
    Real           rotations[] );

BICAPI  void  procrustes(
    int         npoints,
    int         ndim, 
    Real        **Apoints,
    Real        **Bpoints,
    Real        translation[],
    Real        centre_of_rotation[],
    Transform   *rotation_transform,
    Real        *scale_ptr );

BICAPI  BOOLEAN  rotmat_to_ang(
    Transform   *rot_trans,
    Real        *ang );

BICAPI void safe_compute_transform_from_tags(
    int                 npoints, 
    Real                **tag_list1, 
    Real                **tag_list2, 
    Trans_type          trans_type,
    General_transform   *transform );

BICAPI  void  make_scale_transform( 
    Real        sx,
    Real        sy,
    Real        sz,
    Transform   *transform );

BICAPI  void  set_transform_x_and_z_axes(
    Transform   *transform,
    Vector      *x_axis,
    Vector      *z_axis );

BICAPI  void  make_translation_transform(
    Real        x_trans,
    Real        y_trans,
    Real        z_trans,
    Transform   *transform );

BICAPI  void  make_origin_transform(
    Point       *origin,
    Transform   *transform );

BICAPI  void  make_rotation_transform(
    Real       radians,
    int        axis,
    Transform  *transform );

BICAPI  void  make_transform_relative_to_point(
    Point      *point,
    Transform  *transform,
    Transform  *rel_transform );

BICAPI  void  make_transform_in_coordinate_system(
    Point      *origin,
    Vector     *x_axis,
    Vector     *y_axis,
    Vector     *z_axis,
    Transform  *transform,
    Transform  *rel_transform );

BICAPI  void  make_rotation_about_axis(
    Vector     *axis,
    Real       angle,
    Transform  *transform );

BICAPI  void  convert_2d_transform_to_rotation_translation(
    Transform  *transform,
    Real       *degrees_clockwise,
    Real       *x_trans,
    Real       *y_trans );

BICAPI  Real  compute_clockwise_rotation( Real x, Real y );

BICAPI  BOOLEAN  is_transform_right_handed(
    Transform   *transform );

BICAPI  void  make_identity_transform_2d(
    Transform_2d *transform );

BICAPI  void  get_inverse_transform_2d(
    Transform_2d   *transform,
    Transform_2d   *inverse );

BICAPI  void  transform_point_2d(
    Transform_2d   *transform,
    Real           x,
    Real           y,
    Real           *x_trans,
    Real           *y_trans );

BICAPI  void  get_least_squares_transform_2d(
    int           n_points,
    Real          x[],
    Real          y[],
    Real          x_trans[],
    Real          y_trans[],
    Transform_2d  *transform_2d );

BICAPI  Status  write_transform_file(
    STRING       filename,
    STRING       comments,
    Transform    *transform );

BICAPI  Status  read_transform_file(
    STRING       filename,
    Transform    *transform );

#ifdef __cplusplus
}
#endif

#endif
