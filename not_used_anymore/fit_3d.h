/*
   Copyright Alan C. Evans
   Professor of Neurology
   McGill University
*/
#ifndef  _DEF_FIT_3D_H
#define  _DEF_FIT_3D_H

#include  <volume_io/internal_volume_io.h>
#include  <bicpl.h>
#include  <conjugate.h>

#undef  USE_STRETCH_SQRT
#define  USE_STRETCH_SQRT

#undef  USE_CURV_SQRT
#define  USE_CURV_SQRT

#define  THIS_IS_UNIQUE_EDGE( node, neigh )    ((node) < (neigh))

typedef  enum  { BOUNDARY_IS_OUTSIDE, BOUNDARY_IS_INSIDE,
                 BOUNDARY_NOT_FOUND } Boundary_flags ;

typedef  enum  { TOWARDS_LOWER, TOWARDS_HIGHER, ANY_DIRECTION }
               Normal_directions;

struct  clip_struct;

typedef  struct  clip_struct  clip_struct;

typedef  struct  voxel_lin_coef_struct
{
    int                               hash_key;
    Real                              coefs[8];
    struct     voxel_lin_coef_struct  *prev;
    struct     voxel_lin_coef_struct  *next;
}
voxel_lin_coef_struct;

typedef struct
{
    unsigned char          ***voxel_ptr;
    int                    sizes[N_DIMENSIONS];
    Real                   voxel_to_real_values[256];
    int                    offset1;
    int                    offset2;
    int                    offset3;
    int                    offset4;
    int                    offset5;
    int                    offset6;
    int                    offset7;
    Real                   trans00;
    Real                   trans01;
    Real                   trans02;
    Real                   trans03;
    Real                   trans10;
    Real                   trans11;
    Real                   trans12;
    Real                   trans13;
    Real                   trans20;
    Real                   trans21;
    Real                   trans22;
    Real                   trans23;
} voxel_coef_struct;

typedef  struct
{
    int                     n_edges;
    int                     n_points;
    int                     n_polygons;
    int                     *n_neighbours;
    int                     **neighbours;
    int                     *triangles;
    Point                   *points;
    int                     n_midpoints;
    int                     *midpoints;
  // Added by June Sic Kim at 5/12/2002
    Real                    max_weight_value;
    Real                    *weights;
    int                     n_weights;
}
surface_struct;

typedef  struct
{
    Volume                  volume;
    voxel_coef_struct       *voxel_lookup;
    bitlist_3d_struct       done_bits;
    bitlist_3d_struct       surface_bits;
    Real                    threshold;
    Real                    image_weight_in;
    Real                    image_weight_out;
    Normal_directions       normal_direction;
    Real                    max_outward;
    Real                    max_inward;
    Real                    max_dist_threshold;
    Real                    max_dist_weight;
    int                     oversample;
    Real                    offset;
    BOOLEAN                 check_direction_flag;
    BOOLEAN                 normal_direction_only;
    BOOLEAN                 clip_to_surface;
  // Added by June Sic Kim at 5/12/2002
    Real                    max_weight_value;
    Real                    *weights;
    int                     n_weights;
} surface_bound_struct;

typedef  struct
{
    Volume                  volume;
    voxel_coef_struct       *voxel_lookup;
    int                     continuity;
    Real                    threshold;
    Real                    min_diff;
    Real                    max_diff;
    Real                    image_weight;
    Real                    max_diff_weight;
    Real                    differential_ratio;
    Real                    differential_offset;
} gradient_struct;

typedef  struct
{
    Volume                  t1;
    Volume                  cls;
    Real                    search_distance;
    Real                    search_increment;
    Real                    image_weight;
    int                     oversample;
    int                     image_type;
    Real *                  t1grad;
    int                     update;
    int *                   mask;
    STRING                  mask_filename;
} gw_gradient_struct;

typedef  struct
{
    Volume                  volume;
    voxel_coef_struct       *voxel_lookup;
    int                     continuity;
    Real                    threshold;
    Real                    min_diff;
    Real                    max_diff;
    Real                    image_weight;
    Real                    max_diff_weight;
    int                     oversample;
    Real                    differential_ratio;
    Real                    differential_offset;
} surface_value_struct;

typedef   struct
{
    int                n_points;
    int                n_edges;
    int                *n_neighbours;
    int                **neighbours;
    float              *model_lengths;
    Real               stretch_weight;
    Real               max_stretch_weight;
    Real               min_stretch;
    Real               max_stretch;
    Real               differential_ratio;
    Real               differential_offset;
} stretch_struct;

typedef   struct
{
    int                n_points;
    int                *n_neighbours;
    int                **neighbours;
    float              *curvatures;
    Real               curvature_weight;
    Real               max_curvature_weight;
    Real               min_curvature;
    Real               max_curvature;
} curvature_struct;

typedef   struct
{
    int                n_points;
    int                *n_neighbours;
    int                **neighbours;
    float              *model_x;
    float              *model_y;
    Real               bend_weight;
    Real               max_bend_weight;
    Real               min_bend;
    Real               max_bend;
} bend_struct;

typedef struct
{
    int          surface1_point;
    int          surface2_point;
    Real         desired_distance;
    Real         min_distance1;
    Real         min_distance2;
    Real         max_distance1;
    Real         max_distance2;
} connection_struct;

typedef   struct
{
    int                surface_index1;
    int                surface_index2;
    Real               weight;
    Real               max_weight1;
    Real               max_weight2;
    int                n_weight_steps;
    int                n_connections;
    connection_struct  *connections;
    int                oversample;
} inter_surface_struct;

typedef struct
{
    int          surface_point;
    Point        anchor_point;
    Real         min_distance;
    Real         desired_distance;
    Real         max_distance;
  // Added by June Sic Kim at 6/12/2002
    Real         weight;
    Real         max_weight;
}
anchor_point_struct;

typedef   struct
{
    Real                 weight;
    Real                 max_dist_weight;
    int                  n_anchor_points;
    anchor_point_struct  *anchor_points;
  // Added by June Sic Kim at 7/12/2002
    Real                 max_weight_value;
} anchor_struct;


typedef struct
{
    int          n_surface_points;
    int          *surface_points;
    Real         *surface_weights;
    Point        anchor_point;
    Real         min_distance;
    Real         desired_distance;
    Real         max_distance;
} weight_struct;

typedef   struct
{
    Real                 weight;
    Real                 max_dist_weight;
    int                  n_weight_points;
    weight_struct        *weight_points;
} weight_point_struct;


typedef struct
{
    int                    n_weights;
    Real                   *weights;
    Real                   *min_distances;
    BOOLEAN                use_tri_tri_dist;
    BOOLEAN                square_flag;
} self_intersect_struct;

typedef struct
{
    int                    surface_index1;
    int                    surface_index2;
    int                    n_weights;
    Real                   *weights;
    Real                   *min_distances;
} surf_surf_struct;

typedef  struct
{
    Real                   weight;
    Real                   max_weight;
    Real                   adaptive_anchor_ratio;
    Real                   adaptive_boundary_ratio;
    int                    direction; // -1: inward, 1: outward
}volume_info_struct;

typedef  struct
{
    Real                   weight;
    Volume                 volume;
    Volume                 gradient_volume;
    voxel_coef_struct      *voxel_lookup;
    Real                   from_value;
    Real                   to_value;
    Real                   deriv_factor;
    Real                   oversample;
    int                    direction; // -1: inward, 1: outward
    int                    type;      // 0: white, 1: gray
}laplacian_struct;

typedef  struct
{
    object_struct          **objects;
    int                    n_objects;
    int                    *intersected;
    int                    n_intersected;
    Volume                 wm_masked;
    voxel_coef_struct      *voxel_lookup;
    Real                   weight;
    int                    direction; // -1: inward, 1: outward
}intersect_wm_struct;

typedef  struct
{
    BOOLEAN                static_flag;
    surface_struct         surface;
    int                    n_bound;
    surface_bound_struct   *bound;
    int                    n_gradient;
    gradient_struct        *gradient;
    int                    n_gw_gradient;
    gw_gradient_struct     *gw_gradient;
    int                    n_value;
    surface_value_struct   *value;
    int                    n_stretch;
    stretch_struct         *stretch;
    int                    n_curvature;
    curvature_struct       *curvature;
    int                    n_bend;
    bend_struct            *bend;
    int                    n_anchors;
    anchor_struct          *anchors;
    int                    n_weight_points;
    weight_point_struct    *weight_points;
    int                    n_self_intersects;
    self_intersect_struct  *self_intersects;
    int                    n_volume;
    volume_info_struct     *volume;
    int                    n_laplacian;
    laplacian_struct       *laplacian;
} one_surface_struct;

typedef  struct
{
    int                    n_surfaces;
    one_surface_struct     *surfaces;
    int                    n_inter_surfaces;
    inter_surface_struct   *inter_surfaces;
    int                    n_surf_surfs;
    surf_surf_struct       *surf_surfs;
  // Added by June
    intersect_wm_struct    *intersect_wm;
    int                    n_intersect_wm;
} Deform_struct;

typedef struct
{
    Real      boundary_fit;
    Real      gradient_fit;
    Real      gw_gradient_fit;
    Real      value_fit;
    Real      stretch_fit;
    Real      curvature_fit;
    Real      bend_fit;
    Real      self_intersect_fit;
    Real      closest_self_intersect;
    Real      surf_surf_fit;
    Real      closest_surf_surf;
    Real      inter_surface_fit;
    Real      anchor_fit;
    Real      weight_point_fit;
    Real      volume_fit;
    Real      laplacian_fit;
} fit_eval_struct;

typedef  struct {
    float   low_limits[N_DIMENSIONS];
    float   high_limits[N_DIMENSIONS];
    int     p1;
    int     n11;
    int     n12;
} poly_info_struct;

typedef struct
{
    int                    n_pairs;
    int                    n_pairs_alloc;
    int                    *p1s;
    int                    *p2s;
    unsigned char          *cases;
    float                  *min_line_dists;
} self_intersect_lookup_struct;

typedef struct
{
    int                    n_pairs;
    int                    n_pairs_alloc;
    int                    *p1s;
    int                    *p2s;
    unsigned char          *cases;
    float                  *min_line_dists;
} surf_surf_lookup_struct;

typedef struct {
    Real                          closest_dist;
    Real                          interval_size;
    Real                          max_movement;
    int                           point_grid_size;
    int                           *start_parameter;
    BOOLEAN                       self_intersect_present;
    self_intersect_lookup_struct  **si_lookups;
    BOOLEAN                       surf_surf_present;
    surf_surf_lookup_struct       *ss_lookups;
} line_lookup_struct;

#include  <fit_3d_prototypes.h>

#endif
