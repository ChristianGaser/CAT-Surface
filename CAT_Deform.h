#ifndef  DEF_DEFORM
#define  DEF_DEFORM

#include <bicpl.h>

#define  MAX_IN_VOXEL_COEF_LOOKUP  10000
#define  N_DEFORM_HISTOGRAM   7

typedef enum { TOWARDS_LOWER, TOWARDS_HIGHER, ANY_DIRECTION }
               Normal_directions;
typedef enum { VOLUME_DATA }
               Deform_data_types;
typedef enum { FLAT_MODEL, AVERAGE_MODEL,
               PARAMETRIC_MODEL, GENERAL_MODEL }
               Deformation_model_types;

typedef struct {
        Real               min_isovalue;
        Real               max_isovalue;
        Real               gradient_thresh;
        Real               min_dot_product;
        Real               max_dot_product;
        Normal_directions  normal_direction;
        Real               tolerance;
} bound_def_struct;

typedef struct {
        Deform_data_types   type;
        Volume              volume;
        Volume              label_volume;
} deform_data_struct;

typedef struct {
        int                       up_to_n_pts;
        Deformation_model_types   model_type;
        Real                      model_weight;
        object_struct             *model_object;

        int                       n_model_pts;
        Point                     *model_centroids;
        Vector                    *model_normals;
        Point                     *model_pts;

        Real                      min_curv_offset;
        Real                      max_curv_offset;
} deform_mod_struct; 

typedef struct {
        int                       n_models;
        deform_mod_struct         *models;
        BOOLEAN                   position_constrained;
        Real                      max_position_offset;
        Point                     *original_positions;
} deform_model_struct;

typedef struct {
        deform_data_struct            deform_data;
        deform_model_struct           deform_model;
        Real                          fractional_step;
        Real                          max_step;
        Real                          max_search_dist;
        int                           degrees_continuity;
        bound_def_struct              bound_def;
        int                           max_iters;
        Real                          stop_thresh;

        int                           n_movements_alloced;
        double                        *prev_movements;
        Real                          movement_thresh;
} deform_struct;

typedef struct {
        int                            hash_key;
        Real                           coeffs[8];
        struct voxel_lin_coeff_struct  *prev;
        struct voxel_lin_coeff_struct  *next;
} voxel_lin_coeff_struct;

typedef struct {
        hash_table_struct      hash;
        int                    n_in_hash;
        voxel_lin_coeff_struct  *head;
        voxel_lin_coeff_struct  *tail;
} voxel_coeff_struct;

typedef struct {
        Real    average;
        Real    maximum;
        int     n_below[N_DEFORM_HISTOGRAM];
} deform_stats;

typedef struct {
        int                  axis;
        Point                *save_pts;
        Real                 *curv_factors;
        Point                *equil_pts;
        Point                *new_equil_pts;
        Point                *bound_pts;
        Point                *new_bound_pts;
        Real                 temp;
        Real                 temp_factor;
        int                  temp_step;
        int                  min_n_to_move;
        int                  max_n_to_move;
        Real                 max_trans;
        Real                 max_angle_rotation;
        Real                 max_scale_offset;
        int                  stop_criteria;
        int                  try;
        int                  max_tries;
        int                  max_successes;
        int                  n_successes;
        int                  n_pos_successes;
        int                  n_no_moves;
        Real                 min_delta_energy;
        Real                 max_delta_energy;
        Real                 energy;
} anneal_struct;

#include "CAT_DeformPrototypes.h"


#endif
