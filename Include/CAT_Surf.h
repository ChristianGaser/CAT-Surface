/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 *
 */

#ifndef _CAT_SURF_H_
#define _CAT_SURF_H_

#include <stdbool.h>
#include <float.h>
#include "bicpl.h"
#include "CAT_Vol.h"
#include "CAT_Math.h"
#include "CAT_Map.h"
#include "CAT_Smooth.h"
#include "CAT_Resample.h"
#include "CAT_Intersect.h"
#include "CAT_Curvature.h"
#include "CAT_Defect.h"
#include "CAT_Deform.h"
#include "CAT_SurfaceIO.h"
#include "CAT_SurfUtils.h"

#define BW 1024
#define FWHM 30.0
#define DATAFORMAT 1 /* 1 = double data, 0 = complex data */
#define _PI 3.14159265358979323846264338327510

/* edge incidence */
typedef struct { int a,b,c, idx; } keytri_t;
typedef struct {
    int a,b;      /* a<b undirected edge key */
    int tri;      /* triangle index */
    double area;  /* triangle area for ranking */
} edge_occ;

void indices(int *, int *, int *);
/**
 * \brief Public API for apply_warp.
 *
 * This function is part of the CAT-Surface public library interface and is used by command-line tools.
 *
 * \param param (in/out) Parameter of apply_warp.
 * \param param (in/out) Parameter of apply_warp.
 * \param param (in/out) Parameter of apply_warp.
 * \param param (in/out) Parameter of apply_warp.
 * \param int (in/out) Parameter of apply_warp.
 * \return void (no return value).
 */
void apply_warp(polygons_struct *, polygons_struct *, double *, int *, int);
/**
 * \brief Public API for apply_uv_warp.
 *
 * This function is part of the CAT-Surface public library interface and is used by command-line tools.
 *
 * \param param (in/out) Parameter of apply_uv_warp.
 * \param param (in/out) Parameter of apply_uv_warp.
 * \param param (in/out) Parameter of apply_uv_warp.
 * \param param (in/out) Parameter of apply_uv_warp.
 * \param int (in/out) Parameter of apply_uv_warp.
 * \return void (no return value).
 */
void apply_uv_warp(polygons_struct *, polygons_struct *, double *,
              double *, int );
/**
 * \brief Public API for get_surface_ratio.
 *
 * This function is part of the CAT-Surface public library interface and is used by command-line tools.
 *
 * \param double (in/out) Parameter of get_surface_ratio.
 * \param param (in/out) Parameter of get_surface_ratio.
 * \param int (in/out) Parameter of get_surface_ratio.
 * \return Return value of get_surface_ratio.
 */
double * get_surface_ratio(double, polygons_struct *, int);
double get_vertex_areas(polygons_struct *, double *);
/**
 * \brief Public API for get_area_of_points.
 *
 * This function is part of the CAT-Surface public library interface and is used by command-line tools.
 *
 * \param param (in/out) Parameter of get_area_of_points.
 * \param param (in/out) Parameter of get_area_of_points.
 * \return Return value of get_area_of_points.
 */
double get_area_of_points(polygons_struct *, double *);
void get_radius_of_points(polygons_struct *, double *);
void localstat_surface_double(polygons_struct *polygons,
                              double *input,
                              unsigned char *mask,
                              int stat_func,
                              int iters);
void correct_bounds_to_target(polygons_struct *, polygons_struct *);
/**
 * \brief Public API for correct_bounds_to_target_with_scaling.
 *
 * This function is part of the CAT-Surface public library interface and is used by command-line tools.
 *
 * \param param (in/out) Parameter of correct_bounds_to_target_with_scaling.
 * \param param (in/out) Parameter of correct_bounds_to_target_with_scaling.
 * \return void (no return value).
 */
void correct_bounds_to_target_with_scaling(polygons_struct *, polygons_struct *);
double get_sphere_radius(polygons_struct *);
/**
 * \brief Public API for compute_point_hausdorff.
 *
 * This function is part of the CAT-Surface public library interface and is used by command-line tools.
 *
 * \param param (in/out) Parameter of compute_point_hausdorff.
 * \param param (in/out) Parameter of compute_point_hausdorff.
 * \param param (in/out) Parameter of compute_point_hausdorff.
 * \param int (in/out) Parameter of compute_point_hausdorff.
 * \return Return value of compute_point_hausdorff.
 */
double compute_point_hausdorff(polygons_struct *, polygons_struct *, double *, int);
/**
 * \brief Public API for compute_exact_hausdorff.
 *
 * This function is part of the CAT-Surface public library interface and is used by command-line tools.
 *
 * \param param (in/out) Parameter of compute_exact_hausdorff.
 * \param p (in/out) Parameter of compute_exact_hausdorff.
 * \param param (in/out) Parameter of compute_exact_hausdorff.
 * \param int (in/out) Parameter of compute_exact_hausdorff.
 * \return Return value of compute_exact_hausdorff.
 */
double compute_exact_hausdorff(polygons_struct *, polygons_struct *p, double *, int);
/**
 * \brief Public API for compute_point_distance.
 *
 * This function is part of the CAT-Surface public library interface and is used by command-line tools.
 *
 * \param param (in/out) Parameter of compute_point_distance.
 * \param param (in/out) Parameter of compute_point_distance.
 * \param param (in/out) Parameter of compute_point_distance.
 * \param int (in/out) Parameter of compute_point_distance.
 * \return Return value of compute_point_distance.
 */
double compute_point_distance(polygons_struct *, polygons_struct *, double *, int);
/**
 * \brief Public API for compute_point_distance_mean.
 *
 * This function is part of the CAT-Surface public library interface and is used by command-line tools.
 *
 * \param param (in/out) Parameter of compute_point_distance_mean.
 * \param param (in/out) Parameter of compute_point_distance_mean.
 * \param param (in/out) Parameter of compute_point_distance_mean.
 * \param int (in/out) Parameter of compute_point_distance_mean.
 * \return Return value of compute_point_distance_mean.
 */
double compute_point_distance_mean(polygons_struct *, polygons_struct *, double *, int);
/**
 * \brief Public API for euler_characteristic.
 *
 * This function is part of the CAT-Surface public library interface and is used by command-line tools.
 *
 * \param param (in/out) Parameter of euler_characteristic.
 * \param verbose (in/out) Parameter of euler_characteristic.
 * \return Return value of euler_characteristic.
 */
int euler_characteristic(polygons_struct *, int verbose);
/**
 * \brief Public API for convert_ellipsoid_to_sphere_with_surface_area.
 *
 * This function is part of the CAT-Surface public library interface and is used by command-line tools.
 *
 * \param param (in/out) Parameter of convert_ellipsoid_to_sphere_with_surface_area.
 * \param double (in/out) Parameter of convert_ellipsoid_to_sphere_with_surface_area.
 * \return void (no return value).
 */
void convert_ellipsoid_to_sphere_with_surface_area(polygons_struct *, double);
void linear_smoothing(polygons_struct *, double, int, int, int *, int);
/**
 * \brief Public API for areal_smoothing.
 *
 * This function is part of the CAT-Surface public library interface and is used by command-line tools.
 *
 * \param param (in/out) Parameter of areal_smoothing.
 * \param double (in/out) Parameter of areal_smoothing.
 * \param int (in/out) Parameter of areal_smoothing.
 * \param int (in/out) Parameter of areal_smoothing.
 * \param param (in/out) Parameter of areal_smoothing.
 * \param int (in/out) Parameter of areal_smoothing.
 * \return void (no return value).
 */
void areal_smoothing(polygons_struct *, double, int, int, int *, int);
void distance_smoothing(polygons_struct *, double, int, int, int *, int);
void inflate_surface_and_smooth_fingers(polygons_struct *, const int,
                                        const double, const int, const double,
                                        const double, const double, const int);                                        
object_struct ** central_to_new_pial(polygons_struct *, double *, double *, 
                                    int, double, int, int);
void central_to_pial(polygons_struct *, double *, double *, int, double, int, int);
/**
 * \brief Public API for get_area_of_points_central_to_pial.
 *
 * This function is part of the CAT-Surface public library interface and is used by command-line tools.
 *
 * \param param (in/out) Parameter of get_area_of_points_central_to_pial.
 * \param param (in/out) Parameter of get_area_of_points_central_to_pial.
 * \param param (in/out) Parameter of get_area_of_points_central_to_pial.
 * \param double (in/out) Parameter of get_area_of_points_central_to_pial.
 * \return Return value of get_area_of_points_central_to_pial.
 */
double get_area_of_points_central_to_pial(polygons_struct *, double *, double *, double);
int correct_mesh_folding(polygons_struct *, polygons_struct *, float *, 
                         nifti_image *, double);
/**
 * \brief Public API for reduce_mesh_quadrics.
 *
 * This function is part of the CAT-Surface public library interface and is used by command-line tools.
 *
 * \param param (in/out) Parameter of reduce_mesh_quadrics.
 * \param int (in/out) Parameter of reduce_mesh_quadrics.
 * \param double (in/out) Parameter of reduce_mesh_quadrics.
 * \param int (in/out) Parameter of reduce_mesh_quadrics.
 * \param int (in/out) Parameter of reduce_mesh_quadrics.
 * \return Return value of reduce_mesh_quadrics.
 */
int reduce_mesh_quadrics(polygons_struct *, int, double, int, int);

#endif
