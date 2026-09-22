/* Rachel Yotter - rachel.yotter@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 *
 */

#ifndef _CAT_RESAMPLE_H_
#define _CAT_RESAMPLE_H_

#include <bicpl.h>
#include <limits.h>

#include "CAT_Surf.h"
#include "CAT_SurfaceIO.h"

/**
 * \brief Normalize source and target spheres to consistent scale and center.
 *
 * \param source_sphere     (in)  source spherical mesh (unmodified)
 * \param target_sphere     (in)  target spherical mesh (unmodified)
 * \param out_source_sphere (out) scaled and centered source sphere copy
 * \param out_target_sphere (out) scaled and centered target sphere copy
 */
void correct_shift_scale_sphere(polygons_struct *source_sphere,
                                polygons_struct *target_sphere,
                                polygons_struct **out_source_sphere,
                                polygons_struct **out_target_sphere);
/**
 * \brief Resample surface mesh to a tetrahedral sphere using barycentric interpolation.
 *
 * \param surface     (in)  source surface mesh with data
 * \param sphere      (in)  source or reference spherical mesh
 * \param n_triangles (in)  number of triangles for output tetrahedral sphere
 * \param invals      (in)  double[n_points]; surface data values (can be NULL)
 * \param outvals     (out) double[n_points]; interpolated values on output sphere
 * \return                  object_struct** containing resampled surface and sphere
 */
object_struct ** resample_surface(polygons_struct *surface,
                                  polygons_struct *sphere, int n_triangles,
                                  double *invals, double *outvals);
/**
 * \brief Resample surface mesh and its data to a target spherical coordinate system.
 *
 * \param polygons            (in)  source surface mesh with data to map
 * \param polygons_sphere     (in)  source sphere (can be NULL for implicit tetrahedral)
 * \param target_sphere       (in)  target spherical coordinate system
 * \param input_values        (in)  double[n_points]; scalar data on source surface
 * \param output_values       (out) double[n_points]; interpolated values on output surface
 * \param label_interpolation (in)  0/1; use per-triangle label refinement if 1
 * \param areal_interpolation (in)  0/1; use area-weighted interpolation if 1
 * \return                          object_struct** containing resampled surface mesh
 */
object_struct ** resample_surface_to_target_sphere(polygons_struct *polygons,
                                                   polygons_struct *polygons_sphere,
                                                   polygons_struct *target_sphere,
                                                   double *input_values,
                                                   double *output_values,
                                                   int label_interpolation,
                                                   int areal_interpolation);
/**
 * \brief Interpolate scalar values from source sphere vertices to target sphere vertices.
 *
 * \param source_sphere       (in)  source spherical mesh with data
 * \param target_sphere       (in)  target spherical mesh (aligned)
 * \param invals              (in)  double[n_points]; source vertex values
 * \param outvals             (out) double[n_points]; interpolated target vertex values
 * \param areal_interpolation (in)  0/1; apply area-weighted interpolation if 1
 */
void resample_values_sphere_noscale(polygons_struct *source_sphere,
                                    polygons_struct *target_sphere,
                                    double *invals, double *outvals,
                                    int areal_interpolation);
/**
 * \brief Interpolate scalar values from source to target sphere with optional alignment.
 *
 * \param source_sphere        (in)  source spherical mesh
 * \param target_sphere        (in)  target spherical mesh
 * \param invals               (in)  double[n_points]; source vertex values
 * \param outvals              (out) double[n_points]; interpolated target values
 * \param scale_and_shift      (in)  0/1; normalize sphere alignment if 1
 * \param areal_interpolation  (in)  0/1; use area-weighted interpolation if 1
 */
void resample_values_sphere(polygons_struct *source_sphere,
                            polygons_struct *target_sphere, double *invals,
                            double *outvals, int scale_and_shift,
                            int areal_interpolation);
/**
 * \brief Resample surface data onto a tetrahedral sphere mesh with parallel processing.
 *
 * \param polygons         (in)  source surface mesh with data to map
 * \param poly_src_sphere  (in)  source spherical mesh aligned with polygons
 * \param resampled_source (out) output tetrahedral sphere mesh (modified in-place)
 * \param input_values     (in)  double[n_points]; surface data values (can be NULL)
 * \param output_values    (out) double[n_points]; interpolated values on output mesh
 * \param n_triangles      (in)  desired number of triangles in tetrahedral output sphere
 */
void resample_spherical_surface(polygons_struct *polygons,
                                polygons_struct *poly_src_sphere,
                                polygons_struct *resampled_source,
                                double *input_values, double *output_values,
                                int n_triangles);

#endif
