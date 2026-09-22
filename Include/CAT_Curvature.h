/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 *
 */

#ifndef _CAT_CURVATURE_H_
#define _CAT_CURVATURE_H_

#include <bicpl.h>

/**
 * \brief Project a 3D vector onto a 2D plane defined by basis vectors.
 *
 * \param projected (in)  3D vector to project onto plane
 * \param basis     (in)  basis[2]; two orthonormal vectors defining the plane
 * \return               2D projection as Vector with z=0
 */
Vector projectToPlane(Vector projected, Vector basis[2]);

/**
 * \brief Compute vertex curvature type and centroid/normal of vertex neighborhood.
 *
 * \param polygons      (in)  surface mesh
 * \param pidx          (in)  vertex index
 * \param n_neighbours  (in)  number of neighbors
 * \param neighbours    (in)  int[n_neighbours]; neighbor vertex indices
 * \param centroid      (out) centroid of neighborhood
 * \param normal        (out) normal of neighborhood
 * \param baselen       (out) neighborhood size measure
 * \param curvtype      (in)  curvature metric type (1=Gaussian, 2=curvedness, 3=shape index, 4=mean, 6=bending energy, 7=sharpness, 8=folding, 9=min, 10=max)
 * \param curvparameter (out) computed curvature value
 */
void compute_points_centroid_and_normal_cg(polygons_struct *polygons, int pidx,
                                           int n_neighbours, int neighbours[],
                                           Point *centroid, Vector *normal,
                                           double *baselen, int curvtype,
                                           double *curvparameter);
/**
 * \brief Compute specified curvature metric for all vertices of a mesh.
 *
 * \param polygons            (in)  surface mesh
 * \param n_neighbours        (in)  int[n_points]; neighbors per vertex
 * \param neighbours          (in)  int*[n_points]; neighbor indices per vertex
 * \param smoothing_distance  (in)  FWHM for heat kernel smoothing (0=no smoothing)
 * \param curvtype            (in)  1=Gaussian, 2=curvedness, 3=shape index, 4=mean, 5=sulcal, 6=bending, 7=sharpness, 8=folding, 9=min, 10=max, >11=depth potential
 * \param curvatures          (out) double[n_points]; computed curvature values
 */
void get_polygon_vertex_curvatures_cg(polygons_struct *polygons,
                                      int n_neighbours[], int *neighbours[],
                                      double smoothing_distance, int curvtype,
                                      double curvatures[]);
/**
 * \brief Curvature of a surface, heat-kernel smoothed and scaled to [0, 1].
 *
 * \param polygons (in)  surface mesh
 * \param values   (out) n_points curvature values in [0, 1]
 * \param fwhm     (in)  FWHM of the heat-kernel smoothing in mm
 * \param curvtype (in)  curvature type, as in get_polygon_vertex_curvatures_cg()
 */
void get_smoothed_curvatures(polygons_struct *polygons, double *values,
                             double fwhm, int curvtype);
/**
 * \brief Compute sulcal depth using convex hull Euclidean distance.
 *
 * \param surface (in)  cortical surface mesh
 * \param depth   (out) double[n_points]; distance to convex hull for each vertex
 */
void compute_sulcus_depth(polygons_struct *surface, double *depth);
/**
 * \brief Compute FreeSurfer-style sulcal depth via iterative surface inflation.
 *
 * Inflates a copy of the surface and projects each vertex's displacement onto
 * its original surface normal. Provides a signed depth measure comparable to
 * FreeSurfer's sulc file (curvtype 11 in get_polygon_vertex_curvatures_cg).
 *
 * \param surface (in/out) surface mesh; centered in-place at its center of mass
 * \param depth   (out)    double[n_points]; signed displacement along surface normal
 */
void compute_sulcal_depth_inflation(polygons_struct *surface, double *depth);
/**
 * \brief Compute local sharpness metric (maximum angular variation at vertices).
 *
 * \param polygons    (in)  surface mesh
 * \param n_neighbours (in)  int[n_points]; neighbors per vertex
 * \param neighbours   (in)  int*[n_points]; neighbor vertex indices
 * \param sharpness    (out) double[n_points]; sharpness measure per vertex (degrees)
 */
void compute_local_sharpness(polygons_struct *polygons, int n_neighbours[],
                             int *neighbours[], double *sharpness);
/**
 * \brief Compute local convexity as projection of neighborhood vector onto normal.
 *
 * \param polygons     (in)  surface mesh
 * \param n_neighbours (in)  int[n_points]; neighbors per vertex
 * \param neighbours    (in)  int*[n_points]; neighbor vertex indices
 * \param convexity     (out) double[n_points]; convexity measure per vertex
 */
void compute_convexity(polygons_struct *polygons, int n_neighbours[],
                       int *neighbours[], double *convexity);

#endif
