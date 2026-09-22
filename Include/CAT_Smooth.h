/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 *
 */

#ifndef _CAT_SMOOTH_H_
#define _CAT_SMOOTH_H_

#include "CAT_SurfaceIO.h"

#ifndef isnan
#define isnan(a) ((a)!=(a)) 
#endif

#ifdef _MSC_VER
  static const unsigned long __nan[2] = {0xffffffff, 0x7fffffff};
  #define FNAN (*(const float *) __nan)
#else
  #define FNAN 0.0f/0.0f
#endif

/**
 * \brief Extract all unique vertex neighbors for every point on a polygon mesh.
 *
 * \param polygons                (in)  polygon mesh structure
 * \param n_point_neighbours_ptr (out)  int[n_points]; count of neighbors for each point
 * \param point_neighbours_ptr   (out)  int*[n_points]; neighbor vertex indices for each point
 */
void get_all_polygon_point_neighbours(polygons_struct *polygons,
                                      int *n_point_neighbours_ptr[],
                                      int **point_neighbours_ptr[]);
/**
 * \brief Apply weighted averaging for a point using heat kernel distance weighting.
 *
 * \param n_polygon_pts (in)  total number of points on mesh
 * \param polygon_pts   (in)  Point[n_polygon_pts]; 3D coordinates of all mesh points
 * \param values        (in)  double[n_polygon_pts]; scalar values at each point
 * \param n_neighbours  (in)  number of neighbors to consider
 * \param neighbours    (in)  int[n_neighbours]; neighbor vertex indices
 * \param ptidx         (in)  current point index
 * \param sigma         (in)  heat kernel bandwidth
 * \param smooth_point  (out) Point; smoothed position
 * \param value         (out) double; smoothed scalar value
 */
void heatkernel_blur_points(int n_polygon_pts, Point polygon_pts[],
                            double values[], int n_neighbours, int *neighbours,
                            int ptidx, double sigma, Point *smooth_point,
                            double *value);
/**
 * \brief Apply heat diffusion-based smoothing to polygon mesh scalar data.
 *
 * \param polygons  (in/out) polygon mesh structure; coordinates smoothed in-place
 * \param values    (in/out) double[n_points]; scalar values; smoothed in-place
 * \param fwhm      (in)     Gaussian smoothing FWHM in mm (converted to sigma internally)
 */
void smooth_heatkernel(polygons_struct *polygons, double *values, double fwhm);
/**
 * \brief HC Laplacian smoothing of a mesh (Vollmer et al., 1999).
 *
 * \param polygons (in/out) mesh, smoothed in-place
 * \param iter     (in)     number of iterations
 * \param alpha    (in)     weight of the original positions in the correction
 *                              (0..1; 0.1 is typical)
 * \param beta     (in)     weight of a point's own correction against that of
 *                              its neighbours (0..1; 0.5 is typical)
 * \return 0 on success, -1 on an empty mesh or allocation failure
 */
int smooth_laplacian(polygons_struct *polygons, int iter, double alpha,
                     double beta);

#endif
