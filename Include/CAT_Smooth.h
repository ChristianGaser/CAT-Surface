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
 * \brief Public API for get_all_polygon_point_neighbours.
 *
 * This function is part of the CAT-Surface public library interface and is used by command-line tools.
 *
 * \param param (in/out) Parameter of get_all_polygon_point_neighbours.
 * \param param (in/out) Parameter of get_all_polygon_point_neighbours.
 * \param param (in/out) Parameter of get_all_polygon_point_neighbours.
 * \return void (no return value).
 */
void get_all_polygon_point_neighbours(polygons_struct *, int *[], int **[]);
void heatkernel_blur_points(int, Point [], double [], int, int *, int, double,
                            Point *, double *);
/**
 * \brief Public API for smooth_heatkernel.
 *
 * This function is part of the CAT-Surface public library interface and is used by command-line tools.
 *
 * \param param (in/out) Parameter of smooth_heatkernel.
 * \param param (in/out) Parameter of smooth_heatkernel.
 * \param double (in/out) Parameter of smooth_heatkernel.
 * \return void (no return value).
 */
void smooth_heatkernel(polygons_struct *, double *, double);
void smooth_heatkernel_sharpness(polygons_struct *, double *, double, double);
/**
 * \brief Public API for smooth_laplacian.
 *
 * This function is part of the CAT-Surface public library interface and is used by command-line tools.
 *
 * \param polygons (in/out) Parameter of smooth_laplacian.
 * \param iter (in/out) Parameter of smooth_laplacian.
 * \param alpha (in/out) Parameter of smooth_laplacian.
 * \param beta (in/out) Parameter of smooth_laplacian.
 * \return Return value of smooth_laplacian.
 */
int smooth_laplacian(polygons_struct *polygons, int iter, double alpha, double beta);

#endif
