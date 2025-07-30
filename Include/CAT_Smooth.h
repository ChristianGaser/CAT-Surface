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

void get_all_polygon_point_neighbours(polygons_struct *, int *[], int **[]);
void heatkernel_blur_points(int, Point [], double [], int, int *, int, double,
                            Point *, double *);
void smooth_heatkernel(polygons_struct *, double *, double);
void smooth_heatkernel_sharpness(polygons_struct *, double *, double, double);
int smooth_laplacian(polygons_struct *polygons, int iter, double alpha, double beta);

#endif
