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

Vector projectToPlane(Vector projected, Vector basis[2]);
Vector projection(Vector vector, Vector normal);

void leastSquares(const int num, Vector dc[], Vector dn[], double *k1, 
                                      double *k2);
void compute_points_centroid_and_normal_cg(polygons_struct *, int, int, int [],
                                      Point *, Vector *, double *, int,
                                      double *);
void get_polygon_vertex_curvatures_cg(polygons_struct *, int [], int *[],
                                      double, int, double []);
void get_smoothed_curvatures(polygons_struct *, double *,
                                      double, int);
/**
 * \brief Public API for compute_sulcus_depth.
 *
 * This function is part of the CAT-Surface public library interface and is used by command-line tools.
 *
 * \param param (in/out) Parameter of compute_sulcus_depth.
 * \param param (in/out) Parameter of compute_sulcus_depth.
 * \return void (no return value).
 */
void compute_sulcus_depth(polygons_struct *, double *);
/**
 * \brief Public API for compute_local_sharpness.
 *
 * This function is part of the CAT-Surface public library interface and is used by command-line tools.
 *
 * \param param (in/out) Parameter of compute_local_sharpness.
 * \param int (in/out) Parameter of compute_local_sharpness.
 * \param param (in/out) Parameter of compute_local_sharpness.
 * \param param (in/out) Parameter of compute_local_sharpness.
 * \return void (no return value).
 */
void compute_local_sharpness(polygons_struct *, int [], int * [], double *);
/**
 * \brief Public API for compute_convexity.
 *
 * This function is part of the CAT-Surface public library interface and is used by command-line tools.
 *
 * \param param (in/out) Parameter of compute_convexity.
 * \param int (in/out) Parameter of compute_convexity.
 * \param param (in/out) Parameter of compute_convexity.
 * \param param (in/out) Parameter of compute_convexity.
 * \return void (no return value).
 */
void compute_convexity(polygons_struct *, int [], int *[], double *);

#endif
