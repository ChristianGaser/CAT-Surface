/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 *
 */

#ifndef _CAT_DEPTHPOTENTIAL_H_
#define _CAT_DEPTHPOTENTIAL_H_

#include <bicpl.h>

struct csr_matrix {
  int   n;        // size of matrix
  int   nnz;      // number of non-zero coeffs
  int * ia;       // row pointers
  int * ja;       // column pointers
  double * A;       // coefficients
};

/**
 * \brief Compute depth potential on a surface mesh.
 *
 * \param polygons (in)  input surface mesh
 * \param alpha    (in)  regularization weight for the Laplacian system
 * \return Allocated depth potential array (length n_points)
 */
double * compute_depth_potential(polygons_struct *polygons, double alpha);
/**
 * \brief Compute per-vertex mixed Voronoi areas.
 *
 * \param n_points (in)  number of vertices
 * \param coords   (in)  vertex coordinates
 * \param n_ngh    (in)  neighbor counts per vertex
 * \param ngh      (in)  ordered neighbor lists
 * \param lambda   (in)  area mode selector
 * \return Allocated per-vertex area array
 */
double * compute_areas(int n_points, Point coords[], int *n_ngh, int **ngh,
                       int lambda);
/**
 * \brief Solve the depth potential linear system.
 *
 * \param n_points (in)  number of vertices
 * \param coords   (in)  vertex coordinates
 * \param areas    (in)  per-vertex areas
 * \param mat      (in)  cotangent Laplacian matrix
 * \param mc       (in)  mean curvature values
 * \param alpha    (in)  regularization weight
 * \param SOR      (in)  successive over-relaxation factor
 * \return Allocated depth potential array
 */
double * local_depth_potential(int n_points, Point coords[], double *areas,
                               struct csr_matrix *mat, double *mc, double alpha,
                               double SOR);

#endif
