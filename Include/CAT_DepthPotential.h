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
 * \brief Public API for compute_depth_potential.
 *
 * This function is part of the CAT-Surface public library interface and is used by command-line tools.
 *
 * \param param (in/out) Parameter of compute_depth_potential.
 * \param double (in/out) Parameter of compute_depth_potential.
 * \return Return value of compute_depth_potential.
 */
double * compute_depth_potential( polygons_struct *, double);
double * compute_areas( int, Point [], int *, 
                              int **, int );
double * local_depth_potential( int, Point [], double *, struct csr_matrix *, 
                              double *, double, double );
private void gauss_seidel( int, int *, int *, double *, double *, 
                              double *, int, double, double, int );
private void init_csr_matrix( int, int *, int **, struct csr_matrix * );
private void free_csr_matrix( struct csr_matrix * );
private void assemble( double, int, int, struct csr_matrix * );
private void stable_normals( int, Point [], Vector [], int *, int ** );
private double * compute_mean_curvature( int, Point [], double *,
                              Vector [] , struct csr_matrix * );
private void cot_laplacian_operator( int, Point [], struct csr_matrix *,
                              int *, int ** );

#endif
