/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 *
 */

#ifndef _CAT_BMAP_H_
#define _CAT_BMAP_H_

#include <math.h>

/**
 * \brief Bias-corrected tissue classification using EM-like updates.
 *
 * \param src      (in)  input intensity volume
 * \param label    (in/out) initial labels, updated in-place
 * \param prob     (out) class probability maps (stacked by class)
 * \param mean     (in/out) class mean estimates
 * \param n_classes (in) number of tissue classes
 * \param BG       (in)  background label threshold
 * \param niters   (in)  maximum number of iterations
 * \param a        (in)  half-window size along x for bias smoothing
 * \param b        (in)  half-window size along y for bias smoothing
 * \param c        (in)  half-window size along z for bias smoothing
 * \param bias     (in/out) bias field per voxel
 * \param dims     (in)  volume dimensions [nx, ny, nz]
 * \param pve      (in)  enable partial volume estimation if non-zero
 * \param verbose  (in)  non-zero to print progress
 */
void Bmap(float *src, unsigned char *label, unsigned char *prob, double *mean,
          int n_classes, int BG, int niters, int a, int b, int c, float *bias,
          int *dims, int pve, int verbose);

#endif
