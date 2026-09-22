/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 *
 */

#ifndef _CAT_NLM_H_
#define _CAT_NLM_H_

#include <math.h>

/**
 * \brief Optimized blockwise non-local means filter for 3D images.
 *
 * \param ima   (in/out) input image volume, filtered in-place
 * \param v     (in)  search window half-size
 * \param f     (in)  patch window half-size
 * \param h     (in)  filtering parameter
 * \param sigma (in)  noise standard deviation for Rician correction
 * \param dims  (in)  volume dimensions [cols, rows, slices]
 */
void ornlm(float *ima, int v, int f, float h, float sigma, const int *dims);
/**
 * \brief Spatially adaptive non-local means filter for 3D images.
 *
 * \param ima         (in/out) input image volume
 * \param v           (in)  search window half-size
 * \param f           (in)  patch window half-size
 * \param use_rician  (in)  non-zero for Rician correction
 * \param strength    (in)  strength scaling for adaptive weights
 * \param dims        (in)  volume dimensions [x, y, z]
 */
void sanlm(float *ima, int v, int f, int use_rician, double strength,
           const int *dims);

#endif