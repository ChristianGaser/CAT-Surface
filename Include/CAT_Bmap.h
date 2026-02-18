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
 * \brief Public API for Bmap.
 *
 * This function is part of the CAT-Surface public library interface and is used by command-line tools.
 *
 * \param src (in/out) Parameter of Bmap.
 * \param label (in/out) Parameter of Bmap.
 * \param prob (in/out) Parameter of Bmap.
 * \param mean (in/out) Parameter of Bmap.
 * \param n_classes (in/out) Parameter of Bmap.
 * \param BG (in/out) Parameter of Bmap.
 * \param niters (in/out) Parameter of Bmap.
 * \param a (in/out) Parameter of Bmap.
 * \param b (in/out) Parameter of Bmap.
 * \param c (in/out) Parameter of Bmap.
 * \param bias (in/out) Parameter of Bmap.
 * \param dims (in/out) Parameter of Bmap.
 * \param pve (in/out) Parameter of Bmap.
 * \param verbose (in/out) Parameter of Bmap.
 * \return void (no return value).
 */
void Bmap(float *src, unsigned char *label, unsigned char *prob, double *mean, int n_classes, int BG, int niters, int a, int b, int c, float *bias, int *dims, int pve, int verbose);

#endif
