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

void Bmap(float *src, unsigned char *label, unsigned char *prob, int n_classes, int BG, int niters, int nflips, int a, int b, int c, float *bias, float bias_thresh, int *dims, int pve);

#endif