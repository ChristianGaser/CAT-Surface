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

void ornlm(float* ima, int v, int f, float h, float sigma, const int* dims);
void sanlm(float* ima, int v, int f, int is_rician, double strength, const int* dims);

#endif