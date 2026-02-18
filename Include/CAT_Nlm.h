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
 * \brief Public API for ornlm.
 *
 * This function is part of the CAT-Surface public library interface and is used by command-line tools.
 *
 * \param ima (in/out) Parameter of ornlm.
 * \param v (in/out) Parameter of ornlm.
 * \param f (in/out) Parameter of ornlm.
 * \param h (in/out) Parameter of ornlm.
 * \param sigma (in/out) Parameter of ornlm.
 * \param dims (in/out) Parameter of ornlm.
 * \return void (no return value).
 */
void ornlm(float* ima, int v, int f, float h, float sigma, const int* dims);
/**
 * \brief Public API for sanlm.
 *
 * This function is part of the CAT-Surface public library interface and is used by command-line tools.
 *
 * \param ima (in/out) Parameter of sanlm.
 * \param v (in/out) Parameter of sanlm.
 * \param f (in/out) Parameter of sanlm.
 * \param is_rician (in/out) Parameter of sanlm.
 * \param strength (in/out) Parameter of sanlm.
 * \param dims (in/out) Parameter of sanlm.
 * \return void (no return value).
 */
void sanlm(float* ima, int v, int f, int is_rician, double strength, const int* dims);

#endif