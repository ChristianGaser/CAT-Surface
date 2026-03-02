/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 *
 */

#ifndef _CAT_MATH_H_
#define _CAT_MATH_H_

#include <stdio.h>
#include <math.h>
#include <memory.h>
#include <float.h>
#include <limits.h>
#include "CAT_NiftiLib.h"

#define SQRT2PI 2.506628
#define TOLSVD 1e-10
#define EPS 1e-15

enum
{
    F_MEAN,
    F_MIN,
    F_MAX,
    F_STD,
    F_SUM,
    F_MAXABS,
    F_EXP,
    F_MEDIAN,
    F_RANGE,
    F_COUNT,
    F_WAVERAGE,
    F_MULTI,
    F_CLOSE,
    F_OPEN
};

#ifndef isfinite
#define isfinite(x) ((x) * (x) >= 0.) /* check for NaNs */
#endif

#ifndef SQR
#define SQR(x) ((x)*(x))
#endif

#ifndef MAX
#define MAX(A,B) ((A) > (B) ? (A) : (B))
#endif

#ifndef MIN
#define MIN(A,B) ((A) < (B) ? (A) : (B))
#endif

#ifndef ROUND
#define ROUND( x ) ((long) ((x) + ( ((x) >= 0) ? 0.5 : (-0.5) ) ))
#endif

#ifndef isnan
#define isnan(a) ((a)!=(a)) 
#endif

#ifdef _MSC_VER
  static const unsigned long __nan[2] = {0xffffffff, 0x7fffffff};
  #define FNAN (*(const float *) __nan)
#else
  #define FNAN 0.0f/0.0f
#endif

int pinv(int m, int n, double **A, double **Ainv);
void convert_input_type(void *data, double *buffer, int n, int datatype);
void convert_output_type(void *data, double *buffer, int n, int datatype);
void convert_input_type_float(void *data, float *buffer, int n, int datatype);
void convert_output_type_float(void *data, float *buffer, int n, int datatype);
void normalize_double(double *arr, int n);
double get_min_double(double *arr, int n, int mask_zeros);
double get_max_double(double *arr, int n, int mask_zeros);
double get_mean_double(double *arr, int n, int mask_zeros);
/**
 * \brief Public API for get_median_double.
 *
 * This function is part of the CAT-Surface public library interface and is used by command-line tools.
 *
 * \param arr (in/out) Parameter of get_median_double.
 * \param n (in/out) Parameter of get_median_double.
 * \param mask_zeros (in/out) Parameter of get_median_double.
 * \return Return value of get_median_double.
 */
double get_median_double(double *arr, int n, int mask_zeros);
/**
 * \brief Public API for get_std_double.
 *
 * This function is part of the CAT-Surface public library interface and is used by command-line tools.
 *
 * \param arr (in/out) Parameter of get_std_double.
 * \param n (in/out) Parameter of get_std_double.
 * \param mask_zeros (in/out) Parameter of get_std_double.
 * \return Return value of get_std_double.
 */
double get_std_double(double *arr, int n, int mask_zeros);
double get_sum_double(double *arr, int n, int mask_zeros);
/**
 * \brief Public API for get_min.
 *
 * This function is part of the CAT-Surface public library interface and is used by command-line tools.
 *
 * \param arr (in/out) Parameter of get_min.
 * \param n (in/out) Parameter of get_min.
 * \param mask_zeros (in/out) Parameter of get_min.
 * \param datatype (in/out) Parameter of get_min.
 * \return Return value of get_min.
 */
double get_min(void *arr, int n, int mask_zeros, int datatype);
/**
 * \brief Public API for get_max.
 *
 * This function is part of the CAT-Surface public library interface and is used by command-line tools.
 *
 * \param arr (in/out) Parameter of get_max.
 * \param n (in/out) Parameter of get_max.
 * \param mask_zeros (in/out) Parameter of get_max.
 * \param datatype (in/out) Parameter of get_max.
 * \return Return value of get_max.
 */
double get_max(void *arr, int n, int mask_zeros, int datatype);
double get_mean(void *arr, int n, int mask_zeros, int datatype);
double get_median(void *arr, int n, int mask_zeros, int datatype);
double get_std(void *arr, int n, int mask_zeros, int datatype);
double get_sum(void *arr, int n, int mask_zeros, int datatype);
double get_masked_mean_array(void *arr, int n, unsigned char *mask, int datatype);
double get_masked_std_array(void *arr, int n, unsigned char *mask, int datatype);
void get_prctile_double(double *data, int n, double threshold[2],
                        double prctile[2], int exclude_zeros);

/**
 * \brief Public API for get_prctile.
 *
 * This function is part of the CAT-Surface public library interface and is used by command-line tools.
 *
 * \param data (in/out) Parameter of get_prctile.
 * \param n (in/out) Parameter of get_prctile.
 * \param threshold (in/out) Parameter of get_prctile.
 * \param prctile (in/out) Parameter of get_prctile.
 * \param exclude_zeros (in/out) Parameter of get_prctile.
 * \param datatype (in/out) Parameter of get_prctile.
 * \return void (no return value).
 */
void get_prctile(void *data, int n, double threshold[2], double prctile[2], int exclude_zeros, int datatype);
/**
 * \brief Public API for get_corrcoef.
 *
 * This function is part of the CAT-Surface public library interface and is used by command-line tools.
 *
 * \param x (in/out) Parameter of get_corrcoef.
 * \param y (in/out) Parameter of get_corrcoef.
 * \param n (in/out) Parameter of get_corrcoef.
 * \param exclude_zeros (in/out) Parameter of get_corrcoef.
 * \param datatype (in/out) Parameter of get_corrcoef.
 * \return Return value of get_corrcoef.
 */
double get_corrcoef(void *x, void *y, int n, int exclude_zeros, int datatype);
/**
 * \brief Public API for clip_data.
 *
 * This function is part of the CAT-Surface public library interface and is used by command-line tools.
 *
 * \param data (in/out) Parameter of clip_data.
 * \param n (in/out) Parameter of clip_data.
 * \param lower_limit (in/out) Parameter of clip_data.
 * \param upper_limit (in/out) Parameter of clip_data.
 * \param datatype (in/out) Parameter of clip_data.
 * \return void (no return value).
 */
void clip_data(void *data, int n, double lower_limit, double upper_limit, int datatype);

#endif
