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
double get_median_double(double *arr, int n, int mask_zeros);
double get_std_double(double *arr, int n, int mask_zeros);
double get_sum_double(double *arr, int n, int mask_zeros);
int min_inplace_with_index_float(float *A, float *B, int n, signed char *index);
double get_min(void *arr, int n, int mask_zeros, int datatype);
double get_max(void *arr, int n, int mask_zeros, int datatype);
double get_mean(void *arr, int n, int mask_zeros, int datatype);
double get_median(void *arr, int n, int mask_zeros, int datatype);
double get_std(void *arr, int n, int mask_zeros, int datatype);
double get_sum(void *arr, int n, int mask_zeros, int datatype);
double get_masked_mean_array(void *arr, int n, unsigned char *mask, int datatype);
double get_masked_std_array(void *arr, int n, unsigned char *mask, int datatype);
void get_prctile(void *data, int n, double threshold[2], double prctile[2], int exclude_zeros, int datatype);
double get_corrcoef(void *x, void *y, int n, int exclude_zeros, int datatype);
void clip_data(void *data, int n, double lower_limit, double upper_limit, int datatype);

#endif
