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

/**
 * \brief Moore-Penrose pseudo-inverse of an m x n matrix.
 *
 * \param m    (in)  number of rows of A
 * \param n    (in)  number of columns of A
 * \param A    (in)  m x n matrix (bicpl ALLOC2D layout)
 * \param Ainv (out) n x m pseudo-inverse, allocated by the caller
 * \return rank of A
 */
int pinv(int m, int n, double **A, double **Ainv);
/**
 * \brief Build an orthogonal polynomial basis, matching R's poly(x, degree).
 *
 * Centres \p x, forms the Vandermonde matrix [1, x, x^2, ..., x^degree] and
 * orthonormalises its columns with modified Gram-Schmidt; the constant
 * column is dropped (it is normally covered by an intercept).  The result
 * is identical to R's default (non-raw) orthogonal polynomials.
 *
 * \param x      (in)  input vector of length \p n
 * \param n      (in)  number of observations
 * \param degree (in)  polynomial degree (>= 1, and < number of distinct x)
 * \param out    (out) caller-allocated array of n*degree doubles, filled in
 *                     column-major order (column j starts at out + j*n)
 * \return 1 on success, 0 on failure (bad arguments or degenerate data)
 */
int orthogonal_poly(const double *x, int n, int degree, double *out);
/**
 * \brief Convert arbitrary datatype array to double-precision buffer.
 *
 * \param data     (in)  void pointer to input array; interpretation based on datatype
 * \param buffer   (out) double[n]; pre-allocated target array
 * \param n        (in)  number of elements to convert
 * \param datatype (in)  data type code (DT_UINT8, DT_UINT16, DT_FLOAT32, DT_FLOAT64, etc.)
 */
void convert_input_type(void *data, double *buffer, int n, int datatype);
/**
 * \brief Convert double-precision buffer back to arbitrary output datatype.
 *
 * \param data     (out) void pointer to output array; interpretation based on datatype
 * \param buffer   (in)  double[n]; source array with converted values
 * \param n        (in)  number of elements to convert back
 * \param datatype (in)  target data type code (DT_UINT8, DT_UINT16, DT_FLOAT32, etc.)
 */
void convert_output_type(void *data, double *buffer, int n, int datatype);
/**
 * \brief Convert arbitrary datatype array to single-precision float buffer.
 *
 * \param data     (in)  void pointer to input array; interpretation based on datatype
 * \param buffer   (out) float[n]; pre-allocated target array
 * \param n        (in)  number of elements to convert
 * \param datatype (in)  data type code (DT_UINT8, DT_UINT16, DT_FLOAT32, DT_FLOAT64, etc.)
 */
void convert_input_type_float(void *data, float *buffer, int n, int datatype);
/**
 * \brief Convert single-precision float buffer back to arbitrary output datatype.
 *
 * \param data     (out) void pointer to output array; interpretation based on datatype
 * \param buffer   (in)  float[n]; source array with converted values
 * \param n        (in)  number of elements to convert back
 * \param datatype (in)  target data type code (DT_UINT8, DT_UINT16, DT_FLOAT32, etc.)
 */
void convert_output_type_float(void *data, float *buffer, int n, int datatype);
/**
 * \brief Subtract mean from an array of doubles.
 *
 * \param arr Array of doubles.
 * \param n Number of elements in the array.
 */
void normalize_double(double *arr, int n);
/**
 * \brief Get minimum value from double array with optional zero exclusion.
 *
 * \param arr            (in)  double[n]; array to search
 * \param n              (in)  array size
 * \param exclude_zeros  (in)  if non-zero, zero values are ignored (use DBL_MAX as minimum)
 * \return               The minimum value (or minimum of non-zero values if exclude_zeros=1)
 */
double get_min_double(double *arr, int n, int exclude_zeros);
/**
 * \brief Get maximum value from double array with optional zero exclusion.
 *
 * \param arr            (in)  double[n]; array to search
 * \param n              (in)  array size
 * \param exclude_zeros  (in)  if non-zero, zero values are ignored
 * \return               The maximum value (or maximum of non-zero values if exclude_zeros=1)
 */
double get_max_double(double *arr, int n, int exclude_zeros);
/**
 * \brief Get mean value from double array with optional zero exclusion.
 *
 * \param arr            (in)  double[n]; array to compute mean from
 * \param n              (in)  array size
 * \param exclude_zeros  (in)  if non-zero, zero values are excluded from mean calculation
 * \return               The mean value (or mean of non-zero values if exclude_zeros=1)
 */
double get_mean_double(double *arr, int n, int exclude_zeros);
/**
 * \brief Get median value from double array with optional zero exclusion.
 *
 * \param arr            (in/out) double[n]; array to compute median from; sorted in-place
 * \param n              (in)     array size
 * \param exclude_zeros  (in)     if non-zero, zero values are ignored in median calculation
 * \return               The median value (or median of non-zero values if exclude_zeros=1)
 */
double get_median_double(double *arr, int n, int exclude_zeros);
/**
 * \brief Get standard deviation from double array with optional zero exclusion.
 *
 * \param arr            (in)  double[n]; array to compute std dev from
 * \param n              (in)  array size
 * \param exclude_zeros  (in)  if non-zero, zero values are excluded from calculation
 * \return               The standard deviation (or std dev of non-zero values)
 */
double get_std_double(double *arr, int n, int exclude_zeros);
/**
 * \brief Get sum of elements in double array with optional zero exclusion.
 *
 * \param arr            (in)  double[n]; array to sum
 * \param n              (in)  array size
 * \param exclude_zeros  (in)  if non-zero, zero values are excluded from sum
 * \return               The sum (or sum of non-zero values if exclude_zeros=1)
 */
double get_sum_double(double *arr, int n, int exclude_zeros);
/**
 * \brief Get minimum from arbitrary datatype array.
 *
 * \param data            (in)  void pointer to input array
 * \param n               (in)  array size
 * \param exclude_zeros   (in)  if non-zero, zeros are excluded (use DBL_MAX as min)
 * \param datatype        (in)  data type code (DT_UINT8, DT_FLOAT32, etc.)
 * \return                The minimum value
 */
double get_min(void *data, int n, int exclude_zeros, int datatype);
/**
 * \brief Get maximum from arbitrary datatype array.
 *
 * \param data            (in)  void pointer to input array
 * \param n               (in)  array size
 * \param exclude_zeros   (in)  if non-zero, zeros are excluded
 * \param datatype        (in)  data type code (DT_UINT8, DT_FLOAT32, etc.)
 * \return                The maximum value
 */
double get_max(void *data, int n, int exclude_zeros, int datatype);
/**
 * \brief Get mean from arbitrary datatype array.
 *
 * \param data            (in)  void pointer to input array
 * \param n               (in)  array size
 * \param exclude_zeros   (in)  if non-zero, zeros are excluded from mean
 * \param datatype        (in)  data type code (DT_UINT8, DT_FLOAT32, etc.)
 * \return                The mean value
 */
double get_mean(void *data, int n, int exclude_zeros, int datatype);
/**
 * \brief Get median from arbitrary datatype array.
 *
 * \param data            (in)  void pointer to input array
 * \param n               (in)  array size
 * \param exclude_zeros   (in)  if non-zero, zeros are excluded from median
 * \param datatype        (in)  data type code (DT_UINT8, DT_FLOAT32, etc.)
 * \return                The median value
 */
double get_median(void *data, int n, int exclude_zeros, int datatype);
/**
 * \brief Get standard deviation from arbitrary datatype array.
 *
 * \param data            (in)  void pointer to input array
 * \param n               (in)  array size
 * \param exclude_zeros   (in)  if non-zero, zeros are excluded from calculation
 * \param datatype        (in)  data type code (DT_UINT8, DT_FLOAT32, etc.)
 * \return                The standard deviation
 */
double get_std(void *data, int n, int exclude_zeros, int datatype);
/**
 * \brief Get sum from arbitrary datatype array.
 *
 * \param data            (in)  void pointer to input array
 * \param n               (in)  array size
 * \param exclude_zeros   (in)  if non-zero, zeros are excluded from sum
 * \param datatype        (in)  data type code (DT_UINT8, DT_FLOAT32, etc.)
 * \return                The sum value
 */
double get_sum(void *data, int n, int exclude_zeros, int datatype);
/**
 * \brief Mean of an array of any NIfTI datatype, optionally within a mask.
 *
 * \param data     (in)  array of n values of type datatype
 * \param n        (in)  number of values
 * \param mask     (in)  n mask values; only entries > 0 count (NULL: all)
 * \param datatype (in)  NIfTI datatype code of data (DT_FLOAT32, ...)
 * \return mean of the included values, NaN if there are none
 */
double get_masked_mean_array(void *data, int n, unsigned char *mask,
                             int datatype);
/**
 * \brief Standard deviation of an array of any NIfTI datatype, optionally
 *               within a mask.
 *
 * \param data     (in)  array of n values of type datatype
 * \param n        (in)  number of values
 * \param mask     (in)  n mask values; only entries > 0 count (NULL: all)
 * \param datatype (in)  NIfTI datatype code of data (DT_FLOAT32, ...)
 * \return sample standard deviation of the included values
 */
double get_masked_std_array(void *data, int n, unsigned char *mask,
                            int datatype);
/**
 * \brief Calculate percentile-based thresholds.
 *
 * \param data          (in)  array of n values
 * \param n             (in)  number of values
 * \param threshold     (out) the two thresholds, in the order of prctile
 * \param prctile       (in)  the two percentiles, in 0..100
 * \param exclude_zeros (in)  if non-zero, zeros are ignored
 */
void get_prctile_double(double *data, int n, double threshold[2],
                        double prctile[2], int exclude_zeros);

/**
 * \brief Percentile thresholds of an array of any NIfTI datatype.
 *
 * \param data          (in)  array of n values of type datatype
 * \param n             (in)  number of values
 * \param threshold     (out) the two thresholds, in the order of prctile
 * \param prctile       (in)  the two percentiles, in 0..100
 * \param exclude_zeros (in)  if non-zero, zeros are ignored
 * \param datatype      (in)  NIfTI datatype code of data (DT_FLOAT32, ...)
 */
void get_prctile(void *data, int n, double threshold[2], double prctile[2],
                 int exclude_zeros, int datatype);
/**
 * \brief Pearson correlation coefficient of two arrays of any NIfTI datatype.
 *
 * \param x             (in)  array of n values of type datatype
 * \param y             (in)  array of n values of type datatype
 * \param n             (in)  number of values
 * \param exclude_zeros (in)  if non-zero, pairs with a zero in x or y are ignored
 * \param datatype      (in)  NIfTI datatype code of x and y (DT_FLOAT32, ...)
 * \return correlation coefficient
 */
double get_corrcoef(void *x, void *y, int n, int exclude_zeros, int datatype);
/**
 * \brief Clip an array of any NIfTI datatype to [lower_limit, upper_limit].
 *
 * \param data        (in/out) array of n values of type datatype, clipped in-place
 * \param n           (in)     number of values
 * \param lower_limit (in)     values below are set to it
 * \param upper_limit (in)     values above are set to it
 * \param datatype    (in)     NIfTI datatype code of data (DT_FLOAT32, ...)
 */
void clip_data(void *data, int n, double lower_limit, double upper_limit,
               int datatype);

#endif
