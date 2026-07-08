#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <limits.h>
#include <string.h>

/* Comparison function for qsort */
static int compare_doubles(const void *a, const void *b)
{
    double diff = *(const double *)a - *(const double *)b;
    return (diff < 0) ? -1 : (diff > 0) ? 1 : 0;
}

double get_median_double(double *arr, int n, int exclude_zeros)
{
    int i, filtered_count = 0;
    double median;

    if (n <= 0) return NAN;  // Handle empty array

    double *copy = malloc(n * sizeof(double));
    if (copy == NULL) {
        fprintf(stderr, "Memory allocation failed");
        exit(EXIT_FAILURE);
    }

    for (i = 0; i < n; i++) {
        if ((exclude_zeros && arr[i] == 0.0) || isnan(arr[i]) || !isfinite(arr[i]))
            continue;
        copy[filtered_count++] = arr[i];
    }

    if (filtered_count == 0) {
        fprintf(stderr, "Error: No valid data points after filtering.\n");
        free(copy);
        exit(EXIT_FAILURE);
    }

    qsort(copy, filtered_count, sizeof(double), compare_doubles);

    if (filtered_count % 2 == 0) {
        median = (copy[filtered_count / 2 - 1] + copy[filtered_count / 2]) / 2.0;
    } else {
        median = copy[filtered_count / 2];
    }

    free(copy);
    return median;
}

double get_sum_double(double *arr, int n, int exclude_zeros)
{
    int i;
    double sum = 0.0;

    if (n <= 0) return NAN;  // Handle empty array

    for (i = 0; i < n; i++) {
        if ((exclude_zeros && arr[i] == 0.0) || isnan(arr[i]) || !isfinite(arr[i]))
            continue; // Skip invalid values and optionally zeros

        sum += arr[i];
    }

    return sum;
}

double get_mean_double(double *arr, int n, int exclude_zeros)
{
    int i, n0 = 0;
    double sum = 0.0;

    if (n <= 0) return NAN;  // Handle empty array

    for (i = 0; i < n; i++) {
        if ((exclude_zeros && arr[i] == 0.0) || isnan(arr[i]) || !isfinite(arr[i]))
            continue; // Skip zero values if exclusion is enabled
        sum += arr[i];
        n0++;
    }

    if (n0 == 0) return NAN;  // Prevent division by zero

    return sum / (double)n0;
}

double get_std_double(double *arr, int n, int exclude_zeros)
{
    int i, n0 = 0;
    double mean, variance = 0.0;

    if (n <= 0) return NAN;  // Handle empty array

    mean = get_mean_double(arr, n, exclude_zeros);

    if (isnan(mean)) return NAN;  // No valid elements

    for (i = 0; i < n; i++) {
        if ((exclude_zeros && arr[i] == 0.0) || isnan(arr[i]) || !isfinite(arr[i]))
            continue; // Skip zero values if exclusion is enabled

        variance += pow(arr[i] - mean, 2);
        n0++;
    }

    if (n0 <= 1) return NAN;  // Prevent division by zero and invalid sqrt()

    variance /= (double)(n0 - 1);

    return sqrt(variance);
}

/* Copy of Lib/CAT_Math.c orthogonal_poly() for the isolated test build. */
int orthogonal_poly(const double *x, int n, int degree, double *out)
{
    int    i, j, k;
    double mean = 0.0, dot, norm;
    double *v, **q;

    if (degree < 1 || n < 2 || x == NULL || out == NULL)
        return 0;

    v = (double *) malloc(sizeof(double) * n);
    q = (double **) malloc(sizeof(double *) * (degree + 1));
    if (v == NULL || q == NULL) { free(v); free(q); return 0; }

    for (i = 0; i < n; i++) mean += x[i];
    mean /= n;

    for (j = 0; j <= degree; j++) {
        q[j] = (double *) malloc(sizeof(double) * n);
        if (q[j] == NULL) {
            for (k = 0; k < j; k++) free(q[k]);
            free(q); free(v);
            return 0;
        }
        for (i = 0; i < n; i++) v[i] = pow(x[i] - mean, (double) j);
        for (k = 0; k < j; k++) {
            dot = 0.0;
            for (i = 0; i < n; i++) dot += q[k][i] * v[i];
            for (i = 0; i < n; i++) v[i] -= dot * q[k][i];
        }
        norm = 0.0;
        for (i = 0; i < n; i++) norm += v[i] * v[i];
        norm = sqrt(norm);
        if (norm < 1e-10) {
            for (k = 0; k <= j; k++) free(q[k]);
            free(q); free(v);
            return 0;
        }
        for (i = 0; i < n; i++) q[j][i] = v[i] / norm;
    }

    for (j = 1; j <= degree; j++)
        for (i = 0; i < n; i++)
            out[(j-1)*n + i] = q[j][i];

    for (j = 0; j <= degree; j++) free(q[j]);
    free(q); free(v);
    return 1;
}
