/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 *
 */

#include <float.h>
#include <math.h>

#include "CAT_Vol2SurfUtils.h"

double
CAT_Vol2SurfEvaluateFunction(const double *val_array, int n_val, int map_func,
                             const double *kernel, int *index_out)
{
    int i;
    double result = 0.0;

    if (index_out)
        *index_out = 0;

    if (!val_array || n_val <= 0)
        return 0.0;

    switch (map_func)
    {
    case F_MEAN:
        result = get_mean_double((double *)val_array, n_val, 0);
        break;
    case F_MEDIAN:
        result = get_median_double((double *)val_array, n_val, 0);
        break;
    case F_WAVERAGE:
        result = 0.0;
        if (kernel)
        {
            for (i = 0; i < n_val; i++)
                result += val_array[i] * kernel[i];
        }
        break;
    case F_MAXABS:
        result = 0.0;
        for (i = 0; i < n_val; i++)
        {
            if (fabs(val_array[i]) > fabs(result))
            {
                result = val_array[i];
                if (index_out)
                    *index_out = i;
            }
        }
        break;
    case F_MAX:
        result = -FLT_MAX;
        for (i = 0; i < n_val; i++)
        {
            if (val_array[i] > result)
            {
                result = val_array[i];
                if (index_out)
                    *index_out = i;
            }
        }
        break;
    case F_MIN:
        result = FLT_MAX;
        for (i = 0; i < n_val; i++)
        {
            if (val_array[i] < result)
            {
                result = val_array[i];
                if (index_out)
                    *index_out = i;
            }
        }
        break;
    case F_EXP:
        result = 0.0;
        if (kernel)
        {
            for (i = 0; i < n_val; i++)
                result += val_array[i] * kernel[i];
        }
        break;
    case F_SUM:
        result = get_sum_double((double *)val_array, n_val, 0);
        break;
    default:
        result = 0.0;
        break;
    }

    return result;
}

void
CAT_Vol2SurfBuildExpKernel(const double *length_array, int n, double exp_half,
                           double *kernel_out)
{
    int j;
    double sum = 0.0;
    const double log05 = log(0.5);

    if (!length_array || !kernel_out || n <= 0)
        return;

    if (exp_half == 0.0)
    {
        kernel_out[0] = 1.0;
        for (j = 1; j < n; j++)
            kernel_out[j] = 0.0;
        return;
    }

    for (j = 0; j < n; j++)
    {
        kernel_out[j] = exp(log05 / exp_half * length_array[j]);
        sum += kernel_out[j];
    }

    if (sum != 0.0)
    {
        for (j = 0; j < n; j++)
            kernel_out[j] /= sum;
    }
}

void
CAT_Vol2SurfBuildGaussianKernel50(int grid_steps, int grid_steps1,
                                  double *kernel_out)
{
    int i;
    double sum = 0.0;
    const double log05 = log(0.5);
    const double pi2 = 6.28319;
    double sigma;

    if (!kernel_out || grid_steps <= 1 || grid_steps1 <= 0)
        return;

    sigma = -1.0 / (2.0 * log05);

    for (i = 0; i < grid_steps1; i++)
    {
        double x = ((2.0 * (double)i) / ((double)grid_steps - 1.0)) -
                   ((double)grid_steps1 - 1.0) / ((double)grid_steps - 1.0);
        kernel_out[i] = (1.0 / sqrt(pi2 * sigma)) * exp(-(x * x) / (2.0 * sigma));
        sum += kernel_out[i];
    }

    if (sum != 0.0)
    {
        for (i = 0; i < grid_steps1; i++)
            kernel_out[i] /= sum;
    }
}
