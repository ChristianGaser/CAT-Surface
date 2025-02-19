/* Christian Gaser - christian.gaseruni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 *
 */

#include "CAT_Math.h"

/**
 * convert_input_type - Converts various data types to a floating point array.
 *
 * This function is designed to convert a data array of various types into an array
 * of floats. This is useful for standardizing data input types for functions that 
 * are specifically defined to work with floating point data, especially in contexts 
 * like image processing where data might come in various formats.
 *
 * data: Pointer to the input data array. The actual data type of this array is
 *        determined by the 'datatype' parameter.
 *
 * buffer: Pointer to the output float array where the converted data will be stored.
 *          This array should be pre-allocated with enough space to hold 'n' elements.
 *          The function fills this array with the converted float values.
 *
 * n: Integer representing the number of elements in the input data array.
 *
 * datatype: Integer that specifies the type of data in the input array. This parameter
 *            uses predefined constants (e.g., DT_INT8, DT_UINT8, etc.) to represent
 *            different data types like char, unsigned char, short, unsigned short, etc.
 *
 * The function iterates over the input array, converting each element to a float based
 * on the specified datatype, and stores the result in the output float array. This 
 * facilitates the use of functions that require floating point input by providing a 
 * uniform data format.
 *
 */
void convert_input_type(void *data, double *buffer, int n, int datatype)
{
    int i;
    double tmp;
    
    /* check success of memory allocation */
    if (!buffer) {
        printf("Memory allocation error\n");
        exit(EXIT_FAILURE);
    }
   
    for (i = 0; i < n; i++) {
        switch (datatype) {
        case DT_INT8:
            tmp = (double) ((char *)data)[i];
            break;
        case DT_UINT8:
            tmp = (double) ((unsigned char *)data)[i];
            break;
        case DT_INT16:
            tmp = (double) ((short *)data)[i];
            break;
        case DT_UINT16:
            tmp = (double) ((unsigned short *)data)[i];
            break;
        case DT_INT32:
            tmp = (double) ((int *)data)[i];
            break;
        case DT_UINT32:
            tmp = (double) ((unsigned int *)data)[i];
            break;
        case DT_FLOAT32:
            tmp = (double) ((float *)data)[i];
            break;
        case DT_FLOAT64:
            tmp = ((double *)data)[i];
            break;
        default:
            fprintf(stderr, "Data type %d not handled\n", datatype);
            break;
        }
        buffer[i] = tmp;
    }
}

void convert_input_type_float(void *data, float *buffer, int n, int datatype)
{
    int i;
    float tmp;
    
    /* check success of memory allocation */
    if (!buffer) {
        printf("Memory allocation error\n");
        exit(EXIT_FAILURE);
    }
   
    for (i = 0; i < n; i++) {
        switch (datatype) {
        case DT_INT8:
            tmp = (float) ((char *)data)[i];
            break;
        case DT_UINT8:
            tmp = (float) ((unsigned char *)data)[i];
            break;
        case DT_INT16:
            tmp = (float) ((short *)data)[i];
            break;
        case DT_UINT16:
            tmp = (float) ((unsigned short *)data)[i];
            break;
        case DT_INT32:
            tmp = (float) ((int *)data)[i];
            break;
        case DT_UINT32:
            tmp = (float) ((unsigned int *)data)[i];
            break;
        case DT_FLOAT32:
            tmp = (float) ((float *)data)[i];
            break;
        case DT_FLOAT64:
            tmp = (float) ((double *)data)[i];
            break;
        default:
            fprintf(stderr, "Data type %d not handled\n", datatype);
            break;
        }
        buffer[i] = tmp;
    }
}


/**
 * convert_output_type - Converts a floating point array back to various data types.
 *
 * This function reverses the operation performed by `convert_input_type`. It converts
 * an array of floats (typically after some processing) back to a specified data type.
 * This is useful in contexts like image processing where data needs to be restored to
 * its original format after processing.
 *
 * data: Pointer to the output data array where the converted data will be stored.
 *        The actual data type of this array is determined by the 'datatype' parameter.
 *        This array should be pre-allocated with enough space to hold 'n' elements.
 *
 * buffer: Pointer to the input float array containing the data to be converted.
 *          This array contains 'n' elements of type float.
 *
 * n: Integer representing the number of elements in the input float array.
 *
 * datatype: Integer that specifies the desired output data type for the 'data' array.
 *            This parameter uses predefined constants (e.g., DT_INT8, DT_UINT8, etc.)
 *            to represent different data types like char, unsigned char, short, etc.
 *
 * The function iterates over the input float array, converting each float element back 
 * to the specified data type using rounding (via `roundf` function) and stores the 
 * result in the output array. This allows for the processed data to be converted back 
 * to its original or a different format as needed.
 *
 */
void convert_output_type(void *data, double *buffer, int n, int datatype)
{
    int i;

    for (i = 0; i < n; i++) {
        switch (datatype) {
        case DT_INT8:
            ((char*)data)[i] = (char) roundf(buffer[i]);
            break;
        case DT_UINT8:
            ((unsigned char*)data)[i] = (unsigned char) roundf(buffer[i]);
            break;
        case DT_INT16:
            ((short*)data)[i] = (short) roundf(buffer[i]);
            break;
        case DT_UINT16:
            ((unsigned short*)data)[i] = (unsigned short) roundf(buffer[i]);
            break;
        case DT_INT32:
            ((int*)data)[i] = (int) roundf(buffer[i]);
            break;
        case DT_UINT32:
            ((unsigned int *)data)[i] = (unsigned int ) roundf(buffer[i]);
            break;
        case DT_FLOAT32:
            ((float*)data)[i] = (float) buffer[i];
            break;
        case DT_FLOAT64:
            ((double*)data)[i] = (double) buffer[i];
            break;
        default:
            fprintf(stderr, "Data type %d not handled\n", datatype);
            break;
        }
    }
}

void convert_output_type_float(void *data, float *buffer, int n, int datatype)
{
    int i;

    for (i = 0; i < n; i++) {
        switch (datatype) {
        case DT_INT8:
            ((char*)data)[i] = (char) roundf(buffer[i]);
            break;
        case DT_UINT8:
            ((unsigned char*)data)[i] = (unsigned char) roundf(buffer[i]);
            break;
        case DT_INT16:
            ((short*)data)[i] = (short) roundf(buffer[i]);
            break;
        case DT_UINT16:
            ((unsigned short*)data)[i] = (unsigned short) roundf(buffer[i]);
            break;
        case DT_INT32:
            ((int*)data)[i] = (int) roundf(buffer[i]);
            break;
        case DT_UINT32:
            ((unsigned int *)data)[i] = (unsigned int ) roundf(buffer[i]);
            break;
        case DT_FLOAT32:
            ((float*)data)[i] = (float) buffer[i];
            break;
        case DT_FLOAT64:
            ((double*)data)[i] = (double) buffer[i];
            break;
        default:
            fprintf(stderr, "Data type %d not handled\n", datatype);
            break;
        }
    }
}

/**
 * swap - Swap the values of two double variables.
 *
 * This utility function is used in sorting algorithms to swap the values 
 * of two double variables.
 *
 * Parameters:
 *  - a: Pointer to the first double variable.
 *  - b: Pointer to the second double variable.
 */
void swap(double *a, double *b) {
    double temp = *a;
    *a = *b;
    *b = temp;
}

/**
 * quicksort - Sort an array of doubles using the Quick Sort algorithm.
 *
 * This function sorts an array of doubles in place using the Quick Sort algorithm. 
 * It's a recursive algorithm that sorts elements by partitioning the array.
 *
 * Parameters:
 *  - arr: The array of doubles to be sorted.
 *  - start: The starting index for the sorting process.
 *  - end: The ending index (exclusive) for the sorting process.
 */
void quicksort(double *arr, int start, int end) {
    if (end > start + 1) {
        double pivot = arr[start];
        int left = start + 1, right = end;
        while (left < right) {
            if (arr[left] <= pivot) {
                left++;
            } else {
                swap(&arr[left], &arr[--right]);
            }
        }
        swap(&arr[--left], &arr[start]);
        quicksort(arr, start, left);
        quicksort(arr, right, end);
    }
}

/**
 * get_median - Calculate the median of an array of doubles.
 *
 * This function finds the median value in an array of doubles. It sorts the array
 * using quicksort and then calculates the median.
 *
 * Parameters:
 *  - arr: Array of doubles.
 *  - n: Number of elements in the array.
 *  - exclude_zeros: Flag to indicate whether zeros should be excluded from calculations.
 *
 * Returns:
 *  The median value of the array.
 *
 * Note:
 *  This function modifies the original array by sorting it.
 */
double get_median_double(double *arr, int n, int exclude_zeros) {
    quicksort(arr, 0, n);

    // Calculate median
    if (n % 2 != 0) // If n is odd
        return arr[n / 2];
    else
        return (arr[(n - 1) / 2] + arr[n / 2]) / 2.0;
}

/**
 * get_sum - Calculate the sum of an array of doubles.
 *
 * This function calculates the sum of all elements in an array of doubles.
 *
 * Parameters:
 *  - arr: Array of doubles.
 *  - n: Number of elements in the array.
 *  - exclude_zeros: Flag to indicate whether zeros should be excluded from calculations.
 *
 * Returns:
 *  The sum of the array elements.
 */
double get_sum_double(double *arr, int n, int exclude_zeros) {
    int i;
    double sum = 0.0;
    
    for (i = 0; i < n; i++) {
        if ((exclude_zeros) && (arr[i] != 0.0))
            sum += arr[i];
        else
            sum += arr[i];
    }
    
    return sum;
}

/**
 * get_mean - Calculate the mean of an array of doubles.
 *
 * This function calculates the mean (average) value of an array of doubles.
 *
 * Parameters:
 *  - arr: Array of doubles.
 *  - n: Number of elements in the array.
 *  - exclude_zeros: Flag to indicate whether zeros should be excluded from calculations.
 *
 * Returns:
 *  The mean value of the array.
 */
double get_mean_double(double *arr, int n, int exclude_zeros) {
    int i, n0 = 0;
    double sum = 0.0;
    
    for (i = 0; i < n; i++) {
        if ((exclude_zeros) && (arr[i] != 0.0)) {
            sum += arr[i];
            n0++;
        } else {
            sum += arr[i];
            n0++;
        }
    }

    return sum / (double)n0;
}

/**
 * get_std - Calculate the standard deviation of an array of doubles.
 *
 * This function calculates the standard deviation of an array of doubles. It first
 * calculates the mean, then computes the variance, and finally the standard deviation.
 *
 * Parameters:
 *  - arr: Array of doubles.
 *  - n: Number of elements in the array
 *  - exclude_zeros: Flag to indicate whether zeros should be excluded from calculations.
 *
 * Returns: The standard deviation of the array.
*/
double get_std_double(double *arr, int n, int exclude_zeros) {
    int i, n0=0;
    double mean, variance = 0.0;

    mean = get_mean_double(arr,n, exclude_zeros);

    /* Calculate variance */
    for (i = 0; i < n; i++) {
        if ((exclude_zeros) && (arr[i] != 0.0)) {
            variance += pow(arr[i] - mean, 2);
            n0++;
        } else {
            variance += pow(arr[i] - mean, 2);
            n0++;
        }
    }

    variance /= (double)n0;

    /* Calculate standard deviation */
    return sqrt(variance);
}

/**
 * get_min - Find the minimum value in an array of doubles.
 *
 * This function iterates through an array of doubles to find the smallest element.
 *
 * Parameters:
 *  - arr: Array of doubles.
 *  - n: Number of elements in the array.
 *  - exclude_zeros: Flag to indicate whether zeros should be excluded from calculations.
 *
 * Returns:
 *  The minimum value in the array.
 */
 double get_min_double(double *arr, int n, int exclude_zeros) {
    int i;
    double result = FLT_MAX;
    
    for (i = 0; i < n; i++) {
        if (arr[i] < result) {
            if ((exclude_zeros) && (arr[i] != 0.0))
                result = arr[i];
            else
                result = arr[i];
        }
    }

    return result;
}

/**
 * get_max - Find the maximum value in an array of doubles.
 *
 * This function iterates through an array of doubles to find the largest element.
 *
 * Parameters:
 *  - arr: Array of doubles.
 *  - n: Number of elements in the array.
 *  - exclude_zeros: Flag to indicate whether zeros should be excluded from calculations.
 *
 * Returns:
 *  The maximum value in the array.
 */
 double get_max_double(double *arr, int n, int exclude_zeros) {
    int i;
    double result = -FLT_MAX;
    
    for (i = 0; i < n; i++) {
        if (arr[i] > result) {
            if ((exclude_zeros) && (arr[i] != 0.0))
                result = arr[i];
            else
                result = arr[i];
        }
    }

    return result;
}

/**
 * get_masked_mean_array_float - Calculates the mean of an array with an optional mask.
 *
 * This function computes the mean value of elements in a floating point array, optionally
 * considering only those elements that are flagged by a mask array. The mask allows for 
 * selective inclusion of elements in the mean calculation, which can be useful in 
 * situations where only specific parts of data (like certain regions in an image) are of 
 * interest.
 *
 * arr: Pointer to the float array whose mean is to be calculated.
 *       The array should contain 'size' elements.
 *
 * n: Integer representing the number of elements in the 'arr' array.
 *
 * mask: Pointer to an unsigned char array that serves as the mask. If the mask is not
 *        NULL, only elements of 'arr' at positions where 'mask' has a non-zero value are
 *        included in the mean calculation. If 'mask' is NULL, all elements of 'arr' are
 *        considered.
 *
 * The function iterates through the 'arr' array, summing up elements that are not NaN
 * and are flagged by the 'mask' array (if provided). It then calculates and returns the
 * mean of these elements. This is particularly useful in data processing tasks where
 * certain elements need to be excluded from calculations based on some criteria.
 *
 * Return: The function returns the mean of the selected elements as a double. If no
 *         elements are selected (e.g., due to all being NaN or masked out), the behavior
 *         is not defined (potential division by zero).
 */
double get_masked_mean_array_double(double *arr, int n, unsigned char *mask)
{
    double sum = 0.0;
    int i, count = 0;
    
    /* Calculate mean */
    for (i = 0; i < n; i++) {
        if (!isnan(arr[i]) && isfinite(arr[i]) && ((mask && mask[i] > 0) || !mask)) {
            sum += arr[i];
            count++;
        }
    }
    return sum / (double)count;
}

/**
 * get_masked_std_array_float - Calculates the standard deviation of an array with an optional mask.
 *
 * This function computes the standard deviation of elements in a floating point array, 
 * optionally considering only those elements that are flagged by a mask array. The mask 
 * allows for selective inclusion of elements in the standard deviation calculation, which 
 * can be useful in situations where variability of specific parts of data (like certain 
 * regions in an image or dataset) is of interest.
 *
 * arr: Pointer to the float array whose standard deviation is to be calculated.
 *       The array should contain 'size' elements.
 *
 * n: Integer representing the number of elements in the 'arr' array.
 *
 * mask: Pointer to an unsigned char array that serves as the mask. If the mask is not
 *        NULL, only elements of 'arr' at positions where 'mask' has a non-zero value are
 *        included in the standard deviation calculation. If 'mask' is NULL, all elements
 *        of 'arr' are considered.
 *
 * The function first calculates the mean of the selected elements in the array. It then
 * iterates through the array a second time to compute the variance of the selected 
 * elements. Finally, the standard deviation is calculated as the square root of the 
 * variance. This method is useful in statistical analysis and data processing tasks 
 * where understanding the dispersion or variability of data is important.
 *
 * Return: The function returns the standard deviation of the selected elements as a double. 
 *         If no elements are selected (e.g., due to all being NaN or masked out), the 
 *         behavior is not defined (potential division by zero).
 */
double get_masked_std_array_double(double *arr, int n, unsigned char *mask)
{
    double mean = 0.0, variance = 0.0;
    int i, count = 0;

    /* Calculate mean */
    for (i = 0; i < n; i++) {
        if (!isnan(arr[i]) && isfinite(arr[i]) && ((mask && mask[i] > 0) || !mask)) {
            mean += arr[i];
            count++;
        }
    }
    mean = mean / (double)count;

    /* Calculate variance */
    for (i = 0; i < n; i++)
        if (!isnan(arr[i]) && isfinite(arr[i]) && ((mask && mask[i] > 0) || !mask))
            variance += pow(arr[i] - mean, 2);
    variance /= (double)n;

    /* Calculate standard deviation */
    return sqrt(variance);
}

/**
 * get_prctile - Calculate percentile-based thresholds.
 *
 * This function computes two thresholds for a given data (src) based on the 
 * specified percentiles. It can optionally exclude zeros from the calculation.
 * The calculated thresholds are stored in the 'threshold' array.
 *
 * Parameters:
 *  - src: Pointer to the source.
 *  - n_vol: Number of data poits.
 *  - threshold: Array where the calculated threshold values will be stored.
 *  - prctile: Array containing two percentile values for which thresholds are calculated.
 *  - exclude_zeros: Flag to indicate whether zeros should be excluded from calculations.
 *
 * Notes:
 *  - The function uses a histogram-based approach to calculate the thresholds.
 */
void get_prctile_double(double *src, int n, double threshold[2], double prctile[2], int exclude_zeros) {
    double mn_thresh, mx_thresh;
    double min_src = FLT_MAX, max_src = -FLT_MAX;
    long *cumsum, *histo;
    int i, sz_histo = 1000;
    
    cumsum = (long *)malloc(sizeof(long) * sz_histo);
    histo  = (long *)malloc(sizeof(long) * sz_histo);
    
    /* check success of memory allocation */
    if (!cumsum || !histo) {
        printf("Memory allocation error\n");
        exit(EXIT_FAILURE);
    }

    // Find the minimum and maximum values in the source volume
    for (i = 0; i < n; i++) {
        min_src = fmin(src[i], min_src);
        max_src = fmax(src[i], max_src);
    }

    // Initialize and build the histogram
    for (i = 0; i < sz_histo; i++) histo[i] = 0;
    for (i = 0; i < n; i++) {
        if (exclude_zeros && (src[i] == 0)) continue; // Exclude zeros if specified
        int index = (int)round((double)sz_histo * (src[i] - min_src) / (max_src - min_src));
        histo[index]++;
    }

    // Build cumulative sum from histogram
    cumsum[0] = histo[0];
    for (i = 1; i < sz_histo; i++) 
        cumsum[i] = cumsum[i - 1] + histo[i];
    
    for (i = 1; i < sz_histo; i++)
        cumsum[i] = (long)round(100000.0 * (double)cumsum[i] / (double)cumsum[sz_histo - 1]);

    // Normalize cumulative sum and find the lower threshold
    for (i = 0; i < sz_histo; i++) {
        if (cumsum[i] >= (long)round(prctile[0] * 1000.0)) {
            threshold[0] = (double)i / (double)sz_histo * (max_src - min_src) + min_src;
            break;
        }
    }
    
    // Find the upper threshold
    for (i = sz_histo - 1; i >= 0; i--) {
        if (cumsum[i] <= (long)round(prctile[1] * 1000.0)) {
            threshold[1] = (double)i / (double)sz_histo * (max_src - min_src) + min_src;
            break;
        }
    }
    
    // Free allocated memory
    free(cumsum);
    free(histo);
}

/**
 * Generic functions
 *
 * These function additionally provide conversion between different data types
 * and call the respective function for data tpye 'double'.
 *
 */

double get_median(void *data, int n, int exclude_zeros, int datatype) {
    double *buffer, result;
   
    buffer = (double *)malloc(sizeof(double)*n);

    /* check success of memory allocation */
    if (!buffer) {
        printf("Memory allocation error\n");
        exit(EXIT_FAILURE);
    }
   
    convert_input_type(data, buffer, n, datatype);
    result = get_median_double(buffer, n, exclude_zeros);
    
    free(buffer);
    return(result);
}

double get_mean(void *data, int n, int exclude_zeros, int datatype) {
    double *buffer, result;
   
    buffer = (double *)malloc(sizeof(double)*n);

    /* check success of memory allocation */
    if (!buffer) {
        printf("Memory allocation error\n");
        exit(EXIT_FAILURE);
    }
   
    convert_input_type(data, buffer, n, datatype);
    result = get_mean_double(buffer, n, exclude_zeros);
    
    free(buffer);
    return(result);
}

double get_sum(void *data, int n, int exclude_zeros, int datatype) {
    double *buffer, result;
   
    buffer = (double *)malloc(sizeof(double)*n);

    /* check success of memory allocation */
    if (!buffer) {
        printf("Memory allocation error\n");
        exit(EXIT_FAILURE);
    }
   
    convert_input_type(data, buffer, n, datatype);
    result = get_sum_double(buffer, n, exclude_zeros);
    
    free(buffer);
    return(result);
}

double get_min(void *data, int n, int exclude_zeros, int datatype) {
    double *buffer, result;
   
    buffer = (double *)malloc(sizeof(double)*n);

    /* check success of memory allocation */
    if (!buffer) {
        printf("Memory allocation error\n");
        exit(EXIT_FAILURE);
    }
   
    convert_input_type(data, buffer, n, datatype);
    result = get_min_double(buffer, n, exclude_zeros);
    
    free(buffer);
    return(result);
}

double get_max(void *data, int n, int exclude_zeros, int datatype) {
    double *buffer, result;
   
    buffer = (double *)malloc(sizeof(double)*n);

    /* check success of memory allocation */
    if (!buffer) {
        printf("Memory allocation error\n");
        exit(EXIT_FAILURE);
    }
   
    convert_input_type(data, buffer, n, datatype);
    result = get_max_double(buffer, n, exclude_zeros);
    
    free(buffer);
    return(result);
}

double get_std(void *data, int n, int exclude_zeros, int datatype) {
    double *buffer, result;
   
    buffer = (double *)malloc(sizeof(double)*n);

    /* check success of memory allocation */
    if (!buffer) {
        printf("Memory allocation error\n");
        exit(EXIT_FAILURE);
    }
   
    convert_input_type(data, buffer, n, datatype);
    result = get_std_double(buffer, n, exclude_zeros);
    
    free(buffer);
    return(result);
}

double get_masked_mean_array(void *data, int n, unsigned char *mask, int datatype) {
    double *buffer, result;
   
    buffer = (double *)malloc(sizeof(double)*n);

    /* check success of memory allocation */
    if (!buffer) {
        printf("Memory allocation error\n");
        exit(EXIT_FAILURE);
    }
   
    convert_input_type(data, buffer, n, datatype);
    result = get_masked_mean_array_double(buffer, n, mask);
    
    free(buffer);
    return(result);
}

double get_masked_std_array(void *data, int n, unsigned char *mask, int datatype) {
    double *buffer, result;
   
    buffer = (double *)malloc(sizeof(double)*n);

    /* check success of memory allocation */
    if (!buffer) {
        printf("Memory allocation error\n");
        exit(EXIT_FAILURE);
    }
   
    convert_input_type(data, buffer, n, datatype);
    result = get_masked_std_array_double(buffer, n, mask);
    
    free(buffer);
    return(result);
}

void get_prctile(void *data, int n, double threshold[2], double prctile[2], int exclude_zeros, int datatype) {
    double *buffer;
   
    buffer = (double *)malloc(sizeof(double)*n);

    /* check success of memory allocation */
    if (!buffer) {
        printf("Memory allocation error\n");
        exit(EXIT_FAILURE);
    }
   
    convert_input_type(data, buffer, n, datatype);
    get_prctile_double(buffer, n, threshold, prctile, exclude_zeros);
    
    free(buffer);
}
