/* Christian Gaser - christian.gaseruni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 *
 */

#include <bicpl.h>
#include "CAT_Math.h"

/**
 * produces a matrix Ainv of the same dimensions as A', so that
 * A*Ainv*A = Ainv and A*Ainv and Ainv*A are Hermitian. The computation is
 * based and SVD(A) and any singular values less than a tolerance of 1e-10
 * are treated as zero. The rank of the matrix A is returned.
 */
int
pinv(int m, int n, double **A, double **Ainv)
{
    int i, j, k, r;
    double **U, **V, **S, *W, **Ut;
  
    ALLOC2D(U, m, n);
    ALLOC2D(V, n, n);
    ALLOC2D(S,  n, n);
    ALLOC2D(Ut, n, m);
    ALLOC(W, n);

    /*
     * copy matrix A to U, because this matrix will be overwritten by
     * singular_value_decomposition
     */
    for (i = 0; i < m; i++)
        for (j = 0; j < n; j++)
            U[i][j] = A[i][j];

    (void) singular_value_decomposition(m, n, U, W, V);
  
    /*
     * rank of matrix using a somewhat arbitrary value of 1e-10 as
     * tolerance for singular values
     */
    r = 0;
    for (i = 0; i < n; i++)
        if (W[i] > TOLSVD)
          r += 1;

    if (r == 0) {
        for (i = 0; i < n; i++)
          for (j = 0; j < m; j++)
              Ainv[i][j] = 0.0;
    } else {
        for (i = 0; i < r; i++) {
            for (j = 0; j < r; j++) {
                if (i == j) S[i][j] = 1/(W[i] + EPS);
                else S[i][j] = 0.0;
            }
        }

        transpose(m, n, U, Ut);
  
        matrix_multiply(n, n, n, V, S, S);
        matrix_multiply(n, n, m, S, Ut, Ainv);
    }

    FREE(W);
    FREE2D(U);
    FREE2D(S);
    FREE2D(V);
    FREE2D(Ut);
  
    return(r);
}

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

/* Comparison function for qsort */
int compare_doubles(const void *a, const void *b) {
    double diff = *(const double *)a - *(const double *)b;
    return (diff < 0) ? -1 : (diff > 0) ? 1 : 0;
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
 */
double get_median_double(double *arr, int n, int exclude_zeros) {

    int i, filtered_count = 0;
    double median;

    if (n <= 0) return NAN;  // Handle empty array
        
    // Allocate memory for a copy of the data
    double *copy = malloc(n * sizeof(double));
    if (copy == NULL) {
        fprintf(stderr, "Memory allocation failed");
        exit(EXIT_FAILURE);
    }

    // Copy data while optionally excluding zeros
    for (i = 0; i < n; i++) {
        if ((exclude_zeros && arr[i] == 0.0) || isnan(arr[i]) || !isfinite(arr[i]))
            continue;
        copy[filtered_count++] = arr[i];
    }

    // Check if we have valid data points left
    if (filtered_count == 0) {
        fprintf(stderr, "Error: No valid data points after filtering.\n");
        free(copy);
        exit(EXIT_FAILURE);
    }

    // Sort the filtered data
    qsort(copy, filtered_count, sizeof(double), compare_doubles);

    // Compute the median
    if (filtered_count % 2 == 0) {
        // Even number of elements: average of the two middle elements
        median = (copy[filtered_count / 2 - 1] + copy[filtered_count / 2]) / 2.0;
    } else {
        // Odd number of elements: middle element
        median = copy[filtered_count / 2];
    }

    free(copy);
    return median;
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

    if (n <= 0) return NAN;  // Handle empty array
    
    for (i = 0; i < n; i++) {
        if ((exclude_zeros && arr[i] == 0.0) || isnan(arr[i]) || !isfinite(arr[i]))
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
    int i, n0 = 0;
    double mean, variance = 0.0;

    if (n <= 0) return NAN;  // Handle empty array

    mean = get_mean_double(arr, n, exclude_zeros);

    if (isnan(mean)) return NAN;  // No valid elements

    /* Calculate variance */
    for (i = 0; i < n; i++) {
        if ((exclude_zeros && arr[i] == 0.0) || isnan(arr[i]) || !isfinite(arr[i]))
            continue; // Skip zero values if exclusion is enabled
            
        variance += pow(arr[i] - mean, 2);
        n0++;
    }

    if (n0 <= 1) return NAN;  // Prevent division by zero and invalid sqrt()

    variance /= (double)(n0 - 1);

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
    double result = DBL_MAX;

    if (n <= 0) return NAN;  // Handle empty array
    
    for (i = 0; i < n; i++) {
        if ((exclude_zeros && arr[i] == 0.0) || isnan(arr[i]) || !isfinite(arr[i]))
            continue; // Skip zero values if exclusion is enabled
            
        if (arr[i] < result)
            result = arr[i];
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
    double result = -DBL_MAX;

    if (n <= 0) return NAN;  // Handle empty array
    
    for (i = 0; i < n; i++) {
        if ((exclude_zeros && arr[i] == 0.0) || isnan(arr[i]) || !isfinite(arr[i]))
            continue; // Skip zero values if exclusion is enabled
            
        if (arr[i] > result)
            result = arr[i];
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

    if (n <= 0) return NAN;  // Handle empty array
    
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

    if (n <= 0) return NAN;  // Handle empty array

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
 *  - n_vol: Number of data points.
 *  - threshold: Array where the calculated threshold values will be stored.
 *  - prctile: Array containing two percentile values for which thresholds are calculated.
 *  - exclude_zeros: Flag to indicate whether zeros should be excluded from calculations.
 *
 * Notes:
 *  - The function uses a histogram-based approach to calculate the thresholds.
 */

// Function to estimate percentiles for given thresholds (in percent, e.g., 1 and 99)
void get_prctile_double(double *data, int n, double threshold[2], 
                           double prctile[2], int exclude_zeros) {
    int i, filtered_count = 0;
    
    // Create a copy of the data for sorting
    double *copy = malloc(n * sizeof(double));
    if(copy == NULL) {
        perror("Failed to allocate memory");
        exit(EXIT_FAILURE);
    }
    
    // Copy data, optionally excluding zeros
    for (i = 0; i < n; i++) {
        if (!exclude_zeros || data[i] != 0.0) {
            copy[filtered_count++] = data[i];
        }
    }

    // Check if we have enough valid data points
    if (filtered_count == 0) {
        fprintf(stderr, "No valid data points after filtering.\n");
        free(copy);
        exit(EXIT_FAILURE);
    }
            
    // Sort the copy
    qsort(copy, filtered_count, sizeof(double), compare_doubles);
    
    // Calculate indices using the formula: index = round((n - 1) * (P/100))
    int lower_index = (int)round((filtered_count - 1) * prctile[0] / 100.0);
    int upper_index = (int)round((filtered_count - 1) * prctile[1] / 100.0);
  
    threshold[0] = copy[lower_index];
    threshold[1] = copy[upper_index];
        
    free(copy);
}

/**
 * normalize_double - Subtract mean from an array of doubles.
 *
 * This function iterates through an array of doubles to obtain a mean of 0.
 *
 * Parameters:
 *  - arr: Array of doubles.
 *  - n: Number of elements in the array.
 *
 * Returns:
 *  Overwrites the array by the mean-corrected array.
 */
void normalize_double(double *arr, int n) {
    int i;

    double mn = get_mean_double(arr, n, 0);

    for (i = 0; i < n; i++) arr[i] -= mn;

}

/**
 * Generic functions
 *
 * These function additionally provide conversion between different data types
 * and call the respective function for data type 'double'.
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