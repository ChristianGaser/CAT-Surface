/* Christian Gaser - christian.gaseruni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 *
 */

#include "CAT_Vol.h"

/**
 * ind2sub - Convert a linear index to 3D array coordinates.
 *
 * This function calculates the x, y, and z coordinates corresponding to a linear index 
 * in a 3D array. It's useful for converting a 1D array index to 3D coordinates in a volume.
 *
 * Parameters:
 *  - i: The linear index in the array.
 *  - x: Pointer to store the calculated x-coordinate.
 *  - y: Pointer to store the calculated y-coordinate.
 *  - z: Pointer to store the calculated z-coordinate.
 *  - sxy: Product of the dimensions in the x and y directions (sx * sy).
 *  - sx: The dimension in the x direction.
 */
void ind2sub(int i, int *x, int *y, int *z, int sxy, int sx) {
    int j = i % sxy; // Modulo to find position within a single z-plane

    *z = (int)floor((double)i / (double)sxy); // Calculate z-coordinate
    *y = (int)floor((double)j / (double)sx);  // Calculate y-coordinate
    *x = j % sx; // Calculate x-coordinate
}


/**
 * sub2ind - Convert 3D array coordinates to a linear index.
 *
 * This function calculates the linear index corresponding to the x, y, and z coordinates 
 * in a 3D array. It is useful for accessing elements in a linearly stored 3D array.
 *
 * Boundary handling is implemented to ensure the coordinates stay within the array limits.
 *
 * Parameters:
 *  - x: The x-coordinate in the array.
 *  - y: The y-coordinate in the array.
 *  - z: The z-coordinate in the array.
 *  - s: Array containing the dimensions of the 3D array.
 *
 * Returns:
 *  The linear index corresponding to the provided 3D coordinates.
 *
 * See Also:
 *  - ind2sub function for the inverse operation.
 */
int sub2ind(int x, int y, int z, int s[]) {
    // Boundary handling to ensure coordinates are within array limits
    x = (x < 0) ? 0 : (x > s[0] - 1) ? s[0] - 1 : x; 
    y = (y < 0) ? 0 : (y > s[1] - 1) ? s[1] - 1 : y; 
    z = (z < 0) ? 0 : (z > s[2] - 1) ? s[2] - 1 : z; 
  
    // Calculate and return the linear index
    return z * s[0] * s[1] + y * s[0] + x;
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
 *          This array should be pre-allocated with enough space to hold 'nvox' elements.
 *          The function fills this array with the converted float values.
 *
 * nvox: Integer representing the number of elements in the input data array.
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
void convert_input_type(void *data, float *buffer, int nvox, int datatype)
{
    int i;
    float tmp;
    
    /* check success of memory allocation */
    if (!buffer) {
        printf("Memory allocation error\n");
        exit(EXIT_FAILURE);
    }
   
    for (i = 0; i < nvox; i++) {
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
 *        This array should be pre-allocated with enough space to hold 'nvox' elements.
 *
 * buffer: Pointer to the input float array containing the data to be converted.
 *          This array contains 'nvox' elements of type float.
 *
 * nvox: Integer representing the number of elements in the input float array.
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
void convert_output_type(void *data, float *buffer, int nvox, int datatype)
{
    int i;

    for (i = 0; i < nvox; i++) {
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
void quicksort(double arr[], int start, int end) {
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
 *
 * Returns:
 *  The median value of the array.
 *
 * Note:
 *  This function modifies the original array by sorting it.
 */
double get_median(double arr[], int n) {
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
double get_sum(double arr[], int n, int exclude_zeros) {
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

float get_sum_float(float arr[], int n, int exclude_zeros) {
    int i;
    float sum = 0.0;
    
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
double get_mean(double arr[], int n, int exclude_zeros) {
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

float get_mean_float(float arr[], int n, int exclude_zeros) {
    int i, n0 = 0;
    float sum = 0.0;
    
    for (i = 0; i < n; i++) {
        if ((exclude_zeros) && (arr[i] != 0.0)) {
            sum += arr[i];
            n0++;
        } else {
            sum += arr[i];
            n0++;
        }
    }
    return sum / (float)n0;
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
double get_std(double arr[], int n, int exclude_zeros) {
    int i, n0=0;
    double mean, variance = 0.0;

    mean = get_mean(arr,n, exclude_zeros);

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
 double get_min(double arr[], int n, int exclude_zeros) {
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

float get_min_float(float arr[], int n, int exclude_zeros) {
    int i;
    float result = FLT_MAX;
    
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
 double get_max(double arr[], int n, int exclude_zeros) {
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

float get_max_float(float arr[], int n, int exclude_zeros) {
    int i;
    float result = -FLT_MAX;
    
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
 * localstat_float - Calculate local statistics for a 3D float array.
 *
 * This function calculates mean, median, min, max, and standard deviation 
 * within a defined distance from voxel center for each element in a 3D float array. It 
 * optionally uses a Euclidean (instead of block) distance to restrict the search area and 
 * an optional mask to optimize performance.
 *
 * Parameters:
 *  - input: Pointer to the input 3D float array.
 *  - dims: Array representing the dimensions of the input array.
 *  - dist: search distance from voxel center (1..10).
 *  - use_euclidean_dist: Flag to use Euclidean instead of block distance in calculations.
 *  - mask: Optional mask array to optimize calculations.
 *  - stat_func: Function selector for the type of statistic to calculate.
 *
 * Note: The function modifies the input array to store the results.
 */
void localstat_float(float *input, unsigned char mask[], int dims[3], int dist, 
                    int stat_func, int iters, int use_euclidean_dist)
{
    double *arr;
    int i, j, k, ind, ni, x, y, z, n, it;
    float *buffer;
    int nvox = dims[0] * dims[1] * dims[2], size_kernel;
    
    // Check for distance parameter
    if ((dist < 1) || (dist > 10)) {
        printf("Distance parameter should be in the range 1..10.\n");
        exit(EXIT_FAILURE);
    }
        
    size_kernel = (2*dist) + 1;
    
    // Memory allocation
    buffer = (float *)malloc(sizeof(float)*nvox);
    arr    = (double *)malloc(sizeof(double)*size_kernel*size_kernel*size_kernel);

    // Check memory allocation success
    if (!buffer || !arr) {
        printf("Memory allocation error\n");
        exit(EXIT_FAILURE);
    }

    // Initialize buffer to zero
    for (i = 0; i<nvox; i++) buffer[i] = input[i];
    
    // Main filter process
    for (it = 0; it<iters; it++) {
        for (z=0; z<dims[2]; z++) for (y=0; y<dims[1]; y++) for (x=0; x<dims[0]; x++) {
            ind = sub2ind(x,y,z,dims);
            n = 0;
            
            memset(arr,0.0,size_kernel*size_kernel*size_kernel*sizeof(double));  
            // Iterate through kernel
            for (i=-dist; i<=dist; i++) for (j=-dist; j<=dist; j++) for (k=-dist; k<=dist; k++) {
                ni = sub2ind(x+i,y+j,z+k,dims);
                
                // Check for NaNs, Infinities
                if (isnan(input[ni]) || !isfinite(input[ni]))
                    continue;
    
                // Check for optional mask
                if (mask && mask[ni] == 0) 
                    continue;
    
                // Check for Euclidean distance if required
                if (use_euclidean_dist && sqrtf((float)((i * i) + (j * j) + (k * k))) > (float)dist)
                    continue;
                
                arr[n] = (double)input[ni];
                n++;
            }
    
            // Check for NaNs, Infinities
            if (isnan(input[ind]) || !isfinite(input[ind]))
                continue;
    
            // Check for optional mask
            if (mask && mask[ind] == 0) 
                continue;
    
            // Calculate local statistics based on the selected function
            switch (stat_func) {
            case F_MEAN:
                buffer[ind] = (float)get_mean(arr, n, 0);
                break;
            case F_MEDIAN:
                buffer[ind] = (float)get_median(arr, n);
                break;
            case F_STD:
                buffer[ind] = (float)get_std(arr, n, 0);
                break;
            case F_MIN:
                buffer[ind] = (float)get_min(arr, n, 0);
                break;
            case F_MAX:
                buffer[ind] = (float)get_max(arr, n, 0);
                break;
            default:
                fprintf(stderr, "Data Function %d not handled\n", stat_func);
                break;
            }
        }
        // Copy results back to input array
        for (i = 0; i<nvox; i++) input[i] = buffer[i];
    }

    // Free allocated memory
    free(buffer);
    free(arr);
}

/**
 * wrapper to call localstat_float for any data type 
 */
void localstat3(void *data, unsigned char mask[], int dims[3], int dist, 
                    int stat_func, int iters, int use_euclidean_dist, int datatype)
{
    int nvox;
    float *buffer;
   
    // Check for distance parameter
    if ((dist < 1) || (dist > 10)) {
        printf("Distance parameter should be in the range 1..10.\n");
        exit(EXIT_FAILURE);
    }

    nvox = dims[0]*dims[1]*dims[2];
    buffer = (float *)malloc(sizeof(float)*nvox);

    /* check success of memory allocation */
    if (!buffer) {
        printf("Memory allocation error\n");
        exit(EXIT_FAILURE);
    }
   
    convert_input_type(data, buffer, nvox, datatype);
    localstat_float(buffer, mask, dims, dist, stat_func, iters, use_euclidean_dist);
    convert_output_type(data, buffer, nvox, datatype);
    
    free(buffer);
}

/**
 * wrapper to use localstat_float for median calculation for any data type 
 */
void median3(void *data, unsigned char *mask, int dims[3], int iters, int datatype)
{
    int nvox;
    float *buffer;
   
    nvox = dims[0]*dims[1]*dims[2];
    buffer = (float *)malloc(sizeof(float)*nvox);

    /* check success of memory allocation */
    if (!buffer) {
        printf("Memory allocation error\n");
        exit(EXIT_FAILURE);
    }
   
    convert_input_type(data, buffer, nvox, datatype);
    
    /* use kernel 3x3x3 */
    localstat_float(buffer, mask, dims, 1, F_MEDIAN, iters, 0);
    convert_output_type(data, buffer, nvox, datatype);
    
    free(buffer);
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
 * size: Integer representing the number of elements in the 'arr' array.
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
double get_masked_mean_array_float(float arr[], int size, unsigned char mask[])
{
    double sum = 0.0;
    int i, n = 0;
    
    /* Calculate mean */
    for (i = 0; i < size; i++) {
        if (!isnan(arr[i]) && isfinite(arr[i]) && ((mask && mask[i] > 0) || !mask)) {
            sum += arr[i];
            n++;
        }
    }
    return sum / (double)n;
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
 * size: Integer representing the number of elements in the 'arr' array.
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
double get_masked_std_array_float(float arr[], int size, unsigned char mask[])
{
    double mean = 0.0, variance = 0.0;
    int i, n = 0;

    /* Calculate mean */
    for (i = 0; i < size; i++) {
        if (!isnan(arr[i]) && isfinite(arr[i]) && ((mask && mask[i] > 0) || !mask)) {
            mean += arr[i];
            n++;
        }
    }
    mean = mean / (double)n;

    /* Calculate variance */
    for (i = 0; i < size; i++)
        if (!isnan(arr[i]) && isfinite(arr[i]) && ((mask && mask[i] > 0) || !mask))
            variance += pow(arr[i] - mean, 2);
    variance /= (double)n;

    /* Calculate standard deviation */
    return sqrt(variance);
}


/**
 * pmin - Finds the minimum positive value in an array and its index.
 *
 * This function iterates through a given array of floats to find the minimum positive
 * value and the index at which this value occurs. It is specifically designed to ignore
 * non-positive values (i.e., values less than or equal to zero).
 *
 * A: Pointer to the float array in which the minimum positive value is to be searched.
 *     The array should contain 'sA' elements.
 *
 * sA: Integer representing the size of the 'A' array.
 *
 * minimum: Pointer to a float where the minimum positive value found in the array will be stored.
 *           If no positive value is found, it will store FLT_MAX.
 *
 * index: Pointer to an integer where the index of the minimum positive value found in the array
 *         will be stored. If no positive value is found, it will store 0.
 *
 * The function initializes 'minimum' to the maximum float value (FLT_MAX) and 'index' to 0.
 * It then iterates through the array, updating 'minimum' and 'index' whenever it finds a
 * new minimum positive value. This is useful for tasks where identifying the smallest
 * positive element of an array and its position is required.
 *
 */
void pmin(float *A, int sA, float *minimum, int *index)
{
    int i; 
    
    *minimum = FLT_MAX;
    *index = 0;
    
    for (i = 0; i<sA; i++) {
        if ((A[i] > 0.0) && (*minimum > A[i])) { 
            *minimum = A[i]; 
            *index   = i;
        }
    }
}

/**
 * isoval - Calculate the linearly interpolated value in a 3D volume.
 *
 * This function reads out the linear interpolated value of a volume at a specific position
 * in 3D space. If a NIfTI image pointer is provided, it first transforms the coordinates from 
 * world to voxel space. It then performs linear interpolation using the nearest neighbors.
 *
 * Parameters:
 *  - vol: The 3D volume in which to interpolate.
 *  - x, y, z: Coordinates at which to interpolate the value.
 *  - dims: Array containing the dimensions of the volume.
 *  - nii_ptr: Pointer to a NIfTI image for coordinate transformation (optional).
 *
 * Returns:
 *  The interpolated value at the given coordinates. Returns NaN if unable to interpolate.
 *
 * Notes:
 *  - The function uses linear interpolation based on the values of the 8 nearest neighbors.
 *  - Coordinates are in C-notation. See 'ind2sub' for details on coordinate system.
 */
float isoval(float vol[], float x, float y, float z, int dims[], nifti_image *nii_ptr) {
    float seg = 0.0, n = 0.0;
    float world_coords[3] = {x, y, z}; // Define world coordinates (in mm)
    int i;
    
    // Convert from world to voxel space if NIfTI pointer is defined
    if (nii_ptr) {
        mat44 mat = nii_ptr->sto_xyz; // Transformation matrix
        mat44 inverse_mat = nifti_mat44_inverse(mat); // Inverse transformation matrix
    
        // Convert world coordinates to voxel coordinates
        x = world_coords[0] * inverse_mat.m[0][0] + world_coords[1] * inverse_mat.m[0][1] + 
            world_coords[2] * inverse_mat.m[0][2] + inverse_mat.m[0][3];
        y = world_coords[0] * inverse_mat.m[1][0] + world_coords[1] * inverse_mat.m[1][1] + 
            world_coords[2] * inverse_mat.m[1][2] + inverse_mat.m[1][3];
        z = world_coords[0] * inverse_mat.m[2][0] + world_coords[1] * inverse_mat.m[2][1] + 
            world_coords[2] * inverse_mat.m[2][2] + inverse_mat.m[2][3];
    }

    // Floor coordinates for interpolation
    float fx = floor(x), fy = floor(y), fz = floor(z);
    float cx = floor(x+1), cy = floor(y+1), cz = floor(z+1);
    
    // Weight factors for interpolation
    float wfx = cx - x, wfy = cy - y, wfz = cz - z;
    float wcx = x - fx, wcy = y - fy, wcz = z - fz;
    float N[8], W[8]; // Neighbors and their weights
        
    // Calculate value of the 8 neighbors and their distance weight
    for (i = 0; i < 8; i++) {
        int ix = (i & 1) ? cx : fx;
        int iy = (i & 2) ? cy : fy;
        int iz = (i & 4) ? cz : fz;

        N[i] = vol[sub2ind(ix, iy, iz, dims)];
        W[i] = ((i & 1) ? wcx : wfx) * ((i & 2) ? wcy : wfy) * ((i & 4) ? wcz : wfz);

        // Perform interpolation using neighbors and weights
        if (!isnan(N[i]) && isfinite(N[i])) {
            seg += N[i] * W[i];
            n += W[i];
        }
    }

    // Return interpolated value or NaN if unable to interpolate
    return n > 0.0 ? seg / n : FNAN;
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
void get_prctile(float *src, int nvox, double threshold[2], double prctile[2], int exclude_zeros) {
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
    for (i = 0; i < nvox; i++) {
        min_src = fmin((double)src[i], min_src);
        max_src = fmax((double)src[i], max_src);
    }

    // Initialize and build the histogram
    for (i = 0; i < sz_histo; i++) histo[i] = 0;
    for (i = 0; i < nvox; i++) {
        if (exclude_zeros && (src[i] == 0)) continue; // Exclude zeros if specified
        int index = (int)round((double)sz_histo * ((double)src[i] - min_src) / (max_src - min_src));
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
 * convxy_float - Apply 2D convolution to a slice of data.
 *
 * This function applies a slice-wise 2D convolution to a given input data array 
 * using specified filter kernels along x and y dimensions. The result is stored 
 * in the output array.
 *
 * Parameters:
 *  - out: Output array where the convolution result is stored.
 *  - xdim, ydim: Dimensions of the input data slice.
 *  - filtx, filty: Filter kernels for convolution along x and y dimensions.
 *  - fxdim, fydim: Dimensions of the filter kernels.
 *  - xoff, yoff: Offsets for the filter kernels.
 *  - buff: Buffer array for intermediate results.
 *
 * Notes:
 * This is a slightly modified function from spm_conv_vol.c from SPM12.
 */
static void convxy_float(float out[], int xdim, int ydim, double filtx[], double filty[], 
                         int fxdim, int fydim, int xoff, int yoff, float buff[]) {
    int x,y,k;

    for (y=0; y<ydim; y++) {
        for (x=0; x<xdim; x++) {
            buff[x] = out[x+y*xdim];
            if (!isfinite(buff[x]))
                buff[x] = 0.0;
        }
        for (x=0; x<xdim; x++) {
            double sum1 = 0.0;
            int fstart, fend;
            fstart = ((x-xoff >= xdim) ? x-xdim-xoff+1 : 0);
            fend = ((x-(xoff+fxdim) < 0) ? x-xoff+1 : fxdim);

            for (k=fstart; k<fend; k++)
                sum1 += (double)buff[x-xoff-k]*filtx[k];
            out[x+y*xdim] = (float)sum1;
        }
    }
    for (x=0; x<xdim; x++) {
        for (y=0; y<ydim; y++)
            buff[y] = out[x+y*xdim];

        for (y=0; y<ydim; y++) {
            double sum1 = 0.0;
            int fstart, fend;
            fstart = ((y-yoff >= ydim) ? y-ydim-yoff+1 : 0);
            fend = ((y-(yoff+fydim) < 0) ? y-yoff+1 : fydim);

            for (k=fstart; k<fend; k++)
                sum1 += (double)buff[y-yoff-k]*filty[k];
            out[y*xdim+x] = (float)sum1;
        }
    }
}

/**
 * convxyz_float - Apply 3D convolution to a volume.
 *
 * This function applies 3D convolution to a given volume using separate 1D 
 * filter kernels along x, y, and z dimensions. The output is stored in a 
 * separate output volume.
 *
 * Parameters:
 *  - iVol: Input volume for convolution.
 *  - filtx, filty, filtz: Filter kernels for convolution along x, y, and z dimensions.
 *  - fxdim, fydim, fzdim: Dimensions of the filter kernels.
 *  - xoff, yoff, zoff: Offsets for the filter kernels.
 *  - oVol: Output volume where the convolution result is stored.
 *  - dims: Array containing the dimensions of the input volume.
 * 
 * Returns:
 * 0 on successful completion.
 *
 * Notes:
 * The function applies a slice-wise 2D convolution using 'convxy_float'
 * and then combines these results to achieve 3D convolution. It allocates
 * temporary buffers for intermediate results and performs necessary memory
 * management.
 * This is a slightly modified function from spm_conv_vol.c from SPM12
*/
int convxyz_float(float *iVol, double filtx[], double filty[], double filtz[],
                int fxdim, int fydim, int fzdim,
                int xoff, int yoff, int zoff, float *oVol, int dims[3]) {
    float *tmp, *buff, **sortedv, *obuf;
    int xy, z, y, x, k, fstart, fend, startz, endz;
    int xdim, ydim, zdim;

    xdim = dims[0];
    ydim = dims[1];
    zdim = dims[2];

    tmp = (float *)malloc(sizeof(float)*xdim*ydim*fzdim);
    buff = (float *)malloc(sizeof(float)*((ydim>xdim) ? ydim : xdim));
    sortedv = (float **)malloc(sizeof(float *)*fzdim);

    if (!tmp || !buff || !sortedv) {
        fprintf(stderr,"Memory allocation error\n");
        exit(EXIT_FAILURE);
    }

    startz = ((fzdim+zoff-1<0) ? fzdim+zoff-1 : 0);
    endz = zdim+fzdim+zoff-1;

    for (z=startz; z<endz; z++) {
        double sum2 = 0.0;

        if (z >= 0 && z<zdim) {
            for (y=0; y<ydim; y++) for (x=0; x<xdim; x++)
                tmp[((z%fzdim)*xdim*ydim)+(y*xdim)+x] = iVol[(z*xdim*ydim)+(y*xdim)+x];  
                convxy_float(tmp+((z%fzdim)*xdim*ydim),xdim, ydim,
                        filtx, filty, fxdim, fydim, xoff, yoff, buff);
        }
        if (z-fzdim-zoff+1>=0 && z-fzdim-zoff+1<zdim) {
            fstart = ((z >= zdim) ? z-zdim+1 : 0);
            fend = ((z-fzdim < 0) ? z+1 : fzdim);

            for (k=0; k<fzdim; k++) {
                int z1 = (((z-k)%fzdim)+fzdim)%fzdim;
                sortedv[k] = &(tmp[z1*xdim*ydim]);
            }

            for (k=fstart, sum2=0.0; k<fend; k++)
                sum2 += filtz[k];

            obuf = &oVol[(z-fzdim-zoff+1)*ydim*xdim];
            if (sum2) {
                for (xy=0; xy<xdim*ydim; xy++) {
                    double sum1=0.0;
                    for (k=fstart; k<fend; k++)
                        sum1 += filtz[k]*(double)sortedv[k][xy];

                    obuf[xy] = (float)sum1/sum2;
                }
            }
            else
                for (xy=0; xy<xdim*ydim; xy++)
                    obuf[xy] = 0.0;
        }
    }
    free(tmp);
    free(buff);
    free(sortedv);
    return(0);
}

/**
 * convxyz_uint8 - Apply 3D convolution to a volume with type unsigned char.
 *
 * This function applies 3D convolution to a given volume using separate 1D 
 * filter kernels along x, y, and z dimensions. The output is stored in a 
 * separate output volume.
 *
 * Parameters:
 *  - iVol: Input volume for convolution.
 *  - filtx, filty, filtz: Filter kernels for convolution along x, y, and z dimensions.
 *  - fxdim, fydim, fzdim: Dimensions of the filter kernels.
 *  - xoff, yoff, zoff: Offsets for the filter kernels.
 *  - oVol: Output volume where the convolution result is stored.
 *  - dims: Array containing the dimensions of the input volume.
 * 
 * Returns:
 * 0 on successful completion.
 *
 * Notes:
 * The function applies a slice-wise 2D convolution using 'convxy_float'
 * and then combines these results to achieve 3D convolution. It allocates
 * temporary buffers for intermediate results and performs necessary memory
 * management.
 * This is a slightly modified function from spm_conv_vol.c from SPM12
*/
int convxyz_uint8(unsigned char *iVol, double filtx[], double filty[], double filtz[],
    int fxdim, int fydim, int fzdim, int xoff, int yoff, int zoff,
    unsigned char *oVol, int dims[3])
{
    float *tmp, *buff, **sortedv;
    int xy, z, y, x, k, fstart, fend, startz, endz;
    int xdim, ydim, zdim;
    float tmp2;
    unsigned char *obuf;

    xdim = dims[0];
    ydim = dims[1];
    zdim = dims[2];

    tmp = (float *)malloc(sizeof(float)*xdim*ydim*fzdim);
    buff = (float *)malloc(sizeof(float)*((ydim>xdim) ? ydim : xdim));
    sortedv = (float **)malloc(sizeof(float *)*fzdim);

    if (!tmp || !buff || !sortedv) {
        fprintf(stderr,"Memory allocation error\n");
        exit(EXIT_FAILURE);
    }

    startz = ((fzdim+zoff-1<0) ? fzdim+zoff-1 : 0);
    endz = zdim+fzdim+zoff-1;

    for (z=startz; z<endz; z++) {
        double sum2 = 0.0;

        if (z >= 0 && z<zdim) {
            for (y=0; y<ydim; y++) for (x=0; x<xdim; x++)
                tmp[((z%fzdim)*xdim*ydim)+(y*xdim)+x] = (float)iVol[(z*xdim*ydim)+(y*xdim)+x];  
            convxy_float(tmp+((z%fzdim)*xdim*ydim),xdim, ydim,
                filtx, filty, fxdim, fydim, xoff, yoff, buff);
        }
        if (z-fzdim-zoff+1>=0 && z-fzdim-zoff+1<zdim) {
            fstart = ((z >= zdim) ? z-zdim+1 : 0);
            fend = ((z-fzdim < 0) ? z+1 : fzdim);

            for (k=0; k<fzdim; k++) {
                int z1 = (((z-k)%fzdim)+fzdim)%fzdim;
                sortedv[k] = &(tmp[z1*xdim*ydim]);
            }

            for (k=fstart, sum2=0.0; k<fend; k++)
                sum2 += filtz[k];

            obuf = oVol;
            obuf = &obuf[(z-fzdim-zoff+1)*ydim*xdim];
            if (sum2) {
                for (xy=0; xy<xdim*ydim; xy++)
                {
                    float sum1=0.0;
                    for (k=fstart; k<fend; k++)
                        sum1 += filtz[k]*sortedv[k][xy];
                    tmp2 = sum1/sum2;
                    if (tmp2<0.0) tmp2 = 0.0;
                    else if (tmp2>255.0) tmp2 = 255.0;
                    obuf[xy] = (unsigned char)roundf(tmp2);
                }
            }
            else
                for (xy=0; xy<xdim*ydim; xy++)
                    obuf[xy] = 0;
        }
    }
    free(tmp);
    free(buff);
    free(sortedv);
    return(0);
}

/**
 * smooth_float - Apply Gaussian smoothing to a 3D volume.
 *
 * This function performs Gaussian smoothing on a 3D volume. The size of the 
 * Gaussian kernel is defined by the Full Width at Half Maximum (FWHM) parameter. 
 * Optionally, the function can perform smoothing only inside a specified mask 
 * and correct the smoothed values based on the mask.
 *
 * Parameters:
 *  - vol: Pointer to the 3D volume to be smoothed.
 *  - dims: Array containing the dimensions of the volume.
 *  - voxelsize: Array containing the size of each voxel.
 *  - fwhm: Array containing the FWHM values for each dimension.
 *  - use_mask: Flag to indicate whether masked smoothing should be used.
 *
 * Notes:
 *  - The function calculates the Gaussian kernel based on FWHM and voxel size.
 *  - Masked smoothing excludes zero-valued voxels and adjusts the smoothing accordingly.
 */
void smooth_float(float *vol, int dims[3], double voxelsize[3], double fwhm[3], int use_mask) {
    int i;
    double xsum, ysum, zsum;
    double *x, *y, *z, s[3];
    int xyz[3], nvox, sum_mask;
    float *mask;
    unsigned char *mask2;

    nvox = dims[0] * dims[1] * dims[2];

    // Calculate the standard deviation and kernel size based on FWHM and voxel size
    for (i = 0; i < 3; i++) {
        s[i] = fwhm[i] / voxelsize[i];
        if (s[i] < 1.0) s[i] = 1.0;
        s[i] /= sqrt(8.0 * log(2.0));
        xyz[i] = (int)round(6.0 * s[i]);
    }

    // Memory allocation for Gaussian kernels
    x = (double *)malloc(sizeof(double) * ((2 * xyz[0]) + 1));
    y = (double *)malloc(sizeof(double) * ((2 * xyz[1]) + 1));
    z = (double *)malloc(sizeof(double) * ((2 * xyz[2]) + 1));

    if (!x || !y || !z) {
        fprintf(stderr,"Memory allocation error\n");
        exit(EXIT_FAILURE);
    }

    // Initialize and build the mask for masked smoothing
    if (use_mask) {
        mask  = (float *) malloc(sizeof(float)*nvox);
        mask2 = (unsigned char *) malloc(sizeof(unsigned char)*nvox);

        if (!mask || !mask2) {
            fprintf(stderr,"Memory allocation error\n");
            exit(EXIT_FAILURE);
        }
        sum_mask = 0;
        // Build the mask based on the volume data
        for (i = 0; i < nvox; i++) {
            if (vol[i] == 0.0) {
                mask[i]  = 0.0;
                mask2[i] = 0;
            } else {
            mask[i]  = 1.0;
                mask2[i] = 1;
                sum_mask++;
            }
        }
    }

    // Build the Gaussian kernel
    for (i = -xyz[0]; i <= xyz[0]; i++) x[i + xyz[0]] = exp(-pow((double)i, 2) / (2.0 * pow(s[0], 2)));
    for (i = -xyz[1]; i <= xyz[1]; i++) y[i + xyz[1]] = exp(-pow((double)i, 2) / (2.0 * pow(s[1], 2)));
    for (i = -xyz[2]; i <= xyz[2]; i++) z[i + xyz[2]] = exp(-pow((double)i, 2) / (2.0 * pow(s[2], 2)));
    
    // Normalize the Gaussian kernel
    xsum = ysum = zsum = 0.0;
    for (i = 0; i < ((2 * xyz[0]) + 1); i++) xsum += x[i];
    for (i = 0; i < ((2 * xyz[1]) + 1); i++) ysum += y[i];
    for (i = 0; i < ((2 * xyz[2]) + 1); i++) zsum += z[i];
    for (i = 0; i < ((2 * xyz[0]) + 1); i++) x[i] /= xsum;
    for (i = 0; i < ((2 * xyz[1]) + 1); i++) y[i] /= ysum;
    for (i = 0; i < ((2 * xyz[2]) + 1); i++) z[i] /= zsum;
    
    // Apply convolution with the Gaussian kernel
    convxyz_float(vol, x, y, z, ((2 * xyz[0]) + 1), ((2 * xyz[1]) + 1), ((2 * xyz[2]) + 1), -xyz[0], -xyz[1], -xyz[2], vol, dims);
    
    // Apply masked smoothing if selected
    if (use_mask) {
        if (sum_mask > 0) {
            convxyz_float(mask, x, y, z, ((2 * xyz[0]) + 1), ((2 * xyz[1]) + 1), ((2 * xyz[2]) + 1), -xyz[0], -xyz[1], -xyz[2], mask, dims);
            for (i = 0; i < nvox; i++) {
                vol[i] = mask2[i] > 0 ? vol[i] / mask[i] : 0.0;
            }
        }
        free(mask);
        free(mask2);
    }
    
    // Free allocated memory for the Gaussian kernels
    free(x);
    free(y);
    free(z);
}

/**
 * wrapper to call smooth_float for any data type 
 */
void smooth3(void *data, int dims[3], double voxelsize[3], double fwhm[3], int use_mask, int datatype)
{
    int nvox;
    float *buffer;
   
    nvox = dims[0]*dims[1]*dims[2];
    buffer = (float *)malloc(sizeof(float)*nvox);

    /* check success of memory allocation */
    if (!buffer) {
        printf("Memory allocation error\n");
        exit(EXIT_FAILURE);
    }
   
    convert_input_type(data, buffer, nvox, datatype);
    smooth_float(buffer, dims, voxelsize, fwhm, use_mask);
    convert_output_type(data, buffer, nvox, datatype);
    
    free(buffer);
}

/**
 * pmax - Calculate a conditional maximum value from a set of voxels.
 *
 * This function is used in the projection_based_thickness process to find the maximum value among
 * the voxels that are in the range of White Matter Distance (WMD), considering certain constraints.
 * It first finds the pure maximum based on several criteria and then calculates the mean of the 
 * highest values under the same constraints.
 *
 * Parameters:
 *  - GMT: Array of thickness/WMD values of neighbors.
 *  - PPM: Array of projection values.
 *  - SEG: Array of segmentation values.
 *  - ND: Array of Euclidean distances.
 *  - WMD: White Matter Distance for the current voxel.
 *  - SEGI: Segmentation value of the current voxel.
 *  - sA: Size of the arrays (number of elements to consider).
 *
 * Returns:
 *  The calculated maximum value under the specified conditions.
 *
 * Notes:
 *  The function applies several constraints based on segmentation and distance measures to determine
 *  the relevant maximum value. This includes checking the range of projection, upper and lower distance 
 *  boundaries, and segmentation-based conditions.
 */
float pmax(const float GMT[], const float PPM[], const float SEG[], const float ND[], const float WMD, const float SEGI, const int sA) {
    float maximum = WMD;
    int i;

    // Calculate the pure maximum under specified conditions
    for (i = 0; i <= sA; i++) {
        if ((GMT[i] < FLT_MAX) && (maximum < GMT[i]) &&              /* thickness/WMD of neighbors should be larger */
                (SEG[i] >= 1.0) && (SEGI > 1.2 && SEGI <= 2.75) &&       /* projection range */
                (((PPM[i] - ND[i] * 1.2) <= WMD)) &&                 /* upper boundary - maximum distance */
                (((PPM[i] - ND[i] * 0.5) >  WMD) || (SEG[i] < 1.5)) && /* lower boundary - minimum distance - corrected values outside */
                ((((SEGI * MAX(1.0,MIN(1.2, SEGI-1.5))) >= SEG[i])) || (SEG[i] < 1.5))) { /* for high values will project data over sulcal gaps */
            maximum = GMT[i];
        }
    }

    // Calculate the mean of the highest values under the same conditions
    float maximum2 = maximum, m2n = 0.0; 
    for (i = 0; i <= sA; i++) {
        if ((GMT[i] < FLT_MAX) && ((maximum - 1) < GMT[i]) && 
                 (SEG[i] >= 1.0) && (SEGI > 1.2 && SEGI <= 2.75) && 
                 (((PPM[i] - ND[i] * 1.2) <= WMD)) && 
                 (((PPM[i] - ND[i] * 0.5) >  WMD) || (SEG[i] < 1.5)) &&
                 ((((SEGI * MAX(1.0,MIN(1.2, SEGI-1.5))) >= SEG[i])) || (SEG[i] < 1.5))) {
            maximum2 += GMT[i]; 
            m2n++;
        } 
    }
    if (m2n > 0.0)
        maximum = (maximum2 - maximum) / m2n;

    return maximum;
}

/**
 * projection_based_thickness - Calculate the thickness of segmented structures.
 *
 * This function estimates the projection-based thickness of segmented structures 
 * in a 3D volume, typically used in medical imaging for brain tissue analysis. 
 * It requires PVE label images and distance maps for WM and CSF.
 *
 * Parameters:
 *  - SEG: PVE label image with labels for CSF, GM, and WM.
 *  - WMD: White Matter distance map.
 *  - CSFD: Cerebrospinal Fluid distance map.
 *  - GMT: Output thickness image.
 *  - dims: Array containing the dimensions of the volume.
 *  - voxelsize: Array containing the size of each voxel.
 *
 * Note:
 *  - The function assumes specific labels for CSF (1), GM (2), and WM (3).
 */
void projection_based_thickness(float *SEG, float *WMD, float *CSFD, float *GMT, int dims[3], double *voxelsize) 
{
    // Initialization and pre-processing
    const int nvox = dims[0] * dims[1] * dims[2];
    const int x = dims[0], y = dims[1], xy = x * y;
    const float s2 = sqrt(2.0), s3 = sqrt(3.0);
    const int   NI[] = {0, -1, -x+1, -x, -x-1, -xy+1, -xy, -xy-1, -xy+x+1, -xy+x, -xy+x-1, -xy-x+1, -xy-x, -xy-x-1}; // Neighbor index offsets
    const float ND[] = {0.0, 1.0, s2, 1.0, s2, s2, 1.0, s2, s3, s2, s3, s3, s2, s3}; // Neighbor distances
    const int sN = sizeof(NI) / sizeof(NI[0]); // Number of neighbors

    // Variables for processing
    float DN[sN], DI[sN], GMTN[sN], WMDN[sN], SEGN[sN], DNm;
    float GMTi, CSFDi;
    int i, n, ni, u, v, w, nu, nv, nw, count_WM = 0, count_CSF = 0;

    /* GMT should be allocated */
    if (!GMT) {
        GMT = (float *)malloc(sizeof(float)*nvox);
        if (!GMT) {
            fprintf(stderr,"Memory allocation error\n");
            exit(EXIT_FAILURE);
        }
    }    

    // Initial distance checks and assignment
    for (i = 0; i < nvox; i++) {
        // Initial GMT value and WM/CSF counts
        GMT[i] = WMD[i] + 0.0;
        
        // proof distance input
        if (SEG[i] >= GWM) count_WM++;
        if (SEG[i] <= CGM) count_CSF++;
    }

    // Error checks for WM and CSF voxels
    if (count_WM == 0) {
        fprintf(stderr,"ERROR: no WM voxels\n");
        exit(EXIT_FAILURE);
    }    
    if (count_CSF == 0) {
        fprintf(stderr,"ERROR: no CSF voxels\n");
        exit(EXIT_FAILURE);
    }    

    // Forward thickness calculation
    for (i = 0; i < nvox; i++) {
        // Process only GM voxels
        if (SEG[i] > CSF && SEG[i] < WM) {
            // Neighborhood processing
            ind2sub(i, &u, &v, &w, xy, x);

            /* read neighbour values */
            for (n = 0; n < sN; n++) {
                ni = i + NI[n];
                ind2sub(ni, &nu, &nv, &nw, xy, x);

                if ((ni < 0) || (ni >= nvox) || (abs(nu - u) > 1) || (abs(nv - v) > 1) || (abs(nw - w) > 1))
                    ni = i;

                GMTN[n] = GMT[ni];
                WMDN[n] = WMD[ni];
                SEGN[n] = SEG[ni];
            }

            /* find minimum distance within the neighborhood */
            DNm = pmax(GMTN, WMDN, SEGN, ND, WMD[i], SEG[i], sN);
            GMT[i] = DNm;
        }
    }

    // Backward search for thickness correction
    for (i = nvox - 1; i >= 0; i--) {
        // Process only GM voxels
        if (SEG[i] > CSF && SEG[i] < WM) {
            ind2sub(i, &u, &v, &w, xy, x);

            /* read neighbour values */
            for (n = 0; n < sN; n++) {
                ni = i - NI[n];
                ind2sub(ni, &nu, &nv, &nw, xy, x);

                if ((ni < 0) || (ni >= nvox) || (abs(nu - u) > 1) || (abs(nv - v) > 1) || (abs(nw - w) > 1))
                    ni = i;

                GMTN[n] = GMT[ni];
                WMDN[n] = WMD[ni];
                SEGN[n] = SEG[ni];
            }

            /* find minimum distance within the neighborhood */
            DNm = pmax(GMTN, WMDN, SEGN, ND, WMD[i], SEG[i], sN);
            if ((GMT[i] < DNm) && (DNm > 0)) GMT[i] = DNm;
        }
    }

    // Post-processing to refine GMT values
    for (i = 0; i < nvox; i++) {
        if (SEG[i] < CGM || SEG[i] > GWM)
            GMT[i] = 0.0;
    }

    // Final GMT adjustment based on CSFD and WMD
    for (i = 0; i < nvox; i++) {
        if (SEG[i] >= CGM && SEG[i] <= GWM) {
            GMTi = CSFD[i] + WMD[i];
            CSFDi = GMT[i] - WMD[i];

            /*if ( CSFD[i]>CSFDi ) CSFD[i] = CSFDi;          
            else GMT[i]  = GMTi;*/
        }
    }    
}

/**
 * vbdist - Calculate the voxel-wise Euclidean distance to an object in a 3D volume.
 *
 * This function computes the Euclidean distance from each voxel within a given mask 
 * to the nearest surface of an object in a 3D volume. The object is defined by voxels 
 * with values below a threshold (typically 0.5, representing the boundary).
 * The input image is modified and returns the distance measure.
 *
 * Parameters:
 *   V: The input image (float) represented as a 3D volume. Voxels with zero value 
 *       are considered non-elements of the object.
 *   M: An uint16 mask defining the region in which the distance calculations 
 *       are performed. The mask should define a convex hull ensuring direct 
 *       connections between the object and the estimation voxels.
 *       If the mask is NULL the distance calculations are performed for the whole image
 *       without any mask.
 *   dims: An integer array representing the dimensions of the 3D volume (width, height, depth).
 *   voxelsize: An float array representing the size of each voxel in the 3D volume, in millimeters.
 *              If voxelsize==NULL then the default voxelsize is 1mm.
 *   replace: An integer flag indicating whether to replace the values inside the mask with 
 *            their neighboring values. If set to 0, the function returns the voxel 
 *            distance; otherwise, the original values outside the mask are retained, 
 *            and inside the mask, the values are replaced.
 *
 * The function performs a forward and backward pass to calculate the minimum distance
 * from each voxel within the mask to the object's surface. The distances are calculated
 * using the voxel sizes to correct for anisotropy. However, the distance is defined in voxels
 * if the 'voxelsize' is NULL.
 *
 * If the 'replace' parameter is set, the function additionally replaces the values inside 
 * the mask with the values from their nearest neighbor outside the object boundary. This 
 * can be useful in scenarios where the original values inside the object are to be retained 
 * for further analysis.
 *
 */
void vbdist(float *V, unsigned char *M, int dims[3], double *voxelsize, int replace) 
{
    
    /* main information about input data (size, dimensions, ...) */
    const int nvox = dims[0]*dims[1]*dims[2];
    const int x    = dims[0];
    const int y    = dims[1];
    const int xy   = x*y;
    
    float s1, s2, s3;

    // use default voxel size of 1mm to get distance in voxel-space 
    if (!voxelsize) {
        s1 = s2 = s3 = 1.0;
    } else {
        s1 = (float)fabs(voxelsize[0]);
        s2 = (float)fabs(voxelsize[1]);
        s3 = (float)fabs(voxelsize[2]);
    }
    const float s12  = (float) sqrt((double)s1*s1   + s2*s2); /* xy - voxel size */
    const float s13  = (float) sqrt((double)s1*s1   + s3*s3); /* xz - voxel size */
    const float s23  = (float) sqrt((double)s2*s2   + s3*s3); /* yz - voxel size */
    const float s123 = (float) sqrt((double)s12*s12 + s3*s3); /* xyz - voxel size */
    
    /* indices of the neighbor Ni (index distance) and euclidean distance NW */
    const int NI[] = { 0, -1,-x+1, -x,-x-1, -xy+1,-xy,-xy-1, -xy+x+1,-xy+x,-xy+x-1, -xy-x+1,-xy-x,-xy-x-1}; 
    const float  ND[] = {0.0, s1, s12, s2, s12, s13, s3,    s13, s123, s23, s123, s123, s23, s123};
    const int sN = sizeof(NI)/4; /* division by 4 to get from the number of bytes to the number of elements */ 
    float DN[sN];
    float DNm = FLT_MAX;
    int  i, n, ni, DNi;
    int  u,v,w,nu,nv,nw; 
    
    /* data */
    float        *D, *buffer;
    unsigned int *I;
    
    I = (unsigned int *)malloc(sizeof(unsigned int)*nvox);
    D = (float *)malloc(sizeof(float)*nvox);
    
    /* save original input in buffer if we want to replace values */
    if (replace > 0) {
        buffer = (float *)malloc(sizeof(float)*nvox);
        memcpy(buffer,V,nvox*sizeof(float));  
    }
    
    if (!D || !I || ((replace > 0) && !buffer)) {
        fprintf(stderr,"Memory allocation error\n");
        exit(EXIT_FAILURE);
    }

    /* Initiaize mask with ones if not defined */
    if (!M) {
        M = (unsigned char *)malloc(sizeof(unsigned char)*nvox);
        if (!M) {
            fprintf(stderr,"Memory allocation error\n");
            exit(EXIT_FAILURE);
        }
        for (i = 0; i<nvox; i++)
            M[i] = 1;
    }
    
    /* initialisation of D and I */
    for (i = 0; i<nvox; i++) {
        if ((roundf(V[i])<0.5) || isnan(V[i])) D[i] = FLT_MAX; else D[i] = 0.0; 
        I[i] = (unsigned int)i;
    }
    
    /* forward direction that consider all points smaller than i */
    for (i = 0; i<nvox; i++) {
        if ((D[i]>0) && (M[i]>0)) {
            ind2sub(i,&u,&v,&w,xy,x);
            
            /* read neighbor values */
            for (n=0; n<sN; n++) {
                ni = i + NI[n];
                ind2sub(ni,&nu,&nv,&nw,xy,x);
                if ((ni<0) || (ni>=nvox) || (abs(nu-u)>1) || (abs(nv-v)>1) || (abs(nw-w)>1)) ni=i;
                DN[n] = D[ni] + ND[n];
            }
    
            /* find minimum distance within the neighborhood */
            pmin(DN,sN,&DNm,&DNi);
    
            /* update values */
            if (DNi>0) {
                I[i] = (unsigned int) I[i+NI[DNi]];
                D[i] = DNm; 
                ind2sub((int)I[i],&nu,&nv,&nw,xy,x); 
                D[i] = (float)sqrt(pow((double)(u-nu)*s1,2) + pow((double)(v-nv)*s2,2) + pow((double)(w-nw)*s3,2));
            }
         }
    }
    
    /* backward direction that consider all points larger than i */
    for (i=nvox-1;i>=0;i--) {
        if ((D[i]>0) && (M[i]>0)) {
            ind2sub(i,&u,&v,&w,xy,x);
        
            /* read neighbour values */
            for (n=0; n<sN; n++) {
                ni = i - NI[n];
                ind2sub(ni,&nu,&nv,&nw,xy,x);
                if ((ni<0) || (ni>=nvox) || (abs(nu-u)>1) || (abs(nv-v)>1) || (abs(nw-w)>1)) ni=i;
                DN[n] = D[ni] + ND[n];
            }
        
            /* find minimum distance within the neighborhood */
            pmin(DN,sN,&DNm,&DNi);
        
            /* update values */
            if (DNi>0) {
                I[i] = (unsigned int) I[i-NI[DNi]];
                D[i] = DNm; 
                ind2sub((int)I[i],&nu,&nv,&nw,xy,x); 
                D[i] = (float)sqrt(pow((double)(u-nu)*s1,2) + pow((double)(v-nv)*s2,2) + pow((double)(w-nw)*s3,2));
            }
        }
    }

    /* finally return output to original variable V */
    for (i = 0; i<nvox; i++)
        V[i] = ((M[i]==0) || (D[i] == FLT_MAX)) ? 0.0 : D[i];

    /* finally replace values inside mask */
    if (replace > 0) {
        for (i = 0; i < nvox; ++i)
            V[i] = buffer[I[i]];
        free(buffer);
    }
        
    free(D);
    free(I);
}

/**
 * laplace3R - Apply Laplace filter on a 3D volume.
 *
 * This function performs Laplace filtering on a 3D volume. It filters the volume 
 * within the intensity range defined by a mask until the changes are below a 
 * specified threshold.
 *
 * Parameters:
 *  - SEG: 3D single input matrix (volume to be filtered).
 *  - M: 3D volume that defines the filter area (mask).
 *  - dims: Array containing the dimensions of the volume.
 *  - TH: Threshold controlling the number of iterations (maximum change allowed after an iteration).
 *
 * Notes:
 *  - The function iterates until the maximum difference in the filtered volume is 
 *    less than the threshold or until it reaches the maximum number of iterations.
 *  - The function uses the neighboring values to calculate the Laplace filtering.
 */
void laplace3R(float *SEG, unsigned char *M, int dims[3], double TH) {
    const int x = dims[0];
    const int y = dims[1];

    const int z = dims[2];
    const int xy = x * y;
    const int nvox = x * y * z;
    
    // Indices of the neighbor and size of neighbors array
    const int NI[] = { -1, 1, -x, x, -xy, xy }; 
    const int sN = sizeof(NI) / sizeof(NI[0]);
    
    int i, n, u, v, w, nu, nv, nw, ni, iter = 0, maxiter = 2000;
    float Nn, diff, maxdiffi, maxdiff = 1.0;

    // Allocate memory for Laplace calculation
    float *L1 = (float *)malloc(sizeof(float) * nvox);
    float *L2 = (float *)malloc(sizeof(float) * nvox);
    unsigned char *LN = (unsigned char *)malloc(sizeof(unsigned char) * nvox);
    
    // Check for successful memory allocation
    if (!L1 || !L2 || !LN) {
        fprintf(stderr, "Memory allocation error\n");
        exit(EXIT_FAILURE);
    }
    
    // Initialization
    for (i = 0; i < nvox; i++) {
        L1[i] = isnan(SEG[i]) ? FLT_MAX : SEG[i];
        L2[i] = L1[i];
        LN[i] = M[i];
    }
    
    // Iterative Laplace filtering
    while (maxdiff > TH && iter < maxiter) {
        maxdiffi = 0;
        iter++;
        for (i = 0; i < nvox; i++) {
            if (M[i] && LN[i]) {  
                ind2sub(i, &u, &v, &w, xy, x);
    
                // Read neighbor values
                L2[i] = 0.0;
                Nn = 0.0;
                for (n = 0; n < sN; n++) {
                    ni = i + NI[n];
                    ind2sub(ni, &nu, &nv, &nw, xy, x);
                    if (((ni<0) || (ni>=nvox) || (abs(nu-u)>1) || (abs(nv-v)>1) || (abs(nw-w)>1) || (L1[ni]==-FLT_MAX) || (L1[ni]==FLT_MAX))==0) {
                        L2[i] += L1[ni];
                        Nn++;
                    }
                }
    
                L2[i] = Nn > 0 ? L2[i] / (float)Nn : L1[i];
                diff = fabs(L1[i] - L2[i]);
                if (diff > (TH / 10.0)) { 
                    for (n=0;n<sN;n++) {
                        ni = i + NI[n];
                        ind2sub(ni,&nu,&nv,&nw,xy,x);
                        if (((ni<0) || (ni>=nvox) || (abs(nu-u)>1) || (abs(nv-v)>1) || (abs(nw-w)>1) || (L1[ni]==-FLT_MAX) || (L1[ni]==FLT_MAX))==0) 
                            LN[ni] = 1; /* if I change his neigbors it has to be recalculated */
                    }                
                }
                LN[i] = 0;
                if (maxdiffi<diff) maxdiffi = diff; 
            }
        }
        maxdiff = maxdiffi;
        
        // Update L1 with the new values from L2
        for (i = 0; i<nvox; i++) 
            L1[i] = L2[i];
    }
    
    // Copy the final result back into SEG
    for (i = 0; i<nvox; i++) 
        SEG[i] = L1[i];

    // Free allocated memory
    free(L1);
    free(L2);
    free(LN);
}

void distclose_float(float *vol, int dims[3], double voxelsize[3], int niter, double th)
{
    float *buffer;
    int i,x,y,z,j,band,dims2[3];
    float max_vol;
    int nvox2,nvox = dims[0]*dims[1]*dims[2];
    
    if (niter < 1) return;

    for (i = 0; i<nvox; i++) max_vol = MAX(max_vol,vol[i]);
    th *= (double)max_vol;

    /* add band with zeros to image to avoid clipping */    
    band = niter;
    for (i = 0;i<3;i++) dims2[i] = dims[i] + 2*band;
    nvox2 = dims2[0]*dims2[1]*dims2[2];

    buffer = (float *)malloc(sizeof(float)*dims2[0]*dims2[1]*dims2[2]);

    if (!buffer) {
        fprintf(stderr,"Memory allocation error\n");
        exit(EXIT_FAILURE);
    }
    
    memset(buffer,0,sizeof(float)*nvox2);
    
    /* threshold input */
    for (z=0;z<dims[2];z++) for (y=0;y<dims[1];y++) for (x=0;x<dims[0];x++) 
        buffer[sub2ind(x+band,y+band,z+band,dims2)] = (vol[sub2ind(x,y,z,dims)]>(float)th);

    vbdist(buffer, NULL, dims2, voxelsize, 0);
    for (i = 0;i<nvox2;i++)
        buffer[i] = buffer[i] > (float)niter;

    vbdist(buffer, NULL, dims2, voxelsize, 0);
    for (i = 0;i<nvox2;i++)
        buffer[i] = buffer[i] > (float)niter;

    /* return image */
    for (z=0;z<dims[2];z++) for (y=0;y<dims[1];y++) for (x=0;x<dims[0];x++) 
        vol[sub2ind(x,y,z,dims)] = buffer[sub2ind(x+band,y+band,z+band,dims2)];
        
    free(buffer);
}

void distclose(void *data, int dims[3], double voxelsize[3], int niter, double th, int datatype)
{
    int nvox;
    float *buffer;
   
    nvox = dims[0]*dims[1]*dims[2];
    buffer = (float *)malloc(sizeof(float)*nvox);
    
    /* check success of memory allocation */
    if (!buffer) {
        printf("Memory allocation error\n");
        exit(EXIT_FAILURE);
    }
   
    convert_input_type(data, buffer, nvox, datatype);
    distclose_float(buffer, dims, voxelsize, niter, th);
    convert_output_type(data, buffer, nvox, datatype);
    
    free(buffer);
}

void distopen_float(float *vol, int dims[3], double voxelsize[3], double dist, double th)
{
    float *buffer;
    int i, j;
    float max_vol;
    int nvox = dims[0]*dims[1]*dims[2];
    
    if (dist == 0.0) return;

    for (i = 0; i<nvox; i++) max_vol = MAX(max_vol,vol[i]);
    th *= (double)max_vol;
    
    buffer = (float *)malloc(sizeof(float)*nvox);

    if (!buffer) {
        fprintf(stderr,"Memory allocation error\n");
        exit(EXIT_FAILURE);
    }
    
    /* threshold input */
    for (i = 0; i<nvox; i++)
        buffer[i] = 1.0 - ((float)vol[i]>th);

    vbdist(buffer, NULL, dims, voxelsize, 0);
    for (i = 0; i<nvox; i++)
        buffer[i] = buffer[i] > (float)dist;

    vbdist(buffer, NULL, dims, voxelsize, 0);
    for (i = 0; i<nvox; i++)
        buffer[i] = buffer[i] <= (float)dist;

    /* return image */
    for (i = 0; i<nvox; i++)
        vol[i] = buffer[i];

    free(buffer);
}

void distopen(void *data, int dims[3], double voxelsize[3], int niter, double th, int datatype)
{
    int nvox;
    float *buffer;
   
    nvox = dims[0]*dims[1]*dims[2];
    buffer = (float *)malloc(sizeof(float)*nvox);
    
    /* check success of memory allocation */
    if (!buffer) {
        printf("Memory allocation error\n");
        exit(EXIT_FAILURE);
    }
   
    convert_input_type(data, buffer, nvox, datatype);
    distopen_float(buffer, dims, voxelsize, niter, th);
    convert_output_type(data, buffer, nvox, datatype);
    
    free(buffer);
}

void morph_erode_float(float *vol, int dims[3], int niter, double th)
{
    double filt[3]={1,1,1};
    int i,j;
    float max_vol;
    int nvox = dims[0]*dims[1]*dims[2];
    
    if (niter < 1) return;

    for (i = 0; i<nvox; i++) max_vol = MAX(max_vol,vol[i]);
    th *= (double)max_vol;

    /* threshold input */
    for (j=0;j<nvox;j++)
        vol[j] = vol[j]>(float)th;

    for (i = 0;i<niter;i++) {
        convxyz_float(vol,filt,filt,filt,3,3,3,-1,-1,-1,vol,dims);
        for (j=0;j<nvox;j++)
            vol[j] = vol[j]>=9.0;
    }
}

void morph_erode(void *data, int dims[3], int niter, double th, int datatype)
{
    int nvox;
    float *buffer;
   
    nvox = dims[0]*dims[1]*dims[2];
    buffer = (float *)malloc(sizeof(float)*nvox);
    
    /* check success of memory allocation */
    if (!buffer) {
        printf("Memory allocation error\n");
        exit(EXIT_FAILURE);
    }
   
    convert_input_type(data, buffer, nvox, datatype);
    morph_erode_float(buffer, dims, niter, th);
    convert_output_type(data, buffer, nvox, datatype);
    
    free(buffer);
}

void morph_dilate_float(float *vol, int dims[3], int niter, double th)
{
    double filt[3]={1,1,1};
    int i,x,y,z,j,band,dims2[3];
    float max_vol;
    unsigned char *buffer;
    int nvox = dims[0]*dims[1]*dims[2];
    
    if (niter < 1) return;

    for (i = 0; i<nvox; i++) max_vol = MAX(max_vol,vol[i]);
    th *= (double)max_vol;

    /* add band with zeros to image to avoid clipping */    
    band = niter;
    for (i = 0;i<3;i++) dims2[i] = dims[i] + 2*band;

    buffer = (unsigned char *)malloc(sizeof(unsigned char)*dims2[0]*dims2[1]*dims2[2]);

    if (!buffer) {
        fprintf(stderr,"Memory allocation error\n");
        exit(EXIT_FAILURE);
    }
    
    memset(buffer,0,sizeof(unsigned char)*dims2[0]*dims2[1]*dims2[2]);
    
    /* threshold input */
    for (x=0;x<dims[0];x++) for (y=0;y<dims[1];y++) for (z=0;z<dims[2];z++) 
        buffer[sub2ind(x+band,y+band,z+band,dims2)] = (unsigned char)((double)vol[sub2ind(x,y,z,dims)]>th);

    for (i = 0;i<niter;i++) {
        convxyz_uint8(buffer,filt,filt,filt,3,3,3,-1,-1,-1,buffer,dims2);
        for (j=0; j<dims2[0]*dims2[1]*dims2[2]; j++) 
            buffer[j] = buffer[j]>0;
    }

    /* return image */
    for (x=0;x<dims[0];x++) for (y=0;y<dims[1];y++) for (z=0;z<dims[2];z++) 
        vol[sub2ind(x,y,z,dims)] = (float)buffer[sub2ind(x+band,y+band,z+band,dims2)];
        
    free(buffer);
    
}

void morph_dilate(void *data, int dims[3], int niter, double th, int datatype)
{
    int nvox;
    float *buffer;
   
    nvox = dims[0]*dims[1]*dims[2];
    buffer = (float *)malloc(sizeof(float)*nvox);

    /* check success of memory allocation */
    if (!buffer) {
        printf("Memory allocation error\n");
        exit(EXIT_FAILURE);
    }
   
    convert_input_type(data, buffer, nvox, datatype);
    morph_dilate_float(buffer, dims, niter, th);
    convert_output_type(data, buffer, nvox, datatype);
    
    free(buffer);
}

void morph_close_float(float *vol, int dims[3], int niter, double th)
{
    double filt[3]={1,1,1};
    unsigned char *buffer;
    int i,x,y,z,j,band,dims2[3];
    float max_vol;
    int nvox = dims[0]*dims[1]*dims[2];
    
    if (niter < 1) return;

    for (i = 0; i<nvox; i++) max_vol = MAX(max_vol,vol[i]);
    th *= (double)max_vol;

    /* add band with zeros to image to avoid clipping */    
    band = niter;
    for (i = 0;i<3;i++) dims2[i] = dims[i] + 2*band;

    buffer = (unsigned char *)malloc(sizeof(unsigned char)*dims2[0]*dims2[1]*dims2[2]);

    if (!buffer) {
        fprintf(stderr,"Memory allocation error\n");
        exit(EXIT_FAILURE);
    }
    
    memset(buffer,0,sizeof(unsigned char)*dims2[0]*dims2[1]*dims2[2]);
    
    /* threshold input */
    for (x=0;x<dims[0];x++) for (y=0;y<dims[1];y++) for (z=0;z<dims[2];z++) 
        buffer[sub2ind(x+band,y+band,z+band,dims2)] = (unsigned char)((double)vol[sub2ind(x,y,z,dims)]>th);
                
    /* dilate */
    for (i = 0;i<niter;i++) {
        convxyz_uint8(buffer,filt,filt,filt,3,3,3,-1,-1,-1,buffer,dims2);
        for (j=0; j<dims2[0]*dims2[1]*dims2[2]; j++) 
            buffer[j] = (buffer[j]>0);
    }

    /* erode */
    for (i = 0;i<niter;i++) {
        convxyz_uint8(buffer,filt,filt,filt,3,3,3,-1,-1,-1,buffer,dims2);
        for (j=0; j<dims2[0]*dims2[1]*dims2[2]; j++) 
            buffer[j] = (buffer[j]>=9);
    }

    /* return image */
    for (x=0;x<dims[0];x++) for (y=0;y<dims[1];y++) for (z=0;z<dims[2];z++) 
        vol[sub2ind(x,y,z,dims)] = (float)buffer[sub2ind(x+band,y+band,z+band,dims2)];
        
    free(buffer);
}

void morph_close(void *data, int dims[3], int niter, double th, int datatype)
{
    int nvox;
    float *buffer;
   
    nvox = dims[0]*dims[1]*dims[2];
    buffer = (float *)malloc(sizeof(float)*nvox);

    /* check success of memory allocation */
    if (!buffer) {
        printf("Memory allocation error\n");
        exit(EXIT_FAILURE);
    }
   
    convert_input_type(data, buffer, nvox, datatype);
    morph_close_float(buffer, dims, niter, th);
    convert_output_type(data, buffer, nvox, datatype);
    
    free(buffer);
}

void morph_open_float(float *vol, int dims[3], int niter, double th, int keep_values)
{
    unsigned char *buffer;
    double filt[3]={1,1,1};
    int i, j, nvox;
    float max_vol;
    
    if (niter < 1) return;

    nvox = dims[0]*dims[1]*dims[2];
    for (i = 0; i<nvox; i++) max_vol = MAX(max_vol,vol[i]);
    th *= (double)max_vol;

    buffer = (unsigned char *)malloc(sizeof(unsigned char)*nvox);

    if (!buffer) {
        fprintf(stderr,"Memory allocation error\n");
        exit(EXIT_FAILURE);
    }

    /* threshold input */
    for (j=0;j<nvox;j++)
        buffer[j] = (unsigned char)(vol[j]>th);

    for (i = 0;i<niter;i++) {
        convxyz_uint8(buffer,filt,filt,filt,3,3,3,-1,-1,-1,buffer,dims);
        for (j=0;j<nvox;j++)
            buffer[j] = (buffer[j]>=9);
    }

    for (i = 0;i<niter;i++) {
        convxyz_uint8(buffer,filt,filt,filt,3,3,3,-1,-1,-1,buffer,dims);
        for (j=0; j<nvox; j++) 
            buffer[j] = (buffer[j]>0);
    }

    // sometimes it's helpful to keep the original volume, but set areas to zero
    // where opening was successful
    if (keep_values > 0) {
        for (i = 0; i<nvox; i++)
            if (buffer[i] == 0) vol[i] = 0.0;
    } else {
        for (i = 0; i<nvox; i++)
            vol[i] = (float)buffer[i];
    }
          
    free(buffer);
}

void morph_open(void *data, int dims[3], int niter, double th, int datatype)
{
    int nvox;
    float *buffer;
   
    nvox = dims[0]*dims[1]*dims[2];
    buffer = (float *)malloc(sizeof(float)*nvox);

    /* check success of memory allocation */
    if (!buffer) {
        printf("Memory allocation error\n");
        exit(EXIT_FAILURE);
    }
   
    convert_input_type(data, buffer, nvox, datatype);
    morph_open_float(buffer, dims, niter, th, 0);
    convert_output_type(data, buffer, nvox, datatype);
    
    free(buffer);
}

/**
 * Function: subsample_float
 * -------------------------
 * Subsamples a 3D array to a different size using trilinear interpolation.
 * 
 * This function takes an input 3D array and resizes it to a specified dimension
 * using trilinear interpolation. Trilinear interpolation is a method of 
 * interpolating within a 3D grid. It works by interpolating data points on a 
 * linear scale along each of the three axes (X, Y, and Z).
 *
 * The function first checks if the output array is already allocated.
 * If not, it allocates the necessary memory. It then calculates the scaling
 * factors for each dimension based on the input and output dimensions.
 * For each point in the output array, it computes the interpolated value using
 * the surrounding points in the input array.
 * 
 * Parameters:
 *   in - Pointer to the input float array.
 *   out - Pointer to the output float array. If this is NULL, the function will
 *         allocate memory for the output array.
 *   dims - An array of 3 integers representing the dimensions (width, height, depth)
 *          of the input array.
 *   dims_samp - An array of 3 integers representing the dimensions of the output array.
 * 
 * Note:
 *   - If 'out' is NULL, the function will allocate memory, which should be freed
 *     by the caller to avoid memory leaks.
 *   - Points outside the bounds of the input array are assigned a value of 0 in the 
 *     output array.
 */
void subsample_float(float *in, float *out, int dims[3], int dims_samp[3])
{
    int i, x, y, z;
    double k111, k112, k121, k122, k211, k212, k221, k222;
    double dx1, dx2, dy1, dy2, dz1, dz2, xi, yi, zi, samp[3];
    int off1, off2, xcoord, ycoord, zcoord;
    
    // Allocate if not already done
    if (!out) {
        out = (float *)malloc(sizeof(float)*dims_samp[0]*dims_samp[1]*dims_samp[2]);
        if (!out) {
            printf("Memory allocation error\n");
            exit(EXIT_FAILURE);
        }
    }
    
    // Estimate scaling factor for each dimension
    for (i = 0; i < 3; i++) {
        if (dims_samp[i] > dims[i]) samp[i] = ceil((double)dims_samp[i]/(double)dims[i]);
        else samp[i] = 1.0/(ceil((double)dims[i]/(double)dims_samp[i]));
    }
    
    for (z = 0; z < dims_samp[2]; z++) {
        zi = 1.0 + (double)z / samp[2];
        for (y = 0; y < dims_samp[1]; y++) {
            yi = 1.0 + (double)y / samp[1];
            for (x = 0; x < dims_samp[0]; x++) {
                xi = 1.0 + (double)x / samp[0];
                i = z * dims_samp[0] * dims_samp[1] + y * dims_samp[0] + x;

                if (zi >= 0 && zi < dims[2] && yi >= 0 && yi < dims[1] && xi >= 0 && xi < dims[0]) {
                    xcoord = (int)floor(xi); dx1 = xi - (double)xcoord; dx2 = 1.0 - dx1;
                    ycoord = (int)floor(yi); dy1 = yi - (double)ycoord; dy2 = 1.0 - dy1;
                    zcoord = (int)floor(zi); dz1 = zi - (double)zcoord; dz2 = 1.0 - dz1;

                    off1 = xcoord-1 + dims[0]*(ycoord-1 + dims[1]*(zcoord-1));
                    k222 = (double)in[off1]; k122 = (double)in[off1+1]; off2 = off1 + dims[0];
                    k212 = (double)in[off2]; k112 = (double)in[off2+1]; off1+= dims[0] * dims[1];
                    k221 = (double)in[off1]; k121 = (double)in[off1+1]; off2 = off1 + dims[0];
                    k211 = (double)in[off2]; k111 = (double)in[off2+1];

                    out[i] = (float)((((k222*dx2 + k122*dx1)*dy2 + (k212*dx2 + k112*dx1)*dy1))*dz2
                        + (((k221*dx2 + k121*dx1)*dy2 + (k211*dx2 + k111*dx1)*dy1))*dz1);
                                 
                } else out[i] = 0;
            }
        }
    }
}

void subsample3(void *in, void *out, int dims[3], int dims_samp[3], int datatype)
{
    int nvox_samp, nvox;
    float *vol_samp, *buffer;
    
    nvox      = dims[0]*dims[1]*dims[2];
    nvox_samp = dims_samp[0]*dims_samp[1]*dims_samp[2];
    
    vol_samp  = (float *)malloc(sizeof(float)*nvox_samp);
    buffer    = (float *)malloc(sizeof(float)*nvox);

    /* check success of memory allocation */
    if (!buffer || !vol_samp) {
        printf("Memory allocation error\n");
        exit(EXIT_FAILURE);
    }
   
    convert_input_type(in, buffer, nvox, datatype);
    subsample_float(buffer, vol_samp, dims, dims_samp);   
    convert_output_type(out, vol_samp, nvox_samp, datatype);

    free(vol_samp);
    free(buffer);
}

void smooth_subsample_float(float *vol, int dims[3], double voxelsize[3], double s[3], int use_mask, int samp)
{
    int i, nvox_samp, nvox;
    int dims_samp[3];
    float *vol_samp;
    double voxelsize_samp[3];
    
    /* define grid dimensions */
    for (i = 0; i<3; i++) dims_samp[i] = (int) ceil((dims[i]-1)/((double) samp))+1;
    for (i = 0; i<3; i++) voxelsize_samp[i] = voxelsize[i]*((double)dims[i]/(double)dims_samp[i]);

    nvox      = dims[0]*dims[1]*dims[2];
    nvox_samp = dims_samp[0]*dims_samp[1]*dims_samp[2];
    vol_samp  = (float *)malloc(sizeof(float)*nvox_samp);

    if (!vol_samp) {
        fprintf(stderr,"Memory allocation error\n");
        exit(EXIT_FAILURE);
    }

    subsample_float(vol, vol_samp, dims, dims_samp);   
    smooth_float(vol_samp, dims_samp, voxelsize_samp, s, use_mask);
    subsample_float(vol_samp, vol, dims_samp, dims);   

    free(vol_samp);
}

void smooth_subsample3(void *data, int dims[3], double voxelsize[3], double s[3], int use_mask, int samp, int datatype)
{
    int i, nvox_samp, nvox;
    int dims_samp[3];
    float *vol_samp;
    float *buffer;
    double voxelsize_samp[3];
    
    /* define grid dimensions */
    for (i = 0; i<3; i++) dims_samp[i] = (int) ceil((dims[i]-1)/((double) samp))+1;
    for (i = 0; i<3; i++) voxelsize_samp[i] = voxelsize[i]*((double)dims[i]/(double)dims_samp[i]);

    nvox      = dims[0]*dims[1]*dims[2];
    nvox_samp = dims_samp[0]*dims_samp[1]*dims_samp[2];
    vol_samp  = (float *)malloc(sizeof(float)*nvox_samp);
    buffer    = (float *)malloc(sizeof(float)*nvox);

    /* check success of memory allocation */
    if (!buffer || !vol_samp) {
        printf("Memory allocation error\n");
        exit(EXIT_FAILURE);
    }
   
    convert_input_type(data, buffer, nvox, datatype);
    subsample_float(buffer, vol_samp, dims, dims_samp);   
    smooth_float(vol_samp, dims_samp, voxelsize_samp, s, use_mask);
    subsample_float(vol_samp, buffer, dims_samp, dims);
    convert_output_type(data, buffer, nvox, datatype);

    free(vol_samp);
    free(buffer);
}

/**
 * correct_bias_label - Performs bias correction on MRI images based on specified labels.
 *
 * This function applies a bias correction process to an MRI image dataset. The correction
 * is performed based on the label values assigned to each voxel, allowing different
 * treatments for different tissue types (e.g., white matter, gray matter).
 *
 * src: Pointer to the source image data (float array).
 *       This array is modified in place with the bias-corrected values.
 *
 * label: Pointer to the label array (unsigned char array) indicating different
 *         tissue types in the MRI scan. Different values in this array represent
 *         different tissues such as white matter, gray matter, etc.
 *
 * dims: Pointer to an integer array of size 3, indicating the dimensions of the
 *        MRI volume (e.g., [width, height, depth]).
 *
 * voxelsize: Pointer to a double array of size 3, indicating the size of each
 *             voxel in the MRI data (e.g., [size_x, size_y, size_z]).
 *
 * bias_fwhm: A double value indicating the full-width half-maximum (FWHM) of
 *             the Gaussian kernel used for smoothing in the bias correction
 *             process.
 *
 * label_th: An integer specifying the threshold label value for selecting
 *            specific tissues for bias correction. For example, using label_th = 2
 *            may focus the correction on white matter only.
 *
 * The function first calculates the mean intensity for each label class and then
 * estimates a bias field based on the ratio between actual voxel values and their
 * respective label class means. Morphological operations are used to refine the 
 * brain mask for the bias estimation. The bias field is then used to adjust the 
 * voxel intensities in the source image.
 *
 */
void correct_bias_label(float *src, float *biasfield, unsigned char *label, int *dims, double *voxelsize, double bias_fwhm, int label_th)
{
    int i, j, nvox, nvoxr, n[MAX_NC], n_classes = 0;
    int dimsr[3];
    unsigned char *mask, *maskr;
    float *biasfieldr;
    double mean_label[MAX_NC], mean_bias, voxelsizer[3];
    double fwhm[] = {bias_fwhm, bias_fwhm, bias_fwhm};

    // downsample resolution by factor samp to increase speed
    int samp = 2;
    
    int set_bias_to_zero = 1;
    
    // define grid dimensions
    for (i = 0; i<3; i++) {
        dimsr[i] = (int) ceil((dims[i]-1)/((double) samp))+1;
        voxelsizer[i] = voxelsize[i]*samp;
    }

    nvox  = dims[0] *  dims[1]  * dims[2];
    nvoxr = dimsr[0] * dimsr[1] * dimsr[2];
    
    biasfieldr = (float *)malloc(sizeof(float)*nvoxr);
    mask  = (unsigned char *)malloc(sizeof(unsigned char)*nvox);
    maskr = (unsigned char *)malloc(sizeof(unsigned char)*nvoxr);
    
    if (!mask || !maskr || !biasfieldr) {
        printf("Memory allocation error\n");
        exit(EXIT_FAILURE);
    }

    /* get number of classes by checking maximum label value */
    for (i = 0; i < nvox; i++)
        n_classes = MAX((int)label[i], n_classes);

    /* initialize parameters */
    for (i = 0; i < n_classes; i++) {
        n[i] = 0;
        mean_label[i] = 0.0;
    }
    
    /* estimate mean for each label class */
    for (i = 0; i < nvox; i++) {
        if (label[i] == 0) continue;
        n[label[i]-1]++;
        mean_label[label[i]-1] += (double)src[i];
    }
    for (i = 0; i < n_classes; i++) mean_label[i] /= (double)n[i];


    /* only use defined labels (i.e. using label_th) for bias estimation
     * use label_th = 2 for focussing on WM only */
    for (i = 0; i < nvox; i++)
        mask[i] = (label[i] >= label_th) ? 1 : 0;

    /* get bias field by ratio between actual values and respective mean of label class */
    for (i = 0; i < nvox; i++)
        if (mask[i] > 0)
            biasfield[i] += (src[i] / (float)mean_label[label[i]-1]);
        else biasfield[i] = 0.0;

    /* downsample to lower resolution */
    subsample3(biasfield, biasfieldr, dims, dimsr, DT_FLOAT32);
    subsample3(mask, maskr, dims, dimsr, DT_UINT8);
    
    // we have to ensure that thin gyri will be still covered if we use WM only
    // which later makes use of min/max values also necessary
    if (label_th > GM) {
      morph_dilate(maskr, dimsr, 1, 0.5, DT_UINT8);
          
      // use maximum/minimum values to exclude areas from other tissue classes
      // check whether mean WM is larger than GM (i.e. T1w)
      if (mean_label[2] > mean_label[1])
          localstat3(biasfieldr, maskr, dimsr, 2, F_MAX, 4, 0, DT_FLOAT32);
      else
          localstat3(biasfieldr, maskr, dimsr, 2, F_MIN, 4, 0, DT_FLOAT32);
    }
    
    // median and iterative means reduces PVE effects
    localstat3(biasfieldr, maskr, dimsr, 1, F_MEDIAN, 1, 0, DT_FLOAT32);
    localstat3(biasfieldr, maskr, dimsr, 1, F_MEAN, 4, 0, DT_FLOAT32);
    
    /* invert mask because we need to estimate dist outside the original mask */
    for (i = 0; i < nvoxr; i++) {
        if (set_bias_to_zero == 1)
            if (maskr[i] == 0) biasfieldr[i] = 0.0;
        maskr[i] = 1 - maskr[i];
    }

    vbdist(biasfieldr, maskr, dimsr, NULL, 1);
    smooth3(biasfieldr, dimsr, voxelsizer, fwhm, 0, DT_FLOAT32);

    /* upsample to original resolution */
    subsample3(biasfieldr, biasfield, dimsr, dims, DT_FLOAT32);

    /* estimate mean of bias field inside label for mean-correction */
    mean_bias = get_masked_mean_array_float(biasfield, nvox, label);

    for (i = 0; i < nvox; i++) {
        biasfield[i] /= mean_bias;
        if (biasfield[i] != 0)
            src[i] /= biasfield[i];
    }

    free(mask);
    free(maskr);
    free(biasfieldr);
}

/**
 * correct_bias - Applies bias correction to MRI data.
 *
 * This function performs bias correction on MRI images, specifically focusing on
 * white matter (WM) and, if specified, on gray matter (GM) as well. It uses a
 * local adaptive segmentation approach for additional GM correction.
 *
 * src: Pointer to the source image data (float array).
 *       This array is modified in place with the bias-corrected values.
 *
 * label: Pointer to the label array (unsigned char array) indicating different
 *         tissue types in the MRI scan. Typically, different values in this array
 *         represent different tissues such as WM, GM, and CSF.
 *
 * dims: Pointer to an integer array of size 3, indicating the dimensions of the
 *        MRI volume (e.g., [width, height, depth]).
 *
 * voxelsize: Pointer to a double array of size 3, indicating the size of each
 *             voxel in the MRI data (e.g., [size_x, size_y, size_z]).
 *
 * bias_fwhm: A double value indicating the full-width half-maximum (FWHM) of
 *             the Gaussian kernel used for smoothing in the bias correction
 *             process. This parameter is primarily used for WM correction.
 *
 * weight_las: A double value indicating the amount of weighting local adaptive 
 *              segmentation (LAS) that should be applied for additional GM correction. 
 *              If non-zero, LAS is applied (use values 0..1).
 *
 * The function first applies WM bias correction. If `weight_las` is non-zero, it
 * then performs GM correction with a small smoothing factor, creates a distance 
 * map for subcortical regions, and applies a weighted average to optimize the 
 * bias correction in these areas. The function modifies the `src` array in place 
 * with the corrected values.
 *
 */
void correct_bias(float *src, float *biasfield, unsigned char *label, int *dims, double *voxelsize, double bias_fwhm, double weight_las, int square_image)
{
    int i, nvox;
    unsigned char *mask;
    float *src_subcortical, *dist, *buffer;
    float max_dist = -FLT_MAX, mx_image, scl;

    nvox = dims[0]*dims[1]*dims[2];

    /* allocate biasfield if necessary */
    if (!biasfield) {
        biasfield = (float *)malloc(sizeof(float)*nvox);
        if (!biasfield) {
            printf("Memory allocation error\n");
            exit(EXIT_FAILURE);
        }
    }
    
    /* square input image improve segmentation of CSF and GM */
    if (square_image) {
        mx_image = get_max_float(src, nvox, 0);
        for (i = 0; i < nvox; i++) src[i] *= src[i];
        scl = mx_image/get_max_float(src, nvox, 0);
        for (i = 0; i < nvox; i++) src[i] *= scl;
    }
      
    /* apply bias correction for WM only */
    correct_bias_label(src, biasfield, label, dims, voxelsize, bias_fwhm, WM);
    
    /* use local adaptive segmentation (LAS) and apply additional GM correction with
     * very small smoothing */
    if (weight_las > 0.0) {
        src_subcortical = (float *)malloc(sizeof(float)*nvox);
        dist = (float *)malloc(sizeof(float)*nvox);
        buffer = (float *)malloc(sizeof(float)*nvox);
        if (!src_subcortical || !dist || !buffer) {
            printf("Memory allocation error\n");
            exit(EXIT_FAILURE); 
        }

        /* apply bias correction for WM and GM */
        for (i = 0; i < nvox; i++) src_subcortical[i] = src[i];
        bias_fwhm = 3.0;
        correct_bias_label(src_subcortical, buffer, label, dims, voxelsize, bias_fwhm, GM);

        /* weight LAS correction */
        for (i = 0; i < nvox; i++)
            src[i] = weight_las*src_subcortical[i] + (1.0 - weight_las)*src[i];

        /* prepare mask for distance map */
        for (i = 0; i < nvox; i++) dist[i] = (label[i] > 0) ? 0.0 : 1.0;
        
        /* get distance to background */
        vbdist(dist, NULL, dims, NULL, 0);

        /* scale distance values to 0..1  */
        for (i = 0; i < nvox; i++) max_dist = MAX(dist[i], max_dist);
        for (i = 0; i < nvox; i++) dist[i] /= max_dist;
        
        /* use cubic distance weights to weight central regions (i.e. subcortical
           structures) more */
        for (i = 0; i < nvox; i++) dist[i] *= dist[i]*dist[i];

        /* apply weighted average (with squared weights) to maximize weighting 
         * for subcortical regions with large distances, otherwise WM-bias 
         * correction is more weighted */
        for (i = 0; i < nvox; i++)
            src[i] = dist[i]*src_subcortical[i] + (1.0 - dist[i])*src[i];
        
        free(buffer);
        free(dist);
        free(src_subcortical);
    }
}

void vol_approx(float *vol, int dims[3], double voxelsize[3], int samp)
{
    int i, nvoxr, nvox;
    int dimsr[3];
    float *volr, *buffer, *TAr;
    double voxelsizer[3];
    float min_vol = FLT_MAX, max_vol = -FLT_MAX;
    unsigned char *BMr, *BMr2;
    double threshold[2], prctile[2] = {5,95};
        
    /* define grid dimensions */
    for (i = 0; i<3; i++) {
        dimsr[i] = (int) ceil((dims[i]-1)/((double) samp))+1;
        voxelsizer[i] = voxelsize[i]*samp;
    }

    nvox   = dims[0]*dims[1]*dims[2];
    nvoxr  = dimsr[0]*dimsr[1]*dimsr[2];
    volr   = (float *)malloc(sizeof(float)*nvoxr);
    buffer = (float *)malloc(sizeof(float)*nvoxr);
    TAr    = (float *)malloc(sizeof(float)*nvoxr);
    BMr    = (unsigned char *)malloc(sizeof(unsigned char)*nvoxr);
    BMr2   = (unsigned char *)malloc(sizeof(unsigned char)*nvoxr);

    /* find values between 0.1% and 99.9% percentile */
    for (i = 0; i < nvox; ++i)
        vol[i] -= 0;

    for (i = 0; i < nvox; ++i) {
        min_vol = MIN(vol[i], min_vol);
        max_vol = MAX(vol[i], max_vol);
    }
    
    /* only keep values between 5..95% percentiles 
       to remove extremes that occur at edges */
    get_prctile(vol, dims[0]*dims[1]*dims[2], threshold, prctile, 1);  
    for (i = 0; i < nvox; ++i)
        if ((vol[i] < threshold[0]) || (vol[i] > threshold[1])) vol[i] = 0;;

    /* scale vol to 0..1 but keep zero-background */
    for (i = 0; i < nvox; ++i)
        if (vol[i] !=0) min_vol = MIN(vol[i], min_vol);
    for (i = 0; i < nvox; ++i) {
        if (vol[i] != 0) {
            vol[i] -= min_vol;   
            max_vol = MAX(vol[i], max_vol);
        }
        vol[i] /= max_vol;
    }
    
    /* downsample to lower resolution */
    subsample_float(vol, volr, dims, dimsr);

    /* create mask by closing holes */ 
    for (i = 0; i < nvoxr; ++i) BMr[i] = volr[i] > 0;
    //morph_close(BMr, dimsr, 20, 0, DT_UINT8);

    /* vbdist to fill values in background with neighbours */
    memcpy(buffer,volr,nvoxr*sizeof(float));    
    vbdist(buffer, BMr, dimsr, NULL, 1);
    for (i = 0; i < nvoxr; ++i)
        TAr[i] = buffer[i];
    
    /* smooth values outside mask */
    double s[3] = {20,20,20};
    smooth_float(buffer, dimsr, voxelsizer, s, 0);
    for (i = 0; i < nvoxr; ++i)
        if (BMr[i] == 0) TAr[i] = buffer[i];

    /* rescue mask and create new mask that is only defined inside */
    for (i = 0; i < nvoxr; ++i) {
        BMr2[i] = BMr[i];
        BMr[i] = (BMr[i] > 0) && (volr[i] == 0);
    }
    
    laplace3R(TAr, BMr, dimsr, 0.4);
    median3(TAr, NULL, dimsr, 1, DT_FLOAT32);
    laplace3R(TAr, BMr, dimsr, 0.4);

    /* only keep TAr inside (closed) mask */
    for (i = 0; i < nvoxr; ++i)
        if ((BMr2[i] == 0) && (vol[i] == 0)) TAr[i] = 0.0;

    /* again apply vbdist to fill values in background with neighbours */
    vbdist(TAr, NULL, dimsr, NULL, 1);

    for (i = 0; i < 3; ++i) s[i] *= 2.0;
    smooth_float(buffer, dimsr, voxelsizer, s, 0);
    for (i = 0; i < nvoxr; ++i)
        if (BMr[i] == 0) TAr[i] = buffer[i];

    for (i = 0; i < nvoxr; ++i)
        BMr[i] = (BMr2[i] == 0);
    laplace3R(TAr, BMr, dimsr, 0.4);

    median3(TAr, NULL, dimsr, 1, DT_FLOAT32);

    for (i = 0; i < nvoxr; ++i)
        BMr[i] = (volr[i] == 0);
    laplace3R(TAr, BMr, dimsr, 0.4);

    memcpy(volr,TAr,nvoxr*sizeof(float));
    
    subsample_float(volr, vol, dimsr, dims);  

    /* get old range back */
    for (i = 0; i < nvox; ++i)
        vol[i] *= max_vol;

    free(volr);
    free(buffer);
    free(BMr);
    free(BMr2);
    free(TAr);
}

void cleanup_brain(unsigned char *prob, int dims[3], double voxelsize[3], int strength)
{
    double scale = 3.0/(voxelsize[0] + voxelsize[1] + voxelsize[2]);
    int niter, iter, nvox, th, th_erode, th_dilate, i, j, gm_offset, wm_offset, csf_offset;
    float *wm, *gm, *dist;
    float bp, tot, max_dist;
    double filt[3] = {0.75, 1.0, 0.75};
    
    niter     = 32;
    th_erode  = 153;                // initial threshold for erosion 0.6*255.0
    th_dilate = 38 + (strength*13); // threshold for dilation (0.15 + strength*0.05)*255
    th_dilate = 38 + (strength*13); // threshold for dilation (0.15 + strength*0.05)*255
        
    nvox = dims[0]*dims[1]*dims[2];

    // we need the offsets for the tissue classes for prob array
    gm_offset  = (int)(GM-1)*nvox;
    wm_offset  = (int)(WM-1)*nvox;
    csf_offset = (int)(CSF-1)*nvox;

    /* ensure that sum of filt is 1 */
    tot = 0.0;
    for (i = 0; i < 3; ++i)
        tot += filt[i];
    for (i = 0; i < 3; ++i)
        filt[i] /= tot;

    wm   = (float *)malloc(sizeof(float)*nvox);
    gm   = (float *)malloc(sizeof(float)*nvox);
    dist = (float *)malloc(sizeof(float)*nvox);

    if (!wm || !gm || !dist) {
        printf("Memory allocation error\n");
        exit(EXIT_FAILURE);
    }

    /* init with WM and GM */
    for (i = 0; i < nvox; ++i) {
        wm[i] = (float)prob[i + wm_offset];
        gm[i] = (float)prob[i + gm_offset];
    }

    /* prepare mask for distance map */
    for (i = 0; i < nvox; i++)
        dist[i] = ((wm[i] + gm[i] + (float)prob[i + csf_offset]) > 0.0) ? 0.0 : 1.0;

    // prevent that correction is made between the hemispheres by closing holes
    morph_close_float(dist, dims, 3, 0.1);
    
    // get distance to background in order to use this to weight cleanup more
    // to the outer parts of the label mask. This prevents cutting thin structures
    // in the cerebellum or central structures
    vbdist(dist, NULL, dims, NULL, 0);

    /* scale distance to 0..1  */
    for (i = 0; i < nvox; i++) max_dist = MAX(dist[i], max_dist);
    for (i = 0; i < nvox; i++) dist[i] /= max_dist;

    /* use cubic distance weights to weight central regions (i.e. subcortical
       structures) more */
    for (i = 0; i < nvox; i++) dist[i] *= dist[i];

    /* mask computed from gm and wm */
    /* erosions and conditional dilations */
    for (iter = 0; iter < niter; iter++) {
        
        /*  start with 2 iterations of erosions */
        th = (iter < 2) ? th_erode : th_dilate;
        
        /* b = (b>th).*(white+gray) */
        for (i = 0; i < nvox; ++i) {
            bp = (float)prob[i + gm_offset] + (float)prob[i + wm_offset];
            wm[i] = (wm[i] > (float)th) ? MIN(bp, 255.0) : 0.0;
        }

        /* convolve mask with filter width of 3 voxel */
        convxyz_float(wm, filt, filt, filt, 3, 3, 3, -1, -1, -1, wm, dims);           
    }
    
    // apply opening, but keep original values otherwise
    morph_open_float(wm, dims, MAX(1,roundf(scale*strength)), 0.1, 1);

    th = 13; // threshold for final cleanup 0.05*255
    for (i = 0; i < nvox; ++i) {
        bp = (float)prob[i + gm_offset] + (float)prob[i + wm_offset];
        bp = ((bp > (float)th && wm[i] > (float)th)) ? 1.0 : 0.0;
        
        // rescue original segmentation for average weighting
        wm[i] = (float)prob[i + wm_offset];
        gm[i] = (float)prob[i + gm_offset];
        
        // apply cleanup
        if (bp == 0) {
            prob[i + gm_offset] = 0;
            prob[i + wm_offset] = 0;
        }
    }

    // use weighted average w.r.t distance to weight cleaned and original segmentation
    for (i = 0; i < nvox; i++) {
        prob[i + gm_offset] = (unsigned char)roundf(dist[i]*gm[i] + (1.0 - dist[i])*(float)prob[i + gm_offset]);
        prob[i + wm_offset] = (unsigned char)roundf(dist[i]*wm[i] + (1.0 - dist[i])*(float)prob[i + wm_offset]);
    }

    // ensure that overall sum of all tissue classes per voxel is 255 for uint8
    for (i = 0; i < nvox; ++i) {
        tot = 0.0;
        for (j = 0; j < 3; ++j)
            tot += (float)prob[i + j*nvox];
        for (j = 0; j < 3; ++j)
            prob[i + j*nvox] = (unsigned char)(roundf((float)prob[i + j*nvox]/tot*255.0));
    }

    free(gm);
    free(wm);
    free(dist);
    
}

/**
 * keep_largest_cluster_float - Identifies the largest cluster in a 3D volume.
 *
 * This sub-function processes a 3D volume (such as an MRI scan) and identifies the largest
 * cluster of voxels that exceed a specified threshold. The function operates on a volume
 * represented as a linear array and modifies the input data to retain only the largest
 * cluster, setting all other values to zero.
 *
 * inData: Pointer to the float array representing the input 3D volume. The array should
 *          have 'numVoxels' elements. This array is modified in place, with only the largest
 *          cluster's voxels retained.
 *
 * thresh: A double value representing the threshold. Voxels with values equal to or greater
 *          than this threshold are considered part of a cluster.
 *
 * dims: Pointer to an integer array of size 3, indicating the dimensions of the volume
 *        (e.g., [width, height, depth]).
 *
 * min_size: Integer value that defines the minimum cluster size. To ignore this parameter
 *            and only keep the largest cluster, you can set min_size to <= 0.
 *
 * retain_above_th: Integer value that defines whether we set all smaller clusters to zero, 
 *                   but retain all original values (by setting retain_above_th to 1), or we 
 *                   only keep values in larger clusters that are then thresholded.
 *
 * conn: Integer value that sets the connection-scheme to 6, 18 or 26 neighbors.
 *
 * The function first initializes auxiliary arrays to track used voxels and to store the output
 * data. It then iterates through the volume, identifying connected voxels that form clusters
 * and exceed the threshold. The size of each cluster is determined, and the largest one is
 * identified. Finally, the input data is modified so that all smaller clusters are filled
 * with zeros and the original values are retained otherwise.
 *
 */
void keep_largest_cluster_float(float *inData, double thresh, int *dims, int min_size, int retain_above_th, int conn)
{
    float valToAdd;
    float *outData;
    int i, j, k, ti, tj, tk, maxi, maxj, maxk, mini, minj, mink, ind, ind1;
    int maxInd = 0, growingInd, growingCur;
    int numVoxels = dims[0] * dims[1] * dims[2];
    char *flagUsed;
    short *growing;

    if ((conn != 6) && (conn != 18) && (conn != 26)) {
        printf("Values for connectivity can be only 6, 18 or 26.\n");
        exit(EXIT_FAILURE);
    }

    int adiff = 1;   //connectivity 6: faces only
    if (conn == 18)
        adiff = 2;   //connectivity 18: faces+edges
    if (conn == 26)
        adiff = 3;   //connectivity 26: faces+edges+corners

    flagUsed = (char*)  malloc(numVoxels*sizeof(char));
    outData  = (float*) malloc(numVoxels*sizeof(float));
    growing  = (short*) malloc(numVoxels*3*sizeof(short));

    /* check success of memory allocation */
    if (!flagUsed || !outData || !growing) {
        printf("Memory allocation error\n");
        exit(EXIT_FAILURE);
    }

    for (i = 0; i < numVoxels; ++i) {
        flagUsed[i] = 0;
        outData[i] = 0.0;
    }

    for (k = 0; k < dims[2]; ++k) for (j = 0; j < dims[1]; ++j) for (i = 0; i < dims[0]; ++i)
    {
        ind = k*(dims[0]*dims[1]) + (j*dims[0]) + i;

        if (!flagUsed[ind] && inData[ind] >= thresh)
        {
            flagUsed[ind] = 1;
            growingInd = 3;
            growingCur = 0;
            growing[0] = i;
            growing[1] = j;
            growing[2] = k;

            while (growingCur < growingInd)
            {
                maxi = MIN(dims[0], growing[growingCur    ] + 2);
                maxj = MIN(dims[1], growing[growingCur + 1] + 2);
                maxk = MIN(dims[2], growing[growingCur + 2] + 2);

                mini = MAX(0, growing[growingCur    ] - 1);
                minj = MAX(0, growing[growingCur + 1] - 1);
                mink = MAX(0, growing[growingCur + 2] - 1);

                for (tk = mink; tk < maxk; ++tk) for (tj = minj; tj < maxj; ++tj) for (ti = mini; ti < maxi; ++ti)
                {
                    /*  Check for corners and edges depending on defined connectivity */
                    if (abs(ti - growing[growingCur]) + abs(tj - growing[growingCur + 1]) + abs(tk - growing[growingCur + 2]) > adiff)
                        continue;
                    ind1 = tk*(dims[0]*dims[1]) + (tj*dims[0]) + ti;

                    if (!flagUsed[ind1] && inData[ind1] >= thresh)
                    {
                        flagUsed[ind1] = 1;
                        growing[growingInd    ] = ti;
                        growing[growingInd + 1] = tj;
                        growing[growingInd + 2] = tk;
                        growingInd += 3;
                    }   
                }
                growingCur += 3;
            }

            growingCur = 0;

            while (growingCur < growingInd)
            {
                outData[growing[growingCur + 2]*(dims[0]*dims[1])+(growing[growingCur + 1]*dims[0])+growing[growingCur]] += (float)growingInd;
                growingCur += 3;
            }
        }
    }

    /* find maximum value which is the largest cluster or use defined minimum cluster size */
    if (min_size >= 0)
        for (i = 0; i < numVoxels; ++i) maxInd = MAX(outData[i], maxInd);
    else maxInd = min_size;
    
    /* set values with smaller clusters to zero */
    /* Depending on the parameter retain_above_th we either set smaller clusters to zero, 
     * but retain all original values otherwise, or we only keep values in larger clusters
     * that are then thresholded */
    for (i = 0; i < numVoxels; ++i) {
        if (retain_above_th) {
            if ((outData[i] > 0) && (outData[i] < maxInd)) inData[i] = 0.0;
        } else inData[i] = (outData[i] >= maxInd) ? inData[i] : 0;
    }

    free(flagUsed);
    free(growing);
    free(outData);
}

/**
 * keep_largest_cluster - Identifies the largest cluster in a 3D volume.
 *
 * This function calls keep_largest_cluster_float and converts datatypes of
 * input and output
 */
void keep_largest_cluster(void *data, double thresh, int *dims, int datatype, int min_size, int retain_above_th, int conn)
{
    int nvox;
    float *buffer;
   
    nvox = dims[0]*dims[1]*dims[2];
    buffer = (float *)malloc(sizeof(float)*nvox);
    
    /* check success of memory allocation */
    if (!buffer) {
        printf("Memory allocation error\n");
        exit(EXIT_FAILURE);
    }
   
    convert_input_type(data, buffer, nvox, datatype);
    keep_largest_cluster_float(buffer, thresh, dims, min_size, retain_above_th, conn);
    convert_output_type(data, buffer, nvox, datatype);
    
    free(buffer);
}

/**
 * fill_holes - Fill holes in a 3D volume after thresholding.
 *
 * This function is designed to process a 3D volume by filling the holes 
 * created after applying a thresholding operation.
 *
 * The process invoxves several steps:
 * 1. Conversion of input data to a floating-point buffer based on the given datatype.
 * 2. Creation of an inverted mask based on the threshold value, where values below 
 *    the threshold are set to 1 (indicating potential holes) and values above are set to 0.
 * 3. Identification and retention of the largest cluster in the inverted mask, typically 
 *    representing the background, using the `keep_largest_cluster_float` function.
 * 4. Filling the identified holes in the original data buffer by setting values corresponding 
 *    to the holes to the threshold value.
 * 5. Conversion of the processed data back to the original datatype.
 *
 * Parameters:
 *  data: Pointer to the input data buffer. This buffer is modified in-place.
 *  thresh: Threshold value used for identifying holes. Values below this threshold 
 *          are considered potential holes.
 *  dims: Array of 3 integers representing the dimensions of the 3D volume 
 *         (width, height, depth).
 *  datatype: An integer representing the datatype of the input data. This is used 
 *             to correctly interpret and manipulate the data buffer.
 *
 */
void fill_holes(void *data, double thresh, int *dims, int datatype)
{
    int i, nvox;
    float *buffer, *mask_inv;
    double voxelsize[3] = {1.0, 1.0, 1.0};
    unsigned char *mask_fill;
   
    nvox = dims[0]*dims[1]*dims[2];
    buffer = (float *)malloc(sizeof(float)*nvox);
    mask_inv = (float *)malloc(sizeof(float)*nvox);
    mask_fill = (unsigned char *)malloc(sizeof(unsigned char)*nvox);
    
    /* check success of memory allocation */
    if (!buffer || !mask_fill || !mask_inv) {
        printf("Memory allocation error\n");
        exit(EXIT_FAILURE);
    }
   
    convert_input_type(data, buffer, nvox, datatype);
    
    /* get inverted mask after thresholding */
    for (i = 0; i < nvox; ++i)
        mask_inv[i] = (buffer[i] >= thresh) ? 0.0 : 1.0;

    /* retain largest cluster of inverted mask, which should be the background */
    keep_largest_cluster_float(mask_inv, thresh, dims, 0, 1, 18);
    
    /* fill those values (=holes) that were removed by the previous keep_largest_cluster step */
    for (i = 0; i < nvox; ++i)
        mask_fill[i] = ((mask_inv[i] == 0.0) && (buffer[i] < thresh)) ? 1 : 0;
    vbdist(buffer, mask_fill, dims, NULL, 1);
    
    /* ensure a minimum filled value that is the threshold */
    for (i = 0; i < nvox; ++i)
        buffer[i] = ((mask_fill[i] == 1) && (buffer[i] < thresh)) ? thresh : buffer[i];
        
    convert_output_type(data, buffer, nvox, datatype);
    
    free(buffer);
    free(mask_inv);
    free(mask_fill);
}

/**
 * get x-gradient of 3D volume
*/
float gradientX(float *src, int i, int j, int k, int dims[3]) 
{
    if( i > 0 )
    {
        if ( i < dims[0] - 1 )
            return (src[sub2ind(i+1, j, k, dims)] - src[sub2ind(i-1, j, k, dims)])/2.0;
        else
            return src[sub2ind(i, j, k, dims)] - src[sub2ind(i-1, j, k, dims)];
    }
    else
        return src[sub2ind(i+1, j, k, dims)] - src[sub2ind(i, j, k, dims)];
}

/**
 * get y-gradient of 3D volume
*/
float gradientY(float *src, int i, int j, int k, int dims[3]) 
{
    if( j > 0 )
    {
      if ( j < dims[1] - 1 )
          return (src[sub2ind(i, j+1, k, dims)] - src[sub2ind(i, j-1, k, dims)])/2.0;
      else
          return src[sub2ind(i, j, k, dims)] - src[sub2ind(i, j-1, k, dims)];
    }
    else
        return src[sub2ind(i, j+1, k, dims)] - src[sub2ind(i, j, k, dims)];
}

/**
 * get z-gradient of 3D volume
*/
float gradientZ(float *src, int i, int j, int k, int dims[3]) 
{
    if( k > 0 )
    {
      if ( k < dims[2] - 1 )
          return (src[sub2ind(i, j, k+1, dims)] - src[sub2ind(i, j, k-1, dims)])/2.0;
      else
          return src[sub2ind(i, j, k, dims)] - src[sub2ind(i, j, k-1, dims)];
    }
    else
        return src[sub2ind(i, j, k-1, dims)] - src[sub2ind(i, j, k, dims)];
}

/**
 * estimate magnitude of local gradient of 3D volume
*/
void gradient3D_magnitude(float *src, float *grad_mag, int dims[3])
{
    int i, j, k;
    float gradx, grady, gradz;
    
    if (grad_mag == NULL)
        grad_mag = (float *)malloc(sizeof(float)*dims[0]*dims[1]*dims[2]);
    
    for (i = 0; i < dims[0]; ++i) {
        for (j = 0; j < dims[1]; ++j) {
            for (k = 0; k < dims[2]; ++k) {
                gradx = gradientX(src, i, j, k, dims);
                grady = gradientY(src, i, j, k, dims);
                gradz = gradientZ(src, i, j, k, dims);
                grad_mag[sub2ind(i, j, k, dims)] = sqrt(gradx*gradx + grady*grady + gradz*gradz);
            }
        }
    }
}
