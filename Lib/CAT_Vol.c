/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 *
 */

#include "CAT_Vol.h"
#include "CAT_Math.h"

#if defined(_WIN32) || defined(_WIN64)
#include <windows.h>
#include <process.h> /* _beginthreadex, _endthreadex */
#else
#include <pthread.h>
#endif

#define MAX_NTHREADS 4 /* Overhead is otherwise too large */

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
void ind2sub(int i, int *x, int *y, int *z, int sxy, int sx)
{
    int j = i % sxy; // Modulo to find position within a single z-plane

    *z = (int)floor((double)i / (double)sxy); // Calculate z-coordinate
    *y = (int)floor((double)j / (double)sx);  // Calculate y-coordinate
    *x = j % sx;                              // Calculate x-coordinate
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
int sub2ind(int x, int y, int z, int s[3])
{
    // Boundary handling to ensure coordinates are within array limits
    x = (x < 0) ? 0 : (x > s[0] - 1) ? s[0] - 1
                                     : x;
    y = (y < 0) ? 0 : (y > s[1] - 1) ? s[1] - 1
                                     : y;
    z = (z < 0) ? 0 : (z > s[2] - 1) ? s[2] - 1
                                     : z;

    // Calculate and return the linear index
    return z * s[0] * s[1] + y * s[0] + x;
}

/**
 * localstat_double - Calculate local statistics for a 3D float array.
 *
 * This function calculates mean, median, min, max, and standard deviation
 * within a defined distance from voxel center for each element in a 3D double array. It
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
void localstat_double(double *input, unsigned char *mask, int dims[3], int dist,
                      int stat_func, int iters, int use_euclidean_dist)
{
    double *arr;
    int i, j, k, ind, ni, x, y, z, n, it;
    double *buffer;
    int nvox = dims[0] * dims[1] * dims[2], size_kernel;

    // Check for distance parameter
    if ((dist < 1) || (dist > 10))
    {
        printf("Distance parameter should be in the range 1..10.\n");
        exit(EXIT_FAILURE);
    }

    size_kernel = (2 * dist) + 1;

    // Memory allocation
    buffer = (double *)malloc(sizeof(double) * nvox);
    arr = (double *)malloc(sizeof(double) * size_kernel * size_kernel * size_kernel);

    // Check memory allocation success
    if (!buffer || !arr)
    {
        printf("Memory allocation error\n");
        exit(EXIT_FAILURE);
    }

    // Initialize buffer to zero
    for (i = 0; i < nvox; i++)
        buffer[i] = input[i];

    // Main filter process
    for (it = 0; it < iters; it++)
    {
        for (z = 0; z < dims[2]; z++)
            for (y = 0; y < dims[1]; y++)
                for (x = 0; x < dims[0]; x++)
                {
                    ind = sub2ind(x, y, z, dims);
                    n = 0;

                    memset(arr, 0.0, size_kernel * size_kernel * size_kernel * sizeof(double));
                    // Iterate through kernel
                    for (i = -dist; i <= dist; i++)
                        for (j = -dist; j <= dist; j++)
                            for (k = -dist; k <= dist; k++)
                            {
                                ni = sub2ind(x + i, y + j, z + k, dims);

                                // Check for NaNs, Infinities
                                if (isnan(input[ni]) || !isfinite(input[ni]))
                                    continue;

                                // Check for optional mask
                                if (mask && mask[ni] == 0)
                                    continue;

                                // Check for Euclidean distance if required
                                if (use_euclidean_dist && sqrtf((float)((i * i) + (j * j) + (k * k))) > (float)dist)
                                    continue;

                                arr[n] = input[ni];
                                n++;
                            }

                    // Check for NaNs, Infinities
                    if (isnan(input[ind]) || !isfinite(input[ind]))
                        continue;

                    // Check for optional mask
                    if (mask && mask[ind] == 0)
                        continue;

                    // Calculate local statistics based on the selected function
                    switch (stat_func)
                    {
                    case F_MEAN:
                        buffer[ind] = (double)get_mean_double(arr, n, 0);
                        break;
                    case F_MEDIAN:
                        buffer[ind] = (double)get_median_double(arr, n, 0);
                        break;
                    case F_STD:
                        buffer[ind] = (double)get_std_double(arr, n, 0);
                        break;
                    case F_MIN:
                        buffer[ind] = (double)get_min_double(arr, n, 0);
                        break;
                    case F_MAX:
                        buffer[ind] = (double)get_max_double(arr, n, 0);
                        break;
                    default:
                        fprintf(stderr, "Data Function %d not handled\n", stat_func);
                        break;
                    }
                }
        // Copy results back to input array
        for (i = 0; i < nvox; i++)
            input[i] = buffer[i];
    }

    // Free allocated memory
    free(buffer);
    free(arr);
}

/**
 * \brief Public API for localstat3.
 *
 * This function is part of the CAT-Surface public library interface and is used by command-line tools.
 *
 * \param data (in/out) Parameter of localstat3.
 * \param mask (in/out) Parameter of localstat3.
 * \param dims (in/out) Parameter of localstat3.
 * \param dist (in/out) Parameter of localstat3.
 * \param stat_func (in/out) Parameter of localstat3.
 * \param iters (in/out) Parameter of localstat3.
 * \param use_euclidean_dist (in/out) Parameter of localstat3.
 * \param datatype (in/out) Parameter of localstat3.
 * \return void (no return value).
 */
void localstat3(void *data, unsigned char *mask, int dims[3], int dist,
                int stat_func, int iters, int use_euclidean_dist, int datatype)
{
    int nvox;
    double *buffer;

    // Check for distance parameter
    if ((dist < 1) || (dist > 10))
    {
        printf("Distance parameter should be in the range 1..10.\n");
        exit(EXIT_FAILURE);
    }

    nvox = dims[0] * dims[1] * dims[2];
    buffer = (double *)malloc(sizeof(double) * nvox);

    /* check success of memory allocation */
    if (!buffer)
    {
        printf("Memory allocation error\n");
        exit(EXIT_FAILURE);
    }

    convert_input_type(data, buffer, nvox, datatype);
    localstat_double(buffer, mask, dims, dist, stat_func, iters, use_euclidean_dist);
    convert_output_type(data, buffer, nvox, datatype);

    free(buffer);
}

/**
 * \brief Apply median filtering to a 3D volume.
 *
 * Applies median filtering using a 3x3x3 neighborhood with optional masking.
 * Median filtering removes noise while preserving edges better than Gaussian smoothing.
 * The computation is performed on a double-precision buffer regardless of input datatype.
 *
 * \param data      (in/out) void pointer to volume data; type given by datatype parameter
 * \param mask      (in)     optional unsigned char mask (NULL to process entire volume)
 * \param dims      (in)     {nx, ny, nz} dimension array
 * \param iters     (in)     number of median filtering iterations
 * \param datatype  (in)     data type code (DT_UINT8, DT_UINT16, DT_FLOAT32, etc.)
 */
void median3(void *data, unsigned char *mask, int dims[3], int iters, int datatype)
{
    int nvox;
    double *buffer;

    nvox = dims[0] * dims[1] * dims[2];
    buffer = (double *)malloc(sizeof(double) * nvox);

    /* check success of memory allocation */
    if (!buffer)
    {
        printf("Memory allocation error\n");
        exit(EXIT_FAILURE);
    }

    convert_input_type(data, buffer, nvox, datatype);

    /* use kernel 3x3x3 */
    localstat_double(buffer, mask, dims, 1, F_MEDIAN, iters, 0);
    convert_output_type(data, buffer, nvox, datatype);

    free(buffer);
}

/**
 * \brief Find the minimum positive value in a float array.
 *
 * Searches for the smallest positive (> 0) value in an array of floats,
 * along with its array index. Non-positive values are ignored.
 * Useful for finding the nearest positive voxel distance or similar applications.
 *
 * \param A        (in)  float array to search through
 * \param sA       (in)  size of array A
 * \param minimum  (out) pointer to store the minimum positive value found (FLT_MAX if none)
 * \param index    (out) pointer to store the index of the minimum positive value (0 if none)
 */
void pmin(float *A, int sA, float *minimum, int *index)
{
    int i;

    *minimum = FLT_MAX;
    *index = 0;

    for (i = 0; i < sA; i++)
    {
        if ((A[i] > 0.0) && (*minimum > A[i]))
        {
            *minimum = A[i];
            *index = i;
        }
    }
}

/**
 * \brief Public API for isoval.
 *
 * This function is part of the CAT-Surface public library interface and is used by command-line tools.
 *
 * \param vol (in/out) Parameter of isoval.
 * \param x (in/out) Parameter of isoval.
 * \param y (in/out) Parameter of isoval.
 * \param z (in/out) Parameter of isoval.
 * \param dims (in/out) Parameter of isoval.
 * \param nii_ptr (in/out) Parameter of isoval.
 * \return Return value of isoval.
 */
float isoval(float *vol, float x, float y, float z, int dims[3], nifti_image *nii_ptr)
{
    float seg = 0.0, n = 0.0;
    float world_coords[3] = {x, y, z}; // Define world coordinates (in mm)
    int i;

    // Convert from world to voxel space if NIfTI pointer is defined
    if (nii_ptr)
    {
        mat44 mat = nii_ptr->sto_xyz;                 // Transformation matrix
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
    float cx = floor(x + 1), cy = floor(y + 1), cz = floor(z + 1);

    // Weight factors for interpolation
    float wfx = cx - x, wfy = cy - y, wfz = cz - z;
    float wcx = x - fx, wcy = y - fy, wcz = z - fz;
    float N[8], W[8]; // Neighbours and their weights

    // Calculate value of the 8 neighbours and their distance weight
    for (i = 0; i < 8; i++)
    {
        int ix = (i & 1) ? cx : fx;
        int iy = (i & 2) ? cy : fy;
        int iz = (i & 4) ? cz : fz;

        N[i] = vol[sub2ind(ix, iy, iz, dims)];
        W[i] = ((i & 1) ? wcx : wfx) * ((i & 2) ? wcy : wfy) * ((i & 4) ? wcz : wfz);

        // Perform interpolation using neighbours and weights
        if (!isnan(N[i]) && isfinite(N[i]))
        {
            seg += N[i] * W[i];
            n += W[i];
        }
    }

    // Return interpolated value or NaN if unable to interpolate
    return n > 0.0 ? seg / n : FNAN;
}

/* ------------------------ Row-pass worker ------------------------- */
#if defined(_WIN32) || defined(_WIN64)
static unsigned int __stdcall
#else
static void *
#endif
conv_row_worker(void *p)
{
    conv_args_row a = *(conv_args_row *)p;

    /* per-thread temporary buffer for a whole row */
    float *tbuf = (float *)malloc((size_t)a.xdim * sizeof(float));
    int x, y, k;

    if (!tbuf)
    {
#if defined(_WIN32) || defined(_WIN64)
        _endthreadex(0);
        return 0;
#else
        return 0;
#endif
    }

    for (y = a.ini; y < a.fin; ++y)
    {
        /* copy row y to buffer & sanitise NaNs/Infs */
        for (x = 0; x < a.xdim; ++x)
        {
            tbuf[x] = a.out[x + y * a.xdim];
            if (!isfinite(tbuf[x]))
                tbuf[x] = 0.0f;
        }
        /* horizontal (x) convolution */
        for (x = 0; x < a.xdim; ++x)
        {
            double sum1 = 0.0;
            int fstart = ((x - a.xoff >= a.xdim) ? x - a.xdim - a.xoff + 1 : 0);
            int fend = ((x - (a.xoff + a.fxdim) < 0) ? x - a.xoff + 1 : a.fxdim);
            for (k = fstart; k < fend; ++k)
                sum1 += (double)tbuf[x - a.xoff - k] * a.filtx[k];
            a.out[x + y * a.xdim] = (float)sum1;
        }
    }

    free(tbuf);

#if defined(_WIN32) || defined(_WIN64)
    _endthreadex(0);
    return 0;
#else
    return 0;
#endif
}

/* ------------------------ Column-pass worker ---------------------- */
#if defined(_WIN32) || defined(_WIN64)
static unsigned int __stdcall
#else
static void *
#endif
conv_col_worker(void *p)
{
    conv_args_col a = *(conv_args_col *)p;

    /* per-thread temporary buffer for a whole column */
    float *tbuf = (float *)malloc((size_t)a.ydim * sizeof(float));
    int x, y, k;

    if (!tbuf)
    {
#if defined(_WIN32) || defined(_WIN64)
        _endthreadex(0);
        return 0;
#else
        return 0;
#endif
    }

    for (x = a.ini; x < a.fin; ++x)
    {
        /* copy column x to buffer */
        for (y = 0; y < a.ydim; ++y)
            tbuf[y] = a.out[x + y * a.xdim];

        /* vertical (y) convolution writing back to out */
        for (y = 0; y < a.ydim; ++y)
        {
            double sum1 = 0.0;
            int fstart = ((y - a.yoff >= a.ydim) ? y - a.ydim - a.yoff + 1 : 0);
            int fend = ((y - (a.yoff + a.fydim) < 0) ? y - a.yoff + 1 : a.fydim);
            for (k = fstart; k < fend; ++k)
                sum1 += (double)tbuf[y - a.yoff - k] * a.filty[k];
            a.out[y * a.xdim + x] = (float)sum1;
        }
    }

    free(tbuf);

#if defined(_WIN32) || defined(_WIN64)
    _endthreadex(0);
    return 0;
#else
    return 0;
#endif
}

/* ===================== Multithreaded convxy_float ===================== */
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
 *
 * In multithreaded mode, per-thread temporary buffers are used internally.
 * The 'buff' argument is only used in the single-thread fallback.
 */
static void convxy_float(float *out, int xdim, int ydim,
                         double *filtx, double *filty,
                         int fxdim, int fydim,
                         int xoff, int yoff,
                         float *buff)
{
    /* choose thread counts similarly to CAT_Sanlm.c style (cap at 8) */
    int Nthreads_row = (ydim < MAX_NTHREADS) ? ydim : MAX_NTHREADS;
    int Nthreads_col = (xdim < MAX_NTHREADS) ? xdim : MAX_NTHREADS;
    int t;
    if (Nthreads_row < 1)
        Nthreads_row = 1;
    if (Nthreads_col < 1)
        Nthreads_col = 1;

    /* If only one thread would be used overall, keep the original serial path. */
    if (Nthreads_row == 1 && Nthreads_col == 1)
    {
        int x, y, k;
        for (y = 0; y < ydim; y++)
        {
            for (x = 0; x < xdim; x++)
            {
                buff[x] = out[x + y * xdim];
                if (!isfinite(buff[x]))
                    buff[x] = 0.0f;
            }
            for (x = 0; x < xdim; x++)
            {
                double sum1 = 0.0;
                int fstart = ((x - xoff >= xdim) ? x - xdim - xoff + 1 : 0);
                int fend = ((x - (xoff + fxdim) < 0) ? x - xoff + 1 : fxdim);
                for (k = fstart; k < fend; k++)
                    sum1 += (double)buff[x - xoff - k] * filtx[k];
                out[x + y * xdim] = (float)sum1;
            }
        }
        for (x = 0; x < xdim; x++)
        {
            for (y = 0; y < ydim; y++)
                buff[y] = out[x + y * xdim];

            for (y = 0; y < ydim; y++)
            {
                double sum1 = 0.0;
                int fstart = ((y - yoff >= ydim) ? y - ydim - yoff + 1 : 0);
                int fend = ((y - (yoff + fydim) < 0) ? y - yoff + 1 : fydim);
                for (k = fstart; k < fend; k++)
                    sum1 += (double)buff[y - yoff - k] * filty[k];
                out[y * xdim + x] = (float)sum1;
            }
        }
        return;
    }

    /* -------------------- Pass 1: parallel across rows -------------------- */
#if defined(_WIN32) || defined(_WIN64)
    HANDLE *RowThreads = (HANDLE *)malloc((size_t)Nthreads_row * sizeof(HANDLE));
    conv_args_row *RowArgs = (conv_args_row *)malloc((size_t)Nthreads_row * sizeof(conv_args_row));

    for (t = 0; t < Nthreads_row; ++t)
    {
        int ini = (t * ydim) / Nthreads_row;
        int fin = ((t + 1) * ydim) / Nthreads_row;

        RowArgs[t].out = out;
        RowArgs[t].xdim = xdim;
        RowArgs[t].ydim = ydim;
        RowArgs[t].filtx = filtx;
        RowArgs[t].filty = filty;
        RowArgs[t].fxdim = fxdim;
        RowArgs[t].fydim = fydim;
        RowArgs[t].xoff = xoff;
        RowArgs[t].yoff = yoff;
        RowArgs[t].ini = ini;
        RowArgs[t].fin = fin;

        RowThreads[t] = (HANDLE)_beginthreadex(NULL, 0, &conv_row_worker, &RowArgs[t], 0, NULL);
    }
    for (t = 0; t < Nthreads_row; ++t)
        WaitForSingleObject(RowThreads[t], INFINITE);
    for (t = 0; t < Nthreads_row; ++t)
        CloseHandle(RowThreads[t]);
    free(RowThreads);
    free(RowArgs);
#else
    pthread_t *RowThreads = (pthread_t *)calloc((size_t)Nthreads_row, sizeof(pthread_t));
    conv_args_row *RowArgs = (conv_args_row *)calloc((size_t)Nthreads_row, sizeof(conv_args_row));
    for (t = 0; t < Nthreads_row; ++t)
    {
        int ini = (t * ydim) / Nthreads_row;
        int fin = ((t + 1) * ydim) / Nthreads_row;

        RowArgs[t].out = out;
        RowArgs[t].xdim = xdim;
        RowArgs[t].ydim = ydim;
        RowArgs[t].filtx = filtx;
        RowArgs[t].filty = filty;
        RowArgs[t].fxdim = fxdim;
        RowArgs[t].fydim = fydim;
        RowArgs[t].xoff = xoff;
        RowArgs[t].yoff = yoff;
        RowArgs[t].ini = ini;
        RowArgs[t].fin = fin;

        pthread_create(&RowThreads[t], NULL, conv_row_worker, &RowArgs[t]);
    }
    for (t = 0; t < Nthreads_row; ++t)
        pthread_join(RowThreads[t], NULL);
    free(RowThreads);
    free(RowArgs);
#endif

    /* ------------------- Pass 2: parallel across columns ------------------ */
#if defined(_WIN32) || defined(_WIN64)
    HANDLE *ColThreads = (HANDLE *)malloc((size_t)Nthreads_col * sizeof(HANDLE));
    conv_args_col *ColArgs = (conv_args_col *)malloc((size_t)Nthreads_col * sizeof(conv_args_col));

    for (t = 0; t < Nthreads_col; ++t)
    {
        int ini = (t * xdim) / Nthreads_col;
        int fin = ((t + 1) * xdim) / Nthreads_col;

        ColArgs[t].out = out;
        ColArgs[t].xdim = xdim;
        ColArgs[t].ydim = ydim;
        ColArgs[t].filtx = filtx;
        ColArgs[t].filty = filty;
        ColArgs[t].fxdim = fxdim;
        ColArgs[t].fydim = fydim;
        ColArgs[t].xoff = xoff;
        ColArgs[t].yoff = yoff;
        ColArgs[t].ini = ini;
        ColArgs[t].fin = fin;

        ColThreads[t] = (HANDLE)_beginthreadex(NULL, 0, &conv_col_worker, &ColArgs[t], 0, NULL);
    }
    for (t = 0; t < Nthreads_col; ++t)
        WaitForSingleObject(ColThreads[t], INFINITE);
    for (t = 0; t < Nthreads_col; ++t)
        CloseHandle(ColThreads[t]);
    free(ColThreads);
    free(ColArgs);
#else
    pthread_t *ColThreads = (pthread_t *)calloc((size_t)Nthreads_col, sizeof(pthread_t));
    conv_args_col *ColArgs = (conv_args_col *)calloc((size_t)Nthreads_col, sizeof(conv_args_col));
    for (t = 0; t < Nthreads_col; ++t)
    {
        int ini = (t * xdim) / Nthreads_col;
        int fin = ((t + 1) * xdim) / Nthreads_col;

        ColArgs[t].out = out;
        ColArgs[t].xdim = xdim;
        ColArgs[t].ydim = ydim;
        ColArgs[t].filtx = filtx;
        ColArgs[t].filty = filty;
        ColArgs[t].fxdim = fxdim;
        ColArgs[t].fydim = fydim;
        ColArgs[t].xoff = xoff;
        ColArgs[t].yoff = yoff;
        ColArgs[t].ini = ini;
        ColArgs[t].fin = fin;

        pthread_create(&ColThreads[t], NULL, conv_col_worker, &ColArgs[t]);
    }
    for (t = 0; t < Nthreads_col; ++t)
        pthread_join(ColThreads[t], NULL);
    free(ColThreads);
    free(ColArgs);
#endif
}

/* ---------- Stage 1: compute 2D convxy for each z-slice in parallel ---------- */

/* forward decl of convxy_float already in the same file */
static void convxy_float(float *out, int xdim, int ydim,
                         double *filtx, double *filty,
                         int fxdim, int fydim,
                         int xoff, int yoff,
                         float *buff);

#if defined(_WIN32) || defined(_WIN64)
static unsigned int __stdcall
#else
static void *
#endif
convxyz_stage1_worker(void *p)
{
    convxyz_s1_args_t a = *(convxyz_s1_args_t *)p;
    const int xy = a.xdim * a.ydim;

    /* Per-thread temp buffer used by convxy_float in serial path */
    float *buff = (float *)malloc((size_t)((a.ydim > a.xdim) ? a.ydim : a.xdim) * sizeof(float));
    int i, z;

    if (!buff)
    {
#if defined(_WIN32) || defined(_WIN64)
        _endthreadex(0);
        return 0;
#else
        return 0;
#endif
    }

    for (z = a.z_ini; z < a.z_fin; ++z)
    {
        float *dst = a.convxy_vol + (size_t)z * xy;
        const float *src = a.iVol + (size_t)z * xy;

        /* copy slice z into dst; convxy_float works in-place on dst */
        for (i = 0; i < xy; ++i)
            dst[i] = src[i];

        /* 2D convolution (may itself be multi-threaded as implemented) */
        convxy_float(dst, a.xdim, a.ydim,
                     (double *)a.filtx, (double *)a.filty,
                     a.fxdim, a.fydim, a.xoff, a.yoff,
                     buff);
    }

    free(buff);

#if defined(_WIN32) || defined(_WIN64)
    _endthreadex(0);
    return 0;
#else
    return 0;
#endif
}

/* ---------- Stage 2: z-direction 1D convolution per output slice in parallel ---------- */

#if defined(_WIN32) || defined(_WIN64)
static unsigned int __stdcall
#else
static void *
#endif
convxyz_stage2_worker(void *p)
{
    convxyz_s2_args_t a = *(convxyz_s2_args_t *)p;
    const int xy = a.xdim * a.ydim;
    int k, idx, z_out;

    for (z_out = a.z_out_ini; z_out < a.z_out_fin; ++z_out)
    {
        /* Map to loop variable 'z' from original code:
           z = z_out + fzdim + zoff - 1  */
        const int z = z_out + a.fzdim + a.zoff - 1;

        /* original bounds:
           fstart = ((z >= zdim) ? z - zdim + 1 : 0);
           fend   = ((z - fzdim < 0) ? z + 1 : fzdim); */
        const int fstart = (z >= a.zdim) ? (z - a.zdim + 1) : 0;
        const int fend = (z - a.fzdim < 0) ? (z + 1) : a.fzdim;

        /* normaliser sum2 */
        double sum2 = 0.0;
        for (k = fstart; k < fend; ++k)
            sum2 += a.filtz[k];

        float *obuf = a.oVol + (size_t)z_out * xy;

        if (sum2 != 0.0)
        {
            for (idx = 0; idx < xy; ++idx)
            {
                double sum1 = 0.0;
                for (k = fstart; k < fend; ++k)
                {
                    const int z_src = z - k; /* guaranteed 0..zdim-1 */
                    sum1 += a.filtz[k] * (double)a.convxy_vol[(size_t)z_src * xy + idx];
                }
                obuf[idx] = (float)(sum1 / sum2);
            }
        }
        else
        {
            for (idx = 0; idx < xy; ++idx)
                obuf[idx] = 0.0f;
        }
    }

#if defined(_WIN32) || defined(_WIN64)
    _endthreadex(0);
    return 0;
#else
    return 0;
#endif
}

/* ========================= Multithreaded convxyz_float ========================= */
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
 * Two-parallel-stage implementation (slice 2D conv, then z-accumulation).
 *  No locks needed: each thread works on disjoint z (or z_out) ranges.
 *  Memory: uses one temporary volume convxy_vol[z][y][x].
 */
int convxyz_float(float *iVol, double *filtx, double *filty, double *filtz,
                  int fxdim, int fydim, int fzdim,
                  int xoff, int yoff, int zoff, float *oVol, int dims[3])
{
    const int xdim = dims[0];
    const int ydim = dims[1];
    const int zdim = dims[2];
    const int xy = xdim * ydim;
    int i;

    if (xdim <= 0 || ydim <= 0 || zdim <= 0 || fxdim <= 0 || fydim <= 0 || fzdim <= 0)
        return 0;

    /* temp storage for all convxy slices */
    float *convxy_vol = (float *)malloc((size_t)zdim * xy * sizeof(float));
    if (!convxy_vol)
    {
        fprintf(stderr, "Memory allocation error (convxy_vol)\n");
        exit(EXIT_FAILURE);
    }

    /* ---------------- Stage 1: parallel convxy over z ---------------- */
    {
        int Nthreads = (zdim < MAX_NTHREADS) ? zdim : MAX_NTHREADS;
        if (Nthreads < 1)
            Nthreads = 1;

#if defined(_WIN32) || defined(_WIN64)
        HANDLE *ThreadList = (HANDLE *)malloc((size_t)Nthreads * sizeof(HANDLE));
        convxyz_s1_args_t *Args = (convxyz_s1_args_t *)malloc((size_t)Nthreads * sizeof(convxyz_s1_args_t));
        if (!ThreadList || !Args)
        {
            fprintf(stderr, "Memory allocation error (threads stage1)\n");
            exit(EXIT_FAILURE);
        }
        for (i = 0; i < Nthreads; ++i)
        {
            int ini = (i * zdim) / Nthreads;
            int fin = ((i + 1) * zdim) / Nthreads;

            Args[i].iVol = iVol;
            Args[i].xdim = xdim;
            Args[i].ydim = ydim;
            Args[i].zdim = zdim;
            Args[i].filtx = filtx;
            Args[i].filty = filty;
            Args[i].fxdim = fxdim;
            Args[i].fydim = fydim;
            Args[i].xoff = xoff;
            Args[i].yoff = yoff;
            Args[i].convxy_vol = convxy_vol;
            Args[i].z_ini = ini;
            Args[i].z_fin = fin;

            ThreadList[i] = (HANDLE)_beginthreadex(NULL, 0, &convxyz_stage1_worker, &Args[i], 0, NULL);
        }
        for (i = 0; i < Nthreads; ++i)
            WaitForSingleObject(ThreadList[i], INFINITE);
        for (i = 0; i < Nthreads; ++i)
            CloseHandle(ThreadList[i]);
        free(ThreadList);
        free(Args);
#else
        pthread_t *ThreadList = (pthread_t *)calloc((size_t)Nthreads, sizeof(pthread_t));
        convxyz_s1_args_t *Args = (convxyz_s1_args_t *)calloc((size_t)Nthreads, sizeof(convxyz_s1_args_t));
        if (!ThreadList || !Args)
        {
            fprintf(stderr, "Memory allocation error (threads stage1)\n");
            exit(EXIT_FAILURE);
        }
        for (i = 0; i < Nthreads; ++i)
        {
            int ini = (i * zdim) / Nthreads;
            int fin = ((i + 1) * zdim) / Nthreads;

            Args[i].iVol = iVol;
            Args[i].xdim = xdim;
            Args[i].ydim = ydim;
            Args[i].zdim = zdim;
            Args[i].filtx = filtx;
            Args[i].filty = filty;
            Args[i].fxdim = fxdim;
            Args[i].fydim = fydim;
            Args[i].xoff = xoff;
            Args[i].yoff = yoff;
            Args[i].convxy_vol = convxy_vol;
            Args[i].z_ini = ini;
            Args[i].z_fin = fin;

            pthread_create(&ThreadList[i], NULL, convxyz_stage1_worker, &Args[i]);
        }
        for (i = 0; i < Nthreads; ++i)
            pthread_join(ThreadList[i], NULL);
        free(ThreadList);
        free(Args);
#endif
    }

    /* ---------------- Stage 2: parallel z-accumulation ---------------- */
    {
        int Nthreads = (zdim < MAX_NTHREADS) ? zdim : MAX_NTHREADS;
        if (Nthreads < 1)
            Nthreads = 1;

#if defined(_WIN32) || defined(_WIN64)
        HANDLE *ThreadList = (HANDLE *)malloc((size_t)Nthreads * sizeof(HANDLE));
        convxyz_s2_args_t *Args = (convxyz_s2_args_t *)malloc((size_t)Nthreads * sizeof(convxyz_s2_args_t));
        if (!ThreadList || !Args)
        {
            fprintf(stderr, "Memory allocation error (threads stage2)\n");
            exit(EXIT_FAILURE);
        }
        for (i = 0; i < Nthreads; ++i)
        {
            int ini = (i * zdim) / Nthreads;
            int fin = ((i + 1) * zdim) / Nthreads;

            Args[i].convxy_vol = convxy_vol;
            Args[i].oVol = oVol;
            Args[i].xdim = xdim;
            Args[i].ydim = ydim;
            Args[i].zdim = zdim;
            Args[i].filtz = filtz;
            Args[i].fzdim = fzdim;
            Args[i].zoff = zoff;
            Args[i].z_out_ini = ini;
            Args[i].z_out_fin = fin;

            ThreadList[i] = (HANDLE)_beginthreadex(NULL, 0, &convxyz_stage2_worker, &Args[i], 0, NULL);
        }
        for (i = 0; i < Nthreads; ++i)
            WaitForSingleObject(ThreadList[i], INFINITE);
        for (i = 0; i < Nthreads; ++i)
            CloseHandle(ThreadList[i]);
        free(ThreadList);
        free(Args);
#else
        pthread_t *ThreadList = (pthread_t *)calloc((size_t)Nthreads, sizeof(pthread_t));
        convxyz_s2_args_t *Args = (convxyz_s2_args_t *)calloc((size_t)Nthreads, sizeof(convxyz_s2_args_t));
        if (!ThreadList || !Args)
        {
            fprintf(stderr, "Memory allocation error (threads stage2)\n");
            exit(EXIT_FAILURE);
        }
        for (i = 0; i < Nthreads; ++i)
        {
            int ini = (i * zdim) / Nthreads;
            int fin = ((i + 1) * zdim) / Nthreads;

            Args[i].convxy_vol = convxy_vol;
            Args[i].oVol = oVol;
            Args[i].xdim = xdim;
            Args[i].ydim = ydim;
            Args[i].zdim = zdim;
            Args[i].filtz = filtz;
            Args[i].fzdim = fzdim;
            Args[i].zoff = zoff;
            Args[i].z_out_ini = ini;
            Args[i].z_out_fin = fin;

            pthread_create(&ThreadList[i], NULL, convxyz_stage2_worker, &Args[i]);
        }
        for (i = 0; i < Nthreads; ++i)
            pthread_join(ThreadList[i], NULL);
        free(ThreadList);
        free(Args);
#endif
    }

    free(convxy_vol);
    return 0;
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
int convxyz_uint8(unsigned char *iVol, double *filtx, double *filty, double *filtz,
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

    tmp = (float *)malloc(sizeof(float) * xdim * ydim * fzdim);
    buff = (float *)malloc(sizeof(float) * ((ydim > xdim) ? ydim : xdim));
    sortedv = (float **)malloc(sizeof(float *) * fzdim);

    if (!tmp || !buff || !sortedv)
    {
        fprintf(stderr, "Memory allocation error\n");
        exit(EXIT_FAILURE);
    }

    startz = ((fzdim + zoff - 1 < 0) ? fzdim + zoff - 1 : 0);
    endz = zdim + fzdim + zoff - 1;

    for (z = startz; z < endz; z++)
    {
        double sum2 = 0.0;

        if (z >= 0 && z < zdim)
        {
            for (y = 0; y < ydim; y++)
                for (x = 0; x < xdim; x++)
                    tmp[((z % fzdim) * xdim * ydim) + (y * xdim) + x] = (float)iVol[(z * xdim * ydim) + (y * xdim) + x];
            convxy_float(tmp + ((z % fzdim) * xdim * ydim), xdim, ydim,
                         filtx, filty, fxdim, fydim, xoff, yoff, buff);
        }
        if (z - fzdim - zoff + 1 >= 0 && z - fzdim - zoff + 1 < zdim)
        {
            fstart = ((z >= zdim) ? z - zdim + 1 : 0);
            fend = ((z - fzdim < 0) ? z + 1 : fzdim);

            for (k = 0; k < fzdim; k++)
            {
                int z1 = (((z - k) % fzdim) + fzdim) % fzdim;
                sortedv[k] = &(tmp[z1 * xdim * ydim]);
            }

            for (k = fstart, sum2 = 0.0; k < fend; k++)
                sum2 += filtz[k];

            obuf = oVol;
            obuf = &obuf[(z - fzdim - zoff + 1) * ydim * xdim];
            if (sum2)
            {
                for (xy = 0; xy < xdim * ydim; xy++)
                {
                    float sum1 = 0.0;
                    for (k = fstart; k < fend; k++)
                        sum1 += filtz[k] * sortedv[k][xy];
                    tmp2 = sum1 / sum2;
                    if (tmp2 < 0.0)
                        tmp2 = 0.0;
                    else if (tmp2 > 255.0)
                        tmp2 = 255.0;
                    obuf[xy] = (unsigned char)roundf(tmp2);
                }
            }
            else
                for (xy = 0; xy < xdim * ydim; xy++)
                    obuf[xy] = 0;
        }
    }
    free(tmp);
    free(buff);
    free(sortedv);
    return (0);
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
void smooth_float(float *vol, int dims[3], double voxelsize[3], double fwhm[3], int use_mask)
{
    int i;
    double xsum, ysum, zsum;
    double *x, *y, *z, s[3];
    int xyz[3], nvox, sum_mask;
    float *mask;
    unsigned char *mask2;

    nvox = dims[0] * dims[1] * dims[2];

    // Calculate the standard deviation and kernel size based on FWHM and voxel size
    for (i = 0; i < 3; i++)
    {
        s[i] = fwhm[i] / voxelsize[i];
        if (s[i] < 1.0)
            s[i] = 1.0;
        s[i] /= sqrt(8.0 * log(2.0));
        xyz[i] = (int)round(6.0 * s[i]);
    }

    // Memory allocation for Gaussian kernels
    x = (double *)malloc(sizeof(double) * ((2 * xyz[0]) + 1));
    y = (double *)malloc(sizeof(double) * ((2 * xyz[1]) + 1));
    z = (double *)malloc(sizeof(double) * ((2 * xyz[2]) + 1));

    if (!x || !y || !z)
    {
        fprintf(stderr, "Memory allocation error\n");
        exit(EXIT_FAILURE);
    }

    // Initialize and build the mask for masked smoothing
    if (use_mask)
    {
        mask = (float *)malloc(sizeof(float) * nvox);
        mask2 = (unsigned char *)malloc(sizeof(unsigned char) * nvox);

        if (!mask || !mask2)
        {
            fprintf(stderr, "Memory allocation error\n");
            exit(EXIT_FAILURE);
        }
        sum_mask = 0;
        // Build the mask based on the volume data
        for (i = 0; i < nvox; i++)
        {
            if (vol[i] == 0.0)
            {
                mask[i] = 0.0;
                mask2[i] = 0;
            }
            else
            {
                mask[i] = 1.0;
                mask2[i] = 1;
                sum_mask++;
            }
        }
    }

    // Build the Gaussian kernel
    for (i = -xyz[0]; i <= xyz[0]; i++)
        x[i + xyz[0]] = exp(-pow((double)i, 2) / (2.0 * pow(s[0], 2)));
    for (i = -xyz[1]; i <= xyz[1]; i++)
        y[i + xyz[1]] = exp(-pow((double)i, 2) / (2.0 * pow(s[1], 2)));
    for (i = -xyz[2]; i <= xyz[2]; i++)
        z[i + xyz[2]] = exp(-pow((double)i, 2) / (2.0 * pow(s[2], 2)));

    // Normalize the Gaussian kernel
    xsum = ysum = zsum = 0.0;
    for (i = 0; i < ((2 * xyz[0]) + 1); i++)
        xsum += x[i];
    for (i = 0; i < ((2 * xyz[1]) + 1); i++)
        ysum += y[i];
    for (i = 0; i < ((2 * xyz[2]) + 1); i++)
        zsum += z[i];
    for (i = 0; i < ((2 * xyz[0]) + 1); i++)
        x[i] /= xsum;
    for (i = 0; i < ((2 * xyz[1]) + 1); i++)
        y[i] /= ysum;
    for (i = 0; i < ((2 * xyz[2]) + 1); i++)
        z[i] /= zsum;

    // Apply convolution with the Gaussian kernel
    convxyz_float(vol, x, y, z, ((2 * xyz[0]) + 1), ((2 * xyz[1]) + 1), ((2 * xyz[2]) + 1), -xyz[0], -xyz[1], -xyz[2], vol, dims);

    // Apply masked smoothing if selected
    if (use_mask)
    {
        if (sum_mask > 0)
        {
            convxyz_float(mask, x, y, z, ((2 * xyz[0]) + 1), ((2 * xyz[1]) + 1), ((2 * xyz[2]) + 1), -xyz[0], -xyz[1], -xyz[2], mask, dims);
            for (i = 0; i < nvox; i++)
            {
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
 * \brief Smooth a 3D volume with a Gaussian filter.
 *
 * Applies Gaussian smoothing to a 3D volume with specified full-width half-maximum (FWHM).
 * Smoothing is applied independently along x, y, and z axes using separable convolution.
 * Handles automatic data type conversion for any supported datatype.
 *
 * \param data       (in/out) void pointer to volume data; type given by datatype parameter
 * \param dims       (in)     {nx, ny, nz} dimension array
 * \param voxelsize  (in)     voxel spacing in mm; used to scale FWHM to physical units
 * \param fwhm       (in)     {fwhm_x, fwhm_y, fwhm_z} smoothing kernel FWHM in mm
 * \param use_mask   (in)     unused/reserved for compatibility (pass 0)
 * \param datatype   (in)     data type code (DT_UINT8, DT_UINT16, DT_FLOAT32, etc.)
 */
void smooth3(void *data, int dims[3], double voxelsize[3], double fwhm[3], int use_mask, int datatype)
{
    int nvox;
    float *buffer;

    nvox = dims[0] * dims[1] * dims[2];
    buffer = (float *)malloc(sizeof(float) * nvox);

    /* check success of memory allocation */
    if (!buffer)
    {
        printf("Memory allocation error\n");
        exit(EXIT_FAILURE);
    }

    convert_input_type_float(data, buffer, nvox, datatype);
    smooth_float(buffer, dims, voxelsize, fwhm, use_mask);
    convert_output_type_float(data, buffer, nvox, datatype);

    free(buffer);
}

/**
 * \brief Public API for euclidean_distance.
 *
 * This function is part of the CAT-Surface public library interface and is used by command-line tools.
 *
 * \param V (in/out) Parameter of euclidean_distance.
 * \param M (in/out) Parameter of euclidean_distance.
 * \param dims (in/out) Parameter of euclidean_distance.
 * \param voxelsize (in/out) Parameter of euclidean_distance.
 * \param replace (in/out) Parameter of euclidean_distance.
 * \return void (no return value).
 */
void euclidean_distance(float *V, unsigned char *M, int dims[3], double *voxelsize, int replace)
{

    /* main information about input data (size, dimensions, ...) */
    const int nvox = dims[0] * dims[1] * dims[2];
    const int x = dims[0];
    const int y = dims[1];
    const int xy = x * y;

    float s1, s2, s3;

    // use default voxel size of 1mm to get distance in voxel-space
    if (!voxelsize)
    {
        s1 = s2 = s3 = 1.0;
    }
    else
    {
        s1 = (float)fabs(voxelsize[0]);
        s2 = (float)fabs(voxelsize[1]);
        s3 = (float)fabs(voxelsize[2]);
    }
    const float s12 = (float)sqrt((double)s1 * s1 + s2 * s2);    /* xy - voxel size */
    const float s13 = (float)sqrt((double)s1 * s1 + s3 * s3);    /* xz - voxel size */
    const float s23 = (float)sqrt((double)s2 * s2 + s3 * s3);    /* yz - voxel size */
    const float s123 = (float)sqrt((double)s12 * s12 + s3 * s3); /* xyz - voxel size */

    /* indices of the neighbour Ni (index distance) and euclidean distance NW */
    const int NI[14] = {0, -1, -x + 1, -x, -x - 1, -xy + 1, -xy, -xy - 1, -xy + x + 1, -xy + x, -xy + x - 1, -xy - x + 1, -xy - x, -xy - x - 1};
    const float ND[14] = {0.0, s1, s12, s2, s12, s13, s3, s13, s123, s23, s123, s123, s23, s123};
    enum
    {
        sN = (int)(sizeof(NI) / sizeof(NI[0]))
    }; // Number of neighbours
    float DN[sN];
    float DNm = FLT_MAX;
    int i, n, ni, DNi, init_M = 0;
    int u, v, w, nu, nv, nw;

    /* data */
    float *D, *buffer;
    unsigned int *I;

    I = (unsigned int *)malloc(sizeof(unsigned int) * nvox);
    D = (float *)malloc(sizeof(float) * nvox);

    /* save original input in buffer if we want to replace values */
    if (replace > 0)
    {
        buffer = (float *)malloc(sizeof(float) * nvox);
        memcpy(buffer, V, nvox * sizeof(float));
    }

    if (!D || !I || ((replace > 0) && !buffer))
    {
        fprintf(stderr, "Memory allocation error\n");
        exit(EXIT_FAILURE);
    }

    /* Initialize mask with ones if not defined */
    if (!M)
    {
        M = (unsigned char *)malloc(sizeof(unsigned char) * nvox);
        if (!M)
        {
            fprintf(stderr, "Memory allocation error\n");
            exit(EXIT_FAILURE);
        }
        init_M = 1;
        for (i = 0; i < nvox; i++)
            M[i] = 1;
    }

    /* initialisation of D and I */
    for (i = 0; i < nvox; i++)
    {
        /* check for zero or small values that should be filled */
        if ((V[i] <= 0.000001f) || isnan(V[i]))
            D[i] = FLT_MAX;
        else
            D[i] = 0.0;
        I[i] = (unsigned int)i;
    }

    /* forward direction that consider all points smaller than i */
    for (i = 0; i < nvox; i++)
    {
        if ((D[i] > 0) && (M[i] > 0))
        {
            ind2sub(i, &u, &v, &w, xy, x);

            /* read neighbour values */
            for (n = 0; n < sN; n++)
            {
                ni = i + NI[n];
                ind2sub(ni, &nu, &nv, &nw, xy, x);
                if ((ni < 0) || (ni >= nvox) || (abs(nu - u) > 1) || (abs(nv - v) > 1) || (abs(nw - w) > 1))
                    ni = i;
                DN[n] = D[ni] + ND[n];
            }

            /* find minimum distance within the neighbourhood */
            pmin(DN, sN, &DNm, &DNi);

            /* update values */
            if (DNi > 0)
            {
                I[i] = (unsigned int)I[i + NI[DNi]];
                D[i] = DNm;
                ind2sub((int)I[i], &nu, &nv, &nw, xy, x);
                D[i] = (float)sqrt(pow((double)(u - nu) * s1, 2) + pow((double)(v - nv) * s2, 2) + pow((double)(w - nw) * s3, 2));
            }
        }
    }

    /* backward direction that consider all points larger than i */
    for (i = nvox - 1; i >= 0; i--)
    {
        if ((D[i] > 0) && (M[i] > 0))
        {
            ind2sub(i, &u, &v, &w, xy, x);

            /* read neighbour values */
            for (n = 0; n < sN; n++)
            {
                ni = i - NI[n];
                ind2sub(ni, &nu, &nv, &nw, xy, x);
                if ((ni < 0) || (ni >= nvox) || (abs(nu - u) > 1) || (abs(nv - v) > 1) || (abs(nw - w) > 1))
                    ni = i;
                DN[n] = D[ni] + ND[n];
            }

            /* find minimum distance within the neighbourhood */
            pmin(DN, sN, &DNm, &DNi);

            /* update values */
            if (DNi > 0)
            {
                I[i] = (unsigned int)I[i - NI[DNi]];
                D[i] = DNm;
                ind2sub((int)I[i], &nu, &nv, &nw, xy, x);
                D[i] = (float)sqrt(pow((double)(u - nu) * s1, 2) + pow((double)(v - nv) * s2, 2) + pow((double)(w - nw) * s3, 2));
            }
        }
    }

    /* finally return output to original variable V */
    for (i = 0; i < nvox; i++)
        V[i] = ((M[i] == 0) || (D[i] == FLT_MAX)) ? 0.0 : D[i];

    /* finally replace values inside mask */
    if (replace > 0)
    {
        for (i = 0; i < nvox; ++i)
            V[i] = buffer[I[i]];
        free(buffer);
    }

    free(D);
    free(I);
    if (init_M)
        free(M);
}

/**
 * \brief Intensity-limited region growing with distance- and intensity-weighted path cost.
 *
 * Grows labels from an initial seed map into unlabeled voxels while enforcing a
 * monotonic intensity constraint:
 *
 *   intensity(neighbour) <= intensity(current) + limit
 *
 * The propagation cost combines geometric neighbour distance and an intensity term:
 *
 *   cost = w_dist * neighbour_distance + w_int * intensity_cost
 *
 * with weights `dd = {w_dist, w_int}`. The behaviour matches the historical
 * `cat_vol_downcut` MEX implementation, including support for a legacy default
 * intensity cost when `dd == NULL`.
 *
 *
 *
 * \param labels     (in/out) float seed label map; 0 marks unlabeled voxels
 * \param intensity  (in)     float intensity image used for monotonic growth constraint
 * \param dist       (out)    float distance/path-cost map (may be NULL if not needed)
 * \param dims       (in)     volume dimensions {nx, ny, nz}
 * \param limit      (in)     intensity limit for neighbour acceptance
 * \param voxelsize  (in)     voxel spacing {sx, sy, sz}; NULL -> {1,1,1}
 * \param dd         (in)     distance/intensity weights {w_dist, w_int}; NULL -> {0.1, 10}
 */
void downcut_float(float *labels, const float *intensity, float *dist,
                   int dims[3], double limit, double voxelsize[3], double dd[2])
{
    static const int DX[26] = {
        1, -1, 0, 0, 0, 0,
        -1, 1, -1, 1,
        -1, 1, -1, 1,
        0, 0, 0, 0,
        -1, 1, -1, 1,
        -1, 1, -1, 1};
    static const int DY[26] = {
        0, 0, 1, -1, 0, 0,
        -1, -1, 1, 1,
        0, 0, 0, 0,
        -1, 1, -1, 1,
        -1, -1, 1, 1,
        -1, -1, 1, 1};
    static const int DZ[26] = {
        0, 0, 0, 0, 1, -1,
        0, 0, 0, 0,
        -1, -1, 1, 1,
        -1, -1, 1, 1,
        -1, -1, -1, -1,
        1, 1, 1, 1};

    const int nx = dims[0];
    const int ny = dims[1];
    const int nz = dims[2];
    const int xy = nx * ny;
    const int nvox = nx * ny * nz;

    float *seed_labels;
    float *dist_map;
    int i, n;
    int active_count = 0;
    int changed = 1;
    int iter = 0;
    const int maxiter = 2000;
    int offsets[26];
    float ndist[26];
    int own_dist = 0;
    const int use_legacy_default_cost = (dd == NULL);

    float sx = 1.0f, sy = 1.0f, sz = 1.0f;
    float w_dist = 0.1f, w_int = 10.0f;
    const float limit_f = (float)limit;

    if (!labels || !intensity || !dims)
    {
        fprintf(stderr, "Invalid NULL input in downcut_float\n");
        exit(EXIT_FAILURE);
    }

    if (voxelsize)
    {
        sx = (float)fabs(voxelsize[0]);
        sy = (float)fabs(voxelsize[1]);
        sz = (float)fabs(voxelsize[2]);
    }

    if (dd)
    {
        w_dist = (float)dd[0];
        w_int = (float)dd[1];
    }

    seed_labels = (float *)malloc((size_t)nvox * sizeof(float));
    if (!seed_labels)
    {
        fprintf(stderr, "Memory allocation error in downcut_float\n");
        exit(EXIT_FAILURE);
    }

    if (dist)
        dist_map = dist;
    else
    {
        dist_map = (float *)malloc((size_t)nvox * sizeof(float));
        if (!dist_map)
        {
            free(seed_labels);
            fprintf(stderr, "Memory allocation error in downcut_float\n");
            exit(EXIT_FAILURE);
        }
        own_dist = 1;
    }

    const float s12 = (float)sqrt((double)sx * sx + (double)sy * sy);
    const float s13 = (float)sqrt((double)sx * sx + (double)sz * sz);
    const float s23 = (float)sqrt((double)sy * sy + (double)sz * sz);
    const float s123 = (float)sqrt((double)s12 * s12 + (double)sz * sz);

    const float nd_template[26] = {
        sx, sx, sy, sy, sz, sz,
        s12, s12, s12, s12,
        s13, s13, s13, s13,
        s23, s23, s23, s23,
        s123, s123, s123, s123,
        s123, s123, s123, s123};

    for (n = 0; n < 26; ++n)
    {
        offsets[n] = DX[n] + (DY[n] * nx) + (DZ[n] * xy);
        ndist[n] = nd_template[n];
    }

    for (i = 0; i < nvox; ++i)
    {
        float seed = labels[i];

        if (!isfinite(seed))
            seed = 0.0f;

        seed_labels[i] = seed;
        labels[i] = seed;

        if (seed == 0.0f)
            dist_map[i] = FLT_MAX;
        else if (seed == -FLT_MAX)
            dist_map[i] = -FLT_MAX;
        else
        {
            dist_map[i] = 0.0f;
            active_count++;
        }
    }

    while (active_count > 0 && iter < maxiter && changed > 0)
    {
        iter++;
        changed = 0;

        for (i = 0; i < nvox; ++i)
        {
            if (dist_map[i] <= 0.0f && dist_map[i] != -FLT_MAX)
            {
                int x, y, z;
                const int plane_offset = i % xy;

                if (dist_map[i] < 0.0f)
                    dist_map[i] = -dist_map[i];

                active_count--;

                z = i / xy;
                y = plane_offset / nx;
                x = plane_offset % nx;

                for (n = 0; n < 26; ++n)
                {
                    const int xn = x + DX[n];
                    const int yn = y + DY[n];
                    const int zn = z + DZ[n];
                    const int ni = i + offsets[n];

                    if (xn < 0 || xn >= nx || yn < 0 || yn >= ny || zn < 0 || zn >= nz)
                        continue;

                    if ((intensity[i] + limit_f) >= intensity[ni] && seed_labels[ni] == 0.0f)
                    {
                        float intensity_cost;
                        float distn;

                        if (use_legacy_default_cost)
                            intensity_cost = fmaxf(0.0f, 4.0f - fmaxf(0.0f, intensity[ni]));
                        else
                            intensity_cost = fmaxf(0.0f, fminf(1.0f, intensity[ni]));

                        distn = dist_map[i] + (w_dist * ndist[n]) + (w_int * intensity_cost);

                        if (dist_map[ni] != -FLT_MAX && fabsf(dist_map[ni]) > fabsf(distn))
                        {
                            if (dist_map[ni] > 0.0f)
                                active_count++;

                            changed++;
                            dist_map[ni] = -distn;
                            labels[ni] = labels[i];
                        }
                    }
                }

                if (dist_map[i] == 0.0f)
                    dist_map[i] = -FLT_MAX;
            }
        }
    }

    if (own_dist)
        free(dist_map);
    free(seed_labels);
}

/**
 * \brief Datatype-generic wrapper for downcut region growing.
 *
 * Converts `labels` and `intensity` from their input datatypes to floating-point,
 * runs downcut propagation, then converts results back to requested output datatypes
 * using `convert_input_type` and `convert_output_type`.
 *
 * \param labels            (in/out) generic label buffer
 * \param intensity         (in)     generic intensity buffer
 * \param dist              (out)    generic distance/path-cost buffer (NULL allowed)
 * \param dims              (in)     volume dimensions {nx, ny, nz}
 * \param limit             (in)     intensity limit for neighbour acceptance
 * \param voxelsize         (in)     voxel spacing {sx, sy, sz}; NULL -> {1,1,1}
 * \param dd                (in)     distance/intensity weights {w_dist, w_int}; NULL -> defaults
 * \param labels_datatype   (in)     datatype code of labels input/output
 * \param intensity_datatype(in)     datatype code of intensity input
 * \param dist_datatype     (in)     datatype code of dist output (ignored if dist==NULL)
 */
void downcut3(void *labels, void *intensity, void *dist,
              int dims[3], double limit, double voxelsize[3], double dd[2],
              int labels_datatype, int intensity_datatype, int dist_datatype)
{
    const int nvox = dims[0] * dims[1] * dims[2];
    double *labels_d;
    double *intensity_d;
    double *dist_d = NULL;
    float *labels_f;
    float *intensity_f;
    float *dist_f = NULL;
    int i;

    if (!labels || !intensity || !dims)
    {
        fprintf(stderr, "Invalid NULL input in downcut3\n");
        exit(EXIT_FAILURE);
    }

    labels_d = (double *)malloc((size_t)nvox * sizeof(double));
    intensity_d = (double *)malloc((size_t)nvox * sizeof(double));
    labels_f = (float *)malloc((size_t)nvox * sizeof(float));
    intensity_f = (float *)malloc((size_t)nvox * sizeof(float));

    if (dist)
    {
        dist_d = (double *)malloc((size_t)nvox * sizeof(double));
        dist_f = (float *)malloc((size_t)nvox * sizeof(float));
    }

    if (!labels_d || !intensity_d || !labels_f || !intensity_f || (dist && (!dist_d || !dist_f)))
    {
        free(labels_d);
        free(intensity_d);
        free(labels_f);
        free(intensity_f);
        free(dist_d);
        free(dist_f);
        fprintf(stderr, "Memory allocation error in downcut3\n");
        exit(EXIT_FAILURE);
    }

    convert_input_type(labels, labels_d, nvox, labels_datatype);
    convert_input_type(intensity, intensity_d, nvox, intensity_datatype);

    for (i = 0; i < nvox; ++i)
    {
        labels_f[i] = (float)labels_d[i];
        intensity_f[i] = (float)intensity_d[i];
    }

    downcut_float(labels_f, intensity_f, dist_f, dims, limit, voxelsize, dd);

    for (i = 0; i < nvox; ++i)
        labels_d[i] = (double)labels_f[i];
    convert_output_type(labels, labels_d, nvox, labels_datatype);

    if (dist)
    {
        for (i = 0; i < nvox; ++i)
            dist_d[i] = (double)dist_f[i];
        convert_output_type(dist, dist_d, nvox, dist_datatype);
    }

    free(labels_d);
    free(intensity_d);
    free(labels_f);
    free(intensity_f);
    free(dist_d);
    free(dist_f);
}

/**
 * \brief Blood-vessel correction for PVE label maps.
 *
 * Implements blood-vessel correction for a PVE map in the range [0..3].
 * A safe WM seed region is estimated by thresholding and distance-based opening,
 * then expanded with downcut region growing constrained by transformed intensities.
 * Identified vessel-like WM outliers are reset and locally median-smoothed.
 *
 * This implementation uses CAT-Surface library tools only:
 * - `downcut_float()` for constrained region growing
 * - `dist_open_float()` / `dist_dilate_float()` for distance morphology
 * - `median3()` for masked local smoothing
 *
 * \param Yp0     (in/out) float PVE label image in [0..3]
 * \param dims    (in)     dimensions {nx, ny, nz}
 * \param vx_vol  (in)     voxel spacing {sx, sy, sz}; NULL -> {1,1,1}
 */
void blood_vessel_correction_pve_float(float *Yp0, int dims[3], double vx_vol[3])
{
    const int nvox = dims[0] * dims[1] * dims[2];
    double vx_local[3] = {1.0, 1.0, 1.0};
    float *F;
    float *YwmA;
    float *YwmB;
    float *Ywm;
    float *Yd;
    float *Yp0s;
    unsigned char *Ymsk;
    int i;

    if (!Yp0 || !dims)
    {
        fprintf(stderr, "Invalid NULL input in blood_vessel_correction_pve_float\n");
        exit(EXIT_FAILURE);
    }

    if (vx_vol)
    {
        vx_local[0] = vx_vol[0];
        vx_local[1] = vx_vol[1];
        vx_local[2] = vx_vol[2];
    }

    F = (float *)malloc((size_t)nvox * sizeof(float));
    YwmA = (float *)malloc((size_t)nvox * sizeof(float));
    YwmB = (float *)malloc((size_t)nvox * sizeof(float));
    Ywm = (float *)malloc((size_t)nvox * sizeof(float));
    Yd = (float *)malloc((size_t)nvox * sizeof(float));
    Yp0s = (float *)malloc((size_t)nvox * sizeof(float));
    Ymsk = (unsigned char *)malloc((size_t)nvox * sizeof(unsigned char));

    if (!F || !YwmA || !YwmB || !Ywm || !Yd || !Yp0s || !Ymsk)
    {
        free(F);
        free(YwmA);
        free(YwmB);
        free(Ywm);
        free(Yd);
        free(Yp0s);
        free(Ymsk);
        fprintf(stderr, "Memory allocation error in blood_vessel_correction_pve_float\n");
        exit(EXIT_FAILURE);
    }

    /*
     * F = max(0, Yp0 - 1); F(Yp0 <= 1.1) = inf;
     */
    for (i = 0; i < nvox; ++i)
    {
        float fval = Yp0[i] - 1.0f;
        if (fval < 0.0f)
            fval = 0.0f;
        if (Yp0[i] <= 1.1f)
            fval = FLT_MAX;
        F[i] = fval;
    }

    /*
     * Ywm = morph(Yp0>2.5,'ldo',2,vx) | morph(Yp0>2.75,'ldo',1,vx)
     * Approximated with distance-based opening for logical masks.
     */
    for (i = 0; i < nvox; ++i)
    {
        YwmA[i] = (Yp0[i] > 2.5f) ? 1.0f : 0.0f;
        YwmB[i] = (Yp0[i] > 2.75f) ? 1.0f : 0.0f;
    }
    dist_open_float(YwmA, dims, vx_local, 2.0, 0.5);
    dist_open_float(YwmB, dims, vx_local, 1.0, 0.5);
    keep_largest_cluster(YwmA, 0.5, dims, DT_FLOAT32, 0, 1, 18);
    keep_largest_cluster(YwmB, 0.5, dims, DT_FLOAT32, 0, 1, 18);
    for (i = 0; i < nvox; ++i)
        Ywm[i] = (YwmA[i] > 0.0f || YwmB[i] > 0.0f) ? 1.0f : 0.0f;

    /*
     * [~,Yd] = cat_vol_downcut(single(Ywm), F, -0.001)
     */
    downcut_float(Ywm, F, Yd, dims, -0.001, vx_local, NULL);

    /*
     * Ymsk = Yd > 1000000 & Yp0 > 2
     */
    for (i = 0; i < nvox; ++i)
        Ymsk[i] = (Yd[i] > 1000000.0f && Yp0[i] > 2.0f) ? 1 : 0;

    /*
     * Yp0(Ymsk) = 2
     */
    for (i = 0; i < nvox; ++i)
        if (Ymsk[i])
            Yp0[i] = 2.0f;

    /*
     * Ymsk = cat_vol_morph(Ymsk,'dd',1,vx) & Yp0>1.5
     * 'dd' maps to distance-based dilation.
     */
    for (i = 0; i < nvox; ++i)
        Ywm[i] = (float)Ymsk[i];
    dist_dilate_float(Ywm, dims, vx_local, 1.0, 0.5);
    for (i = 0; i < nvox; ++i)
        Ymsk[i] = (Ywm[i] > 0.0f && Yp0[i] > 1.5f) ? 1 : 0;

    /*
     * Yp0s = cat_vol_median3(Yp0, Ymsk)
     * Yp0(Ymsk) = Yp0s(Ymsk)
     */
    memcpy(Yp0s, Yp0, (size_t)nvox * sizeof(float));
    median3(Yp0s, Ymsk, dims, 1, DT_FLOAT32);
    for (i = 0; i < nvox; ++i)
        if (Ymsk[i])
            Yp0[i] = Yp0s[i];

    free(F);
    free(YwmA);
    free(YwmB);
    free(Ywm);
    free(Yd);
    free(Yp0s);
    free(Ymsk);
}

/**
 * \brief Datatype-generic blood-vessel correction wrapper for PVE labels.
 *
 * Converts input PVE labels to float, runs `blood_vessel_correction_pve_float()`, and converts back
 * to the requested datatype.
 *
 * \param data      (in/out) PVE label volume data
 * \param dims      (in)     dimensions {nx, ny, nz}
 * \param vx_vol    (in)     voxel spacing {sx, sy, sz}; NULL -> {1,1,1}
 * \param datatype  (in)     datatype code (DT_UINT8, DT_UINT16, DT_FLOAT32, etc.)
 */
void blood_vessel_correction_pve(void *data, int dims[3], double vx_vol[3], int datatype)
{
    const int nvox = dims[0] * dims[1] * dims[2];
    float *buffer;

    if (!data || !dims)
    {
        fprintf(stderr, "Invalid NULL input in blood_vessel_correction_pve\n");
        exit(EXIT_FAILURE);
    }

    buffer = (float *)malloc((size_t)nvox * sizeof(float));
    if (!buffer)
    {
        fprintf(stderr, "Memory allocation error in blood_vessel_correction_pve\n");
        exit(EXIT_FAILURE);
    }

    convert_input_type_float(data, buffer, nvox, datatype);
    blood_vessel_correction_pve_float(buffer, dims, vx_vol);
    convert_output_type_float(data, buffer, nvox, datatype);

    free(buffer);
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
 *  - The function uses the neighbouring values to calculate the Laplace filtering.
 */
void laplace3R(float *SEG, unsigned char *M, int dims[3], double TH)
{
    const int x = dims[0];
    const int y = dims[1];

    const int z = dims[2];
    const int xy = x * y;
    const int nvox = x * y * z;

    // Indices of the neighbour and size of neighbours array
    const int NI[6] = {-1, 1, -x, x, -xy, xy};
    const int sN = sizeof(NI) / sizeof(NI[0]);

    int i, n, u, v, w, nu, nv, nw, ni, iter = 0, maxiter = 2000;
    float Nn, diff, maxdiffi, maxdiff = 1.0;

    // Allocate memory for Laplace calculation
    float *L1 = (float *)malloc(sizeof(float) * nvox);
    float *L2 = (float *)malloc(sizeof(float) * nvox);
    unsigned char *LN = (unsigned char *)malloc(sizeof(unsigned char) * nvox);

    // Check for successful memory allocation
    if (!L1 || !L2 || !LN)
    {
        fprintf(stderr, "Memory allocation error\n");
        exit(EXIT_FAILURE);
    }

    // Initialization
    for (i = 0; i < nvox; i++)
    {
        L1[i] = isnan(SEG[i]) ? FLT_MAX : SEG[i];
        L2[i] = L1[i];
        LN[i] = M[i];
    }

    // Iterative Laplace filtering
    while (maxdiff > TH && iter < maxiter)
    {
        maxdiffi = 0;
        iter++;
        for (i = 0; i < nvox; i++)
        {
            if (M[i] && LN[i])
            {
                ind2sub(i, &u, &v, &w, xy, x);

                // Read neighbour values
                L2[i] = 0.0;
                Nn = 0.0;
                for (n = 0; n < sN; n++)
                {
                    ni = i + NI[n];
                    ind2sub(ni, &nu, &nv, &nw, xy, x);
                    if (((ni < 0) || (ni >= nvox) || (abs(nu - u) > 1) || (abs(nv - v) > 1) || (abs(nw - w) > 1) || (L1[ni] == -FLT_MAX) || (L1[ni] == FLT_MAX)) == 0)
                    {
                        L2[i] += L1[ni];
                        Nn++;
                    }
                }

                L2[i] = Nn > 0 ? L2[i] / (float)Nn : L1[i];
                diff = fabs(L1[i] - L2[i]);
                if (diff > (TH / 10.0))
                {
                    for (n = 0; n < sN; n++)
                    {
                        ni = i + NI[n];
                        ind2sub(ni, &nu, &nv, &nw, xy, x);
                        if (((ni < 0) || (ni >= nvox) || (abs(nu - u) > 1) || (abs(nv - v) > 1) || (abs(nw - w) > 1) || (L1[ni] == -FLT_MAX) || (L1[ni] == FLT_MAX)) == 0)
                            LN[ni] = 1; /* if I change his neighbours it has to be recalculated */
                    }
                }
                LN[i] = 0;
                if (maxdiffi < diff)
                    maxdiffi = diff;
            }
        }
        maxdiff = maxdiffi;

        // Update L1 with the new values from L2
        for (i = 0; i < nvox; i++)
            L1[i] = L2[i];
    }

    // Copy the final result back into SEG
    for (i = 0; i < nvox; i++)
        SEG[i] = L1[i];

    // Free allocated memory
    free(L1);
    free(L2);
    free(LN);
}

/**
 * \brief Morphological erosion (binary) using the Euclidean distance transform.
 *
 * Thresholds the input to a binary mask M = (vol > th * max(vol)).
 * Erosion by radius 'dist' is implemented as:
 *
 *   E = { x | dist_to_background(x) > dist }.
 *
 * We compute distance to the background by running the distance transform
 * on the inverted mask (background==1), then thresholding:
 *
 *   E = ( EDT( 1 - M ) > dist ).
 *
 * \param vol        (in/out) float[dims[0]*dims[1]*dims[2]]; overwritten with 0/1
 * \param dims       (in)     {nx, ny, nz}
 * \param voxelsize  (in)     voxel spacing in mm (or consistent units)
 * \param dist       (in)     structuring radius in same units as voxelsize (<=0: no-op)
 * \param th         (in)     threshold as fraction of max(vol) in [0,1]
 */
void dist_erode_float(float *vol, int dims[3], double voxelsize[3], double dist, double th)
{
    if (dist <= 0.0)
        return;

    const int nvox = dims[0] * dims[1] * dims[2];
    float *buffer = (float *)malloc(sizeof(float) * nvox);
    int i;

    if (!buffer)
    {
        fprintf(stderr, "Memory allocation error in dist_erode_float\n");
        exit(EXIT_FAILURE);
    }

    /* Threshold to binary mask relative to max */
    float max_vol = get_max(vol, nvox, 0, DT_FLOAT32);
    const float thr = (float)(th * (double)max_vol);

    /* Invert mask to measure distance to background (background = 1, object = 0) */
    for (i = 0; i < nvox; ++i)
        buffer[i] = (vol[i] > thr) ? 0.0f : 1.0f;

    /* Distance to background, then keep voxels farther than dist from background */
    euclidean_distance(buffer, NULL, dims, voxelsize, 0);
    for (i = 0; i < nvox; ++i)
        buffer[i] = (buffer[i] > (float)dist) ? 1.0f : 0.0f;

    /* Write back binary result */
    for (i = 0; i < nvox; ++i)
        vol[i] = buffer[i];

    free(buffer);
}

/**
 * \brief Wrapper for morphological erosion (generic datatype).
 *
 * Applies dist_erode_float() by converting input to float, processing, and converting back.
 * This allows erosion to work with any supported image datatype (uint8, uint16, float, etc.).
 *
 * \param data       (in/out) void pointer to image data; type given by datatype parameter
 * \param dims       (in)     {nx, ny, nz}
 * \param voxelsize  (in)     voxel spacing in mm (or consistent units)
 * \param dist       (in)     structuring radius in same units as voxelsize (<=0: no-op)
 * \param th         (in)     threshold as fraction of max(data) in [0,1]
 * \param datatype   (in)     data type code (DT_UINT8, DT_UINT16, DT_FLOAT32, etc.)
 */
void dist_erode(void *data, int dims[3], double voxelsize[3], double dist, double th, int datatype)
{
    int nvox;
    float *buffer;

    nvox = dims[0] * dims[1] * dims[2];
    buffer = (float *)malloc(sizeof(float) * nvox);

    /* check success of memory allocation */
    if (!buffer)
    {
        printf("Memory allocation error\n");
        exit(EXIT_FAILURE);
    }

    convert_input_type_float(data, buffer, nvox, datatype);
    dist_erode_float(buffer, dims, voxelsize, dist, th);
    convert_output_type_float(data, buffer, nvox, datatype);

    free(buffer);
}

/**
 * \brief Morphological dilation (binary) using the Euclidean distance transform.
 *
 * Thresholds the input to a binary mask M = (vol > th * max(vol)).
 * Dilation by radius 'dist' is implemented as:
 *
 *   D = { x | dist_to_foreground(x) <= dist }.
 *
 * We compute distance to the *foreground* by running the distance transform
 * on the mask itself (foreground==1), then thresholding:
 *
 *   D = ( EDT( M ) <= dist ).
 *
 * To avoid clipping at the volume borders (growth beyond edges), we add a
 * zero-valued band of width floor(dist) around the image, run the transforms,
 * then copy the center region back (as in dist_close_float).
 *
 * \param vol        (in/out) float[dims[0]*dims[1]*dims[2]]; overwritten with 0/1
 * \param dims       (in)     {nx, ny, nz}
 * \param voxelsize  (in)     voxel spacing in mm (or consistent units)
 * \param dist       (in)     structuring radius in same units as voxelsize (<=0: no-op)
 * \param th         (in)     threshold as fraction of max(vol) in [0,1]
 */
void dist_dilate_float(float *vol, int dims[3], double voxelsize[3], double dist, double th)
{
    if (dist <= 0.0)
        return;

    /* pad band to avoid border clipping during dilation */
    const int band = (int)floor(dist);
    int d, i, x, y, z;

    int dims2[3];
    for (d = 0; d < 3; ++d)
        dims2[d] = dims[d] + 2 * band;
    const int nvox2 = dims2[0] * dims2[1] * dims2[2];

    float *buffer = (float *)malloc(sizeof(float) * nvox2);
    if (!buffer)
    {
        fprintf(stderr, "Memory allocation error in dist_dilate_float\n");
        exit(EXIT_FAILURE);
    }
    memset(buffer, 0, sizeof(float) * nvox2);

    const int nx = dims[0], ny = dims[1], nz = dims[2];

    /* Threshold to binary mask relative to max */
    const int nvox = nx * ny * nz;
    float max_vol = get_max(vol, nvox, 0, DT_FLOAT32);
    const float thr = (float)(th * (double)max_vol);

    /* Embed thresholded mask into padded buffer at offset 'band' */
    for (z = 0; z < nz; ++z)
        for (y = 0; y < ny; ++y)
            for (x = 0; x < nx; ++x)
            {
                int dst = sub2ind(x + band, y + band, z + band, dims2);
                int src = sub2ind(x, y, z, dims);
                buffer[dst] = (vol[src] > thr) ? 1.0f : 0.0f;
            }

    /* Distance to foreground, then keep voxels within dist of the original mask */
    euclidean_distance(buffer, NULL, dims2, voxelsize, 0);
    for (i = 0; i < nvox2; ++i)
        buffer[i] = (buffer[i] <= (float)dist) ? 1.0f : 0.0f;

    /* Copy center region back to 'vol' */
    for (z = 0; z < nz; ++z)
        for (y = 0; y < ny; ++y)
            for (x = 0; x < nx; ++x)
            {
                int src = sub2ind(x + band, y + band, z + band, dims2);
                int dst = sub2ind(x, y, z, dims);
                vol[dst] = buffer[src];
            }

    free(buffer);
}

/**
 * \brief Wrapper for morphological dilation (generic datatype).
 *
 * Applies dist_dilate_float() by converting input to float, processing, and converting back.
 * This allows dilation to work with any supported image datatype (uint8, uint16, float, etc.).
 * Growth is computed using the Euclidean distance transform with padding to avoid border clipping.
 *
 * \param data       (in/out) void pointer to image data; type given by datatype parameter
 * \param dims       (in)     {nx, ny, nz}
 * \param voxelsize  (in)     voxel spacing in mm (or consistent units)
 * \param dist       (in)     structuring radius in same units as voxelsize (<=0: no-op)
 * \param th         (in)     threshold as fraction of max(data) in [0,1]
 * \param datatype   (in)     data type code (DT_UINT8, DT_UINT16, DT_FLOAT32, etc.)
 */
void dist_dilate(void *data, int dims[3], double voxelsize[3], double dist, double th, int datatype)
{
    int nvox;
    float *buffer;

    nvox = dims[0] * dims[1] * dims[2];
    buffer = (float *)malloc(sizeof(float) * nvox);

    /* check success of memory allocation */
    if (!buffer)
    {
        printf("Memory allocation error\n");
        exit(EXIT_FAILURE);
    }

    convert_input_type_float(data, buffer, nvox, datatype);
    dist_dilate_float(buffer, dims, voxelsize, dist, th);
    convert_output_type_float(data, buffer, nvox, datatype);

    free(buffer);
}

/**
 * \brief Morphological closing (binary) using the Euclidean distance transform.
 *
 * Closing combines the operations:
 *   1. Dilation: expand foreground by radius 'dist'
 *   2. Erosion: shrink foreground by radius 'dist'
 *
 * This operation fills small holes (background regions surrounded by foreground)
 * while preserving overall shape. Computed via distance transform with padding
 * to prevent border clipping.
 *
 * \param vol        (in/out) float[dims[0]*dims[1]*dims[2]]; overwritten with 0/1
 * \param dims       (in)     {nx, ny, nz}
 * \param voxelsize  (in)     voxel spacing in mm (or consistent units)
 * \param dist       (in)     structuring radius in same units as voxelsize (<=0: no-op)
 * \param th         (in)     threshold as fraction of max(vol) in [0,1]
 */
void dist_close_float(float *vol, int dims[3], double voxelsize[3], double dist, double th)
{
    float *buffer;
    int i, x, y, z, j, band, dims2[3];
    float max_vol;
    int nvox2, nvox = dims[0] * dims[1] * dims[2];

    if (dist <= 0.0)
        return;

    max_vol = get_max(vol, nvox, 0, DT_FLOAT32);
    th *= (double)max_vol;

    /* add band with zeros to image to avoid clipping */
    band = floor(dist);
    for (i = 0; i < 3; i++)
        dims2[i] = dims[i] + 2 * band;
    nvox2 = dims2[0] * dims2[1] * dims2[2];

    buffer = (float *)malloc(sizeof(float) * dims2[0] * dims2[1] * dims2[2]);

    if (!buffer)
    {
        fprintf(stderr, "Memory allocation error\n");
        exit(EXIT_FAILURE);
    }

    memset(buffer, 0, sizeof(float) * nvox2);

    /* threshold input */
    for (z = 0; z < dims[2]; z++)
        for (y = 0; y < dims[1]; y++)
            for (x = 0; x < dims[0]; x++)
                buffer[sub2ind(x + band, y + band, z + band, dims2)] = (vol[sub2ind(x, y, z, dims)] > (float)th);

    euclidean_distance(buffer, NULL, dims2, voxelsize, 0);
    for (i = 0; i < nvox2; i++)
        buffer[i] = buffer[i] > (float)dist;

    euclidean_distance(buffer, NULL, dims2, voxelsize, 0);
    for (i = 0; i < nvox2; i++)
        buffer[i] = buffer[i] > (float)dist;

    /* return image */
    for (z = 0; z < dims[2]; z++)
        for (y = 0; y < dims[1]; y++)
            for (x = 0; x < dims[0]; x++)
                vol[sub2ind(x, y, z, dims)] = buffer[sub2ind(x + band, y + band, z + band, dims2)];

    free(buffer);
}

/**
 * \brief Wrapper for morphological closing (generic datatype).
 *
 * Applies dist_close_float() by converting input to float, processing, and converting back.
 * This allows closing to work with any supported image datatype (uint8, uint16, float, etc.).
 *
 * \param data       (in/out) void pointer to image data; type given by datatype parameter
 * \param dims       (in)     {nx, ny, nz}
 * \param voxelsize  (in)     voxel spacing in mm (or consistent units)
 * \param dist       (in)     structuring radius in same units as voxelsize (<=0: no-op)
 * \param th         (in)     threshold as fraction of max(data) in [0,1]
 * \param datatype   (in)     data type code (DT_UINT8, DT_UINT16, DT_FLOAT32, etc.)
 */
void dist_close(void *data, int dims[3], double voxelsize[3], double dist, double th, int datatype)
{
    int nvox;
    float *buffer;

    nvox = dims[0] * dims[1] * dims[2];
    buffer = (float *)malloc(sizeof(float) * nvox);

    /* check success of memory allocation */
    if (!buffer)
    {
        printf("Memory allocation error\n");
        exit(EXIT_FAILURE);
    }

    convert_input_type_float(data, buffer, nvox, datatype);
    dist_close_float(buffer, dims, voxelsize, dist, th);
    convert_output_type_float(data, buffer, nvox, datatype);

    free(buffer);
}

/**
 * \brief Morphological opening (binary) using the Euclidean distance transform.
 *
 * Opening combines the operations:
 *   1. Erosion: shrink foreground by radius 'dist'
 *   2. Dilation: expand foreground by radius 'dist'
 *
 * This operation removes small foreground objects (noise) while preserving the
 * shape of larger structures. Computed via distance transform on the thresholded
 * input without explicit padding.
 *
 * \param vol        (in/out) float[dims[0]*dims[1]*dims[2]]; overwritten with 0/1
 * \param dims       (in)     {nx, ny, nz}
 * \param voxelsize  (in)     voxel spacing in mm (or consistent units)
 * \param dist       (in)     structuring radius in same units as voxelsize (<=0: no-op)
 * \param th         (in)     threshold as fraction of max(vol) in [0,1]
 */
void dist_open_float(float *vol, int dims[3], double voxelsize[3], double dist, double th)
{
    float *buffer;
    int i, j;
    float max_vol;
    int nvox = dims[0] * dims[1] * dims[2];

    if (dist <= 0.0)
        return;

    max_vol = get_max(vol, nvox, 0, DT_FLOAT32);
    th *= (double)max_vol;

    buffer = (float *)malloc(sizeof(float) * nvox);

    if (!buffer)
    {
        fprintf(stderr, "Memory allocation error\n");
        exit(EXIT_FAILURE);
    }

    /* threshold input */
    for (i = 0; i < nvox; i++)
        buffer[i] = 1.0 - ((float)vol[i] > th);

    euclidean_distance(buffer, NULL, dims, voxelsize, 0);
    for (i = 0; i < nvox; i++)
        buffer[i] = buffer[i] > (float)dist;

    euclidean_distance(buffer, NULL, dims, voxelsize, 0);
    for (i = 0; i < nvox; i++)
        buffer[i] = buffer[i] <= (float)dist;

    /* return image */
    for (i = 0; i < nvox; i++)
        vol[i] = buffer[i];

    free(buffer);
}

/**
 * \brief Wrapper for morphological opening (generic datatype).
 *
 * Applies dist_open_float() by converting input to float, processing, and converting back.
 * This allows opening to work with any supported image datatype (uint8, uint16, float, etc.).
 *
 * \param data       (in/out) void pointer to image data; type given by datatype parameter
 * \param dims       (in)     {nx, ny, nz}
 * \param voxelsize  (in)     voxel spacing in mm (or consistent units)
 * \param dist       (in)     structuring radius in same units as voxelsize (<=0: no-op)
 * \param th         (in)     threshold as fraction of max(data) in [0,1]
 * \param datatype   (in)     data type code (DT_UINT8, DT_UINT16, DT_FLOAT32, etc.)
 */
void dist_open(void *data, int dims[3], double voxelsize[3], double dist, double th, int datatype)
{
    int nvox;
    float *buffer;

    nvox = dims[0] * dims[1] * dims[2];
    buffer = (float *)malloc(sizeof(float) * nvox);

    /* check success of memory allocation */
    if (!buffer)
    {
        printf("Memory allocation error\n");
        exit(EXIT_FAILURE);
    }

    convert_input_type_float(data, buffer, nvox, datatype);
    dist_open_float(buffer, dims, voxelsize, dist, th);
    convert_output_type_float(data, buffer, nvox, datatype);

    free(buffer);
}

/**
 * \brief Binary morphological erosion using separable convolution.
 *
 * Thresholds input to a binary mask M = (vol > th * max(vol)), then applies
 * erosion with a 3x3x3 cube structuring element repeated niter times.
 * Erosion removes foreground voxels on the boundary, effectively shrinking objects.
 *
 * \param vol        (in/out) float[dims[0]*dims[1]*dims[2]]; overwritten with 0/1
 * \param dims       (in)     {nx, ny, nz}
 * \param niter      (in)     number of erosion iterations (<=0: no-op)
 * \param th         (in)     threshold as fraction of max(vol) in [0,1]
 */
void morph_erode_float(float *vol, int dims[3], int niter, double th)
{
    double filt[3] = {1, 1, 1};
    int i, j;
    float max_vol;
    int nvox = dims[0] * dims[1] * dims[2];

    if (niter < 1)
        return;

    max_vol = get_max(vol, nvox, 0, DT_FLOAT32);
    th *= (double)max_vol;

    /* threshold input */
    for (j = 0; j < nvox; j++)
        vol[j] = vol[j] > (float)th;

    for (i = 0; i < niter; i++)
    {
        convxyz_float(vol, filt, filt, filt, 3, 3, 3, -1, -1, -1, vol, dims);
        for (j = 0; j < nvox; j++)
            vol[j] = vol[j] >= 9.0;
    }
}

/**
 * \brief Wrapper for binary morphological erosion (generic datatype).
 *
 * Applies morph_erode_float() by converting input to float, processing, and converting back.
 * Erosion is performed using separable convolution with a 3x3x3 structuring element.
 *
 * \param data       (in/out) void pointer to image data; type given by datatype parameter
 * \param dims       (in)     {nx, ny, nz}
 * \param niter      (in)     number of erosion iterations (<=0: no-op)
 * \param th         (in)     threshold as fraction of max(data) in [0,1]
 * \param datatype   (in)     data type code (DT_UINT8, DT_UINT16, DT_FLOAT32, etc.)
 */
void morph_erode(void *data, int dims[3], int niter, double th, int datatype)
{
    int nvox;
    float *buffer;

    nvox = dims[0] * dims[1] * dims[2];
    buffer = (float *)malloc(sizeof(float) * nvox);

    /* check success of memory allocation */
    if (!buffer)
    {
        printf("Memory allocation error\n");
        exit(EXIT_FAILURE);
    }

    convert_input_type_float(data, buffer, nvox, datatype);
    morph_erode_float(buffer, dims, niter, th);
    convert_output_type_float(data, buffer, nvox, datatype);

    free(buffer);
}

/**
 * \brief Binary morphological dilation using separable convolution.
 *
 * Thresholds input to a binary mask M = (vol > th * max(vol)), then applies
 * dilation with a 3x3x3 cube structuring element repeated niter times.
 * Input image is padded with a zero band to avoid clipping at volume borders.
 * Dilation expands foreground objects by niter voxels in all directions.
 *
 * \param vol        (in/out) float[dims[0]*dims[1]*dims[2]]; overwritten with 0/1
 * \param dims       (in)     {nx, ny, nz}
 * \param niter      (in)     number of dilation iterations (<=0: no-op)
 * \param th         (in)     threshold as fraction of max(vol) in [0,1]
 */
void morph_dilate_float(float *vol, int dims[3], int niter, double th)
{
    double filt[3] = {1, 1, 1};
    int i, x, y, z, j, band, dims2[3];
    float max_vol;
    unsigned char *buffer;
    int nvox = dims[0] * dims[1] * dims[2];

    if (niter < 1)
        return;

    max_vol = get_max(vol, nvox, 0, DT_FLOAT32);
    th *= (double)max_vol;

    /* add band with zeros to image to avoid clipping */
    band = niter;
    for (i = 0; i < 3; i++)
        dims2[i] = dims[i] + 2 * band;

    buffer = (unsigned char *)malloc(sizeof(unsigned char) * dims2[0] * dims2[1] * dims2[2]);

    if (!buffer)
    {
        fprintf(stderr, "Memory allocation error\n");
        exit(EXIT_FAILURE);
    }

    memset(buffer, 0, sizeof(unsigned char) * dims2[0] * dims2[1] * dims2[2]);

    /* threshold input */
    for (x = 0; x < dims[0]; x++)
        for (y = 0; y < dims[1]; y++)
            for (z = 0; z < dims[2]; z++)
                buffer[sub2ind(x + band, y + band, z + band, dims2)] = (unsigned char)((double)vol[sub2ind(x, y, z, dims)] > th);

    for (i = 0; i < niter; i++)
    {
        convxyz_uint8(buffer, filt, filt, filt, 3, 3, 3, -1, -1, -1, buffer, dims2);
        for (j = 0; j < dims2[0] * dims2[1] * dims2[2]; j++)
            buffer[j] = buffer[j] > 0;
    }

    /* return image */
    for (x = 0; x < dims[0]; x++)
        for (y = 0; y < dims[1]; y++)
            for (z = 0; z < dims[2]; z++)
                vol[sub2ind(x, y, z, dims)] = (float)buffer[sub2ind(x + band, y + band, z + band, dims2)];

    free(buffer);
}

/**
 * \brief Wrapper for binary morphological dilation (generic datatype).
 *
 * Applies morph_dilate_float() by converting input to float, processing, and converting back.
 * Dilation is performed using separable convolution with a 3x3x3 structuring element applied iteratively.
 *
 * \param data       (in/out) void pointer to image data; type given by datatype parameter
 * \param dims       (in)     {nx, ny, nz}
 * \param niter      (in)     number of dilation iterations (<=0: no-op)
 * \param th         (in)     threshold as fraction of max(data) in [0,1]
 * \param datatype   (in)     data type code (DT_UINT8, DT_UINT16, DT_FLOAT32, etc.)
 */
void morph_dilate(void *data, int dims[3], int niter, double th, int datatype)
{
    int nvox;
    float *buffer;

    nvox = dims[0] * dims[1] * dims[2];
    buffer = (float *)malloc(sizeof(float) * nvox);

    /* check success of memory allocation */
    if (!buffer)
    {
        printf("Memory allocation error\n");
        exit(EXIT_FAILURE);
    }

    convert_input_type_float(data, buffer, nvox, datatype);
    morph_dilate_float(buffer, dims, niter, th);
    convert_output_type_float(data, buffer, nvox, datatype);

    free(buffer);
}

/**
 * \brief Binary morphological closing (dilation followed by erosion).
 *
 * Thresholds input to a binary mask, then applies:
 *   1. niter dilations (using 3x3x3 structuring element)
 *   2. niter erosions (using 3x3x3 structuring element)
 *
 * Closing fills small holes (background regions) while preserving shape.
 * Input is padded with zeros to avoid border clipping during growth.
 *
 * \param vol        (in/out) float[dims[0]*dims[1]*dims[2]]; overwritten with 0/1
 * \param dims       (in)     {nx, ny, nz}
 * \param niter      (in)     number of iterations for each operation (<=0: no-op)
 * \param th         (in)     threshold as fraction of max(vol) in [0,1]
 */
void morph_close_float(float *vol, int dims[3], int niter, double th)
{
    double filt[3] = {1, 1, 1};
    unsigned char *buffer;
    int i, x, y, z, j, band, dims2[3];
    float max_vol;
    int nvox = dims[0] * dims[1] * dims[2];

    if (niter < 1)
        return;

    max_vol = get_max(vol, nvox, 0, DT_FLOAT32);
    th *= (double)max_vol;

    /* add band with zeros to image to avoid clipping */
    band = niter;
    for (i = 0; i < 3; i++)
        dims2[i] = dims[i] + 2 * band;

    buffer = (unsigned char *)malloc(sizeof(unsigned char) * dims2[0] * dims2[1] * dims2[2]);

    if (!buffer)
    {
        fprintf(stderr, "Memory allocation error\n");
        exit(EXIT_FAILURE);
    }

    memset(buffer, 0, sizeof(unsigned char) * dims2[0] * dims2[1] * dims2[2]);

    /* threshold input */
    for (x = 0; x < dims[0]; x++)
        for (y = 0; y < dims[1]; y++)
            for (z = 0; z < dims[2]; z++)
                buffer[sub2ind(x + band, y + band, z + band, dims2)] = (unsigned char)((double)vol[sub2ind(x, y, z, dims)] > th);

    /* dilate */
    for (i = 0; i < niter; i++)
    {
        convxyz_uint8(buffer, filt, filt, filt, 3, 3, 3, -1, -1, -1, buffer, dims2);
        for (j = 0; j < dims2[0] * dims2[1] * dims2[2]; j++)
            buffer[j] = (buffer[j] > 0);
    }

    /* erode */
    for (i = 0; i < niter; i++)
    {
        convxyz_uint8(buffer, filt, filt, filt, 3, 3, 3, -1, -1, -1, buffer, dims2);
        for (j = 0; j < dims2[0] * dims2[1] * dims2[2]; j++)
            buffer[j] = (buffer[j] >= 9);
    }

    /* return image */
    for (x = 0; x < dims[0]; x++)
        for (y = 0; y < dims[1]; y++)
            for (z = 0; z < dims[2]; z++)
                vol[sub2ind(x, y, z, dims)] = (float)buffer[sub2ind(x + band, y + band, z + band, dims2)];

    free(buffer);
}

/**
 * \brief Wrapper for binary morphological closing (generic datatype).
 *
 * Applies morph_close_float() by converting input to float, processing, and converting back.
 * Closing combines dilation followed by erosion using separable convolution.
 *
 * \param data       (in/out) void pointer to image data; type given by datatype parameter
 * \param dims       (in)     {nx, ny, nz}
 * \param niter      (in)     number of iterations for each operation (<=0: no-op)
 * \param th         (in)     threshold as fraction of max(data) in [0,1]
 * \param datatype   (in)     data type code (DT_UINT8, DT_UINT16, DT_FLOAT32, etc.)
 */
void morph_close(void *data, int dims[3], int niter, double th, int datatype)
{
    int nvox;
    float *buffer;

    nvox = dims[0] * dims[1] * dims[2];
    buffer = (float *)malloc(sizeof(float) * nvox);

    /* check success of memory allocation */
    if (!buffer)
    {
        printf("Memory allocation error\n");
        exit(EXIT_FAILURE);
    }

    convert_input_type_float(data, buffer, nvox, datatype);
    morph_close_float(buffer, dims, niter, th);
    convert_output_type_float(data, buffer, nvox, datatype);

    free(buffer);
}

/**
 * \brief Binary morphological opening (erosion followed by dilation).
 *
 * Thresholds input to a binary mask, then applies:
 *   1. niter erosions (using 3x3x3 structuring element)
 *   2. niter dilations (using 3x3x3 structuring element)
 *
 * Opening removes small foreground objects (noise) while preserving shape of larger structures.
 * If keep_values > 0, output marks only zero regions (areas successfully removed) from input;
 * otherwise outputs the binary result of opening.
 *
 * \param vol         (in/out) float[dims[0]*dims[1]*dims[2]]
 * \param dims        (in)     {nx, ny, nz}
 * \param niter       (in)     number of iterations for each operation (<=0: no-op)
 * \param th          (in)     threshold as fraction of max(vol) in [0,1]
 * \param keep_values (in)     if >0, preserve original foreground values and zero only removed regions
 */
void morph_open_float(float *vol, int dims[3], int niter, double th, int keep_values)
{
    unsigned char *buffer;
    double filt[3] = {1, 1, 1};
    int i, j, nvox;
    float max_vol;

    if (niter < 1)
        return;

    nvox = dims[0] * dims[1] * dims[2];
    max_vol = get_max(vol, nvox, 0, DT_FLOAT32);
    th *= (double)max_vol;

    buffer = (unsigned char *)malloc(sizeof(unsigned char) * nvox);

    if (!buffer)
    {
        fprintf(stderr, "Memory allocation error\n");
        exit(EXIT_FAILURE);
    }

    /* threshold input */
    for (j = 0; j < nvox; j++)
        buffer[j] = (unsigned char)(vol[j] > th);

    for (i = 0; i < niter; i++)
    {
        convxyz_uint8(buffer, filt, filt, filt, 3, 3, 3, -1, -1, -1, buffer, dims);
        for (j = 0; j < nvox; j++)
            buffer[j] = (buffer[j] >= 9);
    }

    for (i = 0; i < niter; i++)
    {
        convxyz_uint8(buffer, filt, filt, filt, 3, 3, 3, -1, -1, -1, buffer, dims);
        for (j = 0; j < nvox; j++)
            buffer[j] = (buffer[j] > 0);
    }

    // sometimes it's helpful to keep the original volume, but set areas to zero
    // where opening was successful
    if (keep_values > 0)
    {
        for (i = 0; i < nvox; i++)
            if (buffer[i] == 0)
                vol[i] = 0.0;
    }
    else
    {
        for (i = 0; i < nvox; i++)
            vol[i] = (float)buffer[i];
    }

    free(buffer);
}

/**
 * \brief Wrapper for binary morphological opening (generic datatype).
 *
 * Applies morph_open_float() by converting input to float, processing, and converting back.
 * Opening combines erosion followed by dilation using separable convolution.
 *
 * \param data        (in/out) void pointer to image data; type given by datatype parameter
 * \param dims        (in)     {nx, ny, nz}
 * \param niter       (in)     number of iterations for each operation (<=0: no-op)
 * \param th          (in)     threshold as fraction of max(data) in [0,1]
 * \param keep_values (in)     if >0, preserve original values and zero only removed regions
 * \param datatype    (in)     data type code (DT_UINT8, DT_UINT16, DT_FLOAT32, etc.)
 */
void morph_open(void *data, int dims[3], int niter, double th, int keep_values, int datatype)
{
    int nvox;
    float *buffer;

    nvox = dims[0] * dims[1] * dims[2];
    buffer = (float *)malloc(sizeof(float) * nvox);

    /* check success of memory allocation */
    if (!buffer)
    {
        printf("Memory allocation error\n");
        exit(EXIT_FAILURE);
    }

    convert_input_type_float(data, buffer, nvox, datatype);
    morph_open_float(buffer, dims, niter, th, keep_values);
    convert_output_type_float(data, buffer, nvox, datatype);

    free(buffer);
}

/**
 * \brief Grey-scale morphological erosion.
 *
 * Applies minimum filtering (gray erosion) to the entire volume using localstat_double().
 * Each voxel is set to the minimum value in its 3x3x3 neighborhood, repeated niter times.
 * Erosion darkens the image and removes bright details (small peaks).
 *
 * \param data       (in/out) void pointer to image data; type given by datatype parameter
 * \param dims       (in)     {nx, ny, nz}
 * \param niter      (in)     number of erosion iterations (<=0: no-op)
 * \param datatype   (in)     data type code (DT_UINT8, DT_UINT16, DT_FLOAT32, etc.)
 */
void grey_erode(void *data, int dims[3], int niter, int datatype)
{
    int nvox;
    double *buffer;

    nvox = dims[0] * dims[1] * dims[2];
    buffer = (double *)malloc(sizeof(double) * nvox);

    /* check success of memory allocation */
    if (!buffer)
    {
        printf("Memory allocation error\n");
        exit(EXIT_FAILURE);
    }

    convert_input_type(data, buffer, nvox, datatype);
    localstat_double(buffer, NULL, dims, 1, F_MIN, niter, 0);
    convert_output_type(data, buffer, nvox, datatype);

    free(buffer);
}

/**
 * \brief Grey-scale morphological dilation.
 *
 * Applies maximum filtering (gray dilation) to the entire volume using localstat_double().
 * Each voxel is set to the maximum value in its 3x3x3 neighborhood, repeated niter times.
 * Dilation brightens the image and removes dark details (small valleys).
 *
 * \param data       (in/out) void pointer to image data; type given by datatype parameter
 * \param dims       (in)     {nx, ny, nz}
 * \param niter      (in)     number of dilation iterations (<=0: no-op)
 * \param datatype   (in)     data type code (DT_UINT8, DT_UINT16, DT_FLOAT32, etc.)
 */
void grey_dilate(void *data, int dims[3], int niter, int datatype)
{
    int nvox;
    double *buffer;

    nvox = dims[0] * dims[1] * dims[2];
    buffer = (double *)malloc(sizeof(double) * nvox);

    /* check success of memory allocation */
    if (!buffer)
    {
        printf("Memory allocation error\n");
        exit(EXIT_FAILURE);
    }

    convert_input_type(data, buffer, nvox, datatype);
    localstat_double(buffer, NULL, dims, 1, F_MAX, niter, 0);
    convert_output_type(data, buffer, nvox, datatype);

    free(buffer);
}

/**
 * \brief Grey-scale morphological opening (erosion followed by dilation).
 *
 * Applies F_MIN (erosion) niter times, then F_MAX (dilation) niter times.
 * Opening removes small bright isolated voxels while preserving overall structure.
 * Implements gray opening = gray_dilate(gray_erode(vol)).
 *
 * \param data       (in/out) void pointer to image data; type given by datatype parameter
 * \param dims       (in)     {nx, ny, nz}
 * \param niter      (in)     number of iterations for each operation (<=0: no-op)
 * \param datatype   (in)     data type code (DT_UINT8, DT_UINT16, DT_FLOAT32, etc.)
 */
void grey_open(void *data, int dims[3], int niter, int datatype)
{
    int nvox;
    double *buffer;

    nvox = dims[0] * dims[1] * dims[2];
    buffer = (double *)malloc(sizeof(double) * nvox);

    /* check success of memory allocation */
    if (!buffer)
    {
        printf("Memory allocation error\n");
        exit(EXIT_FAILURE);
    }

    convert_input_type(data, buffer, nvox, datatype);
    localstat_double(buffer, NULL, dims, 1, F_MIN, niter, 0);
    localstat_double(buffer, NULL, dims, 1, F_MAX, niter, 0);
    convert_output_type(data, buffer, nvox, datatype);

    free(buffer);
}

/**
 * \brief Grey-scale morphological closing (dilation followed by erosion).
 *
 * Applies F_MAX (dilation) niter times, then F_MIN (erosion) niter times.
 * Closing fills small dark isolated voxels while preserving overall structure.
 * Implements gray closing = gray_erode(gray_dilate(vol)).
 *
 * \param data       (in/out) void pointer to image data; type given by datatype parameter
 * \param dims       (in)     {nx, ny, nz}
 * \param niter      (in)     number of iterations for each operation (<=0: no-op)
 * \param datatype   (in)     data type code (DT_UINT8, DT_UINT16, DT_FLOAT32, etc.)
 */
void grey_close(void *data, int dims[3], int niter, int datatype)
{
    int nvox;
    double *buffer;

    nvox = dims[0] * dims[1] * dims[2];
    buffer = (double *)malloc(sizeof(double) * nvox);

    /* check success of memory allocation */
    if (!buffer)
    {
        printf("Memory allocation error\n");
        exit(EXIT_FAILURE);
    }

    convert_input_type(data, buffer, nvox, datatype);
    localstat_double(buffer, NULL, dims, 1, F_MAX, niter, 0);
    localstat_double(buffer, NULL, dims, 1, F_MIN, niter, 0);
    convert_output_type(data, buffer, nvox, datatype);

    free(buffer);
}

void estimate_target_dimensions(int dims[3], double voxelsize[3], double target_voxelsize, int min_target_dim, int dims_samp[3])
{
    int i;

    for (i = 0; i < 3; i++)
    {
        // Compute target dimension
        dims_samp[i] = (int)ceil(dims[i] * (voxelsize[i] / target_voxelsize));

        // Ensure minimum target dimension
        if (dims_samp[i] < min_target_dim)
            dims_samp[i] = min_target_dim;
    }
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
    int x, y, z, x0, y0, z0, i;
    double k111, k112, k121, k122, k211, k212, k221, k222;
    double dx1, dx2, dy1, dy2, dz1, dz2, xi, yi, zi, samp[3];
    int off1, off2;

    // Estimate scaling factor
    for (i = 0; i < 3; i++)
    {
        samp[i] = (double)dims[i] / (double)dims_samp[i];
    }

    for (z = 0; z < dims_samp[2]; z++)
    {
        zi = z * samp[2];
        z0 = (int)floor(zi);
        dz1 = zi - z0;
        dz2 = 1.0 - dz1;
        z0 = z0 < dims[2] - 1 ? z0 : dims[2] - 2;

        for (y = 0; y < dims_samp[1]; y++)
        {
            yi = y * samp[1];
            y0 = (int)floor(yi);
            dy1 = yi - y0;
            dy2 = 1.0 - dy1;
            y0 = y0 < dims[1] - 1 ? y0 : dims[1] - 2;

            for (x = 0; x < dims_samp[0]; x++)
            {
                xi = x * samp[0];
                x0 = (int)floor(xi);
                dx1 = xi - x0;
                dx2 = 1.0 - dx1;
                x0 = x0 < dims[0] - 1 ? x0 : dims[0] - 2;

                i = z * dims_samp[0] * dims_samp[1] + y * dims_samp[0] + x;
                off1 = x0 + dims[0] * (y0 + dims[1] * z0);

                // Trilinear interpolation
                k222 = (double)in[off1];
                k122 = (double)in[off1 + 1];
                off2 = off1 + dims[0];
                k212 = (double)in[off2];
                k112 = (double)in[off2 + 1];
                off1 += dims[0] * dims[1];
                k221 = (double)in[off1];
                k121 = (double)in[off1 + 1];
                off2 = off1 + dims[0];
                k211 = (double)in[off2];
                k111 = (double)in[off2 + 1];

                out[i] = (float)((((k222 * dx2 + k122 * dx1) * dy2 + (k212 * dx2 + k112 * dx1) * dy1) * dz2) +
                                 (((k221 * dx2 + k121 * dx1) * dy2 + (k211 * dx2 + k111 * dx1) * dy1) * dz1));
            }
        }
    }
}

/**
 * \brief Resample a 3D volume to a different size using trilinear interpolation.
 *
 * Resizes a 3D volume from one dimension to another using trilinear interpolation.
 * The routine handles automatic type conversion for different data types.
 *
 * \param in        (in)  input volume data (pointer to any supported datatype)
 * \param out       (out) output volume data (pointer to pre-allocated array)
 * \param dims      (in)  original volume dimensions {nx, ny, nz}
 * \param dims_samp (in)  target volume dimensions {nx_new, ny_new, nz_new}
 * \param datatype  (in)  data type descriptor (e.g., DT_FLOAT32, DT_UINT8)
 */
void subsample3(void *in, void *out, int dims[3], int dims_samp[3], int datatype)
{
    int nvox_samp, nvox;
    float *vol_samp, *buffer;

    nvox = dims[0] * dims[1] * dims[2];
    nvox_samp = dims_samp[0] * dims_samp[1] * dims_samp[2];

    vol_samp = (float *)malloc(sizeof(float) * nvox_samp);
    buffer = (float *)malloc(sizeof(float) * nvox);

    /* check success of memory allocation */
    if (!buffer || !vol_samp)
    {
        printf("Memory allocation error\n");
        exit(EXIT_FAILURE);
    }

    convert_input_type_float(in, buffer, nvox, datatype);
    subsample_float(buffer, vol_samp, dims, dims_samp);
    convert_output_type_float(out, vol_samp, nvox_samp, datatype);

    free(vol_samp);
    free(buffer);
}

void smooth_subsample_float(float *vol, int dims[3], double voxelsize[3],
                            double s[3], int use_mask, double samp_voxelsize)
{
    int i, nvox_samp, nvox;
    int dims_samp[3];
    float *vol_samp;
    double voxelsize_samp[3];

    /* Estimate target dimensions by ensuring a minimum of 32 */
    estimate_target_dimensions(dims, voxelsize, samp_voxelsize, 32, dims_samp);

    /* Define new voxel size */
    for (i = 0; i < 3; i++)
    {
        voxelsize_samp[i] = (double)dims[i] / (double)dims_samp[i] * voxelsize[i];
    }

    nvox = dims[0] * dims[1] * dims[2];
    nvox_samp = dims_samp[0] * dims_samp[1] * dims_samp[2];

    /* Allocate memory for downsampled volume */
    vol_samp = (float *)calloc(nvox_samp, sizeof(float)); // Zero initialization to avoid artifacts
    if (!vol_samp)
    {
        fprintf(stderr, "Memory allocation error\n");
        exit(EXIT_FAILURE);
    }

    /* Subsample to lower resolution */
    subsample_float(vol, vol_samp, dims, dims_samp);

    /* Apply smoothing */
    smooth_float(vol_samp, dims_samp, voxelsize_samp, s, use_mask);

    /* Upsample back to original resolution */
    subsample_float(vol_samp, vol, dims_samp, dims);

    /* Free allocated memory */
    free(vol_samp);
}

void smooth_subsample3(void *data, int dims[3], double voxelsize[3], double s[3],
                       int use_mask, double samp_voxelsize, int datatype)
{
    int i, nvox_samp, nvox;
    int dims_samp[3];
    float *vol_samp;
    float *buffer;
    double voxelsize_samp[3];

    /* Estimate target dimensions by ensuring a minimum of 32 */
    estimate_target_dimensions(dims, voxelsize, samp_voxelsize, 32, dims_samp);

    /* Define new voxel size and compute smoothing size */
    for (i = 0; i < 3; i++)
    {
        voxelsize_samp[i] = (double)dims[i] / (double)dims_samp[i] * voxelsize[i];
    }

    nvox = dims[0] * dims[1] * dims[2];
    nvox_samp = dims_samp[0] * dims_samp[1] * dims_samp[2];
    vol_samp = (float *)calloc(nvox_samp, sizeof(float));
    buffer = (float *)calloc(nvox, sizeof(float));

    /* check success of memory allocation */
    if (!buffer || !vol_samp)
    {
        printf("Memory allocation error\n");
        exit(EXIT_FAILURE);
    }

    convert_input_type_float(data, buffer, nvox, datatype);

    /* Subsample to lower resolution */
    subsample_float(buffer, vol_samp, dims, dims_samp);

    /* Apply smoothing */
    smooth_float(vol_samp, dims_samp, voxelsize_samp, s, use_mask);

    /* Upsample back to original resolution and data type*/
    subsample_float(vol_samp, buffer, dims_samp, dims);
    convert_output_type_float(data, buffer, nvox, datatype);

    /* Free allocated memory */
    free(vol_samp);
    free(buffer);
}
/**
 * \brief Public API for median_subsample3.
 *
 * This function is part of the CAT-Surface public library interface and is used by command-line tools.
 *
 * \param data (in/out) Parameter of median_subsample3.
 * \param dims (in/out) Parameter of median_subsample3.
 * \param voxelsize (in/out) Parameter of median_subsample3.
 * \param niter (in/out) Parameter of median_subsample3.
 * \param samp_voxelsize (in/out) Parameter of median_subsample3.
 * \param datatype (in/out) Parameter of median_subsample3.
 * \return void (no return value).
 */

void median_subsample3(void *data, int dims[3], double voxelsize[3], int niter, double samp_voxelsize, int datatype)
{
    int i, nvox_samp, nvox;
    int dims_samp[3];
    float *vol_samp;
    float *buffer;

    /* Estimate target dimensions by ensuring a minimum of 32 */
    estimate_target_dimensions(dims, voxelsize, samp_voxelsize, 32, dims_samp);

    nvox = dims[0] * dims[1] * dims[2];
    nvox_samp = dims_samp[0] * dims_samp[1] * dims_samp[2];
    vol_samp = (float *)calloc(nvox_samp, sizeof(float));
    buffer = (float *)calloc(nvox, sizeof(float));

    /* check success of memory allocation */
    if (!buffer || !vol_samp)
    {
        printf("Memory allocation error\n");
        exit(EXIT_FAILURE);
    }

    convert_input_type_float(data, buffer, nvox, datatype);

    /* Subsample to lower resolution */
    subsample_float(buffer, vol_samp, dims, dims_samp);

    /* Apply median */
    localstat3(vol_samp, NULL, dims_samp, 1, F_MEDIAN, niter, 0, DT_FLOAT32);

    /* Upsample back to original resolution and data type*/
    subsample_float(vol_samp, buffer, dims_samp, dims);
    convert_output_type_float(data, buffer, nvox, datatype);

    /* Free allocated memory */
    free(vol_samp);
    free(buffer);
}

/**
 * \brief Bias correction on MRI images based on tissue labels.
 *
 * Corrects multiplicative bias field in MRI images using tissue-specific label information.
 * Computes mean intensities for each tissue class, then estimates and removes the bias field
 * based on the ratio between actual and mean values. Different tissues can receive different
 * corrections based on label thresholding.
 *
 * Algorithm:
 *  1. Compute mean intensity for each tissue class from labeled regions
 *  2. Estimate bias field as local ratio of actual intensity to class mean
 *  3. Refine brain mask using morphological operations
 *  4. Smooth bias field with specified FWHM Gaussian kernel
 *  5. Divide source image by the estimated bias field
 *
 * \param src       (in/out) float[nvox]; source image, modified in-place with bias correction
 * \param biasfield (out)    float[nvox]; estimated bias field (can be NULL)
 * \param label     (in)     unsigned char[nvox]; tissue label map (e.g., CSF=1, GM=2, WM=3)
 * \param dims      (in)     {nx, ny, nz} volume dimensions
 * \param voxelsize (in)     {sx, sy, sz} voxel spacing in mm
 * \param bias_fwhm (in)     FWHM of Gaussian smoothing kernel for bias field (mm)
 * \param label_th  (in)     label threshold for selective correction (e.g., 2 for WM only)
 */
void correct_bias_label(float *src, float *biasfield, unsigned char *label, int *dims, double *voxelsize, double bias_fwhm, int label_th)
{
    int i, j, nvox, nvoxr, n[MAX_NC], n_classes = 0;
    int dimsr[3];
    unsigned char *mask, *maskr;
    float *biasfieldr;
    double mean_label[MAX_NC], mean_bias, voxelsizer[3];
    double fwhm[3] = {bias_fwhm, bias_fwhm, bias_fwhm};

    // downsample resolution by factor samp to increase speed
    int samp = 2;

    int set_bias_to_zero = 1;

    // define grid dimensions
    for (i = 0; i < 3; i++)
    {
        dimsr[i] = (int)ceil((dims[i] - 1) / ((double)samp)) + 1;
        voxelsizer[i] = voxelsize[i] * samp;
    }

    nvox = dims[0] * dims[1] * dims[2];
    nvoxr = dimsr[0] * dimsr[1] * dimsr[2];

    biasfieldr = (float *)malloc(sizeof(float) * nvoxr);
    mask = (unsigned char *)malloc(sizeof(unsigned char) * nvox);
    maskr = (unsigned char *)malloc(sizeof(unsigned char) * nvoxr);

    if (!mask || !maskr || !biasfieldr)
    {
        printf("Memory allocation error\n");
        exit(EXIT_FAILURE);
    }

    /* get number of classes by checking maximum label value */
    n_classes = get_max(label, nvox, 0, DT_UINT8);

    /* initialize parameters */
    for (i = 0; i < n_classes; i++)
    {
        n[i] = 0;
        mean_label[i] = 0.0;
    }

    /* estimate mean for each label class */
    for (i = 0; i < nvox; i++)
    {
        if ((label[i] == 0) || isnan(src[i]) || !isfinite(src[i]))
            continue;
        n[label[i] - 1]++;
        mean_label[label[i] - 1] += (double)src[i];
    }
    for (i = 0; i < n_classes; i++)
        mean_label[i] /= (double)n[i];

    /* only use defined labels (i.e. using label_th) for bias estimation
     * use label_th = 2 for focussing on WM only */
    for (i = 0; i < nvox; i++)
        mask[i] = (label[i] >= label_th) ? 1 : 0;

    /* get bias field by ratio between actual values and respective mean of label class */
    for (i = 0; i < nvox; i++)
        if (mask[i] > 0)
            biasfield[i] += (src[i] / (float)mean_label[label[i] - 1]);
        else
            biasfield[i] = 0.0;

    /* downsample to lower resolution */
    subsample3(biasfield, biasfieldr, dims, dimsr, DT_FLOAT32);
    subsample3(mask, maskr, dims, dimsr, DT_UINT8);

    // we have to ensure that thin gyri will be still covered if we use WM only
    // which later makes use of min/max values also necessary
    if (label_th > GM)
    {
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
    for (i = 0; i < nvoxr; i++)
    {
        if (set_bias_to_zero == 1)
            if (maskr[i] == 0)
                biasfieldr[i] = 0.0;
        maskr[i] = 1 - maskr[i];
    }

    vol_approx(biasfieldr, dimsr, voxelsizer);
    //    euclidean_distance(biasfieldr, maskr, dimsr, NULL, 1);
    smooth3(biasfieldr, dimsr, voxelsizer, fwhm, 0, DT_FLOAT32);

    /* upsample to original resolution */
    subsample3(biasfieldr, biasfield, dimsr, dims, DT_FLOAT32);

    /* estimate mean of bias field inside label for mean-correction */
    mean_bias = get_masked_mean_array(biasfield, nvox, label, DT_FLOAT32);

    for (i = 0; i < nvox; i++)
    {
        biasfield[i] /= mean_bias;
        if (biasfield[i] != 0)
            src[i] /= biasfield[i];
    }

    free(mask);
    free(maskr);
    free(biasfieldr);
}

/**
 * \brief Adaptive bias correction for MRI images with optional subcortical refinement.
 *
 * Performs bias field estimation and correction on MRI images with two stages:
 * WM (white matter) correction followed by optional GM (gray matter) correction
 * using local adaptive segmentation for subcortical regions.
 *
 * Algorithm:
 *  1. Estimate and correct WM bias field
 *  2. If weight_las > 0:
 *     - Estimate and correct GM bias field with light smoothing
 *     - Create distance map to subcortical/central brain regions
 *     - Blend WM and GM corrections using distance weighting
 *
 * \param src       (in/out) float[nvox]; source image, modified in-place with bias correction
 * \param biasfield (out)    float[nvox]; estimated bias field (can be NULL)
 * \param label     (in)     unsigned char[nvox]; tissue label map (CSF=1, GM=2, WM=3, etc.)
 * \param dims      (in)     {nx, ny, nz} volume dimensions
 * \param voxelsize (in)     {sx, sy, sz} voxel spacing in mm
 * \param bias_fwhm (in)     FWHM of Gaussian smoothing kernel for WM correction (mm)
 * \param weight_las (in)    weight for local adaptive segmentation GM correction (0..1);
 *                            0 = WM only, >0 = blend WM and GM with distance weighting
 */
void correct_bias(float *src, float *biasfield, unsigned char *label, int *dims, double *voxelsize, double bias_fwhm, double weight_las)
{
    int i, nvox;
    unsigned char *mask;
    float *src_subcortical, *dist, *buffer;
    float max_dist = -FLT_MAX, mx_image, scl;

    nvox = dims[0] * dims[1] * dims[2];

    /* allocate biasfield if necessary */
    if (!biasfield)
    {
        biasfield = (float *)malloc(sizeof(float) * nvox);
        if (!biasfield)
        {
            printf("Memory allocation error\n");
            exit(EXIT_FAILURE);
        }
    }

    /* apply bias correction for WM only */
    correct_bias_label(src, biasfield, label, dims, voxelsize, bias_fwhm, WM);

    /* use local adaptive segmentation (LAS) and apply additional GM correction with
     * very small smoothing */
    if (weight_las > 0.0)
    {
        src_subcortical = (float *)malloc(sizeof(float) * nvox);
        dist = (float *)malloc(sizeof(float) * nvox);
        buffer = (float *)malloc(sizeof(float) * nvox);

        if (!src_subcortical || !dist || !buffer)
        {
            printf("Memory allocation error\n");
            exit(EXIT_FAILURE);
        }

        /* apply bias correction for GM */
        for (i = 0; i < nvox; i++)
            src_subcortical[i] = src[i];
        bias_fwhm = 3.0;
        correct_bias_label(src_subcortical, buffer, label, dims, voxelsize, bias_fwhm, GM);

        /* weight LAS correction */
        for (i = 0; i < nvox; i++)
            src[i] = weight_las * src_subcortical[i] + (1.0 - weight_las) * src[i];

        /* prepare mask for distance map */
        for (i = 0; i < nvox; i++)
            dist[i] = (label[i] > 0) ? 0.0 : 1.0;

        /* get distance to background */
        euclidean_distance(dist, NULL, dims, NULL, 0);

        /* scale distance values to 0..1  */
        max_dist = get_max(dist, nvox, 0, DT_FLOAT32);
        for (i = 0; i < nvox; i++)
            dist[i] /= max_dist;

        /* use cubic distance weights to weight central regions (i.e. subcortical
           structures) more */
        for (i = 0; i < nvox; i++)
            dist[i] *= dist[i] * dist[i];

        /* apply weighted average (with squared weights) to maximize weighting
         * for subcortical regions with large distances, otherwise WM-bias
         * correction is more weighted */
        for (i = 0; i < nvox; i++)
            src[i] = dist[i] * src_subcortical[i] + (1.0 - dist[i]) * src[i];

        free(buffer);
        free(dist);
        free(src_subcortical);
    }
}

/**
 * \brief Approximate missing values in a volume by interpolating from neighbors.
 *
 * This function uses the Laplace approximation to fill zero and negative values
 * by interpolating from their neighboring non-zero voxels. Internally, it applies
 * a smoothing/interpolation procedure via euclidean_distance() and laplace3R() to
 * estimate reasonable values for interior voxels that have been zeroed.
 *
 * \param vol       (in/out) float[dims[0]*dims[1]*dims[2]]; modified in place
 * \param dims      (in)     {nx, ny, nz}
 * \param voxelsize (in)     voxel spacing in mm (or consistent units)
 */
void vol_approx(float *vol, int dims[3], double voxelsize[3])
{
    int i, nvox;
    float *buffer, *TAr;
    float min_vol;
    unsigned char *mask, *mask2;

    nvox = dims[0] * dims[1] * dims[2];
    buffer = (float *)malloc(sizeof(float) * nvox);
    TAr = (float *)malloc(sizeof(float) * nvox);
    mask = (unsigned char *)malloc(sizeof(unsigned char) * nvox);
    mask2 = (unsigned char *)malloc(sizeof(unsigned char) * nvox);

    /* The function approximates only values larger than zeros. The easiest was
       to shift the value and use a mask to redefine the zero values. */
    min_vol = get_min(vol, nvox, 0, DT_FLOAT32);
    if (min_vol < 0)
    {
        for (i = 0; i < nvox; ++i)
            vol[i] = (vol[i] == 0.0) ? 0.0 : vol[i] - min_vol + 1.0;
    }

    /* euclidean_distance to fill values in background with neighbours */
    memcpy(buffer, vol, nvox * sizeof(float));
    euclidean_distance(buffer, NULL, dims, NULL, 1);
    for (i = 0; i < nvox; ++i)
        TAr[i] = (vol[i] > 0) ? vol[i] : buffer[i];

    /* create mask by closing holes */
    for (i = 0; i < nvox; ++i)
        mask[i] = vol[i] > 0;
    dist_close(mask, dims, voxelsize, 50.0, 0.0, DT_UINT8);

    /* smooth values outside mask */
    memcpy(buffer, TAr, nvox * sizeof(float));
    double s[3] = {20, 20, 20};
    smooth_subsample_float(buffer, dims, voxelsize, s, 1, 4.0);
    for (i = 0; i < nvox; ++i)
        if (mask[i] == 0)
            TAr[i] = buffer[i];

    /* rescue mask and create new mask that is only defined inside */
    for (i = 0; i < nvox; ++i)
        mask2[i] = (mask[i] > 0) && (vol[i] == 0);

    double laplace_thresh = 0.4;
    laplace3R(TAr, mask2, dims, laplace_thresh);
    median_subsample3(TAr, dims, voxelsize, 1, 2.0, DT_FLOAT32);
    laplace3R(TAr, mask2, dims, laplace_thresh);

    /* only keep TAr inside (closed) mask */
    for (i = 0; i < nvox; ++i)
        if ((mask[i] == 0) && (vol[i] == 0))
            TAr[i] = 0.0;

    /* again apply euclidean_distance to fill values in background with neighbours */
    euclidean_distance(TAr, NULL, dims, NULL, 1);

    memcpy(buffer, TAr, nvox * sizeof(float));
    for (i = 0; i < 3; ++i)
        s[i] = 40.0;
    smooth_subsample_float(buffer, dims, voxelsize, s, 1, 4.0);
    for (i = 0; i < nvox; ++i)
        if (mask[i] == 0)
            TAr[i] = buffer[i];
    for (i = 0; i < nvox; ++i)
        mask[i] = (mask[i] == 0);
    laplace3R(TAr, mask, dims, laplace_thresh);
    median_subsample3(TAr, dims, voxelsize, 1, 2.0, DT_FLOAT32);

    for (i = 0; i < nvox; ++i)
        mask[i] = (vol[i] == 0);
    laplace3R(TAr, mask, dims, laplace_thresh);

    memcpy(vol, TAr, nvox * sizeof(float));

    free(buffer);
    free(mask);
    free(mask2);
    free(TAr);
}

/**
 * \brief Clean up tissue probability map by morphological refinement.
 *
 * Refines a tissue probability array (containing CSF, GM, WM probabilities) by iteratively
 * applying erosion and conditional dilation operations. This cleanup process improves
 * tissue classification by removing small isolated regions and bridging small gaps,
 * while preserving thin structures (e.g., cerebellum, deep brain structures) using
 * distance-weighted operations.
 *
 * Algorithm:
 *  1. Extract WM and GM probability maps from input array
 *  2. Create distance map to tissue boundary (to weight central regions higher)
 *  3. Initial erosion phase (2 iterations) to remove noise
 *  4. Conditional dilation phase to restore shape while maintaining tissue separation
 *  5. Prevent false gaps between hemispheres using morphological closing
 *  6. Update probability map in-place with cleaned tissue probabilities
 *
 * \param prob     (in/out) unsigned char[3*nvox]; tissue probability array
 *                          [0:nvox-1]=CSF, [nvox:2*nvox-1]=GM, [2*nvox:3*nvox-1]=WM
 *                          Modified in-place by cleanup operations
 * \param dims     (in)     {nx, ny, nz} volume dimensions
 * \param voxelsize (in)     {sx, sy, sz} voxel spacing in mm; used for morphological scaling
 * \param strength  (in)     cleanup strength (0..N); controls dilation threshold
 *                           (higher = more aggressive cleanup)
 */
void cleanup_brain(unsigned char *prob, int dims[3], double voxelsize[3], int strength)
{
    double scale = 3.0 / (voxelsize[0] + voxelsize[1] + voxelsize[2]);
    int niter, iter, nvox, th, th_erode, th_dilate, i, j, gm_offset, wm_offset, csf_offset;
    float *wm, *gm, *dist;
    float bp, tot, max_dist;
    double filt[3] = {0.75, 1.0, 0.75};

    niter = 32;
    th_erode = 153;                   // initial threshold for erosion 0.6*255.0
    th_dilate = 38 + (strength * 13); // threshold for dilation (0.15 + strength*0.05)*255
    th_dilate = 38 + (strength * 13); // threshold for dilation (0.15 + strength*0.05)*255

    nvox = dims[0] * dims[1] * dims[2];

    // we need the offsets for the tissue classes for prob array
    gm_offset = (int)(GM - 1) * nvox;
    wm_offset = (int)(WM - 1) * nvox;
    csf_offset = (int)(CSF - 1) * nvox;

    /* ensure that sum of filt is 1 */
    tot = 0.0;
    for (i = 0; i < 3; ++i)
        tot += filt[i];
    for (i = 0; i < 3; ++i)
        filt[i] /= tot;

    wm = (float *)malloc(sizeof(float) * nvox);
    gm = (float *)malloc(sizeof(float) * nvox);
    dist = (float *)malloc(sizeof(float) * nvox);

    if (!wm || !gm || !dist)
    {
        printf("Memory allocation error\n");
        exit(EXIT_FAILURE);
    }

    /* init with WM and GM */
    for (i = 0; i < nvox; ++i)
    {
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
    euclidean_distance(dist, NULL, dims, NULL, 0);

    /* scale distance to 0..1  */
    max_dist = get_max(dist, nvox, 0, DT_FLOAT32);
    for (i = 0; i < nvox; i++)
        dist[i] /= max_dist;

    /* use cubic distance weights to weight central regions (i.e. subcortical
       structures) more */
    for (i = 0; i < nvox; i++)
        dist[i] *= dist[i];

    /* mask computed from gm and wm */
    /* erosions and conditional dilations */
    for (iter = 0; iter < niter; iter++)
    {

        /*  start with 2 iterations of erosions */
        th = (iter < 2) ? th_erode : th_dilate;

        /* b = (b>th).*(white+gray) */
        for (i = 0; i < nvox; ++i)
        {
            bp = (float)prob[i + gm_offset] + (float)prob[i + wm_offset];
            wm[i] = (wm[i] > (float)th) ? MIN(bp, 255.0) : 0.0;
        }

        /* convolve mask with filter width of 3 voxel */
        convxyz_float(wm, filt, filt, filt, 3, 3, 3, -1, -1, -1, wm, dims);
    }

    // apply opening, but keep original values otherwise
    morph_open_float(wm, dims, MAX(1, roundf(scale * strength)), 0.1, 1);

    th = 13; // threshold for final cleanup 0.05*255
    for (i = 0; i < nvox; ++i)
    {
        bp = (float)prob[i + gm_offset] + (float)prob[i + wm_offset];
        bp = ((bp > (float)th && wm[i] > (float)th)) ? 1.0 : 0.0;

        // rescue original segmentation for average weighting
        wm[i] = (float)prob[i + wm_offset];
        gm[i] = (float)prob[i + gm_offset];

        // apply cleanup
        if (bp == 0)
        {
            prob[i + gm_offset] = 0;
            prob[i + wm_offset] = 0;
        }
    }

    // use weighted average w.r.t distance to weight cleaned and original segmentation
    for (i = 0; i < nvox; i++)
    {
        prob[i + gm_offset] = (unsigned char)roundf(dist[i] * gm[i] + (1.0 - dist[i]) * (float)prob[i + gm_offset]);
        prob[i + wm_offset] = (unsigned char)roundf(dist[i] * wm[i] + (1.0 - dist[i]) * (float)prob[i + wm_offset]);
    }

    // ensure that overall sum of all tissue classes per voxel is 255 for uint8
    for (i = 0; i < nvox; ++i)
    {
        tot = 0.0;
        for (j = 0; j < 3; ++j)
            tot += (float)prob[i + j * nvox];
        for (j = 0; j < 3; ++j)
            prob[i + j * nvox] = (unsigned char)(roundf((float)prob[i + j * nvox] / tot * 255.0));
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
 *          have 'nvox' elements. This array is modified in place, with only the largest
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
 * conn: Integer value that sets the connection-scheme to 6, 18 or 26 neighbours.
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
    int nvox = dims[0] * dims[1] * dims[2];
    char *flagUsed;
    short *growing;

    if ((conn != 6) && (conn != 18) && (conn != 26))
    {
        printf("Values for connectivity can be only 6, 18 or 26.\n");
        exit(EXIT_FAILURE);
    }

    int adiff = 1; // connectivity 6: faces only
    if (conn == 18)
        adiff = 2; // connectivity 18: faces+edges
    if (conn == 26)
        adiff = 3; // connectivity 26: faces+edges+corners

    flagUsed = (char *)malloc(nvox * sizeof(char));
    outData = (float *)malloc(nvox * sizeof(float));
    growing = (short *)malloc(nvox * 3 * sizeof(short));

    /* check success of memory allocation */
    if (!flagUsed || !outData || !growing)
    {
        printf("Memory allocation error\n");
        exit(EXIT_FAILURE);
    }

    for (i = 0; i < nvox; ++i)
    {
        flagUsed[i] = 0;
        outData[i] = 0.0;
    }

    for (k = 0; k < dims[2]; ++k)
        for (j = 0; j < dims[1]; ++j)
            for (i = 0; i < dims[0]; ++i)
            {
                ind = k * (dims[0] * dims[1]) + (j * dims[0]) + i;

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
                        maxi = MIN(dims[0], growing[growingCur] + 2);
                        maxj = MIN(dims[1], growing[growingCur + 1] + 2);
                        maxk = MIN(dims[2], growing[growingCur + 2] + 2);

                        mini = MAX(0, growing[growingCur] - 1);
                        minj = MAX(0, growing[growingCur + 1] - 1);
                        mink = MAX(0, growing[growingCur + 2] - 1);

                        for (tk = mink; tk < maxk; ++tk)
                            for (tj = minj; tj < maxj; ++tj)
                                for (ti = mini; ti < maxi; ++ti)
                                {
                                    /*  Check for corners and edges depending on defined connectivity */
                                    if (abs(ti - growing[growingCur]) + abs(tj - growing[growingCur + 1]) + abs(tk - growing[growingCur + 2]) > adiff)
                                        continue;
                                    ind1 = tk * (dims[0] * dims[1]) + (tj * dims[0]) + ti;

                                    if (!flagUsed[ind1] && inData[ind1] >= thresh)
                                    {
                                        flagUsed[ind1] = 1;
                                        growing[growingInd] = ti;
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
                        outData[growing[growingCur + 2] * (dims[0] * dims[1]) + (growing[growingCur + 1] * dims[0]) + growing[growingCur]] += (float)growingInd;
                        growingCur += 3;
                    }
                }
            }

    /* find maximum value which is the largest cluster or use defined minimum cluster size */
    if (min_size >= 0)
        maxInd = get_max(outData, nvox, 0, DT_FLOAT32);
    else
        maxInd = min_size;

    /* set values with smaller clusters to zero */
    /* Depending on the parameter retain_above_th we either set smaller clusters to zero,
     * but retain all original values otherwise, or we only keep values in larger clusters
     * that are then thresholded */
    for (i = 0; i < nvox; ++i)
    {
        if (retain_above_th)
        {
            if ((outData[i] > 0) && (outData[i] < maxInd))
                inData[i] = 0.0;
        }
        else
            inData[i] = (outData[i] >= maxInd) ? inData[i] : 0;
    }

    free(flagUsed);
    free(growing);
    free(outData);
}

/**
 * keep_largest_cluster - Identifies the largest cluster in a 3D volume.
 *
 * This function calls keep_largest_cluster_float and converts datatypes of
 * input and output. It identifies connected components in a binary or labeled volume
 * and retains only the largest cluster, optionally filtered by size constraints.
 *
 * \param data        (in/out) volume data (pointer to any supported datatype)
 * \param thresh      (in)     threshold for connectivity; voxels above this connect
 * \param dims        (in)     volume dimensions {nx, ny, nz}
 * \param datatype    (in)     data type descriptor (e.g., DT_FLOAT32, DT_UINT8)
 * \param min_size    (in)     minimum cluster size (if <0, retain clusters >= |min_size|)
 * \param retain_above_th (in) if 1, set smaller clusters to 0; if 0, invert logic
 * \param conn        (in)     connectivity: 6 (face), 18 (face+edge), 26 (face+edge+corner)
 */
void keep_largest_cluster(void *data, double thresh, int *dims, int datatype, int min_size, int retain_above_th, int conn)
{
    int nvox;
    float *buffer;

    nvox = dims[0] * dims[1] * dims[2];
    buffer = (float *)malloc(sizeof(float) * nvox);

    /* check success of memory allocation */
    if (!buffer)
    {
        printf("Memory allocation error\n");
        exit(EXIT_FAILURE);
    }

    convert_input_type_float(data, buffer, nvox, datatype);
    keep_largest_cluster_float(buffer, thresh, dims, min_size, retain_above_th, conn);
    convert_output_type_float(data, buffer, nvox, datatype);

    free(buffer);
}

/**
 * \\brief Fill holes in a binary or thresholded volume.
 *
 * Identifies and fills small background regions (holes) in a thresholded volume.
 * Uses connected-component analysis to identify isolated background voxels and
 * replaces them with estimated foreground or specified fill values.
 *
 * Algorithm:
 *  1. Threshold input volume (values < thresh are potential holes)
 *  2. Invert mask to identify connected background components
 *  3. Keep only the largest background cluster (typically representing true background)
 *  4. Fill identified holes with specified value or smooth estimate from neighbors
 *  5. Convert result back to original datatype
 *
 * \\param data      (in/out) void pointer to volume data; type given by datatype parameter
 * \\param dims      (in)     {nx, ny, nz} volume dimensions
 * \\param thresh    (in)     threshold value; voxels < thresh are treated as potential holes
 * \\param fill_value (in)    value to fill holes with;
 *                            if negative, holes are filled with locally estimated values
 *                            if >=0, holes are filled with this fixed value
 * \\param datatype  (in)     data type code (DT_UINT8, DT_UINT16, DT_FLOAT32, etc.)
 */
void fill_holes(void *data, int *dims, double thresh, double fill_value, int datatype)
{
    int i, nvox;
    float *buffer, *mask_inv;
    double voxelsize[3] = {1.0, 1.0, 1.0};
    unsigned char *mask_fill;

    nvox = dims[0] * dims[1] * dims[2];
    buffer = (float *)malloc(sizeof(float) * nvox);
    mask_inv = (float *)malloc(sizeof(float) * nvox);
    mask_fill = (unsigned char *)malloc(sizeof(unsigned char) * nvox);

    /* check success of memory allocation */
    if (!buffer || !mask_fill || !mask_inv)
    {
        printf("Memory allocation error\n");
        exit(EXIT_FAILURE);
    }

    convert_input_type_float(data, buffer, nvox, datatype);

    /* get inverted mask after thresholding */
    for (i = 0; i < nvox; ++i)
        mask_inv[i] = (buffer[i] >= thresh) ? 0.0 : 1.0;

    /* retain largest cluster of inverted mask, which should be the background */
    keep_largest_cluster_float(mask_inv, thresh, dims, 0, 1, 18);

    /* fill those values (=holes) that were removed by the previous keep_largest_cluster step */
    for (i = 0; i < nvox; ++i)
        mask_fill[i] = ((mask_inv[i] == 0.0) && (buffer[i] < thresh)) ? 1 : 0;

    if (fill_value < 0.0)
    {
        euclidean_distance(buffer, mask_fill, dims, NULL, 1);

        /* ensure a minimum filled value that is to the threshold*1.001 */
        for (i = 0; i < nvox; ++i)
            buffer[i] = ((mask_fill[i] == 1) && (buffer[i] <= thresh)) ? thresh * 1.001 : buffer[i];
    }
    else
    {
        for (i = 0; i < nvox; ++i)
            buffer[i] = (mask_fill[i] == 1) ? fill_value : buffer[i];
    }

    convert_output_type_float(data, buffer, nvox, datatype);

    free(buffer);
    free(mask_inv);
    free(mask_fill);
}

/**
 * Compute x-gradient of 3D volume
 */
float gradientX(float *src, int i, int j, int k, int dims[3], double voxelsize[3])
{
    int index = sub2ind(i, j, k, dims);
    float dx = voxelsize[0]; // X-axis voxel size

    if (i > 0 && i < dims[0] - 1)
        return (src[sub2ind(i + 1, j, k, dims)] - src[sub2ind(i - 1, j, k, dims)]) / (2.0 * dx);
    else if (i == 0)
        return (src[sub2ind(i + 1, j, k, dims)] - src[index]) / dx;
    else // i == dims[0] - 1
        return (src[index] - src[sub2ind(i - 1, j, k, dims)]) / dx;
}

/**
 * Compute y-gradient of 3D volume
 */
float gradientY(float *src, int i, int j, int k, int dims[3], double voxelsize[3])
{
    int index = sub2ind(i, j, k, dims);
    float dy = voxelsize[1]; // Y-axis voxel size

    if (j > 0 && j < dims[1] - 1)
        return (src[sub2ind(i, j + 1, k, dims)] - src[sub2ind(i, j - 1, k, dims)]) / (2.0 * dy);
    else if (j == 0)
        return (src[sub2ind(i, j + 1, k, dims)] - src[index]) / dy;
    else // j == dims[1] - 1
        return (src[index] - src[sub2ind(i, j - 1, k, dims)]) / dy;
}

/**
 * Compute z-gradient of 3D volume
 */
float gradientZ(float *src, int i, int j, int k, int dims[3], double voxelsize[3])
{
    int index = sub2ind(i, j, k, dims);
    float dz = voxelsize[2]; // Z-axis voxel size

    if (k > 0 && k < dims[2] - 1)
        return (src[sub2ind(i, j, k + 1, dims)] - src[sub2ind(i, j, k - 1, dims)]) / (2.0 * dz);
    else if (k == 0)
        return (src[sub2ind(i, j, k + 1, dims)] - src[index]) / dz;
    else // k == dims[2] - 1
        return (src[index] - src[sub2ind(i, j, k - 1, dims)]) / dz;
}

/**
 * \brief Compute local gradient magnitude and components for a 3D volume.
 *
 * Calculates the spatial gradient (first derivatives) at each voxel using
 * finite differences, with boundary-aware handling at volume edges.
 * Optionally outputs the gradient magnitude and individual x, y, z components.
 *
 * \param src       (in)  input volume float[dims[0]*dims[1]*dims[2]]
 * \param grad_mag  (out) gradient magnitude; NULL to skip (optional)
 * \param grad_x    (out) x-component of gradient; NULL to skip (optional)
 * \param grad_y    (out) y-component of gradient; NULL to skip (optional)
 * \param grad_z    (out) z-component of gradient; NULL to skip (optional)
 * \param dims      (in)  volume dimensions {nx, ny, nz}
 * \param voxelsize (in)  voxel spacing in mm {dx, dy, dz}
 *
 * \note Boundary voxels use one-sided differences to avoid out-of-bounds access.
 */
void gradient3D(float *src, float *grad_mag, float *grad_x, float *grad_y,
                float *grad_z, int dims[3], double voxelsize[3])
{
    int i, j, k, index;
    float gradx, grady, gradz;

    for (i = 0; i < dims[0]; ++i)
    {
        for (j = 0; j < dims[1]; ++j)
        {
            for (k = 0; k < dims[2]; ++k)
            {
                index = sub2ind(i, j, k, dims);
                gradx = gradientX(src, i, j, k, dims, voxelsize);
                grady = gradientY(src, i, j, k, dims, voxelsize);
                gradz = gradientZ(src, i, j, k, dims, voxelsize);

                /* Also output xyz-gradients */
                if (grad_x)
                    grad_x[index] = gradx;
                if (grad_y)
                    grad_y[index] = grady;
                if (grad_z)
                    grad_z[index] = gradz;

                /* Also output gradient magnitude */
                if (grad_mag)
                    grad_mag[index] = sqrtf(gradx * gradx + grady * grady + gradz * gradz);
            }
        }
    }
}
