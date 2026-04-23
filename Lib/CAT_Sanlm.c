/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 *
 *
 *
 * This code is a modified version of MABONLM3D.c
 * Jose V. Manjon - jmanjon@fis.upv.es
 * Pierrick Coupe - pierrick.coupe@gmail.com
 * Brain Imaging Center, Montreal Neurological Institute.
 * Mc Gill University
 *
 * Copyright (C) 2010 Jose V. Manjon and Pierrick Coupe

 ***************************************************************************
 * Adaptive Non-Local Means Denoising of MR Images
 * With Spatially Varying Noise Levels
 *
 * Jose V. Manjon, Pierrick Coupe, Luis Marti-Bonmati,
 * D. Louis Collins and Montserrat Robles
 ***************************************************************************
 *
 * Details on SANLM filter
 ***************************************************************************
 *    The SANLM filter is described in:
 *
 *    Jose V. Manjon, Pierrick Coupe, Luis Marti-Bonmati, Montserrat Robles
 *    and D. Louis Collins.
 *    Adaptive Non-Local Means Denoising of MR Images with Spatially Varying
 *    Noise Levels. Journal of Magnetic Resonance Imaging, 31,192-203, 2010.
 *
 ***************************************************************************/

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

/* Multithreading stuff */
#if !defined(_WIN32)
#include <pthread.h>
#endif

#define PI 3.1415926535
#define MAX_NTHREADS 16

#if !defined(_WIN32)
pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;
#endif

typedef struct
{
    int rows;
    int cols;
    int slices;
    float *in_image;
    float *means_image;
    float *var_image;
    float *estimate;
    float *bias;
    unsigned char *label;
    int ini;
    int fin;
    int radioB;
    int radioS;
    double strength;
} myargument;

int rician;
double max;

/*Returns the modified Bessel function I0(x) for any real x.*/
/**
 * \brief Compute modified Bessel function I0(x).
 *
 * Uses a polynomial approximation for small x and an asymptotic form
 * for large x.
 *
 * \param x (in) input value
 * \return I0(x)
 */
double bessi0(double x)
{
    double ax, res, a;
    double y;
    if ((ax = fabs(x)) < 3.75)
    {
        y = x / 3.75;
        y *= y;
        res = 1.0 + y * (3.5156229 + y * (3.0899424 + y * (1.2067492 + y * (0.2659732 + y * (0.360768e-1 + y * 0.45813e-2)))));
    }
    else
    {
        y = 3.75 / ax;
        res = (exp(ax) / sqrt(ax));
        a = y * (0.916281e-2 + y * (-0.2057706e-1 + y * (0.2635537e-1 + y * (-0.1647633e-1 + y * 0.392377e-2))));
        res = res * (0.39894228 + y * (0.1328592e-1 + y * (0.225319e-2 + y * (-0.157565e-2 + a))));
    }
    return res;
}

/*Returns the modified Bessel function I1(x) for any real x.*/
/**
 * \brief Compute modified Bessel function I1(x).
 *
 * Uses a polynomial approximation for small x and an asymptotic form
 * for large x.
 *
 * \param x (in) input value
 * \return I1(x)
 */
double bessi1(double x)
{
    double ax, res;
    double y;
    if ((ax = fabs(x)) < 3.75)
    {
        y = x / 3.75;
        y *= y;
        res = ax * (0.5 + y * (0.87890594 + y * (0.51498869 + y * (0.15084934 + y * (0.2658733e-1 + y * (0.301532e-2 + y * 0.32411e-3))))));
    }
    else
    {
        y = 3.75 / ax;
        res = 0.2282967e-1 + y * (-0.2895312e-1 + y * (0.1787654e-1 - y * 0.420059e-2));
        res = 0.39894228 + y * (-0.3988024e-1 + y * (-0.362018e-2 + y * (0.163801e-2 + y * (-0.1031555e-1 + y * res))));
        res *= (exp(ax) / sqrt(ax));
    }
    return x < 0.0 ? -res : res;
}

/**
 * \brief Compute Rician correction factor Epsi for a given SNR.
 *
 * \param snr (in) signal-to-noise ratio
 * \return Rician correction factor
 */
double Epsi(double snr)
{
    double val;
    val = 2 + snr * snr - (PI / 8) * exp(-(snr * snr) / 2) * ((2 + snr * snr) * bessi0((snr * snr) / 4) + (snr * snr) * bessi1((snr * snr) / 4)) * ((2 + snr * snr) * bessi0((snr * snr) / 4) + (snr * snr) * bessi1((snr * snr) / 4));
    if (val < 0.001)
        val = 1;
    if (val > 10)
        val = 1;
    return val;
}

/* Function which compute the weighted average for one block */
/**
 * \brief Accumulate weighted block values for SANLM filtering.
 *
 * Updates the block accumulator with values from a 3D neighborhood.
 *
 * \param ima             (in)  input image volume
 * \param x               (in)  block center x index
 * \param y               (in)  block center y index
 * \param z               (in)  block center z index
 * \param neighborhoodsize (in) half-size of block window
 * \param average         (in/out) block accumulator array
 * \param weight          (in)  weight for this block
 * \param sx              (in)  volume size x
 * \param sy              (in)  volume size y
 * \param sz              (in)  volume size z
 * \return void
 */
void Average_block(float *ima, int x, int y, int z, int neighborhoodsize, float *average, double weight, int sx, int sy, int sz)
{
    int x_pos, y_pos, z_pos;
    int is_outside;
    int a, b, c, ns, sxy, index;

    extern int rician;

    ns = 2 * neighborhoodsize + 1;
    sxy = sx * sy;

    for (c = 0; c < ns; c++)
    {
        for (b = 0; b < ns; b++)
        {
            for (a = 0; a < ns; a++)
            {
                is_outside = 0;
                x_pos = x + a - neighborhoodsize;
                y_pos = y + b - neighborhoodsize;
                z_pos = z + c - neighborhoodsize;

                if ((z_pos < 0) || (z_pos > sz - 1))
                    is_outside = 1;
                if ((y_pos < 0) || (y_pos > sy - 1))
                    is_outside = 1;
                if ((x_pos < 0) || (x_pos > sx - 1))
                    is_outside = 1;

                index = a + ns * b + ns * ns * c;

                if (rician)
                {
                    if (is_outside)
                        average[index] += ima[z * (sxy) + (y * sx) + x] * ima[z * (sxy) + (y * sx) + x] * (float)weight;
                    else
                        average[index] += ima[z_pos * (sxy) + (y_pos * sx) + x_pos] * ima[z_pos * (sxy) + (y_pos * sx) + x_pos] * (float)weight;
                }
                else
                {
                    if (is_outside)
                        average[index] += ima[z * (sxy) + (y * sx) + x] * (float)weight;
                    else
                        average[index] += ima[z_pos * (sxy) + (y_pos * sx) + x_pos] * (float)weight;
                }
            }
        }
    }
}

/* Function which computes the value assigned to each voxel */
/**
 * \brief Update voxel estimates from accumulated SANLM block weights.
 *
 * \param Estimate        (in/out) running estimate accumulator
 * \param Label           (in/out) contribution counts per voxel
 * \param x               (in)  block center x index
 * \param y               (in)  block center y index
 * \param z               (in)  block center z index
 * \param neighborhoodsize (in) half-size of block window
 * \param average         (in)  block accumulator array
 * \param global_sum      (in)  sum of weights for this block
 * \param sx              (in)  volume size x
 * \param sy              (in)  volume size y
 * \param sz              (in)  volume size z
 * \return void
 */
void Value_block(float *Estimate, unsigned char *Label, int x, int y, int z, int neighborhoodsize, float *average, double global_sum, int sx, int sy, int sz)
{
    int x_pos, y_pos, z_pos;
    int is_outside;
    double value = 0.0;
    unsigned char label = 0;
    int count = 0;
    int a, b, c, ns, sxy;

    ns = 2 * neighborhoodsize + 1;
    sxy = sx * sy;

    for (c = 0; c < ns; c++)
    {
        for (b = 0; b < ns; b++)
        {
            for (a = 0; a < ns; a++)
            {
                is_outside = 0;
                x_pos = x + a - neighborhoodsize;
                y_pos = y + b - neighborhoodsize;
                z_pos = z + c - neighborhoodsize;

                if ((z_pos < 0) || (z_pos > sz - 1))
                    is_outside = 1;
                if ((y_pos < 0) || (y_pos > sy - 1))
                    is_outside = 1;
                if ((x_pos < 0) || (x_pos > sx - 1))
                    is_outside = 1;

                if (!is_outside)
                {
                    value = (double)Estimate[z_pos * (sxy) + (y_pos * sx) + x_pos];
                    value += ((double)average[count] / global_sum);

                    label = Label[(x_pos + y_pos * sx + z_pos * sxy)];
                    Estimate[z_pos * (sxy) + (y_pos * sx) + x_pos] = (float)value;
                    Label[(x_pos + y_pos * sx + z_pos * sxy)] = label + 1;
                }
                count++;
            }
        }
    }
}

/**
 * \brief Compute patch distance between two neighborhoods.
 *
 * \param ima (in) input image volume
 * \param x   (in) first patch center x
 * \param y   (in) first patch center y
 * \param z   (in) first patch center z
 * \param nx  (in) second patch center x
 * \param ny  (in) second patch center y
 * \param nz  (in) second patch center z
 * \param f   (in) half-size of patch window
 * \param sx  (in) volume size x
 * \param sy  (in) volume size y
 * \param sz  (in) volume size z
 * \return Average squared distance between patches
 */
double distance(float *ima, int x, int y, int z, int nx, int ny, int nz, int f, int sx, int sy, int sz)
{
    double d, acu, distancetotal;
    int i, j, k, ni1, nj1, ni2, nj2, nk1, nk2;

    distancetotal = 0;
    for (k = -f; k <= f; k++)
    {
        nk1 = z + k;
        nk2 = nz + k;
        if (nk1 < 0)
            nk1 = -nk1;
        if (nk2 < 0)
            nk2 = -nk2;
        if (nk1 >= sz)
            nk1 = 2 * sz - nk1 - 1;
        if (nk2 >= sz)
            nk2 = 2 * sz - nk2 - 1;

        for (j = -f; j <= f; j++)
        {
            nj1 = y + j;
            nj2 = ny + j;
            if (nj1 < 0)
                nj1 = -nj1;
            if (nj2 < 0)
                nj2 = -nj2;
            if (nj1 >= sy)
                nj1 = 2 * sy - nj1 - 1;
            if (nj2 >= sy)
                nj2 = 2 * sy - nj2 - 1;

            for (i = -f; i <= f; i++)
            {
                ni1 = x + i;
                ni2 = nx + i;
                if (ni1 < 0)
                    ni1 = -ni1;
                if (ni2 < 0)
                    ni2 = -ni2;
                if (ni1 >= sx)
                    ni1 = 2 * sx - ni1 - 1;
                if (ni2 >= sx)
                    ni2 = 2 * sx - ni2 - 1;

                d = ((double)ima[nk1 * (sx * sy) + (nj1 * sx) + ni1] - (double)ima[nk2 * (sx * sy) + (nj2 * sx) + ni2]);
                distancetotal += d * d;
            }
        }
    }

    acu = (2 * f + 1) * (2 * f + 1) * (2 * f + 1);
    d = distancetotal / acu;

    return d;
}

/**
 * \brief Compute patch distance after mean removal.
 *
 * Uses the local mean image to compute centered patch distances for
 * adaptive weighting.
 *
 * \param ima   (in) input image volume
 * \param means (in) local mean image
 * \param x     (in) first patch center x
 * \param y     (in) first patch center y
 * \param z     (in) first patch center z
 * \param nx    (in) second patch center x
 * \param ny    (in) second patch center y
 * \param nz    (in) second patch center z
 * \param f     (in) half-size of patch window
 * \param sx    (in) volume size x
 * \param sy    (in) volume size y
 * \param sz    (in) volume size z
 * \return Average squared distance between patches
 */
double distance2(float *ima, float *means, int x, int y, int z, int nx, int ny, int nz, int f, int sx, int sy, int sz)
{
    double d, acu, distancetotal;
    int i, j, k, ni1, nj1, ni2, nj2, nk1, nk2;

    acu = 0;
    distancetotal = 0;

    for (k = -f; k <= f; k++)
    {
        nk1 = z + k;
        nk2 = nz + k;
        if (nk1 < 0)
            nk1 = -nk1;
        if (nk2 < 0)
            nk2 = -nk2;
        if (nk1 >= sz)
            nk1 = 2 * sz - nk1 - 1;
        if (nk2 >= sz)
            nk2 = 2 * sz - nk2 - 1;

        for (j = -f; j <= f; j++)
        {
            nj1 = y + j;
            nj2 = ny + j;
            if (nj1 < 0)
                nj1 = -nj1;
            if (nj2 < 0)
                nj2 = -nj2;
            if (nj1 >= sy)
                nj1 = 2 * sy - nj1 - 1;
            if (nj2 >= sy)
                nj2 = 2 * sy - nj2 - 1;

            for (i = -f; i <= f; i++)
            {
                ni1 = x + i;
                ni2 = nx + i;
                if (ni1 < 0)
                    ni1 = -ni1;
                if (ni2 < 0)
                    ni2 = -ni2;
                if (ni1 >= sx)
                    ni1 = 2 * sx - ni1 - 1;
                if (ni2 >= sx)
                    ni2 = 2 * sx - ni2 - 1;

                d = ((double)ima[nk1 * (sx * sy) + (nj1 * sx) + ni1] - (double)means[nk1 * (sx * sy) + (nj1 * sx) + ni1]) - ((double)ima[nk2 * (sx * sy) + (nj2 * sx) + ni2] - (double)means[nk2 * (sx * sy) + (nj2 * sx) + ni2]);
                distancetotal += d * d;
            }
        }
    }

    acu = (2 * f + 1) * (2 * f + 1) * (2 * f + 1);
    d = distancetotal / acu;

    return d;
}

/**
 * \brief Apply separable box regularization to a volume.
 *
 * Performs a separable averaging with radius r, skipping zero entries.
 *
 * \param in  (in)  input volume
 * \param out (out) output volume
 * \param r   (in)  half-window size
 * \param sx  (in)  volume size x
 * \param sy  (in)  volume size y
 * \param sz  (in)  volume size z
 * \return void
 */
void Regularize(float *in, float *out, int r, int sx, int sy, int sz)
{
    double acu, *temp;
    int ind, i, j, k, ni, nj, nk, ii, jj, kk;

    temp = (double *)malloc(sx * sy * sz * sizeof(double));

    /* separable convolution */
    for (k = 0; k < sz; k++)
        for (j = 0; j < sy; j++)
            for (i = 0; i < sx; i++)
            {
                if (in[k * (sx * sy) + (j * sx) + i] == 0.0)
                    continue;

                acu = 0;
                ind = 0;
                for (ii = -r; ii <= r; ii++)
                {
                    ni = i + ii;
                    if (ni < 0)
                        ni = -ni;
                    if (ni >= sx)
                        ni = 2 * sx - ni - 1;
                    if (in[k * (sx * sy) + (j * sx) + ni] > 0)
                    {
                        acu += (double)in[k * (sx * sy) + (j * sx) + ni];
                        ind++;
                    }
                }
                if (ind == 0)
                    ind = 1;
                out[k * (sx * sy) + (j * sx) + i] = (float)acu / ind;
            }

    for (k = 0; k < sz; k++)
        for (j = 0; j < sy; j++)
            for (i = 0; i < sx; i++)
            {
                if (out[k * (sx * sy) + (j * sx) + i] == 0)
                    continue;

                acu = 0;
                ind = 0;
                for (jj = -r; jj <= r; jj++)
                {
                    nj = j + jj;
                    if (nj < 0)
                        nj = -nj;
                    if (nj >= sy)
                        nj = 2 * sy - nj - 1;
                    if (out[k * (sx * sy) + (nj * sx) + i] > 0)
                    {
                        acu += (double)out[k * (sx * sy) + (nj * sx) + i];
                        ind++;
                    }
                }
                if (ind == 0)
                    ind = 1;
                temp[k * (sx * sy) + (j * sx) + i] = acu / ind;
            }

    for (k = 0; k < sz; k++)
        for (j = 0; j < sy; j++)
            for (i = 0; i < sx; i++)
            {
                if (temp[k * (sx * sy) + (j * sx) + i] == 0)
                    continue;

                acu = 0;
                ind = 0;
                for (kk = -r; kk <= r; kk++)
                {
                    nk = k + kk;
                    if (nk < 0)
                        nk = -nk;
                    if (nk >= sz)
                        nk = 2 * sz - nk - 1;
                    if (temp[nk * (sx * sy) + (j * sx) + i] > 0)
                    {
                        acu += temp[nk * (sx * sy) + (j * sx) + i];
                        ind++;
                    }
                }
                if (ind == 0)
                    ind = 1;
                out[k * (sx * sy) + (j * sx) + i] = (float)acu / ind;
            }

    free(temp);
    return;
}

#if defined(_WIN32)
unsigned int
#else
void *
#endif
/**
 * \brief Thread worker for SANLM filtering.
 *
 * Processes a subset of slices and accumulates block estimates using
 * adaptive non-local means weights.
 *
 * \param pArguments (in) thread argument struct
 * \return Thread exit code
 */
ThreadFunc(void *pArguments)
{
    float *bias, *Estimate, *ima, *means, *variances, *average;
    double epsilon, mu1, var1, totalweight, wmax, t1, t1i, t2, d, w, distanciaminima;
    double strength, scale;
    unsigned char *Label;
    int rows, cols, slices, ini, fin, v, f, i, j, k, l, rc, ii, jj, kk, ni, nj, nk, Ndims;

    extern int rician;
    extern double max;

    myargument arg;
    arg = *(myargument *)pArguments;

    rows = arg.rows;
    cols = arg.cols;
    slices = arg.slices;
    ini = arg.ini;
    fin = arg.fin;
    ima = arg.in_image;
    means = arg.means_image;
    variances = arg.var_image;
    Estimate = arg.estimate;
    bias = arg.bias;
    Label = arg.label;
    v = arg.radioB;
    f = arg.radioS;
    strength = arg.strength;

    /* filter */
    epsilon = 1e-5;
    mu1 = 0.95;
    var1 = 0.5;
    rc = rows * cols;

    Ndims = (2 * f + 1) * (2 * f + 1) * (2 * f + 1);

    average = (float *)malloc(Ndims * sizeof(float));

    wmax = 0.0;

    for (k = ini; k < fin; k += 2)
        for (j = 0; j < rows; j += 2)
            for (i = 0; i < cols; i += 2)
            {
                /* init */
                for (l = 0; l < Ndims; l++)
                    average[l] = 0.0;
                totalweight = 0.0;
                distanciaminima = 1e15;

                if (ima[(k * rc) + (j * cols) + i] > 0 && (means[(k * rc) + (j * cols) + i]) > epsilon && (variances[(k * rc) + (j * cols) + i] > epsilon))
                {
                    /* calculate minimum distance */
                    for (kk = -v; kk <= v; kk++)
                    {
                        nk = k + kk;
                        for (jj = -v; jj <= v; jj++)
                        {
                            nj = j + jj;
                            for (ii = -v; ii <= v; ii++)
                            {
                                ni = i + ii;
                                if (ii == 0 && jj == 0 && kk == 0)
                                    continue;

                                if (ni >= 0 && nj >= 0 && nk >= 0 && ni < cols && nj < rows && nk < slices)
                                {
                                    if (ima[(nk * rc) + (nj * cols) + ni] > 0 && (means[(nk * rc) + (nj * cols) + ni]) > epsilon && (variances[(nk * rc) + (nj * cols) + ni] > epsilon))
                                    {
                                        t1 = ((double)means[(k * rc) + (j * cols) + i]) / ((double)means[(nk * rc) + (nj * cols) + ni]);
                                        t1i = (max - (double)means[(k * rc) + (j * cols) + i]) / (max - (double)means[(nk * rc) + (nj * cols) + ni]);
                                        t2 = ((double)variances[(k * rc) + (j * cols) + i]) / ((double)variances[(nk * rc) + (nj * cols) + ni]);

                                        if ((t1 > mu1 && t1 < (1.0 / mu1)) || ((t1i > mu1 && t1i < (1.0 / mu1)) && t2 > var1 && t2 < (1.0 / var1)))
                                        {
                                            d = distance2(ima, means, i, j, k, ni, nj, nk, f, cols, rows, slices);
                                            if (d < distanciaminima)
                                                distanciaminima = d;
                                        }
                                    }
                                }
                            }
                        }
                    }
                    if (distanciaminima == 0)
                        distanciaminima = 1;

                    /* rician correction */
                    if (rician)
                    {
                        for (kk = -f; kk <= f; kk++)
                        {
                            nk = k + kk;
                            for (ii = -f; ii <= f; ii++)
                            {
                                ni = i + ii;
                                for (jj = -f; jj <= f; jj++)
                                {
                                    nj = j + jj;
                                    if (ni >= 0 && nj >= 0 && nk >= 0 && ni < cols && nj < rows && nk < slices)
                                    {
                                        if (distanciaminima == 1e15)
                                            bias[(nk * rc) + (nj * cols) + ni] = 0.0;
                                        else
                                            bias[(nk * rc) + (nj * cols) + ni] = (float)distanciaminima;
                                    }
                                }
                            }
                        }
                    }

                    /* block filtering */
                    for (kk = -v; kk <= v; kk++)
                    {
                        nk = k + kk;
                        for (jj = -v; jj <= v; jj++)
                        {
                            nj = j + jj;
                            for (ii = -v; ii <= v; ii++)
                            {
                                ni = i + ii;
                                if (ii == 0 && jj == 0 && kk == 0)
                                    continue;

                                if (ni >= 0 && nj >= 0 && nk >= 0 && ni < cols && nj < rows && nk < slices)
                                {
                                    if (ima[(nk * rc) + (nj * cols) + ni] > 0 && (means[(nk * rc) + (nj * cols) + ni]) > epsilon && (variances[(nk * rc) + (nj * cols) + ni] > epsilon))
                                    {
                                        t1 = ((double)means[(k * rc) + (j * cols) + i]) / ((double)means[(nk * rc) + (nj * cols) + ni]);
                                        t1i = (max - (double)means[(k * rc) + (j * cols) + i]) / (max - (double)means[(nk * rc) + (nj * cols) + ni]);
                                        t2 = ((double)variances[(k * rc) + (j * cols) + i]) / ((double)variances[(nk * rc) + (nj * cols) + ni]);

                                        if ((t1 > mu1 && t1 < (1.0 / mu1)) || ((t1i > mu1 && t1i < (1.0 / mu1)) && t2 > var1 && t2 < (1.0 / var1)))
                                        {
                                            d = distance(ima, i, j, k, ni, nj, nk, f, cols, rows, slices);

                                            scale = distanciaminima;
                                            if (strength > 0.0)
                                                scale *= strength;
                                            if (scale <= 0.0)
                                                scale = 1.0;

                                            if (d > 3.0 * scale)
                                                w = 0;
                                            else
                                                w = exp(-d / scale);

                                            if (w > wmax)
                                                wmax = w;

                                            if (w > 0)
                                            {
                                                Average_block(ima, ni, nj, nk, f, average, w, cols, rows, slices);
                                                totalweight = totalweight + w;
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }

                    if (wmax == 0.0)
                        wmax = 1.0;
                    Average_block(ima, i, j, k, f, average, wmax, cols, rows, slices);
                    totalweight = totalweight + wmax;
#if !defined(_WIN32)
                    pthread_mutex_lock(&mutex);
#endif
                    Value_block(Estimate, Label, i, j, k, f, average, totalweight, cols, rows, slices);
#if !defined(_WIN32)
                    pthread_mutex_unlock(&mutex);
#endif
                }
                else
                {
                    wmax = 1.0;
                    Average_block(ima, i, j, k, f, average, wmax, cols, rows, slices);
                    totalweight = totalweight + wmax;
#if !defined(_WIN32)
                    pthread_mutex_lock(&mutex);
#endif
                    Value_block(Estimate, Label, i, j, k, f, average, totalweight, cols, rows, slices);
#if !defined(_WIN32)
                    pthread_mutex_unlock(&mutex);
#endif
                }
            }

#if defined(_WIN32)
    return 0;
#else
    pthread_exit(0);
#endif

    free(average);
    return 0;
}

/**
 * \brief Spatially adaptive non-local means filter for 3D images.
 *
 * Applies SANLM denoising with optional Rician correction using
 * multithreading. The image is filtered in-place.
 *
 * \param ima         (in/out) input image volume
 * \param v           (in)  search window half-size
 * \param f           (in)  patch window half-size
 * \param use_rician  (in)  non-zero for Rician correction
 * \param strength    (in)  strength scaling for adaptive weights
 * \param dims        (in)  volume dimensions [x, y, z]
 * \return void
 */
void sanlm(float *ima, int v, int f, int use_rician, double strength, const int *dims)
{
    float *means, *variances, *Estimate, *bias;
    unsigned char *Label;
    int ndim = 3;
    double SNR, mean, var, estimate, d;
    int vol, slice, label, Ndims, i, j, k, ii, jj, kk, ni, nj, nk, indice, Nthreads, ini, fin, r;

    extern int rician;

    myargument *ThreadArgs;

#if !defined(_WIN32)
    pthread_t *ThreadList;
#endif

    if (strength <= 0.0)
        strength = 1.0;

    Ndims = (int)floor(pow((2.0 * f + 1.0), ndim));
    slice = dims[0] * dims[1];
    vol = dims[0] * dims[1] * dims[2];

    /* Allocate memory */
    means = (float *)malloc(vol * sizeof(float));
    variances = (float *)malloc(vol * sizeof(float));
    Estimate = (float *)malloc(vol * sizeof(float));
    Label = (unsigned char *)malloc(vol * sizeof(unsigned char));

    /* set global parameter */
    if (use_rician)
        rician = 1;

    if (rician)
        bias = (float *)malloc(vol * sizeof(float));

    for (i = 0; i < vol; i++)
    {
        Estimate[i] = 0.0;
        Label[i] = 0;
        if (rician)
            bias[i] = 0.0;
    }

    max = 0.0;
    for (k = 0; k < dims[2]; k++)
    {
        for (j = 0; j < dims[1]; j++)
        {
            for (i = 0; i < dims[0]; i++)
            {
                if (ima[k * (slice) + (j * dims[0]) + i] > max)
                    max = (double)ima[k * (slice) + (j * dims[0]) + i];

                mean = 0.0;
                indice = 0;
                for (ii = -1; ii <= 1; ii++)
                {
                    for (jj = -1; jj <= 1; jj++)
                    {
                        for (kk = -1; kk <= 1; kk++)
                        {
                            ni = i + ii;
                            nj = j + jj;
                            nk = k + kk;

                            if (ni < 0)
                                ni = -ni;
                            if (nj < 0)
                                nj = -nj;
                            if (nk < 0)
                                nk = -nk;
                            if (ni >= dims[0])
                                ni = 2 * dims[0] - ni - 1;
                            if (nj >= dims[1])
                                nj = 2 * dims[1] - nj - 1;
                            if (nk >= dims[2])
                                nk = 2 * dims[2] - nk - 1;

                            mean += (double)ima[nk * (slice) + (nj * dims[0]) + ni];
                            indice++;
                        }
                    }
                }
                mean /= (double)indice;
                means[k * (slice) + (j * dims[0]) + i] = (float)mean;
            }
        }
    }

    for (k = 0; k < dims[2]; k++)
    {
        for (j = 0; j < dims[1]; j++)
        {
            for (i = 0; i < dims[0]; i++)
            {
                var = 0;
                indice = 0;
                for (ii = -1; ii <= 1; ii++)
                {
                    for (jj = -1; jj <= 1; jj++)
                    {
                        for (kk = -1; kk <= 1; kk++)
                        {
                            ni = i + ii;
                            nj = j + jj;
                            nk = k + kk;
                            if (ni >= 0 && nj >= 0 && nk > 0 && ni < dims[0] && nj < dims[1] && nk < dims[2])
                            {
                                d = (double)ima[nk * (slice) + (nj * dims[0]) + ni] - (double)means[k * (slice) + (j * dims[0]) + i];
                                var += d * d;
                                indice++;
                            }
                        }
                    }
                }
                var /= (indice - 1);
                variances[k * (slice) + (j * dims[0]) + i] = (float)var;
            }
        }
    }

    Nthreads = dims[2] < MAX_NTHREADS ? dims[2] : MAX_NTHREADS;
    if (Nthreads < 1)
        Nthreads = 1;

#if defined(_WIN32)
    /* Sequential execution on Windows (no pthread dependency) */
    ThreadArgs = (myargument *)malloc(Nthreads * sizeof(myargument));

    for (i = 0; i < Nthreads; i++)
    {
        /* Make Thread Structure */
        ini = (i * dims[2]) / Nthreads;
        fin = ((i + 1) * dims[2]) / Nthreads;
        ThreadArgs[i].cols = dims[0];
        ThreadArgs[i].rows = dims[1];
        ThreadArgs[i].slices = dims[2];
        ThreadArgs[i].in_image = ima;
        ThreadArgs[i].var_image = variances;
        ThreadArgs[i].means_image = means;
        ThreadArgs[i].estimate = Estimate;
        ThreadArgs[i].bias = bias;
        ThreadArgs[i].label = Label;
        ThreadArgs[i].ini = ini;
        ThreadArgs[i].fin = fin;
        ThreadArgs[i].radioB = v;
        ThreadArgs[i].radioS = f;
        ThreadArgs[i].strength = strength;

        ThreadFunc(&ThreadArgs[i]);
    }
    free(ThreadArgs);

#else

    /* Reserve room for handles of threads in ThreadList*/
    ThreadList = (pthread_t *)calloc(Nthreads, sizeof(pthread_t));
    ThreadArgs = (myargument *)calloc(Nthreads, sizeof(myargument));

    for (i = 0; i < Nthreads; i++)
    {
        /* Make Thread Structure */
        ini = (i * dims[2]) / Nthreads;
        fin = ((i + 1) * dims[2]) / Nthreads;
        ThreadArgs[i].cols = dims[0];
        ThreadArgs[i].rows = dims[1];
        ThreadArgs[i].slices = dims[2];
        ThreadArgs[i].in_image = ima;
        ThreadArgs[i].var_image = variances;
        ThreadArgs[i].means_image = means;
        ThreadArgs[i].estimate = Estimate;
        ThreadArgs[i].bias = bias;
        ThreadArgs[i].label = Label;
        ThreadArgs[i].ini = ini;
        ThreadArgs[i].fin = fin;
        ThreadArgs[i].radioB = v;
        ThreadArgs[i].radioS = f;
        ThreadArgs[i].strength = strength;
    }

    for (i = 0; i < Nthreads; i++)
    {
        if (pthread_create(&ThreadList[i], NULL, ThreadFunc, &ThreadArgs[i]))
        {
            printf("Threads cannot be created\n");
            exit(1);
        }
    }

    for (i = 0; i < Nthreads; i++)
        pthread_join(ThreadList[i], NULL);

#endif

    if (rician)
    {
        r = 5;
        Regularize(bias, variances, r, dims[0], dims[1], dims[2]);
        for (i = 0; i < vol; i++)
        {
            if (variances[i] > 0.0)
            {
                SNR = (double)means[i] / sqrt((double)variances[i]);
                bias[i] = 2 * (variances[i] / (float)Epsi(SNR));
#if defined(_WIN32)
                if (_isnan(bias[i]))
                    bias[i] = 0.0;
#else
                if (isnan(bias[i]))
                    bias[i] = 0.0;
#endif
            }
        }
    }

    /* Aggregation of the estimators (i.e. means computation) */
    label = 0;
    estimate = 0.0;
    for (i = 0; i < vol; i++)
    {
        label = Label[i];
        if (label > 0)
        {
            estimate = (double)Estimate[i];
            estimate /= (double)label;
            if (rician)
            {
                estimate = (estimate - (double)bias[i]) < 0 ? 0 : (estimate - (double)bias[i]);
                ima[i] = (float)sqrt(estimate);
            }
            else
                ima[i] = (float)estimate;
        }
    }

    free(ThreadList);
    free(ThreadArgs);
    free(means);
    free(variances);
    free(Estimate);
    free(Label);
    if (rician)
        free(bias);

    return;
}
