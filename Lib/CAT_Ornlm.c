/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 *
 */

/*
 *
 * This code is a modified version of ornlm.c
 * from Pierrick Coupe and Jose V. Manon
 *
 * Original authors :
 * Pierrick Coupe - pierrick.coupe@gmail.com
 * Jose V. Manjon - jmanjon@fis.upv.es
 * Brain Imaging Center, Montreal Neurological Institute.
 * Mc Gill University
 *
 * Copyright (C) 2008 Pierrick Coupe and Jose V. Manjon
 *
 *
 *                          Details on ONLM filter
 ****************************************************************************
 *  The ONLM filter is described in:
 *
 *  P. Coupe, P. Yger, S. Prima, P. Hellier, C. Kervrann, C. Barillot.
 *  An Optimized Blockwise Non Local Means Denoising Filter for 3D Magnetic
 *  Resonance Images. IEEE Transactions on Medical Imaging, 27(4):425-441,
 *  April 2008
 ****************************************************************************
 *
 *
 *                      Details on Rician adaptation
 ****************************************************************************
 *  The adaptation to Rician noise is described in:
 *
 *  N. Wiest-Daessle, S. Prima, P. Coupe, S.P. Morrissey, C. Barillot.
 *  Rician noise removal by non-local means filtering for low
 *  signal-to-noise ratio MRI: Applications to DT-MRI. In 11th
 *  International Conference on Medical Image Computing and
 *  Computer-Assisted Intervention, MICCAI'2008,
 *  Pages 171-179, New York, USA, September 2008
 ****************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "CAT_Nlm.h"
#include "CAT_SafeAlloc.h"
#include "CAT_Math.h"

#if defined(_WIN32) || defined(_WIN64)
  #include <windows.h>
  #include <process.h>  /* _beginthreadex, _endthreadex */
#else
  #include <pthread.h>
  pthread_mutex_t mutex_ornlm = PTHREAD_MUTEX_INITIALIZER;
#endif

#define MAX_NTHREADS  16

/* Function which compute the weighted average for one block */
void Average_block_ornlm(float *ima, int x, int y, int z, int neighborhoodsize,
                         float *average, float weight, const int* vol_size)
{
    int x_pos, y_pos, z_pos;
    int is_outside;

    int a, b, c;

    int count = 0;

    for (c = 0; c<(2*neighborhoodsize+1); c++) {
        for (b = 0; b<(2*neighborhoodsize+1); b++) {
            for (a = 0; a<(2*neighborhoodsize+1); a++) {

                is_outside = 0;
                x_pos = x+a-neighborhoodsize;
                y_pos = y+b-neighborhoodsize;
                z_pos = z+c-neighborhoodsize;

                if ((z_pos < 0) || (z_pos > vol_size[2]-1)) is_outside = 1;
                if ((y_pos < 0) || (y_pos > vol_size[0]-1)) is_outside = 1;
                if ((x_pos < 0) || (x_pos > vol_size[1]-1)) is_outside = 1;

                if (is_outside)
                    average[count] += ima[z*(vol_size[0]*vol_size[1])+(x*vol_size[0])+y]
                                    * ima[z*(vol_size[0]*vol_size[1])+(x*vol_size[0])+y] * weight;
                else
                    average[count] += ima[z_pos*(vol_size[0]*vol_size[1])+(x_pos*vol_size[0])+y_pos]
                                    * ima[z_pos*(vol_size[0]*vol_size[1])+(x_pos*vol_size[0])+y_pos] * weight;

                count++;
            }
        }
    }
}

/* Function which computes the value assigned to each voxel */
void Value_block_ornlm(float *Estimate, unsigned char *Label, int x, int y, int z,
                       int neighborhoodsize, float *average, float global_sum,
                       const int* vol_size, float hh)
{
    int x_pos, y_pos, z_pos;
    int is_outside;
    float value = 0.0f;
    float denoised_value = 0.0f;
    unsigned char label = 0;
    int count=0 ;
    int a, b, c;

    for (c = 0; c<(2*neighborhoodsize+1); c++) {
        for (b = 0; b<(2*neighborhoodsize+1); b++) {
            for (a = 0; a<(2*neighborhoodsize+1); a++) {

                is_outside = 0;
                x_pos = x+a-neighborhoodsize;
                y_pos = y+b-neighborhoodsize;
                z_pos = z+c-neighborhoodsize;

                if ((z_pos < 0) || (z_pos > vol_size[2]-1)) is_outside = 1;
                if ((y_pos < 0) || (y_pos > vol_size[0]-1)) is_outside = 1;
                if ((x_pos < 0) || (x_pos > vol_size[1]-1)) is_outside = 1;

                if (!is_outside) {
                    value = Estimate[z_pos*(vol_size[0]*vol_size[1])+(x_pos*vol_size[0])+y_pos];
                    denoised_value  = (average[count]/global_sum) - hh;
                    if (denoised_value > 0)
                        denoised_value = sqrt(denoised_value);
                    else denoised_value = 0.0f;

                    value += denoised_value;
                    label = Label[(y_pos + x_pos*vol_size[0] + z_pos *vol_size[0] * vol_size[1])];
                    Estimate[z_pos*(vol_size[0]*vol_size[1])+(x_pos*vol_size[0])+y_pos] = value;
                    Label[(y_pos + x_pos*vol_size[0] + z_pos *vol_size[0] * vol_size[1])] = label + 1;
                }
                count++;
            }
        }
    }
}

float distance_ornlm(float* ima, int x, int y, int z,
                     int nx, int ny, int nz, int f,
                     int sx, int sy, int sz)
{
    float d, acu, distancetotal;
    int i, j, k, ni1, nj1, ni2, nj2, nk1, nk2;

    acu=0;
    distancetotal=0;

    for(k=-f; k<=f; k++) {
        for(i=-f; i<=f; i++) {
            for(j=-f; j<=f; j++) {
                ni1=x+i; nj1=y+j; nk1=z+k;
                ni2=nx+i; nj2=ny+j; nk2=nz+k;

                if(ni1<0) ni1=-ni1;
                if(nj1<0) nj1=-nj1;
                if(ni2<0) ni2=-ni2;
                if(nj2<0) nj2=-nj2;
                if(nk1<0) nk1=-nk1;
                if(nk2<0) nk2=-nk2;

                if(ni1>=sx) ni1=2*sx-ni1-1;
                if(nj1>=sy) nj1=2*sy-nj1-1;
                if(nk1>=sz) nk1=2*sz-nk1-1;
                if(ni2>=sx) ni2=2*sx-ni2-1;
                if(nj2>=sy) nj2=2*sy-nj2-1;
                if(nk2>=sz) nk2=2*sz-nk2-1;

                float diff = ima[nk1*(sx*sy)+(ni1*sy)+nj1] - ima[nk2*(sx*sy)+(ni2*sy)+nj2];
                distancetotal += diff*diff;
                acu += 1;
            }
        }
    }

    d = distancetotal/acu;
    return d;
}

/* ----------------------- Multithreading additions ----------------------- */

typedef struct{
    int rows;
    int cols;
    int slices;
    float* in_image;
    float* means_image;
    float* var_image;
    float* estimate;
    unsigned char* label;
    int ini;
    int fin;
    int radioB;        /* v */
    int radioS;        /* f */
    float hh;          /* 2*h*h */
    float inv_h2;      /* 1/(h*h) */
    const int* dims;   /* dims pointer for Average_block_ornlm */
} myargument;

#if defined(_WIN32) || defined(_WIN64)
unsigned int __stdcall
#else
void *
#endif
ThreadFunc_ornlm( void* pArguments )
{
    myargument arg = *(myargument*)pArguments;

    const int cols   = arg.cols;
    const int rows   = arg.rows;
    const int slices = arg.slices;
    float* ima       = arg.in_image;
    float* means     = arg.means_image;
    float* variances = arg.var_image;
    float* Estimate  = arg.estimate;
    unsigned char* Label = arg.label;
    const int ini    = arg.ini;
    const int fin    = arg.fin;
    const int v      = arg.radioB;
    const int f      = arg.radioS;
    const float hh   = arg.hh;
    const float inv_h2 = arg.inv_h2;
    const int* dims  = arg.dims;

    const float epsilon = 1e-5f;
    const float mu1 = 0.9f;
    const float var1 = 0.5f;

    const int Ndims = (2*f+1)*(2*f+1)*(2*f+1);
    float* average = SAFE_MALLOC(float, Ndims);
    if (!average) {
    #if defined(_WIN32) || defined(_WIN64)
        _endthreadex(0);
        return 0;
    #else
        pthread_exit(0);
        return 0;
    #endif
    }

    int k,i,j,ii,jj,kk,ni,nj,nk,init;
    float totalweight, wmax, d, t1, t2;

    for (k = ini; k < fin; k += 2) {
        for (i = 0; i < rows; i += 2) {
            for (j = 0; j < cols; j += 2) {

                for (init = 0; init < Ndims; ++init) average[init] = 0.0f;
                totalweight = 0.0f;

                if ((means[k*(cols*rows)+(i*cols)+j] > epsilon) &&
                    (variances[k*(cols*rows)+(i*cols)+j] > epsilon)) {

                    wmax = 0.0f;

                    for (kk = -v; kk <= v; kk++) {
                        for (ii = -v; ii <= v; ii++) {
                            for (jj = -v; jj <= v; jj++) {

                                if (ii==0 && jj==0 && kk==0) continue;

                                ni = i + ii;
                                nj = j + jj;
                                nk = k + kk;

                                if (ni>=0 && nj>=0 && nk>=0 &&
                                    ni<rows && nj<cols && nk<slices) {

                                    if ((means[nk*(cols*rows)+(ni*cols)+nj] > epsilon) &&
                                        (variances[nk*(cols*rows)+(ni*cols)+nj] > epsilon)) {

                                        t1 = (means[k*(cols*rows)+(i*cols)+j]) /
                                             (means[nk*(cols*rows)+(ni*cols)+nj]);
                                        t2 = (variances[k*(cols*rows)+(i*cols)+j]) /
                                             (variances[nk*(cols*rows)+(ni*cols)+nj]);

                                        if (t1>mu1 && t1<(1.0f/mu1) && t2>var1 && t2<(1.0f/var1)) {

                                            d = distance_ornlm(ima, i, j, k, ni, nj, nk, f,
                                                               rows, cols, slices);

                                            /* original weight: exp(-d/(h*h)) */
                                            float w = expf(-d * inv_h2);

                                            if (w > wmax) wmax = w;

                                            Average_block_ornlm(ima, ni, nj, nk, f, average, w, dims);
                                            totalweight += w;
                                        }
                                    }
                                }
                            }
                        }
                    }

                    if (wmax == 0.0f) wmax = 1.0f;

                    Average_block_ornlm(ima, i, j, k, f, average, wmax, dims);
                    totalweight += wmax;

                    if (totalweight != 0.0f) {
                    #if !defined(_WIN32) && !defined(_WIN64)
                        pthread_mutex_lock(&mutex_ornlm);
                    #endif
                        Value_block_ornlm(Estimate, Label, i, j, k, f, average, totalweight, dims, hh);
                    #if !defined(_WIN32) && !defined(_WIN64)
                        pthread_mutex_unlock(&mutex_ornlm);
                    #endif
                    }

                } else {
                    wmax = 1.0f;
                    Average_block_ornlm(ima, i, j, k, f, average, wmax, dims);
                    totalweight += wmax;

                #if !defined(_WIN32) && !defined(_WIN64)
                    pthread_mutex_lock(&mutex_ornlm);
                #endif
                    Value_block_ornlm(Estimate, Label, i, j, k, f, average, totalweight, dims, hh);
                #if !defined(_WIN32) && !defined(_WIN64)
                    pthread_mutex_unlock(&mutex_ornlm);
                #endif
                }
            }
        }
    }

    free(average);

#if defined(_WIN32) || defined(_WIN64)
    _endthreadex(0);
    return 0;
#else
    pthread_exit(0);
    return 0;
#endif
}

/* ----------------------- end multithreading additions ------------------- */

void ornlm(float* ima, int v, int f, float h, float sigma, const int* dims)
{
    float *means, *variances, *Estimate, *ima_out;
    unsigned char *Label;
    float mean, var, hh;
    double max_val;

    int vol;
    int i, j, k, ii, jj, kk, ni, nj, nk, indice;

    const float epsilon = 0.00001f;
    
    vol = dims[0]*dims[1]*dims[2];

    /* normalize image to 0..1 to make parameter invariant to scaling */
    max_val = get_max(ima, vol, 0, DT_FLOAT32);
    if (max_val > 0.0) {
        for (i = 0; i < vol; i++) ima[i] /= (float)max_val;
    }

    hh = 2.0f*sigma*sigma;

    ima_out  = SAFE_MALLOC(float, vol);
    means    = SAFE_MALLOC(float, vol);
    variances= SAFE_MALLOC(float, vol);
    Estimate = SAFE_MALLOC(float, vol);
    Label    = SAFE_MALLOC(unsigned char, vol);

    for (i = 0; i < vol; i++) {
        Estimate[i] = 0.0f;
        Label[i]    = 0;
        ima_out[i]  = 0.0f;
    }

    /* local means */
    for(k=0; k<dims[2]; k++) {
        for(i=0; i<dims[1]; i++) {
            for(j=0; j<dims[0]; j++) {
                mean=0.0f;
                indice=0;

                for(ii=-1; ii<=1; ii++) {
                    for(jj=-1; jj<=1; jj++) {
                        for(kk=-1; kk<=1; kk++) {
                            ni=i+ii; nj=j+jj; nk=k+kk;

                            if(ni<0) ni=-ni;
                            if(nj<0) nj=-nj;
                            if(nk<0) nk=-nk;
                            if(ni>=dims[1]) ni=2*dims[1]-ni-1;
                            if(nj>=dims[0]) nj=2*dims[0]-nj-1;
                            if(nk>=dims[2]) nk=2*dims[2]-nk-1;

                            mean += ima[nk*(dims[0]*dims[1])+(ni*dims[0])+nj];
                            indice += 1;
                        }
                    }
                }
                means[k*(dims[0]*dims[1])+(i*dims[0])+j] = mean/(float)indice;
            }
        }
    }

    /* local variances */
    for(k=0; k<dims[2]; k++) {
        for(i=0; i<dims[1]; i++) {
            for(j=0; j<dims[0]; j++) {
                var=0.0f;
                indice=0;

                for(ii=-1; ii<=1; ii++) {
                    for(jj=-1; jj<=1; jj++) {
                        for(kk=-1; kk<=1; kk++) {
                            ni=i+ii; nj=j+jj; nk=k+kk;
                            if(ni>=0 && nj>=0 && nk>0 && ni<dims[1] && nj<dims[0] && nk<dims[2]) {
                                float diff = ima[nk*(dims[0]*dims[1])+(ni*dims[0])+nj] -
                                             means[k*(dims[0]*dims[1])+(i*dims[0])+j];
                                var += diff*diff;
                                indice += 1;
                            }
                        }
                    }
                }
                var /= (float)(indice-1);
                variances[k*(dims[0]*dims[1])+(i*dims[0])+j]=var;
            }
        }
    }

    /* --------------------- multithreaded filter ------------------------ */
    {
        int Nthreads = dims[2] < MAX_NTHREADS ? dims[2] : MAX_NTHREADS;
        if (Nthreads < 1) Nthreads = 1;

        myargument *ThreadArgs;

    #if defined(_WIN32) || defined(_WIN64)
        HANDLE *ThreadList; /* Handles to the worker threads*/
    ThreadList = SAFE_MALLOC(HANDLE, Nthreads);
    ThreadArgs = SAFE_MALLOC(myargument, Nthreads);

        for (i=0; i<Nthreads; i++) {
            int ini = (i*dims[2])/Nthreads;
            int fin = ((i+1)*dims[2])/Nthreads;

            ThreadArgs[i].cols    = dims[0];
            ThreadArgs[i].rows    = dims[1];
            ThreadArgs[i].slices  = dims[2];
            ThreadArgs[i].in_image   = ima;
            ThreadArgs[i].means_image = means;
            ThreadArgs[i].var_image   = variances;
            ThreadArgs[i].estimate    = Estimate;
            ThreadArgs[i].label       = Label;
            ThreadArgs[i].ini = ini;
            ThreadArgs[i].fin = fin;
            ThreadArgs[i].radioB = v;
            ThreadArgs[i].radioS = f;
            ThreadArgs[i].hh     = hh;
            ThreadArgs[i].inv_h2 = 1.0f/(h*h);
            ThreadArgs[i].dims   = dims;

            ThreadList[i] = (HANDLE)_beginthreadex(NULL, 0, &ThreadFunc_ornlm, &ThreadArgs[i], 0, NULL);
        }

        for (i=0; i<Nthreads; i++) { WaitForSingleObject(ThreadList[i], INFINITE); }
        for (i=0; i<Nthreads; i++) { CloseHandle(ThreadList[i]); }

    #else /* POSIX */
        pthread_t *ThreadList;
    ThreadList = SAFE_CALLOC(pthread_t, Nthreads);
    ThreadArgs = SAFE_CALLOC(myargument, Nthreads);

        for (i=0; i<Nthreads; i++) {
            int ini = (i*dims[2])/Nthreads;
            int fin = ((i+1)*dims[2])/Nthreads;

            ThreadArgs[i].cols    = dims[0];
            ThreadArgs[i].rows    = dims[1];
            ThreadArgs[i].slices  = dims[2];
            ThreadArgs[i].in_image   = ima;
            ThreadArgs[i].means_image = means;
            ThreadArgs[i].var_image   = variances;
            ThreadArgs[i].estimate    = Estimate;
            ThreadArgs[i].label       = Label;
            ThreadArgs[i].ini = ini;
            ThreadArgs[i].fin = fin;
            ThreadArgs[i].radioB = v;
            ThreadArgs[i].radioS = f;
            ThreadArgs[i].hh     = hh;
            ThreadArgs[i].inv_h2 = 1.0f/(h*h);
            ThreadArgs[i].dims   = dims;

            if (pthread_create(&ThreadList[i], NULL, ThreadFunc_ornlm, &ThreadArgs[i])) {
                printf("Threads cannot be created\n");
                exit(1);
            }
        }

        for (i=0; i<Nthreads; i++) pthread_join(ThreadList[i], NULL);
    #endif

        free(ThreadList);
        free(ThreadArgs);
    }
    /* ------------------- end multithreaded filter ---------------------- */

    /* Aggregation of the estimators (i.e. means computation) */
    for (k = 0; k < dims[2]; k++ ) {
        for (i = 0; i < dims[1]; i++ ) {
            for (j = 0; j < dims[0]; j++ ) {
                unsigned char label = Label[k*(dims[0]*dims[1])+(i*dims[0])+j];

                if (label == 0) {
                    /* FIX: i*dims[0] (vorher i*dims[1]) */
                    ima_out[k*(dims[0]*dims[1])+(i*dims[0])+j] =
                        ima[k*(dims[0]*dims[1])+(i*dims[0])+j];
                } else {
                    float estimate = Estimate[k*(dims[0]*dims[1])+(i*dims[0])+j];
                    estimate /= (float)label;
                    ima_out[k*(dims[0]*dims[1])+(i*dims[0])+j] = estimate;
                }
            }
        }
    }

    /* Return filtered image to input */
    if (max_val > 0.0) {
        for (i = 0; i < vol; i++) ima[i] = ima_out[i]*(float)max_val;
    } else {
        for (i = 0; i < vol; i++) ima[i] = ima_out[i];
    }

    free(means);
    free(variances);
    free(Estimate);
    free(Label);
    free(ima_out);
}
