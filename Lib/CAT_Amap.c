/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "CAT_Amap.h"
#include "CAT_Vol.h"
#include "CAT_Math.h"

/* calculate the mean and variance for every class on a grid size SUBxSUBxSUB */
#if defined(_WIN32) || defined(_WIN64)
  #include <windows.h>
  #include <process.h>
#else
  #include <pthread.h>
#endif

#define MAX_NTHREADS 16

#if defined(_WIN32) || defined(_WIN64)
static unsigned int __stdcall
#else
static void *
#endif
gmv_accum_worker(void *p)
{
    gmv_accum_args_t a = *(gmv_accum_args_t*)p;
    int k, l, m, x, y, z;

    const int sub = a.sub;
    const int area = a.area;
    const int narea = a.narea;
    const int nvol  = a.nvol;
    const int *dims = a.dims;

    for (k = -sub; k <= sub; ++k)
    for (l = -sub; l <= sub; ++l)
    for (m = -sub; m <= sub; ++m) {

        for (z = a.z_ini; z < a.z_fin; ++z) {
            const int zsub = z*sub + k;
            if (zsub < 0 || zsub >= dims[2]) continue;

            const int zsub2   = zsub * area;
            const int zoffset = z * narea;

            for (y = 0; y < a.niy; ++y) {
                const int ysub = y*sub + l;
                if (ysub < 0 || ysub >= dims[1]) continue;

                const int ysub2  = ysub * dims[0];
                const int yoffset= zoffset + y*a.nix;

                for (x = 0; x < a.nix; ++x) {
                    const int xsub = x*sub + m;
                    if (xsub < 0 || xsub >= dims[0]) continue;

                    const int vox_idx = zsub2 + ysub2 + xsub;
                    int label_value = (int)a.label[vox_idx];
                    int label_value_BG = label_value - 1;
                    if (label_value_BG < 0) continue;

                    double val = (double)a.src[vox_idx];

                    if (val < a.thresh[0]) continue;
                    /* if (val > a.thresh[1]) continue; */

                    const int grid_idx = yoffset + x;          /* j */
                    const int ind = label_value_BG * nvol + grid_idx;

                    struct ipoint *pi = &a.ir[ind];

                    pi->arr[pi->n] = val;
                    pi->n++;
                    pi->s  += val;
                    pi->ss += val*val;
                }
            }
        }
    }

#if defined(_WIN32) || defined(_WIN64)
    _endthreadex(0);
    return 0;
#else
    return 0;
#endif
}

#if defined(_WIN32) || defined(_WIN64)
static unsigned int __stdcall
#else
static void *
#endif
gmv_reduce_worker(void *p)
{
    gmv_reduce_args_t a = *(gmv_reduce_args_t*)p;
    int i, j;

    for (i = 0; i < a.n_classes; ++i) {
        const int base = i * a.nvol;
        for (j = a.j_ini; j < a.j_fin; ++j) {
            const int ind = base + j;
            const struct ipoint *pi = &a.ir[ind];

            if (pi->n > G) {
                if (pi->n == 1) {
                    a.r[ind].var  = 0.0;
                    a.r[ind].mean = pi->arr[0];
                } else {
                    if (a.use_median)
                        a.r[ind].mean = get_median_double(pi->arr, pi->n, 0);
                    else
                        a.r[ind].mean = pi->s / pi->n;
                    a.r[ind].var = (pi->ss - pi->n * SQR(pi->s / pi->n)) / (pi->n - 1);
                }
            } else {
                a.r[ind].mean = 0.0;
            }
        }
    }

#if defined(_WIN32) || defined(_WIN64)
    _endthreadex(0);
    return 0;
#else
    return 0;
#endif
}

/* calculate the mean and variance for every class on a grid size SUBxSUBxSUB */
static void GetMeansVariances(float *src, unsigned char *label, int n_classes,
                              struct point *r, int sub, int *dims, double *thresh, int use_median)
{
    int area, narea, nvol, nix, niy, niz;

    area = dims[0]*dims[1];
    int i, j, t;

    /* grid-size */
    nix = (int)ceil((dims[0]-1)/((double)sub)) + 1;
    niy = (int)ceil((dims[1]-1)/((double)sub)) + 1;
    niz = (int)ceil((dims[2]-1)/((double)sub)) + 1;
    narea = nix*niy;
    nvol  = nix*niy*niz;

    struct ipoint *ir = (struct ipoint*)malloc((size_t)n_classes * (size_t)nvol * sizeof(struct ipoint));
    if (!ir) { printf("Memory allocation error\n"); exit(EXIT_FAILURE); }

    /* maximal size of neighbourhood per grid point */
    const int sz_cube = (2*sub + 1)*(2*sub + 1)*(2*sub + 1);

    /* init accumulators + arr-buffer */
    for (i = 0; i < n_classes; ++i) {
        for (j = 0; j < nvol; ++j) {
            const int ind = i*nvol + j;
            ir[ind].n  = 0;
            ir[ind].s  = 0.0;
            ir[ind].ss = 0.0;
            ir[ind].arr = (double*)malloc((size_t)sz_cube * sizeof(double));
            if (!ir[ind].arr) { printf("Memory allocation error\n"); exit(EXIT_FAILURE); }
        }
    }

    /* ---------------- Stage 1: parallel accumulation ---------------- */
    {
        int Nthreads = (niz < MAX_NTHREADS) ? niz : MAX_NTHREADS;
        if (Nthreads < 1) Nthreads = 1;
        

    #if defined(_WIN32) || defined(_WIN64)
        HANDLE *ThreadList = (HANDLE*)malloc((size_t)Nthreads * sizeof(HANDLE));
        gmv_accum_args_t *Args = (gmv_accum_args_t*)malloc((size_t)Nthreads * sizeof(gmv_accum_args_t));
        if (!ThreadList || !Args) { printf("Memory allocation error\n"); exit(EXIT_FAILURE); }

        for (t = 0; t < Nthreads; ++t) {
            int ini = (t * niz) / Nthreads;
            int fin = ((t + 1) * niz) / Nthreads;

            Args[t].src = src;
            Args[t].label = label;
            Args[t].n_classes = n_classes;
            Args[t].sub = sub;
            Args[t].dims = dims;
            Args[t].thresh = thresh;

            Args[t].nix = nix; Args[t].niy = niy; Args[t].niz = niz;
            Args[t].narea = narea; Args[t].nvol = nvol; Args[t].area = area;

            Args[t].ir = ir;

            Args[t].z_ini = ini; Args[t].z_fin = fin;

            ThreadList[t] = (HANDLE)_beginthreadex(NULL, 0, &gmv_accum_worker, &Args[t], 0, NULL);
        }
        for (t = 0; t < Nthreads; ++t) WaitForSingleObject(ThreadList[t], INFINITE);
        for (t = 0; t < Nthreads; ++t) CloseHandle(ThreadList[t]);
        free(ThreadList); free(Args);
    #else
        pthread_t *ThreadList = (pthread_t*)calloc((size_t)Nthreads, sizeof(pthread_t));
        gmv_accum_args_t *Args = (gmv_accum_args_t*)calloc((size_t)Nthreads, sizeof(gmv_accum_args_t));
        if (!ThreadList || !Args) { printf("Memory allocation error\n"); exit(EXIT_FAILURE); }

        for (t = 0; t < Nthreads; ++t) {
            int ini = (t * niz) / Nthreads;
            int fin = ((t + 1) * niz) / Nthreads;

            Args[t].src = src;
            Args[t].label = label;
            Args[t].n_classes = n_classes;
            Args[t].sub = sub;
            Args[t].dims = dims;
            Args[t].thresh = thresh;

            Args[t].nix = nix; Args[t].niy = niy; Args[t].niz = niz;
            Args[t].narea = narea; Args[t].nvol = nvol; Args[t].area = area;

            Args[t].ir = ir;

            Args[t].z_ini = ini; Args[t].z_fin = fin;

            pthread_create(&ThreadList[t], NULL, gmv_accum_worker, &Args[t]);
        }
        for (t = 0; t < Nthreads; ++t) pthread_join(ThreadList[t], NULL);
        free(ThreadList); free(Args);
    #endif
    }

    /* ---------------- Stage 2: parallel reduction ------------------- */
    {
        int Nthreads = (nvol < 8) ? nvol : 8;
        if (Nthreads < 1) Nthreads = 1;

    #if defined(_WIN32) || defined(_WIN64)
        HANDLE *ThreadList = (HANDLE*)malloc((size_t)Nthreads * sizeof(HANDLE));
        gmv_reduce_args_t *Args = (gmv_reduce_args_t*)malloc((size_t)Nthreads * sizeof(gmv_reduce_args_t));
        if (!ThreadList || !Args) { printf("Memory allocation error\n"); exit(EXIT_FAILURE); }

        for (t = 0; t < Nthreads; ++t) {
            int j_ini = (t * nvol) / Nthreads;
            int j_fin = ((t + 1) * nvol) / Nthreads;

            Args[t].r = r;
            Args[t].ir = ir;
            Args[t].n_classes = n_classes;
            Args[t].nvol = nvol;
            Args[t].use_median = use_median;
            Args[t].j_ini = j_ini; Args[t].j_fin = j_fin;

            ThreadList[t] = (HANDLE)_beginthreadex(NULL, 0, &gmv_reduce_worker, &Args[t], 0, NULL);
        }
        for (t = 0; t < Nthreads; ++t) WaitForSingleObject(ThreadList[t], INFINITE);
        for (t = 0; t < Nthreads; ++t) CloseHandle(ThreadList[t]);
        free(ThreadList); free(Args);
    #else
        pthread_t *ThreadList = (pthread_t*)calloc((size_t)Nthreads, sizeof(pthread_t));
        gmv_reduce_args_t *Args = (gmv_reduce_args_t*)calloc((size_t)Nthreads, sizeof(gmv_reduce_args_t));
        if (!ThreadList || !Args) { printf("Memory allocation error\n"); exit(EXIT_FAILURE); }

        for (t = 0; t < Nthreads; ++t) {
            int j_ini = (t * nvol) / Nthreads;
            int j_fin = ((t + 1) * nvol) / Nthreads;

            Args[t].r = r;
            Args[t].ir = ir;
            Args[t].n_classes = n_classes;
            Args[t].nvol = nvol;
            Args[t].use_median = use_median;
            Args[t].j_ini = j_ini; Args[t].j_fin = j_fin;

            pthread_create(&ThreadList[t], NULL, gmv_reduce_worker, &Args[t]);
        }
        for (t = 0; t < Nthreads; ++t) pthread_join(ThreadList[t], NULL);
        free(ThreadList); free(Args);
    #endif
    }

    /* correct freeing of arrays */
    for (i = 0; i < n_classes; ++i) {
        for (j = 0; j < nvol; ++j) {
            int ind = i*nvol + j;
            free(ir[ind].arr);
        }
    }
    free(ir);
}

/* This PVE calculation is a modified version from
 * the PVE software bundle:
 * Copyright (C) Jussi Tohka, Institute of Signal Processing, Tampere University of
 * Technology, 2002 - 2004.
 * P.O. Box 553, FIN-33101, Finland
 * E-mail: jussi.tohka@tut.fi
 */
void Pve5(float *src, unsigned char *prob, unsigned char *label, double *mean, int *dims)
{
    int i, mxi;
    double w, mx;
    unsigned char new_val[MAX_NC];
    
    int vol = dims[0]*dims[1]*dims[2];
        
    for (i = 0; i < vol; i++) {

        switch(label[i]) {
        case 0: /* BG */
            new_val[CSFLABEL-1] = 0;
            new_val[GMLABEL-1]  = 0;
            new_val[WMLABEL-1]  = 0;
            break;
        case CSFLABEL: /* CSF */
            new_val[CSFLABEL-1] = 255;
            new_val[GMLABEL-1]  = 0;
            new_val[WMLABEL-1]  = 0;
            label[i] = (unsigned char) ROUND(255.0/3.0);
            break;
        case GMLABEL: /* GM */
            new_val[CSFLABEL-1] = 0;
            new_val[GMLABEL-1]  = 255;
            new_val[WMLABEL-1]  = 0;
            label[i] = (unsigned char) ROUND(2.0*255.0/3.0);
            break;
        case WMLABEL: /* WM */
            new_val[CSFLABEL-1] = 0;
            new_val[GMLABEL-1]  = 0;
            new_val[WMLABEL-1]  = 255;
            label[i] = 255;
            break;
        case GMCSFLABEL: /* GMCSF */
            w = ((double)src[i] - mean[CSFLABEL-1])/(mean[GMLABEL-1]-mean[CSFLABEL-1]);
            if (w > 1.0) w = 1.0; if (w < 0.0) w = 0.0;
            new_val[CSFLABEL-1] = (unsigned char) ROUND(255.0*(1-w));
            new_val[GMLABEL-1]  = (unsigned char) ROUND(255.0*w);
            new_val[WMLABEL-1]  = 0;
            label[i] = ROUND(255.0/3.0*(1.0 + w));
            break;
        case WMGMLABEL: /* WMGM */
            w = ((double)src[i] - mean[GMLABEL-1])/(mean[WMLABEL-1]-mean[GMLABEL-1]);
            if (w > 1.0) w = 1.0; if (w < 0.0) w = 0.0;
            new_val[CSFLABEL-1] = 0;
            new_val[GMLABEL-1]  = (unsigned char) ROUND(255.0*(1-w));
            new_val[WMLABEL-1]  = (unsigned char) ROUND(255.0*w);
            label[i] = ROUND(255.0/3.0*(2.0 + w));
            break;
        }

        prob[          i] = new_val[CSFLABEL-1];
        prob[vol     + i] = new_val[GMLABEL-1];
        prob[(2*vol) + i] = new_val[WMLABEL-1];
        
        /* set old probabilities for mixed classes to zero */
        prob[(3*vol) + i] = 0;
        prob[(4*vol) + i] = 0;
    }    
}


/* This code is a substantially modified version of MrfPrior.C 
 * from Jagath C. Rajapakse
 * 
 * Original author : Jagath C. Rajapakse
 *
 * See:
 * Statistical approach to single-channel MR brain scans
 * J. C. Rajapakse, J. N. Giedd, and J. L. Rapoport
 * IEEE Transactions on Medical Imaging, Vol 16, No 2, 1997
 *
 * Comments to raja@cns.mpg.de, 15.10.96
 */

void MrfPrior(unsigned char *label, int n_classes, double *alpha, double *beta, int init, int *dims, int verbose)
{
    int i, j, k, x, y, z;
    int fi, fj;
    long color[MAX_NC][7][7][7][7];
    long area;

    int f[MAX_NC-1], plab, zero, iBG;
    int n, z_area, y_dims;
    double XX, YY, L;
    
    area = dims[0]*dims[1];

    /* initialize configuration counts */
    if (beta != NULL) {
        for (i = 0; i < n_classes; i++)
            for (f[0] = 0; f[0] < 7; f[0]++)
                for (f[1] = 0; f[1] < 7; f[1]++)
                    for (f[2] = 0; f[2] < 7; f[2]++)
                        for (f[3] = 0; f[3] < 7; f[3]++)
                            color[i][f[0]][f[1]][f[2]][f[3]]=0;
    }

    /* calculate configuration counts */
    n = 0;
    for (i = 0; i < n_classes; i++) alpha[i] = 0.0;

    for (z = 1; z < dims[2]-1; z++) {
        z_area=z*area; 
        for (y = 1; y < dims[1]-1; y++) {
            y_dims = y*dims[0];
            for (x = 1; x < dims[0]-1; x++) {
            
                plab = (int)label[z_area + y_dims + x];
                
                zero = plab;
                if (zero < 1) continue;
                n++;
                alpha[zero - 1] += 1.0;
                
                if (beta != NULL) {
                    for (i = 1; i < n_classes; i++) {
                        f[i-1] = 0;       
                        iBG = i+1;
                        if ((int)label[z_area + y_dims + x-1] == iBG)        f[i-1]++;
                        if ((int)label[z_area + y_dims + x+1] == iBG)        f[i-1]++;
                        if ((int)label[z_area + ((y-1)*dims[0]) + x] == iBG) f[i-1]++;
                        if ((int)label[z_area + ((y+1)*dims[0]) + x] == iBG) f[i-1]++;
                        if ((int)label[((z-1)*area) + y_dims + x] == iBG)    f[i-1]++;
                        if ((int)label[((z+1)*area) + y_dims + x] == iBG)    f[i-1]++;
                    }
                    color[zero-1][f[0]][f[1]][f[2]][f[3]]++;
                }
            }
        }
    }

    /* evaluate alphas */
    if (verbose) printf("MRF priors: alpha ");
    for (i = 0; i < n_classes; i++) {
        if (init == 0) alpha[i] /= n; else alpha[i] = 1.0;
        if (verbose) printf("%3.3f ", alpha[i]);
    }

    /* compute beta */
    if (beta != NULL) {
        n = 0;
        XX=0.0, YY=0.0;
        for (f[0] = 0; f[0] < 7; f[0]++)
            for (f[1] = 0; f[1] < 7; f[1]++)
                for (f[2] = 0; f[2] < 7; f[2]++)
                    for (f[3] = 0; f[3] < 7; f[3]++)
                        for (i = 0; i < n_classes; i++)
                            for (j = 0; j < i; j++) {
                                n++;
                                if (color[i][f[0]][f[1]][f[2]][f[3]] < TH_COLOR ||
                                        color[j][f[0]][f[1]][f[2]][f[3]] < TH_COLOR) continue;
                     
                                L = log(((double) color[i][f[0]][f[1]][f[2]][f[3]])/
                                        (double) color[j][f[0]][f[1]][f[2]][f[3]]);
                     
                                if (i == 0) 
                                    fi = 6 - f[0] - f[1] - f[2] - f[3];
                                else fi = f[i-1];
                                                         
                                if (j == 0) 
                                    fj = 6 - f[0] - f[1] - f[2] - f[3];
                                else fj = f[j-1];
    
                                XX += (fi-fj)*(fi-fj);
                                YY += L*(fi-fj);
        }
        /* weighting of beta was empirically estimated using brainweb data with different noise levels
           because old beta estimation was not working */
        beta[0] = XX/YY;
        if (verbose) {
            printf("\t beta %3.3f\n", beta[0]);
            fflush(stdout);
        }
    } else if (verbose) {
        printf("\n");
        fflush(stdout);
    }
}

/* Computes likelihood of value given parameters mean and variance */ 
double ComputeGaussianLikelihood(double value, double mean, double var)

{ 
    return(exp(-(SQR(value - mean))/(2.0 * var))/SQRT2PI/sqrt(var));
}

/* -------------------------------------------------------------------
Computes the likelihoods for the mixed classes. Returns the likelihood.
var1,var2 are the variances of pdfs representing pure classes. measurement_var 
is the measurement noise. So the model for the variable y (representing the 
intensity value) that is composed of t * tissue1 and (1 - t)* tissue2 becomes :

y = t*x1 + (1 - t)*x2 + xm,
x1 ~ N(mean1,var1) , x2 ~ N(mean2,var2) , xm ~ N(0,measurement_var).

Note: The numerical integration routine used by the 
function is primitive , but so is the mankind...
Jussi Tohka
*/

double ComputeMarginalizedLikelihood(double value, double mean1, double mean2, 
    double var1, double var2, unsigned int nof_intervals)

{ 
    double lh, tmean, tvar, delta, step;
    
    step = 1.0 / (double) nof_intervals;
    lh = 0.0;
    
    for (delta = 0.0; delta <= 1.0; delta += step) {
        tmean = delta * mean1 + ( 1 - delta ) * mean2;
        tvar = SQR(delta) * var1 + SQR(1 - delta) * var2;
        lh += ComputeGaussianLikelihood(value, tmean, tvar)*step;
    }
    
    return(lh);
 }


/* Find maximum argument out of the n possibilities */
unsigned char MaxArg(double *val, unsigned char n)
{
    double maximum;
    unsigned char i, index;
    
    maximum = val[0];
    index = 1;
    
    for (i = 1; i < n; i++) {
        if (val[i] > maximum) {
            index = i + 1;
            maximum = val[i];
        }
    }
    return(index);
}

/* Normalize values to an overall sum of 1 */
void Normalize(double* val, char n)
{
    double sum_val = 0.0;
    int i;

    for (i = 0; i < n; i++) 
        sum_val += val[i];
 
    if (fabs(sum_val) > TINY) {    /* To avoid divisions by zero */
        for (i = 0; i < n; i++) {
            val[i] /= sum_val;
        }
    }
}

/* Compute initial PVE labeling based on marginalized likelihood for Amap */
void ComputeInitialPveLabelSub(float *src, unsigned char *label, unsigned char *prob, struct point *r, int n_pure_classes, int sub, int *dims)
{
    int x, y, z, z_area, y_dims, index, label_value;
    int i, ix, iy, iz, ind, ind2, nix, niy, niz, narea, nvol;
    int area, vol;
    double val, sub_1, mean[MAX_NC], var[MAX_NC], d_pve[MAX_NC];
    
    area = dims[0]*dims[1];
    vol = area*dims[2];

    /* find grid point conversion factor */
    sub_1 = 1.0/((double) sub);

    /* define grid dimensions */
    nix = (int) ceil((dims[0]-1)/((double) sub))+1;
    niy = (int) ceil((dims[1]-1)/((double) sub))+1;
    niz = (int) ceil((dims[2]-1)/((double) sub))+1; 

    narea = nix*niy;
    nvol = nix*niy*niz;
    
    /* loop over image points */
    for (z = 1; z < dims[2]-1; z++) {
        z_area=z*area;
        for (y = 1; y < dims[1]-1; y++) {
            y_dims=y*dims[0];
            for (x = 1; x < dims[0]-1; x++)  {
         
                index = x + y_dims + z_area;
                label_value = (int)label[index];
                if (label_value == 0) continue;
                val = (double)src[index];
                    
                /* find the interpolation factors */
                ix = (int)(sub_1*x), iy = (int)(sub_1*y), iz = (int)(sub_1*z);
                ind = iz*narea + iy*nix + ix;
                    
                for (i = 0; i < n_pure_classes; i++) {
                    ind2 = (i*nvol) + ind;                      
                    if (r[ind2].mean > 0.0) {
                        mean[i*2] = r[ind2].mean;
                        var[i*2]  = r[ind2].var;
                    }
                }

                if (fabs(mean[CSFLABEL-1]) > TINY) {
                    d_pve[CSFLABEL-1] = ComputeGaussianLikelihood(val, mean[CSFLABEL-1], var[CSFLABEL-1]);
                } else d_pve[CSFLABEL-1] = HUGE;

                if (fabs(mean[GMLABEL-1]) > TINY) {
                    d_pve[GMLABEL-1] = ComputeGaussianLikelihood(val, mean[GMLABEL-1], var[GMLABEL-1]);
                } else d_pve[GMLABEL-1] = HUGE;

                if (fabs(mean[WMLABEL-1]) > TINY) {
                    d_pve[WMLABEL-1] = ComputeGaussianLikelihood(val, mean[WMLABEL-1], var[WMLABEL-1]);
                } else d_pve[WMLABEL-1] = HUGE;

                if ((fabs(mean[WMLABEL-1]) > TINY) && (fabs(mean[GMLABEL-1]) > TINY)) {
                    d_pve[WMGMLABEL-1] = ComputeMarginalizedLikelihood(val, mean[WMLABEL-1], mean[GMLABEL-1],
                        var[WMLABEL-1], var[GMLABEL-1], 100 );
                } else d_pve[WMGMLABEL-1] = HUGE;
                        
                if ((fabs(mean[CSFLABEL-1]) > TINY) && (fabs(mean[GMLABEL-1]) > TINY)) {
                    d_pve[GMCSFLABEL-1] = ComputeMarginalizedLikelihood(val, mean[GMLABEL-1], mean[CSFLABEL-1],
                        var[GMLABEL-1], var[CSFLABEL-1], 100 );
                } else d_pve[GMCSFLABEL-1] = HUGE;
                
                Normalize(d_pve, n_pure_classes+2);
                
                for (i = 0; i < n_pure_classes+2; i++) 
                    prob[(vol*i) + index] = (unsigned char)ROUND(255.0*d_pve[i]);

                label[index] = (unsigned char) MaxArg(d_pve, n_pure_classes+2);
            }
        }
    }       
} 

void ComputeMrfProbability(double *mrf_probability, double *exponent, unsigned char *label, int x, int y, int z, int *dims,
              int n_classes, double beta, double *voxelsize_squared,
              const double *class_weights)
{
    int i,j,k;
    unsigned char label1, label2;  
    double distance;
    int similarity_value;
    int same = -2;
    int similar = -1;
    int different = 1; 
    
    /* To determine if it's possible to get out of image limits. 
         If not true (as it usually is) this saves the trouble calculating this 27 times */
    for (label1 = 0; label1 < n_classes; label1++)
        exponent[label1] = 0;
    
    for (i = -1; i < 2; i++) for (j = -1; j < 2; j++) for (k = -1; k < 2; k++) 
        if ( i != 0 || j != 0 || k != 0 ) {
                     
            label2 = label[(x+i)+dims[0]*(y+j)+dims[0]*dims[1]*(z+k)];
                             
            for (label1 = 1; label1 < n_classes+1; label1++) { 
                double class_weight = 1.0;
                if (class_weights) class_weight = class_weights[label1-1];
                if (class_weight < 0.0) class_weight = 0.0;
                if (label1 == label2) similarity_value = same;
                else if (abs(label1 - label2) < 2) similarity_value = similar;
                else similarity_value = different;

                distance = sqrt(voxelsize_squared[0] * abs(i) + voxelsize_squared[1] * abs(j) + voxelsize_squared[2] * abs(k));

                exponent[label1-1] += (class_weight * (double)similarity_value) / distance;
            }
        }       

    for (label1 = 0; label1 < n_classes; label1++)
        mrf_probability[label1] = exp(-(beta * exponent[label1])); 
    
} 

/* Iterative conditional mode */
void ICM(unsigned char *prob, unsigned char *label, int n_classes, int *dims, double beta, int iterations,
         double *voxelsize, int verbose, const double *class_weights)
{
    
    int i, iter, x, y, z, z_area, y_dims, index, sum_voxel;
    long area, vol;
    double rel_changed, mrf_probability[MAX_NC], voxelsize_squared[3];
    double exponent[MAX_NC], sum_voxelsize = 0.0;
    unsigned char new_label;
        
    area = dims[0]*dims[1];
    vol = area*dims[2];
    
    /* normalize voxelsize to a sum of 3 and calculate its squared value */
    for (i = 0; i < 3; i++) sum_voxelsize += voxelsize[i];
    for (i = 0; i < 3; i++) voxelsize_squared[i] = SQR(3.0*voxelsize[i]/sum_voxelsize);
        
    for (iter=0; iter < iterations; iter++) {
        sum_voxel = 0;
        rel_changed = 0.0;
        
        /* loop over image points */
        for (z = 1; z < dims[2]-1; z++) {
            z_area=z*area;
            for (y = 1; y < dims[1]-1; y++) {
                y_dims=y*dims[0];
                for (x = 1; x < dims[0]-1; x++)  {
         
                    index = x + y_dims + z_area;
                    if (label[index] == 0) continue;
                    
                    sum_voxel++;
                    ComputeMrfProbability(mrf_probability, exponent, label, x, y, z, dims, n_classes, beta,
                                          voxelsize_squared, class_weights);
                    
                    for (i = 0; i < n_classes; i++)
                        mrf_probability[i] *= (double)prob[index+i*vol];

                    new_label = (unsigned char) MaxArg(mrf_probability, n_classes);
                    if (new_label != label[index]) {
                        rel_changed += 1.0;
                        label[index] = new_label;            
                    }
                }
            }
        }

        rel_changed /= (double)sum_voxel;
#if !defined(_WIN32) && !defined(_WIN64)
        if (verbose) {
            printf("ICM: %d relative change: %2.4f\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b",iter+1, 100.0*rel_changed);
            fflush(stdout);
        }
#endif
        if (rel_changed < TH_CHANGE) break;
    }   
    if (verbose) printf("\n");
} 

void EstimateSegmentation(float *src, unsigned char *label, unsigned char *prob, 
                          struct point *r, struct point *r_large, double *mean, double *var, int n_classes, 
                          int niters, int sub, int sub_large, int *dims, double *voxelsize, double *thresh, 
                          double *beta, int verbose, int use_median)
{
    int i;
    int area, narea, nvol, narea_large, nvol_large, vol, z_area, y_dims, index, ind;
    double sub_1, sub_1_large, dmin, val, sum_voxelsize = 0.0;
    double d[MAX_NC], alpha[MAX_NC], log_alpha[MAX_NC], log_var[MAX_NC], exponent[MAX_NC];
    double pvalue[MAX_NC], mrf_probability[MAX_NC], voxelsize_squared[3], psum;
    int nix, niy, niz, nix_large, niy_large, niz_large, iters, count_change;
    int x, y, z, label_value, xBG;
    int ix, iy, iz, ix_large, iy_large, iz_large, ind2, ind2_large, replace = 1;
    double ll, ll_old, change_ll;
        
    MrfPrior(label, n_classes, alpha, beta, 0, dims, verbose);
    
    area = dims[0]*dims[1];
    vol = area*dims[2];

    /* find grid point conversion factor */
    sub_1 = 1.0/((double) sub);
    sub_1_large = 1.0/((double) sub_large);

    /* define grid dimensions */
    nix = (int) ceil((dims[0]-1)/((double) sub))+1;
    niy = (int) ceil((dims[1]-1)/((double) sub))+1;
    niz = (int) ceil((dims[2]-1)/((double) sub))+1; 

    nix_large = (int) ceil((dims[0]-1)/((double) sub_large))+1;
    niy_large = (int) ceil((dims[1]-1)/((double) sub_large))+1;
    niz_large = (int) ceil((dims[2]-1)/((double) sub_large))+1; 

    narea = nix*niy;
    nvol = nix*niy*niz;
    narea_large = nix_large*niy_large;
    nvol_large = nix_large*niy_large*niz_large;

    for (i = 0; i < n_classes; i++) log_alpha[i] = log(alpha[i]);

    /* normalize voxelsize to a sum of 3 and calculate its squared value */
    for (i = 0; i < 3; i++) sum_voxelsize += voxelsize[i];
    for (i = 0; i < 3; i++) voxelsize_squared[i] = SQR(3.0*voxelsize[i]/sum_voxelsize);
        
    ll_old = HUGE;
    count_change = 0;
    
    for (iters = 0; iters < niters; iters++) {
            
        ll = 0.0;
        
        /* get means for grid points */
        GetMeansVariances(src, label, n_classes, r, sub, dims, thresh, use_median);
        if (r_large && sub_large != sub)
            GetMeansVariances(src, label, n_classes, r_large, sub_large, dims, thresh, use_median);

        /* loop over image points */
        for (z = 1; z < dims[2]-1; z++) {
            z_area=z*area;
            for (y = 1; y < dims[1]-1; y++) {
                y_dims=y*dims[0];
                for (x = 1; x < dims[0]-1; x++)  {
         
                    index = x + y_dims + z_area;
                    label_value = (int) label[index];
                    if (label_value < 1) continue;
                    val = (double)src[index];
                    
                    /* find the interpolation factors */
                    ix = (int)(sub_1*x);
                    iy = (int)(sub_1*y);
                    iz = (int)(sub_1*z);
                    ind = iz*narea + iy*nix + ix;

                    ix_large = (int)(sub_1_large*x);
                    iy_large = (int)(sub_1_large*y);
                    iz_large = (int)(sub_1_large*z);
                    ind2_large = iz_large*narea_large + iy_large*nix_large + ix_large;
                    
                    for (i = 0; i < n_classes; i++) {
                        ind2 = (i*nvol) + ind;
                        if (r[ind2].mean > TINY) {
                            mean[i] = r[ind2].mean;
                            var[i]  = r[ind2].var;
                        } else if (r_large && sub_large != sub) {
                            int ind_large = (i*nvol_large) + ind2_large;
                            if (r_large[ind_large].mean > TINY) {
                                mean[i] = r_large[ind_large].mean;
                                var[i]  = r_large[ind_large].var;
                            }
                        }

                        if (var[i] > TINY)
                            log_var[i] = log(var[i]);
                    }
                    
                    /* compute energy at each point */
                    dmin = HUGE; xBG = 1; 
                    psum = 0.0;

                    for (i = 0; i < n_classes; i++) {
                        if (fabs(mean[i]) > TINY) {

                            d[i] = 0.5*(SQR(val-mean[i])/var[i]+log_var[i])-log_alpha[i];
                            pvalue[i] = exp(-d[i])/SQRT2PI;
                            psum += pvalue[i];
                            
                        } else d[i] = HUGE;
                        
                        if (d[i] < dmin) {
                            dmin = d[i];
                            xBG = i;
                        }
                    }

                    /* scale p-values to a sum of 1 */
                    if (psum > TINY) {
                        for (i = 0; i < n_classes; i++) pvalue[i] /= psum;
                        ll -= log(psum);
                    } else for (i = 0; i < n_classes; i++) pvalue[i] = 0.0;
                 
                    for (i = 0; i < n_classes; i++)
                        prob[(vol*i) + index] = (unsigned char)ROUND(255*pvalue[i]);
                 
                    /* if the class has changed modify the label */
                    if (xBG + 1 != label_value) label[index] = (unsigned char) (xBG + 1); 
                 
                }
            }
        }

        ll /= (double)vol;
        change_ll = (ll_old - ll)/fabs(ll);

        /* break if log-likelihood has not changed significantly two iterations */
        if (change_ll < TH_CHANGE) count_change++;
        if ((count_change > 2) && (iters > 5)) break;

#if !defined(_WIN32) && !defined(_WIN64)
        if (verbose) {
            printf("iters:%3d log-likelihood: %7.5f\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b",iters+1, ll);
            fflush(stdout);
        }
#endif
        ll_old = ll;
    }

    if (verbose) {
        printf("\nFinal Mean*Std: ");
        for (i = 0; i < n_classes; i++) printf("%.3f*%.3f    ",mean[i],sqrt(var[i])); 
        printf("\n"); 
    }
    
}

/* perform adaptive MAP on given src and initial segmentation label */
void Amap(float *src, unsigned char *label, unsigned char *prob, double *mean, 
          int n_classes, int niters, int sub, int *dims, int pve, double weight_MRF, 
          double *voxelsize, int niters_ICM, int verbose, 
          int use_median, const double *mrf_class_weights)
{
    int i, nix, niy, niz;
    int area, nvol, vol;
    int n[MAX_NC], j;
    double var[MAX_NC], mean_voxelsize;
    double thresh[2], beta[1];
    struct point *r = NULL;
    struct point *r_large = NULL;
    int sub_large;
    double class_weights[MAX_NC];
    
    /* we have to make sub independent from voxel size */
    mean_voxelsize = (voxelsize[0] + voxelsize[1] + voxelsize[2])/3.0;
    sub = ROUND((double)sub/mean_voxelsize);
    
    /* ICM is not needed if we skip MRF approach */
    if (weight_MRF == 0)
        niters_ICM = 0;

    sub_large = sub * 2;
    if (sub_large < sub + 1) sub_large = sub + 1;
    sub_large = MIN(sub_large, MAX(dims[0], MAX(dims[1], dims[2])));
            
    area = dims[0]*dims[1];
    vol = area*dims[2];
 
    /* define grid dimensions */
    nix = (int) ceil((dims[0]-1)/((double) sub))+1;
    niy = (int) ceil((dims[1]-1)/((double) sub))+1;
    niz = (int) ceil((dims[2]-1)/((double) sub))+1; 
    nvol  = nix*niy*niz;
    
    r = (struct point*)malloc(sizeof(struct point)*MAX_NC*nvol);
    if (r == NULL) {
        printf("Memory allocation error\n");
        exit(EXIT_FAILURE);
    }

    /* allocate larger-grid stats for adaptive sub */
    {
        int nix_l = (int) ceil((dims[0]-1)/((double) sub_large))+1;
        int niy_l = (int) ceil((dims[1]-1)/((double) sub_large))+1;
        int niz_l = (int) ceil((dims[2]-1)/((double) sub_large))+1; 
        int nvol_l  = nix_l*niy_l*niz_l;

        r_large = (struct point*)malloc(sizeof(struct point)*MAX_NC*nvol_l);
        if (r_large == NULL) {
            printf("Memory allocation error\n");
            exit(EXIT_FAILURE);
        }
    }
    
    double max_label = get_max(label, vol, 0, DT_UINT8);
    if (max_label > n_classes) {
        printf("Label maximum %g exceeds number of tissue classes.\n", max_label);
        exit(EXIT_FAILURE);
    }

    /* estimate 3 classes before PVE */
    EstimateSegmentation(src, label, prob, r, r_large, mean, var, n_classes, niters, sub, sub_large, dims, voxelsize, 
                            thresh, beta, verbose, use_median);

    /* Use marginalized likelihood to estimate initial 5 classes */
    if (pve) {
        ComputeInitialPveLabelSub(src, label, prob, r, n_classes, sub, dims);
        n_classes = 5;

        /* recalculate means for pure and mixed classes */
        for (j = 0; j < n_classes; j++) {
            n[j] = 0;
            mean[j] = 0.0;
        }
        for (i = 0; i < vol; i++) {
            if (label[i] == 0) continue;
            n[label[i]-1]++;
            mean[label[i]-1] += (double)src[i];
        }
        for (j = 0; j < n_classes; j++) mean[j] /= n[j];
    }

    /* use much smaller beta for if no pve is selected */
    if (!pve) beta[0] /= 20.0;
    
    /* Iterative Conditional Mode (run after PVE initialization) */
    if ((niters_ICM > 0) && (weight_MRF > 0.0)) {
        for (i = 0; i < MAX_NC; i++) class_weights[i] = 1.0;
        if (mrf_class_weights) {
            for (i = 0; i < n_classes; i++) class_weights[i] = mrf_class_weights[i];
        }

        beta[0] *= weight_MRF;
        if (verbose) printf("Weighted MRF beta %3.3f\n",beta[0]);
        ICM(prob, label, n_classes, dims, beta[0], niters_ICM, voxelsize, verbose,
            class_weights);
    }
    
    free(r);
    if (r_large) free(r_large);

    return;      
}
