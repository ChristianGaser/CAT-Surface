/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include "CAT_Amap.h"
#include "CAT_Bmap.h"
#include "CAT_Vol.h"
#include "CAT_Math.h"

void get_order(const double mean[3], int order[3]) {
    // Find the order of mean[0], mean[1], mean[2]
    if (mean[0] >= mean[1] && mean[0] >= mean[2]) {
        order[0] = 0;
        if (mean[1] >= mean[2]) { order[1] = 1; order[2] = 2; }
        else { order[1] = 2; order[2] = 1; }
    } else if (mean[1] >= mean[0] && mean[1] >= mean[2]) {
        order[0] = 1;
        if (mean[0] >= mean[2]) { order[1] = 0; order[2] = 2; }
        else { order[1] = 2; order[2] = 0; }
    } else { // mean[2] is largest
        order[0] = 2;
        if (mean[0] >= mean[1]) { order[1] = 0; order[2] = 1; }
        else { order[1] = 1; order[2] = 0; }
    }
}

bool is_order_changed(const double old_mean[3], const double new_mean[3]) {
    int old_order[3], new_order[3], i;
    get_order(old_mean, old_order);
    get_order(new_mean, new_order);

    // Compare orders
    for (i = 0; i < 3; ++i) {
        if (old_order[i] != new_order[i]) return true; // order has changed
    }
    return false; // order is unchanged
}

/* Compute initial PVE labeling based on marginalized likelihood for Bmap*/
void ComputeInitialPveLabel(float *src, unsigned char *label, unsigned char *prob, double *mean, double *var, int n_pure_classes, int *dims)
{
    int x, y, z, z_area, y_dims, index, label_value;
    int i, ix, iy, iz, ind, ind2;
    int area, vol;
    double val, mean_pve[MAX_NC], var_pve[MAX_NC], d_pve[MAX_NC];
    
    area = dims[0]*dims[1];
    vol = area*dims[2];

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
                    
                for (i = 0; i < n_pure_classes; i++) {
                    mean_pve[i*2] = mean[i];
                    var_pve[i*2]  = var[i];
                }

                if (fabs(mean_pve[CSFLABEL-1]) > TINY) {
                    d_pve[CSFLABEL-1] = ComputeGaussianLikelihood(val, mean_pve[CSFLABEL-1], var_pve[CSFLABEL-1]);
                } else d_pve[CSFLABEL-1] = HUGE;

                if (fabs(mean_pve[GMLABEL-1]) > TINY) {
                    d_pve[GMLABEL-1] = ComputeGaussianLikelihood(val, mean_pve[GMLABEL-1], var_pve[GMLABEL-1]);
                } else d_pve[GMLABEL-1] = HUGE;

                if (fabs(mean_pve[WMLABEL-1]) > TINY) {
                    d_pve[WMLABEL-1] = ComputeGaussianLikelihood(val, mean_pve[WMLABEL-1], var_pve[WMLABEL-1]);
                } else d_pve[WMLABEL-1] = HUGE;

                if ((fabs(mean_pve[WMLABEL-1]) > TINY) && (fabs(mean_pve[GMLABEL-1]) > TINY)) {
                    d_pve[WMGMLABEL-1] = ComputeMarginalizedLikelihood(val, mean_pve[WMLABEL-1], mean_pve[GMLABEL-1],
                        var_pve[WMLABEL-1], var_pve[GMLABEL-1], 100 );
                } else d_pve[WMGMLABEL-1] = HUGE;
                        
                if ((fabs(mean_pve[CSFLABEL-1]) > TINY) && (fabs(mean_pve[GMLABEL-1]) > TINY)) {
                    d_pve[GMCSFLABEL-1] = ComputeMarginalizedLikelihood(val, mean_pve[GMLABEL-1], mean_pve[CSFLABEL-1],
                        var_pve[GMLABEL-1], var_pve[CSFLABEL-1], 100 );
                } else d_pve[GMCSFLABEL-1] = HUGE;
                
                Normalize(d_pve, n_pure_classes+2);
                
                for (i = 0; i < n_pure_classes+2; i++) 
                    prob[(vol*i) + index] = (unsigned char)ROUND(255.0*d_pve[i]);

                label[index] = (unsigned char) MaxArg(d_pve, n_pure_classes+2);
            }
        }
    }       
} 

void xaverage(unsigned char *label, float *bias, long *n1, double *bs, int bg, int a, int *dims)
{
    int x,y,z,r;

    int area = dims[0]*dims[1];
    for (z = 0; z < dims[2]; z++) {
        int i = z*area;
        for (y = 0; y < dims[1]; y++) {
            for (x = 0; x < dims[0]; x++) {
                i++; n1[i] = 0;  bs[i] = 0;
                int l = x-a < 0 ? 0 : x-a;
                int h = x+a > dims[0]-1 ? dims[0]-1 : x+a;
                for (r=l; r <=h  ; r++) {
                    int j = i + (r-x);
                    if ((int) label[j]>=bg) {
                        n1[i]++;
                        bs[i] += (double)bias[j];
                    }
                }
            }
         }
    }
}

void yaverage(long *n1, double *bs, long *n2, float *bias, int b, int *dims)
{
    int x,y,z,r;

    int area = dims[0]*dims[1];

    for (z=0; z<dims[2]; z++) {
        int i = z*area;
        for (y = 0; y < dims[1]; y++) {
            for (x = 0; x < dims[0]; x++) {
                i++; n2[i] = 0; bias[i] = 0.0;
                int l = y-b < 0 ? 0 : y-b;
                int h = y+b > dims[1]-1 ? dims[1]-1 : y+b;
                for (r=l; r <=h  ; r++) {
                    int j = i + (r-y)*dims[0];
                    n2[i] += n1[j];
                    bias[i] += (float)bs[j];
                }
            }
         }
    }
}

void zaverage(long *n2, float *bias, long *n1, double *bs, int bg, int c, int *dims)
{
    int x,y,z,r;

    int area = dims[0]*dims[1];

    for (z = 0; z < dims[2]; z++) {
        int i = z*area;
        for (y = 0; y < dims[1]; y++) {
            for (x = 0; x < dims[0]; x++) {
                i++; bs[i] = 0; n1[i] = 0;
                int l = z-c < 0 ? 0 : z-c;
                int h = z+c > dims[2]-1 ? dims[2]-1 : z+c;
                for (r=l; r <=h  ; r++) {
                    int j = i + (r-z)*area;
                    n1[i] += n2[j];
                    bs[i] += (double)bias[j];
                }
            }
         }
    }
}

void movingAverage(unsigned char *label, float *bias, int BG, int a, int b, int c, int *dims)
{
    int i, j, x, y, z, r, l, h;
    double *bs;
    long *n1, *n2;

    long vol = dims[0]*dims[1]*dims[2];
    
    bs = (double *) malloc(sizeof(double)*vol);
    n1 = (long *) malloc(sizeof(long)*vol);
    n2 = (long *) malloc(sizeof(long)*vol);
    
    xaverage(label, bias, n1, bs, BG, a, dims);
    yaverage(n1, bs, n2, bias, b, dims);
    zaverage(n2, bias, n1, bs, BG, c, dims);
        
    for (i=0; i<vol; i++) {
        if ((int) label[i] >= BG) bias[i] =    n1[i] > 0 ? bs[i]/n1[i] : 0.0;
        else bias[i] = 0.0;
    }
    
    free(bs);
    free(n1);
    free(n2);

    double voxel_size[3] = {1.0, 1.0, 1.0};
    double fwhm[3] = {50.0, 50.0, 50.0};

    //smooth_subsample3(bias, dims, voxel_size, fwhm, 0, 8, DT_FLOAT32); 

    return;
}

void Bmap(float *src, unsigned char *label, unsigned char *prob, double *mean, 
        int n_classes, int BG, int niters, int a, int b, int c, 
        float *bias, int *dims, int pve, int verbose)
{
    int i, j, index, x, y, z, iters, count_change;
    int n[MAX_NC];
    double thresh[2], val, d;
    double alpha[n_classes], s[n_classes], ss[n_classes], mean_old[n_classes];
    double var[n_classes], lvar[n_classes], p[n_classes];
    double ll, ll_old, change_ll;
    
    int area = dims[0]*dims[1];
    int vol = area*dims[2];

    double prctile[2] = {1,99};
    get_prctile(src, vol, thresh, prctile, 1, DT_FLOAT32);  

    // initialize means
    for (j=0; j<n_classes; j++) s[j]=ss[j]=0.0;
    for (i=0; i<vol; i++) {
        int lab = (int) label[i];
        if (lab < BG) continue;
        s[lab-BG] += 1.0; ss[lab-BG] += (double)src[i];
    }
    for (j=0; j<n_classes; j++)
        mean[j] = s[j] > 0.0 ? ss[j]/s[j]: 0.0;

    // intitialize standard deviations
    for (j=0; j<n_classes; j++) ss[j] = 0.0;
    for (i=0; i<vol; i++) {
        int lab = (int) label[i];
        if (lab < BG) continue;
        ss[lab-BG] += SQR((double)src[i]-mean[lab-BG]);
    }
    for (j=0; j<n_classes; j++) 
        var[j] = s[j]>1.0 ? ss[j]/(s[j]-1.0): 1.0;

    for (i=0; i<vol; i++) bias[i] = 0.0;
    
    // initialize prior parameters
    MrfPrior(label, n_classes, alpha, NULL, 0, dims, 0);

    if (verbose) {
        fprintf(stdout,"Initial means*vars: "); 
        for (i = 0; i < n_classes; i++) fprintf(stdout,"%.3f*%.3f\t",mean[i],sqrt(var[i]));
        fprintf(stdout,"\n"); 
    }

    // set new variables to speed up
    for (j=0; j<n_classes; j++) {
        lvar[j] = var[j] > 0.0 ? 0.5*log(var[j]) - log(alpha[j]): - log(alpha[j]);
        var[j] = 0.5/var[j];
    }

    ll_old = HUGE;
    count_change = 0;
    
    // iterative condition modes
    for (iters = 0; iters <= niters; iters++)    {
        ll = 0.0;

        // loop over image voxels
        for (z = 1; z < dims[2]-1; z++) 
            for (y = 1; y < dims[1]-1; y++)
                for (x = 1; x < dims[0]-1; x++)  {
        
                    index = x + y*dims[0] + z*area;

                    int lab = (int) label[index];
                    if (lab < BG) continue;
        
                    // loop over all classes
                    int xi = BG; double dmin = HUGE;
                    double psum = 0.0;
                    
                    for (j = 0; j < n_classes; j++) {

                        d = SQR((double)src[index]-(1.0+(double)bias[index])*mean[j])*var[j]+lvar[j];
            
                        p[j] = exp(-d)/(SQRT2PI*sqrt(var[j]));
                        psum += p[j];
            
                        if (d < dmin) {xi = j; dmin = d;}
                    }

                    if (psum > TINY) {
                        for (j = 0; j < n_classes; j++)
                            prob[(vol*j) + index] = (unsigned char)round(255*p[j]/psum);
                        ll -= log(psum);
                    } else for (j = 0; j < n_classes; j++) prob[(vol*j) + index] = 0.0;

                    // change the label
                    if (xi+BG != lab) label[index] = (unsigned char) (xi+BG);
        
                    // find bias values
                    bias[index] = src[index]/mean[xi] - 1.0;
                }

        ll /= (double)vol;
        change_ll = (ll_old - ll)/fabs(ll);
                
        // continue if log-likehihood increases
        if (change_ll <= 0) continue;

        // break if log-likelihood has not changed significantly two iterations
        if (change_ll < TH_CHANGE) count_change++;
        if ((count_change > 2) && (iters > 20)) break;

        // smoothout bias values (remove misclassifications) only in first iterations
        if (iters < 1)
            movingAverage(label, bias, BG, a, b, c, dims);
        
        for (j = 0; j < n_classes; j++) { s[j] = ss[j] = 0.0; }
        
        for (i = 0; i < vol; i++) {
            int lab = (int) label[i];
            if (lab < BG) continue;
            double val = (1.0 + bias[i]) * (double)src[i];
    
            s[lab - BG] += 1.0;
            ss[lab - BG] += val;
        }
        
        int mean_too_close = 0;
        for (j = 0; j < n_classes; j++) {
            mean[j] = (s[j] > 0.0) ? (ss[j] / s[j]) : 0.0;
            if (j > 0) {
                double mean_diff = fabs(mean[j] - mean[j-1]);
                if (mean_diff < 0.05) mean_too_close = 1;
            }
        }
        
        if ((iters > 0) & ((mean_too_close) | (is_order_changed(mean, mean_old)))) {
            for (j = 0; j < n_classes; j++) mean[j] = mean_old[j]; 
            break;
        }
        
        for (j = 0; j < n_classes; j++) { s[j] = ss[j] = 0.0; }
        
        for (i = 0; i < vol; i++) {
            int lab = (int) label[i];
            if (lab < BG) continue;
            double val = (1.0 + bias[i]) * (double)src[i];
    
            s[lab - BG] += 1.0;
            ss[lab - BG] += SQR(val - mean[lab - BG]);
        }
        
        for (j = 0; j < n_classes; j++) {
            var[j] = (s[j] > 0.0) ? (ss[j] / s[j]) : 1.0;
            lvar[j] = var[j] > 0.0 ? 0.5*log(var[j]) - log(alpha[j]): -log(alpha[j]);
            var[j] = 0.5/var[j];
        }

#if !defined(_WIN32) && !defined(_WIN64)
        if (verbose) {
            printf("iters:%2d log-likelihood: %7.5f ", iters, ll);
            for (i = 0; i < n_classes; i++) printf("%.3f*%.3f ",mean[i],sqrt(0.5/var[i]));
            printf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b");
            fflush(stdout);
        }
#endif

        ll_old = ll;
        for (j = 0; j < n_classes; j++) mean_old[j] = mean[j]; 

    }
        
    for (j = 0; j < n_classes; j++) var[j] = 0.5/var[j];

    if (verbose) {
        fprintf(stdout,"\nFinal means*vars: "); 
        for (i = 0; i < n_classes; i++) fprintf(stdout,"%.3f*%.3f\t",mean[i],sqrt(var[i]));
        fprintf(stdout,"\n"); 
    }
    
    /* Use marginalized likelihood to estimate initial 5 classes */
    if (pve) {
        ComputeInitialPveLabel(src, label, prob, mean, var, n_classes, dims);
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

    i = 0;
    for (z=0; z<dims[2]; z++) {
        for (y=0; y<dims[1]; y++) 
            for (x=0; x<dims[0]; x++, i++) {
            int lab = (int) label[i];
            if (lab < BG) continue;
            bias[i] = (double)src[i] -(1.0+bias[i])*mean[lab-BG];
        }
    }

    return;
}    
