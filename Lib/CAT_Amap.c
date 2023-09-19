/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "CAT_Amap.h"

/* This PVE calculation is a modified version from
 * the PVE software bundle:
 * Copyright (C) Jussi Tohka, Institute of Signal Processing, Tampere University of
 * Technology, 2002 - 2004.
 * P.O. Box 553, FIN-33101, Finland
 * E-mail: jussi.tohka@tut.fi
 */

void Pve5(float *src, unsigned char *prob, unsigned char *label, double *mean, int *dims)
{
        int x,y,z,i,z_area,y_dims,ind,mxi;
        double w, mx;
        unsigned char new_val[MAX_NC];
        
        int area = dims[0]*dims[1];
        int vol = area*dims[2];
                
        for (z = 1; z < dims[2]-1; z++) {
                z_area = z*area;
                for (y = 1; y < dims[1]-1; y++) {
                        y_dims = y*dims[0];
                        for (x = 1; x < dims[0]-1; x++) {
                                ind = z_area + y_dims + x;

                                switch(label[ind]) {
                                case 0: /* BG */
                                        new_val[CSFLABEL-1] = 0;
                                        new_val[GMLABEL-1]  = 0;
                                        new_val[WMLABEL-1]  = 0;
                                        break;
                                case CSFLABEL: /* CSF */
                                        new_val[CSFLABEL-1] = 255;
                                        new_val[GMLABEL-1]  = 0;
                                        new_val[WMLABEL-1]  = 0;
                                        label[ind] = (unsigned char) ROUND(255.0/3.0);
                                        break;
                                case GMLABEL: /* GM */
                                        new_val[CSFLABEL-1] = 0;
                                        new_val[GMLABEL-1]  = 255;
                                        new_val[WMLABEL-1]  = 0;
                                        label[ind] = (unsigned char) ROUND(2.0*255.0/3.0);
                                        break;
                                case WMLABEL: /* WM */
                                        new_val[CSFLABEL-1] = 0;
                                        new_val[GMLABEL-1]  = 0;
                                        new_val[WMLABEL-1]  = 255;
                                        label[ind] = 255;
                                        break;
                                case GMCSFLABEL: /* GMCSF */
                                        w = ((double)src[ind] - mean[CSFLABEL-1])/(mean[GMLABEL-1]-mean[CSFLABEL-1]);
                                        if (w > 1.0) w = 1.0; if (w < 0.0) w = 0.0;
                                        new_val[CSFLABEL-1] = (unsigned char) ROUND(255.0*(1-w));
                                        new_val[GMLABEL-1]  = (unsigned char) ROUND(255.0*w);
                                        new_val[WMLABEL-1]  = 0;
                                        label[ind] = ROUND(255.0/3.0*(1.0 + w));
                                        break;
                                case WMGMLABEL: /* WMGM */
                                        w = ((double)src[ind] - mean[GMLABEL-1])/(mean[WMLABEL-1]-mean[GMLABEL-1]);
                                        if (w > 1.0) w = 1.0; if (w < 0.0) w = 0.0;
                                        new_val[CSFLABEL-1] = 0;
                                        new_val[GMLABEL-1]  = (unsigned char) ROUND(255.0*(1-w));
                                        new_val[WMLABEL-1]  = (unsigned char) ROUND(255.0*w);
                                        label[ind] = ROUND(255.0/3.0*(2.0 + w));
                                        break;
                                }

                                prob[          ind] = new_val[CSFLABEL-1];
                                prob[vol     + ind] = new_val[GMLABEL-1];
                                prob[(2*vol) + ind] = new_val[WMLABEL-1];
                                
                                /* set old probabilities for mixed classes to zero */
                                prob[(3*vol) + ind] = 0;
                                prob[(4*vol) + ind] = 0;
                                
                        }
                }
        }        
}

void Pve6(float *src, unsigned char *prob, unsigned char *label, double *mean, int *dims)
{
        int x,y,z,i,z_area,y_dims,ind,mxi;
        double w, mx;
        unsigned char new_val[MAX_NC];
        
        int area = dims[0]*dims[1];
        int vol = area*dims[2];
                
        for (z = 1; z < dims[2]-1; z++) {
                z_area = z*area;
                for (y = 1; y < dims[1]-1; y++) {
                        y_dims = y*dims[0];
                        for (x = 1; x < dims[0]-1; x++) {
                                ind = z_area + y_dims + x;

                                switch(label[ind]) {
                                case 0: /* BG */
                                        new_val[CSFLABEL] = 0;
                                        new_val[GMLABEL]  = 0;
                                        new_val[WMLABEL]  = 0;
                                        break;
                                case CSFLABEL+1: /* CSF */
                                        new_val[CSFLABEL] = 255;
                                        new_val[GMLABEL]  = 0;
                                        new_val[WMLABEL]  = 0;
                                        label[ind] = (unsigned char) ROUND(255.0/3.0);
                                        break;
                                case GMLABEL+1: /* GM */
                                        new_val[CSFLABEL] = 0;
                                        new_val[GMLABEL]  = 255;
                                        new_val[WMLABEL]  = 0;
                                        label[ind] = (unsigned char) ROUND(2.0*255.0/3.0);
                                        break;
                                case WMLABEL+1: /* WM */
                                        new_val[CSFLABEL] = 0;
                                        new_val[GMLABEL]  = 0;
                                        new_val[WMLABEL]  = 255;
                                        label[ind] = 255;
                                        break;
                                case BKGCSFLABEL+1: /* BKGCSF */
                                        w = (double)src[ind]/mean[CSFLABEL];
                                        if (w > 1.0) w = 1.0; if (w < 0.0) w = 0.0;
                                        new_val[CSFLABEL] = (unsigned char) ROUND(255.0*w);
                                        new_val[GMLABEL]  = 0;
                                        new_val[WMLABEL]  = 0;
                                        label[ind] = ROUND(255.0/3.0*w);
                                        break;
                                case GMCSFLABEL+1: /* GMCSF */
                                        w = ((double)src[ind] - mean[CSFLABEL])/(mean[GMLABEL]-mean[CSFLABEL]);
                                        if (w > 1.0) w = 1.0; if (w < 0.0) w = 0.0;
                                        new_val[CSFLABEL] = (unsigned char) ROUND(255.0*(1-w));
                                        new_val[GMLABEL]  = (unsigned char) ROUND(255.0*w);
                                        new_val[WMLABEL]  = 0;
                                        label[ind] = ROUND(255.0/3.0*(1.0 + w));
                                        break;
                                case WMGMLABEL+1: /* WMGM */
                                        w = ((double)src[ind] - mean[GMLABEL])/(mean[WMLABEL]-mean[GMLABEL]);
                                        if (w > 1.0) w = 1.0; if (w < 0.0) w = 0.0;
                                        new_val[CSFLABEL] = 0;
                                        new_val[GMLABEL]  = (unsigned char) ROUND(255.0*(1-w));
                                        new_val[WMLABEL]  = (unsigned char) ROUND(255.0*w);
                                        label[ind] = ROUND(255.0/3.0*(2.0 + w));
                                        break;
                                }

                                prob[          ind] = new_val[CSFLABEL-1];
                                prob[vol     + ind] = new_val[GMLABEL-1];
                                prob[(2*vol) + ind] = new_val[WMLABEL-1];
                                
                                /* set old probabilities for mixed classes to zero */
                                prob[(3*vol) + ind] = 0;
                                prob[(4*vol) + ind] = 0;
                                prob[(5*vol) + ind] = 0;
                                
                        }
                }
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

void MrfPrior(unsigned char *label, int n_classes, double *alpha, double *beta, int init, int *dims)
{
        int i, j, k, x, y, z;
        int fi, fj;
        int verb = 0;
        long color[MAX_NC][7][7][7][7];
        long area;

        int f[MAX_NC-1], plab, zero, iBG;
        int n, z_area, y_dims;
        double XX, YY, L;
        
        area = dims[0]*dims[1];

        /* initialize configuration counts */
        for (i = 0; i < n_classes; i++)
                for (f[0] = 0; f[0] < 7; f[0]++)
                        for (f[1] = 0; f[1] < 7; f[1]++)
                                for (f[2] = 0; f[2] < 7; f[2]++)
                                        for (f[3] = 0; f[3] < 7; f[3]++)
                                                color[i][f[0]][f[1]][f[2]][f[3]]=0;

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

        /* evaluate alphas */
        if (verb) printf("MRF priors: alpha ");
        for (i = 0; i < n_classes; i++) {
                if (init == 0) alpha[i] /= n; else alpha[i] = 1.0;
                if (verb) printf("%3.3f ", alpha[i]);
        }

        /* compute beta */
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
        if (verb) printf("\t beta %3.3f\n", beta[0]);
        fflush(stdout);
}

/* calculate the mean and variance for every class on a grid size SUBxSUBxSUB */
static void GetMeansVariances(float *src, unsigned char *label, int n_classes, struct point *r, int sub, int *dims, double *thresh)
{
        int i, j, ind;
        int area, narea, nvol, zsub, ysub, xsub, yoffset, zoffset;
        int zsub2, ysub2;
        int nix, niy, niz, k, l, m, z, y, x, label_value;
        double val;
        struct ipoint *ir;
        int label_value_BG;

        area = dims[0]*dims[1];

        /* define grid dimensions */
        nix = (int) ceil((dims[0]-1)/((double) sub))+1;
        niy = (int) ceil((dims[1]-1)/((double) sub))+1;
        niz = (int) ceil((dims[2]-1)/((double) sub))+1; 
        narea = nix*niy;
        nvol    = nix*niy*niz;
        
        ir = (struct ipoint*)malloc(sizeof(struct ipoint)*n_classes*nvol);
        if (ir == NULL) {
                printf("Memory allocation error\n");
                exit(EXIT_FAILURE);
        }

        for(i = 0; i < n_classes; i++) {
                for(j = 0; j < nvol; j++) {
                        ind = (i*nvol)+j; 
                        ir[ind].n = 0;
                        ir[ind].s = 0.0;
                        ir[ind].ss = 0.0;
                }
        }
        

        /* loop over neighborhoods of the grid points */
        for(k=-sub; k<=sub; k++) for(l=-sub; l<=sub; l++) for(m=-sub; m<=sub; m++) 
                for(z = 0; z < niz; z++) {
                        zsub = z*sub + k;
                        if ((zsub >= 0) && (zsub < dims[2])) {
                                zsub2 = zsub*area;
                                zoffset = z*narea;
                                for(y = 0; y < niy; y++) {
                                        ysub = y*sub + l;
                                        if ((ysub >= 0) && (ysub < dims[1])) {
                                                ysub2 = ysub*dims[0];
                                                yoffset = zoffset + y*nix;
                                                for(x = 0; x < nix; x++) {
                                                        xsub = x*sub + m;
                                                        if ((xsub >= 0) && (xsub < dims[0])) {
                                                                label_value = (int)label[zsub2 + ysub2 + xsub];
                                                                label_value_BG = label_value - 1;
                                                                if (label_value_BG < 0) continue;
                                                                val = (double)src[zsub2 + ysub2 + xsub];
                                                                                
                                                                /* exclude values out of quartile 1-99% */
                                                                if ((val < thresh[0]) || (val > thresh[1])) continue;
                                                                ind = ((label_value_BG)*nvol)+yoffset+x;
                                                                ir[ind].n++;
                                                                ir[ind].s += val; ir[ind].ss += val*val;
                                                        }
                                                }
                                        }
                                }
                        }
                }


        /* find means and standard deviations */
        for(i = 0; i < n_classes; i++) {
                for(j = 0; j < nvol; j++) {
                        ind = (i*nvol)+j;
                        if (ir[ind].n > G) {
                                r[ind].mean = ir[ind].s/ir[ind].n;
                                if (ir[ind].n == 1)
                                        r[ind].var = 0.0;
                                else            r[ind].var = (ir[ind].ss -ir[ind].n*SQR(r[ind].mean))/(ir[ind].n-1);
                        } else r[ind].mean = 0.0;
                }
        }

        free(ir);

        return;
}

/* Computes likelihood of value given parameters mean and variance */ 
double ComputeGaussianLikelihood(double value, double mean , double var)

{ 
        return(exp(-(SQR(value - mean))/(2 * var))/SQRT2PI/sqrt(var));
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

double ComputeMarginalizedLikelihood(double value, double mean1 , double mean2, 
        double var1, double var2, unsigned int nof_intervals)

{ 
        double lh, tmean, tvar, delta, step;
        
        step = 1.0 / (double) nof_intervals;
        lh = 0;
        
        for(delta = 0.0; delta <= 1.0; delta += step) {
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
        
        for(i = 1; i < n; i++) {
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

        for(i = 0; i < n; i++) 
                sum_val += val[i];
 
        if (fabs(sum_val) > TINY) {      /* To avoid divisions by zero */
                for(i = 0; i < n; i++) {
                        val[i] /= sum_val;
                }
        }
}

/* Compute initial PVE labeling based on marginalized likelihood */
void ComputeInitialPveLabel(float *src, unsigned char *label, unsigned char *prob, struct point *r, int n_pure_classes, int sub, int *dims, int pve)
{
        int x, y, z, z_area, y_dims, index, label_value, off;
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
        
        /* use 5 or 6 classes */
        if (pve == 6) off = 1;
        else             off = 0;
        
        /* loop over image points */
        for(z = 1; z < dims[2]-1; z++) {
                z_area=z*area;
                for(y = 1; y < dims[1]-1; y++) {
                        y_dims=y*dims[0];
                        for(x = 1; x < dims[0]-1; x++)  {
                 
                                index = x + y_dims + z_area;
                                label_value = (int)label[index];
                                if (label_value == 0) continue;
                                val = (double)src[index];
                                        
                                /* find the interpolation factors */
                                ix = (int)(sub_1*x), iy = (int)(sub_1*y), iz = (int)(sub_1*z);
                                ind = iz*narea + iy*nix + ix;
                                        
                                for(i = 0; i < n_pure_classes; i++) {
                                        ind2 = (i*nvol) + ind;                                          
                                        if (r[ind2].mean > 0.0) {
                                                mean[off+i*2] = r[ind2].mean;
                                                var[off+i*2]    = r[ind2].var;
                                        }
                                }

                                if (fabs(mean[CSFLABEL+off-1]) > TINY) {
                                        d_pve[CSFLABEL+off-1] = ComputeGaussianLikelihood(val, mean[CSFLABEL+off-1], var[CSFLABEL+off-1]);
                                } else d_pve[CSFLABEL+off-1] = HUGE;

                                if (fabs(mean[GMLABEL+off-1]) > TINY) {
                                        d_pve[GMLABEL+off-1] = ComputeGaussianLikelihood(val, mean[GMLABEL+off-1], var[GMLABEL+off-1]);
                                } else d_pve[GMLABEL+off-1] = HUGE;

                                if (fabs(mean[WMLABEL+off-1]) > TINY) {
                                        d_pve[WMLABEL+off-1] = ComputeGaussianLikelihood(val, mean[WMLABEL+off-1], var[WMLABEL+off-1]);
                                } else d_pve[WMLABEL+off-1] = HUGE;

                                if ((fabs(mean[WMLABEL+off-1]) > TINY) && (fabs(mean[GMLABEL+off-1]) > TINY)) {
                                        d_pve[WMGMLABEL+off-1] = ComputeMarginalizedLikelihood(val, mean[WMLABEL+off-1], mean[GMLABEL+off-1],
                                                var[WMLABEL+off-1], var[GMLABEL+off-1], 100 );
                                } else d_pve[WMGMLABEL+off-1] = HUGE;
                                                
                                if ((fabs(mean[CSFLABEL+off-1]) > TINY) && (fabs(mean[GMLABEL+off-1]) > TINY)) {
                                        d_pve[GMCSFLABEL+off-1] = ComputeMarginalizedLikelihood(val, mean[GMLABEL+off-1], mean[CSFLABEL+off-1],
                                                var[GMLABEL+off-1], var[CSFLABEL+off-1], 100 );
                                } else d_pve[GMCSFLABEL+off-1] = HUGE;
                                
                                /* BKGCSF only for 6 classes */
                                if (pve == 6) {
                                        if (fabs(mean[CSFLABEL+off-1]) > TINY) {
                                                d_pve[BKGCSFLABEL+off-1] = ComputeMarginalizedLikelihood(val, 0.0, mean[CSFLABEL+off-1],
                                                        0.1*MIN3(var[CSFLABEL+off-1],var[GMLABEL+off-1],var[WMLABEL+off-1]), var[CSFLABEL+off-1], 100 );
                                        } else d_pve[BKGCSFLABEL+off-1] = HUGE;
                                }

                                Normalize(d_pve, n_pure_classes+2+off);
                                
                                for(i = 0; i < n_pure_classes+2+off; i++) 
                                        prob[(vol*i) + index] = (unsigned char)ROUND(255*d_pve[i]);

                                label[index] = (unsigned char) MaxArg(d_pve, n_pure_classes+2+off);
                        }
                }
        }               
} 

void ComputeMrfProbability(double *mrf_probability, double *exponent, unsigned char *label, int x, int y , int z, int *dims,
                          int n_classes, double beta, double *voxelsize_squared)
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
        for(label1 = 0; label1 < n_classes; label1++)
                exponent[label1] = 0;
        
        for(i = -1; i < 2; i++) for(j = -1; j < 2; j++) for(k = -1; k < 2; k++) 
                if ( i != 0 || j != 0 || k != 0 ) {
                                         
                        label2 = label[(x+i)+dims[0]*(y+j)+dims[0]*dims[1]*(z+k)];
                                                         
                        for(label1 = 1; label1 < n_classes+1; label1++) { 
                                if (label1 == label2) similarity_value = same;
                                else if (abs(label1 - label2) < 2) similarity_value = similar;
                                else similarity_value = different;

                                distance = sqrt(voxelsize_squared[0] * abs(i) + voxelsize_squared[1] * abs(j) + voxelsize_squared[2] * abs(k));

                                exponent[label1-1] += (double)similarity_value/distance;                                                 
                        }
                }               

        for(label1 = 0; label1 < n_classes; label1++)
                mrf_probability[label1] = exp(-(beta*exponent[label1])); 
        
} 

/* Iterative conditional mode */
void ICM(unsigned char *prob, unsigned char *label, int n_classes, int *dims, double beta, int iterations, double *voxelsize)
{
        
        int i, iter, x, y, z, z_area, y_dims, index, sum_voxel;
        int area, vol;
        double rel_changed, mrf_probability[MAX_NC], voxelsize_squared[3];
        double exponent[MAX_NC], sum_voxelsize = 0.0;
        unsigned char new_label;
                
        area = dims[0]*dims[1];
        vol = area*dims[2];
        
        /* normalize voxelsize to a sum of 3 and calculate its squared value */
        for(i = 0; i < 3; i++) sum_voxelsize += voxelsize[i];
        for(i = 0; i < 3; i++) voxelsize_squared[i] = SQR(3.0*voxelsize[i]/sum_voxelsize);
        
                
        for(iter=0; iter < iterations; iter++) {
                sum_voxel = 0;
                rel_changed = 0.0;
                
                /* loop over image points */
                for(z = 1; z < dims[2]-1; z++) {
                        z_area=z*area;
                        for(y = 1; y < dims[1]-1; y++) {
                                y_dims=y*dims[0];
                                for(x = 1; x < dims[0]-1; x++)  {
                 
                                        index = x + y_dims + z_area;
                                        if (label[index] == 0) continue;
                                        
                                        sum_voxel++;
                                        ComputeMrfProbability(mrf_probability, exponent, label, x, y, z, dims, n_classes, beta, voxelsize_squared);
                                        
                                        for(i = 0; i < n_classes; i++)
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
#if !defined(_WIN32)
                printf("ICM: %d relative change: %2.4f\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b",iter+1, 100.0*rel_changed);
                fflush(stdout);
#endif
                if (rel_changed < TH_CHANGE) break;
        }               
        printf("\n");
} 

void EstimateSegmentation(float *src, unsigned char *label, unsigned char *prob, struct point *r, double *mean, double *var, int n_classes, int niters, int sub, int *dims, double *thresh, double *beta, double offset)
{
        int i;
        int area, narea, nvol, vol, z_area, y_dims, index, ind;
        double sub_1, dmin, val;
        double d[MAX_NC], alpha[MAX_NC], log_alpha[MAX_NC], log_var[MAX_NC];
        double pvalue[MAX_NC], psum;
        int nix, niy, niz, iters, count_change;
        int x, y, z, label_value, xBG;
        int ix, iy, iz, ind2;
        double ll, ll_old, change_ll;

        MrfPrior(label, n_classes, alpha, beta, 0, dims);                

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

        for(i = 0; i < n_classes; i++) log_alpha[i] = log(alpha[i]);
                
        ll_old = HUGE;
        count_change = 0;
                        
        for(iters = 0; iters < niters; iters++)  {
                        
                ll = 0.0;
                
                /* get means for grid points */
                GetMeansVariances(src, label, n_classes, r, sub, dims, thresh);          

                /* loop over image points */
                for(z = 1; z < dims[2]-1; z++) {
                        z_area=z*area;
                        for(y = 1; y < dims[1]-1; y++) {
                                y_dims=y*dims[0];
                                for(x = 1; x < dims[0]-1; x++)  {
                 
                                        index = x + y_dims + z_area;
                                        label_value = (int) label[index];
                                        if (label_value < 1) continue;
                                        val = (double)src[index];
                                        
                                        /* find the interpolation factors */
                                        ix = (int)(sub_1*x);
                                        iy = (int)(sub_1*y);
                                        iz = (int)(sub_1*z);
                                        ind = iz*narea + iy*nix + ix;
                                        
                                        for(i = 0; i < n_classes; i++) {
                                                ind2 = (i*nvol) + ind;  
                                                if (r[ind2].mean > TINY) {
                                                        mean[i] = r[ind2].mean;
                                                        var[i]  = r[ind2].var;
                                                        log_var[i] = log(var[i]);
                                                }
                                        }
                                        
                                        /* compute energy at each point */
                                        dmin = HUGE; xBG = 1; 
                                        psum = 0.0;

                                        for(i = 0; i < n_classes; i++) {
                                                if (fabs(mean[i]) > TINY) {
                                                        d[i] = 0.5*(SQR(val-mean[i])/var[i]+log_var[i])-log_alpha[i];
                                                        pvalue[i] = exp(-d[i])/SQRT2PI;
                                                        psum += pvalue[i];
                                                } else d[i] = HUGE;
                                                if ( d[i] < dmin) {
                                                        dmin = d[i];
                                                        xBG = i;
                                                }
                                        }
                                                         
                                        /* scale p-values to a sum of 1 */
                                        if (psum > TINY) {
                                                for(i = 0; i < n_classes; i++) pvalue[i] /= psum;
                                                ll -= log(psum);
                                        } else  for(i = 0; i < n_classes; i++) pvalue[i] = 0.0;
                                 
                                        for(i = 0; i < n_classes; i++)
                                                prob[(vol*i) + index] = (unsigned char)ROUND(255*pvalue[i]);
                                 
                                        /* if the class has changed modify the label */
                                        if (xBG + 1 != label_value) label[index] = (unsigned char) (xBG + 1); 
                                 
                                }
                        }
                }

                ll /= (double)vol;
                change_ll = (ll_old - ll)/fabs(ll);
#if !defined(_WIN32)
                printf("iters:%3d log-likelihood: %7.5f\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b",iters+1, ll);
                fflush(stdout);
#endif
                ll_old = ll;
                
                /* break if log-likelihood has not changed significantly two iterations */
                if (change_ll < TH_CHANGE) count_change++;
                if (count_change > 1) break;            
        }

        printf("\nFinal Mean*Std: "); 
        for(i = 0; i < n_classes; i++) printf("%.3f*%.3f        ",mean[i]-offset,sqrt(var[i])); 
        printf("\n"); 

}


/* perform adaptive MAP on given src and initial segmentation label */
void Amap(float *src, unsigned char *label, unsigned char *prob, double *mean, int n_classes, int niters, int sub, int *dims, int pve, double weight_MRF, double *voxelsize, int niters_ICM, double offset)
{
        int i, nix, niy, niz;
        int area, nvol, vol;
        int histo[65536];
        int n[MAX_NC], j;
        double var[MAX_NC];
        double thresh[2], beta[1];
        double min_src = HUGE, max_src = -HUGE;
        int cumsum[65536];
        struct point *r;
        
        /* ICM is not needed if we skip MRF approach */
        if (weight_MRF == 0)
                niters_ICM = 0;
                        
        area = dims[0]*dims[1];
        vol = area*dims[2];
 
        for(i = 0; i < vol; i++) {
                min_src = MIN((double)src[i], min_src);
                max_src = MAX((double)src[i], max_src);
        }
        
        /* build histogram */
        for(i = 0; i < 65536; i++) histo[i] = 0;
        for(i = 0; i < vol; i++) {
                if (label[i] == 0) continue;
                histo[(int)ROUND(65535.0*((double)src[i]-min_src)/(max_src-min_src))]++;
        }

        /* find values between 1% and 99% quartile */
        cumsum[0] = histo[0];
        for(i = 1; i < 65536; i++) cumsum[i] = cumsum[i-1] + histo[i];
        for(i = 0; i < 65536; i++) cumsum[i] = (int) ROUND(1000.0*(double)cumsum[i]/(double)cumsum[65535]);
        for(i = 0; i < 65536; i++) if (cumsum[i] >= 10) break;
        thresh[0] = (double)i/65535.0*(max_src-min_src);
        for(i = 65535; i > 0; i--) if (cumsum[i] <= 990) break;
        thresh[1] = (double)i/65535.0*(max_src-min_src);
 
        /* define grid dimensions */
        nix = (int) ceil((dims[0]-1)/((double) sub))+1;
        niy = (int) ceil((dims[1]-1)/((double) sub))+1;
        niz = (int) ceil((dims[2]-1)/((double) sub))+1; 
        nvol    = nix*niy*niz;

        r = (struct point*)malloc(sizeof(struct point)*MAX_NC*nvol);
        if (r == NULL) {
                printf("Memory allocation error\n");
                exit(EXIT_FAILURE);
        }
                
        /* estimate 3 classes before PVE */
        EstimateSegmentation(src, label, prob, r, mean, var, n_classes, niters, sub, dims, thresh, beta, offset);
        
        /* Use marginalized likelihood to estimate initial 5 or 6 classes */
        if (pve) {

                ComputeInitialPveLabel(src, label, prob, r, n_classes, sub, dims, pve);
                n_classes = pve;
                
                /* recalculate means for pure and mixed classes */
                for(j = 0; j < n_classes; j++) {
                        n[j] = 0;
                        mean[j] = 0.0;
                }
                for(i = 0; i < vol; i++) {
                        if (label[i] == 0) continue;
                        n[label[i]-1]++;
                        mean[label[i]-1] += (double)src[i];
                }
                for(j = 0; j < n_classes; j++) mean[j] /= n[j];
        }
        
        /* use much smaller beta for if no pve is selected */
        if (!pve) beta[0] /= 20.0;
        
        /* Iterative Conditional Mode */
        if (niters_ICM > 0) {
                if (weight_MRF != 1.0) {
                        beta[0] *= weight_MRF;
                        printf("Weighted MRF beta %3.3f\n",beta[0]);
                }
        
                ICM(prob, label, n_classes, dims, beta[0], niters_ICM, voxelsize);
        }
        
        free(r);

        return;          
}



