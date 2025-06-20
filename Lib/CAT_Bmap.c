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
#include "CAT_Bmap.h"
#include "CAT_Vol.h"
#include "CAT_Math.h"

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

void xaverage(unsigned char *label, float *bias, long *n1, double *bs, int bg, int a, float thresh, int *dims)
{
  int x,y,z,r;

  int area = dims[0]*dims[1];
  for (z = 0; z < dims[2]; z++)  {
    int i = z*area;
    for (y = 0; y < dims[1]; y++)  
      for (x = 0; x < dims[0]; x++) {
        i++; n1[i] = 0;  bs[i] = 0;
        int l = x-a < 0 ? 0 : x-a;
        int h = x+a > dims[0]-1 ? dims[0]-1 : x+a;
        for (r=l; r <=h  ; r++) {
          int j = i + (r-x);
          if ((int) label[j]>=bg && fabs(bias[j])<thresh) {
            n1[i]++;
            bs[i] += (double)bias[j];
          }
        }
     }
  }
}

void yaverage(long *n1, double *bs, long *n2, float *bias, int b, int *dims)
{
  int x,y,z,r;

  int area = dims[0]*dims[1];

  for (z=0; z<dims[2]; z++)  {
    int i = z*area;
    for (y = 0; y < dims[1]; y++)  
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

void zaverage(long *n2, float *bias, long *n1, double *bs, int bg, int c, int *dims)
{
  int x,y,z,r;

  int area = dims[0]*dims[1];

  for (z = 0; z < dims[2]; z++)  {
    int i = z*area;
    for (y = 0; y < dims[1]; y++)  
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

void movingAverage(unsigned char *label, float *bias, int BG, int a, int b, int c, float thresh, int *dims)
{
  int i, j, x, y, z, r, l, h;
  double *bs;
  long *n1, *n2;

  long vol = dims[0]*dims[1]*dims[2];
  
  bs = (double *) malloc(sizeof(double)*vol);
  n1 = (long *) malloc(sizeof(long)*vol);
  n2 = (long *) malloc(sizeof(long)*vol);
  
  xaverage(label, bias, n1, bs, BG, a, thresh, dims);
  yaverage(n1, bs, n2, bias, b, dims);
  zaverage(n2, bias, n1, bs, BG, c, dims);
    
  for (i=0; i<vol; i++) {
    if ((int) label[i] >= BG) 
      bias[i] =  n1[i] > 0 ? bs[i]/n1[i] : 0.0;
    else 
      bias[i] = 0.0;
  }
  
  free(bs);
  free(n1);
  free(n2);

  return;
}

void Bmap(float *src, unsigned char *label, unsigned char *prob, int n_classes, int BG, int niters, int nflips, int a, int b, int c, float *bias, float bias_thresh, int *dims, int pve)
{
  int i,j,index,x,y,z,iters;
  int n[MAX_NC];
  double thresh[2], val, d;
  double beta[1],alpha[n_classes],s[n_classes],ss[n_classes];
  double mean[n_classes],var[n_classes],lvar[n_classes],p[n_classes];
  
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

  for (i=0; i<vol; i++) bias[i]=0.0;
  
  // initialize prior parameters
  MrfPrior(label, n_classes, alpha, beta, 0, dims, 0);

  for (i = 0; i < n_classes; i++) fprintf(stderr,"%.3f*%.3f    ",mean[i],sqrt(var[i])); 
  fprintf(stderr,"\n"); 

  // set new variables to speed up
  for (j=0; j<n_classes; j++) {
    lvar[j] = var[j] > 0.0 ? 0.5*log(var[j]) - log(alpha[j]): - log(alpha[j]);
    var[j] = 0.5/var[j];
  }

  int do_mrf = 0;
  int do_bias = 0;
  
  // iterative condition modes
  for (iters=0; iters<=niters; iters++)  {
    int flips = 0;

    // loop over image voxels
    for (z=1; z<dims[2]-1; z++) 
      for (y=1; y<dims[1]-1; y++)
        for (x=1; x<dims[0]-1; x++)  {
    
          index = x + y*dims[0] + z*area;

          int lab = (int) label[index];
          if (lab < BG) continue;
    
          // loop over all classes
          int xi=BG; double dmin = 1e15;
          double psum = 0.0;
          
          for (j=0; j<n_classes; j++) {

            // find the number of first order neighbors in the class
            int first=0;
            if (do_mrf) {
              if (label[index-1] == j+BG) first++;
              if (label[index+1] == j+BG) first++;
              if (label[index-dims[0]] == j+BG) first++;
              if (label[index+dims[0]] == j+BG) first++;
              if (label[index-area] == j+BG) first++;
              if (label[index+area] == j+BG) first++;
            } else first = 5;
      
            d = SQR((double)src[index]-(1.0+(double)bias[index])*mean[j])*var[j]+lvar[j]-beta[0]*first;
      
            p[j] = exp(-d)/(SQRT2PI*sqrt(var[j]));
            psum += p[j];
      
            if (d < dmin) {xi = j; dmin = d;}
          }
          for (j=0; j<n_classes; j++)
            prob[(vol*j) + index] = (unsigned char)round(255*p[j]/psum);

          // if the class has changed increment flips and change the label
          if (xi+BG != lab) {flips++;  label[index] = (unsigned char) (xi+BG); }
    
          // find bias values
          bias[index] = src[index]/mean[xi] - 1.0;
        }
        
    fprintf(stderr,"iters:%2d flips:%6d\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b",iters, flips);

    if ((flips <= nflips) && (do_mrf)) do_mrf = 0;
    if ((flips <= nflips) && (iters > 10)) break;

    // smoothout bias values (remove misclassifications)
    if ((do_bias) && (iters<=1))
      movingAverage(label, bias, BG, a, b, c, bias_thresh, dims);
    
    for (j = 0; j < n_classes; j++) { s[j] = ss[j] = 0.0; }
    for (i = 0; i < vol; i++) {
        int lab = (int) label[i];
        if (lab < BG) continue;
        double val = (1.0 + bias[i]) * (double)src[i];
    
        s[lab - BG] += 1.0;
        ss[lab - BG] += val;
    }
    for (j = 0; j < n_classes; j++) {
        mean[j] = (s[j] > 0.0) ? (ss[j] / s[j]) : 0.0;
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

  }
    
  for (j = 0; j < n_classes; j++) {
      var[j] = 0.5/var[j];
  }

  fprintf(stderr,"\nFinal means*vars: "); 
  for (i = 0; i < n_classes; i++) fprintf(stderr,"%.3f*%.3f    ",mean[i],sqrt(var[i])); 
  fprintf(stderr,"\n"); 

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
  for (z=0; z<dims[2]; z++) 
    for (y=0; y<dims[1]; y++) 
      for (x=0; x<dims[0]; x++, i++)  {
      int lab = (int) label[i];
      if (lab < BG) continue;
      bias[i] = (double)src[i] -(1.0+bias[i])*mean[lab-BG];
    }

  return;
}  
