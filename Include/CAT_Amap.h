/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 *
 */

#ifndef _CAT_AMAP_H_
#define _CAT_AMAP_H_

#define SQRT2PI 2.506628
#define G 6

#define TH_COLOR 1
#define TH_CHANGE 0.001
#ifndef TINY
#define TINY 1e-15 
#endif
#ifndef HUGE
#define HUGE 1e15 
#endif
#ifndef NULL
#define NULL ((void *) 0)
#endif

#define CSFLABEL    1
#define GMCSFLABEL  2
#define GMLABEL     3
#define WMGMLABEL   4
#define WMLABEL     5

#ifndef SQR
#define SQR(x) ((x)*(x))
#endif

#ifndef MAX
#define MAX(A,B) ((A) > (B) ? (A) : (B))
#endif

#ifndef MIN
#define MIN(A,B) ((A) < (B) ? (A) : (B))
#endif

#ifndef ROUND
#define ROUND( x ) ((long) ((x) + ( ((x) >= 0) ? 0.5 : (-0.5) ) ))
#endif

#include <math.h>

void Amap(float *src, unsigned char *label, unsigned char *prob, double *mean, 
                 int nc, int niters, int sub, int *dims, int pve, double weight_MRF, 
                 double *voxelsize, int niters_ICM, int verbose, 
                 int use_median);
void Pve5(float *src, unsigned char *prob, unsigned char *label, double *mean, int *dims);
double ComputeGaussianLikelihood(double value, double mean , double var);
double ComputeMarginalizedLikelihood(double value, double mean1 , double mean2, 
    double var1, double var2, unsigned int nof_intervals);
void MrfPrior(unsigned char *label, int n_classes, double *alpha, double *beta, int init, int *dims, int verbose);
void Normalize(double* val, char n);
unsigned char MaxArg(double *val, unsigned char n);

struct point {
  double median;
  double mean;
  double var;
};

struct ipoint {
  int n;
  double s;
  double ss;
  double *arr;
};

#endif