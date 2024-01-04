/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 *
 */

#ifndef _CAT_VOL_H_
#define _CAT_VOL_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <memory.h>
#include <float.h>
#include <limits.h>
#include "CAT_NiftiLib.h"

#define SQRT2PI 2.506628

#define RINT(A) floor((A)+0.5)
#ifndef isfinite
#define isfinite(x) ((x) * (x) >= 0.) /* check for NaNs */
#endif

#ifndef TINY
#define TINY 1e-15 
#endif

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

#ifndef isnan
#define isnan(a) ((a)!=(a)) 
#endif

#ifdef _MSC_VER
  static const unsigned long __nan[2] = {0xffffffff, 0x7fffffff};
  #define FNAN (*(const float *) __nan)
#else
  #define FNAN 0.0f/0.0f
#endif

#define index(A,B,C,DIM) ((C)*DIM[0]*DIM[1] + (B)*DIM[0] + (A))

#define CSF    1.0
#define CGM    1.5
#define GM     2.0
#define GWM    2.5
#define WM     3.0

#define MAX_NC 6

void correct_bias(float *src, unsigned char *label, int *dims, double *voxelsize, double bias_fwhm, int do_las);
double get_masked_mean_array_float(float arr[], int n, unsigned char mask[]);
double get_masked_std_array_float(float arr[], int n, unsigned char mask[]);
void correct_outer_rim(float *src, unsigned char *label, int *dims, double *voxelsize, int erosion_steps);
void get_prctile(float *src, int dims[3], double threshold[2], double prctile[2], int exclude_zeros);
void morph_erode(void *vol, int dims[3], int niter, double th, int datatype);
void morph_dilate(void *vol, int dims[3], int niter, double th, int datatype);
void morph_close(void *vol, int dims[3], int niter, double th, int datatype);
void morph_open(void *vol, int dims[3], int niter, double th, int datatype);
void distclose(void *vol, int dims[3], double voxelsize[3], int niter, double th, int datatype);
void distopen(void *vol, int dims[3], double voxelsize[3], int niter, double th, int datatype);
void subsample_double(double *in, double *out, int dim_in[3], int dim_out[3], int offset_in, int offset_out);
void subsample_float(float *in, float *out, int dim_in[3], int dim_out[3], int offset_in, int offset_out);
void subsample_uint8(unsigned char *in, float *out, int dim_in[3], int dim_out[3], int offset_in, int offset_out);
void vol_approx(float *vol, int dims[3], double voxelsize[3], int samp);
void smooth_double(double *vol, int dims[3], double voxelsize[3], double s[3], int use_mask);
void smooth_float(float *vol, int dims[3], double voxelsize[3], double s[3], int use_mask);
void smooth_subsample_double(double *vol, int dims[3], double voxelsize[3], double s[3], int use_mask, int samp);
void smooth_subsample_float(float *vol, int dims[3], double voxelsize[3], double s[3], int use_mask, int samp);
void initial_cleanup(unsigned char *probs, unsigned char *mask, int dims[3], double voxelsize[3], int strength, int remove_sinus);
void cleanup(unsigned char *probs, unsigned char *mask, int dims[3], double voxelsize[3], int strength, int gmwm_only);
void median3(void *D, int dims[3], int datatype);
void laplace3(float *SEG, int dims[3], int maxiter);
void laplace3R(float *SEG, unsigned char *M, int dims[3], double TH);
void vbdist(float *V, unsigned char *IO, int dims[3], double *voxelsize, int replace);
void ind2sub(int i,int *x,int *y, int *z, int sxy, int sy);
void projection_based_thickness(float *SEG, float *WMD, float *CSFD, float *GMT, int dims[3], double *voxelsize);
void get_largest_cluster(float *inData, double thresh, const int *dims);
float isoval(float vol[], float x, float y, float z, int s[], nifti_image *nii_ptr);

#endif