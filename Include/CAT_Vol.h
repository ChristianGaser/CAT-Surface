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
#include <math.h>
#include <memory.h>
#include <float.h>
#include <limits.h>
#include "CAT_NiftiLib.h"
#include "CAT_Math.h"

#define index(A,B,C,DIM) ((C)*DIM[0]*DIM[1] + (B)*DIM[0] + (A))

#define CSF    1.0
#define CGM    1.5
#define GM     2.0
#define GWM    2.5
#define WM     3.0

#define MAX_NC 6

enum
{
    F_MEAN,
    F_MIN,
    F_MAX,
    F_STD,
    F_SUM,
    F_MAXABS,
    F_EXP,
    F_MEDIAN,
    F_RANGE,
    F_COUNT,
    F_WAVERAGE,
    F_MULTI
};

void median3(void *D, unsigned char *mask, int dims[3], int iters, int datatype);
void localstat3(void *input, unsigned char mask[], int dims[3], int dist, int stat_func, int iters, int use_euclidean_dist, int datatype);
void laplace3R(float *SEG, unsigned char *M, int dims[3], double TH);
void smooth3(void *vol, int dims[3], double voxelsize[3], double s[3], int use_mask, int datatype);
void smooth_subsample3(void *vol, int dims[3], double voxelsize[3], double s[3], int use_mask, double samp_voxelsize, int datatype);
float isoval(float vol[], float x, float y, float z, int s[], nifti_image *nii_ptr);
void correct_bias(float *src, float *biasfield, unsigned char *label, int *dims, double *voxelsize, double bias_fwhm, double weight_las, int square_image);
void morph_erode(void *vol, int dims[3], int niter, double th, int datatype);
void morph_dilate(void *vol, int dims[3], int niter, double th, int datatype);
void morph_close(void *vol, int dims[3], int niter, double th, int datatype);
void morph_open(void *vol, int dims[3], int niter, double th, int datatype);
void grey_erode(void *data, int dims[3], int niter, int datatype);
void grey_dilate(void *data, int dims[3], int niter, int datatype);
void grey_open(void *data, int dims[3], int niter, int datatype);
void grey_close(void *data, int dims[3], int niter, int datatype);
void distclose(void *vol, int dims[3], double voxelsize[3], double dist, double th, int datatype);
void distopen(void *vol, int dims[3], double voxelsize[3], double dist, double th, int datatype);
void subsample3(void *in, void *out, int dims[3], int dims_samp[3], int datatype);
void vol_approx(float *vol, int dims[3], double voxelsize[3]);
void cleanup_brain(unsigned char *probs, int dims[3], double voxelsize[3], int strength);
void euclidean_distance(float *V, unsigned char *IO, int dims[3], double *voxelsize, int replace);
void ind2sub(int i,int *x,int *y, int *z, int sxy, int sy);
int sub2ind(int x, int y, int z, int s[]);
void projection_based_thickness(float *SEG, float *WMD, float *CSFD, float *GMT, int dims[3], double *voxelsize);
void keep_largest_cluster(void *inData, double thresh, int *dims, int datatype, int min_size, int retain_above_th, int conn);
void fill_holes(void *data, int *dims, double thresh, double fill_value, int datatype);
void gradient3D(float *src, float *grad_mag, float *grad_x, float *grad_y, float *grad_z, int dims[3], double voxelsize[3]);

#endif