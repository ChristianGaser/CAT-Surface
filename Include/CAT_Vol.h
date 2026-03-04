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

#define index(A, B, C, DIM) ((C) * DIM[0] * DIM[1] + (B) * DIM[0] + (A))

#define CSF 1.0
#define CGM 1.5
#define GM 2.0
#define GWM 2.5
#define WM 3.0

#define MAX_NC 6

/* --------------------------- Thread args --------------------------- */
typedef struct
{
    float *out;
    int xdim, ydim;
    const double *filtx, *filty;
    int fxdim, fydim;
    int xoff, yoff;
    int ini, fin; /* range on outer index: rows for row-pass, cols for col-pass */
} conv_args_row;

typedef struct
{
    float *out;
    int xdim, ydim;
    const double *filtx, *filty;
    int fxdim, fydim;
    int xoff, yoff;
    int ini, fin; /* columns range */
} conv_args_col;

typedef struct
{
    /* inputs */
    const float *iVol;
    int xdim, ydim, zdim;
    const double *filtx, *filty;
    int fxdim, fydim;
    int xoff, yoff;
    /* outputs */
    float *convxy_vol; /* [zdim * xdim * ydim] */
    /* range */
    int z_ini, z_fin; /* [z_ini, z_fin) */
} convxyz_s1_args_t;

typedef struct
{
    const float *convxy_vol; /* [zdim * xdim * ydim] */
    float *oVol;             /* [zdim * xdim * ydim] */
    int xdim, ydim, zdim;
    const double *filtz;
    int fzdim;
    int zoff;
    int z_out_ini, z_out_fin; /* [z_out_ini, z_out_fin) */
} convxyz_s2_args_t;

/**
 * \brief Public API for median3.
 *
 * This function is part of the CAT-Surface public library interface and is used by command-line tools.
 *
 * \param D (in/out) Parameter of median3.
 * \param mask (in/out) Parameter of median3.
 * \param dims (in/out) Parameter of median3.
 * \param iters (in/out) Parameter of median3.
 * \param datatype (in/out) Parameter of median3.
 * \return void (no return value).
 */
void median3(void *D, unsigned char *mask, int dims[3], int iters, int datatype);
/**
 * \brief Public API for localstat3.
 *
 * This function is part of the CAT-Surface public library interface and is used by command-line tools.
 *
 * \param input (in/out) Parameter of localstat3.
 * \param mask (in/out) Parameter of localstat3.
 * \param dims (in/out) Parameter of localstat3.
 * \param dist (in/out) Parameter of localstat3.
 * \param stat_func (in/out) Parameter of localstat3.
 * \param iters (in/out) Parameter of localstat3.
 * \param use_euclidean_dist (in/out) Parameter of localstat3.
 * \param datatype (in/out) Parameter of localstat3.
 * \return void (no return value).
 */
void localstat3(void *input, unsigned char mask[], int dims[3], int dist, int stat_func, int iters, int use_euclidean_dist, int datatype);
void laplace3R(float *SEG, unsigned char *M, int dims[3], double TH);
/**
 * \brief Public API for smooth3.
 *
 * This function is part of the CAT-Surface public library interface and is used by command-line tools.
 *
 * \param vol (in/out) Parameter of smooth3.
 * \param dims (in/out) Parameter of smooth3.
 * \param voxelsize (in/out) Parameter of smooth3.
 * \param s (in/out) Parameter of smooth3.
 * \param use_mask (in/out) Parameter of smooth3.
 * \param datatype (in/out) Parameter of smooth3.
 * \return void (no return value).
 */
void smooth3(void *vol, int dims[3], double voxelsize[3], double s[3], int use_mask, int datatype);
void smooth_subsample3(void *vol, int dims[3], double voxelsize[3], double s[3], int use_mask, double samp_voxelsize, int datatype);
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
void median_subsample3(void *data, int dims[3], double voxelsize[3], int niter, double samp_voxelsize, int datatype);
/**
 * \brief Public API for isoval.
 *
 * This function is part of the CAT-Surface public library interface and is used by command-line tools.
 *
 * \param vol (in/out) Parameter of isoval.
 * \param x (in/out) Parameter of isoval.
 * \param y (in/out) Parameter of isoval.
 * \param z (in/out) Parameter of isoval.
 * \param s (in/out) Parameter of isoval.
 * \param nii_ptr (in/out) Parameter of isoval.
 * \return Return value of isoval.
 */
float isoval(float vol[], float x, float y, float z, int s[], nifti_image *nii_ptr);
/**
 * \brief Public API for correct_bias.
 *
 * This function is part of the CAT-Surface public library interface and is used by command-line tools.
 *
 * \param src (in/out) Parameter of correct_bias.
 * \param biasfield (in/out) Parameter of correct_bias.
 * \param label (in/out) Parameter of correct_bias.
 * \param dims (in/out) Parameter of correct_bias.
 * \param voxelsize (in/out) Parameter of correct_bias.
 * \param bias_fwhm (in/out) Parameter of correct_bias.
 * \param weight_las (in/out) Parameter of correct_bias.
 * \return void (no return value).
 */
void correct_bias(float *src, float *biasfield, unsigned char *label, int *dims, double *voxelsize, double bias_fwhm, double weight_las);
void morph_erode(void *vol, int dims[3], int niter, double th, int datatype);
/**
 * \brief Public API for morph_dilate.
 *
 * This function is part of the CAT-Surface public library interface and is used by command-line tools.
 *
 * \param vol (in/out) Parameter of morph_dilate.
 * \param dims (in/out) Parameter of morph_dilate.
 * \param niter (in/out) Parameter of morph_dilate.
 * \param th (in/out) Parameter of morph_dilate.
 * \param datatype (in/out) Parameter of morph_dilate.
 * \return void (no return value).
 */
void morph_dilate(void *vol, int dims[3], int niter, double th, int datatype);
/**
 * \brief Public API for morph_close.
 *
 * This function is part of the CAT-Surface public library interface and is used by command-line tools.
 *
 * \param vol (in/out) Parameter of morph_close.
 * \param dims (in/out) Parameter of morph_close.
 * \param niter (in/out) Parameter of morph_close.
 * \param th (in/out) Parameter of morph_close.
 * \param datatype (in/out) Parameter of morph_close.
 * \return void (no return value).
 */
void morph_close(void *vol, int dims[3], int niter, double th, int datatype);
/**
 * \brief Public API for morph_open.
 *
 * This function is part of the CAT-Surface public library interface and is used by command-line tools.
 *
 * \param vol (in/out) Parameter of morph_open.
 * \param dims (in/out) Parameter of morph_open.
 * \param niter (in/out) Parameter of morph_open.
 * \param th (in/out) Parameter of morph_open.
 * \param keep_values (in/out) Parameter of morph_open.
 * \param datatype (in/out) Parameter of morph_open.
 * \return void (no return value).
 */
void morph_open(void *vol, int dims[3], int niter, double th, int keep_values, int datatype);
void grey_erode(void *data, int dims[3], int niter, int datatype);
void grey_dilate(void *data, int dims[3], int niter, int datatype);
void grey_open(void *data, int dims[3], int niter, int datatype);
void grey_close(void *data, int dims[3], int niter, int datatype);
void dist_close(void *vol, int dims[3], double voxelsize[3], double dist, double th, int datatype);
void dist_close_float(float *vol, int dims[3], double voxelsize[3], double dist, double th, unsigned char *mask);
void dist_open(void *vol, int dims[3], double voxelsize[3], double dist, double th, int datatype);
void dist_open_float(float *vol, int dims[3], double voxelsize[3], double dist, double th, unsigned char *mask);
void dist_erode(void *vol, int dims[3], double voxelsize[3], double dist, double th, int datatype);
void dist_erode_float(float *vol, int dims[3], double voxelsize[3], double dist, double th, unsigned char *mask);
void dist_dilate(void *vol, int dims[3], double voxelsize[3], double dist, double th, int datatype);
void dist_dilate_float(float *vol, int dims[3], double voxelsize[3], double dist, double th, unsigned char *mask);
/**
 * \brief Public API for subsample3.
 *
 * This function is part of the CAT-Surface public library interface and is used by command-line tools.
 *
 * \param in (in/out) Parameter of subsample3.
 * \param out (in/out) Parameter of subsample3.
 * \param dims (in/out) Parameter of subsample3.
 * \param dims_samp (in/out) Parameter of subsample3.
 * \param datatype (in/out) Parameter of subsample3.
 * \return void (no return value).
 */
void subsample3(void *in, void *out, int dims[3], int dims_samp[3], int datatype);
void vol_approx(float *vol, int dims[3], double voxelsize[3]);
/**
 * \brief Public API for cleanup_brain.
 *
 * This function is part of the CAT-Surface public library interface and is used by command-line tools.
 *
 * \param probs (in/out) Parameter of cleanup_brain.
 * \param dims (in/out) Parameter of cleanup_brain.
 * \param voxelsize (in/out) Parameter of cleanup_brain.
 * \param strength (in/out) Parameter of cleanup_brain.
 * \return void (no return value).
 */
void cleanup_brain(unsigned char *probs, int dims[3], double voxelsize[3], int strength);
/**
 * \brief Public API for euclidean_distance.
 *
 * This function is part of the CAT-Surface public library interface and is used by command-line tools.
 *
 * \param V (in/out) Parameter of euclidean_distance.
 * \param IO (in/out) Parameter of euclidean_distance.
 * \param dims (in/out) Parameter of euclidean_distance.
 * \param voxelsize (in/out) Parameter of euclidean_distance.
 * \param replace (in/out) Parameter of euclidean_distance.
 * \return void (no return value).
 */
void euclidean_distance(float *V, unsigned char *IO, int dims[3], double *voxelsize, int replace);
/**
 * \brief Intensity-limited region growing with distance/intensity path cost.
 *
 * Grows labels from seeded voxels into unlabeled voxels under a monotonic
 * intensity constraint and weighted path cost. This is the float core function.
 *
 * \param labels     (in/out) seed label map; 0 means unlabeled
 * \param intensity  (in)     intensity image controlling growth
 * \param dist       (out)    path-cost map (NULL allowed)
 * \param dims       (in)     dimensions {nx, ny, nz}
 * \param limit      (in)     neighbour intensity limit
 * \param voxelsize  (in)     voxel spacing {sx, sy, sz}; NULL -> {1,1,1}
 * \param dd         (in)     weights {distance_weight, intensity_weight}; NULL -> defaults
 */
void downcut_float(float *labels, const float *intensity, float *dist,
                   int dims[3], double limit, double voxelsize[3], double dd[2]);
/**
 * \brief Datatype-generic wrapper for downcut region growing.
 *
 * Converts labels/intensity with convert_input_type(), runs downcut_float(), then
 * converts output buffers with convert_output_type().
 *
 * \param labels             (in/out) labels buffer in labels_datatype
 * \param intensity          (in)     intensity buffer in intensity_datatype
 * \param dist               (out)    distance buffer in dist_datatype (NULL allowed)
 * \param dims               (in)     dimensions {nx, ny, nz}
 * \param limit              (in)     neighbour intensity limit
 * \param voxelsize          (in)     voxel spacing {sx, sy, sz}; NULL -> {1,1,1}
 * \param dd                 (in)     weights {distance_weight, intensity_weight}; NULL -> defaults
 * \param labels_datatype    (in)     datatype code for labels input/output
 * \param intensity_datatype (in)     datatype code for intensity input
 * \param dist_datatype      (in)     datatype code for distance output
 */
void downcut3(void *labels, void *intensity, void *dist,
              int dims[3], double limit, double voxelsize[3], double dd[2],
              int labels_datatype, int intensity_datatype, int dist_datatype);
/**
 * \brief Blood-vessel correction for float PVE label maps.
 *
 * Implements blood-vessel correction on PVE labels in [0..3] using downcut growing,
 * ring-constrained median inpainting, and class-aware value clamping.
 *
 * \param Yp0              (in/out) float PVE label map
 * \param dims             (in)     dimensions {nx, ny, nz}
 * \param vx_vol           (in)     voxel spacing {sx, sy, sz}; NULL -> {1,1,1}
 */
void blood_vessel_correction_pve_float(float *Yp0, int dims[3], double vx_vol[3]);
/**
 * \brief Datatype-generic wrapper for PVE blood-vessel correction.
 *
 * Converts input volume to float, applies `blood_vessel_correction_pve_float()`, then converts back
 * to the specified datatype.
 *
 * \param data             (in/out) PVE label volume
 * \param dims             (in)     dimensions {nx, ny, nz}
 * \param vx_vol           (in)     voxel spacing {sx, sy, sz}; NULL -> {1,1,1}
 * \param datatype         (in)     datatype code (DT_UINT8, DT_UINT16, DT_FLOAT32, etc.)
 */
void blood_vessel_correction_pve(void *data, int dims[3], double vx_vol[3], int datatype);
void ind2sub(int i, int *x, int *y, int *z, int sxy, int sy);
int sub2ind(int x, int y, int z, int s[]);
void keep_largest_cluster(void *inData, double thresh, int *dims, int datatype, int min_size, int retain_above_th, int conn);
/**
 * \brief Public API for fill_holes.
 *
 * This function is part of the CAT-Surface public library interface and is used by command-line tools.
 *
 * \param data (in/out) Parameter of fill_holes.
 * \param dims (in/out) Parameter of fill_holes.
 * \param thresh (in/out) Parameter of fill_holes.
 * \param fill_value (in/out) Parameter of fill_holes.
 * \param datatype (in/out) Parameter of fill_holes.
 * \return void (no return value).
 */
void fill_holes(void *data, int *dims, double thresh, double fill_value, int datatype);
void gradient3D(float *src, float *grad_mag, float *grad_x, float *grad_y, float *grad_z, int dims[3], double voxelsize[3]);
#endif
