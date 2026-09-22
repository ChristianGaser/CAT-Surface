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
 * \brief Apply median filtering to a 3D volume.
 *
 * \param data      (in/out) void pointer to volume data; type given by datatype parameter
 * \param mask      (in)     optional unsigned char mask (NULL to process entire volume)
 * \param dims      (in)     {nx, ny, nz} dimension array
 * \param iters     (in)     number of median filtering iterations
 * \param datatype  (in)     data type code (DT_UINT8, DT_UINT16, DT_FLOAT32, etc.)
 */
void median3(void *data, unsigned char *mask, int dims[3], int iters,
             int datatype);
/**
 * \brief Local statistic over a voxel neighbourhood, for any NIfTI datatype.
 *
 * \param data               (in/out) volume of type datatype, replaced by the result
 * \param mask               (in)     optional mask; voxels with 0 are skipped (NULL: all)
 * \param dims               (in)     volume dimensions {nx, ny, nz}
 * \param dist               (in)     search distance from the voxel centre, 1..10 voxels
 * \param stat_func          (in)     statistic: F_MEAN, F_MIN, F_MAX, F_STD, F_MEDIAN, ...
 *                                       (0=mean, 1=min, 2=max, 3=std, 7=median, 12=close,
 *                                       13=open; see CAT_VolLocalStat)
 * \param iters              (in)     number of iterations
 * \param use_euclidean_dist (in)     non-zero for a Euclidean, zero for a block neighbourhood
 * \param datatype           (in)     NIfTI datatype code of data (DT_FLOAT32, ...)
 */
void localstat3(void *data, unsigned char *mask, int dims[3], int dist,
                int stat_func, int iters, int use_euclidean_dist, int datatype);
/**
 * \brief Apply Laplace filter on a 3D volume.
 *
 * \param SEG 3D single input matrix (volume to be filtered).
 * \param M 3D volume that defines the filter area (mask).
 * \param dims Array containing the dimensions of the volume.
 * \param TH Threshold controlling the number of iterations (maximum change allowed after an iteration).
 */
void laplace3R(float *SEG, unsigned char *M, int dims[3], double TH);
/**
 * \brief Smooth a 3D volume with a Gaussian filter.
 *
 * \param data       (in/out) void pointer to volume data; type given by datatype parameter
 * \param dims       (in)     {nx, ny, nz} dimension array
 * \param voxelsize  (in)     voxel spacing in mm; used to scale FWHM to physical units
 * \param fwhm       (in)     {fwhm_x, fwhm_y, fwhm_z} smoothing kernel FWHM in mm
 * \param use_mask   (in)     unused/reserved for compatibility (pass 0)
 * \param datatype   (in)     data type code (DT_UINT8, DT_UINT16, DT_FLOAT32, etc.)
 */
void smooth3(void *data, int dims[3], double voxelsize[3], double fwhm[3],
             int use_mask, int datatype);
/**
 * \brief Gaussian smoothing on a subsampled grid, for any NIfTI datatype.
 *
 * \param data           (in/out) volume of type datatype, smoothed in-place
 * \param dims           (in)     volume dimensions {nx, ny, nz}
 * \param voxelsize      (in)     voxel size in mm
 * \param s              (in)     FWHM in mm per axis
 * \param use_mask       (in)     non-zero for masked smoothing (zeros are excluded)
 * \param samp_voxelsize (in)     voxel size in mm of the subsampled grid
 * \param datatype       (in)     NIfTI datatype code of data (DT_FLOAT32, ...)
 */
void smooth_subsample3(void *data, int dims[3], double voxelsize[3],
                       double s[3], int use_mask, double samp_voxelsize,
                       int datatype);
/**
 * \brief Median filter on a subsampled grid, for any NIfTI datatype.
 *
 * \param data           (in/out) volume of type datatype, filtered in-place
 * \param dims           (in)     volume dimensions {nx, ny, nz}
 * \param voxelsize      (in)     voxel size in mm
 * \param niter          (in)     number of median iterations
 * \param samp_voxelsize (in)     voxel size in mm of the subsampled grid
 * \param datatype       (in)     NIfTI datatype code of data (DT_FLOAT32, ...)
 */
void median_subsample3(void *data, int dims[3], double voxelsize[3], int niter,
                       double samp_voxelsize, int datatype);
/**
 * \brief Local statistic computed on a subsampled grid, for any NIfTI datatype.
 *
 * \param data               (in/out) volume of type datatype, replaced by the result
 * \param dims               (in)     volume dimensions {nx, ny, nz}
 * \param voxelsize          (in)     voxel size in mm
 * \param dist               (in)     search distance on the subsampled grid, 1..10 voxels
 * \param stat_func          (in)     statistic, as in localstat3()
 * \param niter              (in)     number of iterations
 * \param use_euclidean_dist (in)     non-zero for a Euclidean, zero for a block neighbourhood
 * \param samp_voxelsize     (in)     voxel size in mm of the subsampled grid
 * \param datatype           (in)     NIfTI datatype code of data (DT_FLOAT32, ...)
 */
void localstat_subsample3(void *data, int dims[3], double voxelsize[3],
                          int dist, int stat_func, int niter,
                          int use_euclidean_dist, double samp_voxelsize,
                          int datatype);
/**
 * \brief Trilinearly interpolated volume value at a point.
 *
 * \param vol     (in) float volume
 * \param x       (in) x coordinate
 * \param y       (in) y coordinate
 * \param z       (in) z coordinate
 * \param dims    (in) volume dimensions {nx, ny, nz}
 * \param nii_ptr (in) header providing the world-to-voxel mapping, or NULL
 * \return interpolated value, or NaN if no neighbour is finite
 */
float isoval(float *vol, float x, float y, float z, int dims[3],
             nifti_image *nii_ptr);
/**
 * \brief Adaptive bias correction for MRI images with optional subcortical refinement.
 *
 * \param src       (in/out) float[nvox]; source image, modified in-place with bias correction
 * \param biasfield (out)    float[nvox]; estimated bias field (can be NULL)
 * \param label     (in)     unsigned char[nvox]; tissue label map (CSF=1, GM=2, WM=3, etc.)
 * \param dims      (in)     {nx, ny, nz} volume dimensions
 * \param voxelsize (in)     {sx, sy, sz} voxel spacing in mm
 * \param bias_fwhm (in)     FWHM of Gaussian smoothing kernel for WM correction (mm)
 * \param weight_las (in)    weight for local adaptive segmentation GM correction (0..1);
 *                                0 = WM only, >0 = blend WM and GM with distance weighting
 */
void correct_bias(float *src, float *biasfield, unsigned char *label, int *dims,
                  double *voxelsize, double bias_fwhm, double weight_las);
/**
 * \brief Wrapper for binary morphological erosion (generic datatype).
 *
 * \param data       (in/out) void pointer to image data; type given by datatype parameter
 * \param dims       (in)     {nx, ny, nz}
 * \param niter      (in)     number of erosion iterations (<=0: no-op)
 * \param th         (in)     threshold as fraction of max(data) in [0,1]
 * \param datatype   (in)     data type code (DT_UINT8, DT_UINT16, DT_FLOAT32, etc.)
 */
void morph_erode(void *data, int dims[3], int niter, double th, int datatype);
/**
 * \brief Wrapper for binary morphological dilation (generic datatype).
 *
 * \param data       (in/out) void pointer to image data; type given by datatype parameter
 * \param dims       (in)     {nx, ny, nz}
 * \param niter      (in)     number of dilation iterations (<=0: no-op)
 * \param th         (in)     threshold as fraction of max(data) in [0,1]
 * \param datatype   (in)     data type code (DT_UINT8, DT_UINT16, DT_FLOAT32, etc.)
 */
void morph_dilate(void *data, int dims[3], int niter, double th, int datatype);
/**
 * \brief Wrapper for binary morphological closing (generic datatype).
 *
 * \param data       (in/out) void pointer to image data; type given by datatype parameter
 * \param dims       (in)     {nx, ny, nz}
 * \param niter      (in)     number of iterations for each operation (<=0: no-op)
 * \param th         (in)     threshold as fraction of max(data) in [0,1]
 * \param datatype   (in)     data type code (DT_UINT8, DT_UINT16, DT_FLOAT32, etc.)
 */
void morph_close(void *data, int dims[3], int niter, double th, int datatype);
/**
 * \brief Wrapper for binary morphological opening (generic datatype).
 *
 * \param data        (in/out) void pointer to image data; type given by datatype parameter
 * \param dims        (in)     {nx, ny, nz}
 * \param niter       (in)     number of iterations for each operation (<=0: no-op)
 * \param th          (in)     threshold as fraction of max(data) in [0,1]
 * \param keep_values (in)     if >0, preserve original values and zero only removed regions
 * \param datatype    (in)     data type code (DT_UINT8, DT_UINT16, DT_FLOAT32, etc.)
 */
void morph_open(void *data, int dims[3], int niter, double th, int keep_values,
                int datatype);
/**
 * \brief Grey-scale morphological erosion.
 *
 * \param data       (in/out) void pointer to image data; type given by datatype parameter
 * \param dims       (in)     {nx, ny, nz}
 * \param niter      (in)     number of erosion iterations (<=0: no-op)
 * \param datatype   (in)     data type code (DT_UINT8, DT_UINT16, DT_FLOAT32, etc.)
 */
void grey_erode(void *data, int dims[3], int niter, int datatype);
/**
 * \brief Grey-scale morphological dilation.
 *
 * \param data       (in/out) void pointer to image data; type given by datatype parameter
 * \param dims       (in)     {nx, ny, nz}
 * \param niter      (in)     number of dilation iterations (<=0: no-op)
 * \param datatype   (in)     data type code (DT_UINT8, DT_UINT16, DT_FLOAT32, etc.)
 */
void grey_dilate(void *data, int dims[3], int niter, int datatype);
/**
 * \brief Grey-scale morphological opening (erosion followed by dilation).
 *
 * \param data       (in/out) void pointer to image data; type given by datatype parameter
 * \param dims       (in)     {nx, ny, nz}
 * \param niter      (in)     number of iterations for each operation (<=0: no-op)
 * \param datatype   (in)     data type code (DT_UINT8, DT_UINT16, DT_FLOAT32, etc.)
 */
void grey_open(void *data, int dims[3], int niter, int datatype);
/**
 * \brief Grey-scale morphological closing (dilation followed by erosion).
 *
 * \param data       (in/out) void pointer to image data; type given by datatype parameter
 * \param dims       (in)     {nx, ny, nz}
 * \param niter      (in)     number of iterations for each operation (<=0: no-op)
 * \param datatype   (in)     data type code (DT_UINT8, DT_UINT16, DT_FLOAT32, etc.)
 */
void grey_close(void *data, int dims[3], int niter, int datatype);
/**
 * \brief Wrapper for morphological closing (generic datatype).
 *
 * \param data       (in/out) void pointer to image data; type given by datatype parameter
 * \param dims       (in)     {nx, ny, nz}
 * \param voxelsize  (in)     voxel spacing in mm (or consistent units)
 * \param dist       (in)     structuring radius in same units as voxelsize (<=0: no-op)
 * \param th         (in)     threshold as fraction of max(data) in [0,1]
 * \param datatype   (in)     data type code (DT_UINT8, DT_UINT16, DT_FLOAT32, etc.)
 */
void dist_close(void *data, int dims[3], double voxelsize[3], double dist,
                double th, int datatype);
/**
 * \brief Morphological closing (binary) using the Euclidean distance transform.
 *
 * \param vol        (in/out) float[dims[0]*dims[1]*dims[2]]; overwritten with 0/1
 * \param dims       (in)     {nx, ny, nz}
 * \param voxelsize  (in)     voxel spacing in mm (or consistent units)
 * \param dist       (in)     structuring radius in same units as voxelsize (<=0: no-op)
 * \param th         (in)     threshold as fraction of max(vol) in [0,1]
 * \param mask       (in)     optional uint8 ROI mask (same dims); NULL = full volume.
 *                                When provided, the close is decomposed into a masked
 *                                dilation followed by a masked erosion for efficiency.
 */
void dist_close_float(float *vol, int dims[3], double voxelsize[3], double dist,
                      double th, unsigned char *mask);
/**
 * \brief Wrapper for morphological opening (generic datatype).
 *
 * \param data       (in/out) void pointer to image data; type given by datatype parameter
 * \param dims       (in)     {nx, ny, nz}
 * \param voxelsize  (in)     voxel spacing in mm (or consistent units)
 * \param dist       (in)     structuring radius in same units as voxelsize (<=0: no-op)
 * \param th         (in)     threshold as fraction of max(data) in [0,1]
 * \param datatype   (in)     data type code (DT_UINT8, DT_UINT16, DT_FLOAT32, etc.)
 */
void dist_open(void *data, int dims[3], double voxelsize[3], double dist,
               double th, int datatype);
/**
 * \brief Morphological opening (binary) using the Euclidean distance transform.
 *
 * \param vol        (in/out) float[dims[0]*dims[1]*dims[2]]; overwritten with 0/1
 * \param dims       (in)     {nx, ny, nz}
 * \param voxelsize  (in)     voxel spacing in mm (or consistent units)
 * \param dist       (in)     structuring radius in same units as voxelsize (<=0: no-op)
 * \param th         (in)     threshold as fraction of max(vol) in [0,1]
 * \param mask       (in)     optional uint8 ROI mask (same dims); NULL = full volume.
 *                                When provided, the open is decomposed into a masked
 *                                erosion followed by a masked dilation for efficiency.
 */
void dist_open_float(float *vol, int dims[3], double voxelsize[3], double dist,
                     double th, unsigned char *mask);
/**
 * \brief Wrapper for morphological erosion (generic datatype).
 *
 * \param data       (in/out) void pointer to image data; type given by datatype parameter
 * \param dims       (in)     {nx, ny, nz}
 * \param voxelsize  (in)     voxel spacing in mm (or consistent units)
 * \param dist       (in)     structuring radius in same units as voxelsize (<=0: no-op)
 * \param th         (in)     threshold as fraction of max(data) in [0,1]
 * \param datatype   (in)     data type code (DT_UINT8, DT_UINT16, DT_FLOAT32, etc.)
 */
void dist_erode(void *data, int dims[3], double voxelsize[3], double dist,
                double th, int datatype);
/**
 * \brief Morphological erosion (binary) using the Euclidean distance transform.
 *
 * \param vol        (in/out) float[dims[0]*dims[1]*dims[2]]; overwritten with 0/1
 * \param dims       (in)     {nx, ny, nz}
 * \param voxelsize  (in)     voxel spacing in mm (or consistent units)
 * \param dist       (in)     structuring radius in same units as voxelsize (<=0: no-op)
 * \param th         (in)     threshold as fraction of max(vol) in [0,1]
 * \param mask       (in)     optional uint8 ROI mask (same dims); NULL = full volume.
 *                                Voxels with mask==0 keep their original foreground/background
 *                                classification and are excluded from the EDT sweep, which
 *                                speeds up computation when most of the volume is irrelevant.
 */
void dist_erode_float(float *vol, int dims[3], double voxelsize[3], double dist,
                      double th, unsigned char *mask);
/**
 * \brief Wrapper for morphological dilation (generic datatype).
 *
 * \param data       (in/out) void pointer to image data; type given by datatype parameter
 * \param dims       (in)     {nx, ny, nz}
 * \param voxelsize  (in)     voxel spacing in mm (or consistent units)
 * \param dist       (in)     structuring radius in same units as voxelsize (<=0: no-op)
 * \param th         (in)     threshold as fraction of max(data) in [0,1]
 * \param datatype   (in)     data type code (DT_UINT8, DT_UINT16, DT_FLOAT32, etc.)
 */
void dist_dilate(void *data, int dims[3], double voxelsize[3], double dist,
                 double th, int datatype);
/**
 * \brief Morphological dilation (binary) using the Euclidean distance transform.
 *
 * \param vol        (in/out) float[dims[0]*dims[1]*dims[2]]; overwritten with 0/1
 * \param dims       (in)     {nx, ny, nz}
 * \param voxelsize  (in)     voxel spacing in mm (or consistent units)
 * \param dist       (in)     structuring radius in same units as voxelsize (<=0: no-op)
 * \param th         (in)     threshold as fraction of max(vol) in [0,1]
 * \param mask       (in)     optional uint8 ROI mask (same dims); NULL = full volume.
 *                                Voxels with mask==0 keep their original foreground/background
 *                                classification and are excluded from the EDT sweep.
 */
void dist_dilate_float(float *vol, int dims[3], double voxelsize[3],
                       double dist, double th, unsigned char *mask);
/**
 * \brief Resample a 3D volume to a different size using trilinear interpolation.
 *
 * \param in        (in)  input volume data (pointer to any supported datatype)
 * \param out       (out) output volume data (pointer to pre-allocated array)
 * \param dims      (in)  original volume dimensions {nx, ny, nz}
 * \param dims_samp (in)  target volume dimensions {nx_new, ny_new, nz_new}
 * \param datatype  (in)  data type descriptor (e.g., DT_FLOAT32, DT_UINT8)
 */
void subsample3(void *in, void *out, int dims[3], int dims_samp[3],
                int datatype);
/**
 * \brief Approximate missing values in a volume by interpolating from neighbors.
 *
 * \param vol       (in/out) float[dims[0]*dims[1]*dims[2]]; modified in place
 * \param dims      (in)     {nx, ny, nz}
 * \param voxelsize (in)     voxel spacing in mm (or consistent units)
 */
void vol_approx(float *vol, int dims[3], double voxelsize[3]);
/**
 * \brief Clean up tissue probability map by morphological refinement.
 *
 * \param prob     (in/out) unsigned char[3*nvox]; tissue probability array
 *                              [0:nvox-1]=CSF, [nvox:2*nvox-1]=GM, [2*nvox:3*nvox-1]=WM
 *                              Modified in-place by cleanup operations
 * \param dims     (in)     {nx, ny, nz} volume dimensions
 * \param voxelsize (in)     {sx, sy, sz} voxel spacing in mm; used for morphological scaling
 * \param strength  (in)     cleanup strength (0..N); controls dilation threshold
 *                               (higher = more aggressive cleanup)
 */
void cleanup_brain(unsigned char *prob, int dims[3], double voxelsize[3],
                   int strength);
/**
 * \brief Euclidean distance transform (see euclidean_distance_src()).
 *
 * \param V         (in/out) float volume; positive values are distance sources
 * \param M         (in)     optional uint8 mask (same dims); NULL = all-ones
 * \param dims      (in)     {nx, ny, nz}
 * \param voxelsize (in)     voxel spacing; NULL -> {1,1,1}
 * \param replace   (in)     0 = output distances; >0 = output nearest values
 */
void euclidean_distance(float *V, unsigned char *M, int dims[3],
                        double *voxelsize, int replace);
/**
 * \brief euclidean_distance() with the value at the nearest source voxel.
 *
 * Same as euclidean_distance(), but when `src` and `src_out` are both non-NULL
 * `src_out` additionally receives `src[nearest source voxel]`. The distance is
 * measured centre-to-centre; the partial volume carried in `src_out` tells the
 * caller how far the boundary lies beyond that centre, which is what turns a
 * centre-to-centre distance into a centre-to-boundary distance.
 *
 * \param V (in/out) float volume; positive values are distance sources.
 * \param M (in) optional uint8 mask; NULL = all-ones.
 * \param dims (in) {nx, ny, nz}.
 * \param voxelsize (in) voxel spacing; NULL -> {1,1,1}.
 * \param replace (in) 0 = output distances; >0 = output nearest values.
 * \param src (in) optional array sampled at the nearest source; NULL to skip.
 * \param src_out (out) optional array receiving src[nearest source]; NULL to skip.
 */
void euclidean_distance_src(float *V, unsigned char *M, int dims[3], double *voxelsize,
                            int replace, const float *src, float *src_out);
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
 * \brief Convert a linear index to 3D array coordinates.
 *
 * \param i The linear index in the array.
 * \param x Pointer to store the calculated x-coordinate.
 * \param y Pointer to store the calculated y-coordinate.
 * \param z Pointer to store the calculated z-coordinate.
 * \param sxy Product of the dimensions in the x and y directions (sx * sy).
 * \param sx The dimension in the x direction.
 */
void ind2sub(int i, int *x, int *y, int *z, int sxy, int sx);
/**
 * \brief Convert 3D array coordinates to a linear index.
 *
 * \param x The x-coordinate in the array.
 * \param y The y-coordinate in the array.
 * \param z The z-coordinate in the array.
 * \param s Array containing the dimensions of the 3D array.
 * \return The linear index corresponding to the provided 3D coordinates.
 */
int sub2ind(int x, int y, int z, int s[3]);
/**
 * \brief Connected-component filter on a thresholded volume.
 *
 * Labels the components of `inData >= thresh` and then either keeps only the
 * largest of them or keeps every component above a size floor, depending on the
 * sign of `min_size`.
 *
 * \param data            (in/out) volume data in `datatype`, filtered in place
 * \param thresh          (in)     voxels >= thresh are cluster members
 * \param dims            (in)     {nx, ny, nz}
 * \param datatype        (in)     datatype code of inData (e.g. DT_FLOAT32)
 * \param min_size        (in)     >=0 keeps only the largest cluster and the
 *                                 magnitude is unused; <0 keeps every cluster of
 *                                 at least |min_size| voxels
 * \param retain_above_th (in)     1 = zero the rejected clusters and leave the
 *                                 surviving values untouched; 0 = additionally
 *                                 zero everything below thresh
 * \param conn            (in)     connectivity: 6, 18 or 26
 */
void keep_largest_cluster(void *data, double thresh, int *dims, int datatype, int min_size, int retain_above_th, int conn);
/**
 * \brief Fill holes in a binary or thresholded volume.
 *
 * \param data      (in/out) void pointer to volume data; type given by datatype parameter
 * \param dims      (in)     {nx, ny, nz} volume dimensions
 * \param thresh    (in)     threshold value; voxels < thresh are treated as potential holes
 * \param fill_value (in)    value to fill holes with;
 *                                if negative, holes are filled with locally estimated values
 *                                if >=0, holes are filled with this fixed value
 * \param datatype  (in)     data type code (DT_UINT8, DT_UINT16, DT_FLOAT32, etc.)
 */
void fill_holes(void *data, int *dims, double thresh, double fill_value,
                int datatype);
/**
 * \brief Compute local gradient magnitude and components for a 3D volume.
 *
 * \param src       (in)  input volume float[dims[0]*dims[1]*dims[2]]
 * \param grad_mag  (out) gradient magnitude; NULL to skip (optional)
 * \param grad_x    (out) x-component of gradient; NULL to skip (optional)
 * \param grad_y    (out) y-component of gradient; NULL to skip (optional)
 * \param grad_z    (out) z-component of gradient; NULL to skip (optional)
 * \param dims      (in)  volume dimensions {nx, ny, nz}
 * \param voxelsize (in)  voxel spacing in mm {dx, dy, dz}
 */
void gradient3D(float *src, float *grad_mag, float *grad_x, float *grad_y,
                float *grad_z, int dims[3], double voxelsize[3]);

/**
 * \brief Matrix that maps a gradient3D() gradient into world space.
 *
 * gradient3D() differentiates along the voxel axes and divides by the voxel
 * size.  Surface normals and world positions live in world space, so the
 * gradient has to be rotated before it is compared with them: using it
 * directly flips its sign on every axis stored with a negative direction and
 * ignores the rotation of oblique images.
 *
 * \param nii_ptr (in)  NIfTI header (sto_xyz and voxel size dx, dy, dz)
 * \param M       (out) 3x3 matrix, g_world = M * g_gradient3D
 */
void gradient3D_world_matrix(const nifti_image *nii_ptr, double M[3][3]);
#endif
