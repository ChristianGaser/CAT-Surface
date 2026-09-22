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
/* fewest voxels in a subvolume for which a class mean and variance are estimated */
#define AMAP_MIN_VOXELS 6

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

typedef struct {
    /* inputs */
    const float *src;
    const unsigned char *label;
    int n_classes;
    int sub;
    const int *dims;         /* dims[0]=X, dims[1]=Y, dims[2]=Z */
    const double *thresh;    /* thresh[0] .. (optional thresh[1]) */
    /* grid geometry (precomputed) */
    int nix, niy, niz, narea, nvol, area;
    /* shared output accumulator */
    struct ipoint *ir;       /* size: n_classes * nvol */
    /* work partition on grid-z */
    int z_ini, z_fin;        /* [z_ini, z_fin) in 0..niz */
} gmv_accum_args_t;

typedef struct {
    struct point *r;
    const struct ipoint *ir;
    int n_classes;
    int nvol;
    int use_median;
    int j_ini, j_fin;   /* j-range [j_ini, j_fin) */
} gmv_reduce_args_t;

/**
 * \brief Adaptive Segmentation atlas Mapping (Amap): tissue classification via EM.
 *
 * \param src              (in)  input MRI intensity image
 * \param label            (out) hard tissue classification labels (1=CSF, 3=GM, 5=WM)
 * \param prob             (out) soft tissue probability maps (n_classes*nvol)
 * \param mean             (in/out) class mean intensity estimates; updated in-place
 * \param n_classes        (in)  number of tissue classes (typically 3-5)
 * \param niters           (in)  maximum EM iterations
 * \param sub              (in)  subsampling factor for speed/accuracy trade-off
 * \param dims             (in)  array [nx, ny, nz] volume dimensions
 * \param pve              (in)  1 for partial volume estimation, 0 to skip
 * \param weight_MRF       (in)  MRF regularization strength (0=no smoothing, 1=strong)
 * \param voxelsize        (in)  array [dx, dy, dz] voxel dimensions in mm
 * \param niters_ICM       (in)  iterations of ICM mode refinement per EM step
 * \param verbose          (in)  1 to print progress, 0 for silent
 * \param use_median       (in)  1 to use median in class statistics, 0 for mean only
 * \param mrf_class_weights (in) per-class MRF weights or NULL for uniform
 * \param use_multistep    (in)  1 for multi-resolution coarse-to-fine, 0 for single
 */
void Amap(float *src, unsigned char *label, unsigned char *prob, double *mean,
          int n_classes, int niters, int sub, int *dims, int pve,
          double weight_MRF, double *voxelsize, int niters_ICM, int verbose,
          int use_median, const double *mrf_class_weights, int use_multistep);
/**
 * \brief Convert tissue classification to partial volume estimates (CSF/GM/WM).
 *
 * \param src     (in)  input intensity image
 * \param prob    (out) probability maps (3*nvol length: CSF, GM, WM stacked)
 * \param label   (in/out) tissue labels; updated with PVE intensity estimates
 * \param mean    (in)  tissue class means [CS, GM, WM]
 * \param dims    (in)  array [nx, ny, nz] specifying volume dimensions
 */
void Pve5(float *src, unsigned char *prob, unsigned char *label, double *mean,
          int *dims);
/**
 * \brief Compute Gaussian probability density at a given value.
 *
 * \param value  (in)  measured intensity value
 * \param mean   (in)  tissue class mean intensity
 * \param var    (in)  tissue class variance (spread)
 * \return Probability density p(value | mean, variance)
 */
double ComputeGaussianLikelihood(double value, double mean, double var);
/**
 * \brief Compute likelihood for mixed-tissue voxels via marginalized integration.
 *
 * \param value           (in)  measured intensity value
 * \param mean1           (in)  first tissue class mean
 * \param mean2           (in)  second tissue class mean
 * \param var1            (in)  first tissue class variance
 * \param var2            (in)  second tissue class variance
 * \param nof_intervals   (in)  number of integration steps (higher = more accurate, slower)
 * \return Marginalized likelihood p(value | tissue1, tissue2)
 */
double ComputeMarginalizedLikelihood(double value, double mean1, double mean2,
                                     double var1, double var2,
                                     unsigned int nof_intervals);
/**
 * \brief Estimate Markov Random Field (MRF) prior parameters from label configuration.
 *
 * \param label      (in)  tissue classification labels (hard assignments)
 * \param n_classes  (in)  total number of classes
 * \param alpha      (out) class prevalence/frequency array (length n_classes)
 * \param beta       (out) strength parameter [1] for MRF smoothing; NULL to skip
 * \param init       (in)  if 1, initialize alphas to 1.0; if 0, compute from data
 * \param dims       (in)  array [nx, ny, nz] specifying volume dimensions
 * \param verbose    (in)  1 to print parameters to stdout, 0 for silent
 */
void MrfPrior(unsigned char *label, int n_classes, double *alpha, double *beta,
              int init, int *dims, int verbose);
/**
 * \brief Scale array values to sum to 1.0 (probability normalization).
 *
 * \param val  (in/out) array of values to normalize; modified in-place
 * \param n    (in)  length of array
 */
void Normalize(double *val, char n);
/**
 * \brief Find index of maximum value in array.
 *
 * \param val  (in)  array of likelihood or probability values
 * \param n    (in)  length of array
 * \return 1-indexed class label (tissue class), 1 to n
 */
unsigned char MaxArg(double *val, unsigned char n);

struct point {
  int n;
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
