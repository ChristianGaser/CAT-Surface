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

/**
 * \brief Arguments of one thread that accumulates class statistics on the grid.
 *
 * The volume is divided into sub x sub x sub blocks; every thread takes a
 * disjoint range of grid planes and sums the intensities of each class into its
 * own block of ir, so no locking is needed.
 */
typedef struct {
    /* inputs */
    const float *src;            /**< intensity image */
    const unsigned char *label;  /**< current hard labels, one per voxel */
    int n_classes;               /**< number of classes accumulated */
    int sub;                     /**< block size of the grid in voxels */
    const int *dims;             /**< volume dimensions {nx, ny, nz} */
    const double *thresh;        /**< lower (and optional upper) intensity bound */
    /* grid geometry (precomputed) */
    int nix;                     /**< grid blocks along x */
    int niy;                     /**< grid blocks along y */
    int niz;                     /**< grid blocks along z */
    int narea;                   /**< blocks per grid plane, nix * niy */
    int nvol;                    /**< blocks in the grid, narea * niz */
    int area;                    /**< voxels per volume slice, nx * ny */
    /* shared output accumulator */
    struct ipoint *ir;           /**< per-class block sums, n_classes * nvol entries */
    /* work partition on grid-z */
    int z_ini;                   /**< first grid plane of this thread's range */
    int z_fin;                   /**< one past its last grid plane, at most niz */
} gmv_accum_args_t;

/**
 * \brief Arguments of one thread that turns the accumulated sums into statistics.
 *
 * Converts the sums in ir into the mean and variance (or median) per class and
 * block; every thread takes a disjoint range of blocks.
 */
typedef struct {
    struct point *r;             /**< per-class block statistics, n_classes * nvol */
    const struct ipoint *ir;     /**< accumulated sums from gmv_accum_args_t */
    int n_classes;               /**< number of classes */
    int nvol;                    /**< blocks in the grid */
    int use_median;              /**< non-zero to use the median instead of the mean */
    int j_ini;                   /**< first block of this thread's range */
    int j_fin;                   /**< one past its last block */
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

/**
 * \brief Statistics of one class within one grid block.
 */
struct point {
  int n;          /**< number of voxels of the class in the block, 0 if too few */
  double median;  /**< median intensity (unused unless use_median) */
  double mean;    /**< mean intensity, or the median with use_median */
  double var;     /**< sample variance of the intensities */
};

/**
 * \brief Intensity sums of one class within one grid block.
 *
 * Filled by the accumulation pass and reduced into a struct point afterwards.
 */
struct ipoint {
  int n;        /**< number of voxels accumulated */
  double s;     /**< sum of the intensities */
  double ss;    /**< sum of the squared intensities */
  double *arr;  /**< the intensities themselves, needed for the median */
};

#endif
