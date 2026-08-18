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

/*
 * Anisotropic MRF regularization.
 *
 * The Potts prior used by ICM() penalizes label disagreement equally in every
 * direction. Any such isotropic prior pays a cost proportional to boundary
 * area, and a thin structure has an extreme area-to-volume ratio, so removing
 * it is always the cheaper labelling -- the "shrinking bias". That is why a
 * stronger MRF closes a cerebellar fissure and a thin sulcus at the same rate
 * as it removes genuine noise.
 *
 * Two ways out are offered here, both driven by the Hessian sheetness field of
 * CAT_Sheetness.h and both no-ops where the sheetness is zero, so the
 * behaviour away from thin structures is exactly the behaviour without them:
 *
 *   MRF_ANISO_BETA  local beta:  beta(x) = beta * (1 - strength * s(x)).
 *                   The regularizer simply stops pulling where a thin sheet
 *                   was detected. Cheapest, and needs nothing but s(x).
 *
 *   MRF_ANISO_POTTS direction-dependent weights: a neighbour at offset d is
 *                   weighted by (1 - s) + s * exp(-(dhat . n)^2 / 2 sigma^2),
 *                   with n the local sheet normal. Neighbours *along* the
 *                   sheet keep full influence, neighbours *across* it lose it,
 *                   so the prior still denoises within a bank but can no
 *                   longer merge two banks facing each other.
 *
 * References
 * ----------
 * Weickert, Int J Comput Vis 31(2/3):111-127, 1999 (coherence-enhancing).
 * Ribal et al., IPMU 2020:601-612 (anisotropic MRF neighbourhoods).
 */
#define MRF_ANISO_OFF   0
#define MRF_ANISO_BETA  1
#define MRF_ANISO_POTTS 2

/** \brief Sheetness field steering the anisotropic MRF prior. */
typedef struct {
    int mode;                /**< MRF_ANISO_OFF, _BETA or _POTTS */
    const float *sheetness;  /**< [nvox] in [0,1]; required when mode > 0 */
    const float *normal;     /**< [3*nvox] unit sheet normals; required for _POTTS */
    double strength;         /**< 0..1 how far the prior is relaxed on a sheet */
    double sigma;            /**< angular width for _POTTS, e.g. 0.5 */
} CAT_MrfAnisoField;

void Amap(float *src, unsigned char *label, unsigned char *prob, double *mean, 
                 int nc, int niters, int sub, int *dims, int pve, double weight_MRF, 
                 double *voxelsize, int niters_ICM, int verbose, 
                 int use_median, const double *mrf_class_weights, int use_multistep);
/**
 * \brief Amap() with an anisotropic MRF prior.
 *
 * Identical to Amap() in every respect except that the ICM stage consults
 * `aniso`. Passing NULL, or a field with mode MRF_ANISO_OFF, reproduces
 * Amap() bit for bit.
 *
 * \param src              (in)  input MRI intensity image
 * \param label            (out) hard tissue classification labels
 * \param prob             (out) soft tissue probability maps
 * \param mean             (in/out) class mean intensity estimates
 * \param nc               (in)  number of tissue classes
 * \param niters           (in)  maximum EM iterations
 * \param sub              (in)  subsampling factor
 * \param dims             (in)  {nx, ny, nz}
 * \param pve              (in)  1 for partial volume estimation, 0 to skip
 * \param weight_MRF       (in)  MRF regularization strength
 * \param voxelsize        (in)  voxel dimensions in mm
 * \param niters_ICM       (in)  ICM iterations
 * \param verbose          (in)  1 to print progress
 * \param use_median       (in)  1 to use median in class statistics
 * \param mrf_class_weights (in) per-class MRF weights or NULL
 * \param use_multistep    (in)  1 for multi-resolution coarse-to-fine
 * \param aniso            (in)  anisotropy field, or NULL for the isotropic prior
 * \return void
 */
void AmapAniso(float *src, unsigned char *label, unsigned char *prob, double *mean,
                 int nc, int niters, int sub, int *dims, int pve, double weight_MRF,
                 double *voxelsize, int niters_ICM, int verbose,
                 int use_median, const double *mrf_class_weights, int use_multistep,
                 const CAT_MrfAnisoField *aniso);
/**
 * \brief Public API for Pve5.
 *
 * This function is part of the CAT-Surface public library interface and is used by command-line tools.
 *
 * \param src (in/out) Parameter of Pve5.
 * \param prob (in/out) Parameter of Pve5.
 * \param label (in/out) Parameter of Pve5.
 * \param mean (in/out) Parameter of Pve5.
 * \param dims (in/out) Parameter of Pve5.
 * \return void (no return value).
 */
void Pve5(float *src, unsigned char *prob, unsigned char *label, double *mean, int *dims);
double ComputeGaussianLikelihood(double value, double mean , double var);
double ComputeMarginalizedLikelihood(double value, double mean1 , double mean2, 
                double var1, double var2, unsigned int nof_intervals);
void MrfPrior(unsigned char *label, int n_classes, double *alpha, double *beta, 
                int init, int *dims, int verbose);
void Normalize(double* val, char n);
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
