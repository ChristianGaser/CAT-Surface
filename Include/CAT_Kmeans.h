/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 *
 */

#ifndef _CAT_KMEANS_H_
#define _CAT_KMEANS_H_

#include "CAT_Vol.h"

#ifndef HUGE
#define HUGE 1e15 
#endif

/**
 * \brief Perform K-means clustering on image data for tissue segmentation.
 *
 * \param src           (in)  float[nvox]; source image intensity values
 * \param label         (out) unsigned char[nvox]; cluster assignments (can be NULL)
 * \param mask          (in)  unsigned char[nvox]; voxel inclusion mask (can be NULL)
 * \param NI            (in)  number of K-means iterations
 * \param n_clusters    (in)  maximum number of clusters to try (2..n_clusters)
 * \param mean          (out) double[n_clusters]; cluster mean intensities
 * \param voxelsize     (in)  double[3]; voxel spacing in mm (unused in current impl)
 * \param dims          (in)  int[3]; volume dimensions {nx, ny, nz}
 * \param thresh_mask   (in)  mask threshold; voxels with mask < thresh_mask excluded
 * \param thresh_kmeans (in)  K-means threshold for clustering refinement
 * \return                    Maximum intensity in the source image
 */
double Kmeans(float *src, unsigned char *label, unsigned char *mask, int NI,
              int n_clusters, double *mean, double *voxelsize, int *dims,
              int thresh_mask, int thresh_kmeans);

#endif
