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

double Kmeans(float *src, unsigned char *label, unsigned char *mask, int NI, int n_clusters, double *mean, double *voxelsize, int *dims, int thresh_mask, int thresh_kmeans);

#endif