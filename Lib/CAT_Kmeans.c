/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 *
 */

/* This code is a substantially modified version of Tskmeans.C 
 * from Jagath C. Rajapakse
 * 
 * Original author : Jagath C. Rajapakse
 *
 * See:
 * Statistical approach to single-channel MR brain scans
 * J. C. Rajapakse, J. N. Giedd, and J. L. Rapoport
 * IEEE Transactions on Medical Imaging, Vol 16, No 2, 1997
 *
 * Tree structure k-means algorithm
 *
 * Jagath C. Rajapakse (raja@cns.mpg.de) 23-07-97
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "CAT_Kmeans.h"

/**
 * \brief Estimate K-means cluster centers using histogram-based refinement.
 *
 * Performs iterative K-means refinement on quantized intensity histogram. Initializes
 * cluster assignments based on intensity and iteratively updates cluster means until
 * convergence or maximum iterations reached. Operates on 256-bin histogram for efficiency,
 * then maps final labels back to original voxel data.
 *
 * Algorithm:
 *  1. Build histogram of source intensities (0-255 scale)
 *  2. Exclude histogram tails (< 1% or > 99% cumsum)
 *  3. Iterate: assign histogram bins to nearest cluster mean, update means
 *  4. Build lookup table mapping histogram bins to cluster labels
 *  5. Apply lookup table to source voxels, respecting mask constraints
 *  6. Return accumulated squared error of final classification
 *
 * \param src           (in)  float[nvox]; source image intensity values
 * \param label         (out) unsigned char[nvox]; cluster assignments (can be NULL)
 * \param mask          (in)  unsigned char[nvox]; voxel inclusion mask (can be NULL)
 * \param n_classes     (in)  number of clusters (classes) to refine
 * \param mean          (in/out) double[n_classes]; cluster mean intensities (updated)
 * \param ni            (in)  maximum number of iterations
 * \param dims          (in)  int[3]; volume dimensions {nx, ny, nz}
 * \param thresh_mask   (in)  mask threshold; voxels with mask < thresh_mask get label 0
 * \param thresh_kmeans (in)  threshold for histogram filtering
 * \param max_src       (in)  maximum source intensity (used for 0-255 scaling)
 * \return                    Sum of squared errors for final classification
 */
double EstimateKmeans(float *src, unsigned char *label, unsigned char *mask, int n_classes, double *mean, int ni, int *dims, int thresh_mask, int thresh_kmeans, double max_src)
{
    int i, j, j0, v;
    int count;
    long histo[256], lut[256], cumsum[256], vol;
    double diff, dmin, dx, xnorm, sum;

    vol  = dims[0]*dims[1]*dims[2];

    /* build intensity histogram */
    for (i = 0; i < 256; i++) histo[i] = 0;
    for (i = 0; i < vol; i++) {
        v = (int)ROUND(255.0*(double)src[i]/max_src);
        if (v < 1) continue;
        if (mask)
            if ((thresh_mask > 0) && ((int)mask[i] < thresh_kmeans))
                continue;
        if (v < 0) v = 0;
        if (v > 255) v = 255;   
        histo[v]++;  
    }

    /* use only value in histogram where cumsum is between 1..99% */
    cumsum[0] = histo[0];
    for (i = 1; i < 256; i++) cumsum[i] = cumsum[i-1] + histo[i];
    for (i = 0; i < 256; i++) cumsum[i] = (long) ROUND(1000.0*(double)cumsum[i]/(double)cumsum[255]);
    for (i = 0; i < 256; i++) if ((cumsum[i] <= 10) || (cumsum[i] >= 990)) histo[i] = 0;

    /* loop through */
    diff = HUGE;    count = 0;
    while (diff > 1.0 && count < ni) {

        /* assign class labels */
        for (i = 0; i < 256; i++) {
            dmin = 256.0 * 256.0;
            for (j = 0; j < n_classes; j++) {
                 dx = (double) i - mean[j];
                 dx *= dx;
                 if (dx < dmin) {
                     lut[i] = j;
                     dmin = dx;
                 }
            }
        }

        /* find the new cluster centers */
        diff = 0;
        for (i = 0; i < n_classes; i++) {
            xnorm = 0.0;
            sum = 0.0;
            for (j = 0; j < 256; j++)
             if (lut[j] == i) {
                 xnorm += histo[j];
                 sum +=  j * histo[j];
             }
            if (xnorm > 0) sum /= xnorm;
            else sum = 0.0;
            dx = sum - mean[i];
            mean[i] = sum;
            dx *= dx;
            diff += dx;
        }
        count++;
    }

    /* assign final labels to voxels */
    for (i = 0; i < 256; i++) {
        dmin = HUGE;
        j0 = 0;
        for (j = 0; j < n_classes; j++) {
            if (fabs((double) i - mean[j]) < dmin) {
                 dmin = fabs((double)i - mean[j]);
                 j0 = j;
            }
        }
        lut[i] = j0;
    }
    
    lut[0] = 0;

    /* adjust for the background label */
    diff = 0;
    
    for (i = 0; i < vol; i++) {
        v = (int)ROUND(255.0*(double)src[i]/max_src);
        if (v >= 1) {
            if (v < 0) v = 0;
            if (v > 255) v = 255;
            if (label) label[i] = (unsigned char)(lut[v] + 1); 
            diff += SQR((double)v - mean[lut[v]]);
            if (mask && label)
                if ((thresh_mask > 0) && ((int)mask[i] < thresh_mask))
                    label[i] = 0;   
        }
        else
            if (label) label[i] = 0;  
    }
    
    /* return square error */
    return(diff);
}

/**
 * \brief Perform K-means clustering on image data for tissue segmentation.
 *
 * Implements a K-means clustering algorithm for image segmentation, iteratively
 * refining cluster centroids based on intensity values. Tests cluster counts from
 * 2 up to n_clusters, selecting the optimal number based on variance reduction.
 * Clustering can be constrained by optional mask and label arrays.
 *
 * Algorithm:
 *  1. Find maximum intensity within masked region
 *  2. For each cluster count (2 to n_clusters):
 *  3.   Initialize centroids (Otsu's method for 2 clusters)
 *  4.   Run K-means for NI iterations, updating cluster assignments and means
 *  5.   Compute log-likelihood to select optimal cluster count
 *  6. Return maximum intensity found
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
double Kmeans(float *src, unsigned char *label, unsigned char *mask, int NI, int n_clusters, double *mean, double *voxelsize, int *dims, int thresh_mask, int thresh_kmeans)
{
    int i, j, k;
    double e, emin, eps, val, sum;
    double max_src = -HUGE;
    long n[n_clusters];
    double var[MAX_NC],  mu[MAX_NC],  Mu[MAX_NC];
    int n_classes;
    long vol;

    vol  = dims[0]*dims[1]*dims[2];

    /* find maximum and mean inside mask */
    for (i = 0; i < vol; i++) {
        if (mask) {
            if (mask[i] > 0)
                max_src = MAX((double)src[i], max_src);
        } else
            max_src = MAX((double)src[i], max_src);
    }

    /* go through all sizes of cluster beginning with two clusters */
    for (n_classes=2; n_classes <= n_clusters; n_classes++) {

        if (n_classes == 2) {
            /* initialize for the two cluster case; */
            n[0]=0; mean[0] = 0.0; var[0] = 0.0;

            for (i = 0; i < vol; i++) {
                val = 255.0*(double)src[i]/max_src;
                if (val < 1.0/255.0) continue;
                n[0]++;
                mean[0] += val;
                var[0]  += SQR(val);
            }
            
            Mu[0] = (n[0] != 0) ? mean[0]/n[0]: 0.0;
            var[0] = (n[0] > 1) ? (var[0] - n[0]*Mu[0]*Mu[0])/(n[0] - 1.0) : 1.0;
            eps = 0.5*sqrt(var[0]);
        }
        else {
            /* find the deviant (epsilon) for the node being divided */
            eps = Mu[0];
            for (i = 0; i < n_classes-2; i++)
                if (Mu[i+1] - Mu[i] < eps)
                    eps = Mu[i+1] - Mu[i];
            if (255 - Mu[n_classes-2] < eps)
                eps = 255 - Mu[n_classes-2];
            eps = eps*0.5;
        }


        /* go through low order clustering */
        emin = HUGE;
        for (k = 0; k < n_classes-1; k++) {
            for (i = n_classes-1; i > k+1; i--) mean[i] = Mu[i-1];
            mean[k+1] = Mu[k] + eps;
            mean[k] = Mu[k] - eps;
            for (i = 1; i < k; i++) mean[i] = Mu[i];
            e = EstimateKmeans(src, label, mask, n_classes, mean, NI, dims, thresh_mask, thresh_kmeans, max_src);
            if (e < emin) {
                emin = e;
                for (i = 0; i < n_classes; i++) 
                    mu[i] = mean[i];
            }
        }
        for (i = 0; i < n_classes; i++) Mu[i] = mu[i];       
    }

    e = EstimateKmeans(src, label, mask, n_clusters, mu, NI, dims, thresh_mask, thresh_kmeans, max_src);

    max_src = -HUGE;
    for (i = 0; i < vol; i++)
        max_src = MAX((double)src[i], max_src);

    for (i = 0; i < n_classes; i++) mean[i] = max_src*mu[i]/255.0;       
    
    return(max_src);        
}













