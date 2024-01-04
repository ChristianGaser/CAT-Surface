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

/* perform k-means algorithm given initial mean estimates */      
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
            //sum = (xnorm > 0) ? sum /= xnorm : 0.0;
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
 * Kmeans - Performs K-means clustering on image data.
 *
 * This function implements a K-means clustering algorithm for segmenting image data into
 * a specified number of clusters. The algorithm iteratively refines the cluster centroids
 * based on the image intensity values. The clustering can be constrained by optional label
 * and mask arrays.
 *
 * @src: Pointer to the float array containing the source image data.
 *       The array should have 'vol' elements, where 'vol' is the product of the dimensions.
 *
 * @label: Pointer to an unsigned char array used as a label map. Can be NULL, in which case
 *         it is ignored. If provided, it should have 'vol' elements.
 *
 * @mask: Pointer to an unsigned char array used as a mask. Can be NULL, in which case all
 *        elements are considered. If provided, only elements where the mask has a non-zero
 *        value are considered in the clustering. Should have 'vol' elements.
 *
 * @NI: Integer specifying the number of iterations for the K-means algorithm.
 *
 * @n_clusters: Integer specifying the number of clusters to segment the data into.
 *
 * @mean: Pointer to a double array where the calculated mean values of the clusters will be stored.
 *        This array should have at least 'n_clusters' elements.
 *
 * @voxelsize: Pointer to a double array indicating the size of each voxel in the image data.
 *
 * @dims: Pointer to an integer array of size 3, indicating the dimensions of the image volume.
 *
 * @thresh_mask: Integer threshold value for the mask. Elements in 'src' are considered for clustering
 *               only if the corresponding mask value is greater than or equal to this threshold.
 *
 * @thresh_kmeans: Integer threshold value for K-means clustering. It is used within the algorithm
 *                 to determine the inclusion of elements in clusters.
 *
 * The function starts by finding the maximum intensity within the masked region of the image and
 * then performs K-means clustering for increasing numbers of clusters up to 'n_clusters'. For each
 * number of clusters, it refines the centroids and updates the mean values. The final cluster means
 * are stored in the 'mean' array. The function returns the maximum source intensity found.
 *
 * Return: The maximum intensity value found in the source image data.
 */
double Kmeans(float *src, unsigned char *label, unsigned char *mask, int NI, int n_clusters, double *mean, double *voxelsize, int *dims, int thresh_mask, int thresh_kmeans)
{
    int i, j, k;
    double e, emin, eps, val, sum;
    double max_src = -HUGE;
    long n[n_clusters];
    double var[n_clusters],  mu[n_clusters],  Mu[n_clusters];
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

    for (i = 0; i < n_classes; i++) mean[i] = mu[i];       
    
    return(max_src);        
}













