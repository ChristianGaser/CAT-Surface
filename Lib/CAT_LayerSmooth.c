/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 *
 * Laplace-guided anisotropic smoothing along cortical layers.
 * Inspired by LAYNII's LN_LAYER_SMOOTH but uses continuous depth field
 * instead of discrete layers for smoother transitions.
 *
 */

#include <math.h>
#include <stdlib.h>
#include <memory.h>
#include <stdio.h>
#include <float.h>

#include "CAT_Vol.h"
#include "CAT_LayerSmooth.h"
#include "CAT_NiftiLib.h"

/* Tissue labels from PVE segmentation */
#define CSF_LABEL   1.0
#define CGM_LABEL   1.5   /* CSF/GM boundary */
#define GM_LABEL    2.0
#define GWM_LABEL   2.5   /* GM/WM boundary */
#define WM_LABEL    3.0

#ifndef MAX
#define MAX(A,B) (((A)>(B)) ? (A) : (B))
#endif
#ifndef MIN
#define MIN(A,B) (((A)<(B)) ? (A) : (B))
#endif
#define SQR(x) ((x)*(x))

/**
 * Gaussian weight function
 * Returns exp(-0.5 * (x/sigma)^2)
 */
static double gaussian_weight(double x, double sigma)
{
    if (sigma <= 0.0) return 1.0;
    double t = x / sigma;
    return exp(-0.5 * t * t);
}

/**
 * Convert FWHM to sigma for Gaussian
 * FWHM = 2 * sqrt(2 * ln(2)) * sigma ≈ 2.355 * sigma
 */
static double fwhm_to_sigma(double fwhm)
{
    return fwhm / (2.0 * sqrt(2.0 * log(2.0)));
}

/**
 * compute_cortical_depth - Compute cortical depth field
 *
 * Computes WM and CSF distances, then derives relative depth:
 *   depth = WMD / (WMD + CSFD)
 *
 * The 'extend' parameter allows including voxels slightly outside
 * the strict GM band (CGM < seg < GWM) by this distance in mm.
 */
void compute_cortical_depth(float *seg, float *depth, float *dist_WM_out, float *dist_CSF_out,
                            int dims[3], double voxelsize[3], double extend)
{
    int i, j, nvox;
    float *input, *dist_WM, *dist_CSF, *GMT;
    unsigned char *mask;
    double range = 0.9;  /* Range for masking similar to CAT_VolThicknessPbt */
    int n_avgs = 4;      /* Number of averages for distance estimation */
    double add_value;
    
    nvox = dims[0] * dims[1] * dims[2];
    
    /* Allocate working arrays */
    input = (float *)malloc(sizeof(float) * nvox);
    dist_WM = (float *)malloc(sizeof(float) * nvox);
    dist_CSF = (float *)malloc(sizeof(float) * nvox);
    GMT = (float *)malloc(sizeof(float) * nvox);
    mask = (unsigned char *)malloc(sizeof(unsigned char) * nvox);
    
    if (!input || !dist_WM || !dist_CSF || !mask) {
        fprintf(stderr, "Memory allocation error in compute_cortical_depth\n");
        exit(EXIT_FAILURE);
    }
    
    median3(seg, NULL, dims, 2, DT_FLOAT32);

    /* Initialize distances */
    for (i = 0; i < nvox; i++) {
        dist_CSF[i] = 0.0f;
        dist_WM[i] = 0.0f;
    }
    
    /* Compute averaged distance maps (similar to CAT_VolThicknessPbt approach) */
    for (j = 0; j < n_avgs; j++) {
        /* Estimate value for shifting the border */
        add_value = ((double)j + 1.0) / ((double)n_avgs + 1.0) - 0.5;
        
        /* Prepare map outside CSF boundary */
        for (i = 0; i < nvox; i++) {
            input[i] = (seg[i] < (CGM_LABEL + add_value)) ? 1.0f : 0.0f;
            mask[i] = (seg[i] < GWM_LABEL + range) ? 1 : 0;
        }
        
        /* Compute CSF distance */
        euclidean_distance(input, mask, dims, NULL, 0);
        for (i = 0; i < nvox; i++)
            dist_CSF[i] += input[i];
        
        /* Prepare map outside WM boundary */
        for (i = 0; i < nvox; i++) {
            input[i] = (seg[i] > (GWM_LABEL + add_value)) ? 1.0f : 0.0f;
            mask[i] = (seg[i] > CGM_LABEL - range) ? 1 : 0;
        }
        
        /* Compute WM distance */
        euclidean_distance(input, mask, dims, NULL, 0);
        for (i = 0; i < nvox; i++)
            dist_WM[i] += input[i];
    }
    
    /* Average the distances */
    for (i = 0; i < nvox; i++) {
        dist_CSF[i] /= (float)n_avgs;
        dist_WM[i] /= (float)n_avgs;
    }

    for (i = 0; i < nvox; i++) input[i] = roundf(seg[i]);

    projection_based_thickness(input, dist_WM, dist_CSF, GMT, dims, voxelsize);

    /* Use minimum/maximum to reduce issues with meninges */
    for (i = 0; i < nvox; i++) {
        float total_dist = dist_WM[i] + dist_CSF[i];
        GMT[i] = MIN(total_dist, MAX(0.0, GMT[i] - 0.125*(GMT[i]  < total_dist)));
    }

    /* Re-estimate CSF distance using corrected GM thickness */
    for (i = 0; i < nvox; i++) {
        if ((seg[i] > CGM) && (seg[i] < GWM) && (GMT[i] > 1e-15))
            dist_CSF[i] = MIN(dist_CSF[i], GMT[i] - dist_WM[i]);
    }
    clip_data(dist_CSF, nvox, 0.0, 1E15, DT_FLOAT32);

    /* Estimate percentage position map (PPM) */
    for (i = 0; i < nvox; i++) {
        if ((seg[i] > CGM) && (seg[i] < GWM) && (GMT[i] > 1e-15)) {
            depth[i] = MAX(0.0, GMT[i] - dist_WM[i]) / GMT[i];
        }
    }
    clip_data(depth, nvox, 0.0, 1.0, DT_FLOAT32);
    
    /* extend depth measure and mask it finally */
    euclidean_distance(depth, NULL, dims, NULL, 1);

    for (i = 0; i < nvox; i++)
        depth[i] = (seg[i] > (CGM_LABEL - extend)) && (seg[i] < (GWM_LABEL + extend)) ? depth[i] : 0.0f;
    
    /* Optionally copy distance maps to output */
    if (dist_WM_out) {
        memcpy(dist_WM_out, dist_WM, sizeof(float) * nvox);
    }
    if (dist_CSF_out) {
        memcpy(dist_CSF_out, dist_CSF, sizeof(float) * nvox);
    }
    
    free(input);
    free(dist_WM);
    free(dist_CSF);
    free(GMT);
    free(mask);
}

/**
 * smooth_within_cortex_float - Anisotropic smoothing along cortical layers
 *
 * For each cortical voxel, smooths using a kernel that weights neighbors by:
 * 1. Euclidean distance (standard Gaussian)
 * 2. Depth similarity (penalizes smoothing across layers)
 *
 * The combined weight is:
 *   w = exp(-d_spatial^2 / (2*sigma_spatial^2)) * exp(-d_depth^2 / (2*sigma_depth^2))
 *
 * sigma_depth controls how strictly we stay within layers. A small value
 * means very strict layer preservation; larger values allow more cross-layer smoothing.
 */
void smooth_within_cortex_float(float *data, float *seg, int dims[3], double voxelsize[3],
                                double fwhm, double extend)
{
    int x, y, z, i, ix, iy, iz, ni;
    int nx, ny, nz, nvox, nxy;
    int radius_x, radius_y, radius_z;
    double sigma_spatial, sigma_depth;
    double dx, dy, dz, dist_spatial, dist_depth;
    double weight, weight_sum, value_sum;
    float *depth, *output, *buffer;
    
    nx = dims[0];
    ny = dims[1];
    nz = dims[2];
    nvox = nx * ny * nz;
    nxy = nx * ny;
    
    /* Convert FWHM to sigma */
    sigma_spatial = fwhm_to_sigma(fwhm);
    
    /* Depth sigma controls layer strictness
     * A value of ~0.1 means depths differing by 0.1 (10% of cortical thickness)
     * will have ~60% weight reduction. Adjust this to control strictness. */
    sigma_depth = 0.15;
    
    /* Kernel radius: go out to 3 sigma (covers 99.7% of Gaussian) */
    radius_x = (int)ceil(3.0 * sigma_spatial / voxelsize[0]);
    radius_y = (int)ceil(3.0 * sigma_spatial / voxelsize[1]);
    radius_z = (int)ceil(3.0 * sigma_spatial / voxelsize[2]);
    
    /* Ensure minimum radius of 1 */
    radius_x = MAX(1, radius_x);
    radius_y = MAX(1, radius_y);
    radius_z = MAX(1, radius_z);
    
    /* Allocate arrays */
    depth = (float *)malloc(sizeof(float) * nvox);
    output = (float *)malloc(sizeof(float) * nvox);
    buffer = (float *)malloc(sizeof(float) * nvox);
    
    if (!depth || !output || !buffer) {
        fprintf(stderr, "Memory allocation error in smooth_within_cortex_float\n");
        exit(EXIT_FAILURE);
    }
    
    /* Copy input to buffer for reading */
    memcpy(buffer, data, sizeof(float) * nvox);
    
    /* Compute cortical depth field */
    compute_cortical_depth(seg, depth, NULL, NULL, dims, voxelsize, extend);

    /* Allocate normals */
    float *normals = (float *)calloc(3 * nvox, sizeof(float));
    if (!normals) {
        fprintf(stderr, "Memory allocation error for normals in smooth_within_cortex_float\n");
        exit(EXIT_FAILURE);
    }

    /* Compute gradients of depth field */
    for (z = 0; z < nz; z++) {
        for (y = 0; y < ny; y++) {
            for (x = 0; x < nx; x++) {
                i = z * nxy + y * nx + x;
                if (depth[i] <= 0.0f) continue;
                
                double grad_x = 0, grad_y = 0, grad_z = 0;
                
                if (x > 0 && x < nx - 1) grad_x = (depth[i + 1] - depth[i - 1]) / (2.0 * voxelsize[0]);
                else if (x > 0) grad_x = (depth[i] - depth[i - 1]) / voxelsize[0];
                else if (x < nx - 1) grad_x = (depth[i + 1] - depth[i]) / voxelsize[0];
                
                if (y > 0 && y < ny - 1) grad_y = (depth[i + nx] - depth[i - nx]) / (2.0 * voxelsize[1]);
                else if (y > 0) grad_y = (depth[i] - depth[i - nx]) / voxelsize[1];
                else if (y < ny - 1) grad_y = (depth[i + nx] - depth[i]) / voxelsize[1];
                
                if (z > 0 && z < nz - 1) grad_z = (depth[i + nxy] - depth[i - nxy]) / (2.0 * voxelsize[2]);
                else if (z > 0) grad_z = (depth[i] - depth[i - nxy]) / voxelsize[2];
                else if (z < nz - 1) grad_z = (depth[i + nxy] - depth[i]) / voxelsize[2];
                
                double norm = sqrt(grad_x * grad_x + grad_y * grad_y + grad_z * grad_z);
                if (norm > 1e-6) {
                    normals[3 * i] = (float)(grad_x / norm);
                    normals[3 * i + 1] = (float)(grad_y / norm);
                    normals[3 * i + 2] = (float)(grad_z / norm);
                }
            }
        }
    }
    
    /* Initialize output */
    memcpy(output, buffer, sizeof(float) * nvox);
    
    /* Main smoothing loop */
    for (z = 0; z < nz; z++) {
        for (y = 0; y < ny; y++) {
            for (x = 0; x < nx; x++) {
                i = z * nxy + y * nx + x;
                
                /* Skip voxels outside cortical band */
                if (depth[i] <= 0.0f) {
                    continue;
                }
                
                weight_sum = 0.0;
                value_sum = 0.0;
                
                /* Loop over neighborhood */
                for (iz = MAX(0, z - radius_z); iz <= MIN(nz - 1, z + radius_z); iz++) {
                    for (iy = MAX(0, y - radius_y); iy <= MIN(ny - 1, y + radius_y); iy++) {
                        for (ix = MAX(0, x - radius_x); ix <= MIN(nx - 1, x + radius_x); ix++) {
                            ni = iz * nxy + iy * nx + ix;
                            
                            /* Skip neighbors outside cortical band */
                            if (depth[ni] <= 0.0f) {
                                continue;
                            }
                            
                            /* Compute spatial distance in mm */
                            dx = (double)(ix - x) * voxelsize[0];
                            dy = (double)(iy - y) * voxelsize[1];
                            dz = (double)(iz - z) * voxelsize[2];
                            dist_spatial = sqrt(dx*dx + dy*dy + dz*dz);
                            
                            /* Compute depth difference */
                            dist_depth = fabs((double)(depth[ni] - depth[i]));
                            
                            /* Compute normal similarity (dot product) */
                            double dot = normals[3 * i] * normals[3 * ni] + 
                                         normals[3 * i + 1] * normals[3 * ni + 1] + 
                                         normals[3 * i + 2] * normals[3 * ni + 2];
                            
                            /* Penalize opposing normals (kissing sulci)
                             * Normals in kissing sulci are anti-parallel (dot approx -1)
                             * Use soft threshold */
                            double weight_angle = MAX(0.0, dot);
                            weight_angle = pow(weight_angle,6);

                            /* Combined weight: spatial * depth similarity * angle similarity */
                            weight = gaussian_weight(dist_spatial, sigma_spatial) *
                                     gaussian_weight(dist_depth, sigma_depth) *
                                     weight_angle;
                            
                            weight_sum += weight;
                            value_sum += weight * buffer[ni];
                        }
                    }
                }
                
                /* Normalize and store result */
                if (weight_sum > 1e-10) {
                    output[i] = (float)(value_sum / weight_sum);
                }
            }
        }
    }
    
    /* Copy result back to input */
    memcpy(data, output, sizeof(float) * nvox);
    
    free(normals);
    free(depth);
    free(output);
    free(buffer);
}

/**
 * smooth_within_cortex - Wrapper for any datatype
 */
void smooth_within_cortex(void *data, float *seg, int dims[3], double voxelsize[3],
                          double fwhm, double extend, int datatype)
{
    int nvox;
    float *buffer;
    
    nvox = dims[0] * dims[1] * dims[2];
    buffer = (float *)malloc(sizeof(float) * nvox);
    
    if (!buffer) {
        fprintf(stderr, "Memory allocation error in smooth_within_cortex\n");
        exit(EXIT_FAILURE);
    }
    
    /* Convert input to float */
    convert_input_type_float(data, buffer, nvox, datatype);
    
    /* Perform smoothing */
    smooth_within_cortex_float(buffer, seg, dims, voxelsize, fwhm, extend);
    
    /* Convert back to original datatype */
    convert_output_type_float(data, buffer, nvox, datatype);
    
    free(buffer);
}
