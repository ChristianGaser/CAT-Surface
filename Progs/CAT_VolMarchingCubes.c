/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 *
 */

#include <ParseArgv.h>
#include "CAT_Separate.h"
#include "CAT_NiftiLib.h"
#include "CAT_SurfaceIO.h"
#include "CAT_Surf.h"
#include "CAT_Smooth.h"
#include "CAT_Curvature.h"
#include "CAT_Vol.h"
#include "CAT_Intersect.h"
#include "CAT_MarchingCubes.h"
#include "genus0.h"

/* argument defaults */
double local_smoothing = 10.0;
double min_threshold = 0.5;
double post_fwhm = 2.0;
double pre_fwhm = 2.0;
double dist_morph = FLT_MAX;
int n_median_filter = 2;
int verbose = 0;
int n_iter = 10;

/* the argument table */
static ArgvInfo argTable[] = {
  {"-thresh", ARGV_FLOAT, (char *) TRUE, (char *) &min_threshold,
    "Define the volume threshold, also known as the isovalue.\n\
     This value is crucial for initial image thresholding."},
  
  {"-pre-fwhm", ARGV_FLOAT, (char *) TRUE, (char *) &pre_fwhm,
    "Specify the Full Width Half Maximum (FWHM) for the preprocessing\n\
     smoothing filter. This helps in preserving gyri and sulci by\n\
     creating a weighted average between original and smoothed\n\
     images based on the gradient of the input image. Areas with \n\
     topology artefacts are often characterized by large gradients,\n\
     thus smoothing in these areas tries to prevent these artefacts."},

  {"-post-fwhm", ARGV_FLOAT, (char *) TRUE, (char *) &post_fwhm,
    "Set FWHM for surface smoothing. This aids in correcting the mesh\n\
     in folded areas like gyri and sulci. Note: Do not use smoothing\n\
     size > 3 mm for reliable compensation in these areas."},
  
  {"-dist-morph", ARGV_FLOAT, (char *) TRUE, (char *) &dist_morph,
    "Apply initial morphological open or close step. Close is used\n\
     by a value around 1.0 and open by negative values around -1.0.\n\
     The default automatically estimates the optimal value"},
  
  {"-median-filter", ARGV_INT, (char *) TRUE, (char *) &n_median_filter,
    "Specify the number of iterations to apply a median filter to areas\n\
     where the gradient of the thresholded image indicates larger clusters.\n\
     These clusters may point to potential topology artifacts and regions\n\
     with high local variations. This process helps to smooth these areas, \n\
     improving the quality of the surface reconstruction in subsequent steps."},
  
  {"-iter", ARGV_INT, (char *) TRUE, (char *) &n_iter,
    "Number of iterations."},
  
  {"-local-smoothing", ARGV_FLOAT, (char *) TRUE, (char *) &local_smoothing,
    "Apply local surface smoothing to resulting surface in areas where the distance\n\
     between the surface and a shifted surface is below the expected distance,\n\
     which often happens due to self intersections of the surface."},
  
  {"-verbose", ARGV_CONSTANT, (char *) TRUE, (char *) &verbose,
    "Enable verbose mode for detailed output during processing."},

    {NULL, ARGV_END, NULL, NULL, NULL}
};


private void
usage(
    char *executable)
{
    char *usage_str = "\n\
Usage: CAT_VolMarchingCubes input.nii output_surface_file [change_map.nii thresholded_map.nii]\n\
\n\
    This method generates a mesh with an Euler number of 2 (genus 0) from the\n\
    thresholded volume. The process involves:\n\
    \n\
    1. **Preprocessing with Smoothing Filter:**\n\
       - Apply a smoothing filter to the input image to remove outliers.\n\
       - Use a weighted average of the original and smoothed images to\n\
         preserve gyri and sulci.\n\
       - Weighting is based on the gradient of the input image.\n\
       - Weights range from 0 (areas with low gradient) to 1 (areas\n\
         with large gradient), with intermediate values scaled linearly.\n\
       - Weighting effect is enhanced by squaring the value.\n\
    \n\
    2. **Preprocessing with Median Filter:**\n\
       - Apply an iterative median filter to areas where the gradient of \n\
         the thresholded image indicates larger clusters.\n\
       - Use a weighted average of the original and median filterd images.\n\
       - Weighting is estimated using gradient of the input image and.\n\
         morphological operations to find larger clusters\n\
    \n\
    3. **Morphological Opening:**\n\
       - Apply additional morphological opening, scaled by `scl_open`,\n\
         to prevent gyri fusion and minimize local artifacts.\n\
       - Opening strength is determined by analyzing the impact of\n\
         different `scl_open` values and tracking RMSE changes.\n\
    \n\
    4. **Extraction of the Largest Component:**\n\
       - Extract the largest component for further processing.\n\
    \n\
    5. **Mesh Smoothing:**\n\
       - Smooth the extracted mesh.\n\
    \n\
    6. **Mesh Correction in Folded Areas:**\n\
       - Correct the mesh in areas with folds, particularly in gyri and\n\
         sulci, to counterbalance the averaging effect from smoothing.\n\
       - Use mean curvature average as a folding measure to estimate\n\
         necessary compensation.\n\
       - Compensation degree is auto-calculated based on deviation\n\
         from the defined isovalue.\n\
    7. **Mesh Correction in Areas with Self Intersections:**\n\
       - Apply local surface smoothing to resulting surface in areas where\n\
         the distance between the surface and a shifted surface is below \n\
         the expected distance, which often happens due to self intersections\n\
         of the surface.\n\n";

    fprintf(stderr,"%s\n %s\n",usage_str, executable);
}

#define IDX(x, y, z, nx, ny) ((z) * (nx) * (ny) + (y) * (nx) + (x))

void correct_topology(float *volume, float thresh, int dims[3], int conn_arr[2], int n_loops) {
    int loop, i, x, y, z, iter, val_corr;
    int V, E, F, conn, chi_local, n_errors, cube[8];
    int edges18[18][2] = {
        {0,1}, {1,3}, {3,2}, {2,0}, // Bottom face edges
        {4,5}, {5,7}, {7,6}, {6,4}, // Top face edges
        {0,4}, {1,5}, {2,6}, {3,7}, // Vertical edges
        {0,5}, {1,4}, {2,7}, {3,6}, // Diagonal edges across faces
        {0,7}, {3,4}  // Extended diagonal connections
    };
    int edges26[26][2] = {
        {0,1}, {1,3}, {3,2}, {2,0}, // Bottom face edges
        {4,5}, {5,7}, {7,6}, {6,4}, // Top face edges
        {0,4}, {1,5}, {2,6}, {3,7}, // Vertical edges
        {0,5}, {1,4}, {2,7}, {3,6}, // Diagonal edges across faces
        {0,7}, {3,4}, {2,5}, {1,6}, // Extended diagonal connections
        {0,6}, {1,7}, {2,4}, {3,5}  // Corner-to-corner vertex edges
    };
    // Count faces (6 faces in a cube)
    int faces[6][4] = {
        {0,1,3,2}, {4,5,7,6}, // XY faces
        {0,2,6,4}, {1,3,7,5}, // YZ faces
        {0,1,5,4}, {2,3,7,6}  // XZ faces
    };
    
    int nx = dims[0], ny = dims[1], nz = dims[2];
    int nvol = nx*ny*nz;
    
    double mn = get_min(volume, nvol, 0, DT_FLOAT32);
    double mx = get_max(volume, nvol, 0, DT_FLOAT32);

    float *vol_euler = (float *)malloc(sizeof(float)*nvol);
    float *vol_euler_orig = (float *)malloc(sizeof(float)*nvol);
    unsigned short *vol_bin = (unsigned short *)malloc(sizeof(unsigned short)*nvol);

    // Initialize Euler map
    for (i = 0; i < nvol; i++) {
        vol_euler[i] = (volume[i] >= thresh) ? 1.0 : 0.0;
        vol_euler_orig[i] = vol_euler[i];
    }

    for (iter = 0; iter < 50; iter++) {
        n_errors = 0;
        for (loop = 0; loop < n_loops; loop++) {
    
            // Threshold volume
            for (i = 0; i < nvol; i++)
                vol_bin[i] = (volume[i] > thresh) ? 1.0 - (float)loop : 0.0 + (float)loop;
                
            conn = conn_arr[((loop + iter) % 2)];
            val_corr = (loop == 0) ? 0.0 : 1.0;

            for (z = 0; z < nz - 1; z++) {
                for (y = 0; y < ny - 1; y++) {
                    for (x = 0; x < nx - 1; x++) {
                            
                        // Extract 2x2x2 voxel cube
                        cube[0] = vol_bin[IDX(x, y, z, nx, ny)];
                        cube[1] = vol_bin[IDX(x+1, y, z, nx, ny)];
                        cube[2] = vol_bin[IDX(x, y+1, z, nx, ny)];
                        cube[3] = vol_bin[IDX(x+1, y+1, z, nx, ny)];
                        cube[4] = vol_bin[IDX(x, y, z+1, nx, ny)];
                        cube[5] = vol_bin[IDX(x+1, y, z+1, nx, ny)];
                        cube[6] = vol_bin[IDX(x, y+1, z+1, nx, ny)];
                        cube[7] = vol_bin[IDX(x+1, y+1, z+1, nx, ny)];
                        
                        // Count vertices, edges, faces
                        V = 0; E = 0; F = 0;
        
                        // Count vertices (fully occupied corners)
                        for (i = 0; i < 8; i++) V += cube[i];
                        
                        if (conn == 18) {
                            for (i = 0; i < 18; i++) 
                                if (cube[edges18[i][0]] && cube[edges18[i][1]]) E++;
                        } else if (conn == 26) {
                            for (i = 0; i < 26; i++) 
                                if (cube[edges26[i][0]] && cube[edges26[i][1]]) E++;
                        }
        
                        for (i = 0; i < 6; i++)
                            if (cube[faces[i][0]] && cube[faces[i][1]] && cube[faces[i][2]] && cube[faces[i][3]]) F++;
        
                        // Compute local Euler number
                        chi_local = V - E + F;
        
                        if (vol_bin[IDX(x, y, z, nx, ny)] > 0) {
                            if (((conn == 18) &&  (chi_local ==  2)) || 
                                ((conn == 26) &&  (chi_local == -6))) {
                                vol_euler[IDX(x,   y,   z,   nx, ny)] = val_corr;
                                vol_euler[IDX(x+1, y,   z,   nx, ny)] = val_corr;
                                vol_euler[IDX(x,   y+1, z,   nx, ny)] = val_corr;
                                vol_euler[IDX(x+1, y+1, z,   nx, ny)] = val_corr;
                                vol_euler[IDX(x,   y,   z+1, nx, ny)] = val_corr;
                                vol_euler[IDX(x+1, y,   z+1, nx, ny)] = val_corr;
                                vol_euler[IDX(x,   y+1, z+1, nx, ny)] = val_corr;
                                vol_euler[IDX(x+1, y+1, z+1, nx, ny)] = val_corr;
                                n_errors++;
                            }
                        }
                    }
                }
            }
        }
        
        /* Apply changes to volume */
        for (i = 0; i < nvol; i++)
            volume[i] = (vol_euler[i] < vol_euler_orig[i]) ? mn : volume[i];
            volume[i] = (vol_euler[i] > vol_euler_orig[i]) ? mx : volume[i];
            
        if (n_errors == 0) break;
    }
    
    free(vol_euler);
    free(vol_bin);
    free(vol_euler_orig);

}

/* Function to apply marching cubes and extract polygons */
object_struct *apply_marching_cubes(float *input_float, nifti_image *nii_ptr) {
    double voxelsize[N_DIMENSIONS];
    double best_dist;
    int dims[MAX_DIMENSIONS], i, k, nvol, count_change = 0, best_change_values;
    object_struct **object2;
    Marching_cubes_methods method = (Marching_cubes_methods)1;

    mat44 nii_mat = nii_ptr->sto_xyz;
    dims[0] = nii_ptr->nx;
    dims[1] = nii_ptr->ny;
    dims[2] = nii_ptr->nz;
    voxelsize[0] = nii_ptr->dx;
    voxelsize[1] = nii_ptr->dy;
    voxelsize[2] = nii_ptr->dz;

    nvol = dims[0] * dims[1] * dims[2];

    /* Memory allocation */
    float *vol_changed = (float *)calloc(nvol, sizeof(float));
    float *vol_float = (float *)malloc(nvol * sizeof(float));
    unsigned short *vol_uint16 = (unsigned short *)malloc(nvol * sizeof(unsigned short));
    unsigned char *vol_uint8 = (unsigned char *)malloc(nvol * sizeof(unsigned char));

    if (!vol_float || !vol_uint16 || !vol_uint8) {
        fprintf(stderr, "Memory allocation error\n");
        exit(EXIT_FAILURE);
    }

    /* Normalize input if necessary */
    if (min_threshold > 1.0) {
        double max_value = get_max(input_float, nvol, 0, DT_FLOAT32);
        for (i = 0; i < nvol; i++) input_float[i] /= (float)max_value;
        min_threshold /= max_value;
    }

    /* Preprocessing Step: Smoothing Filter
       - Purpose: To remove outliers in the input image.
       - Method: A weighted average is calculated between the original image and the 
         smoothed image to protect the structural integrity of gyri and sulci.
       - Weight Estimation: Based on the gradient of the input image.
         For values below mean gradient, the weight is 0.
       - Enhanced Weighting: Weights are squared to emphasize larger weightings,
         providing a more robust distinction between regions of interest and outliers.
    */
    if (pre_fwhm) {
        float *grad = (float *)malloc(nvol * sizeof(float));
        if (!grad) {
            fprintf(stderr, "Memory allocation error\n");
            exit(EXIT_FAILURE);
        }
        
        for (i = 0; i < nvol; i++)
            vol_float[i] = input_float[i];

        double s[] = {pre_fwhm, pre_fwhm, pre_fwhm};
        smooth3(vol_float, dims, voxelsize, s, 0, DT_FLOAT32);

        gradient3D(input_float, grad, NULL, NULL, NULL, dims, voxelsize);
        double mean_grad = get_mean(grad, nvol, 1, DT_FLOAT32);
        double max_grad = get_max(grad, nvol, 1, DT_FLOAT32);

        for (i = 0; i < nvol; i++) {
            float weight = (grad[i] - mean_grad) / (max_grad - mean_grad);
            weight = (weight < 0.0) ? 0.0f : weight;
            weight = (weight > 1.0) ? 1.0f : weight;
            weight *= weight;

            input_float[i] = (1.0 - weight) * input_float[i] + weight * vol_float[i];
        }
        free(grad);
    }
    
    /* Apply iterative median filter to strengthen structures */
    unsigned char *mask = (unsigned char *)malloc(nvol * sizeof(unsigned char));
    for (i = 0; i < nvol; i++)
        mask[i] = input_float[i] != 0;
    
    median3(input_float, mask, dims, 3, DT_FLOAT32);
    free(mask);
            
    /* Keep largest cluster and fill holes */
    keep_largest_cluster(input_float, min_threshold, dims, DT_FLOAT32, 0, 1, 18);
    fill_holes(input_float, dims, min_threshold, -1.0, DT_FLOAT32);

    /* Correct topology */
    int conn_arr[2] = {18, 26};
    correct_topology(input_float, min_threshold, dims, conn_arr, 2);

    /* Apply genus-0 correction */
    genus0parameters g0[1];
    genus0init(g0);
    for (i = 0; i < N_DIMENSIONS; i++) g0->dims[i] = dims[i];

    g0->input = vol_uint16;
    g0->connected_component = 1;
    g0->value = 1;
    g0->contour_value = 1;
    g0->any_genus = 0;
    g0->biggest_component = 1;
    g0->pad[0] = g0->pad[1] = g0->pad[2] = 2;
    g0->ijk2ras = NULL;
    g0->verbose = 0;
    g0->return_surface = 0;
    g0->extraijkscale[2] = 1;

    /* Iterative genus-0 correction */
    int count = 0, EC = -1;

    object_struct *object = create_object(POLYGONS);
    polygons_struct *polygons = get_polygons_ptr(object);
    
    /* Find optimal dist-parameter for open/close to minimize the number of voxel that
       have to be changed during topology correction */
    if (dist_morph == FLT_MAX) {
        double dist_values[] = {-1.0, -0.5, 0.0, 0.5, 1.0, 1.5};
        int change_values[] = {0, 0, 0, 0, 0, 0};
        int n_values = sizeof(change_values) / sizeof(change_values[0]);
        for (k = 0; k < n_values; k++) {
    
            for (i = 0; i < nvol; i++)
                vol_uint16[i] = (input_float[i] >= min_threshold) ? 1 : 0;
    
            if (dist_values[k] > 0.0)
                distclose(vol_uint16, dims, voxelsize, dist_values[k], 0.0, DT_UINT16);
            else if (dist_values[k] < 0.0)
                distopen(vol_uint16, dims, voxelsize, -dist_values[k], 0.0, DT_UINT16);
    
            /* call genus0 for the 1st time */
            g0->input = vol_uint16;
            g0->cut_loops = 0;
            g0->connectivity = 6;
            g0->alt_value = 1;
            g0->alt_contour_value = 1;
            if (genus0(g0)) return(NULL); /* check for error */
        
            /* save results as next input */
            for (i = 0; i < nvol; i++)
                vol_uint16[i] = g0->output[i];
    
            /* save changes */
            for (i = 0; i < nvol; i++)
                change_values[k] += (int)fabs((float)g0->output[i] - (float)g0->input[i]);   
    
            /* call genus0 a 2nd time with other parameters */
            g0->input = vol_uint16;
            g0->cut_loops = 1;
            g0->connectivity = 18;
            g0->alt_value = 0;
            g0->alt_contour_value = 0;    
            if (genus0(g0)) return(NULL); 
    
            /* save changes */
            for (i = 0; i < nvol; i++)
                change_values[k] += (int)fabs((float)g0->output[i] - (float)g0->input[i]);
            if (change_values[k] == 0) {
                n_values = k+1;
                break;
            }
        }
        int ind_min_value;
        int min_value = 1E9;
        for (k = 0; k < n_values; k++) {
            if (change_values[k] < min_value) {
                min_value = change_values[k];
                ind_min_value = k;
            }
        }
        best_dist = dist_values[ind_min_value];
        best_change_values = change_values[ind_min_value];
    } else {
        best_dist = dist_morph;
        best_change_values = 1E3;
    }    
    if (verbose) printf("Optimal dist-parameter for morphological operations: %f\n", best_dist);

    /* Convert to binary image */
    for (i = 0; i < nvol; i++)
        vol_uint16[i] = (input_float[i] >= min_threshold) ? 1 : 0;
    
    while ((EC != 2) && (count < n_iter)) {
        
        /* Only move on if topology correction is still necessary */
        if (best_change_values > 0) {
            if (best_dist > 0.0)
                distclose(vol_uint16, dims, voxelsize, best_dist, 0.0, DT_UINT16);
            else if (best_dist < 0.0)
                distopen(vol_uint16, dims, voxelsize, -best_dist, 0.0, DT_UINT16);
    
            /* call genus0 for the 1st time */
            g0->input = vol_uint16;
            g0->cut_loops = 0;
            g0->connectivity = 6;
            g0->alt_value = 1;
            g0->alt_contour_value = 1;
            if (genus0(g0)) return(NULL); /* check for error */
        
            /* save results as next input */
            for (i = 0; i < nvol; i++)
                vol_uint16[i] = g0->output[i];
    
            /* save changes */
            for (i = 0; i < nvol; i++)
                vol_changed[i] += (float)(count + 1)*((float)g0->output[i] - (float)g0->input[i]);   
    
            /* call genus0 a 2nd time with other parameters */
            g0->input = vol_uint16;
            g0->cut_loops = 1;
            g0->connectivity = 18;
            g0->alt_value = 0;
            g0->alt_contour_value = 0;    
            if (genus0(g0)) return(NULL); 
    
            /* save changes */
            count_change = 0;
            for (i = 0; i < nvol; i++) {
                vol_changed[i] += (float)(count + 1)*((float)g0->output[i] - (float)g0->input[i]);
                if (vol_changed[i] != 0.0) count_change++;
            }
                
            if (n_median_filter) {
                /* find areas that were corrected for topology artefacts and dilate them */
                for (i = 0; i < nvol; i++)
                    vol_uint8[i] = (unsigned char)(vol_uint16[i] != g0->output[i]);
                    
                morph_dilate(vol_uint8, dims, 4, 0.5, DT_UINT8);
    
                /* use previous output for filtering */
                for (i = 0; i < nvol; i++)
                    vol_uint16[i] = g0->output[i];
                
                /* Apply iterative median filter */
                median3(vol_uint16, NULL, dims, n_median_filter, DT_UINT16);
    
                /* replace genus0 output with its median filtered version in (dilated)
                   areas with topology artefacts */
                for (i = 0; i < nvol; i++)
                    g0->output[i] = (unsigned short)(((vol_uint8[i] > 0)) ? vol_uint16[i] : g0->output[i]);
            }
        }
        
        keep_largest_cluster(g0->output, min_threshold, dims, DT_UINT16, 0, 1, 18);
        fill_holes(g0->output, dims, min_threshold, -1.0, DT_UINT16);

        for (i = 0; i < nvol; i++)
            vol_float[i] = (float)g0->output[i];

        extract_isosurface(vol_float, dims, 0.0, -1.0, nii_mat, method, FALSE,
                           min_threshold, min_threshold, 0.0, -1.0, polygons, verbose);

        compute_polygon_normals(polygons);
        check_polygons_neighbours_computed(polygons);
        int n_out = separate_polygons(polygons, -1, &object2);          
        triangulate_polygons(get_polygons_ptr(object2[0]), polygons);

        EC = euler_characteristic(polygons);
        count++;
        if (verbose) printf("Euler characteristics after %d iterations is %d (%d voxel changed).\n", count, EC, count_change);
  
        /* save results as next input */
        for (i = 0; i < nvol; i++)
            vol_uint16[i] = g0->output[i];
    }

    /* Mesh Correction in Folded Areas
       - Objective: To compensate for the averaging effect observed in gyri and sulci.
       - Method: Utilization of a folding measure, specifically the mean curvature average, 
         to estimate the necessary compensation in these areas.
       - Compensation Estimation: Automatically calculated based on the difference 
         between the actual mesh curvature and the predefined isovalue. This approach 
         ensures accurate correction in folded regions, maintaining the integrity of 
         gyri and sulci structures.
    */
    if (post_fwhm > 0.0) {
        smooth_heatkernel(polygons, NULL, post_fwhm);
        correct_mesh_folding(polygons, NULL, input_float, nii_ptr, min_threshold);
    }

    /* Mesh Correction in Areas with Self Intersections
       - Objective: To correct areas with self intersections
       - Method: Apply local surface smoothing to resulting surface in areas 
         where the distance between the surface and a shifted surface is below 
         the expected distance, which often happens due to self intersections of 
         the surface.
    */
    if (local_smoothing > 0.0) {
        double *values = (double *)malloc(sizeof(double) * polygons->n_points);
        double *extents = (double *)malloc(sizeof(double) * polygons->n_points);
        polygons_struct *smooth_polygons = (polygons_struct *)malloc(sizeof(polygons_struct));
        copy_polygons(polygons, smooth_polygons);

        for (i = 0; i < polygons->n_points; i++) {
            extents[i] = 0.1;
            values[i] = 3.0;
        }

        object2 = central_to_new_pial(polygons, values, extents, NULL, NULL, 0, 0);
        compute_exact_hausdorff(polygons, get_polygons_ptr(object2[0]), values, 0);
        smooth_heatkernel(polygons, values, 5.0);

        double min_value = 0.25;
        double max_value = 0.3;
        for (i = 0; i < polygons->n_points; i++) {
            /* scale values between 0.25..0.3 and force min=0 and max=1 */
            values[i] = (values[i] - min_value)/(max_value - min_value);
            values[i] = MIN(1.0, values[i]);
            values[i] = MAX(0.0, values[i]);
            values[i] = values[i]*values[i];
        }
        
        /* smooth values to obtain a smooth border */
        smooth_heatkernel(polygons, values, 5.0);
        
        /* obtain smoothed surface */
        smooth_heatkernel(smooth_polygons, NULL, local_smoothing);

        /* use smoothed or original surface w.r.t. weighting */ 
        for (i = 0; i < polygons->n_points; i++)
            Point_x(polygons->points[i]) = values[i]*Point_x(polygons->points[i]) + (1.0 - values[i])*Point_x(smooth_polygons->points[i]);

        free(smooth_polygons);
        free(values);
        free(extents);
    }

    for (i = 0; i < nvol; i++)
        input_float[i] = vol_changed[i];

    free(vol_changed);
    free(vol_uint8);
    free(vol_uint16);
    free(vol_float);
    
    return object;
}

int main(int argc, char *argv[]) {
    float *input_float;
    char out_diff[1024];

    initialize_argument_processing(argc, argv);

    /* get the arguments from the command line */
    if (ParseArgv(&argc, argv, argTable, 0)) {
        usage(argv[0]);
        fprintf(stderr, "     %s -help\n\n", argv[0]);
        return EXIT_FAILURE;
    }

    char *input_filename, *output_filename;
    if (!get_string_argument(NULL, &input_filename) || !get_string_argument(NULL, &output_filename)) {
        usage(argv[0]);
        fprintf(stderr, "Usage: CAT_VolMarchingCubes input.nii output_surface_file\n");
        return EXIT_FAILURE;
    }

    /* Read input volume */
    nifti_image *nii_ptr = read_nifti_float(input_filename, &input_float, 0);
    if (!nii_ptr) {
        fprintf(stderr, "Error reading %s.\n", input_filename);
        return EXIT_FAILURE;
    }

    object_struct *object = apply_marching_cubes(input_float, nii_ptr);
    if (object) {
        output_graphics_any_format(output_filename, ASCII_FORMAT, 1, &object, NULL);
    } else {
        fprintf(stderr, "Error generating surface.\n");
    }

    if(argc > 3) {
        double voxelsize[N_DIMENSIONS];
        int dims[MAX_DIMENSIONS];
        dims[0] = nii_ptr->nx;
        dims[1] = nii_ptr->ny;
        dims[2] = nii_ptr->nz;
        voxelsize[0] = nii_ptr->dx;
        voxelsize[1] = nii_ptr->dy;
        voxelsize[2] = nii_ptr->dz;
        (void) sprintf(out_diff, "%s", argv[3]); 
        if (!write_nifti_float(out_diff, input_float, DT_FLOAT32, 1.0, dims, voxelsize, nii_ptr)) 
            exit(EXIT_FAILURE);
    }

    free(input_float);
    delete_marching_cubes_table();
    return 0;
}
