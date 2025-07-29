/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 *
 */

#include <bicpl/marching.h>
#include "CAT_Vol.h"
#include "genus0.h"
#include "CAT_Separate.h"
#include "CAT_Surf.h"
#include "CAT_Smooth.h"
#include "CAT_Curvature.h"
#include "CAT_Intersect.h"
#include "CAT_MarchingCubes.h"

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

    if (!vol_euler || !vol_euler_orig || !vol_bin) {
        fprintf(stderr, "Memory allocation error\n");
        exit(EXIT_FAILURE);
    }


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
        for (i = 0; i < nvol; i++) {
            volume[i] = (vol_euler[i] < vol_euler_orig[i]) ? mn : volume[i];
            volume[i] = (vol_euler[i] > vol_euler_orig[i]) ? mx : volume[i];
        }
            
        if (n_errors == 0) break;
    }
    
    free(vol_euler);
    free(vol_bin);
    free(vol_euler_orig);

}

private void
clear_slice(
    double **slice,
    int sizes[3])
{
    int    x, y;

    for (x = 0; x < sizes[0]; x++)
    for (y = 0; y < sizes[1]; y++)
    {
        slice[x][y] = 0.0;
    }
}

private void 
input_slice(
    float *vol,
    double **slice,
    int    z,
    int sizes[3])
{
    int    x, y;

    for (x = 0; x < sizes[0]; x++)
    for (y = 0; y < sizes[1]; y++)
    {
        slice[x][y] = (double)vol[sub2ind(x, y, z, sizes)];
    }
}

private double
get_slice_value(
    double ***slices,
    int    x_size,
    int    y_size,
    int    z,
    int    x,
    int    y)
{
    if (x < 0 || x >= x_size || y < 0 || y >= y_size)
        return(0.0);
    else
        return(slices[z][x][y]);
}

private void
clear_points(
    int x_size,
    int y_size,
    int max_edges,
    int ***point_ids)
{
    int x, y, edge;

    for (x = 0; x < x_size+2; x++)
    for (y = 0; y < y_size+2; y++)
    for (edge = 0; edge < max_edges; edge++)
    {
        point_ids[x][y][edge] = -1;
    }
}

private void
get_world_point(
    double slice,
    double x,
    double y,
    mat44 nii_mat,
    Point  *point)
{
    double xw, yw, zw;
    double real_voxel[N_DIMENSIONS];

    real_voxel[0] = slice;
    real_voxel[1] = x;
    real_voxel[2] = y;

    xw = real_voxel[Y] * nii_mat.m[0][0] + 
         real_voxel[Z] * nii_mat.m[0][1] + 
         real_voxel[X] * nii_mat.m[0][2] + nii_mat.m[0][3];
    yw = real_voxel[Y] * nii_mat.m[1][0] + 
         real_voxel[Z] * nii_mat.m[1][1] + 
         real_voxel[X] * nii_mat.m[1][2] + nii_mat.m[1][3];
    zw = real_voxel[Y] * nii_mat.m[2][0] + 
         real_voxel[Z] * nii_mat.m[2][1] + 
         real_voxel[X] * nii_mat.m[2][2] + nii_mat.m[2][3];

    fill_Point(*point, xw, yw, zw);
}

private int
get_point_index(
    int x,
    int y,
    int slice_index,
    int x_size,
    int y_size,
    voxel_point_type *point,
    double corners[2][2][2],
    mat44  nii_mat,
    BOOLEAN  binary_flag,
    double   min_threshold,
    double   max_threshold,
    int    ***point_ids[],
    polygons_struct *polygons)
{
    int    voxel[N_DIMENSIONS], edge, point_index;
    int    edge_voxel[N_DIMENSIONS];
    double v[N_DIMENSIONS];
    Point  world_point;
    Point_classes point_class;

    voxel[X] = x + point->coord[X];
    voxel[Y] = y + point->coord[Y];
    voxel[Z] = point->coord[Z];
    edge = point->edge_intersected;

    point_index = point_ids[voxel[Z]][voxel[X]+1][voxel[Y]+1][edge];
    if (point_index < 0)
    {
        edge_voxel[X] = point->coord[X];
        edge_voxel[Y] = point->coord[Y];
        edge_voxel[Z] = point->coord[Z];
        point_class = get_isosurface_point(corners, edge_voxel, edge,
                        binary_flag,
                        min_threshold, max_threshold, v);

        get_world_point(v[Z] + (float) slice_index,
                        v[X] + (float) x, v[Y] + (float) y,
                        nii_mat, &world_point);

        point_index = polygons->n_points;
        ADD_ELEMENT_TO_ARRAY(polygons->points, polygons->n_points,
                        world_point, CHUNK_SIZE);

        point_ids[voxel[Z]][voxel[X]+1][voxel[Y]+1][edge] = point_index;
    }

    return(point_index);
}

private void
extract_surface(
    Marching_cubes_methods    method,
    BOOLEAN           binary_flag,
    double            min_threshold,
    double            max_threshold,
    double            valid_low,
    double            valid_high,
    int               x_size,
    int               y_size,
    double            ***slices,
    double            min_label,
    double            max_label,
    double            ***label_slices,
    int               slice_index,
    BOOLEAN           right_handed,
    mat44             nii_mat,
    int               ***point_ids[],
    polygons_struct   *polygons)
{
    int         x, y, *sizes, tx, ty, tz, n_polys, ind;
    int         p, point_index, poly, size, start_points, dir;
    voxel_point_type  *points;
    double        corners[2][2][2], label;
    BOOLEAN       valid;

    for (x = -1; x < x_size; x++)
    for (y = -1; y < y_size; y++) {
        valid = TRUE;
        for (tx = 0; tx < 2; tx++)
        for (ty = 0; ty < 2; ty++)
        for (tz = 0; tz < 2; tz++)
        {
            corners[tx][ty][tz] = get_slice_value(slices, x_size, y_size,
                tz, x + tx, y + ty);
            if (valid_low <= valid_high &&
                   (corners[tx][ty][tz] < min_threshold ||
                  corners[tx][ty][tz] > max_threshold) &&
                   (corners[tx][ty][tz] < valid_low ||
                  corners[tx][ty][tz] > valid_high))
                valid = FALSE;

            if (min_label <= max_label)
            {
                label = get_slice_value(label_slices, x_size, y_size,
                          tz, x + tx, y + ty);
                if (label < min_label || label > max_label)
                    corners[tx][ty][tz] = 0.0;
            }
        }

        if (!valid)
            continue;

        n_polys = compute_isosurface_in_voxel(method, x, y, slice_index,
                          corners, binary_flag, min_threshold,
                          max_threshold, &sizes, &points);

        if (n_polys == 0)
            continue;

        if (right_handed)
        {
            start_points = 0;
            dir = 1;
        }
        else
        {
            start_points = sizes[0]-1;
            dir = -1;
        }

        for (poly = 0; poly < n_polys; poly++)
        {
            size = sizes[poly];

            start_new_polygon(polygons);

            /*--- orient polygons properly */

            for (p = 0; p < size; p++)
            {
                ind = start_points + p * dir;
                point_index = get_point_index(x, y, slice_index,
                          x_size, y_size, &points[ind], corners,
                          nii_mat,
                          binary_flag, min_threshold, max_threshold,
                          point_ids, polygons);

                ADD_ELEMENT_TO_ARRAY(polygons->indices,
                          polygons->end_indices[polygons->n_items-1],
                          point_index, CHUNK_SIZE);
            }

            if (right_handed)
                start_points += size;
            else if (poly < n_polys-1)
                start_points += sizes[poly+1];
        }
    }
}

void
extract_isosurface(
    float  *vol,
    int sizes[3],
    double  min_label,
    double  max_label,
    mat44 nii_mat,
    Marching_cubes_methods    method,
    BOOLEAN binary_flag,
    double  min_threshold,
    double  max_threshold,
    double  valid_low,
    double  valid_high,
    polygons_struct *polygons,
    int verbose)
{
    int   n_slices, x_size, y_size, slice;
    int   ***point_ids[2], ***tmp_point_ids;
    int   max_edges;
    double **slices[2], **tmp_slices;
    double **label_slices[2];
    progress_struct progress;
    Surfprop spr;
    Point  point000, point100, point010, point001;
    Vector   v100, v010, v001, perp;
    BOOLEAN  right_handed;

    get_world_point(0.0, 0.0, 0.0, nii_mat, &point000);
    get_world_point(1.0, 0.0, 0.0, nii_mat, &point100);
    get_world_point(0.0, 1.0, 0.0, nii_mat, &point010);
    get_world_point(0.0, 0.0, 1.0, nii_mat, &point001);

    SUB_POINTS(v100, point100, point000);
    SUB_POINTS(v010, point010, point000);
    SUB_POINTS(v001, point001, point000);
    CROSS_VECTORS(perp, v100, v010);

    right_handed = DOT_VECTORS(perp, v001) >= 0.0;

    x_size = sizes[X];
    y_size = sizes[Y];
    n_slices = sizes[Z];

    ALLOC2D(slices[0], x_size, y_size);
    ALLOC2D(slices[1], x_size, y_size);

    max_edges = get_max_marching_edges(method);

    ALLOC3D(point_ids[0], x_size+2, y_size+2, max_edges);
    ALLOC3D(point_ids[1], x_size+2, y_size+2, max_edges);

    clear_slice(slices[1], sizes);

    clear_points(x_size, y_size, max_edges, point_ids[0]);
    clear_points(x_size, y_size, max_edges, point_ids[1]);

    Surfprop_a(spr) = 0.3f;
    Surfprop_d(spr) = 0.6f;
    Surfprop_s(spr) = 0.6f;
    Surfprop_se(spr)= 30.0f;
    Surfprop_t(spr) = 1.0f;
    initialize_polygons(polygons, WHITE, &spr);

    if (verbose) initialize_progress_report(&progress, FALSE, n_slices+1, "Extracting Surface");

    for (slice = 0; slice < n_slices; slice++)
    {
        tmp_slices = slices[0];
        slices[0] = slices[1];
        slices[1] = tmp_slices;
        if (slice < n_slices - 1)
            input_slice(vol, slices[1], slice, sizes);
        else
            clear_slice(slices[1], sizes);

        tmp_point_ids = point_ids[0];
        point_ids[0] = point_ids[1];
        point_ids[1] = tmp_point_ids;
        clear_points(x_size, y_size, max_edges, point_ids[1]);

        extract_surface(method, binary_flag, min_threshold, max_threshold,
            valid_low, valid_high,
            x_size, y_size, slices,
            min_label, max_label, label_slices, slice - 1,
            right_handed, nii_mat,
            point_ids, polygons);

        if (verbose) update_progress_report(&progress, slice+2);
    }

    if (verbose) terminate_progress_report(&progress);

    if (polygons->n_points > 0)
    {
        ALLOC(polygons->normals, polygons->n_points);
        compute_polygon_normals(polygons);
    }

    FREE2D(slices[0]);
    FREE2D(slices[1]);

    FREE3D(point_ids[0]);
    FREE3D(point_ids[1]);
}

/* Function to apply marching cubes and extract polygons */
object_struct *apply_marching_cubes(float *input_float, nifti_image *nii_ptr,
                        double min_threshold, double pre_fwhm, double post_fwhm,
                        double dist_morph, int n_median_filter, int n_iter, int verbose) 
{
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
    for (i = 0; i < nvol; i++) mask[i] = input_float[i] != 0;
    
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

        EC = euler_characteristic(polygons, verbose);
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

    /* Replace input values by change map */
    for (i = 0; i < nvol; i++)
        input_float[i] = vol_changed[i];

    free(vol_changed);
    free(vol_uint8);
    free(vol_uint16);
    free(vol_float);
    
    return object;
}

