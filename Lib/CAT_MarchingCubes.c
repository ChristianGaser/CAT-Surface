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
#include "CAT_VolPbt.h"
#include "genus0.h"
#include "CAT_Separate.h"
#include "CAT_Smooth.h"
#include "CAT_Curvature.h"
#include "CAT_Intersect.h"
#include "CAT_MarchingCubes.h"

int euler_characteristic(polygons_struct *polygons, int verbose);

/*
 * Edge tables for local Euler characteristic computation in 2x2x2 cubes.
 *
 * Corner layout:  0:(x,y,z)     1:(x+1,y,z)   2:(x,y+1,z)   3:(x+1,y+1,z)
 *                 4:(x,y,z+1)   5:(x+1,y,z+1) 6:(x,y+1,z+1) 7:(x+1,y+1,z+1)
 *
 * 18-connectivity: pairs differing in 1 or 2 bit positions (face + edge adjacent)
 *   = 12 face edges + 12 face diagonals = 24 entries
 * 26-connectivity: all pairs = 18-conn + 4 space diagonals (3-bit diff) = 28 entries
 */
static const int g_edges_18conn[24][2] = {
    /* 12 face edges (1-bit difference) */
    {0,1},{1,3},{3,2},{2,0},   /* z=0 face */
    {4,5},{5,7},{7,6},{6,4},   /* z=1 face */
    {0,4},{1,5},{2,6},{3,7},   /* vertical  */
    /* 12 face diagonals (2-bit difference) */
    {0,3},{1,2},               /* z=0 face  */
    {4,7},{5,6},               /* z=1 face  */
    {0,5},{1,4},               /* y=0 face  */
    {2,7},{3,6},               /* y=1 face  */
    {0,6},{2,4},               /* x=0 face  */
    {1,7},{3,5}                /* x=1 face  */
};

static const int g_edges_26conn[28][2] = {
    /* 12 face edges */
    {0,1},{1,3},{3,2},{2,0},
    {4,5},{5,7},{7,6},{6,4},
    {0,4},{1,5},{2,6},{3,7},
    /* 12 face diagonals */
    {0,3},{1,2},
    {4,7},{5,6},
    {0,5},{1,4},
    {2,7},{3,6},
    {0,6},{2,4},
    {1,7},{3,5},
    /* 4 space diagonals (3-bit difference) */
    {0,7},{1,6},{2,5},{3,4}
};

/* 6 axis-aligned quad faces of the 2x2x2 cube */
static const int g_faces[6][4] = {
    {0,1,3,2},{4,5,7,6},   /* z=const */
    {0,2,6,4},{1,3,7,5},   /* x=const */
    {0,1,5,4},{2,3,7,6}    /* y=const */
};

/**
 * \brief Run one topology-correction scan pass over all 2x2x2 cubes.
 *
 * Binarises volume according to the pass polarity, scans every non-empty
 * 2x2x2 cube, and corrects configurations where chi_local equals the
 * bad value for the given connectivity. For each defective cube the
 * single foreground corner whose flip resolves the defect is chosen by
 * proximity of vol_prob[corner] to thresh — the most probability-ambiguous
 * voxel is preferred because it has the least anatomical commitment. Falls
 * back to flipping all foreground corners when no single flip suffices.
 * Changes are accumulated then written back to volume before returning.
 *
 * \param volume        (in/out) current working volume
 * \param vol_prob      (in)     original pre-correction probability snapshot
 * \param thresh        (in)     foreground threshold (>= is foreground)
 * \param mn            (in)     value written when a voxel is removed
 * \param mx            (in)     value written when a voxel is added
 * \param nx,ny,nz      (in)     volume dimensions
 * \param conn          (in)     connectivity: 18 or 26
 * \param vol_euler     (in/out) scratch buffer, nvol floats
 * \param vol_euler_orig(in/out) scratch buffer, nvol floats
 * \param vol_bin       (in/out) scratch buffer, nvol unsigned shorts
 * \return number of defective cubes corrected
 */
static int
run_topology_pass(float *volume,float *vol_changed,  const float *vol_prob, 
                  float thresh, float mn, float mx, int nx, int ny, int nz,
                  int conn, float *vol_euler, float *vol_euler_orig,
                  unsigned short *vol_bin)
{
    int i, j, x, y, z;
    int nvol    = nx * ny * nz;
    int bad_chi = (conn == 18) ? 2 : -6;
    int n_edges = (conn == 18) ? 24 : 28;
    const int (*edges)[2] = (conn == 18) ? g_edges_18conn : g_edges_26conn;

    /* Reset working copies from current volume state */
    for (i = 0; i < nvol; i++)
    {
        vol_euler[i]      = (volume[i] >= thresh) ? 1.0f : 0.0f;
        vol_euler_orig[i] = vol_euler[i];
    }

    /* Binarise according to this pass's polarity */
    for (i = 0; i < nvol; i++)
        vol_bin[i] = (unsigned short)((volume[i] >= thresh) ? 1 : 0);

    int n_errors = 0;

    for (z = 0; z < nz - 1; z++)
    {
        for (y = 0; y < ny - 1; y++)
        {
            for (x = 0; x < nx - 1; x++)
            {
                /* Linear indices for the 8 corners */
                int cidx[8] = {
                    IDX(x,   y,   z,   nx, ny),
                    IDX(x+1, y,   z,   nx, ny),
                    IDX(x,   y+1, z,   nx, ny),
                    IDX(x+1, y+1, z,   nx, ny),
                    IDX(x,   y,   z+1, nx, ny),
                    IDX(x+1, y,   z+1, nx, ny),
                    IDX(x,   y+1, z+1, nx, ny),
                    IDX(x+1, y+1, z+1, nx, ny)
                };

                int cube[8];
                for (i = 0; i < 8; i++)
                    cube[i] = vol_bin[cidx[i]];

                /* Compute chi_local = V - E + F */
                int V = 0, E = 0, F = 0;
                for (i = 0; i < 8; i++) V += cube[i];

                /* Skip completely empty cubes */
                if (V == 0)
                    continue;

                for (i = 0; i < n_edges; i++)
                    if (cube[edges[i][0]] && cube[edges[i][1]]) E++;
                for (i = 0; i < 6; i++)
                    if (cube[g_faces[i][0]] && cube[g_faces[i][1]] &&
                        cube[g_faces[i][2]] && cube[g_faces[i][3]]) F++;

                if (V - E + F != bad_chi)
                    continue;

                /* Find the single foreground corner whose flip resolves the
                 * defect. Prefer the corner with the smallest |prob - thresh|:
                 * that voxel has the lowest anatomical confidence and is the
                 * most "free" to change without distorting anatomy.
                 * vol_prob holds the original pre-correction probabilities so
                 * that this metric is not degraded by earlier iterations.     */
                int   best_corner = -1;
                float best_dist   = FLT_MAX;

                for (i = 0; i < 8; i++)
                {
                    if (!cube[i])
                        continue;

                    int test[8];
                    for (j = 0; j < 8; j++) test[j] = cube[j];
                    test[i] = 0;

                    int tV = 0, tE = 0, tF = 0;
                    for (j = 0; j < 8; j++) tV += test[j];
                    for (j = 0; j < n_edges; j++)
                        if (test[edges[j][0]] && test[edges[j][1]]) tE++;
                    for (j = 0; j < 6; j++)
                        if (test[g_faces[j][0]] && test[g_faces[j][1]] &&
                            test[g_faces[j][2]] && test[g_faces[j][3]]) tF++;

                    if (tV - tE + tF == bad_chi)
                        continue;  /* this flip does not resolve the defect */

                    float dist = (float)fabs((double)vol_prob[cidx[i]] - (double)thresh);
                    if (dist < best_dist)
                    {
                        best_dist   = dist;
                        best_corner = i;
                    }
                }

                if (best_corner >= 0)
                {
                    /* Minimal fix: flip just the most ambiguous corner */
                    vol_euler[cidx[best_corner]] = 0.0f;
                }
                else
                {
                    /* No single flip resolves the defect — flip all foreground
                     * corners of this cube as fallback                        */
                    for (i = 0; i < 8; i++)
                        if (cube[i])
                            vol_euler[cidx[i]] = 0.0f;
                }
                n_errors++;
            }
        }
    }

    /* Apply accumulated changes to volume */
    for (i = 0; i < nvol; i++)
    {
        if (vol_euler[i] < vol_euler_orig[i])
        {
            volume[i] = mn;
            vol_changed[i] = -1.0f;
        } else if (vol_euler[i] > vol_euler_orig[i])
        {
            volume[i] = mx;
            vol_changed[i] = 1.0f;
        }
    }

    return n_errors;
}

/**
 * \brief Correct topological defects in binary volume using Euler characteristic.
 *
 * Eliminates handles and cavities in a segmented volume using local Euler
 * characteristic analysis. Iteratively calls run_topology_pass() with an
 * adaptive connectivity-alternation strategy:
 *
 *   1. The first connectivity sub-pass always runs (alternates 18/26 each
 *      outer iteration to avoid systematic bias).
 *   2. The second connectivity sub-pass runs only if the first made changes.
 *      This avoids oscillation where 18-conn fixes create 26-conn defects and
 *      vice versa — the second pass is only needed when cross-type defects
 *      were actually introduced.
 *
 * Anatomical confidence weighting: a snapshot of the original probability
 * volume is taken once before any iterations. The per-corner selection of
 * which voxel to flip always uses this snapshot so that corrections from
 * early iterations do not degrade the weighting for later ones. Voxels
 * whose probability is closest to thresh (most ambiguous, e.g. in sulcal
 * fundi) are preferred, minimising anatomical distortion.
 *
 * \param volume        (in/out) floating-point volume (probability or binary)
 * \param thresh        (in)     voxels >= thresh are foreground
 * \param dims          (in)     [nx, ny, nz] volume dimensions
 * \param conn_arr      (in)     two connectivity values, e.g. {18, 26}
 * \return void
 */
void correct_topology(float *volume, float *vol_changed, float thresh, int dims[3], int conn_arr[2])
{
    int i, iter;
    int nx   = dims[0], ny = dims[1], nz = dims[2];
    int nvol = nx * ny * nz;

    float mn = (float)get_min(volume, nvol, 1, DT_FLOAT32);
    float mx = (float)get_max(volume, nvol, 0, DT_FLOAT32);

    /* Snapshot of the original probabilities — used throughout all iterations
     * for the per-corner anatomical-confidence weighting. Taking it once here
     * ensures that corrections applied in early iterations do not distort the
     * weighting metric for subsequent iterations.                             */
    float          *vol_prob       = (float *)malloc(sizeof(float) * nvol);
    float          *vol_euler      = (float *)malloc(sizeof(float) * nvol);
    float          *vol_euler_orig = (float *)malloc(sizeof(float) * nvol);
    unsigned short *vol_bin        = (unsigned short *)malloc(sizeof(unsigned short) * nvol);

    if (!vol_prob || !vol_euler || !vol_euler_orig || !vol_bin)
    {
        fprintf(stderr, "Memory allocation error\n");
        exit(EXIT_FAILURE);
    }

    for (i = 0; i < nvol; i++)
        vol_prob[i] = volume[i];

    for (iter = 0; iter < 50; iter++)
    {
        int n_total = 0;
        int n_err;

        /* Alternate which connectivity type leads each outer iteration.
         * This prevents the same connectivity from accumulating a systematic
         * directional correction bias across many iterations.                */
        int ci0 = iter % 2;       /* index of the leading connectivity     */
        int ci1 = 1 - ci0;        /* index of the following connectivity   */

        /* --- Foreground-removal pass, leading connectivity (always runs) --- */
        n_err = run_topology_pass(volume, vol_changed, vol_prob, thresh, mn, mx,
                                  nx, ny, nz, conn_arr[ci0],
                                  vol_euler, vol_euler_orig, vol_bin);
        n_total += n_err;

        /* --- Foreground-removal pass, following connectivity (adaptive) ---
         * Only run if the leading pass made changes. Cross-type defects can
         * only appear if the leading pass altered voxels; skipping this when
         * n_err==0 prevents unnecessary work and avoids the 18↔26 oscillation
         * that can pin the iteration count at the maximum.                   */
        if (n_err > 0)
        {
            n_err = run_topology_pass(volume, vol_changed, vol_prob, thresh, mn, mx,
                                      nx, ny, nz, conn_arr[ci1],
                                      vol_euler, vol_euler_orig, vol_bin);
            n_total += n_err;
        }

        if (n_total == 0)
            break;
    }

    free(vol_prob);
    free(vol_euler);
    free(vol_euler_orig);
    free(vol_bin);
}

private void
clear_slice(
    double **slice,
    int sizes[3])
{
    int x, y;

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
    int z,
    int sizes[3])
{
    int x, y;

    for (x = 0; x < sizes[0]; x++)
        for (y = 0; y < sizes[1]; y++)
        {
            slice[x][y] = (double)vol[sub2ind(x, y, z, sizes)];
        }
}

private double
get_slice_value(
    double ***slices,
    int x_size,
    int y_size,
    int z,
    int x,
    int y)
{
    if (x < 0 || x >= x_size || y < 0 || y >= y_size)
        return (0.0);
    else
        return (slices[z][x][y]);
}

private void
clear_points(
    int x_size,
    int y_size,
    int max_edges,
    int ***point_ids)
{
    int x, y, edge;

    for (x = 0; x < x_size + 2; x++)
        for (y = 0; y < y_size + 2; y++)
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
    Point *point)
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
    mat44 nii_mat,
    BOOLEAN binary_flag,
    double min_threshold,
    double max_threshold,
    int ***point_ids[],
    polygons_struct *polygons)
{
    int voxel[N_DIMENSIONS], edge, point_index;
    int edge_voxel[N_DIMENSIONS];
    double v[N_DIMENSIONS];
    Point world_point;
    Point_classes point_class;

    voxel[X] = x + point->coord[X];
    voxel[Y] = y + point->coord[Y];
    voxel[Z] = point->coord[Z];
    edge = point->edge_intersected;

    point_index = point_ids[voxel[Z]][voxel[X] + 1][voxel[Y] + 1][edge];
    if (point_index < 0)
    {
        edge_voxel[X] = point->coord[X];
        edge_voxel[Y] = point->coord[Y];
        edge_voxel[Z] = point->coord[Z];
        point_class = get_isosurface_point(corners, edge_voxel, edge,
                                           binary_flag,
                                           min_threshold, max_threshold, v);

        get_world_point(v[Z] + (float)slice_index,
                        v[X] + (float)x, v[Y] + (float)y,
                        nii_mat, &world_point);

        point_index = polygons->n_points;
        ADD_ELEMENT_TO_ARRAY(polygons->points, polygons->n_points,
                             world_point, CHUNK_SIZE);

        point_ids[voxel[Z]][voxel[X] + 1][voxel[Y] + 1][edge] = point_index;
    }

    return (point_index);
}

private void
extract_surface(
    Marching_cubes_methods method,
    BOOLEAN binary_flag,
    double min_threshold,
    double max_threshold,
    double valid_low,
    double valid_high,
    int x_size,
    int y_size,
    double ***slices,
    double min_label,
    double max_label,
    double ***label_slices,
    int slice_index,
    BOOLEAN right_handed,
    mat44 nii_mat,
    int ***point_ids[],
    polygons_struct *polygons)
{
    int x, y, *sizes, tx, ty, tz, n_polys, ind;
    int p, point_index, poly, size, start_points, dir;
    voxel_point_type *points;
    double corners[2][2][2], label;
    BOOLEAN valid;

    for (x = -1; x < x_size; x++)
        for (y = -1; y < y_size; y++)
        {
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
                start_points = sizes[0] - 1;
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
                                         polygons->end_indices[polygons->n_items - 1],
                                         point_index, CHUNK_SIZE);
                }

                if (right_handed)
                    start_points += size;
                else if (poly < n_polys - 1)
                    start_points += sizes[poly + 1];
            }
        }
}

void
/**
 * \brief Extract polygonal surface mesh from volumetric data using marching cubes.
 *
 * Generates an isosurface mesh from a 3D volume by sweeping slices through the data
 * and applying the marching cubes algorithm. Handles multi-label images via
 * min/max_label ranges and applies threshold bounding boxes. Stores mesh vertices
 * in world coordinates using the NIfTI affine transformation matrix. Computes
 * surface normals after extraction. Supports various marching cubes methods and
 * hand-edness correction for proper normal orientation.
 *
 * \param vol             (in)  input 3D volume as linearized float array
 * \param sizes           (in)  array [nx, ny, nz] dimensions of volume
 * \param min_label       (in)  minimum label value for extraction (-1 to ignore)
 * \param max_label       (in)  maximum label value for extraction (-1 to ignore)
 * \param nii_mat         (in)  NIfTI affine 4x4 matrix (voxel to world coordinates)
 * \param method          (in)  marching cubes algorithm variant
 * \param binary_flag     (in)  treat volume as binary (0/1) if TRUE
 * \param min_threshold   (in)  minimum intensity threshold for surface
 * \param max_threshold   (in)  maximum intensity threshold for surface
 * \param valid_low       (in)  lowest valid data value
 * \param valid_high      (in)  highest valid data value
 * \param polygons        (out) output mesh structure; allocated and populated by function
 * \param verbose         (in)  1 to print progress, 0 for silent operation
 * \return void
 */
extract_isosurface(
    float *vol,
    int sizes[3],
    double min_label,
    double max_label,
    mat44 nii_mat,
    Marching_cubes_methods method,
    BOOLEAN binary_flag,
    double min_threshold,
    double max_threshold,
    double valid_low,
    double valid_high,
    polygons_struct *polygons,
    int verbose)
{
    int n_slices, x_size, y_size, slice;
    int ***point_ids[2], ***tmp_point_ids;
    int max_edges;
    double **slices[2], **tmp_slices;
    double **label_slices[2];
    progress_struct progress;
    Surfprop spr;
    Point point000, point100, point010, point001;
    Vector v100, v010, v001, perp;
    BOOLEAN right_handed;

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

    ALLOC3D(point_ids[0], x_size + 2, y_size + 2, max_edges);
    ALLOC3D(point_ids[1], x_size + 2, y_size + 2, max_edges);

    clear_slice(slices[1], sizes);

    clear_points(x_size, y_size, max_edges, point_ids[0]);
    clear_points(x_size, y_size, max_edges, point_ids[1]);

    Surfprop_a(spr) = 0.3f;
    Surfprop_d(spr) = 0.6f;
    Surfprop_s(spr) = 0.6f;
    Surfprop_se(spr) = 30.0f;
    Surfprop_t(spr) = 1.0f;
    initialize_polygons(polygons, WHITE, &spr);

    if (verbose)
        initialize_progress_report(&progress, FALSE, n_slices + 1, "Extracting Surface");

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

        if (verbose)
            update_progress_report(&progress, slice + 2);
    }

    if (verbose)
        terminate_progress_report(&progress);

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
/**
 * \brief Extract brain surface mesh with advanced preprocessing and topology correction.
 *
 * Comprehensive surface extraction pipeline including smoothing, median filtering,
 * gyral masking, and topology correction. Applies optional edge-preserving smoothing,
 * median filter to strengthen structures, largest component selection, hole filling,
 * and topology correction using Euler characteristic. Iteratively applies Laplacian
 * smoothing and morphological post-processing. Creates separate inner and outer
 * surfaces with appropriate labeling for CAT12 cortical mesh processing.
 *
 * \param input_float       (in)  input 3D probabilistic tissue segmentation
 * \param nii_ptr           (in)  NIfTI image header with voxel dimensions and affine
 * \param label             (in)  optional tissue label mask (NULL to skip)
 * \param min_threshold     (in)  isosurface threshold value (typically 0.5 for probabilities)
 * \param pre_fwhm          (in)  Gaussian smoothing FWHM in mm (0 to skip)
 * \param iter_laplacian    (in)  number of Laplacian smoothing iterations
 * \param dist_morph        (in)  distance offset for morphological expansion (mm)
 * \param n_median_filter   (in)  iterations of median filtering to apply
 * \param n_iter            (in)  total outer loop iterations
 * \param strength_gyri_mask (in) weighting factor for gyral preservation masking (0-1)
 * \param verbose           (in)  1 to print progress, 0 for silent
 * \return Allocated object_struct containing pial surface polygons; caller must free
 */
object_struct *apply_marching_cubes(float *input_float, nifti_image *nii_ptr,
                                    float *label, double min_threshold, double pre_fwhm,
                                    int iter_laplacian, double dist_morph, int n_median_filter,
                                    int n_iter, double strength_gyri_mask, int verbose)
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

    if (!vol_float || !vol_uint16 || !vol_uint8)
    {
        fprintf(stderr, "Memory allocation error\n");
        exit(EXIT_FAILURE);
    }

    /* Normalize input if necessary */
    if (min_threshold > 1.0)
    {
        double max_value = get_max(input_float, nvol, 0, DT_FLOAT32);
        for (i = 0; i < nvol; i++)
            input_float[i] /= (float)max_value;
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
    if (pre_fwhm)
    {
        float *grad = (float *)malloc(nvol * sizeof(float));
        if (!grad)
        {
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

        for (i = 0; i < nvol; i++)
        {
            float weight = (grad[i] - mean_grad) / (max_grad - mean_grad);
            weight = (weight < 0.0) ? 0.0f : weight;
            weight = (weight > 1.0) ? 1.0f : weight;
            weight *= weight;

            input_float[i] = (1.0 - weight) * input_float[i] + weight * vol_float[i];
        }
        free(grad);
    }

    if ((label != NULL) && (strength_gyri_mask != 0.0))
    {
        /* Estimate smooth gyrus mask */
        smooth_gyri_mask(label, vol_float, dims, voxelsize, 1.5, 8.0);

        /* Scale smooth mask to a range of -strength_gyri_mask..strength_gyri_mask */
        for (i = 0; i < nvol; i++)
            vol_float[i] = vol_float[i] * strength_gyri_mask * 2 - strength_gyri_mask;

        /* And subtract scaled mask to input so that gyri have lower isovalues by 0.1
           (to preserve gyral crowns) and sulci higher isovalues by -0.1 to preserve
           sulcal spaces */
        for (i = 0; i < nvol; i++)
            input_float[i] -= vol_float[i];
    }

    /* Apply iterative median filter to strengthen structures */
    unsigned char *mask = (unsigned char *)malloc(nvol * sizeof(unsigned char));
    for (i = 0; i < nvol; i++)
        mask[i] = input_float[i] != 0;

    median3(input_float, mask, dims, n_median_filter, DT_FLOAT32);
    free(mask);

    /* Keep largest cluster and fill holes */
    keep_largest_cluster(input_float, min_threshold, dims, DT_FLOAT32, 0, 1, 18);
    fill_holes(input_float, dims, min_threshold, -1.0, DT_FLOAT32);

    /* Correct topology */
    int conn_arr[2] = {18, 26};
    correct_topology(input_float, vol_changed, min_threshold, dims, conn_arr);

    /* Apply genus-0 correction */
    genus0parameters g0[1];
    genus0init(g0);
    for (i = 0; i < N_DIMENSIONS; i++)
        g0->dims[i] = dims[i];

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
    if (dist_morph == FLT_MAX)
    {
        double dist_values[] = {-1.0, -0.5, 0.0, 0.5, 1.0, 1.5};
        int change_values[] = {0, 0, 0, 0, 0, 0};
        int n_values = sizeof(change_values) / sizeof(change_values[0]);
        for (k = 0; k < n_values; k++)
        {

            for (i = 0; i < nvol; i++)
                vol_uint16[i] = (input_float[i] >= min_threshold) ? 1 : 0;

            if (dist_values[k] > 0.0)
                dist_close(vol_uint16, dims, voxelsize, dist_values[k], 0.0, DT_UINT16);
            else if (dist_values[k] < 0.0)
                dist_open(vol_uint16, dims, voxelsize, -dist_values[k], 0.0, DT_UINT16);

            /* call genus0 for the 1st time */
            g0->input = vol_uint16;
            g0->cut_loops = 0;
            g0->connectivity = 6;
            g0->alt_value = 1;
            g0->alt_contour_value = 1;
            if (genus0(g0))
                return (NULL); /* check for error */

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
            if (genus0(g0))
                return (NULL);

            /* save changes */
            for (i = 0; i < nvol; i++)
                change_values[k] += (int)fabs((float)g0->output[i] - (float)g0->input[i]);
            if (change_values[k] == 0)
            {
                n_values = k + 1;
                break;
            }
        }
        int ind_min_value;
        int min_value = 1E9;
        for (k = 0; k < n_values; k++)
        {
            if (change_values[k] < min_value)
            {
                min_value = change_values[k];
                ind_min_value = k;
            }
        }
        best_dist = dist_values[ind_min_value];
        best_change_values = change_values[ind_min_value];
    }
    else
    {
        best_dist = dist_morph;
        best_change_values = 1E3;
    }
    if (verbose)
        printf("Optimal dist-parameter for morphological operations: %f\n", best_dist);

    /* Convert to binary image */
    for (i = 0; i < nvol; i++)
        vol_uint16[i] = (input_float[i] >= min_threshold) ? 1 : 0;

    while ((EC != 2) && (count < n_iter))
    {

        /* Only move on if topology correction is still necessary */
        if (best_change_values > 0)
        {
            if (best_dist > 0.0)
                dist_close(vol_uint16, dims, voxelsize, best_dist, 0.0, DT_UINT16);
            else if (best_dist < 0.0)
                dist_open(vol_uint16, dims, voxelsize, -best_dist, 0.0, DT_UINT16);

            /* call genus0 for the 1st time */
            g0->input = vol_uint16;
            g0->cut_loops = 0;
            g0->connectivity = 6;
            g0->alt_value = 1;
            g0->alt_contour_value = 1;
            if (genus0(g0))
                return (NULL); /* check for error */

            /* save results as next input */
            for (i = 0; i < nvol; i++)
                vol_uint16[i] = g0->output[i];

            /* save changes */
            for (i = 0; i < nvol; i++)
                vol_changed[i] += (float)(count + 2) * ((float)g0->output[i] - (float)g0->input[i]);

            /* call genus0 a 2nd time with other parameters */
            g0->input = vol_uint16;
            g0->cut_loops = 1;
            g0->connectivity = 18;
            g0->alt_value = 0;
            g0->alt_contour_value = 0;
            if (genus0(g0))
                return (NULL);

            /* save changes */
            for (i = 0; i < nvol; i++)
            {
                vol_changed[i] += (float)(count + 2) * ((float)g0->output[i] - (float)g0->input[i]);
                if (vol_changed[i] != 0.0)
                    count_change++;
            }

            if (n_median_filter)
            {
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
        if (verbose)
            printf("Euler characteristics after %d iterations: %d (%d voxel changed).\n", 
                count, EC, count_change);

        /* save results as next input */
        for (i = 0; i < nvol; i++)
            vol_uint16[i] = g0->output[i];
    }

    if (!verbose) 
        printf("Euler characteristics after %d iterations: %d (%d voxel changed).\n", 
            count, EC, count_change);

    /* Laplacian smoothing to reduce noise in the mesh */
    if (iter_laplacian > 0)
        smooth_laplacian(polygons, iter_laplacian, 0.0, 0.5);

    /* Replace input values by change map */
    for (i = 0; i < nvol; i++)
        input_float[i] = vol_changed[i];

    free(vol_changed);
    free(vol_uint8);
    free(vol_uint16);
    free(vol_float);

    return object;
}

/* Function to apply marching cubes and extract polygons without any corrections */
/**
 * \brief Fast surface mesh extraction with minimal preprocessing.
 *
 * Streamlined isosurface extraction optimized for speed, skipping topology correction
 * and advanced filtering. Applies only essential preprocessing: largest component
 * selection and genus-0 correction (optional). Suitable for rapid prototyping or
 * when input data is already well-segmented with correct topology. Often used as a
 * fallback when full preprocessing would introduce too many smoothing artifacts.
 *
 * \param input_float     (in)  input 3D volume (typically already binary/thresholded)
 * \param nii_ptr         (in)  NIfTI image header with dimensions and affine transform
 * \param min_threshold   (in)  isosurface threshold value
 * \param iter_laplacian  (in)  number of Laplacian smoothing iterations post-extraction
 * \param verbose         (in)  1 to print progress, 0 for silent operation
 * \return Allocated object_struct containing surface polygons; caller must free
 */
object_struct *apply_marching_cubes_fast(float *input_float, nifti_image *nii_ptr,
                                         double min_threshold, int iter_laplacian, int verbose)
{
    double voxelsize[N_DIMENSIONS];
    double best_dist;
    int dims[MAX_DIMENSIONS], i, k, nvol;
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
    float *vol_float = (float *)malloc(nvol * sizeof(float));
    unsigned short *vol_uint16 = (unsigned short *)malloc(nvol * sizeof(unsigned short));
    unsigned char *vol_uint8 = (unsigned char *)malloc(nvol * sizeof(unsigned char));

    if (!vol_float || !vol_uint16)
    {
        fprintf(stderr, "Memory allocation error\n");
        exit(EXIT_FAILURE);
    }

    /* Normalize input if necessary */
    if (min_threshold > 1.0)
    {
        double max_value = get_max(input_float, nvol, 0, DT_FLOAT32);
        for (i = 0; i < nvol; i++)
            input_float[i] /= (float)max_value;
        min_threshold /= max_value;
    }

    /* Keep largest cluster */
    keep_largest_cluster(input_float, min_threshold, dims, DT_FLOAT32, 0, 1, 18);

    /* Convert to binary image */
    for (i = 0; i < nvol; i++)
        vol_uint16[i] = (input_float[i] >= min_threshold) ? 1 : 0;

    /* Do not apply genus-0 correction */
    genus0parameters g0[1];
    genus0init(g0);
    for (i = 0; i < N_DIMENSIONS; i++)
        g0->dims[i] = dims[i];

    g0->input = vol_uint16;
    g0->connected_component = 1;
    g0->value = 1;
    g0->contour_value = 1;
    g0->any_genus = 1;
    g0->biggest_component = 1;
    g0->pad[0] = g0->pad[1] = g0->pad[2] = 2;
    g0->ijk2ras = NULL;
    g0->verbose = 0;
    g0->return_surface = 0;
    g0->extraijkscale[2] = 1;
    g0->cut_loops = 0;
    g0->connectivity = 6;
    g0->alt_value = 1;
    g0->alt_contour_value = 1;

    if (genus0(g0))
        return (NULL); /* check for error */

    keep_largest_cluster(g0->output, min_threshold, dims, DT_UINT16, 0, 1, 18);

    for (i = 0; i < nvol; i++)
        vol_float[i] = (float)g0->output[i];

    object_struct *object = create_object(POLYGONS);
    polygons_struct *polygons = get_polygons_ptr(object);

    extract_isosurface(vol_float, dims, 0.0, -1.0, nii_mat, method, FALSE,
                       min_threshold, min_threshold, 0.0, -1.0, polygons, verbose);

    compute_polygon_normals(polygons);
    check_polygons_neighbours_computed(polygons);
    int n_out = separate_polygons(polygons, -1, &object2);
    triangulate_polygons(get_polygons_ptr(object2[0]), polygons);

    int EC = euler_characteristic(polygons, verbose);
    if (verbose)
        printf("Euler characteristics: %d\n", EC);

    /* Laplacian smoothing to reduce noise in the mesh */
    if (iter_laplacian > 0)
        smooth_laplacian(polygons, iter_laplacian, 0.0, 0.5);

    free(vol_uint16);
    free(vol_float);

    return object;
}
