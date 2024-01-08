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

#define   CHUNK_SIZE    1000000

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
    polygons_struct *polygons)
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

    initialize_progress_report(&progress, FALSE, n_slices+1, "Extracting Surface");

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

        update_progress_report(&progress, slice+2);
    }

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

