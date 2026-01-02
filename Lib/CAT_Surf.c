/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 *
 */

/*
 * \file CAT_Surf.c
 * \brief Surface utilities for CAT: geometry helpers, smoothing, inflation, metric distances,
 *        central↔pial displacement, and QEM-based mesh reduction.
 *
 * The functions in this file operate on BICPL \c polygons_struct meshes (points, indices,
 * end-indices). Several routines create or rely on polygon neighbourhoods and/or bintrees;
 * see individual function docs for requirements and side-effects.
 *
 * Some code is adapted from caret 5.3 (BrainModelSurface.cxx) and BICPL.
 *
 */

#include <quadric.h>
#include "CAT_Surf.h"
#include "CAT_SafeAlloc.h"

/**
 * \brief Apply mixed boundary conditions to grid indices.
 *
 * Implements circulant (wrap-around) boundary condition in x and a mirrored
 * Neumann condition in y for a 2-D lattice of size \c dm[0]×\c dm[1]. Useful
 * for projecting 3-D vertex coordinates into a fixed voxel lattice.
 *
 * \param i  input x index (can be negative/out of bounds)
 * \param j  input y index (can be negative/out of bounds)
 * \param dm two-element array with lattice extents \c {nx, ny}
 * \return flattened index in \c [0, nx*ny)
 */
int
bound(int i, int j, int dm[])
{
    int i1, j1;
    int m2;

    /* circulant boundary condition for x-coordinates */
    i1 = (i >= 0 ? (i % dm[0]) : (dm[0] + (i % dm[0])) % dm[0]);

    /* Neumann boundary condition for y-coordinates */
    if (dm[1] == 1) {
        j1 = 0;
    } else {
        m2 = dm[1] * 2;
        j = (j < 0) ? (-j - m2*(-j / m2) - 1) : (j - m2*(j / m2));
        if (dm[1] <= j)
            j1 = m2 - j - 1;
        else
            j1 = j;
    }

    return(i1 + dm[0]*j1);
}


/* copy Point values to/from double array of length 3 */
/**
 * \brief Mesh geometry utility.
 *
 * Function: to_array
 *
 * \param p (Point *)
 * \param xyz (double *)
 */
void
to_array(Point *p, double *xyz) {
    int i;

    for (i = 0; i < 3; i++) {
        xyz[i] = (double) Point_coord(*p,i);
    }
}

/**
 * \brief Mesh geometry utility.
 *
 * Function: from_array
 *
 * \param xyz (double *)
 * \param p (Point *)
 */
void
from_array(double *xyz, Point *p) {
    int i;

    for (i = 0; i < 3; i++) {
        Point_coord(*p,i) = xyz[i];
    }
}

double *
get_surface_ratio(double radius, polygons_struct *polygons, int normalize)
{
    int i, j, x, y, z, a, b, c, nan = 0;
    int size, poly_size;
    double *lf, *avol;
    double area, asum, total_area, mean_surface_ratio;
    char str[512];
    Point points[MAX_POINTS_PER_POLYGON];
    
    avol = SAFE_CALLOC(double, 256*256*256);

    lf = SAFE_CALLOC(double, polygons->n_points);
        
    for (i = 0; i < 256*256*256; i++)
        avol[i] = 0.0;

    total_area = 0.0;
    for (i = 0; i < polygons->n_items; i++) {
        size = get_polygon_points(polygons, i, points);
        area = get_polygon_surface_area(size, points);
        if (isnan(area)) area = 0.0;
        total_area += area;

        poly_size = GET_OBJECT_SIZE(*polygons, i);

        for (j = 0; j < poly_size; j++) {
            *points = polygons->points[polygons->indices[
                   POINT_INDEX(polygons->end_indices, i, j)]];
            x = 128 + (int) Point_x(*points);
            y = 128 + (int) Point_y(*points);
            z = 128 + (int) Point_z(*points);

            if (x >= 0 && x < 256 &&
              y >= 0 && y < 256 &&
              z >= 0 && z < 256)
                avol[z*256*256 + y*256 + x] += area;
        }
    }
    
    if (radius < 0) {
        /* this was empirically estimated using global gyrification index and surface ratios using
           different ranges of radii */
        radius = 0.0000471*total_area + 7.456;
    }
    
    if (normalize) {
        /* the individual surface area is related to the total surface area of the template meshes
           ?h.central.Template_T1_IXI555_MNI152_GS.gii from CAT12 */
        radius = radius*sqrt(total_area/90000.0);
        printf("Estimated radius: %3.1f\n", radius);
    }

    mean_surface_ratio = 0.0;
    for (i = 0; i < polygons->n_points; i++) {

        asum = 0;
        for (x = -radius; x <= radius; x++) {
            for (y = -radius; y <= radius; y++) {
                for (z = -radius; z <= radius; z++) {
                    if (x*x + y*y + z*z < radius*radius) {
                        *points = polygons->points[i];
                        a = x + 128 + (int) Point_x(*points);
                        b = y + 128 + (int) Point_y(*points);
                        c = z + 128 + (int) Point_z(*points);
                        if (a >= 0 && a < 256 &&
                          b >= 0 && b < 256 &&
                          c >= 0 && c < 256)
                            asum += avol[c*65536 + b*256 + a];
                    }
                }
            }
        }
        /*
         * The area of a triangle completely inside the sphere will
         * will be added 3 times (one per vertex).   The area of a
         * the area of a triangle partially inside the sphere will
         * will be accounted proportionally to the number of vertices
         * of vertices inside the sphere.
         */
        lf[i] = asum / (_PI*radius*radius) / 3.0;
        
        if (lf[i] != lf[i])
            nan++;

        mean_surface_ratio += lf[i];
    }

    free(avol);
    
    printf("Mean surface ratio: %3.3f\n",mean_surface_ratio/polygons->n_points);
    
    if (nan)
        printf("ERROR: there are %i NaN\n", nan);
    
    return lf;
}

/* assign each point the average of area values for neighboring polygons */
/**
 * \brief Resample per-vertex area values to the sphere and normalize to equal area.
 *
 * Resamples \c area_values from \c polygons to a spherical mesh of identical topology,
 * normalizes locally by the sphere’s area at each vertex, and globally so that the
 * sum of values equals the mesh surface area.
 *
 * \param polygons        source mesh.
 * \param sphere          target spherical mesh (same topology).
 * \param area_values     output array (length \c n_points).
 * \return total surface area of the resampled spherical mesh.
 */
double
get_area_of_points_normalized_to_sphere(polygons_struct *polygons, polygons_struct *sphere, double *area_values)
{
    double sphere_area, native_area, *area_values_resampled, *areas_sphere, area_sum = 0.0;
    int i, n_points;
    object_struct **object;
    object_struct *resampled_sphere_object;
    polygons_struct *resampled, *resampled_sphere;
    Point center;
    
    n_points = 81920;
    
    object = resample_surface(polygons, sphere, n_points, NULL, NULL);
    resampled = get_polygons_ptr(*object);
    
    /* create tetrahedral sphere according to resampled surface */
    resampled_sphere_object = create_object(POLYGONS);
    resampled_sphere = get_polygons_ptr(resampled_sphere_object);
    fill_Point(center, 0.0, 0.0, 0.0);
    create_tetrahedral_sphere(&center, 100.0, 100.0, 100.0,
                  resampled->n_items, resampled_sphere);

    area_values_resampled = SAFE_MALLOC(double, resampled->n_points);
    areas_sphere      = SAFE_MALLOC(double, resampled->n_points);

    /*
     * Use a sum-preserving per-vertex area definition (triangle contributions).
     * This matches the global constraint that vertex areas sum to surface area
     * and avoids bias from the legacy "average-of-neighbour-polygons" metric.
     */
    sphere_area = get_vertex_areas(resampled_sphere, areas_sphere);
    native_area = get_vertex_areas(resampled, area_values_resampled);
    
    /* normalize values to local sphere areas */
    for (i = 0; i < resampled->n_points; i++) 
        area_values_resampled[i] /= areas_sphere[i];

    resample_values_sphere(resampled_sphere, sphere, area_values_resampled, area_values, 1, 0);
    
    /* normalize values to overall area */
    for (i = 0; i < polygons->n_points; i++) 
        area_sum += area_values[i];
    for (i = 0; i < polygons->n_points; i++) 
        area_values[i] *= native_area/area_sum;

    free(area_values_resampled);
    free(areas_sphere);

    delete_object_list(1, object);
    delete_object_list(1, &resampled_sphere_object);
    
    return(native_area);
}

/*
 * Sum-preserving per-vertex surface area: distribute each polygon's area
 * equally across its vertices. The sum of vertex areas equals total area.
 */
double
get_vertex_areas(polygons_struct *polygons, double *vertex_areas)
{
    int ptidx, poly, vertidx, size;
    double area;
    double surface_area = 0.0;
    Point points[MAX_POINTS_PER_POLYGON];

    memset(vertex_areas, 0.0, sizeof(double) * polygons->n_points);

    for (poly = 0; poly < polygons->n_items; poly++) {
        size = get_polygon_points(polygons, poly, points);
        area = get_polygon_surface_area(size, points);
        if (isnan(area))
            area = 0.0;
        surface_area += area;

        for (vertidx = 0; vertidx < size; vertidx++) {
            ptidx = polygons->indices[
                        POINT_INDEX(polygons->end_indices,
                        poly, vertidx)];
            vertex_areas[ptidx] += area / (double) size;
        }
    }

    return surface_area;
}

/* assign each point the average of area values for neighboring polygons */
/**
 * \brief Compute or return a derived quantity from the mesh.
 *
 * Function: get_area_of_points
 *
 * \param polygons (polygons_struct *)
 * \param area_values (double *)
 * \return See function description for return value semantics.
 */
double
get_area_of_points(polygons_struct *polygons, double *area_values)
{
    int *pcount;
    int ptidx, poly, vertidx, size;
    double poly_size, area;
    double surface_area = 0.0;
    Point points[MAX_POINTS_PER_POLYGON];

    pcount = SAFE_MALLOC(int, polygons->n_points);
    memset(pcount, 0, sizeof(int) * polygons->n_points);
    memset(area_values, 0.0, sizeof(double) * polygons->n_points);

    for (poly = 0; poly < polygons->n_items; poly++) {    
        size = get_polygon_points(polygons, poly, points);
        area = get_polygon_surface_area(size, points);
        if (isnan(area)) area = 0.0;
        surface_area += area;

        poly_size = GET_OBJECT_SIZE(*polygons, poly);

        for (vertidx = 0; vertidx < poly_size; vertidx++) {
            ptidx = polygons->indices[
                        POINT_INDEX(polygons->end_indices,
                        poly, vertidx)];
            pcount[ptidx]++;
            area_values[ptidx] += area;
        }
    }

    for (ptidx = 0; ptidx < polygons->n_points; ptidx++) 
        if (pcount[ptidx] > 0)
            area_values[ptidx] /= pcount[ptidx];
  

    free(pcount);
    return(surface_area);
}

/**
 * \brief Compute or return a derived quantity from the mesh.
 *
 * Function: get_area_of_polygons
 *
 * \param polygons (polygons_struct *)
 * \param area_values (double *)
 * \return See function description for return value semantics.
 */
double
get_area_of_polygons(polygons_struct *polygons, double *area_values)
{
    int poly, size;
    double surface_area = 0.0;
    Point points[MAX_POINTS_PER_POLYGON];

    for (poly = 0; poly < polygons->n_items; poly++) {    
        size = get_polygon_points(polygons, poly, points);
        area_values[poly] = get_polygon_surface_area(size, points);
        if (isnan(area_values[poly])) area_values[poly] = 0.0;
        surface_area += area_values[poly];
    }
    return(surface_area);
}

/**
 * localstat_surface_double — Per-vertex local statistics on a surface (fixed 1-ring).
 *
 * Computes mean/median/std/min/max over each vertex’s immediate neighbours
 * (plus the center vertex) as defined by get_all_polygon_point_neighbours().
 * Optional mask (0 = skip) and multi-iteration behaviour matches localstat_double.
 *
 * Parameters
 * ----------
 *  polygons : mesh (BICPL polygons_struct)
 *  input    : in/out, length = polygons->n_points (per-vertex doubles)
 *  mask     : optional uchar mask, length = n_points (can be NULL)
 *  stat_func: F_MEAN, F_MEDIAN, F_STD, F_MIN, F_MAX
 *  iters    : number of iterations (>=1)
 */
void localstat_surface_double(polygons_struct *polygons,
                              double *input,
                              unsigned char *mask,
                              int stat_func,
                              int iters)
{
    if (!polygons || polygons->n_points <= 0) {
        fprintf(stderr, "localstat_surface_double: invalid polygons.\n");
        return;
    }
    if (iters < 1) iters = 1;

    /* Use neighbours to calculate local statistics */
    int *n_neighbours = NULL;
    int **neighbours  = NULL;
    int i, it, v;
    get_all_polygon_point_neighbours(polygons, &n_neighbours, &neighbours);

    double *buffer = SAFE_MALLOC(double, polygons->n_points);
    double *arr    = SAFE_MALLOC(double, polygons->n_points);
    for (i = 0; i < polygons->n_points; ++i) buffer[i] = input[i];

    for (it = 0; it < iters; ++it) {
        for (v = 0; v < polygons->n_points; ++v) {
            /* Skip masked/invalid centers exactly like the volume code. */
            if ((mask && mask[v] == 0) || isnan(input[v]) || !isfinite(input[v])) {
                buffer[v] = input[v];
                continue;
            }

            int n = 0;

            /* Include center value if valid. */
            if (!mask || mask[v] != 0) {
                if (isfinite(input[v]) && !isnan(input[v])) arr[n++] = input[v];
            }

            /* Add neighbours from the precomputed lists. */
            const int nn = n_neighbours[v];
            int *nb = neighbours[v];
            for (i = 0; i < nn; ++i) {
                const int u = nb[i];
                if (mask && mask[u] == 0) continue;
                if (!isfinite(input[u]) || isnan(input[u])) continue;
                arr[n++] = input[u];
            }

            if (n <= 0) { buffer[v] = input[v]; continue; }

            switch (stat_func) {
                case F_MEAN:   buffer[v] = get_mean_double(arr,   n, 0); break;
                case F_MEDIAN: buffer[v] = get_median_double(arr, n, 0); break;
                case F_STD:    buffer[v] = get_std_double(arr,    n, 0); break;
                case F_MIN:    buffer[v] = get_min_double(arr,    n, 0); break;
                case F_MAX:    buffer[v] = get_max_double(arr,    n, 0); break;
                default:
                    fprintf(stderr, "Data Function %d not handled\n", stat_func);
                    buffer[v] = input[v];
                    break;
            }
        }
        /* write back this iteration */
        for (i = 0; i < polygons->n_points; ++i) input[i] = buffer[i];
    }

    /* Cleanup (get_all_polygon_point_neighbours allocates a single flat block for neighbours[0]). */
    free(buffer);
    free(arr);
    if (n_neighbours) free(n_neighbours);
    if (neighbours)   { if (neighbours[0]) free(neighbours[0]); free(neighbours); }
}

/**
 * \brief Mixed boundary condition index mapping for a 2‑D lattice.
 *
 * Function: correct_bounds_to_target
 *
 * Subtracts the offset between mesh and target bounds from all vertices.
 * \param polygons (polygons_struct *)
 * \param target (polygons_struct *)
 */
void
correct_bounds_to_target(polygons_struct *polygons, polygons_struct *target)
{
    int i, j;
    double bounds_dest[6], bounds_src[6];
    Point center;

    /* Calc. sphere center based on bounds of input (correct for shifts) */    
    get_bounds(polygons, bounds_src);
    get_bounds(target,   bounds_dest);
    
    fill_Point(center, bounds_src[0]-bounds_dest[0],
               bounds_src[2]-bounds_dest[2], 
               bounds_src[4]-bounds_dest[4]);

    for (i = 0; i < polygons->n_points; i++) 
        for (j = 0; j < 3; j++)
            Point_coord(polygons->points[i], j) -= Point_coord(center, j);
}

/**
 * \brief Mixed boundary condition index mapping for a 2‑D lattice.
 *
 * Function: correct_bounds_to_target_with_scaling
 *
 * Scales each axis to match target range, recomputes bounds, then translates to align.
 * \param polygons (polygons_struct *)
 * \param target (polygons_struct *)
 */
void
correct_bounds_to_target_with_scaling(polygons_struct *polygons, polygons_struct *target)
{
    int i, j;
    double bounds_dest[6], bounds_src[6];
    Point center;

    /* Calc. sphere center based on bounds of input (correct for shifts) */    
    get_bounds(polygons, bounds_src);
    get_bounds(target,   bounds_dest);
    
    for (i = 0; i < polygons->n_points; i++) 
        for (j = 0; j < 3; j++)
            Point_coord(polygons->points[i], j) /= (bounds_src[2*j+1]-bounds_src[2*j])/(bounds_dest[2*j+1]-bounds_dest[2*j]);

    get_bounds(polygons, bounds_src);
    fill_Point(center, bounds_src[0]-bounds_dest[0],
               bounds_src[2]-bounds_dest[2], 
               bounds_src[4]-bounds_dest[4]);


    for (i = 0; i < polygons->n_points; i++) 
        for (j = 0; j < 3; j++)
            Point_coord(polygons->points[i], j) -= Point_coord(center, j);

}

/**
 * \brief Translate or align a mesh.
 *
 * Function: translate_to_center_of_mass
 *
 * \param polygons (polygons_struct *)
 */
void
translate_to_center_of_mass(polygons_struct *polygons)
{
    int i, j;
    double c[3] = {0.0, 0.0, 0.0};

    for (i = 0; i < polygons->n_points; i++) {
        for (j = 0; j < 3; j++)
            c[j] += Point_coord(polygons->points[i], j);
    }

    for (j = 0; j < 3; j++)
        c[j] /= polygons->n_points;

    for (i = 0; i < polygons->n_points; i++) {
        for (j = 0; j < 3; j++)
            Point_coord(polygons->points[i], j) -= c[j];
    }
}

/**
 * \brief Compute or return a derived quantity from the mesh.
 *
 * Function: get_sphere_radius
 *
 * \param polygons (polygons_struct *)
 * \return See function description for return value semantics.
 */
double
get_sphere_radius(polygons_struct *polygons)
{
    int i, j;
    double radius = 0.0, dist = 0.0;

    for (i = 0; i < polygons->n_points; i++) {
        dist = 0.0;
        for (j = 0; j < 3; j++)
            dist += Point_coord(polygons->points[i], j) *
                Point_coord(polygons->points[i], j);
        radius = MAX(radius, dist);      
    }
    return(sqrt(radius));
}


/**
 * \brief Set or scale a geometric quantity.
 *
 * Function: set_vector_length
 *
 * \param p (Point *)
 * \param newLength (double)
 */
void
set_vector_length(Point *p, double newLength)
{
    int j;
    const double len = sqrt(Point_x(*p)*Point_x(*p) + 
                Point_y(*p)*Point_y(*p) +
                Point_z(*p)*Point_z(*p));
    if (len > 0) {
        const double scale = newLength / len;
        for (j = 0; j < 3; j++)
            Point_coord(*p, j) *= scale;
    }
}

/**
 * \brief Compute or return a derived quantity from the mesh.
 *
 * Function: get_radius_of_points
 *
 * \param polygons (polygons_struct *)
 * \param radius (double *)
 */
void
get_radius_of_points(polygons_struct *polygons, double *radius)
{
    int i;
    Vector xyz;
  
    for (i = 0; i < polygons->n_points; i++) {
        fill_Vector(xyz, Point_x(polygons->points[i]),
                 Point_y(polygons->points[i]),
                 Point_z(polygons->points[i]));
        radius[i] = MAGNITUDE(xyz);
    }

}

/**
 * \brief Compute or return a derived quantity from the mesh.
 *
 * Function: get_bounds
 *
 * \param polygons (polygons_struct *)
 * \param param (double bounds[6)
 */
void
get_bounds(polygons_struct *polygons, double bounds[6])
{
    int i;

    bounds[0] = FLT_MAX;
    bounds[1] = -FLT_MAX;
    bounds[2] = FLT_MAX;
    bounds[3] = -FLT_MAX;
    bounds[4] = FLT_MAX;
    bounds[5] = -FLT_MAX;

    for (i = 0; i < polygons->n_points; i++) {
        bounds[0] = MIN(bounds[0], Point_x(polygons->points[i]));
        bounds[1] = MAX(bounds[1], Point_x(polygons->points[i]));
        bounds[2] = MIN(bounds[2], Point_y(polygons->points[i]));
        bounds[3] = MAX(bounds[3], Point_y(polygons->points[i]));
        bounds[4] = MIN(bounds[4], Point_z(polygons->points[i]));
        bounds[5] = MAX(bounds[5], Point_z(polygons->points[i]));
    }
}

/**
 * \brief Mesh geometry utility.
 *
 * Function: count_edges
 *
 * Duplicated edges (due to non-manifold connectivity) are detected and reported.
 * \param polygons (polygons_struct *)
 * \param n_neighbours (int)
 * \param neighbours (int *)
 * \return See function description for return value semantics.
 */
int
count_edges(polygons_struct *polygons, int n_neighbours[], int *neighbours[], int verbose)
{
    int p, n, nn, n_edges, n_duplicate_edges;

    n_edges = 0;
    n_duplicate_edges = 0;

    for (p = 0; p < polygons->n_points; p++) {
        for (n = 0; n < n_neighbours[p]; n++) {
            for (nn = n+1; nn < n_neighbours[p]; nn++) {
                if (neighbours[p][nn] == neighbours[p][n])
                    break;
            }
            if (nn < n_neighbours[p])
                n_duplicate_edges++;
            else if (p < neighbours[p][n])
                n_edges++;
        }
    }

    if ((n_duplicate_edges > 0) && verbose)
        printf("N duplicate edges: %d\n", n_duplicate_edges);

    return(n_edges);
}

/*
 * Calculate the mean of the closest distances between mesh1 and mesh2 and vice versa (Tfs).
 */
/**
 * \brief Compute distances between meshes or vertices.
 *
 * Function: compute_point_distance_mean
 *
 * For each vertex, computes the average of distances from p1[i]→mesh2 and p2[i]→mesh1,
 * then averages over all vertices.
 * \param p (polygons_struct *)
 * \param p2 (polygons_struct *)
 * \param dist (double *)
 * \param verbose (int)
 * \return See function description for return value semantics.
 */
double
compute_point_distance_mean(polygons_struct *p, polygons_struct *p2, double *dist, int verbose)
{
    int i, poly;
    Point closest, closest2;
    double avg_dist = 0.0, dist_tmp;

    if (p->n_points != p2->n_points) {
        printf("ERROR: Sizes of surfaces don't match\n");
        return(avg_dist);
    }
    
    create_polygons_bintree(p,  ROUND((double)  p->n_items * 0.5));
    create_polygons_bintree(p2, ROUND((double) p2->n_items * 0.5));

    /* walk through the points */
    for (i = 0; i < p->n_points; i++) {
        dist_tmp = 0.0;
        
        /* find closest distance between mesh1 and mesh2 and vice versa */
        poly = find_closest_polygon_point(&p->points[i], p2, &closest);
        poly = find_closest_polygon_point(&p2->points[i], p, &closest2);
        dist_tmp += distance_between_points(&p->points[i],  &closest);
        dist_tmp += distance_between_points(&p2->points[i], &closest2);
        
        /* calculate average of both distances */
        dist[i] = dist_tmp/2.0;
        avg_dist += dist[i];
    }

    avg_dist /= p2->n_points;
    if (verbose) 
        printf("Mean of closest distance: %f\n", avg_dist);
        
    delete_the_bintree(&p->bintree);
    delete_the_bintree(&p2->bintree);
    return(avg_dist);
}

/*
 * Calculate the exact Hausdorff distance.  This assumes that the two
 * input meshes are the same size and of the same brain.
 */
/**
 * \brief Compute Hausdorff distance between meshes.
 *
 * Function: compute_exact_hausdorff
 *
 * Uses direct point-wise distances; returns the maximum (Hausdorff) and prints means if requested.
 * \param p (polygons_struct *)
 * \param p2 (polygons_struct *)
 * \param hd (double *)
 * \param verbose (int)
 * \return See function description for return value semantics.
 */
double
compute_exact_hausdorff(polygons_struct *p, polygons_struct *p2, double *hd, int verbose)
{
    int i;
    double max_hd = 0.0, avg_hd = 0.0;

    /* walk through the points */
    for (i = 0; i < p->n_points; i++) {
        hd[i] = distance_between_points(&p->points[i], &p2->points[i]);

        if (hd[i] > max_hd)
            max_hd = hd[i];
        avg_hd += hd[i];
    }

    if (verbose) {
        avg_hd /= p->n_points;
        printf("Hausdorff distance: %f\n", max_hd);
        printf("Mean distance error: %f\n", avg_hd);
    }

    return(max_hd);
}

/*
 * Calculate the Hausdorff distance using mesh points only.
 */
/**
 * \brief Compute Hausdorff distance between meshes.
 *
 * Function: compute_point_hausdorff
 *
 * Computes directed Hausdorff from \c p→\c p2 using nearest vertices, then the reverse.
 * \param p (polygons_struct *)
 * \param p2 (polygons_struct *)
 * \param hd (double *)
 * \param verbose (int)
 * \return See function description for return value semantics.
 */
double
compute_point_hausdorff(polygons_struct *p, polygons_struct *p2, double *hd, int verbose)
{
    int i, poly;
    double *revhd;
    double max_hd = 0.0, avg_hd = 0.0;
    double max_revhd = 0.0, avg_revhd = 0.0;
    Point closest;

    create_polygons_bintree(p, ROUND((double) p->n_items * 0.5));
    create_polygons_bintree(p2, ROUND((double) p2->n_items * 0.5));

    /* walk through the points */
    for (i = 0; i < p->n_points; i++) {
        poly = find_closest_polygon_point(&p->points[i], p2, &closest);
        hd[i] = distance_between_points(&p->points[i], &closest);

        if (hd[i] > max_hd)
            max_hd = hd[i];
        avg_hd += hd[i];
    }

    avg_hd /= p->n_points;
    if (verbose) {
        printf("Hausdorff distance: %f\n", max_hd);
        printf("Mean distance error: %f\n", avg_hd);
    }
    
    /* Calculate the reverse Hausdorff */

    revhd = SAFE_MALLOC(double, p2->n_points);

    for (i = 0; i < p2->n_points; i++) {
        poly = find_closest_polygon_point(&p2->points[i], p, &closest);
        revhd[i] = distance_between_points(&p2->points[i], &closest);

        if (revhd[i] > max_revhd)
            max_revhd = revhd[i];
        avg_revhd += revhd[i];
    }

    avg_revhd /= p2->n_points;
    if (verbose) {
        printf("Reverse Hausdorff distance: %f\n", max_revhd);
        printf("Mean reverse distance error: %f\n", avg_revhd);
    }

    delete_the_bintree(&p->bintree);
    delete_the_bintree(&p2->bintree);
    free(revhd);
    return(max_hd);
}

/*
 * Calculate the closest distance using mesh points only from mesh1 to mesh2.
 */
/**
 * \brief Compute distances between meshes or vertices.
 *
 * Function: compute_point_distance
 *
 * \param p (polygons_struct *)
 * \param p2 (polygons_struct *)
 * \param hd (double *)
 * \param verbose (int)
 * \return See function description for return value semantics.
 */
double
compute_point_distance(polygons_struct *p, polygons_struct *p2, double *hd, int verbose)
{
    int i, poly;
    double avg_dist = 0.0;
    Point closest;

    create_polygons_bintree(p2, ROUND((double) p2->n_items * 0.5));

    /* walk through the points */
    for (i = 0; i < p->n_points; i++) {
        poly = find_closest_polygon_point(&p->points[i], p2, &closest);
        hd[i] = distance_between_points(&p->points[i], &closest);
        avg_dist += hd[i];
    }

    avg_dist /= p2->n_points;
    if (verbose) 
        printf("Mean of closest distance: %f\n", avg_dist);
        
    return(avg_dist);
}

/**
 * \brief Euler characteristic \f$\chi = F + V - E\f$ computed from point neighbourhoods.
 *
 * Function: euler_characteristic
 *
 * \param polygons (polygons_struct *)
 * \return See function description for return value semantics.
 */
int
euler_characteristic(polygons_struct *polygons, int verbose)
{
    int n_edges, *n_neighbours, **neighbours;

    create_polygon_point_neighbours(polygons, TRUE, &n_neighbours,
                    &neighbours, NULL, NULL);
    n_edges = count_edges(polygons, n_neighbours, neighbours, verbose);
  
    delete_polygon_point_neighbours(polygons, n_neighbours,
                    neighbours, NULL, NULL);

    create_polygon_point_neighbours(polygons, FALSE, &n_neighbours,
                    &neighbours, NULL, NULL);

    n_edges = count_edges(polygons, n_neighbours, neighbours, verbose);

    return(polygons->n_items + polygons->n_points - n_edges);
}

/**
 * \brief Convert or project a mesh to spherical form.
 *
 * Function: convert_ellipsoid_to_sphere_with_surface_area
 *
 * Scales coordinates so that the resulting sphere has area \c desiredSurfaceArea.
 * \param polygons (polygons_struct *)
 * \param desiredSurfaceArea (double)
 */
void
convert_ellipsoid_to_sphere_with_surface_area(polygons_struct *polygons,
                        double desiredSurfaceArea)
{
    int i;
    double radius, A, B, C, A2, B2, C2;
    double t1, t2, t3, f;
    double bounds[6];
    double xyz[3] = { 0.0, 0.0, 0.0 };
   
    /* Find radius of output sphere.  Sphere SA = 4 * PI * r * r */
    radius = sqrt(desiredSurfaceArea / (4.0 * PI));
   
    get_bounds(polygons, bounds); /* Determine lengths of the axes */

    A = (fabs(bounds[0]) + fabs(bounds[1])) * 0.5;
    B = (fabs(bounds[2]) + fabs(bounds[3])) * 0.5;
    C = (fabs(bounds[4]) + fabs(bounds[5])) * 0.5;

    /* Convert the coordinates from ellipsoid to sphere */
    A2 = A * A;
    B2 = B * B;
    C2 = C * C;

    for (i = 0; i < polygons->n_points; i++) {
        to_array(&polygons->points[i], xyz);
        /*
         *  ellipsoidal coordinates
         *
         *  x*x   y*y z*z
         *  --- + --- + --- = 1.0
         *  A*A   B*B C*C
         */
        t1 = (xyz[0]*xyz[0]) / (A2);
        t2 = (xyz[1]*xyz[1]) / (B2);
        t3 = (xyz[2]*xyz[2]) / (C2);
        f = sqrt(t1 + t2 + t3);

        if (f != 0.0) {
            xyz[0] /= f;
            xyz[1] /= f;
            xyz[2] /= f;
        }
     
        /* Push coordinate onto the sphere */
        xyz[0] = (radius * xyz[0]) / A;
        xyz[1] = (radius * xyz[1]) / B;
        xyz[2] = (radius * xyz[2]) / C;
        from_array(xyz, &polygons->points[i]);
    }
}

/**
 * \brief Linear (umbrella) smoothing with optional edge-only passes and periodic spherical projection.
 *
 * \param strength in (0,1]; larger moves points more towards neighbour averages.
 * \param iters    number of iterations.
 * \param smoothEdgesEveryXIters apply smoothing every X iterations (0 disables).
 * \param smoothOnlyTheseNodes   optional mask (length \c n_points) of nodes to smooth when active.
 * \param projectToSphereEveryXIters project to sphere of current radius every X iterations (0 disables).
 */
void
linear_smoothing(polygons_struct *polygons, double strength, int iters,
         int smoothEdgesEveryXIters, int *smoothOnlyTheseNodes,
         int projectToSphereEveryXIters)
{
    int i, j, k, l;
    int *n_neighbours, **neighbours;
    int pidx;
    BOOLEAN smoothSubsetOfNodes = 0;
    const double invstr = 1.0 - strength;
    double  xyz[3], pt[3];
    double radius = get_sphere_radius(polygons);
  
    create_polygon_point_neighbours(polygons, TRUE, &n_neighbours,
                    &neighbours, NULL, NULL);
  
    if (smoothOnlyTheseNodes != NULL) {
        smoothSubsetOfNodes = 1;
    }

    for (k = 1; k < iters; k++) {
        /* Should edges be smoothed? */
        BOOLEAN smoothEdges = 0;
        if (smoothEdgesEveryXIters > 0) {
            if ((k % smoothEdgesEveryXIters) == 0)
                smoothEdges = 1;
        }

        for (i = 0; i < polygons->n_points; i++) {
            BOOLEAN smoothIt = smoothEdges;
            if (smoothIt && smoothSubsetOfNodes) {
                smoothIt = (smoothOnlyTheseNodes)[i];
                if (!smoothIt) continue;
            }
            if (!smoothIt || n_neighbours[i] <= 0)
                continue; /* skip this point */
    
            for (j = 0; j < 3; j++)
                xyz[j] = 0.0;
            for (j = 0; j < n_neighbours[i]; j++) {
                pidx = neighbours[i][j];
                to_array(&polygons->points[pidx], pt);
                xyz[0] += pt[0];
                xyz[1] += pt[1];
                xyz[2] += pt[2];
            }
            /* Update the nodes position */
            to_array(&polygons->points[i], pt);
            for (l = 0; l < 3; l++) {
                pt[l] = (pt[l] * invstr) +
                    (xyz[l] / (double) n_neighbours[i] *
                     strength);
            }
            from_array(pt, &polygons->points[i]);
        }
    
        /* If the surface should be projected to a sphere */
        if (projectToSphereEveryXIters > 0) {
            if ((k % projectToSphereEveryXIters) == 0)
                for (i = 0; i < polygons->n_points; i++) 
                    set_vector_length(&polygons->points[i],
                              radius);
        }
    }
    delete_polygon_point_neighbours(polygons, n_neighbours,
                    neighbours, NULL, NULL);
}

/**
 * \brief Areal smoothing: neighbour influence weighted by local triangle areas around each node.
 *
 * Function: areal_smoothing
 *
 * See \c linear_smoothing for parameters; here, weights are proportional to incident tile areas.
 * \param polygons (polygons_struct *)
 * \param strength (double)
 * \param iters (int)
 * \param smoothEdgesEveryXIters (int)
 * \param smoothOnlyTheseNodes (int *)
 * \param projectToSphereEveryXIters (int)
 */
void
areal_smoothing(polygons_struct *polygons, double strength, int iters,
        int smoothEdgesEveryXIters, int *smoothOnlyTheseNodes,
        int projectToSphereEveryXIters)
{
    int i, j, k, l;
    int n1, n2, next;
    int *n_neighbours, **neighbours;
    double *area_values;
    Point pts[1000];
    double tileAreas[32], tileCenters[32*3];
    double xyz[3], pt1[3], pt2[3], pt3[3];
    double totalArea, weight;
    BOOLEAN smoothSubsetOfNodes = 0;
    BOOLEAN smoothEdges, smoothIt;
    double invstr = 1.0 - strength;
    double radius = get_sphere_radius(polygons);

    create_polygon_point_neighbours(polygons, TRUE, &n_neighbours,
                    &neighbours, NULL, NULL);
  
    if (smoothOnlyTheseNodes != NULL)
        smoothSubsetOfNodes = 1;

    for (k = 1; k < iters; k++) {
        /* should edges be smoothed? */
        smoothEdges = 0;
        if (smoothEdgesEveryXIters > 0) {
            if ((k % smoothEdgesEveryXIters) == 0)
                smoothEdges = 1;
        }

        for (i = 0; i < polygons->n_points; i++) {
            smoothIt = smoothEdges;
            if (smoothIt && smoothSubsetOfNodes)
                smoothIt = (smoothOnlyTheseNodes)[i];
    
            if (!smoothIt || n_neighbours[i] <= 1)
                continue; /* skip this point */

            totalArea = 0.0;  
          
            /* Get 2 consecutive neighbors of this node */
            for (j = 0; j < n_neighbours[i]; j++) {
                n1 = neighbours[i][j];
                next = j + 1;
                if (next >= n_neighbours[i])
                    next = 0;
                n2 = neighbours[i][next];
      
                /* Area of the triangle */
                pts[0] = polygons->points[i];
                pts[1] = polygons->points[n1];
                pts[2] = polygons->points[n2];
                tileAreas[j] = get_polygon_surface_area(3, pts);
                
                if (isnan(tileAreas[j])) tileAreas[j] = 0.0;
                tileAreas[j] = MIN(tileAreas[j], 1.0);
                tileAreas[j] = MAX(tileAreas[j], 0.0);
                totalArea += tileAreas[j];

                /* Save center of this tile */
                to_array(&polygons->points[i], pt1);
                to_array(&polygons->points[n1], pt2);
                to_array(&polygons->points[n2], pt3);
                for (l = 0; l < 3; l++) {
                    tileCenters[j*3+l] = (pt1[l] +
                                pt2[l] +
                                pt3[l]) / 3.0;
                }
            }

            /* Compute the influence of the neighboring nodes */
            for (j = 0; j < 3; j++)
                xyz[j] = 0.0;
            for (j = 0; j < n_neighbours[i]; j++) {
                if (tileAreas[j] > 0.0) {
                    weight = tileAreas[j] / totalArea;
                    for (l = 0; l < 3; l++) {
                        xyz[l] += weight *
                              tileCenters[j*3+l];
                    }
                }
            }
            /* Update the nodes position */
            to_array(&polygons->points[i], pt1);
            for (l = 0; l < 3; l++) {
                pt1[l] = (pt1[l] * invstr) +
                     (xyz[l] * strength);
            }
            from_array(pt1, &polygons->points[i]);
        }
    
        /* If the surface should be projected to a sphere */
        if (projectToSphereEveryXIters > 0) {
            if ((k % projectToSphereEveryXIters) == 0) {
                for (i = 0; i < polygons->n_points; i++) {
                    set_vector_length(&polygons->points[i],
                              radius);
                }
            }
        }
    }
    delete_polygon_point_neighbours(polygons, n_neighbours,
                    neighbours, NULL, NULL);
}

/* Use faster Manhattan distance */
/**
 * \brief Compute distances between meshes or vertices.
 *
 * Function: manhattan_distance_between_points
 *
 * \param p1 (Point *)
 * \param p2 (Point *)
 * \return See function description for return value semantics.
 */
double manhattan_distance_between_points(Point *p1, Point *p2) {
    double dx = fabs(Point_x(*p2) - Point_x(*p1));
    double dy = fabs(Point_y(*p2) - Point_y(*p1));
    double dz = fabs(Point_z(*p2) - Point_z(*p1));
    return dx + dy + dz;
}

/**
 * \brief Surface smoothing on polygon meshes.
 *
 * Function: distance_smoothing
 *
 * \param polygons (polygons_struct *)
 * \param strength (double)
 * \param iters (int)
 * \param smoothEdgesEveryXIters (int)
 * \param smoothOnlyTheseNodes (int *)
 * \param projectToSphereEveryXIters (int)
 */
void
distance_smoothing(polygons_struct *polygons, double strength, int iters,
           int smoothEdgesEveryXIters, int *smoothOnlyTheseNodes,
           int projectToSphereEveryXIters)
{
    int i, j, k, l, pidx;
    int *n_neighbours, **neighbours;
    double tileDist[MAX_POINTS_PER_POLYGON];
    double pt[3];
    double weight;
    double invstr = 1.0 - strength;
    double radius = get_sphere_radius(polygons);
    BOOLEAN smoothSubsetOfNodes = (smoothOnlyTheseNodes != NULL);
  
    create_polygon_point_neighbours(polygons, TRUE, &n_neighbours,
                    &neighbours, NULL, NULL);

    for (k = 1; k < iters; k++) {
        /* should the edges be smoothed? */
        BOOLEAN smoothEdges = (smoothEdgesEveryXIters > 0 && (k % smoothEdgesEveryXIters) == 0);

        for (i = 0; i < polygons->n_points; i++) {
            BOOLEAN smoothIt = smoothEdges;
            if (smoothIt && smoothSubsetOfNodes)
                smoothIt = (smoothOnlyTheseNodes)[i];
    
            if (!smoothIt || n_neighbours[i] <= 1)
                continue; /* skip this point */

            /* Compute the influence of neighboring nodes */
            double xyz[3] = {0.0, 0.0, 0.0};
            double totalDistance = 0.0;
            
            for (j = 0; j < n_neighbours[i]; j++) {  
                tileDist[j] = manhattan_distance_between_points(
                   &polygons->points[i],
                   &polygons->points[neighbours[i][j]]);
                totalDistance += tileDist[j];            
            }

            for (j = 0; j < n_neighbours[i]; j++) {
                if (tileDist[j] <= 0.0)
                    continue;

                pidx = neighbours[i][j];
                to_array(&polygons->points[pidx], pt);
                weight = tileDist[j] / totalDistance;
                for (l = 0; l < 3; l++)
                    xyz[l] += weight * pt[l];
            }

            /* Update the nodes position */
            to_array(&polygons->points[i], pt);
            for (l = 0; l < 3; l++) {
                pt[l] = (pt[l] * invstr) + (xyz[l] * strength);
            }
            from_array(pt, &polygons->points[i]);
        }
      
        /* If the surface should be projected to a sphere */
        if ((projectToSphereEveryXIters > 0) && ((k % projectToSphereEveryXIters) == 0)) {
            for (i = 0; i < polygons->n_points; i++)
                set_vector_length(&polygons->points[i], radius);
        }
    }
    
    delete_polygon_point_neighbours(polygons, n_neighbours,
                    neighbours, NULL, NULL);
}


/**
 * \brief Surface smoothing on polygon meshes.
 *
 * Function: inflate_surface_and_smooth_fingers
 *
 * Alternates displacement (inflation) and smoothing; computes metrics 
 * (linear distortion, areal compression) to target highly distorted nodes.
 * \param polygonsIn (polygons_struct *)
 * \param n_smoothingCycles (const int)
 * \param regSmoothStrength (const double)
 * \param regSmoothIters (const int)
 * \param inflationFactorIn (const double)
 * \param compStretchThresh (const double)
 * \param fingerSmoothStrength (const double)
 * \param fingerSmoothIters (const int)
 */
void
inflate_surface_and_smooth_fingers(polygons_struct *polygonsIn,
                   const int n_smoothingCycles,
                   const double regSmoothStrength,
                   const int regSmoothIters,
                   const double inflationFactorIn,
                   const double compStretchThresh,
                   const double fingerSmoothStrength,
                   const int fingerSmoothIters)
{

    polygons_struct *polygons;
    int i, j, cycle, n;
    int numDistortionAboveThresh;
    int *n_neighbours, **neighbours, *needSmoothing, nidx;
    const double inflationFactor = inflationFactorIn - 1.0;
    double tileArea, tileAreaIn, distort;    
    double bounds[6], xyz[3], nodept[3], nodeptIn[3];
    double *area_values, *area_valuesIn, diff_bound[3];
    double x, y, z, r, k;
    double numNeighbors, neighpt[3], neighptIn[3];
    float *avgCompStretch, *compStretch;
    float *maxLinDistort, *avgArealComp;
    float *stretching;
    float SA, SA_ratio, inflatedSA;
    float maxDistort, minDistort, ratio;
    float dx, dy, dz, dist, distIn;
    object_struct *out_object;
  
    /* Copy the fiducial surface since it will be modified 
       (translated to center of mass) */
    out_object = create_object(POLYGONS);
    polygons = get_polygons_ptr(out_object);
    copy_polygons(polygonsIn, polygons);
  
    /* Get bounds of fiducial surface */
    get_bounds(polygons, bounds);

    for (j = 0; j < 3; j++) 
        diff_bound[j] = bounds[2*j+1] - bounds[2*j];

    /* Translate the fiducial to center of mass */
    translate_to_center_of_mass(polygonsIn);

    SA = get_polygons_surface_area(polygons);
    
    avgCompStretch = SAFE_MALLOC(float, polygons->n_points);
    maxLinDistort  = SAFE_MALLOC(float, polygons->n_points);
    avgArealComp   = SAFE_MALLOC(float, polygons->n_points);
    compStretch    = SAFE_MALLOC(float, polygons->n_points);
    stretching     = SAFE_MALLOC(float, polygons->n_points);
    area_values    = SAFE_MALLOC(double, polygons->n_points);
    area_valuesIn  = SAFE_MALLOC(double, polygons->n_points);
    needSmoothing  = SAFE_MALLOC(int, polygons->n_points);
    SA_ratio = 0.0;

    for (cycle = 0; cycle < (n_smoothingCycles + 1); cycle++) {
        if (cycle < n_smoothingCycles) {
            /* Step 6a: Apply Smoothing to AUX coord */
            distance_smoothing(polygonsIn, regSmoothStrength, 
                    regSmoothIters, 1, NULL, 0);

            /* Step 6b: Incrementally inflate AUX surface by */
            /*      Ellipsoidal Projection  */
            for (i = 0; i < polygons->n_points; i++) {
                to_array(&polygonsIn->points[i], xyz);
        
                x = xyz[0] / diff_bound[0];
                y = xyz[1] / diff_bound[1];
                z = xyz[2] / diff_bound[2];

                r = sqrt(x*x + y*y + z*z);

                k = 1.0 + inflationFactor * (1.0 - r);
                xyz[0] *= k;
                xyz[1] *= k;
                xyz[2] *= k;
                from_array(xyz, &polygonsIn->points[i]);
            }
        }
    
        /* Step 6c: Calculate surface area of this surface */
        inflatedSA = get_polygons_surface_area(polygonsIn);
    
        /* Ratio of inflated and spherical surfaces */
        SA_ratio = inflatedSA / SA;

        create_polygon_point_neighbours(polygons, TRUE, &n_neighbours,
                        &neighbours, NULL, NULL);

        get_area_of_points(polygons, area_values);
        get_area_of_points(polygonsIn, area_valuesIn);

        // Step 6d: Calculate compress/stretched value for each node
        for (i = 0; i < polygons->n_points; i++) {
            /* Get node pos. in both this and fiducial surface */
            to_array(&polygonsIn->points[i], nodeptIn);
            to_array(&polygons->points[i], nodept);
     
            maxLinDistort[i] = 0.0;
            avgArealComp[i] = 0.0;
     
            numNeighbors = 0;
     
            /* Loop through the neighbors */
            for (j = 0; j < n_neighbours[i]; j++) {
                nidx = neighbours[i][j];
                to_array(&polygonsIn->points[nidx], neighptIn); 
                to_array(&polygons->points[nidx], neighpt); 
      
                /* max lin. distortion on aux & ref surface */
                dx = neighptIn[0] - nodeptIn[0];
                dy = neighptIn[1] - nodeptIn[1];
                dz = neighptIn[2] - nodeptIn[2];
                distIn = sqrt(dx*dx + dy*dy + dz*dz);

                dx = neighpt[0] - nodept[0];
                dy = neighpt[1] - nodept[1];
                dz = neighpt[2] - nodept[2];

                dist = sqrt(dx*dx + dy*dy + dz*dz);
                if (dist > 0.0) {
                    ratio = (float)distIn / dist;
                    if (ratio > maxLinDistort[i]) {
                        maxLinDistort[i] = ratio;
                    }
                }

                /* area of tiles on aux and ref surfaces */
                tileAreaIn = area_valuesIn[i];
                tileArea = area_values[i];
                 
                /* avg compression of tiles assoc w/ node */
                distort = 0.0;
                if (tileAreaIn > 0.0) {
                    distort = tileArea / tileAreaIn;
                } else {
                    if (tileArea != 0.0) {
                        /* big since denominator =0 */
                        distort = 10000.0;
                    } else {
                        /* if both 0 then assume */
                        /* same area */
                        distort = 1.0;
                    }
                }

                /* Zero will cause -inf */
                if (distort < 0.00000001) {
                    distort = 0.00000001;
                }

                avgArealComp[i] += (float)distort;
                numNeighbors += 1.0;
            }

            if (numNeighbors > 0) {
                avgArealComp[i] /= numNeighbors;
            }

            /* compressed/stretched for node */
            compStretch[i] = maxLinDistort[i] * avgArealComp[i] *
                     SA_ratio;

            /* stretching for node */
            stretching[i] = maxLinDistort[i] *
                    sqrt(avgArealComp[i] * SA_ratio);
        }

        /* avg comp/stretch for all nodes by averaging w/ neighbors */
        for (i = 0; i < polygons->n_points; i++) {
            avgCompStretch[i] = compStretch[i];

            if (n_neighbours[i] > 0) {
                for (j = 0; j < n_neighbours[i]; j++) {
                    n = neighbours[i][j];
                    avgCompStretch[i] += compStretch[n];
                }
                avgCompStretch[i] /= (float)
                            (n_neighbours[i] + 1);
            }
        }

        /* Step 6e: Flag highly compressed/stretched nodes for 
           targeted smoothing */
        numDistortionAboveThresh = 0;
        maxDistort = -FLT_MAX;
        minDistort =  FLT_MAX;
    
        for (i = 0; i < polygons->n_points; i++) {
            if (avgCompStretch[i] > compStretchThresh) {
                numDistortionAboveThresh++;
                needSmoothing[i] = 1;
            } else needSmoothing[i] = 0;
            if (avgCompStretch[i] > maxDistort) 
                maxDistort = avgCompStretch[i];
            if (avgCompStretch[i] < minDistort) 
                minDistort = avgCompStretch[i];
        }

        if (cycle < n_smoothingCycles) {
            /* Step 6f: Targeted smoothing */
            distance_smoothing(polygonsIn,
                    fingerSmoothStrength,
                    fingerSmoothIters,
                    1, needSmoothing, 0);
        }
    }
    
    compute_polygon_normals(polygonsIn);
    
    delete_object(out_object);
    free(area_values);
    free(area_valuesIn);
    free(avgCompStretch);
    free(maxLinDistort);
    free(avgArealComp);
    free(compStretch);
    free(stretching);
    free(needSmoothing);
    delete_polygon_point_neighbours(polygons, n_neighbours,
                    neighbours, NULL, NULL);
}

/**
 * \brief Multi-stage pipeline to convert a mesh into (increasingly smoothed/inflated) spherical form.
 *
 * Stages: low smooth → inflate → very inflate → high smooth → ellipsoid projection → optional areal smoothing.
 * \param stop_at stage index (1..5) controlling how far to proceed.
 * \param verbose print stage info and iteration scaling for large meshes.
 */
void
surf_to_sphere(polygons_struct *polygons, int stop_at, int verbose)
{
    BOOLEAN enableFingerSmoothing = 1;
    int fingerSmoothingIters, arealSmoothingIters;
    double surfarea, factor;
    
    surfarea = get_polygons_surface_area(polygons);

    /* use more iterations for larger surfaces */
    if (polygons->n_items > 350000) {
        factor = (double)polygons->n_items/350000.0;
        if (verbose)
            fprintf(stderr, "Large number polygons -> Increase # of iterations by factor %g.\n", factor);
    } else factor = 1.0;

    /* low smooth */
    if (verbose) printf("Low Smooth\n");
    inflate_surface_and_smooth_fingers(polygons,
          /*           cycles */ 1,
          /* regular smoothing strength */ 0.2,
          /*  regular smoothing iters */ round(factor*50),
          /*       inflation factor */ 1.0,
          /* finger comp/stretch thresh */ 3.0,
          /*   finger smooth strength */ 1.0,
          /*    finger smooth iters */ 0);

    if (stop_at > 1) {
      if (verbose) printf("Inflate\n");
        /* inflated */
        fingerSmoothingIters = 0;
        if (enableFingerSmoothing)
            fingerSmoothingIters = 30;
        inflate_surface_and_smooth_fingers(polygons,
          /*           cycles */ 2,
          /* regular smoothing strength */ 1.0,
          /*  regular smoothing iters */ round(factor*30),
          /*       inflation factor */ 1.4,
          /* finger comp/stretch thresh */ 3.0,
          /*   finger smooth strength */ 1.0,
          /*    finger smooth iters */ fingerSmoothingIters);
    }                       
  
    if (stop_at > 2) {
      if (verbose) printf("Very inflate\n");
        /* very inflated */
        inflate_surface_and_smooth_fingers(polygons,
          /*           cycles */ 4,
          /* regular smoothing strength */ 1.0,
          /*  regular smoothing iters */ round(factor*30),
          /*       inflation factor */ 1.1,
          /* finger comp/stretch thresh */ 3.0,
          /*   finger smooth strength */ 1.0,
          /*    finger smooth iters */ 0);
    }                       
  
    if (stop_at > 3) {
      if (verbose) printf("High smooth\n");
        /* high smooth */
        fingerSmoothingIters = 0;
        if (enableFingerSmoothing)
            fingerSmoothingIters = 60;
        inflate_surface_and_smooth_fingers(polygons,
          /*           cycles */ 6,
          /* regular smoothing strength */ 1.0,
          /*  regular smoothing iters */ round(factor*60),
          /*       inflation factor */ 1.6,
          /* finger comp/stretch thresh */ 3.0,
          /*   finger smooth strength */ 1.0,
          /*    finger smooth iters */ fingerSmoothingIters);
    }

    if (stop_at > 4) {
      if (verbose) printf("Ellipsoid\n");
        /* ellipsoid */
        inflate_surface_and_smooth_fingers(polygons,
          /*           cycles */ 6,
          /* regular smoothing strength */ 1.0,
          /*  regular smoothing iters */ round(factor*50),
          /*       inflation factor */ 1.4,
          /* finger comp/stretch thresh */ 4.0,
          /*   finger smooth strength */ 1.0,
          /*    finger smooth iters */ fingerSmoothingIters);
        convert_ellipsoid_to_sphere_with_surface_area(polygons,
                                surfarea);
    }
  
    if (stop_at > 5) {
      if (verbose) printf("Areal smoothing\n");
        /* areal smoothing */
        arealSmoothingIters = 1000*(stop_at - 5);
        areal_smoothing(polygons, 1.0, arealSmoothingIters, 1, NULL, 1000);
        convert_ellipsoid_to_sphere_with_surface_area(polygons,
                                surfarea);
    }

    compute_polygon_normals(polygons);
}


/**
 * \brief Check counter‑clockwise neighbour ordering.
 *
 * Function: ccw_neighbours
 *
 * Increments \c point_error for neighbours involved in local ordering inversions.
 * \param centroid (Point *)
 * \param normal (Vector *)
 * \param points (Point)
 * \param n_nb (int)
 * \param neighbours (int)
 * \param point_error (signed char)
 * \return See function description for return value semantics.
 */
BOOLEAN
ccw_neighbours(Point *centroid, Vector *normal, Point points[],
         int n_nb, int neighbours[], signed char point_error[])
{
    Vector to_nb, prev_to_nb, up, offset;
    double len;
    int i;
    BOOLEAN ccw;

    ccw = TRUE;
    fill_Vector(to_nb, 0.0, 0.0, 0.0);

    for (i = 0; i < n_nb + 1; i++) {
        prev_to_nb = to_nb;
        SUB_VECTORS(to_nb, points[neighbours[i % n_nb]], *centroid);
        len = DOT_VECTORS(to_nb, *normal);
        SCALE_VECTOR(offset, *normal, len);
        SUB_VECTORS(to_nb, to_nb, offset);

        if (i != 0) {
            CROSS_VECTORS(up, prev_to_nb, to_nb);
            if (DOT_VECTORS(up, *normal) < 0.0) {
                ++point_error[neighbours[i % n_nb]];
                ++point_error[neighbours[(i-1 + n_nb) % n_nb]];
                ccw = FALSE;
            }
        }
    }

    return(ccw);
}

/**
 * \brief Mesh geometry utility.
 *
 * Function: check_polygons_shape_integrity
 *
 * Uses neighbour ordering checks around each vertex. Vertices implicated in 
 * inversions are snapped to their polygon-point centroid estimate.
 * \param polygons (polygons_struct *)
 * \param new_points (Point)
 */
void
check_polygons_shape_integrity(polygons_struct *polygons, Point new_points[])
{
    int vertidx, ptidx, poly, size;
    int n_nb, neighbours[MAX_NEIGHBOURS];
    int n_errors, n_bad_points;
    double base_length, curv_factor;
    Vector normal;
    BOOLEAN interior_flag;

    signed char *point_done   = SAFE_MALLOC(signed char, polygons->n_points);
    signed char *point_error  = SAFE_MALLOC(signed char, polygons->n_points);
    Point       *centroids    = SAFE_MALLOC(Point, polygons->n_points);

    for (ptidx = 0; ptidx <  polygons->n_points; ptidx++) {
        point_done[ptidx] = FALSE;
        point_error[ptidx] = 0;
    }

    for (poly = 0; poly < polygons->n_items; poly++) {
        size = GET_OBJECT_SIZE(*polygons, poly);

        for (vertidx = 0; vertidx <  size; vertidx++) {
            ptidx = polygons->indices[
              POINT_INDEX(polygons->end_indices, poly, vertidx)];

            if (!point_done[ptidx]) {
                point_done[ptidx] = TRUE;

                compute_polygon_point_centroid(polygons, poly,
                                vertidx, ptidx,
                                &centroids[ptidx],
                                &normal,
                                &base_length,
                                &curv_factor);
                n_nb = get_neighbours_of_point(polygons, poly,
                                 vertidx,
                                 neighbours,
                                 MAX_NEIGHBOURS,
                                 &interior_flag);

#ifdef CHECK_CLOCKWISE_NEIGHBOURS
                if (!ccw_neighbours(&centroids[ptidx], &normal,
                          new_points, n_nb,
                          neighbours, point_error)) {
                    point_error[ptidx] = TRUE;
                    ++n_errors;
                }
#else
                ccw_neighbours(&centroids[ptidx], &normal,
                         new_points, n_nb,
                         neighbours, point_error);
#endif
            }
        }
    }


#ifdef DEBUG
    n_errors = 0;
    n_bad_points = 0;
#endif
    for (ptidx = 0; ptidx < polygons->n_points; ptidx++) {
        if (point_error[ptidx] > 0) {
#ifdef DEBUG
            ++n_errors;
            n_bad_points += point_error[ptidx];
            if (n_errors < 10)
                printf(" %d", ptidx);
#endif
            new_points[ptidx] = centroids[ptidx];
        }
    }

#ifdef DEBUG
    if (n_errors > 0)
        printf(": Shape errors %d/%d\n", n_errors, n_bad_points);
#endif

    free(point_error);
    free(centroids);
    free(point_done);
}

/*
 * Calls central_to_pial, but creates a new surface object and does not modify the original surface
 * The direct implementation into central_to_pial was not working because of issues with the function 
 * check_polygons_shape_integrity.
 */
/**
 * \brief Create a new pial/white surface from a central surface without modifying the input.
 *
 * Thin wrapper around \c central_to_pial that works on a copy (due to integrity checks).
 * \return new object array containing a single POLYGONS object.
 */
object_struct **
central_to_new_pial(polygons_struct *polygons, double *thickness_values, double *extents, int check_intersects, double sigma, int iterations, int verbose)
{
    polygons_struct *polygons_out;
    object_struct **objects_out;
    
    objects_out  = SAFE_MALLOC(object_struct *, 1);
    *objects_out = create_object(POLYGONS);
    polygons_out = get_polygons_ptr(*objects_out);
    
    copy_polygons(polygons, polygons_out);
    central_to_pial(polygons_out, thickness_values, extents, check_intersects, sigma, iterations, verbose);
    
    return(objects_out);
}

/*
 * Estimate pial surface from central surface using cortical thickness values. In order to estimate the pial surface 
 * an extent of 0.5 should be used, while an extent of -0.5 results in the estimation of the white matter surface.
 */
/**
 * \brief Displace a central surface along normals to estimate pial/white surfaces from thickness.
 *
 * Performs 20 incremental steps with optional smoothing of the displacement field; checks and
 * corrects integrity at each step; optionally removes near self-intersections.
 *
 * \param extents signed per-vertex factors; +0.5 ~ pial, −0.5 ~ white.
 * \param sigma   Gaussian smoothing (if >0) for the displacement field.
 * \param iterations number of smoothing iterations (used when \c sigma>0).
 */
void
central_to_pial(polygons_struct *polygons, double *thickness_values, double *extents, 
              int check_intersects, double sigma, int iterations, int verbose)
{
    int *n_neighbours, **neighbours;
    int i, p, n_steps, dims[3];
    double x, y, z, val, length;
    polygons_struct *polygons_out;
    Point *new_pts;
    object_struct **objects_out;

    compute_polygon_normals(polygons);
    check_polygons_neighbours_computed(polygons);

    objects_out  = SAFE_MALLOC(object_struct *, 1);
    *objects_out = create_object(POLYGONS);
    polygons_out = get_polygons_ptr(*objects_out);
            
    copy_polygons(polygons, polygons_out);
    if (((sigma > 0.0) && (iterations > 0)))
        create_polygon_point_neighbours(polygons, TRUE, &n_neighbours, &neighbours, NULL, NULL);

    /* use 20 steps to add thickness values to central surface and check in each step shape integrity and self intersections */
    n_steps = 20;
    length = 1.0;

    /* Allocate displacement field */
    double (*displacements)[3] = (double (*)[3]) SAFE_MALLOC(double, 3 * polygons->n_points);

    for (i = 0; i < n_steps; i++) {
        copy_polygons(polygons, polygons_out);
        
        /* decrease length of change for each step by ratio of 2 to achieve larger changes
           in the first steps which will get smaller from step to step */
        if ((i+1) < n_steps)
            length /= 2.0;

        /* Calculate displacements */
        for (p = 0; p < polygons->n_points; p++) {
            double factor = extents[p] * length * thickness_values[p];
            displacements[p][0] = factor * Point_x(polygons->normals[p]);
            displacements[p][1] = factor * Point_y(polygons->normals[p]);
            displacements[p][2] = factor * Point_z(polygons->normals[p]);
        }

        /* Smooth displacement vectors (optional: adjust iterations and smoothing strength) */
        if ((sigma > 0.0) && (iterations > 0))
            smooth_displacement_field(displacements, polygons, n_neighbours, 
                                      neighbours, iterations, sigma);

        /* Apply smoothed displacements */
        for (p = 0; p < polygons->n_points; p++) {
            Point_x(polygons_out->points[p]) += displacements[p][0];
            Point_y(polygons_out->points[p]) += displacements[p][1];
            Point_z(polygons_out->points[p]) += displacements[p][2];
        }

        new_pts = polygons_out->points;
                
        /* check polygon integrity of new points */
        check_polygons_shape_integrity(polygons, new_pts);
        polygons->points = new_pts;
    }

    compute_polygon_normals(polygons);

    /* final check and correction for self intersections */
    if (check_intersects)
        remove_near_intersections(polygons, 0.75, verbose);
        
    free(displacements);
    if (((sigma > 0.0) && (iterations > 0))) {
        delete_polygon_point_neighbours(polygons, n_neighbours, neighbours, NULL, NULL);
    }
}

/*
 * Estimate area for each point of a surface that is shifted along normals by thickness extent 
 * Extent of 0.5 for a central surface results in pial surface while -0.5 leads to a white surface
 */
/**
 * \brief Compute or return a derived quantity from the mesh.
 *
 * Function: get_area_of_points_central_to_pial
 *
 * \param polygons (polygons_struct *)
 * \param area (double *)
 * \param thickness_values (double *)
 * \param extent (double)
 * \return See function description for return value semantics.
 */
double
get_area_of_points_central_to_pial(polygons_struct *polygons, double *area, double *thickness_values, double extent)
{
    double surface_area, *extents;
    int p;
    polygons_struct *polygons_transformed;
    object_struct **objects_transformed;

    extents = SAFE_MALLOC(double, polygons->n_points);                
    for (p = 0; p < polygons->n_points; p++) extents[p] = extent;

    objects_transformed = central_to_new_pial(polygons, thickness_values, extents, 
        0, 0.0, 0, 0);
    polygons_transformed = get_polygons_ptr(objects_transformed[0]);
    surface_area = get_area_of_points(polygons_transformed, area);
    
    free(extents);

    return(surface_area);
}

/* objective function for finding the minimum of the difference to either a defined 
 * isovalue in a volume or to a reference mesh
*/
/**
 * \brief Compute or return a derived quantity from the mesh.
 *
 * Function: get_distance_mesh_correction
 *
 * Displaces points by \c weight*curvature*normal and measures either squared intensity error to
 * \c isovalue in \c vol, or point-wise distance to \c polygons_reference.
 * \param polygons (polygons_struct *)
 * \param polygons_reference (polygons_struct *)
 * \param vol (float*)
 * \param nii_ptr (nifti_image *)
 * \param isovalue (double)
 * \param curvatures (double *)
 * \param weight (double)
 * \return See function description for return value semantics.
 */
double
get_distance_mesh_correction(polygons_struct *polygons, polygons_struct *polygons_reference, float* vol,
    nifti_image *nii_ptr, double isovalue, double *curvatures, double weight)
{
    int dims[3];
    int p;
    double distance = 0.0, val, x, y, z;
    float *input;
    Point point;

    if (nii_ptr != NULL) {
        dims[0] = nii_ptr->nx;
        dims[1] = nii_ptr->ny;
        dims[2] = nii_ptr->nz;
    }
    
    for (p = 0; p < polygons->n_points; p++) {
        x = Point_x(polygons->points[p]) + weight*curvatures[p]*Point_x(polygons->normals[p]);
        y = Point_y(polygons->points[p]) + weight*curvatures[p]*Point_y(polygons->normals[p]);
        z = Point_z(polygons->points[p]) + weight*curvatures[p]*Point_z(polygons->normals[p]);
                    
        /* if nifti pointer is defined then use squared difference between isovalue and real value 
           at surface border otherwise use distance between the mesh points */
        if (nii_ptr != NULL) {
            val = (double)isoval(vol, x, y, z, dims, nii_ptr);
            if (!isnan(val)) distance += (val-isovalue)*(val-isovalue);
        } else {
            point = polygons->points[p];
            Point_x(point) += weight*curvatures[p]*Point_x(polygons->normals[p]);
            Point_y(point) += weight*curvatures[p]*Point_y(polygons->normals[p]);
            Point_z(point) += weight*curvatures[p]*Point_z(polygons->normals[p]);

            val = distance_between_points(&point,&polygons_reference->points[p]);
            if (!isnan(val)) distance  += val;
        }
    }
    return(distance);
}

/*
 * Correct mesh in folded areas to compensate for the averaging effect in gyri and sulci.
 * We use a folding measure (i.e. mean curvature averaged) to estimate the compensation. The amount
 * of compensation is automatically estimated using the difference to either a defined 
 * isovalue in a volume or to a reference mesh.
*/
/**
 * \brief Correct mesh properties or alignment.
 *
 * Function: correct_mesh_folding
 *
 * Estimates an optimal \c weight that minimizes an intensity/mesh distance objective, smooths curvature,
 * and applies the final displacement along vertex normals.
 * \param polygons (polygons_struct *)
 * \param polygons_reference (polygons_struct *)
 * \param vol (float *)
 * \param nii_ptr (nifti_image *)
 * \param isovalue (double)
 * \return See function description for return value semantics.
 */
int
correct_mesh_folding(polygons_struct *polygons, polygons_struct *polygons_reference, float *vol, 
     nifti_image *nii_ptr, double isovalue)
{
    int curvtype = 0; /* mean curvature averaged over 3mm, in degrees */
    int *n_neighbours, **neighbours, p;
    double distance, eps = 1e-6, fwhm = 2.0;
    double a, b, f1, f2;
    double weight = 0.0, weight1, weight2, avg;
    const double  phi = (1.0 + sqrt(5.0))/2.0; /* Golden ratio constant */
    
    if (!polygons_reference && !nii_ptr) {
      printf("ERROR: You have to define polygons_reference or volume.\n");
      exit(EXIT_FAILURE);
    }
    
    /* for curvtype 0 (mean curvature in degrees) use average around point */
    if (curvtype == 0) avg = 3.0; else avg = 0.0;
    
    double *curvatures = SAFE_MALLOC(double, polygons->n_points);
    compute_polygon_normals(polygons);
    get_all_polygon_point_neighbours(polygons, &n_neighbours, &neighbours);
    get_polygon_vertex_curvatures_cg(polygons, n_neighbours, neighbours,
                 avg, curvtype, curvatures);
    
    /* smooth curvature values to minimize noisy meshes */
    smooth_heatkernel(polygons, curvatures, fwhm);
    
    /* for curvtype 0 (mean curvature in degrees) we have to transform to radians */
    if (curvtype == 0) {
        for (p = 0; p < polygons->n_points; p++)
            curvatures[p] /= 57.3;
    }
    
    /* define lower and upper bound */
    a = -2.0; b = 2.0;
    
    weight1 = b - (b - a) / phi;
    f1 = get_distance_mesh_correction(polygons, polygons_reference, vol, nii_ptr, 
          isovalue, curvatures, weight1);
  
    weight2 = a + (b - a) / phi;
    f2 = get_distance_mesh_correction(polygons, polygons_reference, vol, nii_ptr, 
          isovalue, curvatures, weight2);
  
    /* Golden Ratio method for finding minimum */
    while (fabs(b - a) > eps) {
      if (f1 < f2) {
        b = weight2;
        weight2 = weight1;
        f2 = f1;
        weight1 = b - (b - a) / phi;
        f1 = get_distance_mesh_correction(polygons, polygons_reference, vol, nii_ptr, 
              isovalue, curvatures, weight1);
      } else {
        a = weight1;
        weight1 = weight2;
        f1 = f2;
        weight2 = a + (b - a) / phi;
        f2 = get_distance_mesh_correction(polygons, polygons_reference, vol, nii_ptr, 
              isovalue, curvatures, weight2);
      }
    }
  
    /* estimate weight using the new bounds */
    weight = (a + b)/2.0;          
    
    /* additionally smooth curvature values w.r.t. estimated weight to 
       minimize noisy meshes for larger weights */
    smooth_heatkernel(polygons, curvatures, fabs(2.0*weight));
    
    /* finally apply the optimized weighting to original mesh */
    for (p = 0; p < polygons->n_points; p++) {
        Point_x(polygons->points[p]) += weight*curvatures[p]*Point_x(polygons->normals[p]);
        Point_y(polygons->points[p]) += weight*curvatures[p]*Point_y(polygons->normals[p]);
        Point_z(polygons->points[p]) += weight*curvatures[p]*Point_z(polygons->normals[p]);          
    }

    free(curvatures);
    
    return(EXIT_SUCCESS);
}

/**
 * Reduce a triangle mesh using Quadric Error Metrics (QEM).
 *
 * This function mirrors the calling style used by spm_mesh_reduce.c: you specify
 * a target number of triangles and an aggressiveness value. Internally, the input
 * polygons are treated as triangles; non-triangle faces are fan-triangulated.
 *
 * Parameters
 * ----------
 * polygons        : (in/out) BICPL polygons_struct to be simplified in-place.
 * target_faces    : desired number of triangles after simplification. If <= 0,
 *                   half of the current triangle count is used.
 * aggressiveness  : simplifier aggressiveness (typical ~7.0; larger => stronger/rougher).
 * preserve_sharp  : if non-zero, prevent aggressive collapses across sharp edges.
 * verbose         : if non-zero, print progress and summary to stderr.
 *
 */
/**
 * \brief Reduce mesh complexity using Quadric Error Metrics (QEM).
 *
 * Function: reduce_mesh_quadrics
 *
 * \param polygons (polygons_struct *)
 * \param target_faces (int)
 * \param aggressiveness (double)
 * \param preserve_sharp (int)
 * \param verbose (int)
 * \return See function description for return value semantics.
 */
int reduce_mesh_quadrics(polygons_struct *polygons,
                         int target_faces,
                         double aggressiveness,
                         int preserve_sharp,
                         int verbose)
{
    if (!polygons || polygons->n_points <= 0 || polygons->n_items <= 0) return -1;

    vec3d *V = NULL;
    vec3i *F = NULL;
    int nv = 0, nf = 0, fan = 0;

    /* --- in: BICPL -> (V,F) triangles --- */
    if (polygons_to_tri_arrays(polygons, &V, &F, &nv, &nf, &fan) != 0)
        return -1;

    int target = qem_target(nf, target_faces);

    if (verbose) {
        fprintf(stderr, "[QEM] input: %d verts / %d tris (%s)\n",
                nv, nf, fan ? "fan-triangulated" : "tri");
        fprintf(stderr, "[QEM] target tris: %d, aggressiveness: %.3f, preserve_sharp: %d\n",
                target, aggressiveness, preserve_sharp);
    }

    /* --- simplify --- */
    quadric_simplify_mesh(&V, &F, &nv, &nf,
                          target,
                          aggressiveness,
                          0, /* extract_uv = false */
                          (preserve_sharp != 0) || (aggressiveness < 7.0));

    if (verbose) fprintf(stderr, "[QEM] output: %d verts / %d tris\n", nv, nf);
    if (nv <= 0 || nf <= 0) { free(V); free(F); return -1; }

    /* --- out: (V,F) -> BICPL (triangles, exclusive end_indices) --- */
    const int ok = tri_arrays_to_polygons(polygons, V, F, nv, nf);
    free(V);
    free(F);
    if (ok != 0) return -1;

    if (verbose) {
        double area = get_polygons_surface_area(polygons);
        fprintf(stderr, "[QEM] surface area after: %g\n", area);
    }
    return 0;
}
