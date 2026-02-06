/* Rachel Yotter - rachel.yotter@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 *
 */

#include <bicpl.h>

#include "CAT_Curvature.h"
#include "CAT_Smooth.h"
#include "CAT_Octree.h"
#include "CAT_Defect.h"
#include "CAT_Math.h"
#include "CAT_SurfaceIO.h"
#include "CAT_Intersect.h"

#define GRID_RES 250  // Number of grid cells per axis (adjustable)

typedef struct PointNode {
    int index;
    struct PointNode *next;
} PointNode;

typedef struct GridCell {
    PointNode *points;
} GridCell;

typedef struct SpatialGrid {
    GridCell *cells;
    int res;
    double cell_size;
    Point min, max;
} SpatialGrid;

/**
 * @brief Compute the axis-aligned bounding box (AABB) of a surface.
 *
 * @param polygons The surface structure containing mesh points.
 * @param min Output: the minimum bounding box corner (x_min, y_min, z_min).
 * @param max Output: the maximum bounding box corner (x_max, y_max, z_max).
 */
void get_polygon_bounding_box(polygons_struct *polygons, Point *min, Point *max)
{
    int i;
    double x, y, z;

    if (polygons->n_points == 0)
        return;

    x = Point_x(polygons->points[0]);
    y = Point_y(polygons->points[0]);
    z = Point_z(polygons->points[0]);

    fill_Point(*min, x, y, z);
    fill_Point(*max, x, y, z);

    for (i = 1; i < polygons->n_points; i++) {
        x = Point_x(polygons->points[i]);
        y = Point_y(polygons->points[i]);
        z = Point_z(polygons->points[i]);

        if (x < Point_x(*min)) Point_x(*min) = x;
        if (y < Point_y(*min)) Point_y(*min) = y;
        if (z < Point_z(*min)) Point_z(*min) = z;

        if (x > Point_x(*max)) Point_x(*max) = x;
        if (y > Point_y(*max)) Point_y(*max) = y;
        if (z > Point_z(*max)) Point_z(*max) = z;
    }
}

// Get grid index from 3D coordinate
int get_grid_index(int x, int y, int z, int res) {
    return x + y * res + z * res * res;
}


// Insert a vertex into a grid cell
void insert_into_grid(SpatialGrid *grid, int v, Point *points) {
    int xi = (int)((Point_x(points[v]) - Point_x(grid->min)) / grid->cell_size);
    int yi = (int)((Point_y(points[v]) - Point_y(grid->min)) / grid->cell_size);
    int zi = (int)((Point_z(points[v]) - Point_z(grid->min)) / grid->cell_size);

    if (xi < 0 || yi < 0 || zi < 0 || xi >= grid->res || yi >= grid->res || zi >= grid->res)
        return;

    int index = get_grid_index(xi, yi, zi, grid->res);

    PointNode *node = malloc(sizeof(PointNode));
    node->index = v;
    node->next = grid->cells[index].points;
    grid->cells[index].points = node;
}


// Build spatial grid from surface vertices
SpatialGrid *build_spatial_grid(polygons_struct *polygons, int res) {
    int i;
    SpatialGrid *grid = malloc(sizeof(SpatialGrid));
    grid->res = res;
    grid->cells = calloc(res * res * res, sizeof(GridCell));

    get_polygon_bounding_box(polygons, &grid->min, &grid->max);

    grid->cell_size = fmax(fmax(Point_x(grid->max) - Point_x(grid->min),
                                Point_y(grid->max) - Point_y(grid->min)),
                                Point_z(grid->max) - Point_z(grid->min)) / res;

    for (i = 0; i < polygons->n_points; i++) {
        insert_into_grid(grid, i, polygons->points);
    }

    return grid;
}


// Free memory used by spatial grid
void destroy_spatial_grid(SpatialGrid *grid) {
    int i, total = grid->res * grid->res * grid->res;
    for (i = 0; i < total; i++) {
        PointNode *node = grid->cells[i].points;
        while (node) {
            PointNode *tmp = node;
            node = node->next;
            free(tmp);
        }
    }
    free(grid->cells);
    free(grid);
}


// Estimate average edge length
double estimate_average_edge_length(polygons_struct *polygons, int *n_neighbours, int **neighbours) {
    double total = 0.0;
    int i, j, count = 0;
    for (i = 0; i < polygons->n_points; i++) {
        Point *p1 = &polygons->points[i];
        for (j = 0; j < n_neighbours[i]; j++) {
            int ni = neighbours[i][j];
            Point *p2 = &polygons->points[ni];
            total += distance_between_points(p1, p2);
            count++;
        }
    }
    return (count > 0) ? (total / count) : 1.0;
}


// Main function to find near-self-intersections
int *find_near_self_intersections(polygons_struct *polygons, double threshold_factor, int *n_hits_out) {
    int i, j, dx, dy, dz;
    int *n_neighbours, **neighbours;
    int *flags = calloc(polygons->n_points, sizeof(int));
    int n_hits = 0;

    check_polygons_neighbours_computed(polygons);
    create_polygon_point_neighbours(polygons, TRUE, &n_neighbours, &neighbours, NULL, NULL);

    SpatialGrid *grid = build_spatial_grid(polygons, GRID_RES);
    double threshold = estimate_average_edge_length(polygons, n_neighbours, neighbours) * threshold_factor;

    for (i = 0; i < polygons->n_points; i++) {
        Point p = polygons->points[i];
        int xi = (int)((Point_x(p) - Point_x(grid->min)) / grid->cell_size);
        int yi = (int)((Point_y(p) - Point_y(grid->min)) / grid->cell_size);
        int zi = (int)((Point_z(p) - Point_z(grid->min)) / grid->cell_size);

        double min_dist = DBL_MAX;
        int closest_index = -1;

        for (dx = -1; dx <= 1; dx++) {
            for (dy = -1; dy <= 1; dy++) {
                for (dz = -1; dz <= 1; dz++) {
                    int nx = xi + dx;
                    int ny = yi + dy;
                    int nz = zi + dz;

                    if (nx < 0 || ny < 0 || nz < 0 || nx >= grid->res || ny >= grid->res || nz >= grid->res)
                        continue;

                    int index = get_grid_index(nx, ny, nz, grid->res);
                    PointNode *node = grid->cells[index].points;

                    while (node) {
                        int ni = node->index;
                        if (ni != i) {
                            // Skip direct neighbors
                            int is_neighbor = 0;
                            for (j = 0; j < n_neighbours[i]; j++) {
                                if (neighbours[i][j] == ni) {
                                    is_neighbor = 1;
                                    break;
                                }
                            }
                            if (is_neighbor) {
                                node = node->next;
                                continue;
                            }

                            double d = distance_between_points(&p, &polygons->points[ni]);
                            if (d < min_dist) {
                                min_dist = d;
                                closest_index = ni;
                            }
                        }
                        node = node->next;
                    }
                }
            }
        }

        if (min_dist < threshold) {
            flags[i] = 1;
            n_hits++;
        }
    }

    destroy_spatial_grid(grid);
    delete_polygon_point_neighbours(polygons, n_neighbours, neighbours, NULL, NULL);

    if (n_hits_out)
        *n_hits_out = n_hits;

    return flags;  // Array of size n_points with 1=potential intersection, 0=none
}

int
intersect_poly_poly(int poly0, int poly1, polygons_struct *surface)
{
    int i, t[3], t2[3];

    for (i = 0; i < 3; i++) {
        t[i] = surface->indices[
                 POINT_INDEX(surface->end_indices, poly0, i)];
    }

    for (i = 0; i < 3; i++) {
        t2[i] = surface->indices[
               POINT_INDEX(surface->end_indices, poly1, i)];
    }

    return(intersect_triangle_triangle(t, t2, surface));
}


int
intersect_triangle_triangle(int pidx0[3], int pidx1[3],
              polygons_struct *surface)
{
    int i, result;
    Point pts[4];

    /* test if neighbors... if so, skip */
    for (i = 0; i < 3; i++) {
        if (pidx0[i] == pidx1[0]) return 0;
        if (pidx0[i] == pidx1[1]) return 0;
        if (pidx0[i] == pidx1[2]) return 0;
    }

    for (i = 0; i < 3; i++)
        pts[i] = surface->points[pidx0[i]];

    pts[3] = surface->points[pidx0[0]];
    
    for (i = 1; i < 4; i++) {
        result = intersect_segment_triangle(pts[i-1], pts[i], pidx1,
                          surface);
        if (result == 1) {
            return 1;
        }
    }

    return 0;
}


/*
 * intersect_segment_triangle(): test if a segment intersects a triangle
 *    Input:  two points for the segment p0 p1, triangle indices tpidx[3]
 *    Return: -1 = triangle is degenerate (a segment or point)
 *         0 = disjoint (no intersect)
 *         1 = intersect in unique point I1
 *         2 = are in the same plane
 */
int
intersect_segment_triangle(Point p0, Point p1, int tpidx[3],
               polygons_struct *surface)
{
    Vector u, v, n;       /* triangle vectors */
    Vector dir, w0, w;      /* ray vectors */
    Vector zero;
    float r, a, b;       /* params to calc ray-plane intersect */
    float uu, uv, vv, wu, wv, D;
    float s, t;
    Point pts[3], I;
    int i;

    for (i = 0; i < 3; i++)
        pts[i] = surface->points[tpidx[i]];

    /* get triangle edge vectors and plane normal */
    SUB_POINTS(u, pts[1], pts[0]);
    SUB_POINTS(v, pts[2], pts[0]);

    CROSS_VECTORS(n, u, v);

    fill_Vector(zero, 0.0, 0.0, 0.0);
    if (EQUAL_VECTORS(n, zero))
        return -1; /* triangle is degenerate */

    SUB_POINTS(dir, p1, p0); /* ray direction vector */
    SUB_POINTS(w0, p0, pts[0]);

    a = -DOT_VECTORS(n, w0);
    b = DOT_VECTORS(n, dir);

    if (fabs(b) < 1e-6) { /* ray is parallel to triangle plane */
        if (a == 0) { /* ray lies in triangle plane */
            return 2;
        } else return 0; /* ray disjoint from plane */
    }

    /* get intersect point of ray with triangle plane */
    r = a / b;
    if (r < 0.0 || r > 1.0)  /* no intersect */
        return 0;

    SCALE_VECTOR(dir, dir, r);
    ADD_POINT_VECTOR(I, p0, dir); /* intersect point of ray and plane */

    /* is I inside T? */
    uu = DOT_VECTORS(u, u);
    uv = DOT_VECTORS(u,v);
    vv = DOT_VECTORS(v,v);
    SUB_POINTS(w, I, pts[0]);
    wu = DOT_VECTORS(w,u);
    wv = DOT_VECTORS(w,v);
    D = uv * uv - uu * vv;

    /* get and test parametric coords */
    s = (uv * wv - vv * wu) / D;
    if (s < 0.0 || s > 1.0) /* I is outside T */
        return 0;
    t = (uv * wu - uu * wv) / D;
    if (t < 0.0 || (s + t) > 1.0)  /* I is outside T */
        return 0;

    return 1; /* I is in T */
}


int
find_selfintersections_masked(polygons_struct *polygons, int *defects, int *polydefects)
{
    int n_intersects, p;
    double *curvatures, threshold[2];
    int *n_neighbours, **neighbours;

    curvatures = (double *) malloc(sizeof(double) * polygons->n_points);
  
    get_all_polygon_point_neighbours(polygons, &n_neighbours, &neighbours);
    get_polygon_vertex_curvatures_cg(polygons, n_neighbours, neighbours,
                     0.0, 2, curvatures);

    double prctile[2] = {10, 90};
    get_prctile(curvatures, polygons->n_points, threshold, prctile, 0, DT_FLOAT64);
    /* Limit search for intersections to areas with large neg. curvature by 
       indicating with zeros */
    for (p = 0; p < polygons->n_items; p++)
        defects[p] = (curvatures[p] < threshold[1]) ? -1 : 0;

    update_polydefects(polygons, defects, polydefects);

    n_intersects = find_selfintersections(polygons, defects, polydefects, 0);
    
    free(curvatures);
    return n_intersects;
}

int
find_selfintersections(polygons_struct *polygons, int *defects, int *polydefects, int init)
{
    int n_intersects, size, i, n, p, p2, b;
    progress_struct progress;
    struct octree *tree;
    struct polynode *cur, *node;

    if (init)
        memset(polydefects, 0, sizeof(int) * polygons->n_items);

    tree = build_octree(polygons);

    initialize_progress_report(&progress, FALSE, polygons->n_items,
                   "find_selfintersections");

    n_intersects = 0;
    for (p = 0; p < polygons->n_items; p++) {
        node = tree->nodelist[p];
        
        /* skip check if neg. values in polydefects indicate that */
        if (!init && polydefects[p] < 0)
            continue;

        for (p2 = p+1; p2 < polygons->n_items; p2++)
            tree->polyflag[p2] = 0;
            /* skip check if neg. values in polydefects indicate that */
            if (!init && polydefects[p2] < 0)
                continue;

        for (b = 0; b < NBOXES; b++) {
            if (xintersect(node->bounds, tree->bounds[b]) == 0) {
                b += XINC - 1;
                continue;
            }
            if (yintersect(node->bounds, tree->bounds[b]) == 0) {
                b += YINC - 1;
                continue;
            }
            if (zintersect(node->bounds, tree->bounds[b]) == 0) 
                continue;

            for (cur = tree->nodes[b]; cur != NULL; cur = cur->next) {
            
                if (cur->num <= p) 
                    break;

                if (intersect(node->bounds, cur->bounds) == 0) 
                    continue;

                if (intersect_triangle_triangle(node->pts,
                                cur->pts,
                                polygons) == 0) 
                    continue;

                n_intersects++;
                polydefects[node->num] = n_intersects;
                polydefects[cur->num] = n_intersects;
            }
        }
        update_progress_report(&progress, p);
    }

    terminate_progress_report(&progress);

    update_defects(polygons, polydefects, defects);

    return n_intersects;
}


int
correct_simple_selfintersections(polygons_struct *surface, int *defects, int *polydefects, int *n_neighbours, int **neighbours)
{
    int i, n, p, n_intersections = 1;
    double extent = 0.05;

    for (p = 0; p < surface->n_points; p++)
        if (defects[p] > 0) n_intersections++;
    
    n = 0;
    while (n_intersections != 0) {
        i = 0;
        n++;
        for (p = 0; p < surface->n_points; p++) {
            if (defects[p] == 0) continue;
            i++;
            Point_x(surface->points[p]) += extent*Point_x(surface->normals[p]);
            Point_y(surface->points[p]) += extent*Point_y(surface->normals[p]);
            Point_z(surface->points[p]) += extent*Point_z(surface->normals[p]);
        }
        n_intersections = find_remaining_intersections(surface, defects, polydefects,
                      n_neighbours, neighbours);
        printf("%03d: %04d Remaining self intersections %d\n",n,i,n_intersections);
        if (n==40) n_intersections = 0;
    }

    return 0;
}

/* consolidate neighboring self-intersections */
int
join_intersections(polygons_struct *surface, int *defects, int *polydefects,
           int *n_neighbours, int **neighbours)
{
    int d, old_d;
    int n_intersects = 0, i, n, p, *dmap;

    update_defects(surface, polydefects, defects);

    for (i = 0; i < surface->n_points; i++) {
        if (defects[i] == 0)
            continue; /* skip */

        n_intersects = n_intersects > defects[i] ?
                 n_intersects : defects[i];

        for (n = 0; n < n_neighbours[i]; n++) {
            d = defects[neighbours[i][n]];
            if (d > 0 && d != defects[i]) {
                old_d = d > defects[i] ? d : defects[i];
                d = d < defects[i] ? d : defects[i];

                defects[i] = d;
                defects[neighbours[i][n]] = d;
                for (p = 0; p < surface->n_points; p++) {
                    if (defects[p] == old_d)
                        defects[p] = d;
                }
            }
        }
    }

    dmap = (int *) malloc(sizeof(int) * n_intersects);
    memset(dmap, 0, sizeof(int) * n_intersects);

    n_intersects = 0;
    for (i = 0; i < surface->n_points; i++) {
        if (defects[i] == 0)
            continue; /* skip */

        d = 0;
        for (n = 0; n < n_intersects; n++) {
            if (defects[i] == dmap[n]) {
                d = 1; /* already remapped */
                defects[i] = n + 1;
                break;
            }
        }
        if (d == 0) {
            dmap[n_intersects++] = defects[i];
            defects[i] = n_intersects;
        }
    }

    update_polydefects(surface, defects, polydefects);

    free(dmap);
    return n_intersects;
}


/* From a list of original self-intersections, find the ones that remain */
int
find_remaining_intersections(polygons_struct *surface, int *defects,
               int *polydefects, int *n_neighbours,
               int **neighbours)
{
    int p, p2, idx, i, siflag;

    update_defects(surface, polydefects, defects);

    for (p = 0; p < surface->n_items; p++) {
        if (polydefects[p] == 0)
            continue; /* skip */

        siflag = 0;

        for (p2 = 0; p2 < surface->n_items; p2++) {
            if (p2 == p) continue;

            if (intersect_poly_poly(p, p2, surface) == 1) {
                siflag = 1;
                break;
            }
        }
        if (siflag == 0) /* no intersections */
            polydefects[p] = 0;
    }

    memset(defects, 0, sizeof(int) * surface->n_points);
    for (p = 0; p < surface->n_items; p++) {
        if (polydefects[p] == 0) continue; /* skip */

        for (i = 0; i < 3; i++) {
            idx = surface->indices[
                 POINT_INDEX(surface->end_indices, p, i)];
            defects[idx] = polydefects[p];
        }
    }
    return(join_intersections(surface, defects, polydefects,
                  n_neighbours, neighbours));
}


/* patch "surface" with "patch" based on a defect list */
int
patch_selfintersections(polygons_struct *surface, polygons_struct *patch,
            int *defects, int *polydefects, int n_defects,
            int *n_neighbours, int **neighbours)
{
    int p;

    update_defects(surface, polydefects, defects);

    /* patch self-intersections */
    for (p = 0; p < surface->n_points; p++) {
        if (defects[p] != 0) { /* patch it */
            fill_Point(surface->points[p],
                   Point_x(patch->points[p]),
                   Point_y(patch->points[p]),
                   Point_z(patch->points[p]));
        }
    }

    /* consolidate remaining self-intersections */
    return(find_remaining_intersections(surface, defects, polydefects,
                      n_neighbours, neighbours));

}


/* smooth out self-intersections until they're repaired */
int
smooth_selfintersections(polygons_struct *surface, int *defects,
             int *polydefects, int n_defects,
             int *n_neighbours, int **neighbours, int maxiter)
{
    int p, i, iter, d, n, n2, npts;
    int *edgeflag, n_prev_defects, count;
    Point tp[3];
    double areas[128], centers[384], xyz[3], weight, t_area;
    FILE *fp;

    if (n_defects == 0)
         return(0); /* done! */

    update_defects(surface, polydefects, defects);

    edgeflag = (int *) malloc(sizeof(int) * surface->n_points);

    /* smooth defect areas.. increase defect area every 5th iter */
    iter = 0;
    n_prev_defects = n_defects;
    count = 0;
    
    while (n_defects != 0) {
        iter++;
        if (iter == maxiter)
            break; /* done */

        npts = 0; /* # of smoothed points */

        /* mark defect edges */
        memset(edgeflag, 0, sizeof(int) * surface->n_points);
        for (p = 0; p < surface->n_points; p++) {
            if (defects[p] == 0)
                continue; /* skip */

            for (n = 0; n < n_neighbours[p]; n++) {
                if (defects[neighbours[p][n]] != defects[p]) {
                    edgeflag[p] = defects[p]; /* an edge! */
                    break;
                }
            }
        }

        /* smooth out self-intersections */
        for (p = 0; p < surface->n_points; p++) {
            if (defects[p] == 0 || edgeflag[p] != 0)
                continue; /* skip */

            t_area = 0;
            tp[0] = surface->points[p];
            for (n = 0; n < n_neighbours[p]; n++) {
                n2 = (n + 1) % n_neighbours[p];

                /* area of the triangle */
                tp[1] = surface->points[neighbours[p][n]];
                tp[2] = surface->points[neighbours[p][n2]];
                areas[n] = get_polygon_surface_area(3, tp);

                t_area += areas[n];

                /* Save center of this tile */
                centers[n*3    ] = (Point_x(tp[0]) +
                          Point_x(tp[1]) +
                          Point_x(tp[2])) / 3.0;
                centers[n*3 + 1] = (Point_y(tp[0]) +
                          Point_y(tp[1]) +
                          Point_y(tp[2])) / 3.0;
                centers[n*3 + 2] = (Point_z(tp[0]) +
                          Point_z(tp[1]) +
                          Point_z(tp[2])) / 3.0;
            }
            if (t_area <= 0)
                continue; /* skip */

            /* Area Smoothing */
            xyz[0] = xyz[1] = xyz[2] = 0.0;
            for (n = 0; n <  n_neighbours[p]; n++) {
                weight = areas[n] / t_area;
                for (i = 0; i < 3; i++)
                    xyz[i] += weight * centers[n*3 + i];
            }
            fill_Point(surface->points[p], xyz[0], xyz[1], xyz[2]);
            npts++;
        }

        if (npts == 0) { /* nothing done, expand defects & restart */
            expand_defects(surface, defects, polydefects,
                     0, 1, n_neighbours, neighbours);
            iter--;
            continue;
        }

        /* test if self-intersections repaired */
        for (d = 1; d <= n_defects; d++) {
            if (has_selfintersections(surface, polydefects,
                           d) == 0) {
                /* delete it, it's fixed! */ 
                for (i = 0; i < surface->n_items; i++) {
                    if (polydefects[i] == d)
                        polydefects[i] = 0;
                    else if (polydefects[i] == n_defects)
                        polydefects[i] = d;
                }
                n_defects--;
                d--;
            }
        }
        update_defects(surface, polydefects, defects);

        if (n_defects == 0)
            break; /* all done! */

        if (npts > 100) {
            /* remap defects to limit # of affected points */
            n_defects = find_remaining_intersections(surface,
                                 defects,
                                 polydefects,
                                 n_neighbours,
                                 neighbours);
            /* increase counter if number of defects is unchanged */
            if (n_defects == n_prev_defects) count++; else count = 0;
                        
            /* expand remaining defects if no improvement can be found */
            if (count > 2) 
                expand_defects(surface, defects, polydefects,
                     0, 1, n_neighbours, neighbours);
            /* or stop if expanding does not help either */
            if (count > 3) break; 
        } else if (iter % 5 == 0) {
            /* expand remaining defects every 5th iteration */
            expand_defects(surface, defects, polydefects,
                     0, 1, n_neighbours, neighbours);
        }
        n_prev_defects = n_defects;
    }

    n_defects = find_remaining_intersections(surface, defects, polydefects,
                         n_neighbours, neighbours);

    free(edgeflag);
    return(n_defects);
}


/* determine if an area still has self-intersections (boolean function) */
int
has_selfintersections(polygons_struct *polygons, int *polydefects, int defect)
{
    int p, p2, i, t[3], t2[3];

    for (p = 0; p < polygons->n_items; p++) {
        if (polydefects[p] != defect)
            continue; /* skip */

        for (i = 0; i < 3; i++) {
            t[i] = polygons->indices[
                    POINT_INDEX(polygons->end_indices, p, i)];
        }

        for (p2 = 0; p2 < polygons->n_items; p2++) {
            if (p2 == p) continue;

            for (i = 0; i < 3; i++) {
                t2[i] = polygons->indices[
                   POINT_INDEX(polygons->end_indices, p2, i)];
            }
            if (intersect_triangle_triangle(t, t2, polygons))
                return 1;
        }
    }
    return 0;
}

/* Find and remove intersections */
void
remove_intersections(polygons_struct *polygons, int verbose)
{
    int *defects, *polydefects, n_intersects;
    int *n_neighbours, **neighbours;
    int counter;
    Point *new_pts;

    defects = (int *) malloc(sizeof(int) * polygons->n_points);
    polydefects = (int *) malloc(sizeof(int) * polygons->n_items);
    create_polygon_point_neighbours(polygons, TRUE, &n_neighbours,
                &neighbours, NULL, NULL);

    check_polygons_neighbours_computed(polygons);

    counter = 0;

    n_intersects = find_selfintersections(polygons, defects, polydefects, 1);        
    n_intersects = join_intersections(polygons, defects, polydefects,
              n_neighbours, neighbours);
    do {
        counter++;
        
        if (n_intersects > 0) {
            if (verbose)
                printf("%3d self intersections found that will be corrected.\n", n_intersects);

            n_intersects = smooth_selfintersections(polygons, defects, polydefects,
                   n_intersects, n_neighbours,
                   neighbours, 50);

        }
    } while (n_intersects > 0 && counter < 10);
    
    free(defects);
    free(polydefects);

    compute_polygon_normals(polygons);
        
}

/* Find and remove near self-intersections */
void
remove_near_intersections(polygons_struct *polygons, double threshold, int verbose)
{
    int *polydefects, n_intersects = 0;
    int *n_neighbours, **neighbours;
    int counter;
    Point *new_pts;

    polydefects = (int *) malloc(sizeof(int) * polygons->n_items);
    create_polygon_point_neighbours(polygons, TRUE, &n_neighbours,
                &neighbours, NULL, NULL);

    check_polygons_neighbours_computed(polygons);

    counter = 0;

    int *defects = find_near_self_intersections(polygons, threshold, &n_intersects);
    update_polydefects(polygons, defects, polydefects);

    n_intersects = join_intersections(polygons, defects, polydefects,
              n_neighbours, neighbours);
    do {
        counter++;
        
        if (n_intersects > 1) {
            if (verbose)
                printf("%3d self intersections found that will be corrected.\n", n_intersects);

            n_intersects = smooth_selfintersections(polygons, defects, polydefects,
                   n_intersects, n_neighbours,
                   neighbours, 200);

        }
    } while (n_intersects > 1 && counter < 5);
    
    free(polydefects);

    compute_polygon_normals(polygons);
        
}
