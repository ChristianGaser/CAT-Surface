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

#define GRID_RES 250 // Number of grid cells per axis (adjustable)

typedef struct PointNode
{
    int index;
    struct PointNode *next;
} PointNode;

typedef struct GridCell
{
    PointNode *points;
} GridCell;

typedef struct SpatialGrid
{
    GridCell *cells;
    int res;
    double cell_size;
    Point min, max;
} SpatialGrid;

/**
 * \brief Compute the axis-aligned bounding box (AABB) of a surface.
 *
 * \param polygons The surface structure containing mesh points.
 * \param min Output: the minimum bounding box corner (x_min, y_min, z_min).
 * \param max Output: the maximum bounding box corner (x_max, y_max, z_max).
 */
static void get_polygon_bounding_box(polygons_struct *polygons, Point *min, Point *max)
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

    for (i = 1; i < polygons->n_points; i++)
    {
        x = Point_x(polygons->points[i]);
        y = Point_y(polygons->points[i]);
        z = Point_z(polygons->points[i]);

        if (x < Point_x(*min))
            Point_x(*min) = x;
        if (y < Point_y(*min))
            Point_y(*min) = y;
        if (z < Point_z(*min))
            Point_z(*min) = z;

        if (x > Point_x(*max))
            Point_x(*max) = x;
        if (y > Point_y(*max))
            Point_y(*max) = y;
        if (z > Point_z(*max))
            Point_z(*max) = z;
    }
}

// Get grid index from 3D coordinate
static int get_grid_index(int x, int y, int z, int res)
{
    return x + y * res + z * res * res;
}

// Insert a vertex into a grid cell
static void insert_into_grid(SpatialGrid *grid, int v, Point *points)
{
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
static SpatialGrid *build_spatial_grid(polygons_struct *polygons, int res)
{
    int i;
    SpatialGrid *grid = malloc(sizeof(SpatialGrid));
    grid->res = res;
    grid->cells = calloc(res * res * res, sizeof(GridCell));

    get_polygon_bounding_box(polygons, &grid->min, &grid->max);

    grid->cell_size = fmax(fmax(Point_x(grid->max) - Point_x(grid->min),
                                Point_y(grid->max) - Point_y(grid->min)),
                           Point_z(grid->max) - Point_z(grid->min)) /
                      res;

    for (i = 0; i < polygons->n_points; i++)
    {
        insert_into_grid(grid, i, polygons->points);
    }

    return grid;
}

// Free memory used by spatial grid
static void destroy_spatial_grid(SpatialGrid *grid)
{
    int i, total = grid->res * grid->res * grid->res;
    for (i = 0; i < total; i++)
    {
        PointNode *node = grid->cells[i].points;
        while (node)
        {
            PointNode *tmp = node;
            node = node->next;
            free(tmp);
        }
    }
    free(grid->cells);
    free(grid);
}

// Estimate average edge length
static double estimate_average_edge_length(polygons_struct *polygons, int *n_neighbours, int **neighbours)
{
    double total = 0.0;
    int i, j, count = 0;
    for (i = 0; i < polygons->n_points; i++)
    {
        Point *p1 = &polygons->points[i];
        for (j = 0; j < n_neighbours[i]; j++)
        {
            int ni = neighbours[i][j];
            Point *p2 = &polygons->points[ni];
            total += distance_between_points(p1, p2);
            count++;
        }
    }
    return (count > 0) ? (total / count) : 1.0;
}

/* Shared core of find_near_self_intersections() and
 * find_near_facing_intersections(): flags vertices whose nearest non-neighbour
 * vertex (optionally restricted to opposing normals) is closer than
 * threshold_factor times the mean edge length. */
static int *
find_near_intersections(polygons_struct *polygons, double threshold_factor,
                        int facing_only, double min_opposition, int *n_hits_out)
{
    int i, j, dx, dy, dz;
    int *n_neighbours, **neighbours;
    int *flags = calloc(polygons->n_points, sizeof(int));
    int n_hits = 0;

    check_polygons_neighbours_computed(polygons);
    create_polygon_point_neighbours(polygons, TRUE, &n_neighbours, &neighbours, NULL, NULL);

    SpatialGrid *grid = build_spatial_grid(polygons, GRID_RES);
    double threshold = estimate_average_edge_length(polygons, n_neighbours, neighbours) * threshold_factor;

    for (i = 0; i < polygons->n_points; i++)
    {
        Point p = polygons->points[i];
        int xi = (int)((Point_x(p) - Point_x(grid->min)) / grid->cell_size);
        int yi = (int)((Point_y(p) - Point_y(grid->min)) / grid->cell_size);
        int zi = (int)((Point_z(p) - Point_z(grid->min)) / grid->cell_size);

        double min_dist = DBL_MAX;
        int closest_index = -1;

        for (dx = -1; dx <= 1; dx++)
        {
            for (dy = -1; dy <= 1; dy++)
            {
                for (dz = -1; dz <= 1; dz++)
                {
                    int nx = xi + dx;
                    int ny = yi + dy;
                    int nz = zi + dz;

                    if (nx < 0 || ny < 0 || nz < 0 || nx >= grid->res || ny >= grid->res || nz >= grid->res)
                        continue;

                    int index = get_grid_index(nx, ny, nz, grid->res);
                    PointNode *node = grid->cells[index].points;

                    while (node)
                    {
                        int ni = node->index;
                        if (ni != i)
                        {
                            // Skip direct neighbors
                            int is_neighbor = 0;
                            for (j = 0; j < n_neighbours[i]; j++)
                            {
                                if (neighbours[i][j] == ni)
                                {
                                    is_neighbor = 1;
                                    break;
                                }
                            }
                            if (is_neighbor)
                            {
                                node = node->next;
                                continue;
                            }

                            /* the other side of a sulcus or a blade faces the
                               opposite way; vertices of the same sheet do not */
                            if (facing_only &&
                                Point_x(polygons->normals[i]) * Point_x(polygons->normals[ni]) +
                                        Point_y(polygons->normals[i]) * Point_y(polygons->normals[ni]) +
                                        Point_z(polygons->normals[i]) * Point_z(polygons->normals[ni]) >
                                    -min_opposition)
                            {
                                node = node->next;
                                continue;
                            }

                            double d = distance_between_points(&p, &polygons->points[ni]);
                            if (d < min_dist)
                            {
                                min_dist = d;
                                closest_index = ni;
                            }
                        }
                        node = node->next;
                    }
                }
            }
        }

        if (min_dist < threshold)
        {
            flags[i] = 1;
            n_hits++;
        }
    }

    destroy_spatial_grid(grid);
    delete_polygon_point_neighbours(polygons, n_neighbours, neighbours, NULL, NULL);

    if (n_hits_out)
        *n_hits_out = n_hits;

    return flags; // Array of size n_points with 1=potential intersection, 0=none
}

/**
 * \brief Find vertices closer to a non-adjacent vertex than a distance threshold.
 *
 * Scans a spatial grid for the nearest vertex that is not a direct neighbour and
 * flags the vertex when it is closer than threshold_factor times the mean edge
 * length.  On an irregular (decimated) mesh most hits are 2-ring neighbours of the
 * same sheet rather than contacts; use find_near_facing_intersections() to count
 * only vertices of an opposing sheet.
 *
 * \param polygons         (in)  source 3D polygonal mesh (normals are not used)
 * \param threshold_factor (in)  multiplier for average edge length to define search radius
 * \param n_hits_out       (out) number of flagged vertices; may be NULL
 * \return Allocated array of flags (length = n_points, 1 = near hit), caller must free
 */
int *find_near_self_intersections(polygons_struct *polygons, double threshold_factor, int *n_hits_out)
{
    return find_near_intersections(polygons, threshold_factor, 0, 0.0, n_hits_out);
}

/**
 * \brief Find vertices close to a facing sheet of the same mesh.
 *
 * Like find_near_self_intersections(), but a nearby vertex only counts when its
 * normal opposes the normal of the query vertex (n_i . n_j < -min_opposition), as
 * across a sulcus or a thin blade.  Neighbouring vertices of the same sheet share
 * the normal direction and are ignored, which removes the false positives that a
 * pure distance test produces on meshes with irregular edge lengths: on reduced
 * central surfaces 88-99.9% of the vertices flagged by the distance test are 2-ring
 * neighbours of the same sheet.
 *
 * \param polygons         (in)  source 3D polygonal mesh with current normals
 * \param threshold_factor (in)  multiplier for average edge length to define search radius
 * \param min_opposition   (in)  required opposition of the normals (e.g. 0.3)
 * \param n_hits_out       (out) number of flagged vertices; may be NULL
 * \return Allocated array of flags (length = n_points, 1 = near hit), caller must free
 */
int *find_near_facing_intersections(polygons_struct *polygons, double threshold_factor,
                                    double min_opposition, int *n_hits_out)
{
    return find_near_intersections(polygons, threshold_factor, 1, min_opposition, n_hits_out);
}

/**
 * \brief Test whether one triangle intersects any other triangle of the mesh.
 *
 * Octree-accelerated replacement for scanning all n_items polygons: only triangles
 * whose bounding box shares an octree box with the query triangle are tested.
 * Because two triangles can only intersect geometrically when their bounding boxes
 * overlap, this returns exactly the same answer as the exhaustive scan.
 *
 * \param polygons (in) mesh to query
 * \param tree (in) octree built from the current vertex positions of polygons
 * \param p (in) index of the query triangle
 * \return 1 if triangle p intersects at least one other triangle, 0 otherwise
 */
static int poly_intersects_any(polygons_struct *polygons, struct octree *tree, int p)
{
    struct polynode *node = tree->nodelist[p];
    struct polynode *cur;
    int b;

    for (b = 0; b < NBOXES; b++)
    {
        if (xintersect(node->bounds, tree->bounds[b]) == 0)
        {
            b += XINC - 1;
            continue;
        }
        if (yintersect(node->bounds, tree->bounds[b]) == 0)
        {
            b += YINC - 1;
            continue;
        }
        if (zintersect(node->bounds, tree->bounds[b]) == 0)
            continue;

        for (cur = tree->nodes[b]; cur != NULL; cur = cur->next)
        {
            if (cur->num == p)
                continue;

            if (intersect(node->bounds, cur->bounds) == 0)
                continue;

            if (intersect_triangle_triangle(node->pts, cur->pts, polygons))
                return 1;
        }
    }

    return 0;
}

/**
 * \brief Test geometric intersection between two triangles in 3D space.
 *
 * Determines if two triangles defined by vertex indices intersect. Skips adjacent triangles
 * (sharing vertices) to avoid false positives. Uses segment-triangle intersection tests for
 * all three edges of the first triangle. Proper 3D geometric computation.
 *
 * \param pidx0    (in) int[3]; vertex indices of first triangle
 * \param pidx1    (in) int[3]; vertex indices of second triangle
 * \param surface (in) mesh containing vertex coordinates
 * \return 1 if triangles intersect, 0 if disjoint
 */
int intersect_triangle_triangle(int pidx0[3], int pidx1[3],
                                polygons_struct *surface)
{
    int i, result;
    Point pts[4];

    /* test if neighbors... if so, skip */
    for (i = 0; i < 3; i++)
    {
        if (pidx0[i] == pidx1[0])
            return 0;
        if (pidx0[i] == pidx1[1])
            return 0;
        if (pidx0[i] == pidx1[2])
            return 0;
    }

    for (i = 0; i < 3; i++)
        pts[i] = surface->points[pidx0[i]];

    pts[3] = surface->points[pidx0[0]];

    for (i = 1; i < 4; i++)
    {
        result = intersect_segment_triangle(pts[i - 1], pts[i], pidx1,
                                            surface);
        if (result == 1)
        {
            return 1;
        }
    }

    return 0;
}

/**
 * \brief Test if a line segment intersects a triangle in 3D space.
 *
 * Uses ray-triangle intersection algorithm: computes ray-plane intersection, tests
 * barycentric coordinates. Returns detailed status (-1=degenerate triangle, 0=disjoint,
 * 1=intersect at unique point, 2=coplanar). Based on M ller-Trumbore algorithm variant.
 *
 * \param p0 (in) first endpoint of line segment
 * \param p1 (in) second endpoint of line segment
 * \param tpidx    (in) int[3]; vertex indices of triangle
 * \param surface (in) mesh with vertex coordinates
 * \return -1 degenerate triangle, 0 no intersection, 1 unique intersection, 2 coplanar
 */
int intersect_segment_triangle(Point p0, Point p1, int tpidx[3],
                               polygons_struct *surface)
{
    Vector u, v, n;    /* triangle vectors */
    Vector dir, w0, w; /* ray vectors */
    Vector zero;
    float r, a, b; /* params to calc ray-plane intersect */
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

    if (fabs(b) < 1e-6)
    { /* ray is parallel to triangle plane */
        if (a == 0)
        { /* ray lies in triangle plane */
            return 2;
        }
        else
            return 0; /* ray disjoint from plane */
    }

    /* get intersect point of ray with triangle plane */
    r = a / b;
    if (r < 0.0 || r > 1.0) /* no intersect */
        return 0;

    SCALE_VECTOR(dir, dir, r);
    ADD_POINT_VECTOR(I, p0, dir); /* intersect point of ray and plane */

    /* is I inside T? */
    uu = DOT_VECTORS(u, u);
    uv = DOT_VECTORS(u, v);
    vv = DOT_VECTORS(v, v);
    SUB_POINTS(w, I, pts[0]);
    wu = DOT_VECTORS(w, u);
    wv = DOT_VECTORS(w, v);
    D = uv * uv - uu * vv;

    /* get and test parametric coords */
    s = (uv * wv - vv * wu) / D;
    if (s < 0.0 || s > 1.0) /* I is outside T */
        return 0;
    t = (uv * wu - uu * wv) / D;
    if (t < 0.0 || (s + t) > 1.0) /* I is outside T */
        return 0;

    return 1; /* I is in T */
}

/**
 * \brief Find self-intersecting triangles of a mesh.
 *
 * Tests every pair of triangles whose bounding boxes overlap in an octree and
 * labels both triangles of each intersecting pair with the running pair number.
 * The per-vertex labels are then derived from the per-triangle ones.
 *
 * \param polygons    (in)     triangle mesh
 * \param defects     (out)    per-vertex labels (0 = no intersection)
 * \param polydefects (in/out) per-triangle labels; with init == 0, triangles
 *                             with a negative label are skipped
 * \param init        (in)     non-zero to clear polydefects first
 * \return number of intersecting triangle pairs
 */
int find_selfintersections(polygons_struct *polygons, int *defects, int *polydefects, int init)
{
    int n_intersects, p, b;
    progress_struct progress;
    struct octree *tree;
    struct polynode *cur, *node;

    if (init)
        memset(polydefects, 0, sizeof(int) * polygons->n_items);

    tree = build_octree(polygons);

    initialize_progress_report(&progress, FALSE, polygons->n_items,
                               "find_selfintersections");

    n_intersects = 0;
    for (p = 0; p < polygons->n_items; p++)
    {
        node = tree->nodelist[p];

        /* skip check if neg. values in polydefects indicate that */
        if (!init && polydefects[p] < 0)
            continue;

        for (b = 0; b < NBOXES; b++)
        {
            if (xintersect(node->bounds, tree->bounds[b]) == 0)
            {
                b += XINC - 1;
                continue;
            }
            if (yintersect(node->bounds, tree->bounds[b]) == 0)
            {
                b += YINC - 1;
                continue;
            }
            if (zintersect(node->bounds, tree->bounds[b]) == 0)
                continue;

            for (cur = tree->nodes[b]; cur != NULL; cur = cur->next)
            {

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

    delete_octree(tree);

    update_defects(polygons, polydefects, defects);

    return n_intersects;
}

/**
 * \brief Find the representative of a label in a union-find forest (path halving).
 *
 * \param parent (in/out) union-find parent array
 * \param x (in) label to look up
 * \return representative label of the set containing x
 */
static int find_label(int *parent, int x)
{
    while (parent[x] != x)
    {
        parent[x] = parent[parent[x]]; /* path halving */
        x = parent[x];
    }
    return x;
}

/**
 * \brief Merge the sets of two labels, keeping the smaller label as representative.
 *
 * \param parent (in/out) union-find parent array
 * \param a (in) first label
 * \param b (in) second label
 */
static void union_labels(int *parent, int a, int b)
{
    a = find_label(parent, a);
    b = find_label(parent, b);

    if (a == b)
        return;

    if (a < b)
        parent[b] = a;
    else
        parent[a] = b;
}

/**
 * \brief Consolidate spatially-connected self-intersection regions into components.
 *
 * Groups neighboring defect vertices into shared defect IDs using local connectivity.
 * Merges defect regions that touch, reducing fragmented ID labels into contiguous
 * components. Remaps IDs sequentially (1..n_intersects) for organized processing.
 *
 * find_selfintersections() hands out a fresh label for every intersecting triangle
 * pair, so a single defect region typically arrives here split over hundreds of
 * labels.  Merging is therefore done with a union-find forest over the labels
 * (near-linear) instead of relabelling the whole vertex array per merge, which was
 * quadratic in the number of merges and dominated the runtime on surfaces with many
 * self-intersections.
 *
 * \param surface (in) mesh structure with neighborhood info
 * \param defects (in/out) per-vertex defect labels (remapped and consolidated)
 * \param polydefects (in/out) per-polygon defect labels (updated from consolidated vertex labels)
 * \param n_neighbours (in) neighbor counts per vertex
 * \param neighbours (in) neighbor lists per vertex
 * \return number of consolidated defect components
 */
int join_intersections(polygons_struct *surface, int *defects, int *polydefects,
                       int *n_neighbours, int **neighbours)
{
    int d, max_label = 0;
    int n_intersects = 0, i, n, *parent, *dmap;

    update_defects(surface, polydefects, defects);

    for (i = 0; i < surface->n_points; i++)
    {
        if (defects[i] > max_label)
            max_label = defects[i];
    }

    if (max_label == 0)
    { /* nothing marked */
        update_polydefects(surface, defects, polydefects);
        return 0;
    }

    /* labels that meet at neighbouring vertices belong to the same defect */
    parent = (int *)malloc(sizeof(int) * (max_label + 1));
    for (i = 0; i <= max_label; i++)
        parent[i] = i;

    for (i = 0; i < surface->n_points; i++)
    {
        if (defects[i] == 0)
            continue; /* skip */

        for (n = 0; n < n_neighbours[i]; n++)
        {
            d = defects[neighbours[i][n]];
            if (d > 0 && d != defects[i])
                union_labels(parent, defects[i], d);
        }
    }

    /* remap the surviving representatives to 1..n_intersects in order of
       first appearance */
    dmap = (int *)malloc(sizeof(int) * (max_label + 1));
    memset(dmap, 0, sizeof(int) * (max_label + 1));

    for (i = 0; i < surface->n_points; i++)
    {
        if (defects[i] == 0)
            continue; /* skip */

        d = find_label(parent, defects[i]);
        if (dmap[d] == 0)
            dmap[d] = ++n_intersects;
        defects[i] = dmap[d];
    }

    update_polydefects(surface, defects, polydefects);

    free(parent);
    free(dmap);
    return n_intersects;
}

/**
 * \brief Re-test marked self-intersections to identify which ones persist after correction.
 *
 * Scans all previously-marked intersecting polygons to verify if they still overlap.
 * Removes defect labels from successfully-corrected regions. Re-consolidates remaining
 * defects using join_intersections(). Used iteratively during correction process.
 *
 * \param surface (in) mesh after partial correction attempt
 * \param defects (in/out) per-vertex labels (recomputed from remaining intersections)
 * \param polydefects (in/out) per-polygon labels (cleared for resolved, updated for persistent)
 * \param n_neighbours (in) vertex connectivity
 * \param neighbours (in) neighbor lists
 * \return count of remaining unresolved intersection groups
 */
int find_remaining_intersections(polygons_struct *surface, int *defects,
                                 int *polydefects, int *n_neighbours,
                                 int **neighbours)
{
    int p, idx, i;
    struct octree *tree;

    update_defects(surface, polydefects, defects);

    tree = build_octree(surface);

    for (p = 0; p < surface->n_items; p++)
    {
        if (polydefects[p] == 0)
            continue; /* skip */

        if (poly_intersects_any(surface, tree, p) == 0)
            polydefects[p] = 0; /* no intersections */
    }

    delete_octree(tree);

    memset(defects, 0, sizeof(int) * surface->n_points);
    for (p = 0; p < surface->n_items; p++)
    {
        if (polydefects[p] == 0)
            continue; /* skip */

        for (i = 0; i < 3; i++)
        {
            idx = surface->indices[POINT_INDEX(surface->end_indices, p, i)];
            defects[idx] = polydefects[p];
        }
    }
    return (join_intersections(surface, defects, polydefects,
                               n_neighbours, neighbours));
}

/**
 * \brief Replace intersection regions with corrected coordinates from reference patch surface.
 *
 * Transfer vertices from reference patch surface to defective regions in target mesh.
 * Direct coordinate copy-based repair (vs. iterative smoothing). Fast but requires
 * pre-computed patch with same topology. Checks remaining intersections post-patch.
 *
 * \param surface (in/out) mesh to repair
 * \param patch (in) reference surface with corrected geometry
 * \param defects (in/out) per-vertex defect labels
 * \param polydefects (in/out) per-polygon defect labels
 * \param n_defects (in) number of distinct defects (for context)
 * \param n_neighbours (in) vertex connectivity
 * \param neighbours (in) neighbor lists
 * \return remaining self-intersections after patching
 */
int patch_selfintersections(polygons_struct *surface, polygons_struct *patch,
                            int *defects, int *polydefects, int n_defects,
                            int *n_neighbours, int **neighbours)
{
    int p;

    update_defects(surface, polydefects, defects);

    /* patch self-intersections */
    for (p = 0; p < surface->n_points; p++)
    {
        if (defects[p] != 0)
        { /* patch it */
            fill_Point(surface->points[p],
                       Point_x(patch->points[p]),
                       Point_y(patch->points[p]),
                       Point_z(patch->points[p]));
        }
    }

    /* consolidate remaining self-intersections */
    return (find_remaining_intersections(surface, defects, polydefects,
                                         n_neighbours, neighbours));
}

/**
 * \brief Iteratively smooth defect regions to resolve self-intersections via Laplacian relaxation.
 *
 * Repeatedly applies area-weighted Laplacian smoothing to defective vertices and neighbors.
 * Computes vertex weights based on adjacent triangle areas. Expands defect regions periodically
 * to catch nearby intersections. Continues until resolved or max iterations reached. Most
 * computationally-intensive but most flexible repair strategy.
 *
 * \param surface (in/out) mesh modified by iterative smoothing
 * \param defects (in/out) per-vertex defect labels
 * \param polydefects (in/out) per-polygon defect labels
 * \param n_defects (in) number of distinct defects to track
 * \param n_neighbours (in) vertex neighbor counts
 * \param neighbours (in/out) neighbor lists (may be updated)
 * \param maxiter (in) maximum smoothing iterations (typical 200-500)
 * \return number of defects successfully repaired
 */
int smooth_selfintersections(polygons_struct *surface, int *defects,
                             int *polydefects, int n_defects,
                             int *n_neighbours, int **neighbours, int maxiter)
{
    int p, i, iter, d, n, n2, npts, maxdefects;
    int *edgeflag, *siflags, n_prev_defects, count;
    Point tp[3];
    double areas[128], centers[384], xyz[3], weight, t_area;
    FILE *fp;

    if (n_defects == 0)
        return (0); /* done! */

    update_defects(surface, polydefects, defects);

    edgeflag = (int *)malloc(sizeof(int) * surface->n_points);

    /* find_remaining_intersections() can split a defect into several, so the
       defect count may grow during the loop; one label needs at least one
       polygon and labels are derived from vertex labels */
    maxdefects = surface->n_items > surface->n_points ? surface->n_items : surface->n_points;
    siflags = (int *)malloc(sizeof(int) * (maxdefects + 1));

    /* smooth defect areas.. increase defect area every 5th iter */
    iter = 0;
    n_prev_defects = n_defects;
    count = 0;

    while (n_defects != 0)
    {
        iter++;
        if (iter == maxiter)
            break; /* done */

        npts = 0; /* # of smoothed points */

        /* mark defect edges */
        memset(edgeflag, 0, sizeof(int) * surface->n_points);
        for (p = 0; p < surface->n_points; p++)
        {
            if (defects[p] == 0)
                continue; /* skip */

            for (n = 0; n < n_neighbours[p]; n++)
            {
                if (defects[neighbours[p][n]] != defects[p])
                {
                    edgeflag[p] = defects[p]; /* an edge! */
                    break;
                }
            }
        }

        /* smooth out self-intersections */
        for (p = 0; p < surface->n_points; p++)
        {
            if (defects[p] == 0 || edgeflag[p] != 0)
                continue; /* skip */

            t_area = 0;
            tp[0] = surface->points[p];
            for (n = 0; n < n_neighbours[p]; n++)
            {
                n2 = (n + 1) % n_neighbours[p];

                /* area of the triangle */
                tp[1] = surface->points[neighbours[p][n]];
                tp[2] = surface->points[neighbours[p][n2]];
                areas[n] = get_polygon_surface_area(3, tp);

                t_area += areas[n];

                /* Save center of this tile */
                centers[n * 3] = (Point_x(tp[0]) +
                                  Point_x(tp[1]) +
                                  Point_x(tp[2])) /
                                 3.0;
                centers[n * 3 + 1] = (Point_y(tp[0]) +
                                      Point_y(tp[1]) +
                                      Point_y(tp[2])) /
                                     3.0;
                centers[n * 3 + 2] = (Point_z(tp[0]) +
                                      Point_z(tp[1]) +
                                      Point_z(tp[2])) /
                                     3.0;
            }
            if (t_area <= 0)
                continue; /* skip */

            /* Area Smoothing */
            xyz[0] = xyz[1] = xyz[2] = 0.0;
            for (n = 0; n < n_neighbours[p]; n++)
            {
                weight = areas[n] / t_area;
                for (i = 0; i < 3; i++)
                    xyz[i] += weight * centers[n * 3 + i];
            }
            fill_Point(surface->points[p], xyz[0], xyz[1], xyz[2]);
            npts++;
        }

        if (npts == 0)
        { /* nothing done, expand defects & restart */
            expand_defects(surface, defects, polydefects,
                           0, 1, n_neighbours, neighbours);
            iter--;
            continue;
        }

        /* test if self-intersections repaired - one octree pass for all
           defects instead of one octree per defect */
        find_intersecting_defects(surface, polydefects, n_defects, siflags);

        for (d = 1; d <= n_defects; d++)
        {
            if (siflags[d] == 0)
            {
                /* delete it, it's fixed! */
                for (i = 0; i < surface->n_items; i++)
                {
                    if (polydefects[i] == d)
                        polydefects[i] = 0;
                    else if (polydefects[i] == n_defects)
                        polydefects[i] = d;
                }
                /* the last defect took the slot of the deleted one */
                siflags[d] = siflags[n_defects];
                n_defects--;
                d--;
            }
        }
        update_defects(surface, polydefects, defects);

        if (n_defects == 0)
            break; /* all done! */

        if (npts > 100)
        {
            /* remap defects to limit # of affected points */
            n_defects = find_remaining_intersections(surface,
                                                     defects,
                                                     polydefects,
                                                     n_neighbours,
                                                     neighbours);
            /* increase counter if number of defects is unchanged */
            if (n_defects == n_prev_defects)
                count++;
            else
                count = 0;

            /* expand remaining defects if no improvement can be found */
            if (count > 2)
                expand_defects(surface, defects, polydefects,
                               0, 1, n_neighbours, neighbours);
            /* or stop if expanding does not help either */
            if (count > 3)
                break;
        }
        else if (iter % 5 == 0)
        {
            /* expand remaining defects every 5th iteration */
            expand_defects(surface, defects, polydefects,
                           0, 1, n_neighbours, neighbours);
        }
        n_prev_defects = n_defects;
    }

    n_defects = find_remaining_intersections(surface, defects, polydefects,
                                             n_neighbours, neighbours);

    free(edgeflag);
    free(siflags);
    return (n_defects);
}

/**
 * \brief Test all defect regions for remaining self-intersections in a single pass.
 *
 * Builds the octree once and marks every defect that still contains an
 * intersecting triangle. Testing each defect on its own would rebuild the octree
 * n_defects times, which dominates the runtime of the repair loop as soon as a
 * surface has more than a handful of defects.
 *
 * \param polygons (in) mesh to check
 * \param polydefects (in) per-polygon defect labels
 * \param n_defects (in) largest defect ID in use
 * \param siflags (out) array of n_defects+1 entries; siflags[d] is set to 1 if defect
 *                      d still self-intersects, 0 otherwise (index 0 is unused)
 * \return number of defects that still self-intersect
 */
int find_intersecting_defects(polygons_struct *polygons, int *polydefects,
                              int n_defects, int *siflags)
{
    int p, d, n_remaining = 0;
    struct octree *tree;

    memset(siflags, 0, sizeof(int) * (n_defects + 1));

    tree = build_octree(polygons);

    for (p = 0; p < polygons->n_items; p++)
    {
        d = polydefects[p];
        if (d <= 0 || d > n_defects)
            continue; /* skip */

        if (siflags[d]) /* this defect is already known to intersect */
            continue;

        if (poly_intersects_any(polygons, tree, p))
        {
            siflags[d] = 1;
            n_remaining++;
        }
    }

    delete_octree(tree);
    return n_remaining;
}

/**
 * \brief Remove all self-intersections from a mesh with explicit iteration limits.
 *
 * The two loop limits are exposed because they are what governs the runtime:
 * each pass re-detects the remaining defects and runs up to maxiter smoothing
 * iterations on them.
 *
 * \param polygons (in/out) mesh to repair
 * \param max_passes (in) maximum number of detect/smooth passes (default 10)
 * \param maxiter (in) maximum smoothing iterations per pass (default 50)
 * \param verbose (in) 1 for progress output; 0 for silent
 * \return number of self-intersecting defect regions that remain (0 = fully repaired)
 */
int remove_intersections_iter(polygons_struct *polygons, int max_passes,
                              int maxiter, int verbose)
{
    int *defects, *polydefects, n_intersects = 0;
    int *n_neighbours, **neighbours;
    int counter;

    defects = (int *)malloc(sizeof(int) * polygons->n_points);
    polydefects = (int *)malloc(sizeof(int) * polygons->n_items);
    create_polygon_point_neighbours(polygons, TRUE, &n_neighbours,
                                    &neighbours, NULL, NULL);

    check_polygons_neighbours_computed(polygons);

    /* Detect afresh before every pass.  smooth_selfintersections() only
     * re-tests the triangles it already labelled, but smoothing a defect can
     * make it intersect triangles outside that label.  Carrying the labels over
     * from one pass to the next therefore kept smoothing stale regions while
     * the actual intersections were never seen, and a few defects survived all
     * passes. */
    for (counter = 0; counter < max_passes; counter++)
    {
        find_selfintersections(polygons, defects, polydefects, 1);
        n_intersects = join_intersections(polygons, defects, polydefects,
                                          n_neighbours, neighbours);
        if (n_intersects == 0)
            break;

        if (verbose)
            printf("%3d self intersections found that will be corrected.\n", n_intersects);

        smooth_selfintersections(polygons, defects, polydefects,
                                 n_intersects, n_neighbours,
                                 neighbours, maxiter);
    }

    if (counter >= max_passes)
    {
        /* passes exhausted: report what the last pass left */
        find_selfintersections(polygons, defects, polydefects, 1);
        n_intersects = join_intersections(polygons, defects, polydefects,
                                          n_neighbours, neighbours);
    }

    free(defects);
    free(polydefects);
    delete_polygon_point_neighbours(polygons, n_neighbours, neighbours,
                                    NULL, NULL);

    compute_polygon_normals(polygons);

    return n_intersects;
}

/**
 * \brief Remove self-intersections, retreating stubborn defects towards a reference.
 *
 * Runs remove_intersections_iter() and, where defects survive it, moves their
 * vertices together with CAT_RETREAT_RINGS rings of neighbours by
 * CAT_RETREAT_FRACTION of the way back to the reference positions, then
 * repairs again, for at most CAT_RETREAT_STEPS steps, stopping early when two
 * steps in a row bring no progress -- which is what happens where the
 * reference crosses itself in the same place.
 *
 * The smoothing inside remove_intersections_iter() can only undo crossings it
 * can reach by relaxing the mesh locally.  Two sheets that pass through each
 * other stay crossed however long they are smoothed: the two sides of a thin
 * gyral blade that a deformation drove through each other (the white surface
 * of thin gyri: the target isovalue is never reached inside the blade), or a
 * sheet of a deformed central surface folded onto itself.  Measured on T1Prep
 * surfaces, three times the passes or four times the iterations left 265 of
 * 269 pairs.  Any surface a deformation started from gives the way back: the
 * central surface for the pial and white surfaces, the start mesh for the
 * central one.  One step resolved 86-226 remaining pairs by moving 0.2-0.4%
 * of the vertices, with the mean label error of the surface unchanged to
 * 5e-4.  The reference need not be free of intersections everywhere, only
 * where the surface is repaired.
 *
 * \param polygons   (in/out) mesh to repair
 * \param reference  (in)     reference positions, one per vertex of polygons
 *                            (same topology, e.g. the start of the
 *                            deformation); NULL makes this identical to
 *                            remove_intersections_iter()
 * \param max_passes (in)     detect/smooth passes of each repair (default 10)
 * \param maxiter    (in)     smoothing iterations per pass (default 50)
 * \param verbose    (in)     1 for progress output; 0 for silent
 * \return number of self-intersecting defect regions that remain (0 = fully repaired)
 */
int remove_intersections_ref(polygons_struct *polygons, const Point *reference,
                             int max_passes, int maxiter, int verbose)
{
    int *defects, *polydefects, *n_neighbours, **neighbours;
    int p, step, n_remaining, n_previous, n_stalled = 0;

    n_remaining = remove_intersections_iter(polygons, max_passes, maxiter, verbose);
    if (n_remaining == 0 || reference == NULL)
        return n_remaining;

    defects = (int *)malloc(sizeof(int) * polygons->n_points);
    polydefects = (int *)malloc(sizeof(int) * polygons->n_items);
    if (!defects || !polydefects)
    {
        free(defects);
        free(polydefects);
        return n_remaining;
    }
    create_polygon_point_neighbours(polygons, TRUE, &n_neighbours,
                                    &neighbours, NULL, NULL);

    for (step = 0; step < CAT_RETREAT_STEPS && n_remaining > 0; step++)
    {
        n_previous = n_remaining;

        if (verbose)
            printf("%3d self intersections left, retreating towards the reference.\n",
                   n_remaining);

        /* the vertices of the remaining defects and their neighbourhood */
        find_selfintersections(polygons, defects, polydefects, 1);
        expand_defects(polygons, defects, polydefects, 0, CAT_RETREAT_RINGS,
                       n_neighbours, neighbours);

        for (p = 0; p < polygons->n_points; p++)
        {
            if (defects[p] == 0)
                continue;
            Point_x(polygons->points[p]) += CAT_RETREAT_FRACTION *
                (Point_x(reference[p]) - Point_x(polygons->points[p]));
            Point_y(polygons->points[p]) += CAT_RETREAT_FRACTION *
                (Point_y(reference[p]) - Point_y(polygons->points[p]));
            Point_z(polygons->points[p]) += CAT_RETREAT_FRACTION *
                (Point_z(reference[p]) - Point_z(polygons->points[p]));
        }

        n_remaining = remove_intersections_iter(polygons, max_passes, maxiter,
                                                verbose);

        /* Retreating cannot help where the reference crosses itself too, and
           every step costs a full repair -- give up after two without
           progress. */
        n_stalled = (n_remaining < n_previous) ? 0 : n_stalled + 1;
        if (n_stalled > 1)
        {
            if (verbose)
                printf("%3d self intersections left that the reference cannot "
                       "resolve.\n", n_remaining);
            break;
        }
    }

    free(defects);
    free(polydefects);
    delete_polygon_point_neighbours(polygons, n_neighbours, neighbours,
                                    NULL, NULL);

    return n_remaining;
}

/* Find and remove near self-intersections */
/**
 * \brief Remove near-intersecting vertices by iterative vertex repositioning.
 *
 * Fixes near-intersection problems by finding pairs of non-adjacent vertices
 * within a distance threshold and repositioning them toward the surface
 * centroid. Applies iterative correction until no intersections remain or
 * maximum iterations exceeded. Unlike remove_intersections_iter(), which handles
 * topological self-intersections, this function targets geometric near-collisions
 * that may not cause topological defects but indicate surface quality issues.
 *
 * \param polygons   (in)  source 3D polygonal mesh to be modified in-place
 * \param threshold  (in)  distance threshold for near-intersection detection (typically 0.05-0.20 times edge length)
 * \param verbose    (in)  1 to print progress messages to stdout, 0 for silent operation
 */
void remove_near_intersections(polygons_struct *polygons, double threshold, int verbose)
{
    int *polydefects, n_intersects = 0;
    int *n_neighbours, **neighbours;
    int counter;
    Point *new_pts;

    polydefects = (int *)malloc(sizeof(int) * polygons->n_items);
    create_polygon_point_neighbours(polygons, TRUE, &n_neighbours,
                                    &neighbours, NULL, NULL);

    check_polygons_neighbours_computed(polygons);

    counter = 0;

    int *defects = find_near_self_intersections(polygons, threshold, &n_intersects);
    update_polydefects(polygons, defects, polydefects);

    n_intersects = join_intersections(polygons, defects, polydefects,
                                      n_neighbours, neighbours);
    do
    {
        counter++;

        if (n_intersects > 1)
        {
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
