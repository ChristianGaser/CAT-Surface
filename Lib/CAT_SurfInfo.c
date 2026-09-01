/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 *
 */

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>

#include "CAT_SurfInfo.h"
#include "CAT_Surf.h"
#include "CAT_SafeAlloc.h"

/* A polygon of this area or smaller is treated as degenerate. */
#define DEGENERATE_AREA 1e-12

/** Undirected edge key, a < b, used to enumerate the edges of the mesh. */
typedef struct {
    int a;
    int b;
} surf_edge;

/**
 * \brief Order undirected edges lexicographically so equal edges become adjacent.
 *
 * \param p1 (in) first edge
 * \param p2 (in) second edge
 * \return negative, zero or positive as p1 sorts before, with or after p2
 */
static int
compare_edges(const void *p1, const void *p2)
{
    const surf_edge *e1 = (const surf_edge *) p1;
    const surf_edge *e2 = (const surf_edge *) p2;

    if (e1->a != e2->a)
        return (e1->a < e2->a) ? -1 : 1;
    if (e1->b != e2->b)
        return (e1->b < e2->b) ? -1 : 1;
    return 0;
}

/**
 * \brief Find the representative of a vertex in a union-find forest.
 *
 * \param parent (in/out) union-find parent array
 * \param x      (in)     vertex to look up
 * \return representative of the set containing x
 */
static int
find_root(int *parent, int x)
{
    while (parent[x] != x) {
        parent[x] = parent[parent[x]]; /* path halving */
        x = parent[x];
    }
    return x;
}

/**
 * \brief Merge the sets of two vertices in a union-find forest.
 *
 * \param parent (in/out) union-find parent array
 * \param a      (in)     first vertex
 * \param b      (in)     second vertex
 * \return void
 */
static void
union_roots(int *parent, int a, int b)
{
    a = find_root(parent, a);
    b = find_root(parent, b);

    if (a == b)
        return;

    if (a < b)
        parent[b] = a;
    else
        parent[a] = b;
}

/**
 * \brief Accumulate area, area statistics, enclosed volume and centroid.
 *
 * The volume is the signed integral obtained from the divergence theorem, so
 * it is only meaningful for a closed mesh; its sign reflects the orientation
 * of the polygon normals.  Non-triangular polygons are fan-triangulated.
 *
 * \param polygons (in)  input mesh
 * \param info     (out) area, volume, centroid and area statistics
 * \return void
 */
static void
accumulate_geometry(polygons_struct *polygons, surf_info_struct *info)
{
    int poly, size, k;
    double area, volume = 0.0, area_sum = 0.0;
    double cx = 0.0, cy = 0.0, cz = 0.0;
    Point points[MAX_POINTS_PER_POLYGON];

    info->area_min = 0.0;
    info->area_max = 0.0;
    info->area_mean = 0.0;

    for (poly = 0; poly < polygons->n_items; poly++) {
        size = get_polygon_points(polygons, poly, points);
        if (size < 3)
            continue;

        if (size == 3)
            info->n_triangles++;
        else
            info->n_nontriangular++;

        area = get_polygon_surface_area(size, points);
        if (isnan(area))
            area = 0.0;

        if (area <= DEGENERATE_AREA)
            info->n_degenerate_polygons++;

        if (poly == 0 || area < info->area_min)
            info->area_min = area;
        if (poly == 0 || area > info->area_max)
            info->area_max = area;
        area_sum += area;

        /* fan-triangulate from the first vertex */
        for (k = 1; k < size - 1; k++) {
            double x0 = Point_x(points[0]),   y0 = Point_y(points[0]),   z0 = Point_z(points[0]);
            double x1 = Point_x(points[k]),   y1 = Point_y(points[k]),   z1 = Point_z(points[k]);
            double x2 = Point_x(points[k+1]), y2 = Point_y(points[k+1]), z2 = Point_z(points[k+1]);

            volume += (x0 * (y1 * z2 - z1 * y2)
                     - y0 * (x1 * z2 - z1 * x2)
                     + z0 * (x1 * y2 - y1 * x2)) / 6.0;
        }

        /* area-weighted centroid of the polygon centroids */
        for (k = 0; k < size; k++) {
            cx += area * Point_x(points[k]) / (double) size;
            cy += area * Point_y(points[k]) / (double) size;
            cz += area * Point_z(points[k]) / (double) size;
        }
    }

    info->surface_area = area_sum;
    info->volume = volume;

    if (polygons->n_items > 0)
        info->area_mean = area_sum / (double) polygons->n_items;

    if (area_sum > 0.0) {
        info->centroid[0] = cx / area_sum;
        info->centroid[1] = cy / area_sum;
        info->centroid[2] = cz / area_sum;
    }
}

/**
 * \brief Enumerate the undirected edges and derive edge topology and lengths.
 *
 * Every polygon contributes its border as half-edges; sorting them groups the
 * duplicates, so the multiplicity of an edge is the number of polygons that
 * use it: one means a boundary edge, two a manifold edge, more a non-manifold
 * one.  Connected components are counted with a union-find forest over the
 * same edges.
 *
 * \param polygons (in)  input mesh
 * \param info     (out) edge counts, edge lengths and component count
 * \return void
 */
static void
accumulate_edges(polygons_struct *polygons, surf_info_struct *info)
{
    int poly, k, size, start, i, j, n_half, n_edges = 0;
    int *parent, *used;
    double len, len_sum = 0.0;
    surf_edge *edges;

    n_half = (polygons->n_items > 0) ?
             polygons->end_indices[polygons->n_items - 1] : 0;
    if (n_half == 0 || polygons->n_points == 0)
        return;

    edges = SAFE_MALLOC(surf_edge, n_half);
    parent = SAFE_MALLOC(int, polygons->n_points);
    used = SAFE_CALLOC(int, polygons->n_points);

    for (i = 0; i < polygons->n_points; i++)
        parent[i] = i;

    n_half = 0;
    for (poly = 0; poly < polygons->n_items; poly++) {
        start = (poly == 0) ? 0 : polygons->end_indices[poly - 1];
        size = polygons->end_indices[poly] - start;

        for (k = 0; k < size; k++) {
            int a = polygons->indices[start + k];
            int b = polygons->indices[start + (k + 1) % size];

            used[a] = 1;
            if (a == b)
                continue; /* degenerate half-edge */

            edges[n_half].a = (a < b) ? a : b;
            edges[n_half].b = (a < b) ? b : a;
            n_half++;

            union_roots(parent, a, b);
        }
    }

    qsort(edges, (size_t) n_half, sizeof(surf_edge), compare_edges);

    for (i = 0; i < n_half; i = j) {
        for (j = i + 1; j < n_half && compare_edges(&edges[i], &edges[j]) == 0; j++)
            ;

        n_edges++;
        if (j - i == 1)
            info->n_boundary_edges++;
        else if (j - i > 2)
            info->n_nonmanifold_edges++;

        len = distance_between_points(&polygons->points[edges[i].a],
                                      &polygons->points[edges[i].b]);
        if (n_edges == 1 || len < info->edge_min)
            info->edge_min = len;
        if (n_edges == 1 || len > info->edge_max)
            info->edge_max = len;
        len_sum += len;
    }

    info->n_edges = n_edges;
    if (n_edges > 0)
        info->edge_mean = len_sum / (double) n_edges;

    /* count components over the vertices that are actually used */
    for (i = 0; i < polygons->n_points; i++) {
        if (!used[i]) {
            info->n_unreferenced_points++;
            continue;
        }
        if (find_root(parent, i) == i)
            info->n_components++;
    }

    free(edges);
    free(parent);
    free(used);
}

/**
 * \brief Count self-intersections of the mesh.
 *
 * Runs the same two-stage test as CAT_SurfSelfIntersect: intersecting triangle
 * pairs are found first, then spatially connected ones are merged into one
 * region each.  The hit count inherits the octree double-counting described in
 * CAT_SurfInfo.h; the region and polygon counts do not.
 *
 * \param polygons (in)  input mesh
 * \param info     (out) intersection counts
 * \return void
 */
static void
accumulate_intersections(polygons_struct *polygons, surf_info_struct *info)
{
    int *defects, *polydefects, *n_neighbours, **neighbours;
    int i, n_marked = 0;

    defects = SAFE_MALLOC(int, polygons->n_points);
    polydefects = SAFE_MALLOC(int, polygons->n_items);

    create_polygon_point_neighbours(polygons, TRUE, &n_neighbours,
                                    &neighbours, NULL, NULL);

    info->n_triangle_hits = find_selfintersections(polygons, defects,
                                                   polydefects, 1);
    info->n_self_intersections = join_intersections(polygons, defects,
                                                    polydefects, n_neighbours,
                                                    neighbours);

    for (i = 0; i < polygons->n_items; i++) {
        if (polydefects[i] > 0)
            n_marked++;
    }
    info->n_intersecting_polygons = n_marked;

    delete_polygon_point_neighbours(polygons, n_neighbours, neighbours,
                                    NULL, NULL);
    free(defects);
    free(polydefects);
}

/**
 * \brief Collect geometric and topological information about a mesh.
 *
 * Computes counts, edge topology, Euler characteristic, connected components,
 * surface area, enclosed volume, bounding box and (optionally) the number of
 * self-intersections.  The mesh is not modified.
 *
 * \param polygons            (in)  input mesh
 * \param info                (out) filled summary, see surf_info_struct
 * \param check_intersections (in)  if non-zero, run the self-intersection test
 * \param verbose             (in)  if non-zero, print progress information
 * \return void (no return value).
 */
void
compute_surf_info(polygons_struct *polygons, surf_info_struct *info,
                  int check_intersections, int verbose)
{
    memset(info, 0, sizeof(*info));

    info->n_self_intersections = CAT_SURFINFO_SKIPPED;
    info->n_triangle_hits = CAT_SURFINFO_SKIPPED;
    info->n_intersecting_polygons = CAT_SURFINFO_SKIPPED;

    info->n_points = polygons->n_points;
    info->n_polygons = polygons->n_items;

    if (verbose)
        fprintf(stderr, "[Info] area, volume and bounding box\n");

    accumulate_geometry(polygons, info);
    get_bounds(polygons, info->bounds);

    if (verbose)
        fprintf(stderr, "[Info] edges and connected components\n");

    accumulate_edges(polygons, info);

    if (verbose)
        fprintf(stderr, "[Info] Euler characteristic\n");

    info->euler = euler_characteristic(polygons, verbose);

    info->closed = (info->n_boundary_edges == 0 &&
                    info->n_nonmanifold_edges == 0 &&
                    info->n_polygons > 0);

    /* chi = 2*components - 2*genus holds for closed orientable surfaces */
    if (info->closed)
        info->genus = (2 * info->n_components - info->euler) / 2;

    if (check_intersections) {
        if (verbose)
            fprintf(stderr, "[Info] self-intersections\n");
        accumulate_intersections(polygons, info);
    }
}

/**
 * \brief Print a human-readable report of a surf_info_struct.
 *
 * \param info (in) summary filled by compute_surf_info()
 * \param name (in) name shown in the header, may be NULL
 * \param fp   (in) output stream
 * \return void (no return value).
 */
void
print_surf_info(const surf_info_struct *info, const char *name, FILE *fp)
{
    static const char *rule =
        "---------------------------------------------------------------\n";
    int euler_from_edges;

    fprintf(fp, "\n%s", rule);
    fprintf(fp, " Surface information%s%s\n", name ? ": " : "", name ? name : "");
    fprintf(fp, "%s", rule);

    fprintf(fp, " Vertices                     %d\n", info->n_points);
    fprintf(fp, " Faces (polygons)             %d\n", info->n_polygons);
    fprintf(fp, " Triangular faces             %d\n", info->n_triangles);
    if (info->n_nontriangular > 0)
        fprintf(fp, " Non-triangular faces         %d\n", info->n_nontriangular);
    fprintf(fp, " Edges                        %d\n", info->n_edges);
    fprintf(fp, " Euler number                 %d\n", info->euler);
    if (info->closed)
        fprintf(fp, " Genus                        %d\n", info->genus);
    else
        fprintf(fp, " Genus                        n/a (surface is not closed)\n");
    fprintf(fp, " Connected components         %d\n", info->n_components);
    fprintf(fp, " Closed surface               %s\n", info->closed ? "yes" : "no");
    fprintf(fp, " Boundary edges               %d\n", info->n_boundary_edges);
    fprintf(fp, " Non-manifold edges           %d\n", info->n_nonmanifold_edges);
    fprintf(fp, " Degenerate polygons          %d\n", info->n_degenerate_polygons);
    if (info->n_unreferenced_points > 0)
        fprintf(fp, " Unreferenced vertices        %d\n",
                info->n_unreferenced_points);

    /* V - E + F must agree with the neighbour-based Euler characteristic */
    euler_from_edges = info->n_points - info->n_edges + info->n_polygons;
    if (euler_from_edges != info->euler)
        fprintf(fp, " WARNING: V-E+F = %d disagrees with the Euler number, the\n"
                    "          mesh connectivity is not manifold.\n",
                euler_from_edges);

    fprintf(fp, "%s", rule);
    fprintf(fp, " Surface area                 %.2f mm^2\n", info->surface_area);
    if (info->closed)
        fprintf(fp, " Enclosed volume              %.2f mm^3%s\n",
                fabs(info->volume),
                (info->volume < 0.0) ? "  (inward-facing normals)" : "");
    else
        fprintf(fp, " Enclosed volume              n/a (surface is not closed)\n");
    fprintf(fp, " Polygon area mean/min/max    %.4f / %.4f / %.4f mm^2\n",
            info->area_mean, info->area_min, info->area_max);
    fprintf(fp, " Edge length mean/min/max     %.4f / %.4f / %.4f mm\n",
            info->edge_mean, info->edge_min, info->edge_max);
    fprintf(fp, " Bounding box x               [%.2f, %.2f]\n",
            info->bounds[0], info->bounds[1]);
    fprintf(fp, "              y               [%.2f, %.2f]\n",
            info->bounds[2], info->bounds[3]);
    fprintf(fp, "              z               [%.2f, %.2f]\n",
            info->bounds[4], info->bounds[5]);
    fprintf(fp, " Centroid                     [%.2f, %.2f, %.2f]\n",
            info->centroid[0], info->centroid[1], info->centroid[2]);

    fprintf(fp, "%s", rule);
    if (info->n_self_intersections == CAT_SURFINFO_SKIPPED) {
        fprintf(fp, " Self-intersections           not checked\n");
    } else {
        double percent = (info->n_polygons > 0) ?
            100.0 * (double) info->n_intersecting_polygons /
                    (double) info->n_polygons : 0.0;

        fprintf(fp, " Self-intersections           %d\n",
                info->n_self_intersections);
        fprintf(fp, " Triangle intersections       %d\n",
                info->n_triangle_hits);
        fprintf(fp, " Intersecting polygons        %d (%.3f%%)\n",
                info->n_intersecting_polygons, percent);
    }
    fprintf(fp, "%s", rule);
}

/**
 * \brief Print a surf_info_struct as tab-separated key/value lines.
 *
 * Machine-readable counterpart of print_surf_info(), one "key<TAB>value" pair
 * per line.  Fields that were not computed are omitted.
 *
 * \param info (in) summary filled by compute_surf_info()
 * \param name (in) value of the "file" key, may be NULL
 * \param fp   (in) output stream
 * \return void (no return value).
 */
void
print_surf_info_tabular(const surf_info_struct *info, const char *name,
                        FILE *fp)
{
    if (name != NULL)
        fprintf(fp, "file\t%s\n", name);

    fprintf(fp, "vertices\t%d\n", info->n_points);
    fprintf(fp, "polygons\t%d\n", info->n_polygons);
    fprintf(fp, "faces\t%d\n", info->n_polygons);
    fprintf(fp, "triangular_faces\t%d\n", info->n_triangles);
    fprintf(fp, "nontriangular_faces\t%d\n", info->n_nontriangular);
    fprintf(fp, "edges\t%d\n", info->n_edges);
    fprintf(fp, "euler\t%d\n", info->euler);
    if (info->closed)
        fprintf(fp, "genus\t%d\n", info->genus);
    fprintf(fp, "components\t%d\n", info->n_components);
    fprintf(fp, "closed\t%d\n", info->closed);
    fprintf(fp, "boundary_edges\t%d\n", info->n_boundary_edges);
    fprintf(fp, "nonmanifold_edges\t%d\n", info->n_nonmanifold_edges);
    fprintf(fp, "degenerate_polygons\t%d\n", info->n_degenerate_polygons);
    fprintf(fp, "unreferenced_vertices\t%d\n", info->n_unreferenced_points);
    fprintf(fp, "surface_area\t%.6f\n", info->surface_area);
    if (info->closed)
        fprintf(fp, "volume\t%.6f\n", fabs(info->volume));
    fprintf(fp, "polygon_area_mean\t%.6f\n", info->area_mean);
    fprintf(fp, "polygon_area_min\t%.6f\n", info->area_min);
    fprintf(fp, "polygon_area_max\t%.6f\n", info->area_max);
    fprintf(fp, "edge_length_mean\t%.6f\n", info->edge_mean);
    fprintf(fp, "edge_length_min\t%.6f\n", info->edge_min);
    fprintf(fp, "edge_length_max\t%.6f\n", info->edge_max);
    fprintf(fp, "bbox_xmin\t%.6f\nbbox_xmax\t%.6f\n",
            info->bounds[0], info->bounds[1]);
    fprintf(fp, "bbox_ymin\t%.6f\nbbox_ymax\t%.6f\n",
            info->bounds[2], info->bounds[3]);
    fprintf(fp, "bbox_zmin\t%.6f\nbbox_zmax\t%.6f\n",
            info->bounds[4], info->bounds[5]);
    fprintf(fp, "centroid_x\t%.6f\ncentroid_y\t%.6f\ncentroid_z\t%.6f\n",
            info->centroid[0], info->centroid[1], info->centroid[2]);

    if (info->n_self_intersections != CAT_SURFINFO_SKIPPED) {
        fprintf(fp, "self_intersections\t%d\n", info->n_self_intersections);
        fprintf(fp, "triangle_intersections\t%d\n", info->n_triangle_hits);
        fprintf(fp, "intersecting_polygons\t%d\n",
                info->n_intersecting_polygons);
    }
}

/**
 * \brief Render the shape of a DataArray as "d0 x d1 x ...".
 *
 * \param info (in)  array description
 * \param buf  (out) destination buffer
 * \param len  (in)  size of the destination buffer
 * \return buf
 */
static const char *
format_dims(const gifti_darray_info *info, char *buf, size_t len)
{
    int d;
    size_t used = 0;

    buf[0] = '\0';
    for (d = 0; d < info->num_dim && d < 6; d++) {
        int written = snprintf(buf + used, len - used, "%s%d",
                               (d == 0) ? "" : " x ", info->dims[d]);
        if (written < 0 || (size_t) written >= len - used)
            break;
        used += (size_t) written;
    }
    return buf;
}

/**
 * \brief Print the DataArrays a GIFTI file carries besides its mesh.
 *
 * A .gii holds the mesh in two DataArrays but may embed per-vertex data next
 * to it -- thickness, curvature, labels -- which the mesh itself does not
 * reveal.  Lists every array with its intent, shape, encoding and value range.
 *
 * \param arrays   (in) descriptions from input_gifti_darrays()
 * \param n_arrays (in) number of entries in arrays
 * \param fp       (in) output stream
 * \return void (no return value).
 */
void
print_gifti_darrays(const gifti_darray_info *arrays, int n_arrays, FILE *fp)
{
    static const char *rule =
        "---------------------------------------------------------------\n";
    char dims[64];
    int i, n_data = 0;

    for (i = 0; i < n_arrays; i++) {
        if (arrays[i].intent != NIFTI_INTENT_POINTSET &&
            arrays[i].intent != NIFTI_INTENT_TRIANGLE)
            n_data++;
    }

    fprintf(fp, " GIFTI data arrays            %d (%d besides the mesh)\n",
            n_arrays, n_data);

    for (i = 0; i < n_arrays; i++) {
        const gifti_darray_info *a = &arrays[i];

        fprintf(fp, "   [%d] %-24s %s\n", i, a->intent_name,
                format_dims(a, dims, sizeof(dims)));
        fprintf(fp, "       %-24s %s, %s\n", "type / encoding",
                a->datatype_name, a->encoding_name);
        /* the mesh arrays inherit the file name as their Name metadata,
           which only repeats what the header already says */
        if (a->name[0] != '\0' && a->intent != NIFTI_INTENT_POINTSET &&
            a->intent != NIFTI_INTENT_TRIANGLE)
            fprintf(fp, "       %-24s %s\n", "name", a->name);
        if (a->ext_fname[0] != '\0')
            fprintf(fp, "       %-24s %s\n", "external file", a->ext_fname);
        if (a->has_range)
            fprintf(fp, "       %-24s %.4f / %.4f / %.4f\n",
                    "value mean/min/max", a->mean, a->min, a->max);
        if (a->n_nonfinite > 0)
            fprintf(fp, "       %-24s %lld of %lld\n", "NaN or infinite",
                    a->n_nonfinite, a->n_values);
    }

    fprintf(fp, "%s", rule);
}

/**
 * \brief Print the GIFTI DataArrays as tab-separated key/value lines.
 *
 * \param arrays   (in) descriptions from input_gifti_darrays()
 * \param n_arrays (in) number of entries in arrays
 * \param fp       (in) output stream
 * \return void (no return value).
 */
void
print_gifti_darrays_tabular(const gifti_darray_info *arrays, int n_arrays,
                            FILE *fp)
{
    char dims[64];
    int i, n_data = 0;

    for (i = 0; i < n_arrays; i++) {
        if (arrays[i].intent != NIFTI_INTENT_POINTSET &&
            arrays[i].intent != NIFTI_INTENT_TRIANGLE)
            n_data++;
    }

    fprintf(fp, "gifti_darrays\t%d\n", n_arrays);
    fprintf(fp, "gifti_data_darrays\t%d\n", n_data);

    for (i = 0; i < n_arrays; i++) {
        const gifti_darray_info *a = &arrays[i];

        fprintf(fp, "darray%d_intent\t%s\n", i, a->intent_name);
        fprintf(fp, "darray%d_dims\t%s\n", i,
                format_dims(a, dims, sizeof(dims)));
        fprintf(fp, "darray%d_datatype\t%s\n", i, a->datatype_name);
        fprintf(fp, "darray%d_encoding\t%s\n", i, a->encoding_name);
        if (a->name[0] != '\0' && a->intent != NIFTI_INTENT_POINTSET &&
            a->intent != NIFTI_INTENT_TRIANGLE)
            fprintf(fp, "darray%d_name\t%s\n", i, a->name);
        if (a->ext_fname[0] != '\0')
            fprintf(fp, "darray%d_external_file\t%s\n", i, a->ext_fname);
        if (a->has_range) {
            fprintf(fp, "darray%d_mean\t%.6f\n", i, a->mean);
            fprintf(fp, "darray%d_min\t%.6f\n", i, a->min);
            fprintf(fp, "darray%d_max\t%.6f\n", i, a->max);
        }
        fprintf(fp, "darray%d_nonfinite\t%lld\n", i, a->n_nonfinite);
    }
}
