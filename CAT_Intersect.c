/* Rachel Yotter - rachel.yotter@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id: CAT_SelfIntersections.c 121 2009-07-27 14:42:19Z raytrace $
 *
 */

#include <volume_io/internal_volume_io.h>
#include <volume_io/geometry.h>
#include <bicpl.h>
#include <ParseArgv.h>

#include "CAT_Curvature.h"
#include "CAT_Blur2d.h"
#include "CAT_Octree.h"
#include "CAT_SurfaceIO.h"

#define PINF  1.7976931348623157e+308 /* for doubles */
#define NINF -1.7976931348623157e+308 /* for doubles */


int
intersect_triangle_triangle(int pidx0[3], int pidx1[3],
                            polygons_struct *surface)
{
        int i, result;
        Point pts[3];

        /* test if neighbors... if so, skip */
        for (i = 0; i < 3; i++) {
                if (pidx0[i] == pidx1[0]) return 0;
                if (pidx0[i] == pidx1[1]) return 0;
                if (pidx0[i] == pidx1[2]) return 0;
        }

        for (i = 0; i < 3; i++)
                pts[i] = surface->points[pidx0[i]];

        for (i = 1; i < 3; i++) {
                result = intersect_segment_triangle(pts[i-1], pts[i], pidx1,
                                                    surface);
                if (result > 0)
                        return 1;
        }

        return 0;
}


/*
 * intersect_segment_triangle(): test if a segment intersects a triangle
 *    Input:  two points for the segment p0 p1, triangle indices tpidx[3]
 *    Return: -1 = triangle is degenerate (a segment or point)
 *             0 = disjoint (no intersect)
 *             1 = intersect in unique point I1
 *             2 = are in the same plane
 */
int
intersect_segment_triangle(Point p0, Point p1, int tpidx[3],
                           polygons_struct *surface)
{
        Vector   u, v, n;             // triangle vectors
        Vector   dir, w0, w;          // ray vectors
        Vector   zero;
        float    r, a, b;             // params to calc ray-plane intersect
        float    uu, uv, vv, wu, wv, D;
        float    s, t;
        Point    pts[3], I;
        int      i;

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
                } else return 0;             /* ray disjoint from plane */
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

        // get and test parametric coords
        s = (uv * wv - vv * wu) / D;
        if (s < 0.0 || s > 1.0) /* I is outside T */
                return 0;
        t = (uv * wu - uu * wv) / D;
        if (t < 0.0 || (s + t) > 1.0)  /* I is outside T */
                return 0;

        return 1; /* I is in T */
}

int
find_selfintersections(polygons_struct *surface, int *defects, int *polydefects)
{
        int                n_intersects, size, i, n, p, p2, b;
        progress_struct    progress;
        struct octree      *tree;
        struct polynode    *cur, *node;

        memset(defects, 0, sizeof(int) * surface->n_points);
        memset(polydefects, 0, sizeof(int) * surface->n_items);

        tree = build_octree(surface);

        initialize_progress_report(&progress, FALSE, surface->n_items,
                                   "Self-Intersect Test");

        n_intersects = 0;
        for (p = 0; p < surface->n_items; p++) {
                node = tree->nodelist[p];

                for (p2 = p+1; p2 < surface->n_items; p2++)
                        tree->polyflag[p2] = 0;

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

                        for (cur = tree->nodes[b];
                             cur != NULL; cur = cur->next) {
                                if (tree->polyflag[cur->num] == 1)
                                        continue;

                                if (cur->num <= p)
                                        continue;

                                tree->polyflag[cur->num] = 1;

                                if (intersect(node->bounds, cur->bounds) == 0)
                                        continue;

                                if (intersect_triangle_triangle(node->pts,
                                                                cur->pts,
                                                                surface) == 0)
                                        continue;

                                n_intersects++;

                                polydefects[p] = n_intersects;
                                for (i = 0; i < 3; i++)
                                        defects[node->pts[i]] = n_intersects;

                                n_intersects++;
                                polydefects[cur->num] = n_intersects;
                                for (i = 0; i < 3; i++)
                                        defects[cur->pts[i]] = n_intersects;
                        }

                }
                update_progress_report(&progress, p);
        }

        terminate_progress_report(&progress);

        delete_octree(tree);

        return n_intersects/2;
}

/* consolidate neighboring self-intersections */
int
join_intersections(polygons_struct *surface, int *defects,
                   int *n_neighbours, int **neighbours)
{
        int                d, old_d;
        int                n_intersects = 0, i, n, p, *dmap;

        for (i = 0; i < surface->n_points; i++) {
                if (defects[i] == 0)
                        continue; /* skip */

                for (n = 0; n < n_neighbours[i]; n++) {
                        d = defects[neighbours[i][n]];
                        if (d > 0 && d != defects[i]) {
                                old_d = d > defects[i] ?
                                        d : defects[i];
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

        free(dmap);
        return n_intersects;
}

int
find_topological_defects(polygons_struct *surface, int *defects)
{
        int                *n_neighbours, **neighbours;
        int                d, pts[3];
        int                n_intersects = 0, size, i, n, p;
        int                n_v, n_e, n_f;
        int                *polydefects;

        polydefects = (int *) malloc(sizeof(int) * surface->n_items);

        n_intersects = find_selfintersections(surface, defects, polydefects);

        printf("Raw Triangle Intersects %d\n", n_intersects);

        if (n_intersects == 0) return 0; /* done! */

        /* consolidate intersections */
        create_polygon_point_neighbours(surface, TRUE, &n_neighbours,
                                        &neighbours, NULL, NULL);
        n_intersects = join_intersections(surface, defects,
                                          n_neighbours, neighbours);

        /* find the matching intersecting polygons */
        for (i = 0; i < surface->n_items; i++) {
                if (polydefects[i] == 0)
                        continue;

                pts[0] = surface->indices[
                                      POINT_INDEX(surface->end_indices, i, 0)];
                pts[1] = surface->indices[
                                      POINT_INDEX(surface->end_indices, i, 1)];
                pts[2] = surface->indices[
                                      POINT_INDEX(surface->end_indices, i, 2)];

                if (defects[pts[0]] == defects[pts[1]] &&
                    defects[pts[0]] == defects[pts[2]]) {
                        polydefects[i] = defects[pts[0]];
                } else
                        polydefects[i] = 0;
        }

        /* calculate the Euler number of each defect, delete if = 1 */
        for (d = 1; d <= n_intersects; d++) {
                n_v = 0; n_e = 0; n_f = 0;
                for (p = 0; p < surface->n_points; p++) {
                        if (defects[p] != d)
                                continue;

                        n_v++; /* found a vertex */

                        for (n = 0; n < n_neighbours[p]; n++) {
                                if (defects[neighbours[p][n]] == d &&
                                    p < neighbours[p][n]) {
                                        n_e++; /* found an edge */
                                }
                        }
                }
                for (p = 0; p < surface->n_items; p++) {
                        if (polydefects[p] == d)
                                n_f++; /* found a face */
                }

                if (n_v + n_f - n_e == 1) { /* not a defect, delete it */
                        for (p = 0; p < surface->n_points; p++) {
                                if (defects[p] == d)
                                        defects[p] = 0;
                                if (defects[p] == n_intersects)
                                        defects[p] = d;
                        }
                        for (p = 0; p < surface->n_items; p++) {
                                if (polydefects[p] == d)
                                        polydefects[p] = 0;
                                if (polydefects[p] == n_intersects)
                                        polydefects[p] = d;
                        }
                        n_intersects--;
                        d--;
                }
        }

        free(polydefects);

        return n_intersects;
}

/* expand defects to include neighboring points.  defect = 0 for all defects,
 * or can just expand one defect.
 */
void
expand_defects(polygons_struct *surface, int *defects, int defect, int level)
{
        int *n_neighbours, **neighbours;
        int *buffer, p, n;

        create_polygon_point_neighbours(surface, TRUE, &n_neighbours,
                                        &neighbours, NULL, NULL);

        buffer = (int *) malloc(sizeof(int) * surface->n_points);

        while (level > 0) {
                for (p = 0; p < surface->n_points; p++)
                        buffer[p] = defects[p];

                for (p = 0; p < surface->n_points; p++) {
                        if (buffer[p] == 0)
                                continue;

                        if (defect > 0 && buffer[p] != defect)
                                continue;

                        for (n = 0; n < n_neighbours[p]; n++) {
                                defects[neighbours[p][n]] = defects[p];
                        }
                }
                level--;
        }

        free(buffer);
}
