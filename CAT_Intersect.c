/* Rachel Yotter - rachel.yotter@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
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
#include "CAT_Intersect.h"


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
 *             0 = disjoint (no intersect)
 *             1 = intersect in unique point I1
 *             2 = are in the same plane
 */
int
intersect_segment_triangle(Point p0, Point p1, int tpidx[3],
                           polygons_struct *surface)
{
        Vector   u, v, n;             /* triangle vectors */
        Vector   dir, w0, w;          /* ray vectors */
        Vector   zero;
        float    r, a, b;             /* params to calc ray-plane intersect */
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
find_selfintersections(polygons_struct *surface, int *defects, int *polydefects)
{
        int                n_intersects, size, i, n, p, p2, b;
        progress_struct    progress;
        struct octree      *tree;
        struct polynode    *cur, *node;

        memset(polydefects, 0, sizeof(int) * surface->n_items);

        tree = build_octree(surface);

        initialize_progress_report(&progress, FALSE, surface->n_items,
                                   "find_selfintersections");

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
                                polydefects[node->num] = n_intersects;
                                polydefects[cur->num] = n_intersects;
                        }

                }
                update_progress_report(&progress, p);
        }

        terminate_progress_report(&progress);

        delete_octree(tree);

        update_defects(surface, polydefects, defects);

        return n_intersects;
}


/* consolidate neighboring self-intersections */
int
join_intersections(polygons_struct *surface, int *defects, int *polydefects,
                   int *n_neighbours, int **neighbours)
{
        int                d, old_d;
        int                n_intersects = 0, i, n, p, *dmap;

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
        int *edgeflag;
        Point tp[3];
        double areas[128], centers[384], xyz[3], weight, t_area;
        FILE *fp;

        if (n_defects == 0)
                 return; /* done! */

        update_defects(surface, polydefects, defects);

        edgeflag = (int *) malloc(sizeof(int) * surface->n_points);

        /* smooth defect areas.. increase defect area every 5th iter */
        iter = 0;
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
                printf("iter %3d: %3d smoothed pts, %2d intersection(s) remaining\n",
                       iter, npts, n_defects);

                if (n_defects == 0)
                        break; /* all done! */

                if (npts > 100) {
                        /* remap defects to limit # of affected points */
                        n_defects = find_remaining_intersections(surface,
                                                                 defects,
                                                                 polydefects,
                                                                 n_neighbours,
                                                                 neighbours);
                        printf("defects remapped, %d defect(s) remaining\n",
                               n_defects);
                } else if (iter % 5 == 0) {
                        /* expand remaining defects every 5th iteration */
                        expand_defects(surface, defects, polydefects,
                                       0, 1, n_neighbours, neighbours);
                }
        }

        n_defects = find_remaining_intersections(surface, defects, polydefects,
                                                 n_neighbours, neighbours);

        free(edgeflag);
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
