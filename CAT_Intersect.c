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
#include "CAT_Intersect.h"


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
find_selfintersections(polygons_struct *surface, int *defects)
{
        int                n_intersects, size, i, n, p, p2, b;
        progress_struct    progress;
        struct octree      *tree;
        struct polynode    *cur, *node;

        memset(defects, 0, sizeof(int) * surface->n_points);

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
                                for (i = 0; i < 3; i++)
                                        defects[node->pts[i]] = n_intersects;

                                n_intersects++;
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

/* consolidate neighboring self-intersections */
int
join_intersections_poly(polygons_struct *surface, int *defects,
                        int *polydefects, int *n_neighbours, int **neighbours)
{
        join_intersections(surface, defects, n_neighbours, neighbours);
        update_polydefects(surface, defects, polydefects);
}

int
find_topological_defects(polygons_struct *surface, int *defects,
                         int *n_neighbours, int **neighbours)
{
        int                d, pts[3];
        int                n_intersects = 0, size, i, n, p;
        int                n_v, n_e, n_f;
        int                *polydefects;

        n_intersects = find_selfintersections(surface, defects);
        if (n_intersects == 0)
                return 0; /* done! */

        /* consolidate intersections */
        n_intersects = join_intersections(surface, defects,
                                          n_neighbours, neighbours);

        printf("%d Topological Defects\n", n_intersects);

        polydefects = (int *) malloc(sizeof(int) * surface->n_items);
        update_polydefects(surface, defects, polydefects);

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
expand_defects(polygons_struct *surface, int *defects, int defect, int level,
               int *n_neighbours, int **neighbours)
{
        int *buffer, p, n;

        buffer = (int *) malloc(sizeof(int) * surface->n_points);

        while (level > 0) {
                for (p = 0; p < surface->n_points; p++)
                        buffer[p] = defects[p];

                for (p = 0; p < surface->n_points; p++) {
                        if (buffer[p] == 0)
                                continue; /* skip */

                        if (defect > 0 && buffer[p] != defect)
                                continue; /* skip */

                        for (n = 0; n < n_neighbours[p]; n++) {
                                defects[neighbours[p][n]] = defects[p];
                        }
                }
                level--;
        }

        free(buffer);
}

/* expand defects to include neighboring points.  defect = 0 for all defects,
 * or can just expand one defect.
 */
void
expand_defects_poly(polygons_struct *surface, int *defects, int *polydefects,
                    int defect, int level, int *n_neighbours, int **neighbours)
{
        expand_defects(surface, defects, defect, level,
                       n_neighbours, neighbours);
        update_polydefects(surface, defects, polydefects);
}

/* Update the polygon defect list from the point defect list */
void
update_polydefects(polygons_struct *surface, int *defects, int *polydefects)
{
        int p;

        /* update polygon defect list */
        for (p = 0; p < surface->n_items; p++) {
                polydefects[p] = defects[surface->indices[
                                    POINT_INDEX(surface->end_indices, p, 0)]];
                if (polydefects[p] == 0) continue;
                polydefects[p] = defects[surface->indices[
                                    POINT_INDEX(surface->end_indices, p, 1)]];
                if (polydefects[p] == 0) continue;
                polydefects[p] = defects[surface->indices[
                                    POINT_INDEX(surface->end_indices, p, 2)]];
        }
}


/* holes = 1, handles = 2, in between = 3 */
/* this code doesn't work correctly */
void
identify_holes_handles(polygons_struct *surface, int *defects, int n_defects,
                       double *curvatures, int *n_neighbours, int **neighbours)
{
        int                d, p, *count;
        double             *curv_score;

        if (n_defects == 0) return;

        get_polygon_vertex_curvatures_cg(surface, n_neighbours, neighbours,
                                         3, 4, curvatures); /* mean curv */

        curv_score = (double *) malloc(sizeof(double) * n_defects);
        count = (int *) malloc(sizeof(int) * n_defects);
        memset(curv_score, 0, sizeof(double) * n_defects);
        memset(count, 0, sizeof(int) * n_defects);

        for (p = 0; p < surface->n_points; p++) {
                if (defects[p] > 0) {
                        curv_score[defects[p]-1] += curvatures[p];
                        count[defects[p]-1]++;
                }
        }

        for (d = 0; d < n_defects; d++) {
                curv_score[d] /= count[d];
                if (curv_score[d] > 5)
                        curv_score[d] = 1; /* hole */
                else if (curv_score[d] < 3)
                        curv_score[d] = 2; /* handle */
                else 
                        curv_score[d] = 3;
        }


        for (p = 0; p < surface->n_points; p++) {
                if (defects[p] > 0) {
                        curvatures[p] = curv_score[defects[p]-1];
                } else curvatures[p] = 0;
        }

        free(curv_score);
        free(count);
}

/* remap topological defects from the original spherical mapping to a new
 * spherical mapping
 */
void
remap_defect_points(polygons_struct *sphere, int *defects,
                    polygons_struct *remap, int *remap_defects)
{
        int p, size, i, idx, poly;
        Point closest_pt;
        progress_struct progress;

        if (remap->bintree == NULL) {
                create_polygons_bintree(remap, round((double) remap->n_items *
                                              BINTREE_FACTOR));
        }

        initialize_progress_report(&progress, FALSE, remap->n_points,
                                   "remap_defect");

        for (p = 0; p < sphere->n_points; p++) {
                if (defects[p] == 0) continue; /* skip */

                poly = find_closest_polygon_point(&sphere->points[p], remap,
                                                  &closest_pt);
                size = GET_OBJECT_SIZE(*remap, poly);
                for (i = 0; i < size; i++) {
                        idx = remap->indices[POINT_INDEX(remap->end_indices,
                                                  poly, i)];
                        remap_defects[idx] = defects[p];
                }
                update_progress_report(&progress, p);
        }

        terminate_progress_report(&progress);
}

void
repair_selfintersections(polygons_struct *surface, int *n_neighbours,
                      int **neighbours)
{
        int p, i, iter, d, n, n2, npts;
        int *defects, *polydefects, *edgeflag;
        int n_defects;
        Point *pts, tp[3];
        double areas[128], centers[128], xyz[3], weight, t_area;
        FILE *fp;

        defects = (int *) malloc(sizeof(int) * surface->n_points);

        n_defects = find_selfintersections(surface, defects);
        n_defects = join_intersections(surface, defects,
                                       n_neighbours, neighbours);
        printf("%d defect(s) to repair\n", n_defects);

        /* optionally dump self intersections */ /*
        if (open_file("si.txt", WRITE_FILE, ASCII_FORMAT, &fp) != OK) {
                exit(0);
        }

        for (p = 0; p < surface->n_points; p++)
                fprintf(fp, " %d.0\n", defects[p]);
        fclose(fp);
        */

        if (n_defects == 0) {
                 free(defects);
                 return;
        }

        polydefects = (int *) malloc(sizeof(int) * surface->n_items);
        update_polydefects(surface, defects, polydefects);

        edgeflag = (int *) malloc(sizeof(int) * surface->n_points);
        pts = (Point *) malloc(sizeof(Point) * surface->n_points);
        for (p = 0; p < surface->n_points; p++)
                pts[p] = surface->points[p];

        /* smooth defect areas.. increase defect area every 3rd iter */
        iter = 0;
        while (n_defects != 0) {
                iter++;
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

                /* test if self-intersections repaired */
                for (d = 1; d <= n_defects; d++) {
                        if (defect_has_intersections(surface, polydefects,
                                                     d) == 0) {
                                /* delete it, it's fixed! */
                                for (i = 0; i < surface->n_points; i++) {
                                        if (defects[i] == d)
                                                defects[i] = 0;
                                        else if (defects[i] == n_defects)
                                                defects[i] = d;
                                }
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
                printf("iter %2d: %3d smoothed pts, %2d defect(s) remaining\n",
                       iter, npts, n_defects);

                if (n_defects == 0)
                        break; /* all done! */

                /* expand remaining defects every 3rd iteration */
                if (iter % 3 == 0) {
                        expand_defects_poly(surface, defects, polydefects,
                                            0, 1, n_neighbours, neighbours);
                }
        }

        free(defects);
        free(polydefects);
        free(edgeflag);
        free(pts);
}

int
defect_has_intersections(polygons_struct *polygons, int *polydefects,
                         int defect)
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
