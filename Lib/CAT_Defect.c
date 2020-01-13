/* Rachel Yotter - rachel.yotter@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 *
 */


#include "CAT_SPH.h"
#include "CAT_Intersect.h"
#include "CAT_Defect.h"
#include "CAT_Surf.h"
#include "CAT_Curvature.h"
#include "CAT_SurfaceIO.h"

#define DUMP_FILES 0

Vector
defect_direction(polygons_struct *surface, int *defects, int defect)
{
        Vector dir;
        int p;

        /* get the direction of the defect */
        fill_Vector(dir, 0.0, 0.0, 0.0);
        for (p = 0; p < surface->n_points; p++) {
                if (defects[p] == defect)
                        ADD_POINT_VECTOR(dir, dir, surface->normals[p]);
        }
        NORMALIZE_VECTOR(dir, dir);
        return dir;
}

int
defect_euler(polygons_struct *surface, int *defects, int *polydefects,
             int defect, int *n_neighbours, int **neighbours)
{
        int p, n, nn, n_v, n_e, n_f, n_dup_e;
        int free_polydefects = 0;

        if (polydefects == NULL) {
                polydefects = (int *) malloc(sizeof(int) * surface->n_items);
                update_polydefects(surface, defects, polydefects);
                free_polydefects = 1;
        }

        n_v = 0; n_e = 0; n_f = 0; n_dup_e = 0;
        for (p = 0; p < surface->n_points; p++) {
                if (defects[p] != defect)
                        continue;

                n_v++; /* found a vertex */

                for (n = 0; n < n_neighbours[p]; n++) {
                        if (defects[neighbours[p][n]] != defect)
                                continue;
                        for (nn = n+1; nn < n_neighbours[p]; nn++) {
                                if (neighbours[p][nn] == neighbours[p][n])
                                        break;
                        }
                        if (nn < n_neighbours[p])
                                n_dup_e++;
                        else if (p < neighbours[p][n])
                                n_e++;
                }
        }
        for (p = 0; p < surface->n_items; p++) {
                if (polydefects[p] == defect)
                        n_f++; /* found a face */
        }

        if (free_polydefects)
                free(polydefects);

        return (n_v + n_f - n_e);
}

int
isedge(polygons_struct *surface, int *defects, int *polydefects,
       int *n_neighbours, int **neighbours, int polygon)
{
        int i, p, n;

        for (i = 0; i < 3; i++) {
                p = surface->indices[POINT_INDEX(surface->end_indices,
                                                 polygon, i)];

                for (n = 0; n < n_neighbours[p]; n++) {
                        if (defects[neighbours[p][n]] != defects[p]) {
                                return 1; /* an edge! */
                        }
                }
        }
        return 0;
}

int
find_topological_defects(polygons_struct *surface, polygons_struct *sphere,
                         int *defects, int *n_neighbours, int **neighbours)
{
        int d, pts[3], n_defects, n_euler, n_intersects;
        int size, i, n, p;
        int euler, *polydefects;
        struct looptree *ltree;

        polydefects = (int *) malloc(sizeof(int) * sphere->n_items);

        n_intersects = find_selfintersections(sphere, defects, polydefects);
        if (n_intersects == 0)
                return 0;

        n_intersects = join_intersections(sphere, defects, polydefects,
                                          n_neighbours, neighbours);

        /* remove simple self-intersections */
        for (d = 1; d <= n_intersects; d++) {
                euler = defect_euler(surface, defects, polydefects, d,
                                     n_neighbours, neighbours);

                if (euler == 1) { /* not a defect, delete it */
                        for (p = 0; p < sphere->n_points; p++) {
                                if (defects[p] == d)
                                        defects[p] = 0;
                                if (defects[p] == n_intersects)
                                        defects[p] = d;
                        }
                        for (p = 0; p < sphere->n_items; p++) {
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
expand_defects(polygons_struct *surface, int *defects, int *polydefects,
              int defect, int level, int *n_neighbours, int **neighbours)
{
        int *buffer, p, n;

        buffer = (int *) malloc(sizeof(int) * surface->n_points);
        update_defects(surface, polydefects, defects);

        while (level > 0) {
                for (p = 0; p < surface->n_points; p++)
                        buffer[p] = defects[p];

                for (p = 0; p < surface->n_points; p++) {
                        if (buffer[p] == 0)
                                continue; /* skip */

                        if (defect > 0 && buffer[p] != defect)
                                continue; /* skip */

                        for (n = 0; n < n_neighbours[p]; n++) {
                                defects[neighbours[p][n]] = buffer[p];
                        }
                }
                level--;
        }
        update_polydefects(surface, defects, polydefects);

        free(buffer);
}


/* Update the defect list from the polygon defect list */
void
update_defects(polygons_struct *surface, int *polydefects, int *defects)
{
        int p;

        memset(defects, 0, sizeof(int) * surface->n_points);

        /* update point defect list */
        for (p = 0; p < surface->n_items; p++) {
                if (polydefects[p] == 0) continue;
                defects[surface->indices[POINT_INDEX(surface->end_indices,
                                                     p, 0)]] = polydefects[p];
                defects[surface->indices[POINT_INDEX(surface->end_indices,
                                                     p, 1)]] = polydefects[p];
                defects[surface->indices[POINT_INDEX(surface->end_indices,
                                                     p, 2)]] = polydefects[p];
        }
}


/* Update the polygon defect list from the point defect list */
void
update_polydefects(polygons_struct *surface, int *defects, int *polydefects)
{
        int p, d;

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


/* returns the bounds of a defect */
void
get_defect_bounds(polygons_struct *surface, int *defects, int defect,
                  double bounds[6])
{
        int p;

        bounds[0] = bounds[2] = bounds[4] = PINF;
        bounds[1] = bounds[3] = bounds[5] = NINF;

        for (p = 0; p < surface->n_points; p++) {
                if (defects[p] != defect)
                        continue;

                if (Point_x(surface->points[p]) < bounds[0])
                        bounds[0] = Point_x(surface->points[p]);
                if (Point_x(surface->points[p]) > bounds[1])
                        bounds[1] = Point_x(surface->points[p]);
                if (Point_y(surface->points[p]) < bounds[2])
                        bounds[2] = Point_y(surface->points[p]);
                if (Point_y(surface->points[p]) > bounds[3])
                        bounds[3] = Point_y(surface->points[p]);
                if (Point_z(surface->points[p]) < bounds[4])
                        bounds[4] = Point_z(surface->points[p]);
                if (Point_z(surface->points[p]) > bounds[5])
                        bounds[5] = Point_z(surface->points[p]);
        }
}


/* holes = 1, handles = 2, large errors = 3, ventricles = 4 */

/* returns the center point of a defect */
Point
get_defect_center(polygons_struct *surface, int *defects, int defect)
{
        Point center;
        int p, npts = 0;

        fill_Point(center, 0.0, 0.0, 0.0);

        for (p = 0; p < surface->n_points; p++) {
                if (defects[p] != defect)
                        continue;

                ADD_POINTS(center, center, surface->points[p]);
                npts++;
        }
        Point_x(center) = Point_x(center) / npts;
        Point_y(center) = Point_y(center) / npts;
        Point_z(center) = Point_z(center) / npts;

        return center;
}

/* returns size of detect in relation to overall surface size */
void
get_defect_size(polygons_struct *surface, int *defects, int n_defects, 
                  double *defect_size)
{
        int i, p, d;
        double *size;

        size = (double *) malloc(sizeof(double) * n_defects);

        memset(defect_size, 0, sizeof(double) * surface->n_points);

        /* get size of defect */
        for (d = 1; d <= n_defects; d++) {
                size[d] = 0.0;
                for (p = 0; p < surface->n_points; p++)
                        if (defects[p] == d)
                                size[d]++;
        }

        for (d = 1; d <= n_defects; d++)
                size[d] /= (double)surface->n_points;

        for (d = 1; d <= n_defects; d++)
                for (p = 0; p < surface->n_points; p++)
                        if (defects[p] == d) 
                                    defect_size[p] = size[d];

        free(size);
}

/* cut defects in half, saving top half of holes and bottom half of handles based on sulcal depth of inflated surface */
void
bisect_defects(polygons_struct *surface, polygons_struct *sphere, int *defects, int n_defects,
               int *holes, int *bisected, int detect_euler)
{
        int p, d, i, fillflag, euler;
        int npts, maxpts, *pts, *polydefects;
        int *n_neighbours, **neighbours;
        Vector dir;
        double *depth, avg, sum;
        polygons_struct *inflated_surface;
        object_struct **objects;

        memset(bisected, 0, sizeof(int) * surface->n_points);

        pts = (int *) malloc(sizeof(int) * surface->n_points);
        depth = (double *) malloc(sizeof(double) * surface->n_points);

        /* inflate surface roughly to an ellopsoid where holes and handles are still visible and have
           different height/depth on the surface */
        objects = (object_struct **) malloc(sizeof(object_struct *));
        *objects = create_object(POLYGONS);
        inflated_surface = get_polygons_ptr(*objects);
        copy_polygons(surface, inflated_surface);
        inflate_surface_with_topology_defects(inflated_surface);

        /* estimate sulcal depth of inflated surface to identify holes and handles due to their different
           height/depth */
        memset(depth, 0, sizeof(double) * surface->n_points);
        compute_sulcus_depth(inflated_surface, depth);

        if (DUMP_FILES) {
                output_values_any_format("depth.txt", surface->n_points, 
                                           depth, TYPE_DOUBLE);
                output_graphics_any_format("inflated.obj", ASCII_FORMAT,
                                           1, objects, NULL);
        }
        
        /* estimate overall depth to check whether estimation was correct */
        sum = 0.0;
        for (p = 0; p < surface->n_points; p++)  
                sum += depth[p];
                
        /* use Hausdorff distance if sulcal depth estimation fails */
        if (sum == 0.0) {
                printf("WARNING: Hausdorff distance is used for bisectioning the defects because estimation of sulcal depth failed.\n");
                compute_point_hausdorff(inflated_surface, sphere, depth, 0);
        }

        /* prepare calculation of euler number to decide whether we have a handle or hole */
        if (detect_euler) {
                create_polygon_point_neighbours(sphere, FALSE, &n_neighbours,
                                        &neighbours, NULL, NULL);

                polydefects = (int *) malloc(sizeof(int) * sphere->n_items);
                update_polydefects(sphere, defects, polydefects);
        }

        for (d = 1; d <= n_defects; d++) {
                npts = 0;
                for (p = 0; p < surface->n_points; p++) {
                        if (defects[p] == d) {
                                pts[npts] = p;
                                npts++;
                        }
                }

                /* estimate mean sulcal depth inside defect */
                avg = 0.0;
                for (p = 0; p < npts; p++)  avg += depth[pts[p]];
                avg /= (double)npts;

                /* we can use euler number to automatically decide whether we have a handle or hole */
                if (detect_euler) {
                        euler = defect_euler(surface, defects, polydefects, d,
                                     n_neighbours, neighbours);
                        fillflag = (euler < 1) ? HOLE : HANDLE;
                } else fillflag = holes[pts[0]]; /* otherwise use predefined decision */

                /* indicate "bottom" (ground) of defect only as either hole or handle */
                for (p = 0; p < npts; p++)
                        if (((fillflag == HOLE) & (depth[pts[p]] >= avg)) | ((fillflag == HANDLE) & (depth[pts[p]] < avg)))
                                bisected[pts[p]] = fillflag;
                       
        }
        free(pts);
        free(depth);
        free(objects);
        if (detect_euler) free(polydefects);

}

void
inflate_surface_with_topology_defects(polygons_struct *polygons)
{
        int              i;
        double           factor;
        

        /* use more iterations for larger surfaces */
        if (polygons->n_items > 350000) {
                factor = (double)polygons->n_items/350000.0;
                printf("Large number polygons -> Increase # of iterations by factor %g.\n",
                            factor);
        } else factor = 1.0;

        /* low smooth */
        inflate_surface_and_smooth_fingers(polygons,
                  /*                     cycles */ 1,
                  /* regular smoothing strength */ 0.2,
                  /*    regular smoothing iters */ round(factor*50),
                  /*           inflation factor */ 1.0,
                  /* finger comp/stretch thresh */ 3.0,
                  /*     finger smooth strength */ 1.0,
                  /*        finger smooth iters */ 0);

                /* inflated */
        inflate_surface_and_smooth_fingers(polygons,
                  /*                     cycles */ 2,
                  /* regular smoothing strength */ 1.0,
                  /*    regular smoothing iters */ round(factor*30),
                  /*           inflation factor */ 1.4,
                  /* finger comp/stretch thresh */ 3.0,
                  /*     finger smooth strength */ 1.0,
                  /*        finger smooth iters */ 0);

                /* very inflated */
        inflate_surface_and_smooth_fingers(polygons,
                  /*                     cycles */ 4,
                  /* regular smoothing strength */ 1.0,
                  /*    regular smoothing iters */ round(factor*30),
                  /*           inflation factor */ 1.1,
                  /* finger comp/stretch thresh */ 3.0,
                  /*     finger smooth strength */ 1.0,
                  /*        finger smooth iters */ 0);

               /* high smooth */
        inflate_surface_and_smooth_fingers(polygons,
                  /*                     cycles */ 3,
                  /* regular smoothing strength */ 1.0,
                  /*    regular smoothing iters */ round(factor*40),
                  /*           inflation factor */ 1.4,
                  /* finger comp/stretch thresh */ 3.0,
                  /*     finger smooth strength */ 1.0,
                  /*        finger smooth iters */ 60);

        compute_polygon_normals(polygons);
}

/* remap topological defects from the original spherical mapping to a new
 * spherical mapping
 */
void
remap_defect(polygons_struct *sphere, int *defects, int *polydefects,
             polygons_struct *remap, int *remap_defects, int *remap_polydefects)
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

        update_polydefects(remap, remap_defects, remap_polydefects);
}

int
find_artifacts(polygons_struct *surface, polygons_struct *sph,
               int *artifacts, int *n_neighbours, int **neighbours, double dist)
{
        int d, euler, aval, aval2;
        int size, i, n, p, poly, done;
        int n_artifacts, *polyart;
        double hd;
        Point closest;

        /* estimate the Hausdorff distance */
        n_artifacts = 0;
        create_polygons_bintree(sph, ROUND((double) sph->n_items * 0.5));

        for (i = 0; i < surface->n_points; i++) {
                poly = find_closest_polygon_point(&surface->points[i],
                                                  sph, &closest);
                hd = distance_between_points(&surface->points[i], &closest);

                if (hd > dist) {
                        artifacts[i] = (++n_artifacts);
                } else {
                        artifacts[i] = 0;
                }
        }

        /* consolidate artifact points */
        done = 0;
        while (!done) {
                done = 1;
                for (i = 0; i < surface->n_points; i++) {
                        if (artifacts[i] == 0) continue;
                        aval = artifacts[i];

                        for (n = 0; n < n_neighbours[i]; n++) {
                                 aval2 = artifacts[neighbours[i][n]];
                                 if (aval2 == 0 || aval2 == aval) continue;

                                 done = 0;
                                 if (aval2 < aval)
                                         artifacts[i] = aval2;
                                 else
                                         artifacts[neighbours[i][n]] = aval;
                        }
                }
        }

        /* renumber */
        n_artifacts++;
        for (i = 0; i < surface->n_points; i++) {
                if (artifacts[i] == 0 || artifacts[i] < n_artifacts) continue;
                aval = (++n_artifacts);

                for (p = 0; p < surface->n_points; p++) {
                         if (artifacts[p] == artifacts[i])
                                 artifacts[p] = n_artifacts;
                }
                artifacts[i] = n_artifacts;
        }

        /* delete if < 4 points */
        for (i = 1; i <= n_artifacts; i++) {
                size = 0;

                for (p = 0; p < surface->n_points; p++) {
                        if (artifacts[p] == i) {
                                size++;
                        }
                }
                if (size <= 4) { /* delete it */
                        for (p = 0; p < surface->n_points; p++) {
                                if (artifacts[p] == i) {
                                        artifacts[p] = 0;
                                } else if (artifacts[p] == n_artifacts) {
                                        artifacts[p] = i;
                                }
                        }
                        n_artifacts--;
                        i--;
                }
        }

        polyart = (int *) malloc(sizeof(int) * surface->n_items);
        update_polydefects(surface, artifacts, polyart);
        expand_defects(surface, artifacts, polyart, 0, 3,
                       n_neighbours, neighbours);

        for (i = 1; i <= n_artifacts; i++) {
                euler = defect_euler(surface, artifacts, polyart, i,
                                     n_neighbours, neighbours);

                if (euler < 0) { /* not an artifact, delete it */
                        for (p = 0; p < surface->n_points; p++) {
                                if (artifacts[p] == i) {
                                        artifacts[p] = 0;
                                } else if (artifacts[p] == n_artifacts) {
                                        artifacts[p] = i;
                                }
                        }

                        for (p = 0; p < surface->n_items; p++) {
                                if (polyart[p] == i)
                                        polyart[p] = 0;
                                if (polyart[p] == n_artifacts)
                                        polyart[p] = i;
                        }
                        n_artifacts--;
                        i--;
                }
        }

        free(polyart);
        return n_artifacts;
}
