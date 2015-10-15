/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id: CAT_SurfParameters.c 333 2015-01-27 10:46:13Z gaser $
 *
 */


#include <bicpl.h>
#include <float.h>

#include "CAT_Surf.h"
#include "CAT_ConvexHull.h"

#define _PI 3.14159265358979323846264338327510
#define PI2 0.6366197724    /* 2/pi */
#define  MAX_NEIGHBOURS   1000

void
compute_sulcus_depth(polygons_struct *surface, polygons_struct *sphere, double *depth)
{
        polygons_struct *convex;
        object_struct **object;
        Point closest;
        int i, poly;

        object = (object_struct **) malloc(sizeof(object_struct *));
        *object = create_object(POLYGONS);
        
        /* get convex hull */
        object = surface_get_convex_hull(surface, NULL);
        convex = get_polygons_ptr(*object);

        if (convex->bintree == NULL) 
                create_polygons_bintree(convex, round((double) convex->n_items * 0.5));

        /* find closest (euclidean) distance between convex hull and surface */
        for (i = 0; i < surface->n_points; i++) {
                poly  = find_closest_polygon_point(&surface->points[i], convex, &closest);
                depth[i] = distance_between_points(&surface->points[i], &closest);
        }

        delete_the_bintree(&convex->bintree);

}

void
compute_local_sharpness(polygons_struct *polygons, int n_neighbours[],
                        int *neighbours[], double *sharpness)
{
        int              n, n2, p;
        Point            pts[3];
        Vector           norms[MAX_NEIGHBOURS];
        double           max_radians, r;

        for (p = 0; p < polygons->n_points; p++) {
                if (n_neighbours[p] <= 1)
                        continue; /* skip this point */

                /* Get 2 consecutive neighbors of this node */
                pts[0] = polygons->points[p];
                for (n = 0; n < n_neighbours[p]; n++) {
                        n2 = (n + 1) % n_neighbours[p];

                        /* normal of the triangle */
                        pts[1] = polygons->points[neighbours[p][n]];
                        pts[2] = polygons->points[neighbours[p][n2]];
                        find_polygon_normal(3, pts, &norms[n]);
                }

                max_radians = 0;
                for (n = 0; n < n_neighbours[p]; n++) {
                        n2 = (n + 1) % n_neighbours[p];
                        r = acos(DOT_VECTORS(norms[n], norms[n2]));
                        if (r > max_radians)
                                max_radians = r;
                }
                sharpness[p] = PI2 * 90 * max_radians;
        }
}

void
compute_convexity(polygons_struct *polygons, int n_neighbours[],
               int *neighbours[], double l_convex, double *convexity)
{
        int              p, n;
        double nx, ny, nz, x, y, z, sx, sy, sz, nc;
        Point npt;

        compute_polygon_normals(polygons);

        for (p = 0; p < polygons->n_points; p++) {
            nx = Point_x(polygons->normals[p]);
            ny = Point_y(polygons->normals[p]);
            nz = Point_z(polygons->normals[p]);
            x = Point_x(polygons->points[p]);
            y = Point_y(polygons->points[p]);
            z = Point_z(polygons->points[p]);

            sx = sy = sz = 0.0;
            for (n = 0; n < n_neighbours[p]; n++) {
                npt = polygons->points[neighbours[p][n]];
                sx += Point_x(npt) - x;
                sy += Point_y(npt) - y;
                sz += Point_z(npt) - z;
            }
            if (n > 0) {
                sx /= n; sy /= n; sz /= n;
            }
            convexity[p] = sx*nx + sy*ny + sz*nz;   /* projection onto normal */
        }

}