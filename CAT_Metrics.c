/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id: CAT_Metrics.c 218 2011-06-21 12:13:03Z raytrace $
 *
 */

#include <volume_io/internal_volume_io.h>
#include <volume_io/geometry.h>
#include <bicpl.h>
#include <float.h>

#include "CAT_Metrics.h"


void
calc_convexity(polygons_struct *polygons, int n_neighbours[],
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
            if (convexity[p] < 0)
                convexity[p] = 0;

            sx = convexity[p]*nx;              /* move in normal direction */
            sy = convexity[p]*ny;
            sz = convexity[p]*nz;

            //Point_x(depths[p]) += l_convex * sx;
            //Point_y(depths[p]) += l_convex * sy;
            //Point_z(depths[p]) += l_convex * sz;
    }


}

