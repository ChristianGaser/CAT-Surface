/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 *
 */

/* Some of the code is used from caret 5.3 (BrainModelSurface.cxx)  */

#include <volume_io/internal_volume_io.h>
#include <bicpl.h>
#include <float.h>

#include "CAT_Surf.h"
#include "CAT_Blur2d.h"

#define _PI 3.14159265358979323846264338327510

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

/* helper functions -- copy Point values to/from double array of length 3 */
void
to_array(Point *p, double *xyz) {
        int i;

        for (i = 0; i < 3; i++) {
                xyz[i] = (double) Point_coord(*p,i);
        }
}
void
from_array(double *xyz, Point *p) {
        int i;

        for (i = 0; i < 3; i++) {
                Point_coord(*p,i) = xyz[i];
        }
}

double *
get_surface_ratio(double r, polygons_struct *polygons)
{
        int      i, j, x, y, z, a, b, c, nan = 0;
        double   *lf, *avol;
        int      size, poly_size;
        double   area, asum;
        char     str[512];
        Point    points[MAX_POINTS_PER_POLYGON];
        
        avol = (double *) calloc(256*256*256, sizeof(double));

        lf = (double *) calloc(polygons->n_points, sizeof(double));
                
        for (i = 0; i < 256*256*256; i++)
                avol[i] = 0.0;

        for (i = 0; i < polygons->n_items; i++) {
                size = get_polygon_points(polygons, i, points);
                area = get_polygon_surface_area(size, points);

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
        
        for (i = 0; i < polygons->n_points; i++) {
                if (i % 100 == 0) {
                        sprintf(str, "%i/%i", i, polygons->n_points);
                        printf("%s",str);
                        for (j = 0; j < strlen(str); j++)
                                printf("\b");
                        fflush(stdout);
                }

                asum = 0;
                for (x = -r; x <= r; x++) {
                        for (y = -r; y <= r; y++) {
                                for (z = -r; z <= r; z++) {
                                        if (x*x + y*y + z*z < r*r) {
                                                *points = polygons->points[i];
                                                a = x + 128 +
                                                    (int) Point_x(*points);
                                                b = y + 128 +
                                                    (int) Point_y(*points);
                                                c = z + 128 +
                                                    (int) Point_z(*points);
                                                if (a >= 0 && a < 256 &&
                                                    b >= 0 && b < 256 &&
                                                    c >= 0 && c < 256)
                                                        asum += avol[c*65536 +
                                                                     b*256 + a];
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
                lf[i] = asum / (_PI*r*r) / 3.0;
                
                if (lf[i] != lf[i])
                        nan++;
        }
        free(avol);
        
        if (nan)
                printf("ERROR: there are %i NaN\n", nan);
        
        return lf;
}

/* assign each point the average of area values for neighboring polygons */
double
get_area_of_points(polygons_struct *polygons, double *area_values)
{
        int             *pcount;
        double          poly_size, area;
        Point           points[MAX_POINTS_PER_POLYGON];
        int             ptidx, poly, vertidx, size;
        double          surface_area = 0.0;

        pcount = (int *) malloc(sizeof(int) * polygons->n_points);
        memset(pcount, 0, sizeof(int) * polygons->n_points);
        memset(area_values, 0.0, sizeof(double) * polygons->n_points);

        for (poly = 0; poly < polygons->n_items; poly++) {    
                size = get_polygon_points(polygons, poly, points);
                area = get_polygon_surface_area(size, points);
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

        for (ptidx = 0; ptidx < polygons->n_points; ptidx++) {
                if (pcount[ptidx] > 0)
                        area_values[ptidx] /= pcount[ptidx];
        }

        FREE(pcount);
        return(surface_area);
}

double
get_area_of_polygons(polygons_struct *polygons, double *area_values)
{
        int             poly, size;
        double          surface_area = 0.0;
        Point           points[MAX_POINTS_PER_POLYGON];

        for (poly = 0; poly < polygons->n_items; poly++) {    
                size = get_polygon_points(polygons, poly, points);
                area_values[poly] = get_polygon_surface_area(size, points);
                surface_area += area_values[poly];
        }
        return(surface_area);
}


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

double
get_largest_dist(polygons_struct *polygons)
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


void
set_vector_length(Point *p, double newLength)
{
        int   j;
        const double len = sqrt(Point_x(*p)*Point_x(*p) + 
                                Point_y(*p)*Point_y(*p) +
                                Point_z(*p)*Point_z(*p));
        if (len > 0) {
                const double scale = newLength / len;
                for (j = 0; j < 3; j++)
                        Point_coord(*p, j) *= scale;
        }
}

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

int
count_edges(polygons_struct *polygons, int n_neighbours[], int *neighbours[])
{
        int    p, n, nn, n_edges, n_duplicate_edges;

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

        if (n_duplicate_edges > 0)
                printf("N duplicate edges: %d\n", n_duplicate_edges);

        return(n_edges);
}

void
apply_warp(polygons_struct *polygons, polygons_struct *sphere, double *flow, int *size_map, int inverse)
{
        Point             centre, unit_point, *new_points, trans_point;
        polygons_struct   unit_sphere;
        double            inflow_x, inflow_y, u, v, x, y, z, ux, vy;
        double            indx, indy;
        int               i, p, ind;

        if (sphere == NULL) {
                /* create unit sphere with same number of triangles as skin surface */
                fill_Point(centre, 0.0, 0.0, 0.0);
                create_tetrahedral_sphere(&centre, 1.0, 1.0, 1.0,
                                  polygons->n_items, &unit_sphere);
        } else {
                copy_polygons(sphere, &unit_sphere);
                /* set radius to 1 */
                for (i = 0; i < unit_sphere.n_points; i++) 
                        set_vector_length(&unit_sphere.points[i], 1.0);
        }

        create_polygons_bintree(polygons, round((double) polygons->n_items *
                                                BINTREE_FACTOR));

        create_polygons_bintree(&unit_sphere,
                                round((double) unit_sphere.n_items *
                                      BINTREE_FACTOR));

        ALLOC(new_points, polygons->n_points);
  
        inflow_x = (double)size_map[0] - 1.0;
        inflow_y = (double)size_map[1] - 1.0;

        for (p = 0; p < polygons->n_points; p++) {
                map_point_to_unit_sphere(polygons, &polygons->points[p],
                                         &unit_sphere, &unit_point);

                point_to_uv(&unit_point, &u, &v);

                indx = u*inflow_x;
                indy = v*inflow_y;
                ind  = (int) round(indx) + (size_map[0] * (int) round(indy));

                ux = (flow[ind] - 1.0 - indx) / inflow_x;
                vy = (flow[ind + size_map[0]*size_map[1]] - 1.0 - indy) / inflow_y;
                      
                if (inverse) {
                        u -= ux;
                        v -= vy;
                } else {
                        u += ux;
                        v += vy;
                }

                /* wrap borders */
                while (u < 0.0)  u += 1.0;
                while (u >= 1.0) u -= 1.0;
                if (v < 0.0)     v = 0.0;
                if (v > 1.0)     v = 1.0;

                uv_to_point(u, v, &unit_point);

                x = Point_x(unit_point);
                y = Point_y(unit_point);
                z = Point_z(unit_point);

                fill_Point(trans_point, x, y, z);

                map_unit_sphere_to_point(&unit_sphere, &trans_point,
                                         polygons, &new_points[p]);

        }

        for (p = 0; p < polygons->n_points; p++)
                polygons->points[p] = new_points[p];

        compute_polygon_normals(polygons);
}

int
euler_characteristic(polygons_struct *polygons)
{
        int n_edges, *n_neighbours, **neighbours;

        create_polygon_point_neighbours(polygons, TRUE, &n_neighbours,
                                        &neighbours, NULL, NULL);
        n_edges = count_edges(polygons, n_neighbours, neighbours);
    
        delete_polygon_point_neighbours(polygons, n_neighbours,
                                        neighbours, NULL, NULL);

        create_polygon_point_neighbours(polygons, FALSE, &n_neighbours,
                                        &neighbours, NULL, NULL);

        n_edges = count_edges(polygons, n_neighbours, neighbours);

        return(polygons->n_items + polygons->n_points - n_edges);
}

void
convert_ellipsoid_to_sphere_with_surface_area(polygons_struct *polygons,
                                              double desiredSurfaceArea)
{
        int     i;
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
                 *  x*x   y*y   z*z
                 *  --- + --- + --- = 1.0
                 *  A*A   B*B   C*C
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

void
linear_smoothing(polygons_struct *polygons, double strength, int iters,
                 int smoothEdgesEveryXIters, int *smoothOnlyTheseNodes,
                 int projectToSphereEveryXIters)
{
        int     i, j, k, l;
        int     *n_neighbours, **neighbours;
        int     pidx;
        double  xyz[3], pt[3];

        BOOLEAN smoothSubsetOfNodes = 0;
        const double invstr = 1.0 - strength;
        double radius = get_largest_dist(polygons);
    
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
}

void
areal_smoothing(polygons_struct *polygons, double strength, int iters,
                int smoothEdgesEveryXIters, int *smoothOnlyTheseNodes,
                int projectToSphereEveryXIters)
{
        int     i, j, k, l;
        int     n1, n2, next;
        int     *n_neighbours, **neighbours;
        double  *area_values;
        Point   pts[1000];
        double tileAreas[32], tileCenters[32*3];
        double xyz[3], pt1[3], pt2[3], pt3[3];
        double totalArea, weight;

        BOOLEAN smoothSubsetOfNodes = 0;
        BOOLEAN smoothEdges, smoothIt;
        double invstr = 1.0 - strength;
        double radius = get_largest_dist(polygons);
    
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
        
                // If the surface should be projected to a sphere
                if (projectToSphereEveryXIters > 0) {
                        if ((k % projectToSphereEveryXIters) == 0) {
                                for (i = 0; i < polygons->n_points; i++) {
                                        set_vector_length(&polygons->points[i],
                                                          radius);
                                }
                        }
                }
        }
}

        
void
distance_smoothing(polygons_struct *polygons, double strength, int iters,
                   int smoothEdgesEveryXIters, int *smoothOnlyTheseNodes,
                   int projectToSphereEveryXIters)
{
        int     i, j, k, l, pidx;
        int     *n_neighbours, **neighbours;
        double  totalDistance, tileDist[MAX_POINTS_PER_POLYGON];
        double  pt[3], xyz[3];
        double  weight, invstr, radius, div;
        BOOLEAN smoothSubsetOfNodes = 0;

        invstr = 1.0 - strength;
        radius = get_largest_dist(polygons);
    
        create_polygon_point_neighbours(polygons, TRUE, &n_neighbours,
                                        &neighbours, NULL, NULL);

        if (smoothOnlyTheseNodes != NULL)
                 smoothSubsetOfNodes = 1;

        for (k = 1; k < iters; k++) {
                /* should the edges be smoothed? */
                BOOLEAN smoothEdges = 0;
                if (smoothEdgesEveryXIters > 0) {
                        if ((k % smoothEdgesEveryXIters) == 0)
                                smoothEdges = 1;
                }

                for (i = 0; i < polygons->n_points; i++) {
                        BOOLEAN smoothIt = smoothEdges;
                        if (smoothIt && smoothSubsetOfNodes)
                                smoothIt = (smoothOnlyTheseNodes)[i];
        
                        if (!smoothIt || n_neighbours[i] <= 1)
                                continue; /* skip this point */

                        /* use distance weighting only for a limited number */
                        /* of neighbours (usually tetras have 5-7 neighbours) */
                        if (n_neighbours[i] < 20) { 
                                totalDistance = 0.0;    
                                for (j = 0; j < n_neighbours[i]; j++) {  
                                        tileDist[j] = distance_between_points(
                                           &polygons->points[i],
                                           &polygons->points[neighbours[i][j]]);
                                        totalDistance += tileDist[j];                      
                                }

                                /* Compute the influence of neighboring nodes */
                                for (j = 0; j < 3; j++)
                                        xyz[j] = 0.0;
                                for (j = 0; j < n_neighbours[i]; j++) {
                                        if (tileDist[j] <= 0.0)
                                                continue;

                                        pidx = neighbours[i][j];
                                        to_array(&polygons->points[pidx], pt);
                                        weight = tileDist[j] / totalDistance;
                                        for (l = 0; l < 3; l++)
                                                xyz[l] += weight * pt[l];
                                }
                                div = 1;

                        /* use linear smoothing if number of */
                        /* neighbours is too large */
                        } else {
                                for (j = 0; j < 3; j++)
                                        xyz[j] = 0.0;
                                for (j = 0; j < n_neighbours[i]; j++) {
                                        pidx = neighbours[i][j];
                                        to_array(&polygons->points[pidx], pt);
                                        for (l = 0; l < 3; l++)
                                                xyz[l] += pt[l];
                                }
                                div = (double) n_neighbours[i];
                        }
                        /* Update the nodes position */
                        to_array(&polygons->points[i], pt);
                        for (l = 0; l < 3; l++) {
                                pt[l] = (pt[l] * invstr) +
                                        (xyz[l] / div * strength);
                        }
                        from_array(pt, &polygons->points[i]);
                }
            
                // If the surface should be projected to a sphere
                if (projectToSphereEveryXIters > 0) {
                        if ((k % projectToSphereEveryXIters) == 0) {
                                for (i = 0; i < polygons->n_points; i++) {
                                        set_vector_length(&polygons->points[i],
                                                          radius);
                                }
                        }
                }
        }
}


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

        polygons_struct     *polygons;
        int                 i, j, cycle, n;
        const double        inflationFactor = inflationFactorIn - 1.0;
        double              *avgCompStretch, *compStretch;
        double              *maxLinDistort, *avgArealComp;
        double              *stretching, *area_values, *area_valuesIn;
        int                 *n_neighbours, **neighbours, *needSmoothing, nidx;
        object_struct       *out_object;
        double              bounds[6], xyz[3], nodept[3], nodeptIn[3];
        double              xdiff, ydiff, zdiff;
        double              dx, dy, dz, dist, distIn, ratio;
        double              x, y, z, r, k;
        double              SA, SA_ratio, inflatedSA;
        double              numNeighbors, neighpt[3], neighptIn[3];
        double              tileArea, tileAreaIn;
        int                 numDistortionAboveThresh;
        double              maxDistort, minDistort, distort;
    
        /* Copy the fiducial surface since it will be modified */
        /* (translated to center of mass) */
        out_object = create_object(POLYGONS);
        polygons = get_polygons_ptr(out_object);
        copy_polygons(polygonsIn, polygons);
    
        /* Translate the fiducial to center of mass */
        translate_to_center_of_mass(polygons);
        translate_to_center_of_mass(polygonsIn);

        /* Get bounds of fiducial surface */
        get_bounds(polygons, bounds);

        xdiff = bounds[1] - bounds[0];
        ydiff = bounds[3] - bounds[2];
        zdiff = bounds[5] - bounds[4];

        SA = get_polygons_surface_area(polygons);

        ALLOC(avgCompStretch, polygons->n_points);
        ALLOC(maxLinDistort, polygons->n_points);
        ALLOC(avgArealComp, polygons->n_points);
        ALLOC(compStretch, polygons->n_points);
        ALLOC(stretching, polygons->n_points);

        SA_ratio = 0.0;

        for (cycle = 0; cycle < (n_smoothingCycles + 1); cycle++) {
                if (cycle < n_smoothingCycles) {
                        /* Step 6a: Apply Smoothing to AUX coord */
                        areal_smoothing(polygonsIn, regSmoothStrength, 
                                        regSmoothIters, 1, NULL, 0);

                        /* Step 6b: Incrementally inflate AUX surface by */
                        /*          Ellipsoidal Projection  */
                        for (i = 0; i < polygons->n_points; i++) {
                                to_array(&polygonsIn->points[i], xyz);
                
                                x = xyz[0] / xdiff;
                                y = xyz[1] / ydiff;
                                z = xyz[2] / zdiff;

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

                ALLOC(area_values, polygons->n_points);
                get_area_of_points(polygons, area_values);
                ALLOC( area_valuesIn, polygonsIn->n_points);
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
                                        ratio = distIn / dist;
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

                                avgArealComp[i] += distort;
                                /* arealComp[i] += distort; */
                                /* optionally,  log(distort) / log2; */

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
                                avgCompStretch[i] /= (double)
                                                      (n_neighbours[i] + 1);
                        }
                }

                /* Step 6e: Flag highly compressed/stretched nodes for */
                /* targeted smoothing */
                numDistortionAboveThresh = 0;
                maxDistort = -FLT_MAX;
                minDistort =  FLT_MAX;
                ALLOC(needSmoothing, polygons->n_points);
      
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
                        areal_smoothing(polygonsIn,
                                        fingerSmoothStrength,
                                        fingerSmoothIters,
                                        1, needSmoothing, 0);
                }
        }
      
        compute_polygon_normals(polygonsIn);

        FREE(area_values);
        FREE(area_valuesIn);
        FREE(avgCompStretch);
        FREE(maxLinDistort);
        FREE(avgArealComp);
        FREE(compStretch);
        FREE(stretching);
        FREE(needSmoothing);
}

