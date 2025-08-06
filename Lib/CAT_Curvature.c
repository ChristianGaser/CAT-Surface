/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 *
 */

/* Some of the code is used from caret 5.3 (BrainModelSurfaceCurvature.cxx)  */

#include <bicpl.h>
#include <float.h>

#include "CAT_Curvature.h"
#include "CAT_Smooth.h"
#include "CAT_Surf.h"
#include "CAT_ConvexHull.h"
#include "CAT_DepthPotential.h"

#define PI2 0.6366197724  /* 2/pi */

Vector
projectToPlane(Vector projected, Vector basis[2])
{
    Vector xyz;

    Vector_x(xyz) = DOT_VECTORS(projected,basis[0]);
    Vector_y(xyz) = DOT_VECTORS(projected,basis[1]);
    Vector_z(xyz) = 0.0;
    return(xyz);
}


Vector
projectionVector(Vector vector, Vector normal)
{
    Vector xyz;
    const double t2 = DOT_VECTORS(vector, normal);
 
    Vector_x(xyz) = Vector_x(vector) - (t2*Vector_x(normal));
    Vector_y(xyz) = Vector_y(vector) - (t2*Vector_y(normal));
    Vector_z(xyz) = Vector_z(vector) - (t2*Vector_z(normal));
    return(xyz);
}


void
leastSquares_curv(const int num, Vector dc[], Vector dn[], double *k1, double *k2)
{
    int i;
    double sum1 = 0.0, sum2 = 0.0, sum3 = 0.0;
    double wx = 0.0, wy = 0.0, wxy = 0.0;
    double a = 0.0, b = 0.0, c = 0.0;
    double wx2, wy2, wxy2, t1;
    double trC, detC, temp;
    double deltaPlus, deltaMinus;
   
    for (i = 0; i < num; i++) {
        sum1 += (Vector_x(dc[i]) * Vector_x(dn[i]));
        sum2 += ((Vector_x(dc[i]) * Vector_y(dn[i])) +
             (Vector_y(dc[i]) * Vector_x(dn[i])));
        sum3 += (Vector_y(dc[i]) * Vector_y(dn[i]));
        wx   += (Vector_x(dc[i]) * Vector_x(dc[i]));
        wy   += (Vector_y(dc[i]) * Vector_y(dc[i]));
        wxy  += (Vector_x(dc[i]) * Vector_y(dc[i]));
    }
   
    wx2  = wx  * wx;
    wy2  = wy  * wy;
    wxy2 = wxy * wxy;
    t1 = (wx + wy) * (-wxy2 + wx * wy);

    if (t1 > 0.0) {
        a = (sum3 * wxy2 - sum2 * wxy * wy +
           sum1 * (-wxy2 + wx * wy + wy2)) / t1;
        b = (-(sum3 * wx * wxy) + sum2 * wx * wy -
           sum1 * wxy * wy) / t1;
        c = (-(sum2 * wx * wxy) + sum1 * wxy2 +
           sum3 * (wx2 - wxy2 + wx * wy)) / t1;
    }
   
    trC = a + c;
    detC = a * c - b * b;
    temp = trC * trC - 4 * detC;
    *k1 = 0.0;
    *k2 = 0.0;

    if (temp > 0.0) {
        deltaPlus = sqrt(temp);
        deltaMinus = -deltaPlus;
        *k1 = (trC + deltaPlus)  / 2.0;
        *k2 = (trC + deltaMinus) / 2.0;
    }
}


/* ----------------------------- MNI Header -----------------------------------
@NAME   : compute_points_centroid_and_normal_cg
@INPUT    : polygons
        pidx
        n_neighbours
        neighbours
@OUTPUT   : centroid
        normal
        baselen
        curvature
@RETURNS  : 
@DESCRIPTION: Computes the centroid and normal of the neighbours of a vertex,
        as well as a measure of the size of the polygon defined by the
        neighbours, and the relative curvature of the surface at the
        vertex.
@METHOD   : 
@GLOBALS  : 
@CALLS    : 
@CREATED  :     1993    David MacDonald
@MODIFIED : 
---------------------------------------------------------------------------- */

void
compute_points_centroid_and_normal_cg(polygons_struct *polygons,
                    int pidx, int n_neighbours,
                    int neighbours[], Point *centroid,
                    Vector *normal, double *baselen,
                    int curvtype, double *curvparameter)
{
    int    i;
    Point  neigh_pts[MAX_NEIGHBOURS];
    Vector   deltaNormal[MAX_NEIGHBOURS],deltaCoord[MAX_NEIGHBOURS];
    Vector   basis[2], t1, dn[MAX_NEIGHBOURS], dc[MAX_NEIGHBOURS];
    Vector   projected;
    double     k1, k2;

    if (n_neighbours > 2) {
        for (i = 0; i < n_neighbours; i++)
            neigh_pts[i] = polygons->points[neighbours[i]];

        get_points_centroid(n_neighbours, neigh_pts, centroid);

        find_polygon_normal(n_neighbours, neigh_pts, normal);
    
        for (i = 0; i < n_neighbours; i++) {
            SUB_VECTORS(deltaNormal[i],
                  polygons->normals[neighbours[i]],
                  polygons->normals[pidx]);
            SUB_VECTORS(deltaCoord[i],
                  polygons->points[neighbours[i]],
                  polygons->points[pidx]);
        }
        basis[0] = projectionVector(deltaCoord[0],
                    polygons->normals[pidx]);
        NORMALIZE_VECTOR(basis[0], basis[0]);
        fill_Vector(t1, 0.0, 0.0, 0.0);
    
        SUB_VECTORS(t1, t1, basis[0]);
        CROSS_VECTORS(basis[1], t1, polygons->normals[pidx]);
        NORMALIZE_VECTOR(basis[1], basis[1]);

        for (i = 0; i < n_neighbours; i++) {
            projected = projectionVector(deltaNormal[i],
                         polygons->normals[pidx]);
            dn[i] = projectToPlane(projected, basis);
      
            projected = projectionVector(deltaCoord[i],
                         polygons->normals[pidx]);
            dc[i] = projectToPlane(projected, basis);
        }

        leastSquares_curv(n_neighbours, dc, dn, &k1, &k2);

        switch (curvtype) {
        case 1: /* gauss curvature */
            *curvparameter = k1 * k2;
            break;
        case 2: /* curvedness */
            *curvparameter = sqrt((k1*k1 + k2*k2) / 2.0);
            break;
        case 3: /* shape index */
            *curvparameter = -PI2 * atan((k1 + k2) / (k1 - k2));
            break;
        case 4: /* mean curvature */
            *curvparameter = (k1 + k2) / 2.0;
            break;
        case 6: /* bending enegery */
            *curvparameter = (k1*k1 + k2*k2);
            break;
        case 7: /* sharpness */
            *curvparameter = (k1 - k2)*(k1 - k2);
            break;
        case 8: /* folding index */
            *curvparameter = fabs(k1)*(fabs(k1)- fabs(k2));
            break;
        case 9: /* minimum curvature */
            *curvparameter = k2;
            break;
        case 10: /* maximum curvature */
            *curvparameter = k1;
            break;
        }
    } else {
        *centroid = polygons->points[pidx];
        fill_Vector(*normal, 0.0, 0.0, 0.0);
        *baselen = 1.0;
        *curvparameter = 0.0;
    }
}


/* ----------------------------- MNI Header -----------------------------------
@NAME   : get_polygon_vertex_curvatures_cg
@INPUT    : polygons
        smoothing_distance
        low_threshold
@OUTPUT   : curvatures
@RETURNS  : 
@DESCRIPTION: Computes the curvatures at each vertex of the polygon, using
        1 of two methods.  If smoothing distance is zero, computes
        instantaneous curvature in terms of a fractional relative
        curvature.  If non-zero, returns +/- angle in degrees of the
        smoothed curvature.
@METHOD   : 
@GLOBALS  : 
@CALLS    : 
@CREATED  :     1994    David MacDonald
@MODIFIED : 
---------------------------------------------------------------------------- */

void
get_polygon_vertex_curvatures_cg(polygons_struct *polygons, int n_neighbours[],
                 int *neighbours[], double smoothing_distance,
                 int curvtype, double curvatures[])
{
    int        size, pidx, vidx, p;
    double       curvature, *values, baselen, alpha;
    signed char    *point_done;
    Point      centroid;
    Vector       normal;
    float      *distances;
    BOOLEAN      initialized;
    progress_struct  progress;
    polygons_struct  *polygonsIn;

    compute_polygon_normals(polygons);

    ALLOC(point_done, polygons->n_points);

    for (pidx = 0; pidx < polygons->n_points; pidx++)
        point_done[pidx] = FALSE;

    if (smoothing_distance > 0.0) {
        ALLOC(distances, polygons->n_points);
        initialized = FALSE;
    }
    
    /* depth potential */
    if (curvtype > 11) {
        /* looks weird, but we need a way to define the small values alpha << 1 */
        alpha = 1/(double)curvtype;
        values = compute_depth_potential( polygons, alpha);
        for (pidx = 0; pidx < polygons->n_points; pidx++) curvatures[pidx] = values[pidx];
    /* sulcal depth like estimator */
    } else if (curvtype == 5) {
        polygonsIn = (polygons_struct *) malloc(sizeof(polygons_struct));
        copy_polygons(polygons, polygonsIn);

        /* use smoothing with FWHM of 50mm */
        smooth_heatkernel(polygonsIn, NULL, 50.0);

        compute_polygon_normals(polygonsIn);

        for (p = 0; p < polygons->n_points; p++) {
            curvatures[p] = ((Point_x(polygonsIn->points[p]) - Point_x(polygons->points[p]))*Point_x(polygonsIn->normals[p]) +
                     (Point_y(polygonsIn->points[p]) - Point_y(polygons->points[p]))*Point_y(polygonsIn->normals[p]) +
                     (Point_z(polygonsIn->points[p]) - Point_z(polygons->points[p]))*Point_z(polygonsIn->normals[p]));
        }
        free(polygonsIn);
    } else {

        for (p = 0; p < polygons->n_items; p++) {
            size = GET_OBJECT_SIZE(*polygons, p);

            for (vidx = 0; vidx < size; vidx++) {
                pidx = polygons->indices[
                       POINT_INDEX(polygons->end_indices, p, vidx)];

                if (!point_done[pidx]) {
                    point_done[pidx] = TRUE;
                    /* if smoothing_distance is > 0 mean curvature
                     * will be calculated and averaged over
                     * smoothing_distance */
                    if (smoothing_distance <= 0.0) {
                        compute_points_centroid_and_normal_cg(polygons,
                          pidx, n_neighbours[pidx],
                          neighbours[pidx],
                          &centroid, &normal, &baselen,
                          curvtype, &curvature);

                    } else {
                        curvature = get_smooth_surface_curvature(polygons,
                            n_neighbours,
                            neighbours,
                            p, vidx,
                            initialized,
                            distances,
                            smoothing_distance);

                        initialized = TRUE;
                    }
                    curvatures[pidx] = curvature;
                }
            }
        }
    }

    if (smoothing_distance > 0.0)
        FREE(distances);
        
    FREE(point_done);
}


void
get_smoothed_curvatures(polygons_struct *polygons,
            double *values, double fwhm, int curvtype)
{
    double distance, mn, mx;
    int *n_neighbours, **neighbours;
    int i;

    get_all_polygon_point_neighbours(polygons, &n_neighbours, &neighbours);

    if (curvtype == 0)
        distance = 3.0;
    else distance = 0.0;

    get_polygon_vertex_curvatures_cg(polygons, n_neighbours, neighbours,
                     distance, curvtype, values);

    smooth_heatkernel(polygons, values, fwhm);

    /* scale data to uint8 range */
    mn = FLT_MAX; mx = -FLT_MAX;
    for (i = 0; i < polygons->n_points; i++) {
      if (values[i] > mx) mx = values[i];
      if (values[i] < mn) mn = values[i];
    }

    for (i = 0; i < polygons->n_points; i++)
        values[i] = (values[i] - mn)/(mx - mn);

}

void
compute_sulcus_depth(polygons_struct *surface, double *depth)
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

    delete_object_list(1, object);
    delete_the_bintree(&convex->bintree);

}

void
compute_local_sharpness(polygons_struct *polygons, int n_neighbours[],
            int *neighbours[], double *sharpness)
{
    int        n, n2, p;
    Point      pts[3];
    Vector       norms[MAX_NEIGHBOURS];
    double       max_radians, r;

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
         int *neighbours[], double *convexity)
{
    int        p, n;
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
      convexity[p] = sx*nx + sy*ny + sz*nz; /* projection onto normal */
    }

}
