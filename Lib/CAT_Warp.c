/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 *
 */

#include <bicpl.h>
#include <float.h>

#include "CAT_Warp.h"
#include "CAT_Map.h"
#include "CAT_Surf.h"
#include "CAT_Curvature.h"
#include "CAT_Smooth.h"
#include "CAT_Resample.h"
#include "CAT_SafeAlloc.h"

void
rotate_polygons(polygons_struct *polygons, polygons_struct *rotated_polygons,
        double *rotation_matrix)
{
    int i;
    double x, y, z;
    
    if (rotated_polygons != NULL)
        copy_polygons(polygons, rotated_polygons);
    
    for (i = 0; i < polygons->n_points; i++) {
        x = Point_x(polygons->points[i])*rotation_matrix[0] 
          + Point_y(polygons->points[i])*rotation_matrix[1]
          + Point_z(polygons->points[i])*rotation_matrix[2];
        y = Point_x(polygons->points[i])*rotation_matrix[3] 
          + Point_y(polygons->points[i])*rotation_matrix[4]
          + Point_z(polygons->points[i])*rotation_matrix[5];
        z = Point_x(polygons->points[i])*rotation_matrix[6] 
          + Point_y(polygons->points[i])*rotation_matrix[7]
          + Point_z(polygons->points[i])*rotation_matrix[8];
        if (rotated_polygons != NULL) {
            fill_Point(rotated_polygons->points[i], x, y, z);
        } else  fill_Point(polygons->points[i], x, y, z);
    }
}

void
rotation_to_matrix(double *rotation_matrix, double alpha, double beta,
           double gamma)
{
    int i, j, k;    
    double sum, rot[9];
    
    /* rotation matrices */
    double rot_x[9] = {1.0, 0.0,     0.0,
               0.0, cos(alpha),  sin(alpha),
               0.0, -sin(alpha), cos(alpha)}; 
    double rot_y[9] = {cos(beta),  0.0, sin(beta),
               0.0,      1.0, 0.0,
               -sin(beta), 0.0, cos(beta)}; 
    double rot_z[9] = {cos(gamma),  sin(gamma), 0.0,
               -sin(gamma), cos(gamma), 0.0,
               0.0,     0.0,    1.0}; 

    /* combine x and y rotation */
    for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++) {
            sum = 0.0;
            for (k = 0; k < 3; k++)
                sum += rot_y[i + 3*k] * rot_x[k + 3*j];
            rot[i + 3*j] = sum;
        }
    }

    /* combine with z rotation */
    for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++) {
            sum = 0.0;
            for (k = 0; k < 3; k++) 
                sum += rot_z[i + 3*k] * rot[k + 3*j];
            rotation_matrix[i + 3*j] = sum;
        }
    }
}

void
apply_warp(polygons_struct *polygons, polygons_struct *sphere, double *deform,
       int *dm, int inverse)
{
    Point centre, unit_point, *new_points;
    polygons_struct unit_sphere;
    double xm, ym, xp, yp, x0, x1, y0, y1, weight;
    double *udeform, *vdeform, u, v, *ux, *vy;
    int i, j, p, x, y, m = dm[0]*dm[1];

    if (sphere == NULL) {
        /* create unit sphere w/ same # of triangles as skin surf */
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

    udeform = SAFE_MALLOC(double, m);
    vdeform = SAFE_MALLOC(double, m);

    for (i = 0; i < dm[0]; i++) {
        for (j = 0; j < dm[1]; j++) {
            double theta;
            p = i + dm[0]*j;
            /* Use sin(theta) weighting consistent with equirectangular projection
             * distortion correction. This smoothly tapers from 1 at equator to 0
             * at poles, matching the area distortion of the projection. */
            v = ((double)j + 0.5) / (double)dm[1];
            theta = v * PI;
            weight = sin(theta);

            udeform[p] = (deform[p  ] - (double)i - 1.0) / (double)dm[0];
            vdeform[p] = (deform[p+m] - (double)j - 1.0) / (double)dm[1];
            if (udeform[p] >=  1.0) udeform[p] -= floor(udeform[p]);
            if (udeform[p] <= -1.0) udeform[p] += floor(-udeform[p]);
            if (udeform[p] >=  0.5) udeform[p] -= 1.0;
            if (udeform[p] <= -0.5) udeform[p] += 1.0;
            if (vdeform[p] >=  1.0) vdeform[p] -= floor(vdeform[p]);
            if (vdeform[p] <= -1.0) vdeform[p] += floor(-vdeform[p]);
            udeform[p] *= weight;
            vdeform[p] *= weight;
        }
    }

    ux = SAFE_MALLOC(double, polygons->n_points);
    vy = SAFE_MALLOC(double, polygons->n_points);

    for (p = 0; p < polygons->n_points; p++) {
        map_point_to_unit_sphere(polygons, &polygons->points[p],
                     &unit_sphere, &unit_point);

        if (isnan(Point_x(unit_point)))
            fill_Point(unit_point, Point_x(unit_sphere.points[p]), 
                                   Point_y(unit_sphere.points[p]), 
                                   Point_z(unit_sphere.points[p]));

        point_to_uv(&unit_point, &u, &v);

        xp = u*((double)dm[0]) - 0.5;
        yp = v*((double)dm[1]) - 0.5;

        x = (int) floor(xp); xp -= x; xm = 1.0 - xp;
        y = (int) floor(yp); yp -= y; ym = 1.0 - yp;

        x0 = udeform[bound(x,  y,  dm)];
        x1 = udeform[bound(x+1,y,  dm)];
        y0 = udeform[bound(x,  y+1,dm)];
        y1 = udeform[bound(x+1,y+1,dm)];

        ux[p] = ((xm*x0 + xp*x1)*ym + (xm*y0 + xp*y1)*yp);
        if (ux[p] >=  1.0) ux[p] -= floor(ux[p]);
        if (ux[p] <= -1.0) ux[p] += floor(-ux[p]);
        if (ux[p] <  -0.5) ux[p] += 1.0;
        if (ux[p] >   0.5) ux[p] -= 1.0;

        x0 = vdeform[bound(x,  y,  dm)];
        x1 = vdeform[bound(x+1,y,  dm)];
        y0 = vdeform[bound(x,  y+1,dm)];
        y1 = vdeform[bound(x+1,y+1,dm)];
        vy[p] = ((xm*x0 + xp*x1)*ym + (xm*y0 + xp*y1)*yp);
        if (vy[p] >=  1.0) vy[p] -= floor(vy[p]);
        if (vy[p] <= -1.0) vy[p] += floor(-vy[p]);

        if (inverse) {
            ux[p] = -ux[p];
            vy[p] = -vy[p];
        }

        u += ux[p];
        v += vy[p];

        /* wrap borders */
        if (v < 0.0) {
            v = -v;
            u += 0.5;
        }
        if (v > 1.0) {
            v = 2 - v;
            u += 0.5;
        }
        while (u < 0.0)  u += 1.0;
        while (u >= 1.0) u -= 1.0;

        uv_to_point(u, v, &new_points[p]);
        set_vector_length(&new_points[p], 1.0);
    }
    for (p = 0; p < polygons->n_points; p++)
        polygons->points[p] = new_points[p];

    compute_polygon_normals(polygons);
    free(new_points);
    delete_the_bintree(&polygons->bintree);
    delete_the_bintree(&unit_sphere.bintree);
}

void
apply_uv_warp(polygons_struct *polygons, polygons_struct *sphere, double *ux,
        double *vy, int inverse)
{
    Point centre, unit_point, *new_points, trans_point;
    polygons_struct unit_sphere;
    double u, v, x, y, z;
    double indx, indy;
    double xo, yo, zo;
    int i, p, ind;

    copy_polygons(sphere, &unit_sphere);
    
    create_polygons_bintree(polygons, round((double) polygons->n_items *
                        BINTREE_FACTOR));
    create_polygons_bintree(&unit_sphere,
                round((double) unit_sphere.n_items *
                    BINTREE_FACTOR));

    ALLOC(new_points, sphere->n_points);
  
    for (p = 0; p < polygons->n_points; p++) {
        xo = Point_x(sphere->points[p]);
        yo = Point_y(sphere->points[p]);
        zo = Point_z(sphere->points[p]);

        if (inverse) {
            u = -ux[p];
            v = -vy[p];
        } else {
            u = ux[p];
            v = vy[p];
        }

        x = xo * cos(u) + xo * cos(v)
          + yo * sin(v) + zo * sin(u) - xo;
        y = xo * -sin(u) * sin(u) + yo * cos(u)
          + zo * cos(u) * sin(u) + xo * -sin(v)
          + yo * cos(v) - yo;
        z = xo * -sin(u) * cos(u)
          + yo * -sin(u) + zo * cos(u) * cos(u);
          
        fill_Point(trans_point, x, y, z);

        map_unit_sphere_to_point(&unit_sphere, &trans_point,
                     polygons, &new_points[p]);

    }

    for (p = 0; p < polygons->n_points; p++)
        polygons->points[p] = new_points[p];

    /* set radius to 1 */
    for (i = 0; i < unit_sphere.n_points; i++) 
        set_vector_length(&unit_sphere.points[i], 1.0);

    compute_polygon_normals(polygons);
    free(new_points);
}

void
apply_poly_warp(polygons_struct *polygons, polygons_struct *sphere,
        double *flow, int inverse)
{
    Point centre, unit_point, *new_points, trans_point;
    polygons_struct unit_sphere;
    double u, v, x, y, z, ux, vy;
    double indx, indy;
    int i, p, ind;

    copy_polygons(sphere, &unit_sphere);
    /* set radius to 1 */
    for (i = 0; i < unit_sphere.n_points; i++) 
        set_vector_length(&unit_sphere.points[i], 1.0);

    create_polygons_bintree(polygons, round((double) polygons->n_items *
                        BINTREE_FACTOR));
    create_polygons_bintree(&unit_sphere,
                round((double) unit_sphere.n_items *
                    BINTREE_FACTOR));

    ALLOC(new_points, polygons->n_points);
  
    for (p = 0; p < polygons->n_points; p++) {
        map_point_to_unit_sphere(polygons, &polygons->points[p],
                     &unit_sphere, &unit_point);

        point_to_uv(&unit_point, &u, &v);

        ux = flow[p] - u;
        vy = flow[p + sphere->n_points] - v;

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
        if (v < 0.0)   v = 0.0;
        if (v > 1.0)   v = 1.0;

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
    free(new_points);
}

double compute_cost(double *angles, void *params) {
    OptimizationParams *opt_params = (OptimizationParams *)params;
    polygons_struct rot_src_sphere;
    double rotation_tmp[9];
    double sum_sq = 0.0;
    double d, weight;
    int i;
    double theta;
    Point unit_pt;

    // Rotate source sphere
    rotation_to_matrix(rotation_tmp, angles[0], angles[1], angles[2]);
    rotate_polygons(opt_params->src_sphere, &rot_src_sphere, rotation_tmp);

    // Resample values
    resample_values_sphere(opt_params->trg_sphere, &rot_src_sphere, opt_params->orig_trg, opt_params->map_trg, 0, 0);

    // Compute squared difference
    for (i = 0; i < opt_params->src->n_points; i++) {
        d = opt_params->map_src[i] - opt_params->map_trg[i];
        sum_sq += d * d;
    }

    delete_polygons(&rot_src_sphere);
    return sum_sq;
}

void nelder_mead(double **simplex, double *f_values, int n, int max_iter, double tol, OptimizationParams *params, double *optimal_params, int verbose) {
    int i, j, iter;
    int highest, second_highest, lowest;
    double centroid[n];
    double reflected[n], expanded[n], contracted[n];
    double f_reflected, f_expanded, f_contracted;

    for (iter = 0; iter < max_iter; iter++) {
        // Identify the lowest, highest, and second-highest points
        highest = 0;
        lowest = 0;
        second_highest = 1;
        for (i = 0; i <= n; i++) {
            if (f_values[i] > f_values[highest]) {
                second_highest = highest;
                highest = i;
            } else if (f_values[i] > f_values[second_highest] && i != highest) {
                second_highest = i;
            }
            if (f_values[i] < f_values[lowest]) {
                lowest = i;
            }
        }

        // Compute the centroid
        for (j = 0; j < n; j++) {
            centroid[j] = 0.0;
            for (i = 0; i <= n; i++) {
                if (i != highest) {
                    centroid[j] += simplex[i][j];
                }
            }
            centroid[j] /= n;
        }

        // Reflection
        for (j = 0; j < n; j++) {
            reflected[j] = centroid[j] + ALPHA * (centroid[j] - simplex[highest][j]);
        }
        f_reflected = compute_cost(reflected, params);

        if (f_reflected < f_values[lowest]) {
            // Expansion
            for (j = 0; j < n; j++) {
                expanded[j] = centroid[j] + GAMMA * (reflected[j] - centroid[j]);
            }
            f_expanded = compute_cost(expanded, params);

            if (f_expanded < f_reflected) {
                for (j = 0; j < n; j++) {
                    simplex[highest][j] = expanded[j];
                }
                f_values[highest] = f_expanded;
            } else {
                for (j = 0; j < n; j++) {
                    simplex[highest][j] = reflected[j];
                }
                f_values[highest] = f_reflected;
            }
        } else if (f_reflected < f_values[second_highest]) {
            for (j = 0; j < n; j++) {
                simplex[highest][j] = reflected[j];
            }
            f_values[highest] = f_reflected;
        } else {
            // Contraction
            for (j = 0; j < n; j++) {
                contracted[j] = centroid[j] + RHO * (simplex[highest][j] - centroid[j]);
            }
            f_contracted = compute_cost(contracted, params);

            if (f_contracted < f_values[highest]) {
                for (j = 0; j < n; j++) {
                    simplex[highest][j] = contracted[j];
                }
                f_values[highest] = f_contracted;
            } else {
                // Shrink the simplex
                for (i = 0; i <= n; i++) {
                    if (i != lowest) {
                        for (j = 0; j < n; j++) {
                            simplex[i][j] = simplex[lowest][j] + SIGMA * (simplex[i][j] - simplex[lowest][j]);
                        }
                        f_values[i] = compute_cost(simplex[i], params);
                    }
                }
            }
        }

        // Check for convergence
        double max_diff = 0.0;
        for (i = 0; i <= n; i++) {
            double diff = fabs(f_values[i] - f_values[lowest]);
            if (diff > max_diff) {
                max_diff = diff;
            }
        }
        if (max_diff < tol) {
            break;
        }
    }

    // Copy the optimal parameters from the lowest point
    for (j = 0; j < n; j++) {
        optimal_params[j] = simplex[lowest][j];
    }

    if (verbose) {
        printf("Optimization completed in %d iterations\n", iter);
        printf("Minimum found at:\n");
        for (j = 0; j < n; j++) {
            printf("Angle[%d] = %.6f\n", j, optimal_params[j]);
        }
        printf("Minimum squared difference: %.6f\n", f_values[lowest]);
    }
}

/* input 2 surfaces w/ weighting along x- or z-axes, output weighted average */
void
average_xz_surf(polygons_struct *xsurf, polygons_struct *zsurf,
        polygons_struct *surface)
{
    double xx, xy, xz, zx, zy, zz, phi, *wx, *wz, wtot;
    int p;
    int *n_neighbours, **neighbours;

    copy_polygons(zsurf, surface);
    compute_polygon_normals(xsurf);
    compute_polygon_normals(zsurf);

    wx = SAFE_MALLOC(double, surface->n_points);
    wz = SAFE_MALLOC(double, surface->n_points);

    create_polygon_point_neighbours(surface, TRUE, &n_neighbours,
                    &neighbours, NULL, NULL);

    for (p = 0; p < surface->n_points; p++) {
        phi = acos(Point_x(xsurf->points[p])) / PI;
        wx[p] = exp(-(pow(2.0*phi - 1.0, 2.0)/0.1));
        if (wx[p] <= 0.0) wx[p] = 1e-19;
        
        phi = acos(Point_z(zsurf->points[p])) / PI;
        wz[p] = exp(-(pow(2.0*phi - 1.0, 2.0)/0.1));
        if (wz[p] <= 0.0) wz[p] = 1e-19;
    }

    for (p = 0; p < surface->n_points; p++) {
        xx = Point_x(xsurf->points[p]);
        xy = Point_y(xsurf->points[p]);
        xz = Point_z(xsurf->points[p]);
        zx = Point_x(zsurf->points[p]);
        zy = Point_y(zsurf->points[p]);
        zz = Point_z(zsurf->points[p]);

        wtot = wx[p] + wz[p];
        if (wtot != 0.0) {
            wx[p] /= wtot;
            wz[p] /= wtot;
        }


        fill_Point(surface->points[p], wx[p]*xx + wz[p]*zx,
               wx[p]*xy + wz[p]*zy, wx[p]*xz + wz[p]*zz);
        set_vector_length(&surface->points[p], 1.0);
    }
    free(wx);
    free(wz);
}


/* This function find the optimal rotation parameters to minimize differences in 
   the curvature maps of target and source using the Nelder-Mead (Downhill) approach 
   Optionally applies distortion correction weighting by sin(theta) for better
   2D-to-sphere registration on spherical surfaces */
void
rotate_polygons_to_atlas(polygons_struct *src, polygons_struct *src_sphere,
             polygons_struct *trg, polygons_struct *trg_sphere,
             double fwhm, int curvtype, double *rot, int verbose)
{
    int i, n;
    int n_angles;
    double alpha, beta, gamma, sum_sq, min_sum_sq;
    double degrees, delta, d;
    double curr_alpha = 0.0, curr_beta = 0.0, curr_gamma = 0.0;
    double best_alpha, best_beta, best_gamma;
    double rotation_tmp[9], min_degrees, max_degrees;
    double *orig_trg, *map_trg, *map_src;
    polygons_struct rot_src_sphere;
        
    min_degrees = RADIANS(1.0);
    max_degrees = RADIANS(32.0);
    degrees = max_degrees;
    n_angles  = 4;
    n = 3;
    
    min_sum_sq = 1e15;

    orig_trg = SAFE_MALLOC(double, trg->n_points);
    map_trg  = SAFE_MALLOC(double, src->n_points);
    map_src  = SAFE_MALLOC(double, src->n_points);

    get_smoothed_curvatures(trg, orig_trg, fwhm, curvtype);
    get_smoothed_curvatures(src, map_src, fwhm, curvtype);

    // Initialize optimization parameters
    OptimizationParams params = {src, src_sphere, trg_sphere, orig_trg, map_trg, map_src};

    // Initial simplex
    double simplex[4][3] = {
        {0.0, 0.0, 0.0},
        {0.1, 0.0, 0.0},
        {0.0, 0.1, 0.0},
        {0.0, 0.0, 0.1}
    };
    double *simplex_ptrs[4] = {simplex[0], simplex[1], simplex[2], simplex[3]};
    double f_values[4];

    // Evaluate the cost function at each vertex
    for (i = 0; i <= n; i++) {
        f_values[i] = compute_cost(simplex_ptrs[i], &params);
    }

    // Run optimization
    nelder_mead(simplex_ptrs, f_values, n, MAX_ITER, TOL, &params, rot, verbose);

    // Free memory
    free(orig_trg);
    free(map_trg);
    free(map_src);

}
