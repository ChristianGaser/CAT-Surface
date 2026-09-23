/* Helpers of the deprecated DARTEL surface registration.
 *
 * These lived in Lib/CAT_Warp.c and Lib/CAT_Map.c until CAT_SurfWarp,
 * CAT_SurfApplyWarp and CAT_SurfApplyWarpValues were deprecated; they had no
 * other caller, so they were moved here with the tools that used them.
 */

#include <bicpl.h>

#include "CAT_Surf.h"
#include "CAT_Map.h"
#include "CAT_Curvature.h"
#include "CAT_SurfaceIO.h"

/**
 * \brief Apply 2D deformation warp to surface using displacement field on sphere.
 *
 * Applies parameterized warping: deformation defined on 2D sphere parameter space
 * (u,v coordinates) is transferred to 3D mesh via barycentric interpolation within
 * spherical triangles. Supports inverse warping (reverse direction). Handles optional
 * internal sphere creation if reference not provided.
 *
 * \param polygons (in/out) source mesh modified by warp
 * \param sphere (in) reference spherical mapping (NULL creates unit sphere)
 * \param deform (in) 2D deformation field (interleaved ux/vy values)
 * \param dm (in) deformation field dimensions [width, height]
 * \param inverse (in) 1 for inverse warp direction; 0 for forward
 */
void apply_warp(polygons_struct *polygons, polygons_struct *sphere, double *deform,
                int *dm, int inverse)
{
    Point centre, unit_point, *new_points;
    polygons_struct unit_sphere;
    double xm, ym, xp, yp, x0, x1, y0, y1, weight;
    double *udeform, *vdeform, u, v, *ux, *vy;
    int i, j, p, x, y, m = dm[0] * dm[1];

    if (sphere == NULL)
    {
        /* create unit sphere w/ same # of triangles as skin surf */
        fill_Point(centre, 0.0, 0.0, 0.0);
        create_tetrahedral_sphere(&centre, 1.0, 1.0, 1.0,
                                  polygons->n_items, &unit_sphere);
    }
    else
    {
        copy_polygons(sphere, &unit_sphere);
        /* set radius to 1 */
        for (i = 0; i < unit_sphere.n_points; i++)
            set_vector_length(&unit_sphere.points[i], 1.0);
    }

    create_polygons_bintree(polygons, round((double)polygons->n_items *
                                            BINTREE_FACTOR));

    create_polygons_bintree(&unit_sphere,
                            round((double)unit_sphere.n_items *
                                  BINTREE_FACTOR));

    ALLOC(new_points, polygons->n_points);

    udeform = SAFE_MALLOC(double, m);
    vdeform = SAFE_MALLOC(double, m);

    for (i = 0; i < dm[0]; i++)
    {
        for (j = 0; j < dm[1]; j++)
        {
            double theta;
            p = i + dm[0] * j;
            /* Use sin(theta) weighting consistent with equirectangular projection
             * distortion correction. This smoothly tapers from 1 at equator to 0
             * at poles, matching the area distortion of the projection. */
            v = ((double)j + 0.5) / (double)dm[1];
            theta = v * PI;
            weight = sin(theta);

            udeform[p] = (deform[p] - (double)i - 1.0) / (double)dm[0];
            vdeform[p] = (deform[p + m] - (double)j - 1.0) / (double)dm[1];
            if (udeform[p] >= 1.0)
                udeform[p] -= floor(udeform[p]);
            if (udeform[p] <= -1.0)
                udeform[p] += floor(-udeform[p]);
            if (udeform[p] >= 0.5)
                udeform[p] -= 1.0;
            if (udeform[p] <= -0.5)
                udeform[p] += 1.0;
            if (vdeform[p] >= 1.0)
                vdeform[p] -= floor(vdeform[p]);
            if (vdeform[p] <= -1.0)
                vdeform[p] += floor(-vdeform[p]);
            udeform[p] *= weight;
            vdeform[p] *= weight;
        }
    }

    ux = SAFE_MALLOC(double, polygons->n_points);
    vy = SAFE_MALLOC(double, polygons->n_points);

    for (p = 0; p < polygons->n_points; p++)
    {
        map_point_to_unit_sphere(polygons, &polygons->points[p],
                                 &unit_sphere, &unit_point);

        if (isnan(Point_x(unit_point)))
            fill_Point(unit_point, Point_x(unit_sphere.points[p]),
                       Point_y(unit_sphere.points[p]),
                       Point_z(unit_sphere.points[p]));

        point_to_uv(&unit_point, &u, &v);

        xp = u * ((double)dm[0]) - 0.5;
        yp = v * ((double)dm[1]) - 0.5;

        x = (int)floor(xp);
        xp -= x;
        xm = 1.0 - xp;
        y = (int)floor(yp);
        yp -= y;
        ym = 1.0 - yp;

        x0 = udeform[bound(x, y, dm)];
        x1 = udeform[bound(x + 1, y, dm)];
        y0 = udeform[bound(x, y + 1, dm)];
        y1 = udeform[bound(x + 1, y + 1, dm)];

        ux[p] = ((xm * x0 + xp * x1) * ym + (xm * y0 + xp * y1) * yp);
        if (ux[p] >= 1.0)
            ux[p] -= floor(ux[p]);
        if (ux[p] <= -1.0)
            ux[p] += floor(-ux[p]);
        if (ux[p] < -0.5)
            ux[p] += 1.0;
        if (ux[p] > 0.5)
            ux[p] -= 1.0;

        x0 = vdeform[bound(x, y, dm)];
        x1 = vdeform[bound(x + 1, y, dm)];
        y0 = vdeform[bound(x, y + 1, dm)];
        y1 = vdeform[bound(x + 1, y + 1, dm)];
        vy[p] = ((xm * x0 + xp * x1) * ym + (xm * y0 + xp * y1) * yp);
        if (vy[p] >= 1.0)
            vy[p] -= floor(vy[p]);
        if (vy[p] <= -1.0)
            vy[p] += floor(-vy[p]);

        if (inverse)
        {
            ux[p] = -ux[p];
            vy[p] = -vy[p];
        }

        u += ux[p];
        v += vy[p];

        /* wrap borders */
        if (v < 0.0)
        {
            v = -v;
            u += 0.5;
        }
        if (v > 1.0)
        {
            v = 2 - v;
            u += 0.5;
        }
        while (u < 0.0)
            u += 1.0;
        while (u >= 1.0)
            u -= 1.0;

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

/* input 2 surfaces w/ weighting along x- or z-axes, output weighted average */
/**
 * \brief Average geometry between two surfaces storing result in second argument.
 *
 * Point-wise averaging: surface = (xsurf + surface) / 2. Used iteratively for
 * surface registration to compute intermediate geometry. Modifies target surface in-place.
 *
 * \param xsurf (in) first surface
 * \param zsurf (in/out) second surface (result stored here)
 * \param surface (in) third surface parameter (unused)
 */
void average_xz_surf(polygons_struct *xsurf, polygons_struct *zsurf,
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

    for (p = 0; p < surface->n_points; p++)
    {
        phi = acos(Point_x(xsurf->points[p])) / PI;
        wx[p] = exp(-(pow(2.0 * phi - 1.0, 2.0) / 0.1));
        if (wx[p] <= 0.0)
            wx[p] = 1e-19;

        phi = acos(Point_z(zsurf->points[p])) / PI;
        wz[p] = exp(-(pow(2.0 * phi - 1.0, 2.0) / 0.1));
        if (wz[p] <= 0.0)
            wz[p] = 1e-19;
    }

    for (p = 0; p < surface->n_points; p++)
    {
        xx = Point_x(xsurf->points[p]);
        xy = Point_y(xsurf->points[p]);
        xz = Point_z(xsurf->points[p]);
        zx = Point_x(zsurf->points[p]);
        zy = Point_y(zsurf->points[p]);
        zz = Point_z(zsurf->points[p]);

        wtot = wx[p] + wz[p];
        if (wtot != 0.0)
        {
            wx[p] /= wtot;
            wz[p] /= wtot;
        }

        fill_Point(surface->points[p], wx[p] * xx + wz[p] * zx,
                   wx[p] * xy + wz[p] * zy, wx[p] * xz + wz[p] * zz);
        set_vector_length(&surface->points[p], 1.0);
    }
    free(wx);
    free(wz);
}

/**
 * \brief Find the sphere rotation that best aligns a source surface with a template.
 *
 * Minimizes the difference between the smoothed curvature maps of source and
 * target over three rotation angles: a coarse multi-start grid (wider along the
 * anterior-posterior axis, where the one-sulcus-off ambiguity is strongest)
 * picks a seed, which the Nelder-Mead simplex method then refines.
 *
 * \param src        (in)  source surface
 * \param src_sphere (in)  its spherical mapping
 * \param trg        (in)  template surface
 * \param trg_sphere (in)  its spherical mapping
 * \param fwhm       (in)  FWHM of the curvature smoothing in mm
 * \param curvtype   (in)  curvature type, as in get_polygon_vertex_curvatures_cg()
 * \param rot        (out) the three rotation angles in radians (see
 *                         rotation_to_matrix())
 * \param verbose    (in)  non-zero to print progress
 */
void rotate_polygons_to_atlas(polygons_struct *src, polygons_struct *src_sphere,
                              polygons_struct *trg, polygons_struct *trg_sphere,
                              double fwhm, int curvtype, double *rot, int verbose)
{
    int i, n;
    double *orig_trg, *map_trg, *map_src;

    n = 3;

    orig_trg = SAFE_MALLOC(double, trg->n_points);
    map_trg = SAFE_MALLOC(double, src->n_points);
    map_src = SAFE_MALLOC(double, src->n_points);

    get_smoothed_curvatures(trg, orig_trg, fwhm, curvtype);
    get_smoothed_curvatures(src, map_src, fwhm, curvtype);

    // Initialize optimization parameters
    OptimizationParams params = {src, src_sphere, trg_sphere, orig_trg, map_trg, map_src, NULL};

    /* Coarse multi-start grid over seed rotations. angles[0] (anterior-posterior)
     * is sampled over a wider range than the other two axes because the
     * one-sulcus-off ambiguity is strongest there. The identity is always among
     * the candidates via best_seed's initial cost. */
    static const double seed0[]  = {-0.6, -0.3, 0.0, 0.3, 0.6}; /* AP axis */
    static const double seed12[] = {-0.3,  0.0, 0.3};           /* other two axes */
    const int n0  = (int)(sizeof(seed0)  / sizeof(seed0[0]));
    const int n12 = (int)(sizeof(seed12) / sizeof(seed12[0]));
    double best_seed[3] = {0.0, 0.0, 0.0};
    double best_cost;
    int a, b, c;

    best_cost = compute_cost(best_seed, &params);
    for (a = 0; a < n0; a++)
        for (b = 0; b < n12; b++)
            for (c = 0; c < n12; c++)
            {
                double seed[3] = {seed0[a], seed12[b], seed12[c]};
                double cost = compute_cost(seed, &params);
                if (cost < best_cost)
                {
                    best_cost = cost;
                    best_seed[0] = seed[0];
                    best_seed[1] = seed[1];
                    best_seed[2] = seed[2];
                }
            }

    if (verbose)
        fprintf(stdout, "Rotation seed: %.3f %.3f %.3f (cost %.4g)\n",
                best_seed[0], best_seed[1], best_seed[2], best_cost);

    /* Refine from the best seed: build the initial simplex around it. */
    double simplex[4][3] = {
        {best_seed[0],       best_seed[1],       best_seed[2]},
        {best_seed[0] + 0.1, best_seed[1],       best_seed[2]},
        {best_seed[0],       best_seed[1] + 0.1, best_seed[2]},
        {best_seed[0],       best_seed[1],       best_seed[2] + 0.1}};
    double *simplex_ptrs[4] = {simplex[0], simplex[1], simplex[2], simplex[3]};
    double f_values[4];

    // Evaluate the cost function at each vertex
    for (i = 0; i <= n; i++)
        f_values[i] = compute_cost(simplex_ptrs[i], &params);

    // Run optimization
    nelder_mead(simplex_ptrs, f_values, n, MAX_ITER, TOL, &params, rot, verbose);

    // Free memory
    free(orig_trg);
    free(map_trg);
    free(map_src);
}

/**
 * \brief Upsample a 2D flow field by factor of 2 using bilinear interpolation.
 *
 * Enlarges a 2D flow field (e.g., displacement field, optical flow) from
 * src_dm dimensions to dst_dm dimensions (typically 2x in each direction).
 * Uses bilinear interpolation to estimate smooth flow values at new grid locations.
 * Commonly used in multi-resolution surface deformation algorithms.
 *
 * \param src_flow (in)  source flow field (size src_dm[0]*src_dm[1])
 * \param src_dm   (in)  source dimensions [width, height]
 * \param dst_flow (out) destination flow field (size dst_dm[0]*dst_dm[1])
 * \param dst_dm   (in)  destination dimensions [width, height]
 */
void upsample_flow_field(double *src_flow, int *src_dm, double *dst_flow, int *dst_dm)
{
    int src_m = src_dm[0] * src_dm[1];
    int dst_m = dst_dm[0] * dst_dm[1];
    int i, j, idx;
    double x_src, y_src;
    double fx, fy;
    int x0, y0, x1, y1;
    double dx, dy;
    double f00, f01, f10, f11;
    double scale_x = (double)src_dm[0] / (double)dst_dm[0];
    double scale_y = (double)src_dm[1] / (double)dst_dm[1];

    for (j = 0; j < dst_dm[1]; j++)
    {
        for (i = 0; i < dst_dm[0]; i++)
        {
            idx = i + j * dst_dm[0];

            /* Map destination coords to source coords */
            x_src = (i + 0.5) * scale_x - 0.5;
            y_src = (j + 0.5) * scale_y - 0.5;

            /* Get integer and fractional parts */
            x0 = (int)floor(x_src);
            y0 = (int)floor(y_src);
            dx = x_src - x0;
            dy = y_src - y0;

            /* Handle boundary with wrapping for x and clamping for y */
            x1 = x0 + 1;
            y1 = y0 + 1;

            /* Wrap x (periodic) */
            if (x0 < 0)
                x0 += src_dm[0];
            if (x0 >= src_dm[0])
                x0 -= src_dm[0];
            if (x1 < 0)
                x1 += src_dm[0];
            if (x1 >= src_dm[0])
                x1 -= src_dm[0];

            /* Clamp y (poles) */
            if (y0 < 0)
                y0 = 0;
            if (y0 >= src_dm[1])
                y0 = src_dm[1] - 1;
            if (y1 < 0)
                y1 = 0;
            if (y1 >= src_dm[1])
                y1 = src_dm[1] - 1;

            /* Bilinear interpolation for u-component */
            f00 = src_flow[x0 + y0 * src_dm[0]];
            f10 = src_flow[x1 + y0 * src_dm[0]];
            f01 = src_flow[x0 + y1 * src_dm[0]];
            f11 = src_flow[x1 + y1 * src_dm[0]];

            /* Scale by 2 since grid spacing is halved */
            dst_flow[idx] = 2.0 * ((1 - dx) * (1 - dy) * f00 +
                                   dx * (1 - dy) * f10 +
                                   (1 - dx) * dy * f01 +
                                   dx * dy * f11);

            /* Bilinear interpolation for v-component */
            f00 = src_flow[src_m + x0 + y0 * src_dm[0]];
            f10 = src_flow[src_m + x1 + y0 * src_dm[0]];
            f01 = src_flow[src_m + x0 + y1 * src_dm[0]];
            f11 = src_flow[src_m + x1 + y1 * src_dm[0]];

            dst_flow[dst_m + idx] = 2.0 * ((1 - dx) * (1 - dy) * f00 +
                                           dx * (1 - dy) * f10 +
                                           (1 - dx) * dy * f01 +
                                           dx * dy * f11);
        }
    }
}

/**
 * \brief Downsample a 2D image by factor of 2 using area averaging.
 *
 * Reduces a 2D image from src_dm to dst_dm dimensions (typically 0.5x in each direction)
 * using neighborhood averaging. Each destination pixel is the mean of the 2x2 source
 * region. Useful for creating multi-resolution pyramids or reducing noise.
 *
 * \param src (in)  source image (size src_dm[0]*src_dm[1])
 * \param src_dm (in) source dimensions [width, height]
 * \param dst (out) destination image (size dst_dm[0]*dst_dm[1])
 * \param dst_dm (in) destination dimensions [width, height]
 */
void downsample_image(double *src, int *src_dm, double *dst, int *dst_dm)
{
    int i, j, di, dj;
    double sum;

    for (dj = 0; dj < dst_dm[1]; dj++)
    {
        for (di = 0; di < dst_dm[0]; di++)
        {
            /* Average 2x2 block */
            i = di * 2;
            j = dj * 2;

            sum = src[i + j * src_dm[0]];
            sum += src[(i + 1) % src_dm[0] + j * src_dm[0]];

            if (j + 1 < src_dm[1])
            {
                sum += src[i + (j + 1) * src_dm[0]];
                sum += src[(i + 1) % src_dm[0] + (j + 1) * src_dm[0]];
                sum /= 4.0;
            }
            else
            {
                sum /= 2.0;
            }

            dst[di + dj * dst_dm[0]] = sum;
        }
    }
}
