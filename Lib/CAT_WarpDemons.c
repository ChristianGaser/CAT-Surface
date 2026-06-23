/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 *
 * Demon-based spherical surface registration. See CAT_WarpDemons.h for the
 * public interface. Method 4 (Spherical Demons, Yeo et al. IEEE TMI 2010)
 * integrates each Gauss-Newton velocity update with a scaling-and-squaring
 * exponential map and composes it onto the running warp, mirroring
 * SD_SphericalExpMap.m / SD_registerAtlas2Sphere.m from the reference toolbox.
 */

#include <bicpl.h>
#include <float.h>
#include <math.h>

#include "CAT_WarpDemons.h"
#include "CAT_Warp.h"
#include "CAT_Surf.h"
#include "CAT_Curvature.h"
#include "CAT_SurfaceIO.h"
#include "CAT_Interpolate.h"
#include "CAT_Resample.h"
#include "CAT_Smooth.h"
#include "CAT_Math.h"
#include "dartel.h"

#define SPHERE_RADIUS 100.0
#define EXP_MAX_COMPOSITIONS 12

/* The running warp is integrated by composing exp(v) with inverse direction,
 * matching the convention of the original additive method-4 update. */
#define EXP_INVERSE 1

/**
 * \brief Free the arrays allocated by init_dartel_poly so it can be re-built.
 *
 * \param dpoly (in/out) dartel helper whose sample-point arrays are released
 * \return void
 */
static void
free_dartel_poly(struct dartel_poly *dpoly)
{
    if (dpoly->u)      { free(dpoly->u);      dpoly->u = NULL; }
    if (dpoly->v)      { free(dpoly->v);      dpoly->v = NULL; }
    if (dpoly->ntheta) { free(dpoly->ntheta); dpoly->ntheta = NULL; }
    if (dpoly->ptheta) { free(dpoly->ptheta); dpoly->ptheta = NULL; }
    if (dpoly->nphi)   { free(dpoly->nphi);   dpoly->nphi = NULL; }
    if (dpoly->pphi)   { free(dpoly->pphi);   dpoly->pphi = NULL; }
}

/**
 * \brief Rescale every sphere vertex to a constant radius.
 *
 * Keeps vertices on the sphere after a warp displaces them off the surface.
 *
 * \param sphere (in/out) spherical mesh whose vertices are reprojected
 * \param radius (in)     target radius
 * \return void
 */
static void
normalize_sphere_radius(polygons_struct *sphere, double radius)
{
    int i;
    double len;

    for (i = 0; i < sphere->n_points; i++) {
        len = sqrt(Point_x(sphere->points[i]) * Point_x(sphere->points[i]) +
                   Point_y(sphere->points[i]) * Point_y(sphere->points[i]) +
                   Point_z(sphere->points[i]) * Point_z(sphere->points[i]));
        if (len > 1e-10) {
            double scale = radius / len;
            Point_x(sphere->points[i]) *= scale;
            Point_y(sphere->points[i]) *= scale;
            Point_z(sphere->points[i]) *= scale;
        }
    }
}

/**
 * \brief Smallest great-circle angle to any neighbour, over all vertices.
 *
 * Used to size the scaling-and-squaring of the exponential map so each
 * sub-step stays well below the inter-vertex spacing (keeps the map
 * diffeomorphic).
 *
 * \param sphere       (in) spherical mesh
 * \param n_neighbours (in) per-vertex neighbour counts
 * \param neighbours   (in) per-vertex neighbour index lists
 * \return minimum neighbour angle in radians (clamped to a small positive value)
 */
static double
min_neighbour_angle(polygons_struct *sphere, int *n_neighbours, int **neighbours)
{
    int i, j;
    double min_angle = PI;

    for (i = 0; i < sphere->n_points; i++) {
        double xi = Point_x(sphere->points[i]);
        double yi = Point_y(sphere->points[i]);
        double zi = Point_z(sphere->points[i]);
        double ri = sqrt(xi*xi + yi*yi + zi*zi);

        if (ri < 1e-10) continue;

        for (j = 0; j < n_neighbours[i]; j++) {
            int nb = neighbours[i][j];
            double xj = Point_x(sphere->points[nb]);
            double yj = Point_y(sphere->points[nb]);
            double zj = Point_z(sphere->points[nb]);
            double rj = sqrt(xj*xj + yj*yj + zj*zj);
            double dot, ang;

            if (rj < 1e-10) continue;
            dot = (xi*xj + yi*yj + zi*zj) / (ri * rj);
            if (dot > 1.0) dot = 1.0;
            if (dot < -1.0) dot = -1.0;
            ang = acos(dot);
            if (ang > 1e-10 && ang < min_angle) min_angle = ang;
        }
    }

    if (min_angle < 1e-6) min_angle = 1e-6;
    return min_angle;
}

/**
 * \brief Central-difference gradient of a scalar field in (theta, phi).
 *
 * Samples the field at the precomputed +/-1 degree neighbour points of each
 * vertex (dartel chart) and forms the symmetric difference.
 *
 * \param polygons (in)  spherical mesh the field lives on
 * \param dpoly    (in)  dartel helper with theta/phi sample points
 * \param f        (in)  double[n_points]; scalar field
 * \param dtheta   (out) double[n_points]; d/dtheta component
 * \param dphi     (out) double[n_points]; d/dphi component
 * \return void
 */
static void
gradient_poly(polygons_struct *polygons, struct dartel_poly *dpoly,
              double f[], double dtheta[], double dphi[])
{
    int i, mm = polygons->n_points;
    double kxm, kxp, kym, kyp;

    for (i = 0; i < mm; i++) {
        kxm = interp_point_unit_sphere(polygons, f, dpoly->ntheta[i]);
        kxp = interp_point_unit_sphere(polygons, f, dpoly->ptheta[i]);
        kym = interp_point_unit_sphere(polygons, f, dpoly->nphi[i]);
        kyp = interp_point_unit_sphere(polygons, f, dpoly->pphi[i]);
        dtheta[i] = (kxp - kxm) / 2.0;
        dphi[i]   = (kyp - kym) / 2.0;
    }
    if (polygons->bintree != NULL) delete_the_bintree(&polygons->bintree);
}

/**
 * \brief Standardize a feature vector (z-score with soft outlier clamping).
 *
 * Subtracts the median, divides by the standard deviation, and softly clamps
 * values beyond +/-3 standard deviations so correlation is not dominated by
 * outliers.
 *
 * \param data   (in/out) double[length]; standardized in place
 * \param length (in)     number of values
 * \return void
 */
static void
normalizeVector(double *data, int length)
{
    double median_data, stdev_data;
    int i;

    median_data = get_median_double(data, length, 0);
    for (i = 0; i < length; i++)
        data[i] -= median_data;

    stdev_data = get_std_double(data, length, 0);
    if (stdev_data < 1e-20) stdev_data = 1.0;

    for (i = 0; i < length; i++)
        data[i] /= stdev_data;

    for (i = 0; i < length; i++) {
        if (data[i] < -3.0)
            data[i] = -3.0 - (1.0 - exp(3.0 - fabs(data[i])));
        if (data[i] > 3.0)
            data[i] = 3.0 + (1.0 - exp(3.0 - fabs(data[i])));
    }
}

/**
 * \brief Pearson correlation coefficient between two vectors.
 *
 * \param x (in) double[N]
 * \param y (in) double[N]
 * \param N (in) number of values
 * \return correlation coefficient in [-1, 1]
 */
static double
Correlation(double *x, double *y, int N)
{
    double EX, EY, EXY, EX2, EY2, denom;
    int i;

    EX = EY = EXY = EX2 = EY2 = 0.0;
    for (i = 0; i < N; i++) {
        EX  += x[i];
        EY  += y[i];
        EXY += x[i] * y[i];
        EX2 += x[i] * x[i];
        EY2 += y[i] * y[i];
    }

    denom = sqrt((N*EX2 - EX*EX) * (N*EY2 - EY*EY));
    if (denom < 1e-20) return 0.0;
    return (N*EXY - EX*EY) / denom;
}

/**
 * \brief Diffeomorphic exponential map of a tangent velocity field.
 *
 * Integrates the velocity field (du, dv), expressed in the dartel theta/phi
 * chart, into a displacement of the reference sphere's vertices using scaling
 * and squaring. The field is scaled by 2^-N, applied once with a first-order
 * step, then composed with itself N times (each composition doubles the
 * integration time). Mirrors SD_SphericalExpMap.m.
 *
 * \param ref_sphere (in)  reference sphere the field is defined on (radius)
 * \param du         (in)  double[n_points]; theta-component velocity (radians)
 * \param dv         (in)  double[n_points]; phi-component velocity (radians)
 * \param radius     (in)  sphere radius used for reprojection
 * \param min_angle  (in)  smallest neighbour angle (radians) for step sizing
 * \param out_points (out) Point[n_points]; integrated vertex positions
 * \return void
 */
static void
spherical_exp_map(polygons_struct *ref_sphere, double *du, double *dv,
                  double radius, double min_angle, Point *out_points)
{
    int i, k, N, n = ref_sphere->n_points;
    double maxnorm = 0.0, scale;
    double *us, *vs;
    polygons_struct def_sphere;

    for (i = 0; i < n; i++) {
        double nrm = sqrt(du[i]*du[i] + dv[i]*dv[i]);
        if (nrm > maxnorm) maxnorm = nrm;
    }

    if (maxnorm < 1e-12 || min_angle < 1e-12) {
        N = 0;
    } else {
        N = (int) ceil(log2(maxnorm / min_angle) + 3.0);
        if (N < 0) N = 0;
        if (N > EXP_MAX_COMPOSITIONS) N = EXP_MAX_COMPOSITIONS;
    }

    us = (double *) malloc(sizeof(double) * n);
    vs = (double *) malloc(sizeof(double) * n);
    scale = pow(2.0, -(double) N);
    for (i = 0; i < n; i++) {
        us[i] = du[i] * scale;
        vs[i] = dv[i] * scale;
    }

    /* first-order step: phi_0 = exp(v / 2^N) applied to the reference grid */
    copy_polygons(ref_sphere, &def_sphere);
    apply_uv_warp(&def_sphere, ref_sphere, us, vs, EXP_INVERSE);
    normalize_sphere_radius(&def_sphere, radius);

    if (N > 0) {
        double *dx, *dy, *dz, *nx, *ny, *nz;
        polygons_struct query_sphere;

        dx = (double *) malloc(sizeof(double) * n);
        dy = (double *) malloc(sizeof(double) * n);
        dz = (double *) malloc(sizeof(double) * n);
        nx = (double *) malloc(sizeof(double) * n);
        ny = (double *) malloc(sizeof(double) * n);
        nz = (double *) malloc(sizeof(double) * n);
        copy_polygons(ref_sphere, &query_sphere);

        /* square N times: phi_{k+1} = phi_k o phi_k */
        for (k = 0; k < N; k++) {
            for (i = 0; i < n; i++) {
                dx[i] = Point_x(def_sphere.points[i]);
                dy[i] = Point_y(def_sphere.points[i]);
                dz[i] = Point_z(def_sphere.points[i]);
                query_sphere.points[i] = def_sphere.points[i];
            }
            /* evaluate the current displacement field at its own image */
            resample_values_sphere(ref_sphere, &query_sphere, dx, nx, 0, 0);
            resample_values_sphere(ref_sphere, &query_sphere, dy, ny, 0, 0);
            resample_values_sphere(ref_sphere, &query_sphere, dz, nz, 0, 0);
            for (i = 0; i < n; i++)
                fill_Point(def_sphere.points[i], nx[i], ny[i], nz[i]);
            normalize_sphere_radius(&def_sphere, radius);
        }

        free(dx); free(dy); free(dz);
        free(nx); free(ny); free(nz);
        delete_polygons(&query_sphere);
    }

    for (i = 0; i < n; i++)
        out_points[i] = def_sphere.points[i];

    delete_polygons(&def_sphere);
    free(us);
    free(vs);
}

/**
 * \brief Run one demon registration stage on a single feature.
 *
 * Computes the chosen curvature feature for both surfaces, then iteratively
 * estimates and applies an update to the source sphere. Method 4 with
 * use_expmap composes a diffeomorphic exponential-map update onto the running
 * warp; all other methods accumulate an additive theta/phi flow.
 *
 * \param src               (in)  source surface (feature + smoothing support)
 * \param src_sphere        (in)  current reference sphere for the source
 * \param orig_sphere       (in)  original source sphere feature values live on
 * \param trg               (in)  template surface
 * \param trg_sphere        (in)  template sphere
 * \param warped_src_sphere (out) deformed source sphere for this stage
 * \param dpoly_src         (in)  dartel helper built on src_sphere
 * \param dpoly_trg         (in)  dartel helper built on trg_sphere
 * \param type              (in)  curvature type for this stage
 * \param opt               (in)  registration options
 * \param fwhm_flow_start   (in)  initial velocity-smoothing FWHM for this stage
 * \return void
 */
static void
warp_demon(polygons_struct *src, polygons_struct *src_sphere,
           polygons_struct *orig_sphere, polygons_struct *trg,
           polygons_struct *trg_sphere, polygons_struct *warped_src_sphere,
           struct dartel_poly *dpoly_src, struct dartel_poly *dpoly_trg,
           int type, const CAT_WarpDemonsOptions *opt, double fwhm_flow_start)
{
    int    *n_neighbours, **neighbours;
    int    i, it, count_break;
    int    diffeo = (opt->method == 4 && opt->use_expmap);
    double idiff, denom_trg, denom_src, sum_diff2;
    double *curv_trg0, *curv_trg, *curv_src0, *curv_src, cc, old_cc, distance;
    double *dtheta_trg, *dphi_trg, *dtheta_src, *dphi_src, *u, *v, *Utheta, *Uphi;
    double fwhm_flow = fwhm_flow_start;
    double min_angle = 0.0;
    double step_factor = opt->step_factor;
    double alpha0 = opt->alpha0;
    double sigma_x = opt->sigma_x;
    /* diffeomorphic-path scratch */
    double *cx = NULL, *cy = NULL, *cz = NULL, *nx = NULL, *ny = NULL, *nz = NULL;
    Point  *inc_points = NULL;
    polygons_struct query_sphere, warped_trg_sphere;
    int n = src->n_points;

    if (src->n_points != trg->n_points) {
        fprintf(stderr, "Source and target have different size!\n");
        return;
    }

    curv_src0  = (double *) malloc(sizeof(double) * n);
    curv_src   = (double *) malloc(sizeof(double) * n);
    curv_trg0  = (double *) malloc(sizeof(double) * n);
    curv_trg   = (double *) malloc(sizeof(double) * n);
    dtheta_trg = (double *) malloc(sizeof(double) * n);
    dphi_trg   = (double *) malloc(sizeof(double) * n);
    dtheta_src = (double *) malloc(sizeof(double) * n);
    dphi_src   = (double *) malloc(sizeof(double) * n);
    u          = (double *) malloc(sizeof(double) * n);
    v          = (double *) malloc(sizeof(double) * n);
    Utheta     = (double *) malloc(sizeof(double) * n);
    Uphi       = (double *) malloc(sizeof(double) * n);

    if (type == 0)
        distance = 3.0;
    else
        distance = 0.0;

    get_all_polygon_point_neighbours(trg, &n_neighbours, &neighbours);
    get_polygon_vertex_curvatures_cg(trg, n_neighbours, neighbours,
                                     distance, type, curv_trg0);
    /* get_all_polygon_point_neighbours allocates a single flat block for
     * neighbours[0]; free it directly (delete_polygon_point_neighbours assumes
     * per-vertex allocations and would corrupt the heap here). */
    free(n_neighbours);
    if (neighbours) { free(neighbours[0]); free(neighbours); }

    get_all_polygon_point_neighbours(src, &n_neighbours, &neighbours);
    get_polygon_vertex_curvatures_cg(src, n_neighbours, neighbours,
                                     distance, type, curv_src0);

    if (diffeo)
        min_angle = min_neighbour_angle(src_sphere, n_neighbours, neighbours);
    free(n_neighbours);
    if (neighbours) { free(neighbours[0]); free(neighbours); }

    normalizeVector(curv_trg0, n);
    normalizeVector(curv_src0, n);

    for (i = 0; i < n; i++) {
        curv_trg[i] = curv_trg0[i];
        u[i] = 0.0;
        v[i] = 0.0;
    }

    /* If src_sphere already carries a warp from a previous stage, pull the
     * source feature through it before starting. */
    if (src_sphere != orig_sphere) {
        resample_values_sphere(orig_sphere, src_sphere, curv_src0, curv_src, 0, 0);
        normalizeVector(curv_src, n);
    } else {
        for (i = 0; i < n; i++)
            curv_src[i] = curv_src0[i];
    }

    if (opt->debug) {
        output_values_any_format("curv.txt", n, curv_src, TYPE_DOUBLE);
        output_values_any_format("curv_trg.txt", n, curv_trg, TYPE_DOUBLE);
    }

    /* gradient of the static (template) feature */
    gradient_poly(trg_sphere, dpoly_trg, curv_trg, dtheta_trg, dphi_trg);

    if (diffeo) {
        cx = (double *) malloc(sizeof(double) * n);
        cy = (double *) malloc(sizeof(double) * n);
        cz = (double *) malloc(sizeof(double) * n);
        nx = (double *) malloc(sizeof(double) * n);
        ny = (double *) malloc(sizeof(double) * n);
        nz = (double *) malloc(sizeof(double) * n);
        inc_points = (Point *) malloc(sizeof(Point) * n);
        copy_polygons(src_sphere, &query_sphere);
    }

    count_break = 0;
    old_cc = -FLT_MAX;

    copy_polygons(src_sphere, warped_src_sphere);
    if (opt->method == 3)
        copy_polygons(trg_sphere, &warped_trg_sphere);

    for (it = 0; it < opt->iters; it++) {

        /* gradient of the moving (source) feature */
        if (opt->method > 1)
            gradient_poly(src_sphere, dpoly_src, curv_src, dtheta_src, dphi_src);

        if (opt->method == 3)
            gradient_poly(trg_sphere, dpoly_trg, curv_trg, dtheta_trg, dphi_trg);

        sum_diff2 = 0.0;
        for (i = 0; i < n; i++) {
            idiff = curv_src[i] - curv_trg[i];
            double idiff2 = idiff * idiff;
            sum_diff2 += idiff2;

            if (opt->method == 3) {
                denom_trg = (dtheta_trg[i]*dtheta_trg[i] + dphi_trg[i]*dphi_trg[i] +
                             dtheta_src[i]*dtheta_src[i] + dphi_src[i]*dphi_src[i]) +
                            alpha0*alpha0*idiff2;
                if (denom_trg == 0.0) {
                    Utheta[i] = 0.0;
                    Uphi[i] = 0.0;
                } else {
                    Utheta[i] = idiff*(dtheta_trg[i]+dtheta_src[i])/denom_trg;
                    Uphi[i]   = idiff*(dphi_trg[i]+dphi_src[i])/denom_trg;
                }
            } else if (opt->method == 4) {
                /* Spherical Demons: per-vertex Gauss-Newton with 2x2 Hessian */
                double g1 = dtheta_trg[i] + dtheta_src[i];
                double g2 = dphi_trg[i] + dphi_src[i];
                double G2 = g1*g1 + g2*g2;

                if (opt->use_hessian && G2 > 1e-9) {
                    double sigma_x_sq = sigma_x * sigma_x;
                    double reg = idiff2 / sigma_x_sq + 0.1;
                    double h11 = g1*g1 + reg;
                    double h12 = g1*g2;
                    double h22 = g2*g2 + reg;
                    double det = h11*h22 - h12*h12;
                    double r1 = idiff * g1;
                    double r2 = idiff * g2;

                    if (fabs(det) > 1e-12) {
                        Utheta[i] = (h22*r1 - h12*r2) / det;
                        Uphi[i]   = (-h12*r1 + h11*r2) / det;
                    } else {
                        double denom = G2 + reg;
                        Utheta[i] = idiff * g1 / denom;
                        Uphi[i]   = idiff * g2 / denom;
                    }
                } else {
                    double sigma_x_sq = sigma_x * sigma_x;
                    double denom = G2 + idiff2/sigma_x_sq + 0.01;
                    Utheta[i] = idiff * g1 / denom;
                    Uphi[i]   = idiff * g2 / denom;
                }
            } else {
                denom_trg = (dtheta_trg[i]*dtheta_trg[i] + dphi_trg[i]*dphi_trg[i]) +
                            alpha0*alpha0*idiff2;
                if (opt->method == 2)
                    denom_src = (dtheta_src[i]*dtheta_src[i] + dphi_src[i]*dphi_src[i]) +
                                alpha0*alpha0*idiff2;
                else
                    denom_src = 1.0;

                if ((denom_trg == 0.0) || (denom_src == 0.0)) {
                    Utheta[i] = 0.0;
                    Uphi[i] = 0.0;
                } else {
                    Utheta[i] = idiff*dtheta_trg[i]/denom_trg;
                    Uphi[i]   = idiff*dphi_trg[i]/denom_trg;
                    if (opt->method == 2) {
                        Utheta[i] += idiff*dtheta_src[i]/denom_src;
                        Uphi[i]   += idiff*dphi_src[i]/denom_src;
                    }
                }
            }
        }

        /* low-pass the velocity update (fluid prior) */
        if (opt->smooth_velocity) {
            smooth_heatkernel(src, Utheta, fwhm_flow);
            smooth_heatkernel(src, Uphi, fwhm_flow);
        }

        /* clamp per-vertex step to limit overshoot/folding */
        if (opt->max_step_deg > 0.0) {
            double max_step = opt->max_step_deg * (PI / 180.0);
            for (i = 0; i < n; i++) {
                double step = sqrt(Utheta[i]*Utheta[i] + Uphi[i]*Uphi[i]);
                if (step > max_step && step > 0.0) {
                    double sc = max_step / step;
                    Utheta[i] *= sc;
                    Uphi[i]   *= sc;
                }
            }
        }

        /* scale the velocity update into the theta/phi chart */
        for (i = 0; i < n; i++) {
            Utheta[i] *= THETA * step_factor;
            Uphi[i]   *= PHI * step_factor;
        }

        if (diffeo) {
            /* integrate exp(v) and compose: warped <- warped o exp(v) */
            spherical_exp_map(src_sphere, Utheta, Uphi, SPHERE_RADIUS,
                              min_angle, inc_points);

            for (i = 0; i < n; i++) {
                cx[i] = Point_x(warped_src_sphere->points[i]);
                cy[i] = Point_y(warped_src_sphere->points[i]);
                cz[i] = Point_z(warped_src_sphere->points[i]);
                query_sphere.points[i] = inc_points[i];
            }
            resample_values_sphere(src_sphere, &query_sphere, cx, nx, 0, 0);
            resample_values_sphere(src_sphere, &query_sphere, cy, ny, 0, 0);
            resample_values_sphere(src_sphere, &query_sphere, cz, nz, 0, 0);
            for (i = 0; i < n; i++)
                fill_Point(warped_src_sphere->points[i], nx[i], ny[i], nz[i]);
            normalize_sphere_radius(warped_src_sphere, SPHERE_RADIUS);

            /* smooth the accumulated displacement field (elastic prior) */
            if (opt->smooth_displacement && opt->fwhm_disp > 0.0) {
                for (i = 0; i < n; i++) {
                    cx[i] = Point_x(warped_src_sphere->points[i]) - Point_x(src_sphere->points[i]);
                    cy[i] = Point_y(warped_src_sphere->points[i]) - Point_y(src_sphere->points[i]);
                    cz[i] = Point_z(warped_src_sphere->points[i]) - Point_z(src_sphere->points[i]);
                }
                smooth_heatkernel(src, cx, opt->fwhm_disp);
                smooth_heatkernel(src, cy, opt->fwhm_disp);
                smooth_heatkernel(src, cz, opt->fwhm_disp);
                for (i = 0; i < n; i++)
                    fill_Point(warped_src_sphere->points[i],
                               Point_x(src_sphere->points[i]) + cx[i],
                               Point_y(src_sphere->points[i]) + cy[i],
                               Point_z(src_sphere->points[i]) + cz[i]);
                normalize_sphere_radius(warped_src_sphere, SPHERE_RADIUS);
            }
        } else {
            /* additive theta/phi flow (methods 1-3 / -no-expmap fallback) */
            for (i = 0; i < n; i++) {
                u[i] += Utheta[i];
                v[i] += Uphi[i];
            }

            if (opt->smooth_displacement && opt->fwhm_disp > 0.0) {
                smooth_heatkernel(src, u, opt->fwhm_disp);
                smooth_heatkernel(src, v, opt->fwhm_disp);
            }

            copy_polygons(src_sphere, warped_src_sphere);
            apply_uv_warp(warped_src_sphere, warped_src_sphere, u, v, 1);
            normalize_sphere_radius(warped_src_sphere, SPHERE_RADIUS);
        }

        /* pull the source feature through the accumulated warp */
        resample_values_sphere(orig_sphere, warped_src_sphere, curv_src0, curv_src, 0, 0);
        normalizeVector(curv_src, n);

        if (opt->method == 3) {
            copy_polygons(trg_sphere, &warped_trg_sphere);
            apply_uv_warp(&warped_trg_sphere, &warped_trg_sphere, u, v, 0);
            normalize_sphere_radius(&warped_trg_sphere, SPHERE_RADIUS);
            resample_values_sphere(trg_sphere, &warped_trg_sphere, curv_trg0, curv_trg, 0, 0);
            normalizeVector(curv_trg, n);
        }

        cc = Correlation(curv_src, curv_trg, n);
        if (opt->verbose)
            printf("%02d: CC=%5.4f diff=%g sigma_x=%g fwhm-flow=%g step=%g\n",
                   it+1, cc, sum_diff2, sigma_x, fwhm_flow, step_factor);

        /* adaptive step size / convergence check */
        if (it >= 3) {
            if (cc < old_cc) {
                count_break++;
                if (opt->use_line_search && step_factor > 0.25) {
                    step_factor *= 0.5;
                    if (opt->verbose)
                        printf("  -> Reducing step factor to %g\n", step_factor);
                }
            } else if (cc/old_cc < 1.00025) {
                count_break++;
            } else {
                count_break = 0;
                if (step_factor < opt->step_factor) step_factor *= 1.1;
            }

            if (count_break > 2) break;
        }

        old_cc = cc;
        if (opt->rate != 0.0 && opt->rate != 1.0)
            fwhm_flow *= opt->rate;
    }

    if (opt->debug) {
        output_values_any_format("u.txt", n, u, TYPE_DOUBLE);
        output_values_any_format("v.txt", n, v, TYPE_DOUBLE);
        output_values_any_format("warped_curv.txt", n, curv_src, TYPE_DOUBLE);
    }

    /* invert deformation for methods that need the inverse transform */
    if (opt->method != 4)
        apply_uv_warp(src_sphere, warped_src_sphere, u, v, 0);

    if (diffeo) {
        delete_polygons(&query_sphere);
        free(cx); free(cy); free(cz);
        free(nx); free(ny); free(nz);
        free(inc_points);
    }
    if (opt->method == 3)
        delete_polygons(&warped_trg_sphere);

    free(curv_src0);
    free(curv_src);
    free(curv_trg0);
    free(curv_trg);
    free(dtheta_trg);
    free(dphi_trg);
    free(dtheta_src);
    free(dphi_src);
    free(u);
    free(v);
    free(Utheta);
    free(Uphi);
}

void
CAT_WarpDemonsDefaults(CAT_WarpDemonsOptions *opt)
{
    opt->n_points            = 20480;
    opt->method              = 4;
    opt->n_steps             = 2;
    opt->curvtype[0]         = 5;   /* sulcal-depth-like */
    opt->curvtype[1]         = 0;   /* mean curvature (averaged over 3mm) */
    opt->curvtype[2]         = 0;
    opt->iters               = 100;
    opt->rotate              = 1;
    opt->smooth_velocity     = 1;
    opt->smooth_displacement = 1;
    opt->use_hessian         = 1;
    opt->use_line_search     = 1;
    opt->use_expmap          = 1;
    opt->fwhm_flow           = 30.0;
    opt->fwhm_curv           = 6.0;
    opt->fwhm_disp           = 5.0;
    opt->rate                = 0.97;
    opt->alpha0              = 0.5;
    opt->max_step_deg        = 10.0;
    opt->sigma_x             = 1.0;
    opt->step_factor         = 1.0;
    opt->verbose             = 0;
    opt->debug               = 0;
}

Status
CAT_WarpDemonsRegister(polygons_struct *src, polygons_struct *src_sphere,
                       polygons_struct *trg, polygons_struct *trg_sphere,
                       polygons_struct *warped_src_sphere,
                       const CAT_WarpDemonsOptions *opt)
{
    polygons_struct *sm_src, *sm_trg, *sm_src_sphere, *sm_trg_sphere;
    polygons_struct *orig_sm_src_sphere;
    struct dartel_poly *dpoly_src, *dpoly_trg;
    double rotation_matrix[9], rot[3];
    int    i, step, n_points = opt->n_points;
    object_struct **objects;
    int dpoly_src_built = 0;

    if (src->n_points != src_sphere->n_points ||
        trg->n_points != trg_sphere->n_points) {
        fprintf(stderr, "Surface and sphere must have matching size.\n");
        return ERROR;
    }

    translate_to_center_of_mass(src_sphere);
    for (i = 0; i < src_sphere->n_points; i++)
        set_vector_length(&src_sphere->points[i], SPHERE_RADIUS);
    translate_to_center_of_mass(trg_sphere);
    for (i = 0; i < trg_sphere->n_points; i++)
        set_vector_length(&trg_sphere->points[i], SPHERE_RADIUS);

    dpoly_src          = (struct dartel_poly *) malloc(sizeof(struct dartel_poly));
    dpoly_trg          = (struct dartel_poly *) malloc(sizeof(struct dartel_poly));
    sm_src             = (polygons_struct *) malloc(sizeof(polygons_struct));
    sm_trg             = (polygons_struct *) malloc(sizeof(polygons_struct));
    sm_src_sphere      = (polygons_struct *) malloc(sizeof(polygons_struct));
    sm_trg_sphere      = (polygons_struct *) malloc(sizeof(polygons_struct));
    orig_sm_src_sphere = (polygons_struct *) malloc(sizeof(polygons_struct));

    resample_spherical_surface(trg, trg_sphere, sm_trg, NULL, NULL, n_points);
    resample_spherical_surface(src_sphere, src_sphere, sm_src_sphere, NULL, NULL, n_points);
    resample_spherical_surface(trg_sphere, trg_sphere, sm_trg_sphere, NULL, NULL, n_points);
    resample_spherical_surface(src, src_sphere, sm_src, NULL, NULL, n_points);
    init_dartel_poly(sm_trg_sphere, dpoly_trg);

    if (opt->verbose)
        printf("Resample surfaces to %d points\n", sm_src->n_points);

    copy_polygons(sm_src_sphere, orig_sm_src_sphere);

    for (step = 0; step < opt->n_steps; step++) {
        int ctype = opt->curvtype[step < CAT_WARP_DEMONS_MAX_STEPS ?
                                  step : CAT_WARP_DEMONS_MAX_STEPS - 1];

        if (step == 0) {
            smooth_heatkernel(sm_src, NULL, opt->fwhm_curv);

            if (opt->rotate) {
                rotate_polygons_to_atlas(sm_src, sm_src_sphere, sm_trg,
                                         sm_trg_sphere, 10.0, 5, rot, opt->verbose);
                rotation_to_matrix(rotation_matrix, rot[0], rot[1], rot[2]);
                rotate_polygons(src_sphere, NULL, rotation_matrix);
                resample_spherical_surface(src_sphere, src_sphere, sm_src_sphere,
                                           NULL, NULL, n_points);
            }
        }

        /* (re)build the source dartel helper on the current reference sphere */
        if (dpoly_src_built)
            free_dartel_poly(dpoly_src);
        init_dartel_poly(sm_src_sphere, dpoly_src);
        dpoly_src_built = 1;

        if (opt->verbose)
            printf("Stage %d/%d: curvature type %d\n", step+1, opt->n_steps, ctype);

        warp_demon(sm_src, sm_src_sphere, orig_sm_src_sphere, sm_trg,
                   sm_trg_sphere, warped_src_sphere, dpoly_src, dpoly_trg,
                   ctype, opt, opt->fwhm_flow);

        /* carry the accumulated warp forward to the next stage */
        if (step < opt->n_steps - 1)
            copy_polygons(warped_src_sphere, sm_src_sphere);
    }
    if (opt->verbose)
        printf("\n");

    /* resample the warp back to the full input sphere resolution */
    objects = resample_surface_to_target_sphere(orig_sm_src_sphere,
                                                 warped_src_sphere, src_sphere,
                                                 NULL, NULL, 0, 0);
    copy_polygons(get_polygons_ptr(objects[0]), warped_src_sphere);
    compute_polygon_normals(warped_src_sphere);

    delete_object_list(1, objects);
    if (dpoly_src_built) free_dartel_poly(dpoly_src);
    free_dartel_poly(dpoly_trg);
    delete_polygons(sm_src);
    delete_polygons(sm_trg);
    delete_polygons(sm_src_sphere);
    delete_polygons(sm_trg_sphere);
    delete_polygons(orig_sm_src_sphere);
    free(sm_src);
    free(sm_trg);
    free(sm_src_sphere);
    free(sm_trg_sphere);
    free(orig_sm_src_sphere);
    free(dpoly_src);
    free(dpoly_trg);

    return OK;
}
