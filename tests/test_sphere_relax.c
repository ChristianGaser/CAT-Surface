#include "minunit.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <bicpl.h>
#include "CAT_SphereRelax.h"

/*
 * Metric relaxation of a spherical map (CAT_SphereRelaxMetric).
 *
 * The surface is a regular sphere mesh, so its metric is known exactly, and
 * the map is the same mesh with every vertex pushed towards one pole: area is
 * compressed there and expanded at the other pole, at the scale of the whole
 * sphere, as the inflation compresses the frontal and occipital poles of a
 * brain.  Relaxed, the map has to recover the even area, stay on the sphere
 * and stay fold-free.
 */

#define RADIUS 100.0
#define N_TRIANGLES 20480

/** \brief Regular sphere mesh of RADIUS around the origin. */
static void
make_sphere(polygons_struct *p)
{
    Point centre;

    fill_Point(centre, 0.0, 0.0, 0.0);
    create_tetrahedral_sphere(&centre, RADIUS, RADIUS, RADIUS, N_TRIANGLES, p);
}

/** \brief Push every vertex towards the south pole: theta -> pi (theta/pi)^0.6. */
static void
distort(polygons_struct *p)
{
    int i;

    for (i = 0; i < p->n_points; i++)
    {
        double x = Point_x(p->points[i]), y = Point_y(p->points[i]);
        double z = Point_z(p->points[i]), r = sqrt(x * x + y * y + z * z);
        double theta = acos(z / r), phi = atan2(y, x);
        double t = PI * pow(theta / PI, 0.6);
        fill_Point(p->points[i], r * sin(t) * cos(phi), r * sin(t) * sin(phi),
                   r * cos(t));
    }
}

static double
tri_area(const polygons_struct *p, int t, int *flipped)
{
    Point *a = &p->points[p->indices[POINT_INDEX(p->end_indices, t, 0)]];
    Point *b = &p->points[p->indices[POINT_INDEX(p->end_indices, t, 1)]];
    Point *c = &p->points[p->indices[POINT_INDEX(p->end_indices, t, 2)]];
    double u[3], v[3], w[3], g[3];
    int k;

    for (k = 0; k < 3; k++)
    {
        u[k] = Point_coord(*b, k) - Point_coord(*a, k);
        v[k] = Point_coord(*c, k) - Point_coord(*a, k);
        g[k] = Point_coord(*a, k) + Point_coord(*b, k) + Point_coord(*c, k);
    }
    w[0] = u[1] * v[2] - u[2] * v[1];
    w[1] = u[2] * v[0] - u[0] * v[2];
    w[2] = u[0] * v[1] - u[1] * v[0];
    if (flipped)
        *flipped = (w[0] * g[0] + w[1] * g[1] + w[2] * g[2]) < 0.0;
    return 0.5 * sqrt(w[0] * w[0] + w[1] * w[1] + w[2] * w[2]);
}

/** \brief SD of log2(area share of the map / area share of the surface); counts folds. */
static double
area_distortion(const polygons_struct *surface, const polygons_struct *map,
                int *folds)
{
    double sa = 0.0, sm = 0.0, s1 = 0.0, s2 = 0.0;
    double *as = (double *)malloc(sizeof(double) * surface->n_items);
    double *am = (double *)malloc(sizeof(double) * surface->n_items);
    int t, f;

    *folds = 0;
    for (t = 0; t < surface->n_items; t++)
    {
        as[t] = tri_area(surface, t, NULL);
        am[t] = tri_area(map, t, &f);
        *folds += f;
        sa += as[t];
        sm += am[t];
    }
    for (t = 0; t < surface->n_items; t++)
    {
        double l = log2((am[t] / sm) / (as[t] / sa));
        s1 += l;
        s2 += l * l;
    }
    free(as);
    free(am);
    s1 /= surface->n_items;
    return sqrt(s2 / surface->n_items - s1 * s1);
}

static double
max_radius_error(const polygons_struct *p)
{
    double worst = 0.0;
    int i;

    for (i = 0; i < p->n_points; i++)
    {
        double x = Point_x(p->points[i]), y = Point_y(p->points[i]);
        double z = Point_z(p->points[i]);
        worst = fmax(worst, fabs(sqrt(x * x + y * y + z * z) - RADIUS));
    }
    return worst;
}

/** \brief A map compressed towards one pole recovers the even area of the surface. */
static void
test_relax_recovers_area(void)
{
    polygons_struct surface, map;
    double before, after;
    int folds_before, folds_after;

    make_sphere(&surface);
    make_sphere(&map);
    distort(&map);
    before = area_distortion(&surface, &map, &folds_before);

    MU_ASSERT("relaxation succeeds", CAT_SphereRelaxMetric(&surface, &map, NULL) == 0);
    after = area_distortion(&surface, &map, &folds_after);
    if (after > 0.5 * before)
        printf("  area distortion %.3f -> %.3f\n", before, after);
    MU_ASSERT("distorted map: the test starts fold-free with area distortion",
              folds_before == 0 && before > 0.3);
    MU_ASSERT("distorted map: area distortion at least halved", after < 0.5 * before);
    MU_ASSERT("distorted map: no folded triangle", folds_after == 0);
    MU_ASSERT("distorted map: every vertex stays on the sphere",
              max_radius_error(&map) < 1e-6 * RADIUS);

    delete_polygons(&surface);
    delete_polygons(&map);
}

/** \brief A map that already has the surface's metric is left where it is. */
static void
test_relax_keeps_isometric_map(void)
{
    polygons_struct surface, map;
    double worst = 0.0;
    int i;

    make_sphere(&surface);
    make_sphere(&map);
    MU_ASSERT("relaxation succeeds", CAT_SphereRelaxMetric(&surface, &map, NULL) == 0);
    for (i = 0; i < map.n_points; i++)
    {
        double dx = Point_x(map.points[i]) - Point_x(surface.points[i]);
        double dy = Point_y(map.points[i]) - Point_y(surface.points[i]);
        double dz = Point_z(map.points[i]) - Point_z(surface.points[i]);
        worst = fmax(worst, sqrt(dx * dx + dy * dy + dz * dz));
    }
    if (worst >= 0.05)
        printf("  isometric map: a vertex moved %.4f mm\n", worst);
    MU_ASSERT("isometric map: no vertex moves by 0.05 mm", worst < 0.05);

    delete_polygons(&surface);
    delete_polygons(&map);
}

/** \brief Meshes that do not belong together are refused. */
static void
test_relax_rejects_mismatch(void)
{
    polygons_struct surface, map;
    Point centre;

    make_sphere(&surface);
    fill_Point(centre, 0.0, 0.0, 0.0);
    create_tetrahedral_sphere(&centre, RADIUS, RADIUS, RADIUS, N_TRIANGLES / 4, &map);
    MU_ASSERT("different meshes are refused",
              CAT_SphereRelaxMetric(&surface, &map, NULL) == -1);
    delete_polygons(&surface);
    delete_polygons(&map);
}

int main(void)
{
    MU_RUN_TEST(test_relax_recovers_area);
    MU_RUN_TEST(test_relax_keeps_isometric_map);
    MU_RUN_TEST(test_relax_rejects_mismatch);
    printf("%d tests run, %d failed\n", tests_run, tests_failed);
    return tests_failed ? 1 : 0;
}
