#include "minunit.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <bicpl.h>
#include "CAT_SphereUntangle.h"

/*
 * Removal of the folds of a spherical map (CAT_SphereRemoveFolds).
 *
 * The surface is a regular sphere mesh and the map is the same mesh with a
 * radial crease in a cap around the north pole: within the cap the polar
 * angle becomes theta + (cap / pi) sin(2 pi theta / cap), which runs
 * backwards between a third and two thirds of the cap, so that band is
 * mirrored and the map covers it three times -- the flaps that the inflation
 * folds into real spheres.  In the second map the creased cap is also
 * compressed a thousandfold, as the folds of an inflated sphere sit where it
 * is compressed most.  The removal has to leave no bad triangle and move
 * nothing outside the cap; a map without folds has to come back unchanged.
 * A surface with a handle, which no map unfolds, keeps the folds around it
 * without an error.
 */

#define RADIUS 100.0
#define N_TRIANGLES 20480
#define CAP (30.0 * PI / 180.0)

/** \brief Regular sphere mesh of RADIUS around the origin. */
static void
make_sphere(polygons_struct *p)
{
    Point centre;

    fill_Point(centre, 0.0, 0.0, 0.0);
    create_tetrahedral_sphere(&centre, RADIUS, RADIUS, RADIUS, N_TRIANGLES, p);
}

/**
 * \brief Crease the cap within CAP of the north pole; with squeeze < 1 the
 * cap is also compressed by that factor, and the band out to 2 CAP stretched
 * to make up for it, so that the map stays continuous.
 */
static void
crease(polygons_struct *p, double squeeze)
{
    int i;

    for (i = 0; i < p->n_points; i++)
    {
        double x = Point_x(p->points[i]), y = Point_y(p->points[i]);
        double z = Point_z(p->points[i]), r = sqrt(x * x + y * y + z * z);
        double theta = acos(z / r), phi = atan2(y, x), t = theta;

        if (theta < CAP)
            t = squeeze * (theta + CAP / PI * sin(2.0 * PI * theta / CAP));
        else if (theta < 2.0 * CAP)
            t = squeeze * CAP + (theta - CAP) * (2.0 * CAP - squeeze * CAP) / CAP;
        fill_Point(p->points[i], r * sin(t) * cos(phi), r * sin(t) * sin(phi),
                   r * cos(t));
    }
}

/** \brief Signed height of a triangle (2 area / longest edge), positive if it
 *  runs counter-clockwise seen from outside. */
static double
signed_height(const polygons_struct *p, int t)
{
    double a[3], b[3], c[3], u[3], v[3], w[3], m[3], e1, e2, e3, emax, nm;
    int k;

    for (k = 0; k < 3; k++)
    {
        a[k] = Point_coord(p->points[p->indices[POINT_INDEX(p->end_indices, t, 0)]], k);
        b[k] = Point_coord(p->points[p->indices[POINT_INDEX(p->end_indices, t, 1)]], k);
        c[k] = Point_coord(p->points[p->indices[POINT_INDEX(p->end_indices, t, 2)]], k);
        u[k] = b[k] - a[k];
        v[k] = c[k] - a[k];
        m[k] = a[k] + b[k] + c[k];
    }
    w[0] = u[1] * v[2] - u[2] * v[1];
    w[1] = u[2] * v[0] - u[0] * v[2];
    w[2] = u[0] * v[1] - u[1] * v[0];
    nm = sqrt(m[0] * m[0] + m[1] * m[1] + m[2] * m[2]);
    e1 = u[0] * u[0] + u[1] * u[1] + u[2] * u[2];
    e2 = v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
    e3 = (c[0] - b[0]) * (c[0] - b[0]) + (c[1] - b[1]) * (c[1] - b[1]) +
         (c[2] - b[2]) * (c[2] - b[2]);
    emax = fmax(e1, fmax(e2, e3));
    return (w[0] * m[0] + w[1] * m[1] + w[2] * m[2]) / nm / sqrt(emax);
}

/** \brief Triangles below a signed height of 1e-6 RADIUS: folded or degenerate. */
static int
count_bad(const polygons_struct *p)
{
    int t, n = 0;

    for (t = 0; t < p->n_items; t++)
        n += signed_height(p, t) < 1e-6 * RADIUS;
    return n;
}

/**
 * \brief Remove the folds of a creased map; checks that nothing outside the
 * creased area moves and that the moved vertices stay on the sphere.
 */
static void
check_creased(double squeeze, const char *name)
{
    polygons_struct surface, map, before;
    double worst_radius = 0.0;
    int i, bad_before, bad_after, left, moved_outside = 0;
    char msg[256];

    make_sphere(&surface);
    make_sphere(&map);
    crease(&map, squeeze);
    copy_polygons(&map, &before);
    bad_before = count_bad(&map);

    left = CAT_SphereRemoveFolds(&surface, &map, 0);
    bad_after = count_bad(&map);
    for (i = 0; i < map.n_points; i++)
    {
        double x = Point_x(map.points[i]), y = Point_y(map.points[i]);
        double z = Point_z(map.points[i]);
        double theta = acos(Point_z(surface.points[i]) / RADIUS);

        worst_radius = fmax(worst_radius, fabs(sqrt(x * x + y * y + z * z) - RADIUS));
        if (theta > 2.0 * CAP &&
            memcmp(&map.points[i], &before.points[i], sizeof(Point)) != 0)
            moved_outside++;
    }
    if (left != 0 || bad_after != 0)
        printf("  %s: %d bad triangles before, %d after (returned %d)\n", name,
               bad_before, bad_after, left);

    sprintf(msg, "%s: the test starts with folds", name);
    MU_ASSERT(msg, bad_before > 0);
    sprintf(msg, "%s: no bad triangle is left", name);
    MU_ASSERT(msg, left == 0 && bad_after == 0);
    sprintf(msg, "%s: nothing outside the creased area moves", name);
    MU_ASSERT(msg, moved_outside == 0);
    sprintf(msg, "%s: every vertex stays on the sphere", name);
    MU_ASSERT(msg, worst_radius < 1e-4 * RADIUS);

    delete_polygons(&surface);
    delete_polygons(&map);
    delete_polygons(&before);
}

/** \brief A crease at the scale of the mesh is removed. */
static void
test_untangle_crease(void)
{
    check_creased(1.0, "crease");
}

/** \brief A crease in a cap compressed a thousandfold is removed. */
static void
test_untangle_compressed_crease(void)
{
    check_creased(1e-3, "compressed crease");
}

/** \brief A map without folds is returned bit for bit. */
static void
test_untangle_keeps_foldfree_map(void)
{
    polygons_struct surface, map, before;
    int left;

    make_sphere(&surface);
    make_sphere(&map);
    /* compress the cap without creasing it: distorted but fold-free */
    {
        int i;
        for (i = 0; i < map.n_points; i++)
        {
            double x = Point_x(map.points[i]), y = Point_y(map.points[i]);
            double z = Point_z(map.points[i]), r = sqrt(x * x + y * y + z * z);
            double theta = acos(z / r), phi = atan2(y, x);
            double t = PI * pow(theta / PI, 0.6);
            fill_Point(map.points[i], r * sin(t) * cos(phi), r * sin(t) * sin(phi),
                       r * cos(t));
        }
    }
    copy_polygons(&map, &before);
    MU_ASSERT("fold-free map: the test starts without bad triangles", count_bad(&map) == 0);
    left = CAT_SphereRemoveFolds(&surface, &map, 0);
    MU_ASSERT("fold-free map: nothing to remove", left == 0);
    MU_ASSERT("fold-free map: returned bit for bit",
              memcmp(map.points, before.points, sizeof(Point) * map.n_points) == 0);
    delete_polygons(&surface);
    delete_polygons(&map);
    delete_polygons(&before);
}

/**
 * \brief Sphere mesh with a handle: triangles t1 and t2, a few rings apart,
 * are replaced by a prism tube between their vertices (genus 1, Euler
 * characteristic 0).  Returns the first triangle of the tube.
 */
static int
make_sphere_with_handle(polygons_struct *p, int *t2_out)
{
    polygons_struct s;
    double c1[3], best = 1e30;
    int a[3], b[3], t, k, i, t2 = -1, n = 0, *idx;

    make_sphere(&s);
    for (k = 0; k < 3; k++)
    {
        a[k] = s.indices[POINT_INDEX(s.end_indices, 0, k)];
        c1[k] = 0.0;
    }
    for (i = 0; i < 3; i++)
        for (k = 0; k < 3; k++)
            c1[k] += Point_coord(s.points[a[i]], k) / 3.0;
    /* t2: the triangle sharing no vertex with t1 whose centre lies closest to
       15 mm (about four edges) from that of t1 */
    for (t = 1; t < s.n_items; t++)
    {
        double d = 0.0, c2[3] = {0.0, 0.0, 0.0};
        int shared = 0;

        for (i = 0; i < 3; i++)
        {
            int v = s.indices[POINT_INDEX(s.end_indices, t, i)];
            shared += (v == a[0] || v == a[1] || v == a[2]);
            for (k = 0; k < 3; k++)
                c2[k] += Point_coord(s.points[v], k) / 3.0;
        }
        for (k = 0; k < 3; k++)
            d += (c2[k] - c1[k]) * (c2[k] - c1[k]);
        d = fabs(sqrt(d) - 15.0);
        if (!shared && d < best)
        {
            best = d;
            t2 = t;
        }
    }
    for (k = 0; k < 3; k++)
        b[k] = s.indices[POINT_INDEX(s.end_indices, t2, k)];

    /* the quad between edge (a_i, a_i+1) of t1 and edge (b_j+1, b_j) of t2,
       with b paired in reverse (a_i with b_-i), keeps the orientation */
    copy_polygons(&s, p);
    p->n_items = s.n_items - 2 + 6;
    idx = (int *)malloc(sizeof(int) * 3 * p->n_items);
    for (t = 0; t < s.n_items; t++)
        if (t != 0 && t != t2)
            for (k = 0; k < 3; k++)
                idx[n++] = s.indices[POINT_INDEX(s.end_indices, t, k)];
    for (i = 0; i < 3; i++)
    {
        int ai = a[i], aj = a[(i + 1) % 3];
        int bi = b[(3 - i) % 3], bj = b[(3 - (i + 1)) % 3];

        idx[n++] = ai;
        idx[n++] = aj;
        idx[n++] = bj;
        idx[n++] = ai;
        idx[n++] = bj;
        idx[n++] = bi;
    }
    free(p->indices);
    free(p->end_indices);
    p->indices = idx;
    p->end_indices = (int *)malloc(sizeof(int) * p->n_items);
    for (t = 0; t < p->n_items; t++)
        p->end_indices[t] = 3 * (t + 1);
    delete_polygons(&s);
    *t2_out = t2;
    return p->n_items - 6;
}

/**
 * \brief A surface that is not genus 0 keeps the folds around its handle
 * (no map of it is fold-free) but is no error: the patches around the handle
 * are no disks, and the removal goes on with the next cluster or returns.
 */
static void
test_untangle_handle(void)
{
    polygons_struct surface, map;
    int tube, t2, left, bad;

    make_sphere_with_handle(&surface, &t2);
    tube = make_sphere_with_handle(&map, &t2);
    MU_ASSERT("handle: Euler characteristic 0",
              surface.n_points - surface.n_items / 2 == 0);
    bad = count_bad(&map);
    MU_ASSERT("handle: the tube folds on the sphere", bad > 0 && tube > 0);
    left = CAT_SphereRemoveFolds(&surface, &map, 0);
    if (left <= 0)
        printf("  handle: %d bad triangles before, returned %d\n", bad, left);
    MU_ASSERT("handle: no error", left >= 0);
    MU_ASSERT("handle: the folds around the handle stay", left > 0);
    delete_polygons(&surface);
    delete_polygons(&map);
}

/** \brief Meshes that do not belong together are refused. */
static void
test_untangle_rejects_mismatch(void)
{
    polygons_struct surface, map;
    Point centre;

    make_sphere(&surface);
    fill_Point(centre, 0.0, 0.0, 0.0);
    create_tetrahedral_sphere(&centre, RADIUS, RADIUS, RADIUS, N_TRIANGLES / 4, &map);
    MU_ASSERT("different meshes are refused",
              CAT_SphereRemoveFolds(&surface, &map, 0) == -1);
    delete_polygons(&surface);
    delete_polygons(&map);
}

int main(void)
{
    MU_RUN_TEST(test_untangle_crease);
    MU_RUN_TEST(test_untangle_compressed_crease);
    MU_RUN_TEST(test_untangle_keeps_foldfree_map);
    MU_RUN_TEST(test_untangle_handle);
    MU_RUN_TEST(test_untangle_rejects_mismatch);
    printf("%d tests run, %d failed\n", tests_run, tests_failed);
    return tests_failed ? 1 : 0;
}
