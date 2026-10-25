#include "minunit.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <bicpl.h>
#include "CAT_Curvature.h"

/*
 * Sulcal depth is the distance to the convex hull, and three defects sat on
 * that path.  compute_sulcus_depth() read the hull's bintree after deleting
 * the hull.  The neighbour lists of the hull code gave every vertex as many
 * slots as triangles but always wrote two, so a vertex that no triangle used
 * wrote into the next vertex's slots, and past the end of the array if it was
 * the last one; an unused vertex 0 was also taken for 1-based numbering, and
 * every index of the caller's mesh was shifted down by one.  A mesh of quads
 * failed the triangle check, whose error was ignored, and the neighbour lists
 * were read uninitialized.
 */

#define SPHERE_R 50.0
#define FOLD_AMP 4.0

/* longitudinal folds: crests where fold() > 0, valleys where it is < 0 */
static double
fold(double x, double y, double z)
{
    return cos(8.0 * atan2(y, x)) * (1.0 - z * z);
}

static double *
make_folded_sphere(polygons_struct *p)
{
    Point centre;
    double *f;
    int i;

    fill_Point(centre, 0.0, 0.0, 0.0);
    create_tetrahedral_sphere(&centre, 1.0, 1.0, 1.0, 20480, p);
    f = (double *)malloc(sizeof(double) * p->n_points);

    for (i = 0; i < p->n_points; i++)
    {
        double x = Point_x(p->points[i]);
        double y = Point_y(p->points[i]);
        double z = Point_z(p->points[i]);
        double len = sqrt(x * x + y * y + z * z);
        double r;

        x /= len;
        y /= len;
        z /= len;
        f[i] = fold(x, y, z);
        r = SPHERE_R + FOLD_AMP * f[i];
        fill_Point(p->points[i], r * x, r * y, r * z);
    }
    compute_polygon_normals(p);
    return f;
}

/**
 * \brief Copy of p with one vertex that no triangle uses, placed far outside.
 *
 * first: the extra vertex becomes vertex 0 and every index moves up by one;
 * otherwise it is appended as the last vertex.
 */
static void
add_unused_vertex(polygons_struct *p, polygons_struct *q, int first)
{
    int i, shift = first ? 1 : 0;

    copy_polygons(p, q);
    q->n_points = p->n_points + 1;
    REALLOC(q->points, q->n_points);
    REALLOC(q->normals, q->n_points);
    if (first)
    {
        for (i = p->n_points; i > 0; i--)
        {
            q->points[i] = q->points[i - 1];
            q->normals[i] = q->normals[i - 1];
        }
    }
    fill_Point(q->points[first ? 0 : p->n_points], 1000.0, 0.0, 0.0);
    fill_Vector(q->normals[first ? 0 : p->n_points], 0.0, 0.0, 0.0);
    for (i = 0; i < NUMBER_INDICES(*q); i++)
        q->indices[i] += shift;
}

static void test_depth_of_folded_sphere(void)
{
    polygons_struct p;
    double *depth, *f;
    double crest = 0.0, valley = 0.0, min_depth = 1e10;
    int i, n, n_crest = 0, n_valley = 0;

    initialize_polygons(&p, WHITE, NULL);
    f = make_folded_sphere(&p);
    n = p.n_points;
    depth = (double *)malloc(sizeof(double) * n);

    MU_ASSERT("sulcal depth of a triangle mesh succeeds",
              compute_sulcus_depth(&p, depth) == 0);

    for (i = 0; i < n; i++)
    {
        min_depth = MIN(min_depth, depth[i]);
        if (f[i] > 0.5)
        {
            crest += depth[i];
            n_crest++;
        }
        else if (f[i] < -0.5)
        {
            valley += depth[i];
            n_valley++;
        }
    }
    crest /= n_crest;
    valley /= n_valley;
    printf("depth: crests %.3f, valleys %.3f\n", crest, valley);
    MU_ASSERT("depth is never negative", min_depth >= 0.0);
    MU_ASSERT("crests lie on the hull", crest < 0.25 * FOLD_AMP);
    /* the hull bridges a valley by a chord, below the crests' radius */
    MU_ASSERT("valleys lie deeper", valley > crest + 0.5 * FOLD_AMP);

    free(depth);
    free(f);
    delete_polygons(&p);
}

static void test_unused_vertices(void)
{
    polygons_struct p, q;
    double *depth, *depth_q;
    int *indices;
    int i, n, first, same, unchanged;

    initialize_polygons(&p, WHITE, NULL);
    free(make_folded_sphere(&p));
    n = p.n_points;
    depth = (double *)malloc(sizeof(double) * n);
    depth_q = (double *)malloc(sizeof(double) * (n + 1));
    compute_sulcus_depth(&p, depth);

    for (first = 0; first <= 1; first++)
    {
        add_unused_vertex(&p, &q, first);
        indices = (int *)malloc(sizeof(int) * NUMBER_INDICES(q));
        memcpy(indices, q.indices, sizeof(int) * NUMBER_INDICES(q));

        MU_ASSERT("a mesh with an unused vertex is handled",
                  compute_sulcus_depth(&q, depth_q) == 0);

        same = 1;
        for (i = 0; i < n; i++)
            if (depth_q[i + first] != depth[i])
                same = 0;
        unchanged = memcmp(indices, q.indices,
                           sizeof(int) * NUMBER_INDICES(q)) == 0;
        if (first)
        {
            MU_ASSERT("an unused vertex 0 changes no depth", same);
            MU_ASSERT("and leaves the caller's indices alone", unchanged);
        }
        else
        {
            MU_ASSERT("an unused last vertex changes no depth", same);
            MU_ASSERT("and leaves the caller's indices alone", unchanged);
        }

        free(indices);
        delete_polygons(&q);
    }

    free(depth);
    free(depth_q);
    delete_polygons(&p);
}

static void test_quads_are_an_error(void)
{
    static const int quads[24] = {
        0, 2, 3, 1,  4, 5, 7, 6,  0, 1, 5, 4,
        2, 6, 7, 3,  0, 4, 6, 2,  1, 3, 7, 5};
    polygons_struct p;
    double depth[8];
    int i, untouched = 1;

    initialize_polygons(&p, WHITE, NULL);
    p.n_points = 8;
    ALLOC(p.points, 8);
    ALLOC(p.normals, 8);
    for (i = 0; i < 8; i++)
    {
        fill_Point(p.points[i], i & 1, (i >> 1) & 1, i >> 2);
        fill_Vector(p.normals[i], 0.0, 0.0, 0.0);
    }
    p.n_items = 6;
    ALLOC(p.end_indices, 6);
    ALLOC(p.indices, 24);
    for (i = 0; i < 6; i++)
        p.end_indices[i] = 4 * (i + 1);
    for (i = 0; i < 24; i++)
        p.indices[i] = quads[i];
    for (i = 0; i < 8; i++)
        depth[i] = 7.0;

    MU_ASSERT("sulcal depth of a quad mesh is an error",
              compute_sulcus_depth(&p, depth) == -1);
    for (i = 0; i < 8; i++)
        if (depth[i] != 7.0)
            untouched = 0;
    MU_ASSERT("and leaves the depth alone", untouched);

    delete_polygons(&p);
}

int main(void)
{
    MU_RUN_TEST(test_depth_of_folded_sphere);
    MU_RUN_TEST(test_unused_vertices);
    MU_RUN_TEST(test_quads_are_an_error);
    printf("%d tests run, %d failed\n", tests_run, tests_failed);
    return tests_failed ? 1 : 0;
}
