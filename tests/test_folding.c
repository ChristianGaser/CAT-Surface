#include "minunit.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <bicpl.h>
#include "CAT_Math.h"
#include "CAT_CorrectThicknessFolding.h"

/*
 * The folding correction is restricted to positive mean curvature, which with
 * outward normals is convex, i.e. gyral, cortex.  That has to be the sign of
 * the geometry: the curvature used to be centred first, so a handful of
 * degenerate vertices with |H| ~ 1e4 could move the mean and flip the
 * selection for almost the whole hemisphere.  A sphere is the cleanest case --
 * convex everywhere, or concave everywhere once turned inside out, while the
 * centred curvature splits either of them in half.
 */

#define SPHERE_R 30.0

/* smooth, low-order shape and thickness pattern on the unit sphere */
static double
pattern(double x, double y, double z)
{
    return x * y + 0.5 * z * z;
}

/* sphere with a gentle l = 2 bump: convex everywhere, curvature not constant */
static void
make_bumpy_sphere(polygons_struct *p, int inside_out)
{
    Point centre;
    int i;

    fill_Point(centre, 0.0, 0.0, 0.0);
    create_tetrahedral_sphere(&centre, 1.0, 1.0, 1.0, 20480, p);

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
        r = SPHERE_R * (1.0 + 0.03 * pattern(x, y, z));
        fill_Point(p->points[i], r * x, r * y, r * z);
    }

    /* reversing every triangle turns the normals inward, which makes the
       same shape concave everywhere */
    if (inside_out)
    {
        for (i = 0; i < p->n_items; i++)
        {
            int k = POINT_INDEX(p->end_indices, i, 1);
            int tmp = p->indices[k];
            p->indices[k] = p->indices[k + 1];
            p->indices[k + 1] = tmp;
        }
    }
    compute_polygon_normals(p);
}

static double *
make_thickness(const polygons_struct *p)
{
    double *t = (double *)malloc(sizeof(double) * p->n_points);
    int i;

    for (i = 0; i < p->n_points; i++)
    {
        double x = Point_x(p->points[i]) / SPHERE_R;
        double y = Point_y(p->points[i]) / SPHERE_R;
        double z = Point_z(p->points[i]) / SPHERE_R;
        t[i] = 2.5 + 0.5 * pattern(x, y, z);
    }
    return t;
}

/* fraction of vertices the correction changed */
static double
changed_fraction(int inside_out)
{
    object_struct *obj = create_object(POLYGONS);
    polygons_struct *p = get_polygons_ptr(obj);
    double *t, *t0;
    int i, n_changed = 0;

    make_bumpy_sphere(p, inside_out);
    t = make_thickness(p);
    t0 = make_thickness(p);

    MU_ASSERT("folding correction succeeds",
              CAT_CorrectThicknessFoldingWeighted(p, p->n_points, t, 1.0) == OK);

    for (i = 0; i < p->n_points; i++)
        if (fabs(t[i] - t0[i]) > 1e-9)
            n_changed++;

    free(t);
    free(t0);
    i = p->n_points;
    delete_object(obj);
    return (double)n_changed / (double)i;
}

static void test_convex_surface_corrected(void)
{
    double f = changed_fraction(0);

    /* the centred curvature put about half of the sphere on either side */
    MU_ASSERT("a convex surface is corrected everywhere", f > 0.99);
}

static void test_concave_surface_untouched(void)
{
    double f = changed_fraction(1);

    MU_ASSERT("a concave surface has no convex vertices to correct", f == 0.0);
}

/* pinv of a rank-deficient matrix must satisfy the Moore-Penrose conditions */
static void test_pinv_rank_deficient(void)
{
    const int m = 6, n = 3;
    double **A, **Ainv, **AAi, **AAiA, **AiA, **AiAAi;
    double err1 = 0.0, err2 = 0.0;
    int i, j, rank, round;

    ALLOC2D(A, m, n);
    for (i = 0; i < m; i++)
    {
        A[i][0] = 1.0;
        A[i][1] = (double)i;
        A[i][2] = 2.0 + 3.0 * (double)i; /* = 2*col0 + 3*col1 */
    }

    /* the unused part of S used to be left uninitialized, so poison what the
       allocator may hand back before each attempt */
    for (round = 0; round < 4; round++)
    {
        double **junk;
        ALLOC2D(junk, n, n);
        for (i = 0; i < n; i++)
            for (j = 0; j < n; j++)
                junk[i][j] = 1e30;
        FREE2D(junk);

        ALLOC2D(Ainv, n, m);
        ALLOC2D(AAi, m, m);
        ALLOC2D(AAiA, m, n);
        ALLOC2D(AiA, n, n);
        ALLOC2D(AiAAi, n, m);

        rank = pinv(m, n, A, Ainv);
        MU_ASSERT("rank of the deficient matrix is 2", rank == 2);

        matrix_multiply(m, n, m, A, Ainv, AAi);
        matrix_multiply(m, m, n, AAi, A, AAiA);
        matrix_multiply(n, m, n, Ainv, A, AiA);
        matrix_multiply(n, n, m, AiA, Ainv, AiAAi);

        for (i = 0; i < m; i++)
            for (j = 0; j < n; j++)
            {
                err1 = fmax(err1, fabs(AAiA[i][j] - A[i][j]));
                err2 = fmax(err2, fabs(AiAAi[j][i] - Ainv[j][i]));
            }

        FREE2D(Ainv);
        FREE2D(AAi);
        FREE2D(AAiA);
        FREE2D(AiA);
        FREE2D(AiAAi);
    }
    FREE2D(A);

    MU_ASSERT("A * pinv(A) * A == A", err1 < 1e-8);
    MU_ASSERT("pinv(A) * A * pinv(A) == pinv(A)", err2 < 1e-8);
}

int main(void)
{
    MU_RUN_TEST(test_pinv_rank_deficient);
    MU_RUN_TEST(test_convex_surface_corrected);
    MU_RUN_TEST(test_concave_surface_untouched);
    printf("%d tests run, %d failed\n", tests_run, tests_failed);
    return tests_failed ? 1 : 0;
}
