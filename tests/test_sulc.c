#include "minunit.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <bicpl.h>
#include "CAT_Curvature.h"

/*
 * curvtype 11 is FreeSurfer's average convexity: smoothing moves sulcal
 * vertices outwards and gyral vertices inwards, so the depth is positive in
 * the valleys of a folded surface and negative on its crests.  It used to be
 * the net displacement of an inflation that also expanded the surface,
 * projected onto the original normals -- dominated by the expansion, which
 * grows with the distance from the centre and so ranked crests above valleys.
 */

#define SPHERE_R 50.0
#define FOLD_AMP 4.0

/* off-centre, so that re-centring the input would show */
static const double offset[3] = {10.0, -5.0, 20.0};

/* longitudinal folds: crests where fold() > 0, valleys where it is < 0 */
static double
fold(double x, double y, double z)
{
    return cos(8.0 * atan2(y, x)) * (1.0 - z * z);
}

static void
make_folded_sphere(polygons_struct *p)
{
    Point centre;
    int i;

    fill_Point(centre, 0.0, 0.0, 0.0);
    create_tetrahedral_sphere(&centre, 1.0, 1.0, 1.0, 81920, p);

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
        r = SPHERE_R + FOLD_AMP * fold(x, y, z);
        fill_Point(p->points[i], r * x + offset[0], r * y + offset[1],
                   r * z + offset[2]);
    }
    compute_polygon_normals(p);
}

static void test_sulc_positive_in_valleys(void)
{
    polygons_struct p;
    Point *before;
    double *depth, *f;
    double mean = 0.0, mf = 0.0, sxy = 0.0, sxx = 0.0, syy = 0.0;
    int i, n, moved = 0;

    make_folded_sphere(&p);
    n = p.n_points;
    depth = (double *)malloc(n * sizeof(double));
    f = (double *)malloc(n * sizeof(double));
    before = (Point *)malloc(n * sizeof(Point));
    for (i = 0; i < n; i++)
    {
        double x = Point_x(p.points[i]) - offset[0];
        double y = Point_y(p.points[i]) - offset[1];
        double z = Point_z(p.points[i]) - offset[2];

        before[i] = p.points[i];
        f[i] = fold(x, y, z / sqrt(x * x + y * y + z * z));
    }

    compute_sulcal_depth_inflation(&p, depth);

    for (i = 0; i < n; i++)
    {
        mean += depth[i];
        mf += f[i];
        if (Point_x(before[i]) != Point_x(p.points[i]) ||
            Point_y(before[i]) != Point_y(p.points[i]) ||
            Point_z(before[i]) != Point_z(p.points[i]))
            moved++;
    }
    mean /= n;
    mf /= n;
    for (i = 0; i < n; i++)
    {
        sxy += (depth[i] - mean) * (f[i] - mf);
        sxx += (depth[i] - mean) * (depth[i] - mean);
        syy += (f[i] - mf) * (f[i] - mf);
    }

    MU_ASSERT("sulcal depth is centred on zero", fabs(mean) < 1e-6);
    MU_ASSERT("sulcal depth is positive in valleys and negative on crests",
              sxy / sqrt(sxx * syy) < -0.9);
    MU_ASSERT("the input surface is not modified", moved == 0);

    free(depth);
    free(f);
    free(before);
    delete_polygons(&p);
}

int main(void)
{
    MU_RUN_TEST(test_sulc_positive_in_valleys);
    printf("%d tests run, %d failed\n", tests_run, tests_failed);
    return tests_failed ? 1 : 0;
}
