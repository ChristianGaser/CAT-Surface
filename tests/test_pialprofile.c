#include "minunit.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "CAT_SurfPialProfile.h"

#include "test_phantoms.h"

/* GM slab |x| < 1.4 between two CSF regions */
static double
slab_phantom(double x, double y, double z)
{
    (void)y;
    (void)z;
    return 1.0 + ramp(1.4 - fabs(x));
}

/* ------------------------------------------------------------------ */

static void test_reaches_boundary_from_far_inside(void)
{
    int n[3] = {64, 64, 64};
    nifti_image *nim = make_nim(n, -1);
    float *vol = fill_volume(nim, sphere_phantom);
    object_struct *obj = create_object(POLYGONS);
    polygons_struct *p = get_polygons_ptr(obj);
    CAT_PialProfileOptions opts;
    double mean, min, max;

    /* 1.5 mm inside the boundary: needs more than a 1 mm search */
    make_sphere(p, 7.5);
    CAT_PialProfileOptionsInit(&opts);
    MU_ASSERT("profile placement succeeds",
              CAT_SurfDeformPialProfile(p, vol, nim, &opts) == 0);
    radius_stats(p, &mean, &min, &max);
    MU_ASSERT("sphere reaches the 1.5 boundary at r = 9", fabs(mean - 9.0) < 0.1);
    MU_ASSERT("no vertex left far from the boundary", min > 8.75 && max < 9.25);

    delete_object(obj);
    free(vol);
    nifti_image_free(nim);
}

static void test_holds_without_boundary_in_range(void)
{
    int n[3] = {64, 64, 64};
    nifti_image *nim = make_nim(n, -1);
    float *vol = fill_volume(nim, sphere_phantom);
    object_struct *obj = create_object(POLYGONS);
    polygons_struct *p = get_polygons_ptr(obj);
    CAT_PialProfileOptions opts;
    double mean, min, max;

    /* boundary 1.5 mm away but only 1 mm searched: grey matter plateau */
    make_sphere(p, 7.5);
    CAT_PialProfileOptionsInit(&opts);
    opts.search_out = 1.0;
    CAT_SurfDeformPialProfile(p, vol, nim, &opts);
    radius_stats(p, &mean, &min, &max);
    MU_ASSERT("plateau vertices are held, not ballooned", fabs(mean - 7.5) < 0.1);

    delete_object(obj);
    free(vol);
    nifti_image_free(nim);
}

static void test_stops_at_valley_bottom(void)
{
    int n[3] = {64, 64, 64};
    nifti_image *nim = make_nim(n, -1);
    float *vol = fill_volume(nim, valley_phantom);
    object_struct *obj = create_object(POLYGONS);
    polygons_struct *p = get_polygons_ptr(obj);
    CAT_PialProfileOptions opts;
    double mean, min, max;

    make_sphere(p, 8.0);
    CAT_PialProfileOptionsInit(&opts);
    CAT_SurfDeformPialProfile(p, vol, nim, &opts);
    radius_stats(p, &mean, &min, &max);
    MU_ASSERT("glued sulcus: stops at the valley bottom r = 9", fabs(mean - 9.0) < 0.15);
    MU_ASSERT("glued sulcus: does not run on to the outer boundary", max < 9.5);

    delete_object(obj);
    free(vol);
    nifti_image_free(nim);
}

/* Two parallel sheets in the y/z plane at x = -0.5 (normal +x) and x = +0.5
 * (normal -x), 17 x 17 vertices each at 0.5 mm spacing. */
static void
make_facing_sheets(polygons_struct *p)
{
    const int m = 17;
    int s, i, j, t = 0;

    initialize_polygons(p, WHITE, NULL);
    p->n_points = 2 * m * m;
    p->points = (Point *)malloc(sizeof(Point) * p->n_points);
    p->normals = (Vector *)malloc(sizeof(Vector) * p->n_points);
    p->n_items = 2 * 2 * (m - 1) * (m - 1);
    p->end_indices = (int *)malloc(sizeof(int) * p->n_items);
    p->indices = (int *)malloc(sizeof(int) * 3 * p->n_items);

    for (s = 0; s < 2; s++)
        for (j = 0; j < m; j++)
            for (i = 0; i < m; i++)
                fill_Point(p->points[s * m * m + i + m * j], s ? 0.5 : -0.5,
                           -4.0 + 0.5 * i, -4.0 + 0.5 * j);

    for (s = 0; s < 2; s++)
        for (j = 0; j < m - 1; j++)
            for (i = 0; i < m - 1; i++)
            {
                int o = s * m * m;
                int a = o + i + m * j, b = a + 1, c = a + m, d = c + 1;
                /* (+y) x (+z) = +x for the first sheet, reversed for the second */
                int tri[2][3] = {{a, b, c}, {b, d, c}};
                int q, k;
                for (q = 0; q < 2; q++, t++)
                {
                    for (k = 0; k < 3; k++)
                        p->indices[3 * t + k] = s ? tri[q][2 - k] : tri[q][k];
                    p->end_indices[t] = 3 * (t + 1);
                }
            }
    compute_polygon_normals(p);
}

static void test_facing_walls_meet_in_the_middle(void)
{
    int n[3] = {13, 25, 25};
    nifti_image *nim = make_nim(n, 1);
    float *vol = fill_volume(nim, slab_phantom);
    polygons_struct p;
    CAT_PialProfileOptions opts;
    const int m = 17;
    int i, j, ok_order = 1, ok_middle = 1, normals_face = 1;
    double gap_min = 1e30, gap_max = -1e30;

    make_facing_sheets(&p);
    normals_face = Point_x(p.normals[0]) > 0.5 && Point_x(p.normals[m * m]) < -0.5;
    MU_ASSERT("sheet normals face each other", normals_face);

    /* each sheet alone would move 1.9 mm, through the other one */
    CAT_PialProfileOptionsInit(&opts);
    CAT_SurfDeformPialProfile(&p, vol, nim, &opts);

    for (j = 4; j < m - 4; j++)
        for (i = 4; i < m - 4; i++)
        {
            double xa = Point_x(p.points[i + m * j]);
            double xb = Point_x(p.points[m * m + i + m * j]);
            if (xa >= xb)
                ok_order = 0;
            if (xa < -0.3 || xb > 0.3)
                ok_middle = 0;
            gap_min = fmin(gap_min, xb - xa);
            gap_max = fmax(gap_max, xb - xa);
        }
    MU_ASSERT("facing walls do not cross", ok_order);
    MU_ASSERT("facing walls both move to the middle", ok_middle);
    MU_ASSERT("facing walls keep a small gap", gap_min > 0.0 && gap_max < 0.4);

    free(p.points);
    free(p.normals);
    free(p.end_indices);
    free(p.indices);
    free(vol);
    nifti_image_free(nim);
}

int main(void)
{
    MU_RUN_TEST(test_reaches_boundary_from_far_inside);
    MU_RUN_TEST(test_holds_without_boundary_in_range);
    MU_RUN_TEST(test_stops_at_valley_bottom);
    MU_RUN_TEST(test_facing_walls_meet_in_the_middle);
    printf("%d tests run, %d failed\n", tests_run, tests_failed);
    return tests_failed ? 1 : 0;
}
