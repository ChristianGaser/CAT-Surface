#include "minunit.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <bicpl.h>
#include "CAT_Intersect.h"

/*
 * Local smoothing cannot separate two sheets that were driven through each
 * other -- the two sides of a thin gyral blade after a deformation.  The
 * phantom does the same to a sphere: its northern cap is pushed through the
 * southern hemisphere, so two sheets that are far apart on the mesh cross.
 * The sphere it came from is the way back.
 */

#define RADIUS 10.0

static void
make_sphere10(polygons_struct *p)
{
    Point centre;

    fill_Point(centre, 0.0, 0.0, 0.0);
    create_tetrahedral_sphere(&centre, RADIUS, RADIUS, RADIUS, 5120, p);
    compute_polygon_normals(p);
}

static void
push_cap_through(polygons_struct *p)
{
    int i;

    for (i = 0; i < p->n_points; i++)
        if (Point_z(p->points[i]) > 0.7 * RADIUS)
            Point_z(p->points[i]) -= 1.9 * RADIUS;
    compute_polygon_normals(p);
}

static int
count_pairs(polygons_struct *p)
{
    int *defects = (int *)calloc(p->n_points, sizeof(int));
    int *polydefects = (int *)calloc(p->n_items, sizeof(int));
    int n = find_selfintersections(p, defects, polydefects, 1);

    free(defects);
    free(polydefects);
    return n;
}

static void test_crossed_sheets_retreat_to_reference(void)
{
    object_struct *o_ref = create_object(POLYGONS);
    object_struct *o_plain = create_object(POLYGONS);
    object_struct *o_fix = create_object(POLYGONS);
    polygons_struct *ref = get_polygons_ptr(o_ref);
    polygons_struct *plain = get_polygons_ptr(o_plain);
    polygons_struct *fix = get_polygons_ptr(o_fix);
    int i, n_before, left_plain, left_ref, equator_moved = 0;

    make_sphere10(ref);
    make_sphere10(plain);
    make_sphere10(fix);
    push_cap_through(plain);
    push_cap_through(fix);

    n_before = count_pairs(fix);
    MU_ASSERT("the pushed cap crosses the sphere", n_before > 100);

    /* the phantom has to defeat the smoothing, or it proves nothing */
    left_plain = remove_intersections_iter(plain, 10, 50, 0);
    MU_ASSERT("smoothing alone leaves the crossed sheets",
              left_plain > 0 && count_pairs(plain) > 0);

    left_ref = remove_intersections_ref(fix, ref->points, 10, 50, 0);
    MU_ASSERT("retreating to the reference removes them",
              left_ref == 0 && count_pairs(fix) == 0);

    /* and only the defect and its neighbourhood moved */
    for (i = 0; i < fix->n_points; i++)
        if (fabs(Point_z(ref->points[i])) < 0.3 * RADIUS &&
            (Point_x(fix->points[i]) != Point_x(ref->points[i]) ||
             Point_y(fix->points[i]) != Point_y(ref->points[i]) ||
             Point_z(fix->points[i]) != Point_z(ref->points[i])))
            equator_moved++;
    MU_ASSERT("vertices far from the defect are untouched", equator_moved == 0);

    delete_object(o_ref);
    delete_object(o_plain);
    delete_object(o_fix);
}

static void test_no_reference_is_plain_repair(void)
{
    object_struct *o_a = create_object(POLYGONS);
    object_struct *o_b = create_object(POLYGONS);
    polygons_struct *a = get_polygons_ptr(o_a);
    polygons_struct *b = get_polygons_ptr(o_b);
    int i, same = 1, ra, rb;

    make_sphere10(a);
    make_sphere10(b);
    push_cap_through(a);
    push_cap_through(b);

    ra = remove_intersections_iter(a, 10, 50, 0);
    rb = remove_intersections_ref(b, NULL, 10, 50, 0);
    for (i = 0; i < a->n_points; i++)
        if (memcmp(&a->points[i], &b->points[i], sizeof(Point)) != 0)
            same = 0;
    MU_ASSERT("without a reference the repair is unchanged", ra == rb && same);

    delete_object(o_a);
    delete_object(o_b);
}

static void test_clean_surface_untouched(void)
{
    object_struct *o_a = create_object(POLYGONS);
    object_struct *o_orig = create_object(POLYGONS);
    object_struct *o_ref = create_object(POLYGONS);
    polygons_struct *a = get_polygons_ptr(o_a);
    polygons_struct *orig = get_polygons_ptr(o_orig);
    polygons_struct *ref = get_polygons_ptr(o_ref);
    int i, same = 1;

    make_sphere10(a);
    make_sphere10(orig);
    make_sphere10(ref);
    /* a reference that differs everywhere must not be used at all */
    for (i = 0; i < ref->n_points; i++)
        Point_x(ref->points[i]) += 1.0f;

    MU_ASSERT("a clean surface reports no defects",
              remove_intersections_ref(a, ref->points, 10, 50, 0) == 0);
    for (i = 0; i < a->n_points; i++)
        if (memcmp(&a->points[i], &orig->points[i], sizeof(Point)) != 0)
            same = 0;
    MU_ASSERT("and is left where it was", same);

    delete_object(o_a);
    delete_object(o_orig);
    delete_object(o_ref);
}

int main(void)
{
    MU_RUN_TEST(test_crossed_sheets_retreat_to_reference);
    MU_RUN_TEST(test_no_reference_is_plain_repair);
    MU_RUN_TEST(test_clean_surface_untouched);
    printf("%d tests run, %d failed\n", tests_run, tests_failed);
    return tests_failed ? 1 : 0;
}
