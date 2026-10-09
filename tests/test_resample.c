#include "minunit.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <bicpl.h>
#include "CAT_Resample.h"
#include "CAT_Surf.h"

/*
 * Areal interpolation has to preserve the sum of an extensive measure such as
 * the local surface area.  It used to drop the source vertices that no target
 * vertex referenced, which happens wherever the target mesh is coarser than
 * the source mesh: resampling native surfaces to the 32k mesh lost 5-25% of
 * their area, the more the larger the brain.  The source spheres are rotated,
 * so that their vertices do not coincide with those of the target spheres.
 */

#define RADIUS 100.0

static object_struct *
make_sphere(int n_triangles, double angle)
{
    object_struct *obj = create_object(POLYGONS);
    polygons_struct *p = get_polygons_ptr(obj);
    double ca = cos(angle), sa = sin(angle);
    double cb = cos(0.7 * angle), sb = sin(0.7 * angle);
    Point centre;
    int i;

    fill_Point(centre, 0.0, 0.0, 0.0);
    create_tetrahedral_sphere(&centre, RADIUS, RADIUS, RADIUS, n_triangles, p);

    /* rotate about z, then about x */
    for (i = 0; i < p->n_points; i++)
    {
        double x = Point_x(p->points[i]);
        double y = Point_y(p->points[i]);
        double z = Point_z(p->points[i]);
        double x1 = ca * x - sa * y, y1 = sa * x + ca * y;

        fill_Point(p->points[i], x1, cb * y1 - sb * z, sb * y1 + cb * z);
    }
    compute_polygon_normals(p);
    return obj;
}

/* vertex areas of the sphere, modulated by a smooth pattern */
static double *
make_values(polygons_struct *p, int modulate)
{
    double *v = (double *)malloc(sizeof(double) * p->n_points);
    int i;

    get_area_of_points(p, v);
    if (modulate)
        for (i = 0; i < p->n_points; i++)
            v[i] *= 1.0 + 0.5 * Point_x(p->points[i]) * Point_y(p->points[i]) / (RADIUS * RADIUS);
    return v;
}

static double
sum(const double *v, int n)
{
    double s = 0.0;
    int i;

    for (i = 0; i < n; i++)
        s += v[i];
    return s;
}

/* resample the source values to the target, return the relative change of
   their sum and optionally the relative RMS error to reference values */
static double
areal_resample(int n_src, int n_trg, double angle, int modulate,
               int use_values_sphere, double *rrmse)
{
    object_struct *src_obj = make_sphere(n_src, angle);
    object_struct *trg_obj = make_sphere(n_trg, 0.0);
    polygons_struct *src = get_polygons_ptr(src_obj);
    polygons_struct *trg = get_polygons_ptr(trg_obj);
    double *in = make_values(src, modulate);
    double *out = (double *)malloc(sizeof(double) * trg->n_points);
    double s_in, s_out;
    int i;

    if (use_values_sphere)
        resample_values_sphere(src, trg, in, out, 1, 1);
    else
    {
        object_struct **res = resample_surface_to_target_sphere(src, src, trg, in, out, 0, 1);
        delete_object_list(1, res);
    }

    s_in = sum(in, src->n_points);
    s_out = sum(out, trg->n_points);

    /* a uniform area density has to give the target vertex areas */
    if (rrmse != NULL)
    {
        double *ref = make_values(trg, 0);
        double se = 0.0;

        for (i = 0; i < trg->n_points; i++)
            se += (out[i] - ref[i]) * (out[i] - ref[i]);
        *rrmse = sqrt(se / trg->n_points) / (sum(ref, trg->n_points) / trg->n_points);
        free(ref);
    }

    free(in);
    free(out);
    delete_object(src_obj);
    delete_object(trg_obj);
    return fabs(s_out - s_in) / s_in;
}

static void test_areal_coarser_target(void)
{
    double rrmse;
    double d = areal_resample(81920, 20480, 0.3, 1, 0, NULL);

    MU_ASSERT("areal resampling to a coarser mesh preserves the sum", d < 1e-9);

    d = areal_resample(81920, 20480, 0.3, 0, 0, &rrmse);
    MU_ASSERT("areal resampling to a coarser mesh preserves the sum (uniform)", d < 1e-9);
    printf("coarser target: relative RMS error to target vertex areas %.4f\n", rrmse);
    MU_ASSERT("areal resampling to a coarser mesh stays local", rrmse < 0.1);
}

static void test_areal_finer_target(void)
{
    double rrmse;
    double d = areal_resample(20480, 81920, 0.3, 1, 0, NULL);

    MU_ASSERT("areal resampling to a finer mesh preserves the sum", d < 1e-9);

    d = areal_resample(20480, 81920, 0.3, 0, 0, &rrmse);
    MU_ASSERT("areal resampling to a finer mesh preserves the sum (uniform)", d < 1e-9);
    printf("finer target: relative RMS error to target vertex areas %.4f\n", rrmse);
    MU_ASSERT("areal resampling to a finer mesh stays local", rrmse < 0.1);
}

static void test_areal_same_mesh(void)
{
    object_struct *obj = make_sphere(20480, 0.0);
    polygons_struct *p = get_polygons_ptr(obj);
    double *in = make_values(p, 1);
    double *out = (double *)malloc(sizeof(double) * p->n_points);
    double err = 0.0;
    object_struct **res;
    int i;

    res = resample_surface_to_target_sphere(p, p, p, in, out, 0, 1);
    for (i = 0; i < p->n_points; i++)
        err = fmax(err, fabs(out[i] - in[i]) / in[i]);
    MU_ASSERT("areal resampling to the same mesh keeps the values", err < 1e-6);

    delete_object_list(1, res);
    free(in);
    free(out);
    delete_object(obj);
}

static void test_areal_values_sphere(void)
{
    double d = areal_resample(81920, 20480, 0.3, 1, 1, NULL);

    MU_ASSERT("resample_values_sphere with areal interpolation preserves the sum", d < 1e-9);
}

int main(void)
{
    MU_RUN_TEST(test_areal_coarser_target);
    MU_RUN_TEST(test_areal_finer_target);
    MU_RUN_TEST(test_areal_same_mesh);
    MU_RUN_TEST(test_areal_values_sphere);
    printf("%d tests run, %d failed\n", tests_run, tests_failed);
    return tests_failed ? 1 : 0;
}
