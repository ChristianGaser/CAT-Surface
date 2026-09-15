#include "minunit.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <bicpl.h>
#include "CAT_NiftiLib.h"
#include "CAT_Deform.h"
#include "CAT_SurfLaplacian.h"
#include "test_phantoms.h"

/*
 * The gradient forces of the mesh deformations must be defined in world space:
 * the same phantom stored with RAS and with all axes flipped has to give the
 * same surface.  Before, the voxel-axis gradient was dotted with the world
 * normal, which flips the sign of the gradient term per negatively stored
 * axis.
 */

static double
max_point_diff(const Point *a, const Point *b, int n_points)
{
    double max_diff = 0.0;
    int i;

    for (i = 0; i < n_points; i++)
    {
        double d = fabs(Point_x(a[i]) - Point_x(b[i])) + fabs(Point_y(a[i]) - Point_y(b[i])) +
                   fabs(Point_z(a[i]) - Point_z(b[i]));
        max_diff = d > max_diff ? d : max_diff;
    }
    return max_diff;
}

/* surf_deform must not depend on how the image axes are stored */
static void test_deform_orientation_invariant(void)
{
    int n[3] = {64, 64, 64};
    int sgn, i, n_points = 0;
    double w[3] = {0.1, 0.1, 1.0}, max_diff = 0.0;
    Point *result[2] = {NULL, NULL};

    for (sgn = 0; sgn < 2; sgn++)
    {
        nifti_image *nim = make_nim(n, sgn ? 1 : -1);
        float *vol = fill_volume(nim, sphere_phantom);
        object_struct *obj = create_object(POLYGONS);
        polygons_struct *p = get_polygons_ptr(obj);

        /* few iterations: later on the discrete near-intersection revert
         * amplifies rounding differences between the two storages */
        make_sphere(p, 8.0);
        n_points = p->n_points;
        surf_deform(p, vol, nim, w, 0.2, 1.5f, 3, 0, 0);

        result[sgn] = (Point *)malloc(sizeof(Point) * n_points);
        memcpy(result[sgn], p->points, sizeof(Point) * n_points);

        delete_object(obj);
        free(vol);
        nifti_image_free(nim);
    }

    max_diff = max_point_diff(result[0], result[1], n_points);
    MU_ASSERT("surf_deform: RAS and LAS storage give the same surface", max_diff < 1e-3);

    free(result[0]);
    free(result[1]);
}

/* surf_deform_dual must not depend on how the image axes are stored */
static void test_deform_dual_orientation_invariant(void)
{
    int n[3] = {64, 64, 64};
    int sgn, i, n_points = 0;
    double w[3] = {0.05, 0.05, 0.05}, max_diff = 0.0;
    Point *result[2] = {NULL, NULL};

    for (sgn = 0; sgn < 2; sgn++)
    {
        nifti_image *nim = make_nim(n, sgn ? 1 : -1);
        float *vol = fill_volume(nim, sphere_phantom);
        object_struct *o_pial = create_object(POLYGONS);
        object_struct *o_orig = create_object(POLYGONS);
        polygons_struct *pial = get_polygons_ptr(o_pial);
        polygons_struct *orig = get_polygons_ptr(o_orig);
        double *thick;

        make_sphere(pial, 8.0);
        make_sphere(orig, 7.5);
        n_points = pial->n_points;
        thick = (double *)malloc(sizeof(double) * n_points);
        for (i = 0; i < n_points; i++)
            thick[i] = 3.0;

        surf_deform_dual(pial, NULL, orig, vol, nim, w, 0.2, 1.25f, 2.7f, thick, 10, 0);

        result[sgn] = (Point *)malloc(sizeof(Point) * n_points);
        memcpy(result[sgn], pial->points, sizeof(Point) * n_points);

        free(thick);
        delete_object(o_pial);
        delete_object(o_orig);
        free(vol);
        nifti_image_free(nim);
    }

    max_diff = max_point_diff(result[0], result[1], n_points);
    MU_ASSERT("surf_deform_dual: RAS and LAS storage give the same surface", max_diff < 1e-3);

    free(result[0]);
    free(result[1]);
}

/* ADE streamlines are traced in world space: both storages must converge.
 * They stop at phi >= 0.999 / <= 0.001, i.e. at the far end of the one-voxel
 * partial-volume ramp, which puts the surfaces at about r = 9.35 / 5.65 around
 * the boundaries at 9 / 6.  A thickness of 2 mm puts the fallback (central
 * surface +/- half the thickness) at r = 8.5 / 6.5, so streamlines that fail
 * cannot pass. */
static void test_ade_orientation_invariant(void)
{
    int n[3] = {64, 64, 64};
    int sgn, i, n_points = 0;
    double mean_r[2][2];
    Point *white[2] = {NULL, NULL};

    for (sgn = 0; sgn < 2; sgn++)
    {
        nifti_image *nim = make_nim(n, sgn ? 1 : -1);
        float *vol = fill_volume(nim, sphere_phantom);
        object_struct *o_c = create_object(POLYGONS);
        object_struct *o_p = create_object(POLYGONS);
        object_struct *o_w = create_object(POLYGONS);
        polygons_struct *central = get_polygons_ptr(o_c);
        double *thick, min, max;

        make_sphere(central, 7.5);
        n_points = central->n_points;
        thick = (double *)malloc(sizeof(double) * n_points);
        for (i = 0; i < n_points; i++)
            thick[i] = 2.0;

        surf_ade_pial_white(central, vol, nim, 1.5f, 2.5f, thick,
                            get_polygons_ptr(o_p), get_polygons_ptr(o_w), 0);
        radius_stats(get_polygons_ptr(o_p), &mean_r[sgn][0], &min, &max);
        radius_stats(get_polygons_ptr(o_w), &mean_r[sgn][1], &min, &max);

        white[sgn] = (Point *)malloc(sizeof(Point) * n_points);
        memcpy(white[sgn], get_polygons_ptr(o_w)->points, sizeof(Point) * n_points);

        free(thick);
        delete_object(o_c);
        delete_object(o_p);
        delete_object(o_w);
        free(vol);
        nifti_image_free(nim);
    }

    MU_ASSERT("ADE: white streamlines reach the GM/WM boundary when stored LAS",
              mean_r[0][1] > 5.4 && mean_r[0][1] < 6.1);
    MU_ASSERT("ADE: pial streamlines reach the CSF/GM boundary when stored LAS",
              mean_r[0][0] > 8.9 && mean_r[0][0] < 9.6);
    MU_ASSERT("ADE: RAS and LAS storage give the same mean radii",
              fabs(mean_r[0][0] - mean_r[1][0]) < 0.01 && fabs(mean_r[0][1] - mean_r[1][1]) < 0.01);
    /* streamlines may stop one 0.1 mm integration step apart; max_point_diff()
     * sums |dx| + |dy| + |dz|, which is at most sqrt(3) times that */
    MU_ASSERT("ADE: RAS and LAS storage give the same white surface",
              max_point_diff(white[0], white[1], n_points) < 0.18);

    free(white[0]);
    free(white[1]);
}

int main(void)
{
    MU_RUN_TEST(test_deform_orientation_invariant);
    MU_RUN_TEST(test_deform_dual_orientation_invariant);
    MU_RUN_TEST(test_ade_orientation_invariant);
    printf("%d tests run, %d failed\n", tests_run, tests_failed);
    return tests_failed ? 1 : 0;
}
