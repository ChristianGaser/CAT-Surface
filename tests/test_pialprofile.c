#include "minunit.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "CAT_SurfPialProfile.h"
#include "CAT_Deform.h"

/*
 * Analytic phantoms on a 0.5 mm grid.  Tissue boundaries are linear ramps one
 * voxel wide, so a label of 1.5 marks the CSF/GM boundary exactly and every
 * expected surface position has a closed-form value.
 */

#define VX 0.5

/**
 * \brief Linear partial-volume ramp: 0 outside, 1 inside, 0.5 at d = 0.
 *
 * \param d (in) signed distance to the boundary in mm, positive inside
 * \return ramp value in [0, 1]
 */
static double
ramp(double d)
{
    double v = d / VX + 0.5;
    return v < 0.0 ? 0.0 : (v > 1.0 ? 1.0 : v);
}

/**
 * \brief Create a header-only NIfTI image with a diagonal sform.
 *
 * \param n    (in) dimensions
 * \param sign (in) +1 stores the axes RAS, -1 stores them LAS-like (all flipped)
 * \return new image; free with nifti_image_free()
 */
static nifti_image *
make_nim(const int n[3], int sign)
{
    nifti_image *nim = nifti_simple_init_nim();
    int k;

    nim->nx = nim->dim[1] = n[0];
    nim->ny = nim->dim[2] = n[1];
    nim->nz = nim->dim[3] = n[2];
    nim->nvox = (size_t)n[0] * n[1] * n[2];
    nim->dx = nim->pixdim[1] = VX;
    nim->dy = nim->pixdim[2] = VX;
    nim->dz = nim->pixdim[3] = VX;
    memset(&nim->sto_xyz, 0, sizeof(mat44));
    for (k = 0; k < 3; k++)
    {
        /* centre the grid on the origin */
        nim->sto_xyz.m[k][k] = (float)(sign * VX);
        nim->sto_xyz.m[k][3] = (float)(-sign * VX * (n[k] - 1) / 2.0);
    }
    nim->sto_xyz.m[3][3] = 1.0f;
    nim->sform_code = 1;
    return nim;
}

/**
 * \brief Fill a volume by evaluating f at the world position of each voxel.
 *
 * \param nim (in)  header defining the grid
 * \param f   (in)  phantom function of the world position
 * \return newly allocated volume
 */
static float *
fill_volume(const nifti_image *nim, double (*f)(double, double, double))
{
    float *vol = (float *)malloc(sizeof(float) * nim->nvox);
    int i, j, k;

    for (k = 0; k < nim->nz; k++)
        for (j = 0; j < nim->ny; j++)
            for (i = 0; i < nim->nx; i++)
            {
                double x = nim->sto_xyz.m[0][0] * i + nim->sto_xyz.m[0][3];
                double y = nim->sto_xyz.m[1][1] * j + nim->sto_xyz.m[1][3];
                double z = nim->sto_xyz.m[2][2] * k + nim->sto_xyz.m[2][3];
                vol[i + nim->nx * (j + nim->ny * k)] = (float)f(x, y, z);
            }
    return vol;
}

/* WM inside r = 6, pial boundary at r = 9 */
static double
sphere_phantom(double x, double y, double z)
{
    double r = sqrt(x * x + y * y + z * z);
    return 1.0 + ramp(9.0 - r) + ramp(6.0 - r);
}

/* Glued sulcus: GM continues to r = 11, with a valley down to 1.8 at r = 9 */
static double
valley_phantom(double x, double y, double z)
{
    double r = sqrt(x * x + y * y + z * z);
    return 1.0 + ramp(11.0 - r) + ramp(6.0 - r) - 0.2 * exp(-(r - 9.0) * (r - 9.0) / (2.0 * 0.16));
}

/* GM slab |x| < 1.4 between two CSF regions */
static double
slab_phantom(double x, double y, double z)
{
    (void)y;
    (void)z;
    return 1.0 + ramp(1.4 - fabs(x));
}

static void
make_sphere(polygons_struct *p, double radius)
{
    Point centre;
    fill_Point(centre, 0.0, 0.0, 0.0);
    create_tetrahedral_sphere(&centre, radius, radius, radius, 5120, p);
    compute_polygon_normals(p);
}

static void
radius_stats(const polygons_struct *p, double *mean, double *min, double *max)
{
    int i;
    *mean = 0.0;
    *min = 1e30;
    *max = -1e30;
    for (i = 0; i < p->n_points; i++)
    {
        double r = sqrt(Point_x(p->points[i]) * Point_x(p->points[i]) +
                        Point_y(p->points[i]) * Point_y(p->points[i]) +
                        Point_z(p->points[i]) * Point_z(p->points[i]));
        *mean += r;
        *min = r < *min ? r : *min;
        *max = r > *max ? r : *max;
    }
    *mean /= p->n_points;
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

    for (i = 0; i < n_points; i++)
    {
        double d = fabs(Point_x(result[0][i]) - Point_x(result[1][i])) +
                   fabs(Point_y(result[0][i]) - Point_y(result[1][i])) +
                   fabs(Point_z(result[0][i]) - Point_z(result[1][i]));
        max_diff = d > max_diff ? d : max_diff;
    }
    MU_ASSERT("surf_deform_dual: RAS and LAS storage give the same surface", max_diff < 1e-3);

    free(result[0]);
    free(result[1]);
}

int main(void)
{
    MU_RUN_TEST(test_reaches_boundary_from_far_inside);
    MU_RUN_TEST(test_holds_without_boundary_in_range);
    MU_RUN_TEST(test_stops_at_valley_bottom);
    MU_RUN_TEST(test_facing_walls_meet_in_the_middle);
    MU_RUN_TEST(test_deform_dual_orientation_invariant);
    printf("%d tests run, %d failed\n", tests_run, tests_failed);
    return tests_failed ? 1 : 0;
}
