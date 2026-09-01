#include "minunit.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "CAT_SurfInfo.h"
#include "CAT_SurfaceIO.h"

/*
 * The reference bodies are built by hand so that every reported quantity has a
 * closed-form value: a unit octahedron (V=6, F=8, E=12, chi=2, genus 0), an
 * open square patch (two triangles, four boundary edges) and two disjoint
 * octahedra (chi=4, two components).
 */

static const double octa_pts[6][3] = {
    { 1, 0, 0}, {-1, 0, 0}, {0,  1, 0},
    {0, -1, 0}, { 0, 0, 1}, {0, 0, -1}
};

static const int octa_tri[8][3] = {
    {0, 2, 4}, {2, 1, 4}, {1, 3, 4}, {3, 0, 4},
    {2, 0, 5}, {1, 2, 5}, {3, 1, 5}, {0, 3, 5}
};

/**
 * \brief Build a mesh from vertex and triangle tables.
 *
 * \param polygons (out) mesh to fill; free with free_mesh()
 * \param pts      (in)  n_pts x 3 vertex coordinates
 * \param n_pts    (in)  number of vertices
 * \param tri      (in)  n_tri x 3 vertex indices
 * \param n_tri    (in)  number of triangles
 * \return void
 */
static void
make_mesh(polygons_struct *polygons, const double (*pts)[3], int n_pts,
          const int (*tri)[3], int n_tri)
{
    int i, k;

    initialize_polygons(polygons, WHITE, NULL);

    polygons->n_points = n_pts;
    polygons->points = (Point *) malloc(sizeof(Point) * n_pts);
    polygons->normals = (Vector *) malloc(sizeof(Vector) * n_pts);
    for (i = 0; i < n_pts; i++) {
        fill_Point(polygons->points[i], pts[i][0], pts[i][1], pts[i][2]);
        fill_Vector(polygons->normals[i], 0.0, 0.0, 1.0);
    }

    polygons->n_items = n_tri;
    polygons->end_indices = (int *) malloc(sizeof(int) * n_tri);
    polygons->indices = (int *) malloc(sizeof(int) * n_tri * 3);
    for (i = 0; i < n_tri; i++) {
        polygons->end_indices[i] = 3 * (i + 1);
        for (k = 0; k < 3; k++)
            polygons->indices[3 * i + k] = tri[i][k];
    }
}

/**
 * \brief Release a mesh built by make_mesh().
 *
 * \param polygons (in/out) mesh to free
 * \return void
 */
static void
free_mesh(polygons_struct *polygons)
{
    free(polygons->points);
    free(polygons->normals);
    free(polygons->end_indices);
    free(polygons->indices);
    polygons->points = NULL;
    polygons->normals = NULL;
    polygons->end_indices = NULL;
    polygons->indices = NULL;
    polygons->n_points = 0;
    polygons->n_items = 0;
}

static void
test_octahedron_counts(void)
{
    polygons_struct polygons;
    surf_info_struct info;

    make_mesh(&polygons, octa_pts, 6, octa_tri, 8);
    compute_surf_info(&polygons, &info, 1, 0);

    MU_ASSERT("6 vertices", info.n_points == 6);
    MU_ASSERT("8 faces", info.n_polygons == 8);
    MU_ASSERT("8 triangular faces", info.n_triangles == 8);
    MU_ASSERT("no non-triangular faces", info.n_nontriangular == 0);
    MU_ASSERT("12 edges", info.n_edges == 12);
    MU_ASSERT("V - E + F = 2", info.euler == 2);
    MU_ASSERT("genus 0", info.genus == 0);
    MU_ASSERT("one component", info.n_components == 1);
    MU_ASSERT("closed", info.closed == 1);
    MU_ASSERT("no boundary edges", info.n_boundary_edges == 0);
    MU_ASSERT("no non-manifold edges", info.n_nonmanifold_edges == 0);
    MU_ASSERT("no degenerate faces", info.n_degenerate_polygons == 0);
    MU_ASSERT("no unreferenced vertices", info.n_unreferenced_points == 0);
    MU_ASSERT("no self-intersections", info.n_self_intersections == 0);
    MU_ASSERT("no triangle intersections", info.n_triangle_hits == 0);

    free_mesh(&polygons);
}

static void
test_octahedron_geometry(void)
{
    polygons_struct polygons;
    surf_info_struct info;

    make_mesh(&polygons, octa_pts, 6, octa_tri, 8);
    compute_surf_info(&polygons, &info, 0, 0);

    /* eight equilateral triangles of edge sqrt(2): area = 8 * sqrt(3) / 2 */
    MU_ASSERT("surface area 4*sqrt(3)",
              fabs(info.surface_area - 4.0 * sqrt(3.0)) < 1e-5);
    /* volume of a regular octahedron with circumradius 1 is 4/3 */
    MU_ASSERT("volume 4/3", fabs(fabs(info.volume) - 4.0 / 3.0) < 1e-5);
    MU_ASSERT("all faces have the same area",
              fabs(info.area_max - info.area_min) < 1e-6);
    MU_ASSERT("all edges have length sqrt(2)",
              fabs(info.edge_mean - sqrt(2.0)) < 1e-5 &&
              fabs(info.edge_min - info.edge_max) < 1e-6);
    MU_ASSERT("bounding box is the unit cube",
              fabs(info.bounds[0] + 1.0) < 1e-6 &&
              fabs(info.bounds[5] - 1.0) < 1e-6);
    MU_ASSERT("centroid at the origin",
              fabs(info.centroid[0]) < 1e-6 &&
              fabs(info.centroid[1]) < 1e-6 &&
              fabs(info.centroid[2]) < 1e-6);

    free_mesh(&polygons);
}

static void
test_open_patch(void)
{
    /* unit square split into two triangles */
    static const double pts[4][3] = {
        {0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}
    };
    static const int tri[2][3] = { {0, 1, 2}, {0, 2, 3} };
    polygons_struct polygons;
    surf_info_struct info;

    make_mesh(&polygons, pts, 4, tri, 2);
    compute_surf_info(&polygons, &info, 0, 0);

    MU_ASSERT("5 edges", info.n_edges == 5);
    MU_ASSERT("chi of a disc is 1", info.euler == 1);
    MU_ASSERT("4 boundary edges", info.n_boundary_edges == 4);
    MU_ASSERT("not closed", info.closed == 0);
    MU_ASSERT("area of the unit square", fabs(info.surface_area - 1.0) < 1e-6);

    free_mesh(&polygons);
}

static void
test_two_components(void)
{
    double pts[12][3];
    int tri[16][3];
    polygons_struct polygons;
    surf_info_struct info;
    int i, k;

    for (i = 0; i < 6; i++) {
        for (k = 0; k < 3; k++) {
            pts[i][k] = octa_pts[i][k];
            pts[i + 6][k] = octa_pts[i][k] + (k == 0 ? 10.0 : 0.0);
        }
    }
    for (i = 0; i < 8; i++) {
        for (k = 0; k < 3; k++) {
            tri[i][k] = octa_tri[i][k];
            tri[i + 8][k] = octa_tri[i][k] + 6;
        }
    }

    make_mesh(&polygons, (const double (*)[3]) pts, 12,
              (const int (*)[3]) tri, 16);
    compute_surf_info(&polygons, &info, 0, 0);

    MU_ASSERT("two components", info.n_components == 2);
    MU_ASSERT("chi is additive over components", info.euler == 4);
    MU_ASSERT("both octahedra are closed", info.closed == 1);
    MU_ASSERT("genus stays 0", info.genus == 0);
    MU_ASSERT("area of two octahedra",
              fabs(info.surface_area - 8.0 * sqrt(3.0)) < 1e-5);

    free_mesh(&polygons);
}

static void
test_self_intersection_detected(void)
{
    /* two triangles crossing each other in the middle */
    static const double pts[6][3] = {
        {-1, 0, 0}, {1, 0, 0}, {0, 2, 0},
        {0, 1, -1}, {0, 1, 1}, {0, -1, 0}
    };
    static const int tri[2][3] = { {0, 1, 2}, {3, 4, 5} };
    polygons_struct polygons;
    surf_info_struct info;

    make_mesh(&polygons, pts, 6, tri, 2);
    compute_surf_info(&polygons, &info, 1, 0);

    /* the raw hit count double-counts a pair once per shared octree cell, so
       only its sign is asserted here; the derived counts are exact */
    MU_ASSERT("the crossing pair is found", info.n_triangle_hits >= 1);
    MU_ASSERT("it is one region", info.n_self_intersections == 1);
    MU_ASSERT("both triangles are marked", info.n_intersecting_polygons == 2);

    free_mesh(&polygons);
}

static void
test_skipped_intersections(void)
{
    polygons_struct polygons;
    surf_info_struct info;

    make_mesh(&polygons, octa_pts, 6, octa_tri, 8);
    compute_surf_info(&polygons, &info, 0, 0);

    MU_ASSERT("skipped test is flagged, not reported as zero",
              info.n_self_intersections == CAT_SURFINFO_SKIPPED &&
              info.n_triangle_hits == CAT_SURFINFO_SKIPPED &&
              info.n_intersecting_polygons == CAT_SURFINFO_SKIPPED);

    free_mesh(&polygons);
}

static void
test_gifti_data_arrays(void)
{
    /* a .gii carries the mesh in two DataArrays but may embed per-vertex data
       next to it, which is what input_gifti_darrays() has to surface */
    char file[] = "test_surfinfo_tmp.gii";
    /* a stack object, so make_mesh keeps ownership of the arrays it
       allocated with malloc and free_mesh can release them */
    object_struct obj;
    object_struct *object = &obj;
    polygons_struct *polygons;
    gifti_darray_info *arrays = NULL;
    double values[6];
    int n_arrays = 0, i;
    int pointset = -1, triangle = -1, shape = -1;

    obj.object_type = POLYGONS;
    obj.visibility = TRUE;
    polygons = get_polygons_ptr(object);
    make_mesh(polygons, octa_pts, 6, octa_tri, 8);

    for (i = 0; i < 6; i++)
        values[i] = (double) (i + 1);
    values[5] = 0.0 / 0.0;  /* a masked-out vertex */

    if (output_graphics_any_format(file, ASCII_FORMAT, 1, &object,
                                   values) != OK) {
        MU_ASSERT("the test file could be written", 0);
        return;
    }

    MU_ASSERT("the data arrays could be listed",
              input_gifti_darrays(file, &n_arrays, &arrays) == OK);
    MU_ASSERT("mesh plus one data array", n_arrays == 3);

    for (i = 0; i < n_arrays; i++) {
        if (arrays[i].intent == NIFTI_INTENT_POINTSET)
            pointset = i;
        else if (arrays[i].intent == NIFTI_INTENT_TRIANGLE)
            triangle = i;
        else if (arrays[i].intent == NIFTI_INTENT_SHAPE)
            shape = i;
    }

    MU_ASSERT("the three intents are found",
              pointset >= 0 && triangle >= 0 && shape >= 0);

    if (pointset >= 0)
        MU_ASSERT("the coordinates are 6 x 3",
                  arrays[pointset].dims[0] == 6 && arrays[pointset].dims[1] == 3);
    if (triangle >= 0)
        MU_ASSERT("the triangles are 8 x 3",
                  arrays[triangle].dims[0] == 8 && arrays[triangle].dims[1] == 3);

    if (shape >= 0) {
        MU_ASSERT("the embedded data has one value per vertex",
                  arrays[shape].n_values == 6);
        MU_ASSERT("the NaN is counted", arrays[shape].n_nonfinite == 1);
        MU_ASSERT("statistics are computed", arrays[shape].has_range == 1);
        /* the NaN must be left out rather than poison the summary */
        MU_ASSERT("mean of the finite values",
                  fabs(arrays[shape].mean - 3.0) < 1e-5);
        MU_ASSERT("min of the finite values",
                  fabs(arrays[shape].min - 1.0) < 1e-5);
        MU_ASSERT("max of the finite values",
                  fabs(arrays[shape].max - 5.0) < 1e-5);
    }

    /* the mesh arrays carry no separate data, so they get no statistics */
    if (pointset >= 0)
        MU_ASSERT("the mesh arrays are not summarized",
                  arrays[pointset].has_range == 0);

    free(arrays);
    free_mesh(polygons);
    remove(file);
}

int main(void)
{
    MU_RUN_TEST(test_octahedron_counts);
    MU_RUN_TEST(test_octahedron_geometry);
    MU_RUN_TEST(test_open_patch);
    MU_RUN_TEST(test_two_components);
    MU_RUN_TEST(test_self_intersection_detected);
    MU_RUN_TEST(test_skipped_intersections);
    MU_RUN_TEST(test_gifti_data_arrays);
    printf("%d tests run, %d failed\n", tests_run, tests_failed);
    return tests_failed ? 1 : 0;
}
