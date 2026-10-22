#include "minunit.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <bicpl.h>
#include "CAT_Separate.h"

/*
 * Clusters of positive values, as CAT_SurfSeparateClusters writes them.  The
 * tool used to test the neighbour's slot number instead of the neighbour and
 * never looked past the first ring around a seed, so a cluster wider than
 * that came out in pieces, numbered arbitrarily.  Clusters are now grown by a
 * breadth-first search and numbered by size.
 */

static void
make_sphere(polygons_struct *p)
{
    Point centre;

    fill_Point(centre, 0.0, 0.0, 0.0);
    create_tetrahedral_sphere(&centre, 1.0, 1.0, 1.0, 20480, p);
    compute_polygon_normals(p);
}

static int
count_label(const int *labels, int n, int label)
{
    int i, count = 0;

    for (i = 0; i < n; i++)
        if (labels[i] == label)
            count++;
    return count;
}

static void test_wide_cluster_is_one(void)
{
    polygons_struct p;
    double *values;
    int *labels, i, n, n_clusters, n_wrong = 0;

    make_sphere(&p);
    n = p.n_points;
    values = (double *)malloc(sizeof(double) * n);
    labels = (int *)malloc(sizeof(int) * n);
    for (i = 0; i < n; i++)
        values[i] = (Point_z(p.points[i]) > 0.5) ? 1.0 : -1.0;

    n_clusters = label_value_clusters(&p, values, labels);
    for (i = 0; i < n; i++)
        if (labels[i] != ((values[i] > 0.0) ? 1 : 0))
            n_wrong++;

    MU_ASSERT("a cap many rings wide is one cluster", n_clusters == 1);
    MU_ASSERT("the cluster is the cap, no more and no less", n_wrong == 0);

    free(values);
    free(labels);
    delete_polygons(&p);
}

static void test_clusters_numbered_by_size(void)
{
    polygons_struct p;
    double *values;
    int *labels, i, n, n_clusters, n_small = 0, n_large = 0, south = 0;

    make_sphere(&p);
    n = p.n_points;
    values = (double *)malloc(sizeof(double) * n);
    labels = (int *)malloc(sizeof(int) * n);

    /* a small cap around the north pole, found first, and a larger one in
       the south; the band between them is exactly 0 */
    for (i = 0; i < n; i++)
    {
        double z = Point_z(p.points[i]);

        values[i] = 0.0;
        if (z > 0.7)
        {
            values[i] = 2.0;
            n_small++;
        }
        else if (z < -0.4)
        {
            values[i] = 0.5;
            n_large++;
        }
        if (z < Point_z(p.points[south]))
            south = i;
    }
    /* NaN inside the large cap is outside, but does not split it */
    values[south] = NAN;

    n_clusters = label_value_clusters(&p, values, labels);

    MU_ASSERT("two caps are two clusters", n_clusters == 2);
    MU_ASSERT("the larger cap is cluster 1", count_label(labels, n, 1) == n_large - 1);
    MU_ASSERT("the smaller cap is cluster 2", count_label(labels, n, 2) == n_small);
    MU_ASSERT("NaN is outside every cluster", labels[south] == 0);
    MU_ASSERT("values of 0 are outside every cluster",
              count_label(labels, n, 0) == n - n_large - n_small + 1);

    free(values);
    free(labels);
    delete_polygons(&p);
}

static void test_ties_keep_vertex_order(void)
{
    polygons_struct p;
    double *values;
    int *labels, i, n, north = 0, south = 0, low, high, n_clusters;

    make_sphere(&p);
    n = p.n_points;
    values = (double *)malloc(sizeof(double) * n);
    labels = (int *)malloc(sizeof(int) * n);
    for (i = 0; i < n; i++)
    {
        values[i] = -1.0;
        if (Point_z(p.points[i]) > Point_z(p.points[north]))
            north = i;
        if (Point_z(p.points[i]) < Point_z(p.points[south]))
            south = i;
    }
    values[north] = 1.0;
    values[south] = 1.0;
    low = (north < south) ? north : south;
    high = (north < south) ? south : north;

    n_clusters = label_value_clusters(&p, values, labels);

    MU_ASSERT("two isolated vertices are two clusters", n_clusters == 2);
    MU_ASSERT("of equal clusters the one with the lower vertex comes first",
              labels[low] == 1 && labels[high] == 2);

    free(values);
    free(labels);
    delete_polygons(&p);
}

static void test_extract_cluster_polygons(void)
{
    polygons_struct p, out;
    double *values;
    int *labels, *used, i, j, n, size, poly, ind, n_items, n_expected = 0;
    int n_used = 0, n_mismatch = 0;

    make_sphere(&p);
    n = p.n_points;
    values = (double *)malloc(sizeof(double) * n);
    labels = (int *)malloc(sizeof(int) * n);
    used = (int *)calloc(n, sizeof(int));
    for (i = 0; i < n; i++)
        values[i] = (Point_z(p.points[i]) > 0.5) ? 1.0 : -1.0;
    label_value_clusters(&p, values, labels);

    /* the polygons inside the cap, counted independently */
    for (poly = 0; poly < p.n_items; poly++)
    {
        int inside = 1;

        size = GET_OBJECT_SIZE(p, poly);
        for (j = 0; j < size; j++)
            if (labels[p.indices[POINT_INDEX(p.end_indices, poly, j)]] != 1)
                inside = 0;
        if (!inside)
            continue;
        n_expected++;
        for (j = 0; j < size; j++)
            used[p.indices[POINT_INDEX(p.end_indices, poly, j)]] = 1;
    }
    for (i = 0; i < n; i++)
        n_used += used[i];

    n_items = extract_label_polygons(&p, labels, 1, &out);
    MU_ASSERT("every polygon inside the cluster and no other", n_items == n_expected &&
              out.n_items == n_expected);
    MU_ASSERT("the vertices of those polygons", out.n_points == n_used);

    /* the vertices keep their order, and each polygon its vertices */
    for (i = 0, j = 0; i < n; i++)
    {
        if (!used[i])
            continue;
        if (!EQUAL_POINTS(out.points[j], p.points[i]))
            n_mismatch++;
        j++;
    }
    ind = 0;
    for (poly = 0; poly < p.n_items; poly++)
    {
        int inside = 1;

        size = GET_OBJECT_SIZE(p, poly);
        for (j = 0; j < size; j++)
            if (labels[p.indices[POINT_INDEX(p.end_indices, poly, j)]] != 1)
                inside = 0;
        if (!inside)
            continue;
        for (j = 0; j < size; j++, ind++)
            if (!EQUAL_POINTS(out.points[out.indices[ind]],
                              p.points[p.indices[POINT_INDEX(p.end_indices, poly, j)]]))
                n_mismatch++;
    }
    MU_ASSERT("the extracted surface is the cap of the input", n_mismatch == 0);
    delete_polygons(&out);

    /* a single vertex has no polygon of its own */
    for (i = 0; i < n; i++)
        values[i] = -1.0;
    values[0] = 1.0;
    label_value_clusters(&p, values, labels);
    n_items = extract_label_polygons(&p, labels, 1, &out);
    MU_ASSERT("a cluster without a complete polygon gives an empty surface",
              n_items == 0 && out.n_items == 0 && out.n_points == 0);
    delete_polygons(&out);

    free(values);
    free(labels);
    free(used);
    delete_polygons(&p);
}

int main(void)
{
    MU_RUN_TEST(test_wide_cluster_is_one);
    MU_RUN_TEST(test_clusters_numbered_by_size);
    MU_RUN_TEST(test_ties_keep_vertex_order);
    MU_RUN_TEST(test_extract_cluster_polygons);
    printf("%d tests run, %d failed\n", tests_run, tests_failed);
    return tests_failed ? 1 : 0;
}
