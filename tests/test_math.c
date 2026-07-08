#include "minunit.h"
#include "math_subset.c"
#include <math.h>

static void test_get_median(void) {
    double arr[] = {3.0, 1.0, 2.0};
    double m = get_median_double(arr, 3, 0);
    MU_ASSERT("median should be 2", fabs(m - 2.0) < 1e-6);
}

static void test_get_sum(void) {
    double arr[] = {1.0, 2.0, 3.0};
    double s = get_sum_double(arr, 3, 0);
    MU_ASSERT("sum should be 6", fabs(s - 6.0) < 1e-6);
}

static void test_get_sum_exclude_zeros(void) {
    double arr[] = {1.0, 0.0, 2.0};
    double s = get_sum_double(arr, 3, 1);
    MU_ASSERT("sum should be 3 when excluding zeros", fabs(s - 3.0) < 1e-6);
}

static void test_get_mean(void) {
    double arr[] = {2.0, 4.0, 6.0};
    double m = get_mean_double(arr, 3, 0);
    MU_ASSERT("mean should be 4", fabs(m - 4.0) < 1e-6);
}

static void test_get_std(void) {
    double arr[] = {1.0, 3.0, 5.0};
    double s = get_std_double(arr, 3, 0);
    MU_ASSERT("std should be 2", fabs(s - 2.0) < 1e-6);
}

static void test_orthogonal_poly(void) {
    /* R: poly(1:4, 2) */
    double x[] = {1.0, 2.0, 3.0, 4.0};
    double p1[] = {-0.6708204, -0.2236068, 0.2236068, 0.6708204};
    double p2[] = { 0.5, -0.5, -0.5, 0.5 };
    double out[8];
    int i, ok, match = 1;
    ok = orthogonal_poly(x, 4, 2, out);
    MU_ASSERT("orthogonal_poly should succeed", ok == 1);
    for (i = 0; i < 4; i++) {
        if (fabs(out[i] - p1[i]) > 1e-5) match = 0;       /* column 0 */
        if (fabs(out[4 + i] - p2[i]) > 1e-5) match = 0;   /* column 1 */
    }
    MU_ASSERT("orthogonal_poly should match R poly(1:4,2)", match);
    /* degree >= #distinct points is degenerate -> failure */
    MU_ASSERT("orthogonal_poly rejects degenerate degree",
              orthogonal_poly(x, 4, 4, out) == 0);
}

int main(void) {
    MU_RUN_TEST(test_get_median);
    MU_RUN_TEST(test_get_sum);
    MU_RUN_TEST(test_get_sum_exclude_zeros);
    MU_RUN_TEST(test_get_mean);
    MU_RUN_TEST(test_get_std);
    MU_RUN_TEST(test_orthogonal_poly);
    printf("%d tests run, %d failed\n", tests_run, tests_failed);
    return tests_failed ? 1 : 0;
}
