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

int main(void) {
    MU_RUN_TEST(test_get_median);
    MU_RUN_TEST(test_get_sum);
    MU_RUN_TEST(test_get_sum_exclude_zeros);
    MU_RUN_TEST(test_get_mean);
    MU_RUN_TEST(test_get_std);
    printf("%d tests run, %d failed\n", tests_run, tests_failed);
    return tests_failed ? 1 : 0;
}
