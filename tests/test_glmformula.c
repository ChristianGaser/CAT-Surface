#include "minunit.h"
#include "../Lib/CAT_GlmFormula.c"
#include <stdio.h>
#include <string.h>
#include <math.h>

/* write a small term file, one value per line */
static void wf(const char *name, const char *content)
{
    FILE *f = fopen(name, "w");
    fputs(content, f);
    fclose(f);
}

static void setup_files(void)
{
    wf("grp.csv", "a\na\nb\nb\n");     /* factor, 2 levels          */
    wf("cov.csv", "10\n20\n30\n40\n"); /* covariate                 */
    wf("num.csv", "1\n1\n2\n2\n");     /* numeric, forced to factor */
    wf("fa.csv",  "a\na\nb\nb\n");
    wf("fb.csv",  "p\nq\np\nq\n");
    wf("fc.csv",  "m\nm\nn\nn\n");
    wf("bad.csv", "a\nb\n");           /* only 2 rows               */
}

static void cleanup_files(void)
{
    remove("grp.csv"); remove("cov.csv"); remove("num.csv");
    remove("fa.csv");  remove("fb.csv");  remove("fc.csv");
    remove("bad.csv");
}

/* 1. simple group comparison: intercept + treatment contrast */
static void test_group_comparison(void)
{
    GlmDesign d;
    int ok = glm_build_design("~ grp.csv", 4, &d);
    MU_ASSERT("group: build should succeed", ok == 1);
    MU_ASSERT("group: 2 columns", d.n_beta == 2);
    MU_ASSERT("group: col0 intercept",
              strcmp(d.colnames[0], "(Intercept)") == 0);
    MU_ASSERT("group: col1 treatment label (ref = a)",
              strcmp(d.colnames[1], "grp[b]") == 0);
    MU_ASSERT("group: level a coded 0", d.G[0][1] == 0.0);
    MU_ASSERT("group: level b coded 1", d.G[2][1] == 1.0);
    MU_ASSERT("group: intercept all ones", d.G[0][0] == 1.0 && d.G[3][0] == 1.0);
    glm_free_design(&d);
}

/* 2. cell-means coding when the intercept is dropped */
static void test_cell_means(void)
{
    GlmDesign d;
    int ok = glm_build_design("~ 0 + grp.csv", 4, &d);
    MU_ASSERT("cellmeans: build should succeed", ok == 1);
    MU_ASSERT("cellmeans: 2 columns (both levels)", d.n_beta == 2);
    MU_ASSERT("cellmeans: col0 label", strcmp(d.colnames[0], "grp[a]") == 0);
    MU_ASSERT("cellmeans: col1 label", strcmp(d.colnames[1], "grp[b]") == 0);
    MU_ASSERT("cellmeans: no intercept, a=(1,0)",
              d.G[0][0] == 1.0 && d.G[0][1] == 0.0);
    MU_ASSERT("cellmeans: b=(0,1)",
              d.G[2][0] == 0.0 && d.G[2][1] == 1.0);
    glm_free_design(&d);
}

/* 3. regression: numeric file becomes a covariate column */
static void test_regression(void)
{
    GlmDesign d;
    int ok = glm_build_design("~ cov.csv", 4, &d);
    MU_ASSERT("reg: build should succeed", ok == 1);
    MU_ASSERT("reg: 2 columns", d.n_beta == 2);
    MU_ASSERT("reg: col1 label", strcmp(d.colnames[1], "cov") == 0);
    MU_ASSERT("reg: covariate values passed through",
              d.G[2][1] == 30.0 && d.G[3][1] == 40.0);
    glm_free_design(&d);
}

/* 4. interaction: '*' expands to main effects plus the product */
static void test_interaction(void)
{
    GlmDesign d;
    int ok = glm_build_design("~ grp.csv * cov.csv", 4, &d);
    MU_ASSERT("inter: build should succeed", ok == 1);
    MU_ASSERT("inter: 4 columns", d.n_beta == 4);
    MU_ASSERT("inter: labels in order",
              strcmp(d.colnames[0], "(Intercept)") == 0 &&
              strcmp(d.colnames[1], "grp[b]") == 0 &&
              strcmp(d.colnames[2], "cov") == 0 &&
              strcmp(d.colnames[3], "grp[b]:cov") == 0);
    MU_ASSERT("inter: product is grp[b]*cov",
              d.G[0][3] == 0.0 && d.G[3][3] == 40.0);
    glm_free_design(&d);
}

/* 5. factor() forces a numeric file to be categorical */
static void test_forced_factor(void)
{
    GlmDesign d;
    int ok = glm_build_design("~ factor(num.csv)", 4, &d);
    MU_ASSERT("factor: build should succeed", ok == 1);
    MU_ASSERT("factor: 2 columns", d.n_beta == 2);
    MU_ASSERT("factor: treated as factor, not covariate",
              strcmp(d.colnames[1], "num[2]") == 0);
    glm_free_design(&d);
}

/* 6. full three-way factorial gives 8 columns for 2x2x2 */
static void test_three_way(void)
{
    GlmDesign d;
    int ok = glm_build_design("~ fa.csv * fb.csv * fc.csv", 4, &d);
    MU_ASSERT("3way: build should succeed", ok == 1);
    MU_ASSERT("3way: 8 columns (2x2x2 full rank)", d.n_beta == 8);
    MU_ASSERT("3way: three-way interaction column last",
              strcmp(d.colnames[7], "fa[b]:fb[q]:fc[n]") == 0);
    glm_free_design(&d);
}

/* 7. wrong number of rows must fail cleanly */
static void test_row_mismatch(void)
{
    GlmDesign d;
    int ok = glm_build_design("~ bad.csv", 4, &d);
    MU_ASSERT("mismatch: build should fail", ok == 0);
}

/* 8. mixing '*' and ':' in one term is rejected */
static void test_mixed_operators(void)
{
    GlmDesign d;
    int ok = glm_build_design("~ grp.csv : cov.csv * fa.csv", 4, &d);
    MU_ASSERT("mixed: build should fail", ok == 0);
}

int main(void)
{
    setup_files();
    MU_RUN_TEST(test_group_comparison);
    MU_RUN_TEST(test_cell_means);
    MU_RUN_TEST(test_regression);
    MU_RUN_TEST(test_interaction);
    MU_RUN_TEST(test_forced_factor);
    MU_RUN_TEST(test_three_way);
    MU_RUN_TEST(test_row_mismatch);
    MU_RUN_TEST(test_mixed_operators);
    cleanup_files();
    printf("%d tests run, %d failed\n", tests_run, tests_failed);
    return tests_failed ? 1 : 0;
}
