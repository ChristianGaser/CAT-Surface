#include "minunit.h"
#include "CAT_Calc.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define EPS 1e-12

/* Evaluate an expression once over the given per-voxel input values. */
static double
eval(const char *expression, const double *values, int n_img)
{
    CAT_CalcExpr *expr;
    CAT_CalcCtx ctx;
    double *scratch = NULL, result;
    char err[256];
    int n_scratch;

    expr = CAT_CalcParse(expression, n_img, err, sizeof(err));
    if (expr == NULL) {
        printf("unexpected parse error in \"%s\": %s\n", expression, err);
        return NAN;
    }

    n_scratch = CAT_CalcScratchSize(expr);
    if (n_scratch > 0)
        scratch = (double *) malloc(sizeof(double) * (size_t) n_scratch);

    ctx.values = values;
    ctx.n_img = n_img;
    ctx.scratch = scratch;

    result = CAT_CalcEval(expr, &ctx);

    free(scratch);
    CAT_CalcFree(expr);

    return result;
}

/* Evaluate an expression that needs no image inputs. */
static double
eval0(const char *expression)
{
    double v[1] = {0.0};

    return eval(expression, v, 1);
}

static int
parse_fails(const char *expression, int n_img)
{
    char err[256];
    CAT_CalcExpr *expr = CAT_CalcParse(expression, n_img, err, sizeof(err));

    if (expr != NULL) {
        CAT_CalcFree(expr);
        return 0;
    }
    return err[0] != '\0';   /* a rejection must come with a diagnostic */
}

/* ------------------------------------------------------------------ */
/* arithmetic and precedence                                           */
/* ------------------------------------------------------------------ */

static void test_arithmetic(void)
{
    MU_ASSERT("addition and multiplication", fabs(eval0("1+2*3") - 7.0) < EPS);
    MU_ASSERT("parentheses override", fabs(eval0("(1+2)*3") - 9.0) < EPS);
    MU_ASSERT("subtraction is left associative",
              fabs(eval0("10-3-2") - 5.0) < EPS);
    MU_ASSERT("division is left associative",
              fabs(eval0("12/3/2") - 2.0) < EPS);
    MU_ASSERT("decimal literal", fabs(eval0("0.5+.25") - 0.75) < EPS);
    MU_ASSERT("exponent literal", fabs(eval0("2e3") - 2000.0) < EPS);
    MU_ASSERT("pi", fabs(eval0("pi") - 3.14159265358979) < 1e-12);
}

static void test_precedence(void)
{
    /* '^' binds tighter than unary minus, and is right associative */
    MU_ASSERT("-2^2 is -(2^2)", fabs(eval0("-2^2") + 4.0) < EPS);
    MU_ASSERT("2^3^2 is 2^(3^2)", fabs(eval0("2^3^2") - 512.0) < EPS);
    MU_ASSERT("negative exponent", fabs(eval0("2^-1") - 0.5) < EPS);
    MU_ASSERT("unary plus is a no-op", fabs(eval0("+3") - 3.0) < EPS);
    MU_ASSERT("double negation", fabs(eval0("--3") - 3.0) < EPS);
    /* comparison binds looser than arithmetic, logical looser still */
    MU_ASSERT("1+1==2 is true", fabs(eval0("1+1==2") - 1.0) < EPS);
    MU_ASSERT("and binds tighter than or",
              fabs(eval0("0 & 0 | 1") - 1.0) < EPS);
}

static void test_comparison_and_logic(void)
{
    MU_ASSERT("less than", fabs(eval0("1<2") - 1.0) < EPS);
    MU_ASSERT("greater or equal", fabs(eval0("2>=2") - 1.0) < EPS);
    MU_ASSERT("not equal, C spelling", fabs(eval0("1!=2") - 1.0) < EPS);
    MU_ASSERT("not equal, MATLAB spelling", fabs(eval0("1~=2") - 1.0) < EPS);
    MU_ASSERT("equality is false", fabs(eval0("1==2")) < EPS);
    MU_ASSERT("logical not", fabs(eval0("!0") - 1.0) < EPS);
    MU_ASSERT("logical not, MATLAB spelling", fabs(eval0("~5")) < EPS);
    MU_ASSERT("and", fabs(eval0("1&0")) < EPS);
    MU_ASSERT("doubled and", fabs(eval0("1&&1") - 1.0) < EPS);
    MU_ASSERT("doubled or", fabs(eval0("0||1") - 1.0) < EPS);
}

/* ------------------------------------------------------------------ */
/* image variables                                                     */
/* ------------------------------------------------------------------ */

static void test_image_variables(void)
{
    double v[3] = {2.0, 3.0, 10.0};

    MU_ASSERT("i1", fabs(eval("i1", v, 3) - 2.0) < EPS);
    MU_ASSERT("i3", fabs(eval("i3", v, 3) - 10.0) < EPS);
    MU_ASSERT("difference", fabs(eval("i1-i2", v, 3) + 1.0) < EPS);
    MU_ASSERT("mixed", fabs(eval("i1*i2+i3", v, 3) - 16.0) < EPS);
}

static void test_matlab_operator_spellings(void)
{
    double v[2] = {3.0, 4.0};

    /* everything here is element-wise already, so .* ./ .^ are synonyms */
    MU_ASSERT("dot-star", fabs(eval("i1.*i2", v, 2) - 12.0) < EPS);
    MU_ASSERT("dot-slash", fabs(eval("i2./i1", v, 2) - 4.0 / 3.0) < EPS);
    MU_ASSERT("dot-caret", fabs(eval("i1.^2", v, 2) - 9.0) < EPS);
}

static void test_masking_idiom(void)
{
    double keep[2] = {1.0, 42.0};
    double drop[2] = {0.0, 42.0};

    MU_ASSERT("mask keeps the voxel",
              fabs(eval("(i1>0.5).*i2", keep, 2) - 42.0) < EPS);
    MU_ASSERT("mask drops the voxel",
              fabs(eval("(i1>0.5).*i2", drop, 2)) < EPS);
}

/* ------------------------------------------------------------------ */
/* functions                                                           */
/* ------------------------------------------------------------------ */

static void test_functions(void)
{
    double v[2] = {-4.0, 9.0};

    MU_ASSERT("abs", fabs(eval("abs(i1)", v, 2) - 4.0) < EPS);
    MU_ASSERT("sqrt", fabs(eval("sqrt(i2)", v, 2) - 3.0) < EPS);
    MU_ASSERT("exp/log round trip", fabs(eval0("log(exp(2))") - 2.0) < EPS);
    MU_ASSERT("log2", fabs(eval0("log2(8)") - 3.0) < EPS);
    MU_ASSERT("log10", fabs(eval0("log10(1000)") - 3.0) < EPS);
    MU_ASSERT("floor", fabs(eval0("floor(2.7)") - 2.0) < EPS);
    MU_ASSERT("ceil", fabs(eval0("ceil(2.1)") - 3.0) < EPS);
    MU_ASSERT("round", fabs(eval0("round(2.5)") - 3.0) < EPS);
    MU_ASSERT("sign of a negative", fabs(eval0("sign(-7)") + 1.0) < EPS);
    MU_ASSERT("sign of zero", fabs(eval0("sign(0)")) < EPS);
    MU_ASSERT("sign of NaN is NaN", isnan(eval0("sign(nan)")));
    MU_ASSERT("isnan detects nan", fabs(eval0("isnan(nan)") - 1.0) < EPS);
    MU_ASSERT("isinf detects inf", fabs(eval0("isinf(inf)") - 1.0) < EPS);
    MU_ASSERT("isfinite rejects inf", fabs(eval0("isfinite(inf)")) < EPS);
    MU_ASSERT("min", fabs(eval("min(i1,i2)", v, 2) + 4.0) < EPS);
    MU_ASSERT("max", fabs(eval("max(i1,i2)", v, 2) - 9.0) < EPS);
    MU_ASSERT("pow", fabs(eval0("pow(2,10)") - 1024.0) < EPS);
    MU_ASSERT("atan2", fabs(eval0("atan2(0,-1)") - 3.14159265358979) < 1e-12);
    /* mod() takes the sign of the divisor, as in MATLAB */
    MU_ASSERT("mod of a negative", fabs(eval0("mod(-1,3)") - 2.0) < EPS);
    MU_ASSERT("mod of a positive", fabs(eval0("mod(7,3)") - 1.0) < EPS);
    MU_ASSERT("nested calls", fabs(eval0("max(min(5,3),2)") - 3.0) < EPS);
}

/* ------------------------------------------------------------------ */
/* matrix mode                                                         */
/* ------------------------------------------------------------------ */

static void test_matrix_reductions(void)
{
    double odd[3] = {1.0, 2.0, 3.0};
    double even[4] = {4.0, 1.0, 3.0, 2.0};

    MU_ASSERT("mean across images", fabs(eval("mean(X)", odd, 3) - 2.0) < EPS);
    MU_ASSERT("median of an odd count",
              fabs(eval("median(X)", odd, 3) - 2.0) < EPS);
    /* the median must not depend on the order the files were given in */
    MU_ASSERT("median of an even count averages the middle pair",
              fabs(eval("median(X)", even, 4) - 2.5) < EPS);
    MU_ASSERT("sample standard deviation",
              fabs(eval("std(X)", odd, 3) - 1.0) < EPS);
    MU_ASSERT("mean of an even count", fabs(eval("mean(X)", even, 4) - 2.5) < EPS);
}

static void test_matrix_broadcasting(void)
{
    double v[3] = {2.0, 4.0, 6.0};

    /* X is an ordinary value, so a reduction can take a computed argument */
    MU_ASSERT("reduction of a shifted vector",
              fabs(eval("std(X-i1)", v, 3) - 2.0) < EPS);
    MU_ASSERT("scalar broadcasts over the vector",
              fabs(eval("mean(X*2)", v, 3) - 8.0) < EPS);
    MU_ASSERT("fraction above a threshold",
              fabs(eval("mean(X>3)", v, 3) - 2.0 / 3.0) < EPS);
    MU_ASSERT("vector against vector",
              fabs(eval("mean(X-X)", v, 3)) < EPS);
    MU_ASSERT("reduction combines with scalars",
              fabs(eval("i1/mean(X)", v, 3) - 0.5) < EPS);
    MU_ASSERT("unary minus over a vector",
              fabs(eval("mean(-X)", v, 3) + 4.0) < EPS);
    MU_ASSERT("element-wise function over a vector",
              fabs(eval("mean(abs(X-4))", v, 3) - 4.0 / 3.0) < EPS);
    /* a bare vector at the root collapses to its first element */
    MU_ASSERT("bare X yields the first image", fabs(eval("X", v, 3) - 2.0) < EPS);
}

static void test_reductions_skip_nonfinite(void)
{
    double with_nan[3] = {1.0, NAN, 3.0};
    double one_finite[2] = {5.0, NAN};
    double none[2] = {NAN, NAN};

    /* one subject with a NaN must not blank the group result at that voxel */
    MU_ASSERT("mean skips NaN", fabs(eval("mean(X)", with_nan, 3) - 2.0) < EPS);
    MU_ASSERT("median skips NaN",
              fabs(eval("median(X)", with_nan, 3) - 2.0) < EPS);
    MU_ASSERT("std skips NaN",
              fabs(eval("std(X)", with_nan, 3) - sqrt(2.0)) < EPS);
    MU_ASSERT("mean of a single finite value",
              fabs(eval("mean(X)", one_finite, 2) - 5.0) < EPS);
    MU_ASSERT("std of a single finite value is undefined",
              isnan(eval("std(X)", one_finite, 2)));
    MU_ASSERT("mean of nothing finite is NaN",
              isnan(eval("mean(X)", none, 2)));
    MU_ASSERT("median of nothing finite is NaN",
              isnan(eval("median(X)", none, 2)));
}

static void test_matrix_introspection(void)
{
    CAT_CalcExpr *plain, *matrix;
    char err[256];

    plain = CAT_CalcParse("i1+i2", 2, err, sizeof(err));
    matrix = CAT_CalcParse("mean(X)", 2, err, sizeof(err));

    MU_ASSERT("plain expression parses", plain != NULL);
    MU_ASSERT("matrix expression parses", matrix != NULL);
    MU_ASSERT("X is not reported for a plain expression",
              CAT_CalcUsesMatrix(plain) == 0);
    MU_ASSERT("X is reported for a matrix expression",
              CAT_CalcUsesMatrix(matrix) == 1);
    MU_ASSERT("scratch covers every slot a vector could land in",
              CAT_CalcScratchSize(matrix) >= 2);

    CAT_CalcFree(plain);
    CAT_CalcFree(matrix);
}

/* ------------------------------------------------------------------ */
/* whole-volume application                                            */
/* ------------------------------------------------------------------ */

static void test_apply_volume(void)
{
    double a[4] = {1.0, 2.0, 3.0, 4.0};
    double b[4] = {10.0, 20.0, 30.0, 40.0};
    double *images[2];
    double out[4];
    CAT_CalcExpr *expr;
    char err[256];
    int ok = 1, i;

    images[0] = a;
    images[1] = b;

    expr = CAT_CalcParse("i1+i2", 2, err, sizeof(err));
    MU_ASSERT("volume expression parses", expr != NULL);
    MU_ASSERT("apply succeeds", CAT_CalcApply(expr, images, 2, 4, out) == 0);

    for (i = 0; i < 4; i++) {
        if (fabs(out[i] - (a[i] + b[i])) > EPS)
            ok = 0;
    }
    MU_ASSERT("every voxel is the sum of its inputs", ok);
    CAT_CalcFree(expr);

    /* the same driver, this time reducing across images */
    expr = CAT_CalcParse("mean(X)", 2, err, sizeof(err));
    MU_ASSERT("matrix volume expression parses", expr != NULL);
    MU_ASSERT("matrix apply succeeds",
              CAT_CalcApply(expr, images, 2, 4, out) == 0);

    ok = 1;
    for (i = 0; i < 4; i++) {
        if (fabs(out[i] - 0.5 * (a[i] + b[i])) > EPS)
            ok = 0;
    }
    MU_ASSERT("every voxel is the mean of its inputs", ok);

    /* a count mismatch must be refused, not read out of bounds */
    MU_ASSERT("apply rejects the wrong image count",
              CAT_CalcApply(expr, images, 1, 4, out) != 0);
    CAT_CalcFree(expr);
}

/* ------------------------------------------------------------------ */
/* diagnostics                                                         */
/* ------------------------------------------------------------------ */

static void test_parse_errors(void)
{
    MU_ASSERT("empty expression", parse_fails("", 2));
    MU_ASSERT("whitespace only", parse_fails("   ", 2));
    MU_ASSERT("unknown name", parse_fails("foo(i1)", 2));
    MU_ASSERT("unbalanced parenthesis", parse_fails("(i1+i2", 2));
    MU_ASSERT("missing operand", parse_fails("i1+", 2));
    MU_ASSERT("trailing input", parse_fails("i1 i2", 2));
    MU_ASSERT("assignment is not comparison", parse_fails("i1=i2", 2));
    MU_ASSERT("stray dot", parse_fails("i1.i2", 2));
    MU_ASSERT("unexpected character", parse_fails("i1 # i2", 2));
    MU_ASSERT("too few arguments", parse_fails("min(i1)", 2));
    MU_ASSERT("too many arguments", parse_fails("abs(i1,i2)", 2));
    MU_ASSERT("missing parentheses after a function", parse_fails("sqrt i1", 2));

    /* referring past the end of the file list is the likeliest mistake */
    MU_ASSERT("image index above the file count", parse_fails("i3", 2));
    MU_ASSERT("image index of zero", parse_fails("i0", 2));
    MU_ASSERT("no images at all", parse_fails("i1", 0));

    MU_ASSERT("a valid expression still parses", !parse_fails("i1+i2", 2));
}

static void test_error_message_is_useful(void)
{
    char err[256];
    CAT_CalcExpr *expr = CAT_CalcParse("i1 + i9", 2, err, sizeof(err));

    MU_ASSERT("out-of-range image is rejected", expr == NULL);
    MU_ASSERT("the message names the offending variable",
              strstr(err, "i9") != NULL);
    MU_ASSERT("the message locates the error",
              strstr(err, "character") != NULL);

    CAT_CalcFree(expr);
}

static void test_null_arguments(void)
{
    char err[256];

    MU_ASSERT("a NULL expression is refused",
              CAT_CalcParse(NULL, 2, err, sizeof(err)) == NULL);
    MU_ASSERT("freeing NULL is safe", (CAT_CalcFree(NULL), 1));
    MU_ASSERT("scratch size of NULL is zero", CAT_CalcScratchSize(NULL) == 0);
    MU_ASSERT("matrix use of NULL is zero", CAT_CalcUsesMatrix(NULL) == 0);
    MU_ASSERT("evaluating NULL yields NaN", isnan(CAT_CalcEval(NULL, NULL)));
}

int main(void)
{
    MU_RUN_TEST(test_arithmetic);
    MU_RUN_TEST(test_precedence);
    MU_RUN_TEST(test_comparison_and_logic);
    MU_RUN_TEST(test_image_variables);
    MU_RUN_TEST(test_matlab_operator_spellings);
    MU_RUN_TEST(test_masking_idiom);
    MU_RUN_TEST(test_functions);
    MU_RUN_TEST(test_matrix_reductions);
    MU_RUN_TEST(test_matrix_broadcasting);
    MU_RUN_TEST(test_reductions_skip_nonfinite);
    MU_RUN_TEST(test_matrix_introspection);
    MU_RUN_TEST(test_apply_volume);
    MU_RUN_TEST(test_parse_errors);
    MU_RUN_TEST(test_error_message_is_useful);
    MU_RUN_TEST(test_null_arguments);
    printf("%d tests run, %d failed\n", tests_run, tests_failed);
    return tests_failed ? 1 : 0;
}
