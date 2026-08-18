/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 *
 */

#ifndef _CAT_CALC_H_
#define _CAT_CALC_H_

#include <stddef.h>

/**
 * \brief A parsed voxel-wise image expression.
 *
 * The expression is parsed once into a small abstract syntax tree and then
 * evaluated at every voxel, so the textual formula is never re-scanned inside
 * the voxel loop.  The structure is opaque; build it with CAT_CalcParse() and
 * release it with CAT_CalcFree().
 */
typedef struct CAT_CalcExpr CAT_CalcExpr;

/**
 * \brief Per-voxel evaluation context.
 *
 * Holds everything the evaluator needs that changes from voxel to voxel, so
 * that the parsed expression itself stays read-only.  One context per thread
 * is therefore enough to evaluate the same expression concurrently.
 */
typedef struct {
    const double *values;   /* n_img values: i1..iN at the current voxel   */
    int           n_img;    /* number of input images                      */
    double       *scratch;  /* CAT_CalcScratchSize() doubles of workspace  */
} CAT_CalcCtx;

/**
 * \brief Parse a voxel-wise image expression.
 *
 * Compiles an spm_imcalc-style formula into an expression tree.  The grammar
 * is a standard precedence-climbing arithmetic grammar:
 *
 *   - operands: numeric literals, the image variables \c i1 ... \c iN, the
 *     matrix variable \c X, and the constants \c pi, \c inf and \c nan;
 *   - operators, weakest binding first: \c | \c || , \c & \c && ,
 *     \c == \c ~= \c != , \c < \c <= \c > \c >= , \c + \c - , \c * \c / ,
 *     unary \c - \c + \c ! \c ~ , and \c ^ (right associative);
 *   - element-wise spellings \c .* \c ./ \c .^ are accepted as synonyms of
 *     \c * \c / \c ^ so that expressions can be copied over from MATLAB
 *     unchanged -- every operator here is element-wise already;
 *   - element-wise functions: \c abs, \c sqrt, \c exp, \c log, \c log2,
 *     \c log10, \c sin, \c cos, \c tan, \c asin, \c acos, \c atan, \c floor,
 *     \c ceil, \c round, \c sign, \c isnan, \c isinf, \c isfinite, and the
 *     two-argument \c min, \c max, \c pow, \c atan2, \c mod;
 *   - reductions across images: \c mean, \c median, \c std.
 *
 * \c X is the matrix variable: at each voxel it stands for the vector of all
 * \p n_img input values, so \c mean(X) is the voxel-wise mean across images
 * and \c std(X-i1) reduces an expression that was itself evaluated per image.
 * Comparisons and logical operators yield 1.0 or 0.0, which makes
 * \c "(i1>0.5).*i2" the idiomatic way to mask.
 *
 * Reductions ignore non-finite values, matching get_mean_double() and
 * get_std_double() in CAT_Math; a voxel at which every input is non-finite
 * reduces to NaN, as does \c std of fewer than two finite values.  \c std is
 * the sample standard deviation (normalised by n-1).  In practice the
 * non-finite values reaching a reduction are the ones the formula itself
 * produced: nifti_read_buffer() in the vendored nifticlib already rewrites
 * every non-finite value in a file to zero as it loads, so a NaN on disk
 * arrives as a zero and is counted as one.
 *
 * \param expression (in)  the formula to parse
 * \param n_img      (in)  number of input images; \c i1..i<n_img> are valid
 * \param err        (out) buffer receiving a diagnostic on failure; may be NULL
 * \param err_len    (in)  size of \p err in bytes
 * \return the parsed expression, or NULL on error (with a message in \p err)
 */
CAT_CalcExpr *CAT_CalcParse(const char *expression, int n_img,
                            char *err, int err_len);

/**
 * \brief Release a parsed expression.
 *
 * \param expr (in/out) expression to free; NULL is ignored
 */
void CAT_CalcFree(CAT_CalcExpr *expr);

/**
 * \brief Number of doubles of scratch space the evaluator needs.
 *
 * Intermediate results of an expression that mentions \c X are vectors of
 * length n_img, and each node that can produce one owns a slice of the
 * caller-supplied scratch buffer.  Allocate this many doubles and point
 * CAT_CalcCtx::scratch at them.
 *
 * \param expr (in) parsed expression
 * \return required scratch length in doubles (may be 0)
 */
int CAT_CalcScratchSize(const CAT_CalcExpr *expr);

/**
 * \brief Report whether the expression uses the matrix variable \c X.
 *
 * \param expr (in) parsed expression
 * \return 1 if \c X occurs in the expression, 0 otherwise
 */
int CAT_CalcUsesMatrix(const CAT_CalcExpr *expr);

/**
 * \brief Evaluate a parsed expression at a single voxel.
 *
 * \p ctx supplies the input values for the voxel and the scratch buffer; the
 * expression itself is not modified, so separate contexts may evaluate the
 * same expression in parallel.  A result that is still a vector at the root
 * (as in the expression \c "X") is reduced to its first element, since a
 * single output image can hold only one value per voxel.
 *
 * \param expr (in)     parsed expression
 * \param ctx  (in/out) evaluation context; scratch is overwritten
 * \return the value of the expression at this voxel
 */
double CAT_CalcEval(const CAT_CalcExpr *expr, CAT_CalcCtx *ctx);

/**
 * \brief Evaluate a parsed expression over a whole volume.
 *
 * Gathers the \p n_img inputs at each voxel, evaluates \p expr and stores the
 * result in \p out.  Scratch space is allocated internally.
 *
 * \param expr   (in)  parsed expression
 * \param images (in)  \p n_img volumes of \p nvox voxels each
 * \param n_img  (in)  number of input volumes; must match the parse-time count
 * \param nvox   (in)  number of voxels per volume
 * \param out    (out) double[nvox] result volume
 * \return 0 on success, -1 on error (with a message on stderr)
 */
int CAT_CalcApply(const CAT_CalcExpr *expr, double *const *images,
                  int n_img, size_t nvox, double *out);

#endif
