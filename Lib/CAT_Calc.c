/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 *
 */

/*
 * Voxel-wise image calculator, in the spirit of SPM's spm_imcalc.
 *
 * The formula is parsed once into a small tree and then evaluated at every
 * voxel.  Two things make the tree worth having over a naive re-scan of the
 * string: the per-voxel cost drops to one traversal, and the matrix variable
 * X falls out for free.  X stands for the whole vector of input values at the
 * current voxel, so a node evaluates either to a scalar or to a vector of
 * length n_img, and the ordinary operators broadcast between the two exactly
 * as MATLAB would.  That is what makes mean(X), median(X) and std(X) work on
 * an arbitrary number of inputs, and what lets a reduction take a computed
 * argument such as std(X - i1) rather than only the bare variable.
 *
 * Every node that can produce a vector owns a slice of a caller-supplied
 * scratch buffer.  Because slots are handed out per node, the slices used by
 * the two operands of a binary node are always disjoint from each other and
 * from the node's own slice, so no intermediate can be clobbered while it is
 * still live.  Keeping the scratch outside the tree also keeps the tree
 * read-only during evaluation, which is what allows one parsed expression to
 * be shared by several threads, each with its own context.
 */

#include "CAT_Calc.h"

#include <ctype.h>
#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* ------------------------------------------------------------------ */
/* operator and function codes                                         */
/* ------------------------------------------------------------------ */

enum {
    OP_ADD = 1, OP_SUB, OP_MUL, OP_DIV, OP_POW,
    OP_LT, OP_LE, OP_GT, OP_GE, OP_EQ, OP_NE,
    OP_AND, OP_OR, OP_NEG, OP_NOT
};

enum {
    FN_ABS = 1, FN_SQRT, FN_EXP, FN_LOG, FN_LOG2, FN_LOG10,
    FN_SIN, FN_COS, FN_TAN, FN_ASIN, FN_ACOS, FN_ATAN,
    FN_FLOOR, FN_CEIL, FN_ROUND, FN_SIGN,
    FN_ISNAN, FN_ISINF, FN_ISFINITE,
    FN_MIN, FN_MAX, FN_POW, FN_ATAN2, FN_MOD,
    FN_MEAN, FN_MEDIAN, FN_STD
};

typedef struct {
    const char *name;
    int         id;
    int         nargs;   /* number of arguments, 1 or 2                  */
    int         reduce;  /* 1 if the result is a scalar reduction of a
                            vector argument, 0 if the function is applied
                            element by element                           */
} FuncDef;

static const FuncDef functions[] = {
    {"abs",      FN_ABS,      1, 0},
    {"sqrt",     FN_SQRT,     1, 0},
    {"exp",      FN_EXP,      1, 0},
    {"log",      FN_LOG,      1, 0},
    {"log2",     FN_LOG2,     1, 0},
    {"log10",    FN_LOG10,    1, 0},
    {"sin",      FN_SIN,      1, 0},
    {"cos",      FN_COS,      1, 0},
    {"tan",      FN_TAN,      1, 0},
    {"asin",     FN_ASIN,     1, 0},
    {"acos",     FN_ACOS,     1, 0},
    {"atan",     FN_ATAN,     1, 0},
    {"floor",    FN_FLOOR,    1, 0},
    {"ceil",     FN_CEIL,     1, 0},
    {"round",    FN_ROUND,    1, 0},
    {"sign",     FN_SIGN,     1, 0},
    {"isnan",    FN_ISNAN,    1, 0},
    {"isinf",    FN_ISINF,    1, 0},
    {"isfinite", FN_ISFINITE, 1, 0},
    {"min",      FN_MIN,      2, 0},
    {"max",      FN_MAX,      2, 0},
    {"pow",      FN_POW,      2, 0},
    {"atan2",    FN_ATAN2,    2, 0},
    {"mod",      FN_MOD,      2, 0},
    {"mean",     FN_MEAN,     1, 1},
    {"median",   FN_MEDIAN,   1, 1},
    {"std",      FN_STD,      1, 1},
    {NULL,       0,           0, 0}
};

/* ------------------------------------------------------------------ */
/* expression tree                                                     */
/* ------------------------------------------------------------------ */

typedef enum {
    ND_NUM,      /* numeric literal or named constant                   */
    ND_IMG,      /* i<k>: the k-th input image at the current voxel      */
    ND_MTX,      /* X: the whole vector of inputs at the current voxel   */
    ND_UNARY,
    ND_BINARY,
    ND_FUNC
} NodeType;

typedef struct Node {
    NodeType     type;
    int          op;      /* operator code, or function id for ND_FUNC   */
    int          reduce;  /* ND_FUNC: function reduces a vector          */
    double       num;     /* ND_NUM                                      */
    int          idx;     /* ND_IMG: zero-based image index              */
    int          slot;    /* scratch slot, or -1 if the node has none    */
    struct Node *a, *b;
} Node;

struct CAT_CalcExpr {
    Node *root;
    int   n_slots;
    int   n_img;
    int   uses_matrix;
};

/* ------------------------------------------------------------------ */
/* lexer                                                               */
/* ------------------------------------------------------------------ */

enum { T_END = 0, T_NUM, T_IDENT, T_LPAR, T_RPAR, T_COMMA, T_OP };

typedef struct {
    const char *begin;      /* start of the expression, for error offsets */
    const char *p;          /* current scan position                      */
    const char *tok;        /* start of the current token                 */
    int         type;       /* token type                                 */
    double      num;        /* T_NUM                                      */
    char        ident[64];  /* T_IDENT                                    */
    int         op;         /* T_OP                                       */
    int         n_img;      /* number of images i1..iN may refer to       */
    int         n_slots;    /* scratch slots handed out so far            */
    int         uses_matrix;
    int         failed;
    char       *err;
    int         err_len;
} Parser;

static void
parse_error(Parser *ps, const char *fmt, ...)
{
    va_list ap;
    char msg[192];

    if (ps->failed)
        return;     /* keep the first, most specific diagnostic */
    ps->failed = 1;

    va_start(ap, fmt);
    vsnprintf(msg, sizeof(msg), fmt, ap);
    va_end(ap);

    if (ps->err != NULL && ps->err_len > 0)
        snprintf(ps->err, (size_t) ps->err_len, "%s (at character %d)",
                 msg, (int) (ps->tok - ps->begin) + 1);
}

static void
next_token(Parser *ps)
{
    const char *p = ps->p;

    while (isspace((unsigned char) *p))
        p++;

    ps->tok = p;

    if (*p == '\0') {
        ps->type = T_END;
        ps->p = p;
        return;
    }

    if (isdigit((unsigned char) *p) ||
        (*p == '.' && isdigit((unsigned char) p[1]))) {
        char *end;

        ps->num = strtod(p, &end);
        ps->p = end;
        ps->type = T_NUM;
        return;
    }

    if (isalpha((unsigned char) *p) || *p == '_') {
        size_t n = 0;

        while (isalnum((unsigned char) *p) || *p == '_') {
            if (n < sizeof(ps->ident) - 1)
                ps->ident[n++] = *p;
            p++;
        }
        ps->ident[n] = '\0';
        ps->p = p;
        ps->type = T_IDENT;
        return;
    }

    ps->type = T_OP;
    ps->p = p + 1;

    switch (*p) {
    case '(':
        ps->type = T_LPAR;
        return;
    case ')':
        ps->type = T_RPAR;
        return;
    case ',':
        ps->type = T_COMMA;
        return;
    case '+':
        ps->op = OP_ADD;
        return;
    case '-':
        ps->op = OP_SUB;
        return;
    case '*':
        ps->op = OP_MUL;
        return;
    case '/':
        ps->op = OP_DIV;
        return;
    case '^':
        ps->op = OP_POW;
        return;
    case '<':
        if (p[1] == '=') {
            ps->op = OP_LE;
            ps->p = p + 2;
        } else
            ps->op = OP_LT;
        return;
    case '>':
        if (p[1] == '=') {
            ps->op = OP_GE;
            ps->p = p + 2;
        } else
            ps->op = OP_GT;
        return;
    case '=':
        if (p[1] == '=') {
            ps->op = OP_EQ;
            ps->p = p + 2;
            return;
        }
        parse_error(ps, "'=' is not an operator, use '==' to compare");
        ps->type = T_END;
        return;
    case '!':
    case '~':
        if (p[1] == '=') {
            ps->op = OP_NE;
            ps->p = p + 2;
        } else
            ps->op = OP_NOT;
        return;
    case '&':
        if (p[1] == '&')
            ps->p = p + 2;
        ps->op = OP_AND;
        return;
    case '|':
        if (p[1] == '|')
            ps->p = p + 2;
        ps->op = OP_OR;
        return;
    case '.':
        /* MATLAB spells the element-wise operators with a leading dot; every
           operator here is element-wise anyway, so accept both spellings and
           let expressions be copied over from spm_imcalc unchanged. */
        if (p[1] == '*' || p[1] == '/' || p[1] == '^') {
            ps->op = (p[1] == '*') ? OP_MUL : (p[1] == '/') ? OP_DIV : OP_POW;
            ps->p = p + 2;
            return;
        }
        parse_error(ps, "'.' must be followed by '*', '/' or '^'");
        ps->type = T_END;
        return;
    default:
        parse_error(ps, "unexpected character '%c'", *p);
        ps->type = T_END;
        return;
    }
}

/* ------------------------------------------------------------------ */
/* tree construction                                                   */
/* ------------------------------------------------------------------ */

static void
free_node(Node *node)
{
    if (node == NULL)
        return;
    free_node(node->a);
    free_node(node->b);
    free(node);
}

static Node *
new_node(Parser *ps, NodeType type)
{
    Node *node = (Node *) calloc(1, sizeof(Node));

    if (node == NULL) {
        parse_error(ps, "out of memory");
        return NULL;
    }
    node->type = type;

    /* only nodes that compute something need somewhere to put a vector */
    node->slot = (type == ND_UNARY || type == ND_BINARY || type == ND_FUNC)
             ? ps->n_slots++ : -1;

    return node;
}

static const FuncDef *
lookup_function(const char *name)
{
    int i;

    for (i = 0; functions[i].name != NULL; i++) {
        if (strcmp(functions[i].name, name) == 0)
            return &functions[i];
    }
    return NULL;
}

static Node *parse_or(Parser *ps);

static Node *
parse_call(Parser *ps, const FuncDef *fn)
{
    Node *node;

    next_token(ps);
    if (ps->type != T_LPAR) {
        parse_error(ps, "expected '(' after '%s'", fn->name);
        return NULL;
    }
    next_token(ps);

    node = new_node(ps, ND_FUNC);
    if (node == NULL)
        return NULL;
    node->op = fn->id;
    node->reduce = fn->reduce;

    node->a = parse_or(ps);

    if (!ps->failed && fn->nargs == 2) {
        if (ps->type != T_COMMA)
            parse_error(ps, "'%s' takes two arguments", fn->name);
        else {
            next_token(ps);
            node->b = parse_or(ps);
        }
    }

    if (!ps->failed && ps->type == T_COMMA)
        parse_error(ps, "too many arguments for '%s'", fn->name);
    if (!ps->failed && ps->type != T_RPAR)
        parse_error(ps, "expected ')' after the arguments of '%s'", fn->name);

    if (ps->failed) {
        free_node(node);
        return NULL;
    }
    next_token(ps);

    return node;
}

static Node *
parse_ident(Parser *ps)
{
    char name[sizeof(ps->ident)];
    const FuncDef *fn;
    Node *node;

    strcpy(name, ps->ident);

    /* the matrix variable: the whole vector of inputs at this voxel */
    if (strcmp(name, "X") == 0) {
        next_token(ps);
        node = new_node(ps, ND_MTX);
        if (node != NULL)
            ps->uses_matrix = 1;
        return node;
    }

    if (strcmp(name, "pi") == 0 || strcmp(name, "inf") == 0 ||
        strcmp(name, "nan") == 0) {
        double value = (name[0] == 'p') ? 4.0 * atan(1.0)
                     : (name[0] == 'i') ? (double) INFINITY : (double) NAN;

        next_token(ps);
        node = new_node(ps, ND_NUM);
        if (node != NULL)
            node->num = value;
        return node;
    }

    /* image variables i1 .. iN */
    if (name[0] == 'i' && isdigit((unsigned char) name[1])) {
        int k = 0, j;

        for (j = 1; name[j] != '\0'; j++) {
            if (!isdigit((unsigned char) name[j])) {
                parse_error(ps, "malformed image variable '%s'", name);
                return NULL;
            }
            if (k < 1000000)
                k = k * 10 + (name[j] - '0');
        }
        if (k < 1 || k > ps->n_img) {
            parse_error(ps, "'%s' refers to image %d, but only %d image%s given",
                        name, k, ps->n_img, ps->n_img == 1 ? " was" : "s were");
            return NULL;
        }
        next_token(ps);
        node = new_node(ps, ND_IMG);
        if (node != NULL)
            node->idx = k - 1;
        return node;
    }

    fn = lookup_function(name);
    if (fn != NULL)
        return parse_call(ps, fn);

    parse_error(ps, "unknown name '%s'", name);

    return NULL;
}

static Node *
parse_primary(Parser *ps)
{
    Node *node;

    if (ps->failed)
        return NULL;

    if (ps->type == T_NUM) {
        node = new_node(ps, ND_NUM);
        if (node == NULL)
            return NULL;
        node->num = ps->num;
        next_token(ps);
        return node;
    }

    if (ps->type == T_LPAR) {
        next_token(ps);
        node = parse_or(ps);
        if (!ps->failed && ps->type != T_RPAR)
            parse_error(ps, "expected ')'");
        if (ps->failed) {
            free_node(node);
            return NULL;
        }
        next_token(ps);
        return node;
    }

    if (ps->type == T_IDENT)
        return parse_ident(ps);

    if (ps->type == T_END)
        parse_error(ps, "unexpected end of expression");
    else
        parse_error(ps, "expected a value");

    return NULL;
}

static Node *parse_unary(Parser *ps);

/* '^' binds tighter than unary minus, so -2^2 is -(2^2), and it is right
   associative, so 2^3^2 is 2^(3^2). */
static Node *
parse_power(Parser *ps)
{
    Node *base = parse_primary(ps);
    Node *node;

    if (ps->failed || ps->type != T_OP || ps->op != OP_POW) {
        if (ps->failed) {
            free_node(base);
            return NULL;
        }
        return base;
    }

    next_token(ps);
    node = new_node(ps, ND_BINARY);
    if (node == NULL) {
        free_node(base);
        return NULL;
    }
    node->op = OP_POW;
    node->a = base;
    node->b = parse_unary(ps);

    if (ps->failed) {
        free_node(node);
        return NULL;
    }
    return node;
}

static Node *
parse_unary(Parser *ps)
{
    Node *node;
    int op;

    if (ps->failed)
        return NULL;

    if (ps->type != T_OP ||
        (ps->op != OP_SUB && ps->op != OP_ADD && ps->op != OP_NOT))
        return parse_power(ps);

    op = ps->op;
    next_token(ps);

    if (op == OP_ADD)
        return parse_unary(ps);     /* unary plus changes nothing */

    node = new_node(ps, ND_UNARY);
    if (node == NULL)
        return NULL;
    node->op = (op == OP_SUB) ? OP_NEG : OP_NOT;
    node->a = parse_unary(ps);

    if (ps->failed) {
        free_node(node);
        return NULL;
    }
    return node;
}

/*
 * One level of left-associative binary operators.  Each level accepts the
 * operators in \p ops and otherwise defers to the next tighter-binding level.
 */
static Node *
parse_binary_level(Parser *ps, const int *ops, int n_ops,
                   Node *(*next_level)(Parser *))
{
    Node *lhs = next_level(ps);

    for (;;) {
        Node *node;
        int i, op = 0;

        if (ps->failed || ps->type != T_OP)
            break;

        for (i = 0; i < n_ops; i++) {
            if (ps->op == ops[i]) {
                op = ops[i];
                break;
            }
        }
        if (op == 0)
            break;

        next_token(ps);
        node = new_node(ps, ND_BINARY);
        if (node == NULL) {
            free_node(lhs);
            return NULL;
        }
        node->op = op;
        node->a = lhs;
        node->b = next_level(ps);
        lhs = node;
    }

    if (ps->failed) {
        free_node(lhs);
        return NULL;
    }
    return lhs;
}

static Node *
parse_multiplicative(Parser *ps)
{
    static const int ops[] = {OP_MUL, OP_DIV};

    return parse_binary_level(ps, ops, 2, parse_unary);
}

static Node *
parse_additive(Parser *ps)
{
    static const int ops[] = {OP_ADD, OP_SUB};

    return parse_binary_level(ps, ops, 2, parse_multiplicative);
}

static Node *
parse_relational(Parser *ps)
{
    static const int ops[] = {OP_LT, OP_LE, OP_GT, OP_GE};

    return parse_binary_level(ps, ops, 4, parse_additive);
}

static Node *
parse_equality(Parser *ps)
{
    static const int ops[] = {OP_EQ, OP_NE};

    return parse_binary_level(ps, ops, 2, parse_relational);
}

static Node *
parse_and(Parser *ps)
{
    static const int ops[] = {OP_AND};

    return parse_binary_level(ps, ops, 1, parse_equality);
}

static Node *
parse_or(Parser *ps)
{
    static const int ops[] = {OP_OR};

    return parse_binary_level(ps, ops, 1, parse_and);
}

/* ------------------------------------------------------------------ */
/* evaluation                                                          */
/* ------------------------------------------------------------------ */

/*
 * A value is either a scalar or a vector of n_img elements.  Vectors arise
 * only from X, so an expression that never mentions X evaluates entirely in
 * scalars.
 */
typedef struct {
    double        s;
    const double *v;
    int           is_vec;
} Val;

#define VAL_AT(val, k) ((val).is_vec ? (val).v[k] : (val).s)

static Val
scalar_val(double s)
{
    Val out;

    out.s = s;
    out.v = NULL;
    out.is_vec = 0;

    return out;
}

static Val
vector_val(const double *v)
{
    Val out;

    out.s = 0.0;
    out.v = v;
    out.is_vec = 1;

    return out;
}

static double *
slot_of(const Node *node, CAT_CalcCtx *ctx)
{
    return ctx->scratch + (size_t) node->slot * (size_t) ctx->n_img;
}

static double
apply_unary(int op, double x)
{
    switch (op) {
    case OP_NEG: return -x;
    case OP_NOT: return (x == 0.0) ? 1.0 : 0.0;
    }
    return NAN;
}

static double
apply_binary(int op, double x, double y)
{
    switch (op) {
    case OP_ADD: return x + y;
    case OP_SUB: return x - y;
    case OP_MUL: return x * y;
    case OP_DIV: return x / y;
    case OP_POW: return pow(x, y);
    case OP_LT:  return (x <  y) ? 1.0 : 0.0;
    case OP_LE:  return (x <= y) ? 1.0 : 0.0;
    case OP_GT:  return (x >  y) ? 1.0 : 0.0;
    case OP_GE:  return (x >= y) ? 1.0 : 0.0;
    case OP_EQ:  return (x == y) ? 1.0 : 0.0;
    case OP_NE:  return (x != y) ? 1.0 : 0.0;
    case OP_AND: return (x != 0.0 && y != 0.0) ? 1.0 : 0.0;
    case OP_OR:  return (x != 0.0 || y != 0.0) ? 1.0 : 0.0;
    }
    return NAN;
}

static double
apply_func1(int id, double x)
{
    switch (id) {
    case FN_ABS:      return fabs(x);
    case FN_SQRT:     return sqrt(x);
    case FN_EXP:      return exp(x);
    case FN_LOG:      return log(x);
    case FN_LOG2:     return log2(x);
    case FN_LOG10:    return log10(x);
    case FN_SIN:      return sin(x);
    case FN_COS:      return cos(x);
    case FN_TAN:      return tan(x);
    case FN_ASIN:     return asin(x);
    case FN_ACOS:     return acos(x);
    case FN_ATAN:     return atan(x);
    case FN_FLOOR:    return floor(x);
    case FN_CEIL:     return ceil(x);
    case FN_ROUND:    return round(x);
    case FN_SIGN:     return isnan(x) ? NAN
                           : (x > 0.0) ? 1.0 : (x < 0.0) ? -1.0 : 0.0;
    case FN_ISNAN:    return isnan(x) ? 1.0 : 0.0;
    case FN_ISINF:    return isinf(x) ? 1.0 : 0.0;
    case FN_ISFINITE: return isfinite(x) ? 1.0 : 0.0;
    }
    return NAN;
}

static double
apply_func2(int id, double x, double y)
{
    double r;

    switch (id) {
    case FN_MIN:   return fmin(x, y);
    case FN_MAX:   return fmax(x, y);
    case FN_POW:   return pow(x, y);
    case FN_ATAN2: return atan2(x, y);
    case FN_MOD:
        /* MATLAB's mod(), i.e. the result takes the sign of the divisor */
        if (y == 0.0)
            return x;
        r = fmod(x, y);
        if (r != 0.0 && ((r < 0.0) != (y < 0.0)))
            r += y;
        return r;
    }
    return NAN;
}

static int
compare_double(const void *a, const void *b)
{
    double x = *(const double *) a;
    double y = *(const double *) b;

    return (x > y) - (x < y);
}

/*
 * Reduce a vector to a scalar.  Non-finite inputs are skipped rather than
 * propagated, which matches get_mean_double() and get_std_double() in
 * CAT_Math and is what multi-subject averaging wants: one subject with a NaN
 * at a voxel should not blank the whole group result there.  \p work must
 * hold at least \p n doubles; the median sorts its filtered copy in place.
 */
static double
reduce_values(int id, const double *v, int n, double *work)
{
    double sum = 0.0, mean, var = 0.0;
    int i, m = 0;

    for (i = 0; i < n; i++) {
        if (isfinite(v[i]))
            work[m++] = v[i];
    }
    if (m == 0)
        return NAN;

    if (id == FN_MEDIAN) {
        qsort(work, (size_t) m, sizeof(double), compare_double);
        if (m % 2 != 0)
            return work[m / 2];
        return 0.5 * (work[m / 2 - 1] + work[m / 2]);
    }

    for (i = 0; i < m; i++)
        sum += work[i];
    mean = sum / (double) m;

    if (id == FN_MEAN)
        return mean;

    /* FN_STD: sample standard deviation, normalised by n-1 as in
       get_std_double(); undefined for a single observation */
    if (m < 2)
        return NAN;
    for (i = 0; i < m; i++)
        var += (work[i] - mean) * (work[i] - mean);

    return sqrt(var / (double) (m - 1));
}

static Val eval_node(const Node *node, CAT_CalcCtx *ctx);

static Val
eval_func(const Node *node, CAT_CalcCtx *ctx)
{
    Val x = eval_node(node->a, ctx);
    double *dst = slot_of(node, ctx);
    Val y;
    int k;

    if (node->reduce) {
        if (x.is_vec)
            return scalar_val(reduce_values(node->op, x.v, ctx->n_img, dst));
        return scalar_val(reduce_values(node->op, &x.s, 1, dst));
    }

    if (node->b == NULL) {
        if (!x.is_vec)
            return scalar_val(apply_func1(node->op, x.s));
        for (k = 0; k < ctx->n_img; k++)
            dst[k] = apply_func1(node->op, x.v[k]);
        return vector_val(dst);
    }

    y = eval_node(node->b, ctx);
    if (!x.is_vec && !y.is_vec)
        return scalar_val(apply_func2(node->op, x.s, y.s));
    for (k = 0; k < ctx->n_img; k++)
        dst[k] = apply_func2(node->op, VAL_AT(x, k), VAL_AT(y, k));

    return vector_val(dst);
}

static Val
eval_node(const Node *node, CAT_CalcCtx *ctx)
{
    Val x, y;
    double *dst;
    int k;

    switch (node->type) {
    case ND_NUM:
        return scalar_val(node->num);

    case ND_IMG:
        return scalar_val(ctx->values[node->idx]);

    case ND_MTX:
        return vector_val(ctx->values);

    case ND_UNARY:
        x = eval_node(node->a, ctx);
        if (!x.is_vec)
            return scalar_val(apply_unary(node->op, x.s));
        dst = slot_of(node, ctx);
        for (k = 0; k < ctx->n_img; k++)
            dst[k] = apply_unary(node->op, x.v[k]);
        return vector_val(dst);

    case ND_BINARY:
        x = eval_node(node->a, ctx);
        y = eval_node(node->b, ctx);
        if (!x.is_vec && !y.is_vec)
            return scalar_val(apply_binary(node->op, x.s, y.s));
        dst = slot_of(node, ctx);
        for (k = 0; k < ctx->n_img; k++)
            dst[k] = apply_binary(node->op, VAL_AT(x, k), VAL_AT(y, k));
        return vector_val(dst);

    case ND_FUNC:
        return eval_func(node, ctx);
    }

    return scalar_val(NAN);
}

/* ------------------------------------------------------------------ */
/* public interface                                                    */
/* ------------------------------------------------------------------ */

/**
 * \brief Parse a voxel-wise image expression.
 *
 * Compiles an spm_imcalc-style formula into an expression tree; see the
 * declaration in CAT_Calc.h for the full grammar, the meaning of the
 * matrix variable \c X and the NaN conventions of the reductions.
 *
 * \param expression (in)  the formula to parse
 * \param n_img      (in)  number of input images; \c i1..i<n_img> are valid
 * \param err        (out) buffer receiving a diagnostic on failure; may be NULL
 * \param err_len    (in)  size of \p err in bytes
 * \return the parsed expression, or NULL on error (with a message in \p err)
 */
CAT_CalcExpr *
CAT_CalcParse(const char *expression, int n_img, char *err, int err_len)
{
    CAT_CalcExpr *expr;
    Parser ps;
    Node *root;

    if (err != NULL && err_len > 0)
        err[0] = '\0';

    memset(&ps, 0, sizeof(ps));
    ps.begin = (expression != NULL) ? expression : "";
    ps.p = ps.begin;
    ps.tok = ps.begin;
    ps.n_img = n_img;
    ps.err = err;
    ps.err_len = err_len;

    if (expression == NULL) {
        parse_error(&ps, "no expression given");
        return NULL;
    }
    if (n_img < 1) {
        parse_error(&ps, "at least one input image is required");
        return NULL;
    }

    next_token(&ps);
    if (!ps.failed && ps.type == T_END)
        parse_error(&ps, "the expression is empty");

    root = parse_or(&ps);

    if (!ps.failed && ps.type != T_END)
        parse_error(&ps, "unexpected trailing input");

    if (ps.failed) {
        free_node(root);
        return NULL;
    }

    expr = (CAT_CalcExpr *) calloc(1, sizeof(*expr));
    if (expr == NULL) {
        parse_error(&ps, "out of memory");
        free_node(root);
        return NULL;
    }
    expr->root = root;
    expr->n_slots = ps.n_slots;
    expr->n_img = n_img;
    expr->uses_matrix = ps.uses_matrix;

    return expr;
}

/**
 * \brief Release a parsed expression.
 *
 * \param expr (in/out) expression to free; NULL is ignored
 */
void
CAT_CalcFree(CAT_CalcExpr *expr)
{
    if (expr == NULL)
        return;
    free_node(expr->root);
    free(expr);
}

/**
 * \brief Number of doubles of scratch space the evaluator needs.
 *
 * One slice of n_img doubles per node that can produce a vector.
 *
 * \param expr (in) parsed expression
 * \return required scratch length in doubles (may be 0)
 */
int
CAT_CalcScratchSize(const CAT_CalcExpr *expr)
{
    if (expr == NULL)
        return 0;

    return expr->n_slots * expr->n_img;
}

/**
 * \brief Report whether the expression uses the matrix variable \c X.
 *
 * \param expr (in) parsed expression
 * \return 1 if \c X occurs in the expression, 0 otherwise
 */
int
CAT_CalcUsesMatrix(const CAT_CalcExpr *expr)
{
    return (expr != NULL) ? expr->uses_matrix : 0;
}

/**
 * \brief Evaluate a parsed expression at a single voxel.
 *
 * The expression is not modified, so separate contexts may evaluate the same
 * expression concurrently.  A result that is still a vector at the root (as
 * in the expression \c "X") is reduced to its first element, since a single
 * output image can hold only one value per voxel.
 *
 * \param expr (in)     parsed expression
 * \param ctx  (in/out) evaluation context; scratch is overwritten
 * \return the value of the expression at this voxel
 */
double
CAT_CalcEval(const CAT_CalcExpr *expr, CAT_CalcCtx *ctx)
{
    Val v;

    if (expr == NULL || ctx == NULL || ctx->values == NULL)
        return NAN;

    v = eval_node(expr->root, ctx);

    return v.is_vec ? v.v[0] : v.s;
}

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
int
CAT_CalcApply(const CAT_CalcExpr *expr, double *const *images,
              int n_img, size_t nvox, double *out)
{
    double *values = NULL, *scratch = NULL;
    CAT_CalcCtx ctx;
    size_t j;
    int i, n_scratch;

    if (expr == NULL || images == NULL || out == NULL) {
        fprintf(stderr, "CAT_CalcApply: invalid argument.\n");
        return -1;
    }
    if (n_img != expr->n_img) {
        fprintf(stderr, "CAT_CalcApply: expression was parsed for %d image%s "
                "but %d were supplied.\n",
                expr->n_img, expr->n_img == 1 ? "" : "s", n_img);
        return -1;
    }

    n_scratch = CAT_CalcScratchSize(expr);

    values = (double *) malloc(sizeof(double) * (size_t) n_img);
    if (n_scratch > 0)
        scratch = (double *) malloc(sizeof(double) * (size_t) n_scratch);

    if (values == NULL || (n_scratch > 0 && scratch == NULL)) {
        fprintf(stderr, "CAT_CalcApply: memory allocation error.\n");
        free(values);
        free(scratch);
        return -1;
    }

    ctx.values = values;
    ctx.n_img = n_img;
    ctx.scratch = scratch;

    for (j = 0; j < nvox; j++) {
        for (i = 0; i < n_img; i++)
            values[i] = images[i][j];
        out[j] = CAT_CalcEval(expr, &ctx);
    }

    free(values);
    free(scratch);

    return 0;
}
