/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 *
 */

/*
 * Minimal R-style model-formula parser and design-matrix builder.
 *
 * Every variable on the right-hand side of the formula names a plain text
 * file holding one value per scan (row).  Numeric files become covariates,
 * text files become factors.  Supported operators are '+', ':' and '*',
 * plus intercept control via a leading '0'/'-1' or '1'.  This is a
 * deliberately small subset of R's formula language; see CAT_GlmFormula.h.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "CAT_GlmFormula.h"

#define GLM_MAXVARS   64
#define GLM_MAXTERMS  512
#define GLM_MAXCOLS   1024
#define GLM_NAME      256
#define GLM_PATH      1024

/* ----------------------------------------------------------------------- */
/*  internal data structures */
/* ----------------------------------------------------------------------- */

typedef struct {
    char    name[GLM_NAME];   /* short label derived from the filename    */
    char    filename[GLM_PATH];
    int     forced_factor;    /* set via factor(...)                      */
    int     is_factor;
    int     n_obs;
    double *num;              /* covariate values (is_factor == 0)        */
    int     n_levels;         /* factor levels (is_factor == 1)           */
    char  **levels;           /* sorted level labels                      */
    int    *code;             /* per-obs level index into levels[]        */
} Var;

typedef struct {
    int n;                    /* number of variables in this term         */
    int v[GLM_MAXVARS];       /* indices into ctx->vars                    */
} Term;

typedef struct {
    int  intercept;
    int  n_vars;
    Var  vars[GLM_MAXVARS];
    int  n_terms;
    Term terms[GLM_MAXTERMS];
} Ctx;

/* A single expanded design column while it is being assembled. */
typedef struct {
    double *val;
    char    lab[GLM_NAME];
} Col;

/* ----------------------------------------------------------------------- */
/*  small helpers */
/* ----------------------------------------------------------------------- */

/* Derive a short label from a path: drop directory and everything from the
 * first '.', e.g. "data/site.nii.gz" -> "site". */
static void
base_label(const char *path, char *out)
{
    const char *slash = strrchr(path, '/');
    const char *base = slash ? slash + 1 : path;
    const char *dot = strchr(base, '.');
    size_t len = dot ? (size_t)(dot - base) : strlen(base);

    if (len >= GLM_NAME)
        len = GLM_NAME - 1;
    memcpy(out, base, len);
    out[len] = '\0';
}

static int
cmp_str(const void *a, const void *b)
{
    return strcmp(*(const char *const *)a, *(const char *const *)b);
}

/* Is the whole token a finite number? */
static int
is_number(const char *tok)
{
    char *end;
    if (tok[0] == '\0')
        return 0;
    (void) strtod(tok, &end);
    return (*end == '\0');
}

/* ----------------------------------------------------------------------- */
/*  read a variable file into a Var (covariate or factor) */
/* ----------------------------------------------------------------------- */
static int
read_var(Var *var, int n_obs)
{
    FILE  *fp;
    char   tok[GLM_NAME];
    char **raw;
    int    i, l, count = 0, extra;

    fp = fopen(var->filename, "r");
    if (fp == NULL) {
        fprintf(stderr, "glm_build_design: cannot open file %s\n",
                var->filename);
        return 0;
    }

    raw = (char **) malloc(sizeof(char *) * n_obs);
    while (count < n_obs && fscanf(fp, "%255s", tok) == 1)
        raw[count++] = strdup(tok);
    extra = (fscanf(fp, "%255s", tok) == 1);   /* is there a surplus row? */
    fclose(fp);

    if (count != n_obs || extra) {
        fprintf(stderr,
            "glm_build_design: file %s has %d%s values, expected %d\n",
            var->filename, count, extra ? "+" : "", n_obs);
        for (i = 0; i < count; i++) free(raw[i]);
        free(raw);
        return 0;
    }

    var->n_obs = n_obs;

    /* covariate unless forced to a factor or any value is non-numeric */
    var->is_factor = var->forced_factor;
    if (!var->is_factor) {
        for (i = 0; i < n_obs; i++) {
            if (!is_number(raw[i])) { var->is_factor = 1; break; }
        }
    }

    if (!var->is_factor) {
        var->num = (double *) malloc(sizeof(double) * n_obs);
        for (i = 0; i < n_obs; i++)
            var->num[i] = strtod(raw[i], NULL);
    } else {
        /* collect sorted unique levels */
        char **uniq = (char **) malloc(sizeof(char *) * n_obs);
        int    nuniq = 0;
        for (i = 0; i < n_obs; i++) {
            int seen = 0;
            for (l = 0; l < nuniq; l++)
                if (strcmp(uniq[l], raw[i]) == 0) { seen = 1; break; }
            if (!seen)
                uniq[nuniq++] = raw[i];   /* borrow pointer from raw */
        }
        qsort(uniq, nuniq, sizeof(char *), cmp_str);

        var->n_levels = nuniq;
        var->levels = (char **) malloc(sizeof(char *) * nuniq);
        for (l = 0; l < nuniq; l++)
            var->levels[l] = strdup(uniq[l]);

        var->code = (int *) malloc(sizeof(int) * n_obs);
        for (i = 0; i < n_obs; i++) {
            for (l = 0; l < nuniq; l++)
                if (strcmp(var->levels[l], raw[i]) == 0) {
                    var->code[i] = l; break;
                }
        }
        free(uniq);
    }

    for (i = 0; i < count; i++) free(raw[i]);
    free(raw);
    return 1;
}

/* Find or register a variable by filename (+ forced flag); read it once. */
static int
get_var(Ctx *ctx, const char *filename, int forced, int n_obs)
{
    int i;
    for (i = 0; i < ctx->n_vars; i++) {
        if (strcmp(ctx->vars[i].filename, filename) == 0 &&
            ctx->vars[i].forced_factor == forced)
            return i;
    }
    if (ctx->n_vars >= GLM_MAXVARS) {
        fprintf(stderr, "glm_build_design: too many variables\n");
        return -1;
    }
    Var *v = &ctx->vars[ctx->n_vars];
    memset(v, 0, sizeof(*v));
    strncpy(v->filename, filename, GLM_PATH - 1);
    base_label(filename, v->name);
    v->forced_factor = forced;
    if (!read_var(v, n_obs))
        return -1;
    return ctx->n_vars++;
}

/* ----------------------------------------------------------------------- */
/*  tokenizer */
/* ----------------------------------------------------------------------- */

enum { T_VAR, T_PLUS, T_STAR, T_COLON, T_LP, T_RP, T_INT0, T_INT1, T_END };

typedef struct { int type; char text[GLM_PATH]; } Token;

/* Split the right-hand side of the formula into tokens. */
static int
tokenize(const char *rhs, Token *toks, int max)
{
    int n = 0;
    const char *p = rhs;

    while (*p) {
        if (isspace((unsigned char) *p)) { p++; continue; }
        if (n >= max - 1) {
            fprintf(stderr, "glm_build_design: formula too long\n");
            return -1;
        }
        if (*p == '+') { toks[n++].type = T_PLUS;  p++; continue; }
        if (*p == '*') { toks[n++].type = T_STAR;  p++; continue; }
        if (*p == ':') { toks[n++].type = T_COLON; p++; continue; }
        if (*p == '(') { toks[n++].type = T_LP;    p++; continue; }
        if (*p == ')') { toks[n++].type = T_RP;    p++; continue; }

        /* read a word up to the next space or operator */
        {
            char *w = toks[n].text;
            int   len = 0;
            while (*p && !isspace((unsigned char) *p) &&
                   *p != '+' && *p != '*' && *p != ':' &&
                   *p != '(' && *p != ')') {
                if (len < GLM_PATH - 1) w[len++] = *p;
                p++;
            }
            w[len] = '\0';

            if (strcmp(w, "0") == 0 || strcmp(w, "-1") == 0)
                toks[n++].type = T_INT0;
            else if (strcmp(w, "1") == 0)
                toks[n++].type = T_INT1;
            else
                toks[n++].type = T_VAR;
        }
    }
    toks[n].type = T_END;
    return n;
}

/* ----------------------------------------------------------------------- */
/*  parse tokens into an intercept flag and a list of terms */
/* ----------------------------------------------------------------------- */

/* Add every non-empty subset of the m atoms as a term (for '*' expansion),
 * ordered by increasing interaction degree. */
static void
add_subset_terms(Ctx *ctx, const int *atoms, int m)
{
    int size, mask;
    for (size = 1; size <= m; size++) {
        for (mask = 1; mask < (1 << m); mask++) {
            int bits = 0, i;
            for (i = 0; i < m; i++) if (mask & (1 << i)) bits++;
            if (bits != size) continue;
            if (ctx->n_terms >= GLM_MAXTERMS) return;
            Term *t = &ctx->terms[ctx->n_terms++];
            t->n = 0;
            for (i = 0; i < m; i++)
                if (mask & (1 << i)) t->v[t->n++] = atoms[i];
        }
    }
}

static int
parse_formula(Ctx *ctx, const char *formula, int n_obs)
{
    Token toks[GLM_MAXTERMS];
    int   ntok, i;

    /* skip an optional leading "~" */
    const char *rhs = strchr(formula, '~');
    rhs = rhs ? rhs + 1 : formula;

    ntok = tokenize(rhs, toks, GLM_MAXTERMS);
    if (ntok < 0)
        return 0;

    ctx->intercept = 1;

    /* walk chunks separated by '+' */
    i = 0;
    while (i < ntok) {
        int atoms[GLM_MAXVARS], m = 0;
        int has_star = 0, has_colon = 0;

        /* collect one chunk */
        while (i < ntok && toks[i].type != T_PLUS) {
            switch (toks[i].type) {
            case T_INT0: ctx->intercept = 0; i++; break;
            case T_INT1: ctx->intercept = 1; i++; break;
            case T_STAR: has_star = 1; i++; break;
            case T_COLON: has_colon = 1; i++; break;
            case T_VAR: {
                int forced = 0;
                const char *fname = toks[i].text;
                /* factor(<file>) */
                if (strcmp(toks[i].text, "factor") == 0 &&
                    i + 3 < ntok + 1 &&
                    toks[i+1].type == T_LP && toks[i+2].type == T_VAR &&
                    toks[i+3].type == T_RP) {
                    forced = 1;
                    fname = toks[i+2].text;
                    i += 4;
                } else {
                    i++;
                }
                if (m >= GLM_MAXVARS) {
                    fprintf(stderr, "glm_build_design: term too large\n");
                    return 0;
                }
                atoms[m] = get_var(ctx, fname, forced, n_obs);
                if (atoms[m] < 0)
                    return 0;
                m++;
                break;
            }
            default: i++; break;
            }
        }
        if (i < ntok && toks[i].type == T_PLUS) i++;  /* consume '+' */

        if (m == 0)
            continue;                                 /* e.g. bare "0"/"1" */
        if (has_star && has_colon) {
            fprintf(stderr, "glm_build_design: mixing '*' and ':' in one "
                            "term is not supported\n");
            return 0;
        }
        if (m == 1) {
            Term *t = &ctx->terms[ctx->n_terms++];
            t->n = 1; t->v[0] = atoms[0];
        } else if (has_star) {
            add_subset_terms(ctx, atoms, m);
        } else {                                      /* ':' or juxtaposition */
            Term *t = &ctx->terms[ctx->n_terms++];
            t->n = m;
            memcpy(t->v, atoms, sizeof(int) * m);
        }
    }
    return 1;
}

/* ----------------------------------------------------------------------- */
/*  build the design columns from the parsed terms */
/* ----------------------------------------------------------------------- */

/* Expand a single variable into its contrast columns for a given term.
 * Returns the number of columns written to out[] (each a Col with an
 * allocated val vector). */
static int
expand_var(const Var *var, int treatment, int n_obs, Col *out)
{
    int l, i, start, nc = 0;

    if (!var->is_factor) {
        out[0].val = (double *) malloc(sizeof(double) * n_obs);
        for (i = 0; i < n_obs; i++) out[0].val[i] = var->num[i];
        strncpy(out[0].lab, var->name, GLM_NAME - 1);
        return 1;
    }

    /* factor: treatment contrasts drop the reference (first) level */
    start = treatment ? 1 : 0;
    for (l = start; l < var->n_levels; l++) {
        out[nc].val = (double *) malloc(sizeof(double) * n_obs);
        for (i = 0; i < n_obs; i++)
            out[nc].val[i] = (var->code[i] == l) ? 1.0 : 0.0;
        snprintf(out[nc].lab, GLM_NAME, "%s[%s]", var->name, var->levels[l]);
        nc++;
    }
    return nc;
}

static int
build_design(Ctx *ctx, GlmDesign *design)
{
    Col cols[GLM_MAXCOLS];
    int ncol = 0, n_obs = -1, i, j, t;

    /* every variable shares the same number of observations */
    if (ctx->n_vars > 0)
        n_obs = ctx->vars[0].n_obs;

    if (ctx->intercept) {
        cols[ncol].val = (double *) malloc(sizeof(double) * design->n_obs);
        for (i = 0; i < design->n_obs; i++) cols[ncol].val[i] = 1.0;
        strcpy(cols[ncol].lab, "(Intercept)");
        ncol++;
    }
    n_obs = design->n_obs;

    for (t = 0; t < ctx->n_terms; t++) {
        Term *term = &ctx->terms[t];
        /* accumulator of the cartesian product, seeded with a "ones" column */
        Col acc[GLM_MAXCOLS];
        int nacc = 1;
        acc[0].val = (double *) malloc(sizeof(double) * n_obs);
        for (i = 0; i < n_obs; i++) acc[0].val[i] = 1.0;
        acc[0].lab[0] = '\0';

        for (j = 0; j < term->n; j++) {
            Var *var = &ctx->vars[term->v[j]];
            /* full cell-means only for a bare main-effect factor w/o intercept */
            int treatment = !(term->n == 1 && !ctx->intercept);
            Col set[GLM_MAXCOLS];
            int nset = expand_var(var, treatment, n_obs, set);

            Col next[GLM_MAXCOLS];
            int nnext = 0, a, s;
            for (a = 0; a < nacc; a++) {
                for (s = 0; s < nset; s++) {
                    next[nnext].val = (double *) malloc(sizeof(double) * n_obs);
                    for (i = 0; i < n_obs; i++)
                        next[nnext].val[i] = acc[a].val[i] * set[s].val[i];
                    if (acc[a].lab[0] == '\0')
                        snprintf(next[nnext].lab, GLM_NAME, "%s", set[s].lab);
                    else
                        snprintf(next[nnext].lab, GLM_NAME, "%s:%s",
                                 acc[a].lab, set[s].lab);
                    nnext++;
                }
            }
            for (a = 0; a < nacc; a++) free(acc[a].val);
            for (s = 0; s < nset; s++) free(set[s].val);
            memcpy(acc, next, sizeof(Col) * nnext);
            nacc = nnext;
        }

        for (i = 0; i < nacc; i++) {
            if (ncol >= GLM_MAXCOLS) {
                fprintf(stderr, "glm_build_design: too many columns\n");
                return 0;
            }
            cols[ncol++] = acc[i];   /* transfer ownership */
        }
    }

    /* materialise into the GlmDesign (row-major G[obs][col]) */
    design->n_beta = ncol;
    design->G = (double **) malloc(sizeof(double *) * n_obs);
    for (i = 0; i < n_obs; i++) {
        design->G[i] = (double *) malloc(sizeof(double) * ncol);
        for (j = 0; j < ncol; j++)
            design->G[i][j] = cols[j].val[i];
    }
    design->colnames = (char **) malloc(sizeof(char *) * ncol);
    for (j = 0; j < ncol; j++) {
        design->colnames[j] = strdup(cols[j].lab);
        free(cols[j].val);
    }
    return 1;
}

/* ----------------------------------------------------------------------- */
/*  public API */
/* ----------------------------------------------------------------------- */
static void
free_ctx(Ctx *ctx)
{
    int i, l;
    for (i = 0; i < ctx->n_vars; i++) {
        Var *v = &ctx->vars[i];
        free(v->num);
        free(v->code);
        if (v->levels) {
            for (l = 0; l < v->n_levels; l++) free(v->levels[l]);
            free(v->levels);
        }
    }
}

int
glm_build_design(const char *formula, int n_obs, GlmDesign *design)
{
    Ctx ctx;
    int ok;

    memset(&ctx, 0, sizeof(ctx));
    memset(design, 0, sizeof(*design));
    design->n_obs = n_obs;

    if (!parse_formula(&ctx, formula, n_obs)) {
        free_ctx(&ctx);
        return 0;
    }
    ok = build_design(&ctx, design);
    free_ctx(&ctx);
    if (!ok) {
        glm_free_design(design);
        return 0;
    }
    return 1;
}

void
glm_free_design(GlmDesign *design)
{
    int i;
    if (design->G) {
        for (i = 0; i < design->n_obs; i++) free(design->G[i]);
        free(design->G);
        design->G = NULL;
    }
    if (design->colnames) {
        for (i = 0; i < design->n_beta; i++) free(design->colnames[i]);
        free(design->colnames);
        design->colnames = NULL;
    }
}

void
glm_print_design(FILE *fp, const GlmDesign *design, char **row_labels)
{
    int i, j;

    fprintf(fp, "%-20s", "");
    for (j = 0; j < design->n_beta; j++)
        fprintf(fp, " %12s", design->colnames[j]);
    fprintf(fp, "\n");

    for (i = 0; i < design->n_obs; i++) {
        fprintf(fp, "%-20s", row_labels ? row_labels[i] : "");
        for (j = 0; j < design->n_beta; j++)
            fprintf(fp, " %12.4f", design->G[i][j]);
        fprintf(fp, "\n");
    }
}
