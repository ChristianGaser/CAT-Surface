/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 *
 */

#ifndef _CAT_GLMFORMULA_H_
#define _CAT_GLMFORMULA_H_

#include <stdio.h>

/**
 * \brief A design matrix built from an R-style model formula.
 *
 * The observations (rows) correspond, in order, to the input scans passed
 * to the GLM tool.  The columns are the terms of the model formula after
 * contrast coding.  Column labels follow R conventions, e.g.
 * "(Intercept)", "age", "group[patient]" or "group[patient]:age".
 */
typedef struct {
    int      n_obs;      /* number of observations (== number of scans)   */
    int      n_beta;     /* number of design-matrix columns               */
    double **G;          /* n_obs x n_beta design matrix (bicpl ALLOC2D)  */
    char   **colnames;   /* n_beta column labels (heap strings)           */
} GlmDesign;

/**
 * \brief Build a design matrix from an R-style model formula.
 *
 * Parses a formula such as "~ group.csv * age.csv" and builds the
 * corresponding design matrix.  Every variable on the right-hand side is a
 * file with exactly \p n_obs rows (one value per scan); files whose values
 * are all numeric become covariates, otherwise they are treated as factors
 * (categorical).  A numeric file can be forced to a factor with factor(...).
 *
 * Supported operators: '+' (add term), ':' (interaction), '*' (main effects
 * plus interaction) and intercept control via a leading '0'/'-1' (drop) or
 * '1' (keep, the default).  A leading '~' is optional.
 *
 * Factor levels are sorted lexicographically; with an intercept a factor is
 * treatment-coded against its first (reference) level, otherwise a bare
 * main-effect factor is coded as full cell means.  Interaction operands are
 * always treatment-coded.
 *
 * \param formula (in)  the model formula string
 * \param n_obs   (in)  number of observations; each variable file must match
 * \param design  (out) design matrix, filled on success
 * \return 1 on success, 0 on error (with a message on stderr)
 */
int  glm_build_design(const char *formula, int n_obs, GlmDesign *design);

/**
 * \brief Release all memory held by a GlmDesign.
 *
 * \param design (in/out) design to free; may be partially built
 */
void glm_free_design(GlmDesign *design);

/**
 * \brief Print a human-readable listing of the design matrix.
 *
 * Writes the column labels followed by one row per observation.  Used both
 * for the glm.log record and for verification output.
 *
 * \param fp         (in) open stream to write to
 * \param design     (in) design matrix to print
 * \param row_labels (in) optional per-row labels (e.g. scan names); may be NULL
 */
void glm_print_design(FILE *fp, const GlmDesign *design, char **row_labels);

#endif
