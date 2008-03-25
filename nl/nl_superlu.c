/*
 *  OpenNL: Numerical Library
 *  Copyright (C) 2004 Bruno Levy
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 *
 *  If you modify this software, you should include a notice giving the
 *  name of the person performing the modification, the date of modification,
 *  and the reason for such modification.
 *
 *  Contact: Bruno Levy
 *
 *     levy@loria.fr
 *
 *     ISA Project
 *     LORIA, INRIA Lorraine, 
 *     Campus Scientifique, BP 239
 *     54506 VANDOEUVRE LES NANCY CEDEX 
 *     FRANCE
 *
 *  Note that the GNU General Public License does not permit incorporating
 *  the Software into proprietary programs. 
 */

#include <nl/nl_superlu.h>
#include <nl/nl_context.h>

/************************************************************************/
/* SuperLU wrapper */

#ifdef NL_USE_SUPERLU

/* SuperLU includes */
#include <dsp_defs.h>
#include <util.h>

/* Note: SuperLU is difficult to call, but it is worth it.    */
/* Here is a driver inspired by A. Sheffer's "cow flattener". */
NLboolean nlSolve_SUPERLU() {

    /* OpenNL Context */
    NLSparseMatrix* M  = &(nlCurrentContext->M) ;
    NLdouble* b          = nlCurrentContext->b ;
    NLdouble* x          = nlCurrentContext->x ;

    /* Compressed Row Storage matrix representation */
    NLuint    n      = nlCurrentContext->n ;
    NLuint    nnz    = nlSparseMatrixNNZ(M) ; /* Number of Non-Zero coeffs */
    NLint*    xa     = NL_NEW_ARRAY(NLint, n+1) ;
    NLdouble* rhs    = NL_NEW_ARRAY(NLdouble, n) ;
    NLdouble* a      = NL_NEW_ARRAY(NLdouble, nnz) ;
    NLint*    asub   = NL_NEW_ARRAY(NLint, nnz) ;

    /* Permutation vector */
    NLint*    perm_r  = NL_NEW_ARRAY(NLint, n) ;
    NLint*    perm    = NL_NEW_ARRAY(NLint, n) ;

    /* SuperLU variables */
    SuperMatrix A, B ; /* System       */
    SuperMatrix L, U ; /* Inverse of A */
    NLint info ;       /* status code  */
    DNformat *vals = NULL ; /* access to result */
    double *rvals  = NULL ; /* access to result */

    /* SuperLU options and stats */
    superlu_options_t options ;
    SuperLUStat_t     stat ;

    /* Temporary variables */
    NLRowColumn* Ri = NULL ;
    NLuint         i,jj,count ;
    
    /* Sanity checks */
    nl_assert(!(M->storage & NL_MATRIX_STORE_SYMMETRIC)) ;
    nl_assert(M->storage & NL_MATRIX_STORE_ROWS) ;
    nl_assert(M->m == M->n) ;
    
    /*
     * Step 1: convert matrix M into SuperLU compressed column 
     *   representation.
     * -------------------------------------------------------
     */

    count = 0 ;
    for(i=0; i<n; i++) {
        Ri = &(M->row[i]) ;
        xa[i] = count ;
        for(jj=0; jj<Ri->size; jj++) {
            a[count]    = Ri->coeff[jj].value ;
            asub[count] = Ri->coeff[jj].index ;
            count++ ;
        }
    }
    xa[n] = nnz ;

    /* Save memory for SuperLU */
    nlSparseMatrixClear(M) ;


    /*
     * Rem: SuperLU does not support symmetric storage.
     * In fact, for symmetric matrix, what we need 
     * is a SuperLLt algorithm (SuperNodal sparse Cholesky),
     * but it does not exist, anybody wants to implement it ?
     * However, this is not a big problem (SuperLU is just
     * a superset of what we really need.
     */
    dCreate_CompCol_Matrix(
        &A, n, n, nnz, a, asub, xa, 
        SLU_NR,              /* Row_wise, no supernode */
        SLU_D,               /* doubles                */ 
        SLU_GE               /* general storage        */
    );

    /* Step 2: create vector */
    dCreate_Dense_Matrix(
        &B, n, 1, b, n, 
        SLU_DN, /* Fortran-type column-wise storage */
        SLU_D,  /* doubles                          */
        SLU_GE  /* general                          */
    );
            

    /* Step 3: set SuperLU options 
     * ------------------------------
     */

    set_default_options(&options) ;

    switch(nlCurrentContext->solver) {
    case NL_SUPERLU_EXT: {
        options.ColPerm = NATURAL ;
    } break ;
    case NL_PERM_SUPERLU_EXT: {
        options.ColPerm = COLAMD ;
    } break ;
    case NL_SYMMETRIC_SUPERLU_EXT: {
        options.ColPerm = MMD_AT_PLUS_A ;
        options.SymmetricMode = YES ;
    } break ;
    default: {
        nl_assert_not_reached ;
    } break ;
    }

    StatInit(&stat) ;

    /* Step 4: call SuperLU main routine
     * ---------------------------------
     */

    dgssv(&options, &A, perm, perm_r, &L, &U, &B, &stat, &info);

    /* Step 5: get the solution
     * ------------------------
     * Fortran-type column-wise storage
     */
    vals = (DNformat*)B.Store;
    rvals = (double*)(vals->nzval);
    if(info == 0) {
        for(i = 0; i <  n; i++){
            x[i] = rvals[i];
        }
    } else {
        nlError("nlSolve()", "SuperLU failed") ;
    }

    /* Step 6: cleanup
     * ---------------
     */

    /*
     *  For these two ones, only the "store" structure
     * needs to be deallocated (the arrays have been allocated
     * by us).
     */
    Destroy_SuperMatrix_Store(&A) ;
    Destroy_SuperMatrix_Store(&B) ;

    /*
     *   These ones need to be fully deallocated (they have been
     * allocated by SuperLU).
     */
    Destroy_SuperNode_Matrix(&L);
    Destroy_CompCol_Matrix(&U);

    /* There are some dynamically allocated vectors in the stats */
    StatFree(&stat) ;

    NL_DELETE_ARRAY(xa) ;
    NL_DELETE_ARRAY(rhs) ;
    NL_DELETE_ARRAY(a) ;
    NL_DELETE_ARRAY(asub) ;
    NL_DELETE_ARRAY(perm_r) ;
    NL_DELETE_ARRAY(perm) ;

    return (info == 0) ;
}

#else

NLboolean nlSolve_SUPERLU() {
    nl_assert_not_reached ;
}

#endif

