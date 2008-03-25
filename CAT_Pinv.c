/* Christian Gaser - christian.gaser@uni-jena.de                             */
/* Department of Psychiatry                                                  */
/* University of Jena                                                        */
/*                                                                           */
/* Copyright Christian Gaser, University of Jena.                            */

#include <volume_io/internal_volume_io.h>

/* produces a matrix Ainv of the same dimensions as A', so that A*Ainv*A = Ainv and A*Ainv and Ainv*A are Hermitian. The computation is based and SVD(A) and any singular values less than a tolerance of 1e-10 are treated as zero. The rank of the matrix A is returned.
*/

#define TOL 1e-10
#define EPS 1e-15

public int pinv(
    int    m,
    int    n,
    Real   **A,
    Real   **Ainv)
{
    int         i, j, k, r;
    Real        **U, **V, **S, *W, **Ut;
    
    ALLOC2D(U, m, n);
    ALLOC2D(V, n, n);
    ALLOC2D(S,  n, n);
    ALLOC2D(Ut, n, m);
    ALLOC(W, n);

    /* copy matrix A to U, because this matrix will be overwritten by singular_value_decomposition */
    for_less(i, 0, m)
        for_less(j, 0, n)
            U[i][j] = A[i][j];

    (void) singular_value_decomposition(m, n, U, W, V );
    
    /* rank of matrix using a somewhat arbitrary value of 1e-10 as tolerance for singular values */
    r = 0;
    for_less(i, 0, n)
        if( W[i] > TOL )
	    r += 1;

    if( r == 0) {
        for_less(i, 0, n)
	    for_less(j, 0, m)
	        Ainv[i][j] = 0.0;
    }
    else {
        for_less(i, 0, r) {
            for_less(j, 0, r) {
				if( i==j ) S[i][j] = 1/(W[i] + EPS);
				else S[i][j] = 0.0;
			}
		}
		    
        transpose(m, n, U, Ut);
    
        matrix_multiply(n, n, n, V, S, S);
        matrix_multiply(n, n, m, S, Ut, Ainv);
    }

    FREE(W);
    FREE2D(U);
    FREE2D(S);
    FREE2D(V);
    FREE2D(Ut);
    
    return( r );
}
