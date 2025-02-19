/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 *
 *  based on depth_potential.c by Maxime Boucher
 *
 *   References:
 *
 *  M. Boucher, S. Whitesides, and A. Evans, "Depth potential function for 
 *      folding pattern representation, registration and analysis", Medical
 *      Image Analysis, Volume 13, Issue 2, pp 203-214, April 2009.
 */


#include <bicpl.h>
#include <float.h>

#include "CAT_Surf.h"
#include "CAT_DepthPotential.h"

#define vec_sub( a, b, c ) { \
    c.coords[0] = a.coords[0] - b.coords[0]; \
    c.coords[1] = a.coords[1] - b.coords[1]; \
    c.coords[2] = a.coords[2] - b.coords[2]; \
}

#define vec_dot_product( a, b ) \
    ( a.coords[0] * b.coords[0] + a.coords[1] * b.coords[1] + a.coords[2] * b.coords[2] )

#define vec_normalize( a ) { \
    double aa = sqrt( a.coords[0] * a.coords[0] + a.coords[1] * a.coords[1] + \
                                    a.coords[2] * a.coords[2] ); \
    a.coords[0] /= aa; \
    a.coords[1] /= aa; \
    a.coords[2] /= aa; \
}

#define vec_cross( a, b, c ) { \
    c.coords[0] = a.coords[1] * b.coords[2] - a.coords[2] * b.coords[1]; \
    c.coords[1] = a.coords[2] * b.coords[0] - a.coords[0] * b.coords[2]; \
    c.coords[2] = a.coords[0] * b.coords[1] - a.coords[1] * b.coords[0]; \
}

double *
compute_depth_potential( polygons_struct *polygons, double alpha) 
{

    Point    * coords = NULL;            // coordinates
    Vector * normals = NULL;         // normal vectors
    double * areas = NULL;           // Voronoi area at nodes
    double * mc = NULL;                  // mean curvature
    double * dp = NULL;                  // depth potential parameter
    int      * n_neighbours = NULL;             // node neighbours (inverse connectivity)
    int     ** neighbours = NULL;
    double SOR = 1.90;           // this value works best with alpha=0.0015
    int lambda = 2, i;

    struct csr_matrix laplacian; // the laplacian matrix operator

    // Defaults

    normals = (Vector *) malloc(sizeof(Vector) * polygons->n_points);

    create_polygon_point_neighbours(polygons, TRUE, &n_neighbours, &neighbours, NULL, NULL);

    stable_normals( polygons->n_points, polygons->points, normals, n_neighbours, neighbours );
    areas = compute_areas( polygons->n_points, polygons->points, n_neighbours, neighbours, lambda );
    init_csr_matrix( polygons->n_points, n_neighbours, neighbours, &laplacian );
    cot_laplacian_operator( polygons->n_points, polygons->points, &laplacian, n_neighbours, neighbours );
    mc = compute_mean_curvature( polygons->n_points, polygons->points, areas, normals, &laplacian );
    dp = local_depth_potential( polygons->n_points, polygons->points, areas, &laplacian, mc, 
                                                                    alpha, SOR );

    if( normals ) FREE( normals );
    if( n_neighbours ) FREE( n_neighbours );
    if( neighbours ) {
        FREE( neighbours[0] );   // this is ngh_array
        FREE( neighbours );
    }
    if( areas ) free( areas );
    if( mc ) free( mc );
    
    return(dp);
}

/* -------------------------------------------------------------------
 Compute Voronoi area around nodes.

    if lambda ==1
            area_mat = sparse([1:nb_tri,1:nb_tri,1:nb_tri],[tri(1,:),tri(2,:),tri(3,:)],
                                                (1/6)*[area_tri,area_tri,area_tri]);
    else
            %mixed area
            good = (dot_product(v1,v2)<=0).*(dot_product(v2,v3)<=0).*(dot_product(v3,v1)<=0);
            lgood = find(good);
            lbad = find(1-good);
            area_mat = sparse([lgood,lgood,lgood],[tri(1,lgood),tri(2,lgood),tri(3,lgood)],
                                -0.125*[(norms(3,lgood).*dot_product(v1(:,lgood),v2(:,lgood))+
                                                 norms(1,lgood).*dot_product(v2(:,lgood),v3(:,lgood)))./
                                                area_tri(lgood),
                                                (norms(2,lgood).*dot_product(v1(:,lgood),v3(:,lgood))+
                                                 norms(1,lgood).*dot_product(v2(:,lgood),v3(:,lgood)))./
                                                area_tri(lgood),
                                                (norms(2,lgood).*dot_product(v1(:,lgood),v3(:,lgood))+
                                                 norms(3,lgood).*dot_product(v1(:,lgood),v2(:,lgood)))
                                                 ./area_tri(lgood)],nb_tri,nb_points);
            area_mixed = sparse([lbad,lbad,lbad],[tri(1,lbad),tri(2,lbad),tri(3,lbad)],
                                     0.125*[area_tri(lbad),area_tri(lbad),area_tri(lbad)],nb_tri,nb_points);
            area_mat = area_mat + area_mixed;
            bad1 = find(dot_product(v3,v1)>0);
            bad2 = find(dot_product(v1,v2)>0);
            bad3 = find(dot_product(v2,v3)>0);
            area_mixed = sparse([bad1,bad2,bad3],[tri(1,bad1),tri(2,bad2),tri(3,bad3)],
                                     0.125*[area_tri(bad1),area_tri(bad2),area_tri(bad3)],nb_tri,nb_points);
            area_mat = area_mat + area_mixed;
    end
    areas = full(sum(area_mat));
*/ 
double * 
compute_areas( int n_points, Point coords[], int * n_ngh, 
               int ** ngh, int lambda ) {

    int         i, n1, n2, nn;
    double      area_tri;
    Vector  v1, v2, v3, cross;

    double * areas = (double*)malloc( n_points * sizeof( double ) );
    for( i = 0; i < n_points; i++ ) {
        areas[i] = 0.0;
    }

    /* Loop over the triangles to compute the area around nodes. */

    for( i = 0; i < n_points; i++ ) {
        for( nn = 0; nn < n_ngh[i]; nn++ ) {
            n1 = ngh[i][nn];
            n2 = ngh[i][(nn+1)%n_ngh[i]];
            if( i < n1 && i < n2 ) {
                vec_sub( coords[n1], coords[i], v1 );
                vec_sub( coords[n2], coords[n1], v2 );
                vec_sub( coords[i], coords[n2], v3 );
                vec_cross( v1, v2, cross );
                area_tri = sqrt( vec_dot_product( cross, cross ) );
                if( lambda == 1 ) {
                    /* equal weights */
                    areas[i] += area_tri / 6.0;
                    areas[n1] += area_tri / 6.0;
                    areas[n2] += area_tri / 6.0;
                } else {
                    /* Voronoi */
                    int bad1 = vec_dot_product( v3, v1 ) > 0;
                    int bad2 = vec_dot_product( v1, v2 ) > 0;
                    int bad3 = vec_dot_product( v2, v3 ) > 0;
                    int bad = bad1 || bad2 || bad3;
                    if( !bad ) {
                        double w1 = vec_dot_product( v3, v3 ) * vec_dot_product( v1, v2 ) +
                                            vec_dot_product( v1, v1 ) * vec_dot_product( v2, v3 );
                        double w2 = vec_dot_product( v2, v2 ) * vec_dot_product( v1, v3 ) +
                                            vec_dot_product( v1, v1 ) * vec_dot_product( v2, v3 );
                        double w3 = vec_dot_product( v2, v2 ) * vec_dot_product( v1, v3 ) +
                                            vec_dot_product( v3, v3 ) * vec_dot_product( v1, v2 );
                        areas[i] -= 0.125 * w1 / area_tri;
                        areas[n1] -= 0.125 * w2 / area_tri;
                        areas[n2] -= 0.125 * w3 / area_tri;
                    } else {
                        areas[i] += 0.125 * area_tri;
                        areas[n1] += 0.125 * area_tri;
                        areas[n2] += 0.125 * area_tri;
                        if( bad1 ) areas[i] += 0.125 * area_tri;
                        if( bad2 ) areas[n1] += 0.125 * area_tri;
                        if( bad3 ) areas[n2] += 0.125 * area_tri;
                    }
                }
            }
        }
        if ((isnan(areas[i])) || (areas[i] == 0.0)) areas[i] = 1.0;
    }
    return( areas );
}

/* -------------------------------------------------------------------
 Compute mean curvature.

 mc obtained using: mean_curvature white.obj mc.dat
 which is: [mat,area] = cot_laplacian_operator(coords,tri);
                n = stable_normals(coords,tri);
                output_vv(output_file,0.5*dot_product(coords*mat,n)./area);
*/
private double * 
compute_mean_curvature( int n_points, Point coords[], double * areas,
                         Vector normals[], struct csr_matrix * mat ) {

    int         i, j;
    Vector  vv;

    double * mc = (double*)malloc( n_points * sizeof( double ) );
    for( i = 0; i < n_points; i++ ) {
        mc[i] = 0.0;
    }

    /* Loop over the triangles to compute the mean curvature. */

    for( i = 0; i < n_points; i++ ) {
        vv.coords[0] = 0.0;
        vv.coords[1] = 0.0;
        vv.coords[2] = 0.0;
        for( j = mat->ia[i]; j < mat->ia[i+1]; j++ ) {
            vv.coords[0] += mat->A[j] * coords[mat->ja[j]].coords[0];
            vv.coords[1] += mat->A[j] * coords[mat->ja[j]].coords[1];
            vv.coords[2] += mat->A[j] * coords[mat->ja[j]].coords[2];
        }
        mc[i] = 0.5 * vec_dot_product( vv, normals[i] ) / areas[i];
        if (isnan(mc[i])) mc[i] = 0.0;
    }
    return( mc );
}

/* -------------------------------------------------------------------
 Creation of the matrix for the Laplacian operator.

 function [mat,areas,mc,area_mat,gc] = cot_laplacian_operator(coords,tri,lambda)
 
 nb_points = size(coords,2);
 nb_tri = size(tri,2);
 v1 = coords(:,tri(2,:)) - coords(:,tri(1,:));
 v2 = coords(:,tri(3,:)) - coords(:,tri(2,:));
 v3 = coords(:,tri(1,:)) - coords(:,tri(3,:));
 
 norms(1,:) = norm_vec(v1).^2;
 norms(2,:) = norm_vec(v2).^2;
 norms(3,:) = norm_vec(v3).^2;
 
 area_tri = norm_vec(cross(v1,v2));
 mat = sparse([1:nb_points,tri(2,:),tri(3,:),tri(1,:)],
                            [1:nb_points,tri(1,:),tri(2,:),tri(3,:)],
                            [ones(1,nb_points),dot_product(v2,v3)./area_tri,
                             dot_product(v1,v3)./area_tri,dot_product(v1,v2)./area_tri]);
 mat = 0.5*(mat + mat');
 mat = mat - sparse(1:nb_points,1:nb_points,sum(mat));
*/
private void 
cot_laplacian_operator( int n_points, Point coords[], struct csr_matrix * mat,
                        int * n_ngh, int ** ngh ) {

    int         i, n1, n2, nn;
    double      area_tri;
    Vector  v1, v2, v3, cross;

    /* Loop over the triangles to compute the coefficients of the matrix. */

    for( i = 0; i < n_points; i++ ) {
        for( nn = 0; nn < n_ngh[i]; nn++ ) {
            n1 = ngh[i][nn];
            n2 = ngh[i][(nn+1)%n_ngh[i]];
            if( i < n1 && i < n2 ) {
                vec_sub( coords[n1], coords[i], v1 );
                vec_sub( coords[n2], coords[n1], v2 );
                vec_sub( coords[i], coords[n2], v3 );

                vec_cross( v1, v2, cross );
                area_tri = sqrt( vec_dot_product( cross, cross ) );
                
                if (area_tri == 0) area_tri = 1.0;
                double dot23 = vec_dot_product( v2, v3 ) / area_tri;
                double dot13 = vec_dot_product( v1, v3 ) / area_tri;
                double dot12 = vec_dot_product( v1, v2 ) / area_tri;

                /* Assemble dot23 at (i,n1) and (n1,i); dot13 at (n1,n2) and (n2,n1);
                     dot12 at (i,n2) and (n2,i). assemble() is symmetric. */
                assemble( dot23, i, n1, mat );
                assemble( dot13, n1, n2, mat );
                assemble( dot12, i, n2, mat );

                assemble( -0.5*(dot23+dot12), i, i, mat );
                assemble( -0.5*(dot13+dot23), n1, n1, mat );
                assemble( -0.5*(dot12+dot13), n2, n2, mat );
            }
        }       
    }
}

/* -------------------------------------------------------------------
 Compute depth potential.

 Linear system is:
 [mat,area] = cot_laplacian_operator(coords,triangles);
 mat = 0.5*mat;
 mc = mean_curvature(coords,triangles);
 c = (mc*area-area*sum(mc*area)/sum(area));
 mat = mat + alpha*sparse(1:nb_coords,1:nb_coords,area);
 
 [R]=cholinc(mat,1e-3);
 res = pcg(mat,c',1e-12,10000,R',R)';
*/
double * 
local_depth_potential( int n_points, Point coords[], double * areas,
                      struct csr_matrix * mat, double * mc, 
                      double alpha, double SOR ) {

    int i, j;

    double sum_area = 0.0;
    double sum_area_mc = 0.0;
    for( i = 0; i < n_points; i++ ) {
        sum_area += areas[i];
        sum_area_mc += mc[i] * areas[i];
    }

    double * rhs = (double*)malloc( n_points * sizeof( double ) );
    double factor = sum_area_mc / sum_area;
    for( i = 0; i < n_points; i++ ) {
        rhs[i] = ( mc[i] - factor ) * areas[i];
    }

    double * dp_mat_array = (double*)malloc( mat->nnz * sizeof( double ) );

    for( i = 0; i < n_points; i++ ) {
        for( j = mat->ia[i]; j < mat->ia[i+1]; j++ ) {
            dp_mat_array[j] = 0.50 * mat->A[j];
            if( mat->ja[j] == i ) {
                dp_mat_array[j] += alpha * areas[i];
            }
        }
    }

    /* Solution by SOR Gauss Seidel */

    double * dp = (double*)malloc( n_points * sizeof( double ) );

    gauss_seidel( n_points, mat->ia, mat->ja, dp_mat_array, rhs, dp, 1000, 1.0e-10, SOR, 0 );

    free( dp_mat_array );
    free( rhs );

    return( dp );
}

/* -------------------------------------------------------------------
 Solve a linear system with SOR Gauss-Seidel. The matrix is in CSR
 format. The convergence threshold is a relative tolerance from the
 initial residual, up to a maximum number of iterations, whichever 
 comes first.
*/
private void 
gauss_seidel( int nnode, int * ia, int * ja, double * mat, 
             double * rhs, double * soln, int max_iter, 
             double tol, double SOR, int verbose ) {

    int i, j, iter;
    double res, res0, val, diag;

    for( i = 0; i < nnode; i++ ) {
        soln[i] = 0;
    }
 
    for( iter = 1; iter <= max_iter; iter++ ) { 
        res = 0.0;
        for( i = 0; i < nnode; i++ ) {
            val = rhs[i];
            for( j = ia[i]; j < ia[i+1]; j++ ) {
                if( i != ja[j] ) {
                    val -= mat[j] * soln[ja[j]];
                } else {
                    diag = mat[j];
                }
            }
            val /= diag;
            val = soln[i] + SOR * ( val - soln[i] );
            res += fabs( val - soln[i] );
            soln[i] = val;
        }
        res /= (double)nnode;
        if( iter == 1 ) res0 = res;
        res /= res0;
        if( verbose && iter%10 == 1 ) printf( "Iteration %d  Res %g\n", iter, res );
        if( res < tol ) break;
    }
    if( verbose ) printf( "Last iteration %d    Res %g\n", iter, res );
}

/* -------------------------------------------------------------------
 Initialize the data structures for the sparse matrix in CSR format
 based on the connectivity (assumed to be triangles). The matrix is
 initialized to zero.
*/
private void 
init_csr_matrix( int n_points, int * n_ngh, int ** ngh, 
                 struct csr_matrix * mat ) {

    int i, j, nnz;

    mat->nnz = n_points;
    for( i = 0; i < n_points; i++ ) {
        mat->nnz += n_ngh[i];
    }

    mat->n = n_points;
    mat->ia = (int *)malloc( ( n_points + 1 ) * sizeof( int ) );
    mat->ja = (int *)malloc( mat->nnz * sizeof( int ) );
    mat->A = (double *)malloc( mat->nnz * sizeof( double ) );

    for( i = 0; i < mat->nnz; i++ ) {
        mat->A[i] = 0.0;
    }

    nnz = 0;
    mat->ia[0] = 0;
    for( i = 0; i < n_points; i++ ) {
        mat->ia[i+1] = mat->ia[i] + n_ngh[i] + 1;
        mat->ja[nnz] = i;
        nnz++;
        for( j = 0; j < n_ngh[i]; j++ ) {
            mat->ja[nnz] = ngh[i][j];
            nnz++;
        }
    }
}

/* -------------------------------------------------------------------
     Free the memory for the data structures of a sparse matrix in CSR 
     format.
*/
private void 
free_csr_matrix( struct csr_matrix * mat ) {

    free( mat->ia );
    free( mat->ja );
    free( mat->A );

    mat->nnz = 0;
    mat->n = 0;
    mat->ia = NULL;
    mat->ja = NULL;
    mat->A = NULL;
}

/* -------------------------------------------------------------------
     Add a coefficient to a symmetric sparse matrix in CSR format, at
     (row,col) and (col,row).
*/
private void 
assemble( double val, int row, int col, struct csr_matrix * mat ) {

    int j;

    for( j = mat->ia[row]; j < mat->ia[row+1]; j++ ) {
        if( col == mat->ja[j] ) mat->A[j] += val;
    }
    for( j = mat->ia[col]; j < mat->ia[col+1]; j++ ) {
        if( row == mat->ja[j] ) mat->A[j] += val;
    }
}

/* -------------------------------------------------------------------
 Compute the normal vectors.

 function res = stable_normals(coords,tri)
 
 nb_points = size(coords,2);
 nb_tri = size(tri,2);
 v1 = normalize(coords(:,tri(2,:)) - coords(:,tri(1,:)));
 v2 = normalize(coords(:,tri(3,:)) - coords(:,tri(2,:)));
 v3 = normalize(coords(:,tri(1,:)) - coords(:,tri(3,:)));
 
 n = normalize(cross(v1,v2));
 
 w(1,:) = (1-dot_product(v3,v1))./norm_vec(cross(v3,v1));
 w(2,:) = (1-dot_product(v1,v2))./norm_vec(cross(v1,v2));
 w(3,:) = (1-dot_product(v2,v3))./norm_vec(cross(v2,v3));
 
 res(1,:) = sum(sparse([1:nb_tri,1:nb_tri,1:nb_tri],
                                             [tri(1,:),tri(2,:),tri(3,:)],
                                             [n(1,:).*w(1,:),n(1,:).*w(2,:),n(1,:).*w(3,:)]));
 res(2,:) = sum(sparse([1:nb_tri,1:nb_tri,1:nb_tri],
                                             [tri(1,:),tri(2,:),tri(3,:)],
                                             [n(2,:).*w(1,:),n(2,:).*w(2,:),n(2,:).*w(3,:)]));
 res(3,:) = sum(sparse([1:nb_tri,1:nb_tri,1:nb_tri],
                                             [tri(1,:),tri(2,:),tri(3,:)],
                                             [n(3,:).*w(1,:),n(3,:).*w(2,:),n(3,:).*w(3,:)]));
 res = normalize(full(res));
*/
private void 
stable_normals( int n_points, Point coords[], Vector normals[], 
                int * n_ngh, int ** ngh) {

    int         i, j, n1, n2, nn;
    double      mag1, mag2, mag3, mag3x1, mag1x2, mag2x3, w1, w2, w3;
    Vector  v1, v2, v3, cross;

    for( i = 0; i < n_points; i++ ) {
        normals[i].coords[0] = 0.0;
        normals[i].coords[1] = 0.0;
        normals[i].coords[2] = 0.0;
    }

    for( i = 0; i < n_points; i++ ) {
        for( nn = 0; nn < n_ngh[i]; nn++ ) {

            n1 = ngh[i][nn];
            n2 = ngh[i][(nn+1)%n_ngh[i]];
            if( i < n1 && i < n2 ) {

                vec_sub( coords[n1], coords[i], v1 );
                vec_sub( coords[n2], coords[n1], v2 );
                vec_sub( coords[i], coords[n2], v3 );
                vec_normalize( v1 );
                vec_normalize( v2 );
                vec_normalize( v3 );

                /* weights */
                vec_cross( v3, v1, cross );
                mag3x1 = sqrt( vec_dot_product( cross, cross ) );
                vec_cross( v2, v3, cross );
                mag2x3 = sqrt( vec_dot_product( cross, cross ) );
                vec_cross( v1, v2, cross );
                mag1x2 = sqrt( vec_dot_product( cross, cross ) );

                w1 = ( 1.0 - vec_dot_product( v3, v1 ) ) / mag3x1;
                w2 = ( 1.0 - vec_dot_product( v1, v2 ) ) / mag1x2;
                w3 = ( 1.0 - vec_dot_product( v2, v3 ) ) / mag2x3;

                /* contribution of tri_norm to each node of triangle */
                for( j = 0; j < 3; j++ ) {
                    normals[i].coords[j] += w1 * cross.coords[j] / mag1x2;
                    normals[n1].coords[j] += w2 * cross.coords[j] / mag1x2;
                    normals[n2].coords[j] += w3 * cross.coords[j] / mag1x2;
                }
            }
        }
    }

    /* normalization of normals */
    for( i = 0; i < n_points; i++ ) {
        vec_normalize( normals[i] );
    }
}
