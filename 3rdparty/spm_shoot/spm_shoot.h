/*
 * Public interface of the CAT-Surface port of the SPM26 geodesic-shooting
 * primitives (3rdparty/spm_shoot): what spm_diffeo, spm_bsplinc and
 * spm_bsplins compute in MATLAB, without MATLAB.
 *
 * Include this header, not the ported shoot_*.h: those rename SPM's generic
 * names (bound, norm, grad, kernel, wrap, ...) by macros, this one declares
 * the spm_ prefixed functions directly and defines no such macro.
 *
 * Conventions (as in MATLAB; all arrays column-major, x fastest):
 *  - A field of n components on a dm[0] x dm[1] x dm[2] grid is n volumes
 *    one after the other: component c of voxel i is at [i + c*m], with
 *    m = dm[0]*dm[1]*dm[2].  Velocities, momenta and deformations have 3.
 *  - Deformations hold SPM's 1-based voxel coordinates: the identity maps
 *    voxel (i,j,k) (0-based) to (i+1, j+1, k+1).
 *  - Jacobian fields (9 components): element (r,c) = d y_r / d x_c is
 *    component r + 3*c, i.e. MATLAB's J(:,:,:,r+1,c+1).
 *  - Hessians (6 components, symmetric): H11 H22 H33 H12 H13 H23.
 *  - Regularisation parameters s[8] of vel2mom, kernel, relax, Atimesp,
 *    trapprox, fmg3 and cgs3: s[0..2] = 1/voxel size (the spm_diffeo
 *    gateways take voxel sizes and invert them: spm_shoot_param()),
 *    s[3] = absolute displacement, s[4] = membrane, s[5] = bending,
 *    s[6] = mu, s[7] = lambda (linear elasticity).  spm_diffeo('fmg')
 *    adds a small penalty to s[3] (spm_shoot_fmg_param()).
 *  - The boundary condition is global: spm_set_bound(0) (circulant, the
 *    default and what shooting uses) or 1 (Neumann).
 *  - Outputs that are accumulated (+=) must be zeroed by the caller, as
 *    mxCreateNumericArray does in MATLAB: kernel, push, pushc, pushc_grads.
 *
 * Threads: spm_shoot_set_threads(0) (the default) runs SPM's serial code.
 * With n >= 1 the loops SPM runs with OpenMP, the per-voxel loops of the
 * deformation code, pushc_grads, invdef, the multigrid transfer and the
 * B-splines run on up to n pthreads.  Every function then gives exactly
 * the serial result, for any n, except trapprox (and the internal sumsq),
 * whose sums are formed per slice and added in slice order: identical for
 * every n >= 1, but not bit-identical to n = 0.  push, pushc, cgs3, norm,
 * grad, divergence, minmax_div, bracket and kernel stay serial.  Calls must
 * not run concurrently from several threads (the boundary condition, the
 * error message and invdef's tables are global).
 *
 * Copyright (C) 2026 Christian Gaser, University of Jena.
 * Licence: GPL-2.0-or-later, as the SPM code it accompanies.
 */

#ifndef SPM_SHOOT_H
#define SPM_SHOOT_H

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

#ifndef SPM_SHOOT_MWTYPES
#define SPM_SHOOT_MWTYPES
typedef size_t    mwSize;
typedef size_t    mwIndex;
typedef ptrdiff_t mwSignedIndex;
#endif

/* ---------------------------------------------------------------------- */
/* Threads and errors                                                     */

/* 0 = SPM's serial code (bit-exact with the mex), n >= 1 = up to n threads */
void   spm_shoot_set_threads(int n);
int    spm_shoot_get_threads(void);
/* Minimum work per thread in ns (default 1e5); smaller jobs use fewer
   threads.  0 = always use all n threads (for tests).  Results do not
   depend on it. */
void   spm_shoot_set_grain(double ns);
double spm_shoot_get_grain(void);

/* SPM's mex functions abort with an error message; here the first message
   since the last spm_shoot_clear_error() is kept (NULL if none) and the
   call returns.  Raised by spm_set_bound (bad code), spm_invdef (a
   tetrahedron with more than 16383 voxels: the output is then incomplete),
   spm_bsplinc3, spm_shoot_fmg, spm_shoot_mat. */
const char *spm_shoot_get_error(void);
void        spm_shoot_clear_error(void);

/* ---------------------------------------------------------------------- */
/* spm_diffeo gateway equivalents                                         */

/* param = [1./prm(1:3) prm(4:8)], as spm_diffeo's vel2mom, kernel,
   trapprox and mom2vel gateways */
void spm_shoot_param(const double prm[8], double param[8]);

/* As spm_shoot_param, plus the penalty spm_diffeo('fmg') adds to param[3]:
   (wx000+wy000+wz000)*1.1921e-7*1e-4, the centres of the difference
   operators without the absolute term */
void spm_shoot_fmg_param(const double prm[8], double param[8]);

/* v = spm_diffeo('fmg', H, g, [prm cyc nit] [, v0]): solves (L+H) v = g
   with full multigrid (v0 NULL or v0[0] == 0) or with V-cycles from v0.
   dm: grid, H: 6-component Hessian (NULL: none), g: 3 components,
   prm: [vx vy vz lam0 lam1 lam2 mu lam] (voxel sizes in mm).  Allocates
   and zeroes the scratch space.  Returns 0, or -1 if out of memory. */
int  spm_shoot_fmg(const mwSize dm[3], const float *H, const float *g,
                   const double prm[8], int cyc, int nit, const float *v0,
                   float *v);

/* get_mat of spm_diffeo('invdef'): MATLAB 4x4 affine M (column-major, 16
   values) to the float M[4][3] of spm_invdef, Mf[col][row].  Returns -1
   (and sets the error) if the last row is not [0 0 0 1]. */
int  spm_shoot_mat(const double M[16], float Mf[4][3]);

/* ---------------------------------------------------------------------- */
/* Boundary condition (shoot_boundary.c)                                  */

void spm_set_bound(int t);         /* 0 circulant, 1 Neumann */
int  spm_get_bound(void);

/* ---------------------------------------------------------------------- */
/* Regularisation and multigrid (shoot_regularisers.c, shoot_optim3d.c,   */
/* shoot_multiscale.c); dm has 3 elements, s[8] as described above        */

/* g = L f (momentum from velocity); f, g: 3 components */
void   spm_vel2mom(mwSize dm[], float f[], double s[], float g[]);

/* Convolution kernel of L at the origin, wrapped around with the current
   boundary condition (spm_shoot_greens sets 0 before calling it).  f must be
   zeroed: m floats if s[6] == s[7] == 0 (scalar kernel, multiply by
   1/s[c]^2 for component c), else 9*m floats, block (r,c) at f + (r+3c)*m
   (MATLAB F(:,:,:,r+1,c+1), symmetric) */
void   spm_kernel(mwSize dm[], double s[], float f[]);

/* nit Gauss-Seidel sweeps of (L+H) u = b, in place; a: Hessian or NULL */
void   spm_relax(mwSize dm[], float a[], float b[], double s[], int nit,
                 float u[]);

/* Ap = (L + A) p; A: Hessian or NULL */
void   spm_Atimesp(mwSize dm[], float A[], double s[], float p[], float Ap[]);

/* approximate trace term used by spm_diffeo('trapprox'); a: Hessian */
double spm_trapprox(mwSize dm[], float a[], double s[]);

/* Full multigrid: solves (L+H) u = b.  n0: grid, a0: Hessian or NULL,
   b0: 3 components (not modified), param0: s[8], c: cycles, nit:
   relaxation iterations; u0: 3 components, the starting estimate for
   V-cycles, or full multigrid if u0[0] == 0 (only that element is tested;
   u0 is then overwritten, as the mex does with its zeros).  scratch:
   fmg3_scratchsize(n0, a0 != NULL) floats, zeroed before every call (the
   full-multigrid branch reads its coarsest level before writing it). */
void   spm_fmg3(mwSize n0[], float *a0, float *b0, double param0[], int c,
                int nit, float *u0, float *scratch);
mwSize spm_fmg3_scratchsize(mwSize n0[], int use_hessian);

/* Conjugate gradients (spm_diffeo('cgs')): x is the starting estimate;
   r, p, Ap are 3-component scratch fields */
void   spm_cgs3(mwSize dm[], float A[], float b[], double param[], double tol,
                int nit, float x[], float r[], float p[], float Ap[]);
double spm_norm(mwSize m, float a[]);

/* Multigrid transfer of one volume: restrict (linear, na -> nc =
   ceil(na/2)) and resize (2nd degree B-spline, any nc).  b: scratch of
   4*nc[0]*nc[1] + na[0]*nc[1] floats */
void   spm_restrict_vol(mwSize na[], float *a, mwSize nc[], float *c, float *b);
void   spm_resize_vol(mwSize na[], float *a, mwSize nc[], float *c, float *b);

/* ---------------------------------------------------------------------- */
/* Deformations (shoot_diffeo3d.c, shoot_expm3.c, shoot_invdef.c)         */

/* F1 = F0(Psi): n volumes F0 on grid dm0 sampled trilinearly at the m1
   points of the deformation Psi (3 x m1, 1-based).  pullc: circulant
   (any finite coordinate, NaN for non-finite); pull: NaN outside
   [1-0.05, dm+0.05] */
void spm_pullc(mwSize dm0[], mwSize m1, mwSize n, float Psi[], float F0[],
               float F1[]);
void spm_pull(mwSize dm0[], mwSize m1, mwSize n, float Psi[], float F0[],
              float F1[]);

/* Adjoint of pullc/pull: F0 (grid dm0, n volumes, zeroed) += F1 pushed
   through Psi (m1 points); S0 (dm0, zeroed, or NULL) counts the weights */
void spm_pushc(mwSize dm0[], mwSize m1, mwSize n, float Psi[], float F1[],
               float F0[], float S0[]);
void spm_push(mwSize dm0[], mwSize m1, mwSize n, float Psi[], float F1[],
              float F0[], float S0[]);

/* spm_diffeo('pushg', pf, Psi [, Jpsi] [, dmo]): pushes the 3-component
   field pf (grid dm) through Psi, multiplied by the inverse transposed
   Jacobian at each voxel (Jpsi, 9 components, or NULL: computed from Psi
   by central differences).  po: 3 components on grid dmo, zeroed */
void spm_pushc_grads(mwSize dmo[], mwSize dm[], float Psi[], float Jpsi[],
                     float pf[], float po[]);

/* C = B(A): B (3 components on grid dm) sampled trilinearly, circulant,
   at the mm points of A; C has mm points.  spm_diffeo('comp', B, A) =
   spm_composition(size(B), numel(A)/3, B, A, C).  Non-finite points of
   A give NaN.  C must not alias B. */
void spm_composition(mwSize dm[], mwSize mm, float *B, float *A, float *C);

/* As above with Jacobians: JC = JB(A) * JA (matrix product per voxel).
   spm_diffeo('comp', B, A, JB, JA) =
   spm_composition_jacobian(size(B), numel(A)/3, B, JB, A, JA, C, JC) */
void spm_composition_jacobian(mwSize dm[], mwSize mm, float *B, float *JB,
                              float *A, float *JA, float *C, float *JC);

/* As above with Jacobian determinants (1 component each): JC = JB(A)*JA */
void spm_composition_jacdet(mwSize dm[], mwSize mm, float *B, float *JB,
                            float *A, float *JA, float *C, float *JC);

/* Psi = Id + sc*V; _jac: Jpsi = I + sc*D V (central differences); _jac1:
   Jpsi = expm(sc*D V) (Pade), as spm_diffeo('smalldef', V, sc) with one
   and two outputs */
void spm_smalldef(mwSize dm[], float sc, float V[], float Psi[]);
void spm_smalldef_jac(mwSize dm[], float sc, float V[], float Psi[],
                      float Jpsi[]);
void spm_smalldef_jac1(mwSize dm[], float sc, float V[], float Psi[],
                       float Jpsi[]);

/* d = det(Jpsi) per voxel (computed in double) */
void spm_determinant(mwSize dm[], float Jpsi[], float d[]);

/* Jacobian determinants (def2det) or matrices (def2jac) of a deformation by
   central differences, circulant.  s = -1: whole volume (always use this;
   for a single slice s, SPM writes to that slice's offset of a one-slice
   output, out of bounds for s > 1) */
void spm_def2det(mwSize dm[], float *Psi, float *Jpsi, mwSignedIndex s);
void spm_def2jac(mwSize dm[], float *Psi, float *Jpsi, mwSignedIndex s);

/* iy = inverse of y (spm_diffeo('invdef', y, dim_iy, M1, M2)): y on grid
   dim_y, iy (fully written, NaN where undefined) on grid dim_iy; M1 maps
   the values of y to voxels of iy, M2 the result to its values.  M1, M2
   from spm_shoot_mat() */
void spm_invdef(mwSize dim_y[3], float y[], mwSize dim_iy[3], float iy[],
                float M1[4][3], float M2[4][3]);

/* Gradients, divergence, Lie bracket and 3x3 matrix exponentials */
void spm_grad(mwSize dm[], mwSize N, float F[], float D[]);
void spm_divergence(mwSize dm[], float V[], float dv[]);
void spm_minmax_div(mwSize dm[], float V[], double mnmx[]);
void spm_bracket(mwSize dm[], float *A, float *B, float *C);
void spm_expm22(float a[], float l[]);
void spm_expm33(float a[], float l[]);

/* ---------------------------------------------------------------------- */
/* B-splines (bsplines.c, double precision)                               */

/* c = spm_bsplinc(vol, [deg wrapflags]): copy of vol (non-finite -> 0),
   deconvolved along x, y, z.  deg: 0..7 per dimension, wrapflags: 1 =
   periodic, 0 = mirrored.  coef: dim[0]*dim[1]*dim[2] doubles.
   Returns 0, -1 for a bad degree, 1 if dim[1] or dim[2] exceed 10240. */
int  spm_bsplinc3(const float *vol, const mwSize dim[3], const int deg[3],
                  const int wrapflags[3], double *coef);

/* [f,d0,d1,d2] = spm_bsplins(coef, x0, x1, x2, [deg wrapflags]) at n
   points with 1-based coordinates.  In mirrored dimensions points outside
   [1-0.05, dim+0.05] give f = NaN (derivatives 0); in periodic ones every
   point is sampled (non-finite or |x| > 2^30: f = 0, derivatives 0, which
   is what the arm64 mex returns).  d0, d1, d2 may be NULL (all NULL: no
   derivatives). */
void spm_bsplins3(const double *coef, const mwSize dim[3], mwSize n,
                  const double *x0, const double *x1, const double *x2,
                  const int deg[3], const int wrapflags[3],
                  double *f, double *d0, double *d1, double *d2);

#ifdef __cplusplus
}
#endif

#endif /* SPM_SHOOT_H */
