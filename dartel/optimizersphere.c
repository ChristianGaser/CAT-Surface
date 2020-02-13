/* $Id$ 
   largely modified version of optimizer2d.c from John Ashburner 
   (c) John Ashburner (2007) */

#include<math.h>
#include "dartel.h"
#include <stdio.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#define S 1.0000000001


static int
bound(int i, int j, int dm[])
{
        int i1, j1;
        int m2;

        /* circulant boundary condition for x-coordinates */
        i1 = (i >= 0 ? (i % dm[0]) : (dm[0] + (i % dm[0])) % dm[0]);

        /* Neumann boundary condition for y-coordinates */ 
        if (dm[1] == 1) {
                j1 = 0;
        } else {
                m2 = dm[1] * 2;
                j = (j < 0) ? (-j - m2*(-j / m2) - 1) : (j - m2*(j / m2));
                if (dm[1] <= j)
                        j1 = m2 - j - 1;
                else
                        j1 = j;
        } 
        return(i1 + dm[0]*j1);
}

/* Neumann boundary conditions */
static int
neumann(int i, int m)
{
	if (m==1)
		return(0);
	else {
		int m2 = m*2;
		i = (i < 0) ? (-i - m2*((-i)/m2) - 1) : (i - m2*(i/m2));
		if (m <= i)
			return(m2 - i - 1);
		else
			return(i);
	}
}
#    define BOUNDY(i,m) neumann(i,m)
#    define BOUNDX(i,m) (((i)>=0) ? (i)%(m) : ((m)+(i)%(m))%m)

double
sumsq_le(int dm[], double a[], double b[], double s[], double u[])
{
    double ss = 0.0;
    int i, j;
    double mu = s[2], lam = s[3], id = s[4];
    double wx0, wx1, wx2, wy0, wy1, wy2, wxy;
    double *pux, *puy, *pbx, *pby, *paxx, *paxy, *payy;

    wx0 = mu*(2*s[1]*s[1] + 4*s[0]*s[0]) + 2*lam*s[0]*s[0] + id;
    wy0 = mu*(2*s[0]*s[0] + 4*s[1]*s[1]) + 2*lam*s[1]*s[1] + id;
    wx1 = -(2*mu + lam)*s[0]*s[0];
    wy1 = -(2*mu + lam)*s[1]*s[1];
    wx2 = -mu*s[0]*s[0];
    wy2 = -mu*s[1]*s[1];
    wxy = (0.25*lam + 0.25*mu)*s[0]*s[1];

    pux  = u;
    puy  = u + dm[0]*dm[1];
    pbx  = b;
    pby  = b + dm[0]*dm[1];
    paxx = a;
    payy = a + dm[0]*dm[1];
    paxy = a + dm[0]*dm[1]*2;

    for (j = 0; j < dm[1]; j++) {
        for (i = 0; i < dm[0]; i++) {
            int ij;
            double tmp;

            ij  = bound(i  ,j  ,dm);

            tmp = (wx0 + paxx[ij])*pux[ij] + paxy[ij]*puy[ij]
                 + wy2*(pux[bound(i,j-1,dm)] + pux[bound(i,j+1,dm)])
                 + wx1*(pux[bound(i-1,j,dm)] + pux[bound(i+1,j,dm)])
                 + wxy*(puy[bound(i-1,j-1,dm)] - puy[bound(i-1,j+1,dm)] -
                        puy[bound(i+1,j-1,dm)] + puy[bound(i+1,j+1,dm)])
                 - pbx[ij];
            ss += tmp*tmp;

            tmp = (wy0 + payy[ij])*puy[ij] + paxy[ij]*pux[ij]
                 + wy1*(puy[bound(i,j-1,dm)] + puy[bound(i,j+1,dm)])
                 + wx2*(puy[bound(i-1,j,dm)] + puy[bound(i+1,j,dm)])
                 + wxy*(pux[bound(i-1,j-1,dm)] - pux[bound(i-1,j+1,dm)] -
                        pux[bound(i+1,j-1,dm)] + pux[bound(i+1,j+1,dm)])
                 - pby[ij];
            ss += tmp*tmp;
        }
    }
    return(ss);
}

void
LtLf_le(int dm[], double f[], double s[], double g[])
{
    int i, j, jm1,jp1, im1,ip1;
    double *pgx, *pgy, *pfx, *pfy;
    double mu = s[2], lam = s[3], id = s[4];
    double wx0, wx1, wx2, wy0, wy1, wy2, wxy;

    wx0 = mu*(2*s[1]*s[1] + 4*s[0]*s[0]) + 2*lam*s[0]*s[0] + id;
    wy0 = mu*(2*s[0]*s[0] + 4*s[1]*s[1]) + 2*lam*s[1]*s[1] + id;
    wx1 = -(2*mu + lam)*s[0]*s[0];
    wy1 = -(2*mu + lam)*s[1]*s[1];
    wx2 = -mu*s[0]*s[0];
    wy2 = -mu*s[1]*s[1];
    wxy = (0.25*lam + 0.25*mu)*s[0]*s[1];

    pgx = g;
    pgy = g + dm[0]*dm[1];
    pfx = f;
    pfy = f + dm[0]*dm[1];

#   pragma omp parallel for collapse(2) private(i)
    for (j = 0; j < dm[1]; j++) {
        for (i = 0; i < dm[0]; i++) {
            int ij;

            ij = bound(i,j,dm);

            pgx[ij] = wx0*pfx[ij] +
                      wy2*(pfx[bound(i,j-1,dm)] + pfx[bound(i,j+1,dm)]) +
                      wx1*(pfx[bound(i-1,j,dm)] + pfx[bound(i+1,j,dm)]) +
                      wxy*(pfy[bound(i-1,j-1,dm)] - pfy[bound(i-1,j+1,dm)] -
                           pfy[bound(i+1,j-1,dm)] + pfy[bound(i+1,j+1,dm)]);
            pgy[ij] = wy0*pfy[ij] +
                      wy1*(pfy[bound(i,j-1,dm)] + pfy[bound(i,j+1,dm)]) +
                      wx2*(pfy[bound(i-1,j,dm)] + pfy[bound(i+1,j,dm)]) +
                      wxy*(pfx[bound(i-1,j-1,dm)] - pfx[bound(i-1,j+1,dm)] -
                           pfx[bound(i+1,j-1,dm)] + pfx[bound(i+1,j+1,dm)]);
        }
    }
}

void
relax_le(int dm[], double a[], double b[], double s[], int nit, double u[])
{
    int it;
    double mu = s[2], lam = s[3], id = s[4];
    double wx0, wx1, wx2, wy0, wy1, wy2, wxy;
    double regx, regy;

    wx0 = mu*(2*s[1]*s[1] + 4*s[0]*s[0]) + 2*lam*s[0]*s[0] + id;
    wy0 = mu*(2*s[0]*s[0] + 4*s[1]*s[1]) + 2*lam*s[1]*s[1] + id;
    wx1 = -(2*mu + lam)*s[0]*s[0];
    wy1 = -(2*mu + lam)*s[1]*s[1];
    wx2 = -mu*s[0]*s[0];
    wy2 = -mu*s[1]*s[1];
    wxy = (0.25*lam + 0.25*mu)*s[0]*s[1];

    regx = (4.0*wxy - 2.0*(wx1+wy2)) - wx0; if (regx<0) regx = 0.0;
    regy = (4.0*wxy - 2.0*(wx2+wy1)) - wy0; if (regy<0) regy = 0.0;

#ifdef VERBOSE
    for (it = 0; it < 10 - (int)ceil(1.44269504088896*log((double)dm[0])); it++)
            printf("  ");
    printf("%dx%d: ", dm[0],dm[1]);
#endif

    for (it = 0; it < 2*nit; it++) {
        int j, jstart;
        int i, istart;

#ifdef VERBOSE
        printf(" %g", sumsq_le(dm, a, b, s, u));
#endif

        jstart = it%2;
        for (j = 0; j < dm[1]; j++) {
            double *pux, *puy, *pbx, *pby, *paxx, *paxy, *payy;
            int jm1, jp1, im1, ip1;

            pux  = u + dm[0]*j;
            puy  = u + dm[0]*(j + dm[1]);
            pbx  = b + dm[0]*j;
            pby  = b + dm[0]*(j + dm[1]);
            paxx = a + dm[0]*j;
            payy = a + dm[0]*(j + dm[1]);
            paxy = a + dm[0]*(j + dm[1]*2);

            jm1 = (BOUNDY(j-1,dm[1]) - j)*dm[0];
    	    jp1 = (BOUNDY(j+1,dm[1]) - j)*dm[0];

            istart = (jstart == (j%2));

            for (i = istart; i < dm[0]; i+=2) {
                double sux, suy, axx, ayy, axy, idt;
                double *px = pux+i, *py = puy+i;

                im1 = BOUNDX(i-1,dm[0]) - i;
    	        ip1 = BOUNDX(i+1,dm[0]) - i;
	    
                sux = pbx[i] - ((wx0 + paxx[i])*px[0] + paxy[i]*py[0] +
                                wy2*(px[jm1] + px[jp1]) +
                                wx1*(px[im1] + px[ip1]) +
                                wxy*(py[jp1+im1] - py[jp1+ip1] -
                                     py[jm1+im1] + py[jm1+ip1]));

                suy = pby[i] - (paxy[i]*px[0] + (wy0 + payy[i])*py[0] +
                                wy1*(py[jm1] + py[jp1]) +
                                wx2*(py[im1] + py[ip1]) +
                                wxy*(px[jp1+im1] - px[jp1+ip1] -
                                     px[jm1+im1] + px[jm1+ip1]));

                /*
                   syms axx ayy axy sux suy
                   A = [axx axy; axy ayy];
                   su = [sux ; suy]
                   inv(A)*su
                */
                axx  = paxx[i] + wx0 + regx;
                ayy  = payy[i] + wy0 + regy;
                axy  = paxy[i];
                idt  = 1.0/(axx*ayy - axy*axy);

                *px += idt*( ayy*sux - axy*suy);
                *py += idt*(-axy*sux + axx*suy);
            }
        }
    }
#ifdef VERBOSE
    printf(" %g\n", sumsq_le(dm, a, b, s, u));
#endif
}


void
Atimesp_le(int dm[], double A[], double s[], double p[], double Ap[])
{
    int i, m = dm[0]*dm[1];

    LtLf_le(dm, p, s, Ap);
#   pragma omp parallel for
    for (i = 0; i < m; i++) {
        Ap[i  ] += A[i  ]*p[i  ] + A[i+2*m]*p[i+m];
        Ap[i+m] += A[i+m]*p[i+m] + A[i+2*m]*p[i  ];
    }
}

double
sumsq_me(int dm[], double a[], double b[], double s[], double u[])
{
    double w00,w01, w10;
    double ss = 0.0;
    int i, j, m = dm[0]*dm[1];

    w00 = s[2]*(2*s[0]*s[0] + 2*s[1]*s[1]) + s[4];
    w01 = s[2]*(-s[1]*s[1]);
    w10 = s[2]*(-s[0]*s[0]);

    for (j = 0; j < dm[1]; j++) {
        double *pux, *puy, *pbx, *pby, *paxx, *paxy, *payy;
        int jm1,jp1,im1,ip1;

        pux  = u;
        puy  = u + dm[0]*dm[1];
        pbx  = b;
        pby  = b + dm[0]*dm[1];
        paxx = a;
        payy = a + dm[0]*dm[1];
        paxy = a + dm[0]*dm[1]*2;

        for (i = 0; i < dm[0]; i++) {
            double tmp;
            int ij;

            ij = bound(i,j,dm);

            tmp = (w00 + paxx[ij])*pux[ij] + paxy[ij]*puy[ij] +
                  w01*(pux[bound(i,j-1,dm)] + pux[bound(i,j+1,dm)]) +
                  w10*(pux[bound(i-1,j,dm)] + pux[bound(i+1,j,dm)]) -
                  pbx[ij];
            ss += tmp*tmp;

            tmp = (w00 + payy[ij])*puy[ij] + paxy[ij]*pux[ij] +
                  w01*(puy[bound(i,j-1,dm)] + puy[bound(i,j+1,dm)]) +
                  w10*(puy[bound(i-1,j,dm)] + puy[bound(i+1,j,dm)]) -
                  pby[ij];
            ss += tmp*tmp;
        }
    }
    return(ss);
}

#define _PI 3.14159265358979323846264338327510

void
LtLf_me(int dm[], double f[], double s[], double g[])
{
    int i, j, m = dm[0]*dm[1];
    double w00,w01,w10;

    w00 = s[2]*(2*s[0]*s[0] + 2*s[1]*s[1]) + s[4];
    w01 = s[2]*(-s[1]*s[1]);
    w10 = s[2]*(-s[0]*s[0]);

#   pragma omp parallel for collapse(2) private(i)
    for (j = 0; j < dm[1]; j++) {
        for (i = 0; i < dm[0]; i++) {
            int ij;

            ij = bound(i,j,dm);

            g[ij  ] = w00* f[ij] +
                      w01*(f[bound(i,j-1,dm)  ] + f[bound(i,j+1,dm)  ]) +
                      w10*(f[bound(i-1,j,dm)  ] + f[bound(i+1,j,dm)  ]);
            g[ij+m] = w00* f[ij+m] +
                      w01*(f[bound(i,j-1,dm)+m] + f[bound(i,j+1,dm)+m]) +
                      w10*(f[bound(i-1,j,dm)+m] + f[bound(i+1,j,dm)+m]);
        }
    }
}

void
relax_me(int dm[], double a[], double b[], double s[], int nit, double u[])
{
    int it;
    double w00,w01, w10;
    double *pux, *puy, *pbx, *pby, *paxx, *paxy, *payy;

    w00 = s[2]*(2*s[0]*s[0] + 2*s[1]*s[1]) + s[4];
    w01 = s[2]*(-s[1]*s[1]);
    w10 = s[2]*(-s[0]*s[0]);

#ifdef VERBOSE
    for (it=0; it < 10 - (int)ceil(1.44269504088896*log((double)dm[0])); it++)
        printf("  ");
    printf("%dx%d: ", dm[0],dm[1]);
#endif

    for (it = 0; it < 2*nit; it++) {
        int j, jstart;
        int i, istart;

#ifdef VERBOSE
        printf(" %g", sumsq_me(dm, a, b, s, u));
#endif

        pux  = u;
        puy  = u + dm[0]*dm[1];
        pbx  = b;
        pby  = b + dm[0]*dm[1];
        paxx = a;
        payy = a + dm[0]*dm[1];
        paxy = a + dm[0]*dm[1]*2;

        jstart = it%2;
        for (j = 0; j < dm[1]; j++) {

            istart = (jstart == (j%2));
            for (i = istart; i < dm[0]; i+=2) {
                double sux, suy, axx, ayy, axy, idt;
                int ij, jm1, jp1, im1, ip1;

                im1 = bound(i-1,j  ,dm);
    	        ip1 = bound(i+1,j  ,dm);
                jm1 = bound(i  ,j-1,dm);
    	        jp1 = bound(i  ,j+1,dm);
    	        ij = bound(i,j,dm);

                sux = pbx[ij] -
                      (w01*(pux[jm1] + pux[jp1]) +
                       w10*(pux[im1] + pux[ip1]));
                suy = pby[ij] -
                      (w01*(puy[jm1] + puy[jp1]) +
                       w10*(puy[im1] + puy[ip1]));

                /*
                   syms axx ayy axy sux suy
                   A = [axx axy; axy ayy];
                   su = [sux ; suy]
                   inv(A)*su
                */
                axx = paxx[ij] + w00;
                ayy = payy[ij] + w00;
                axy = paxy[ij];
                idt = 1.0/(axx*ayy*S - axy*axy);

                pux[ij] = idt*( ayy*sux - axy*suy);
                puy[ij] = idt*(-axy*sux + axx*suy);
            }
        }
    }
#ifdef VERBOSE
    printf(" %g\n", sumsq_me(dm, a, b, s, u));
#endif
}


void
Atimesp_me(int dm[], double A[], double s[], double p[], double Ap[])
{
    int i, m = dm[0]*dm[1];

    LtLf_me(dm, p, s, Ap);

#   pragma omp parallel for
    for (i = 0; i < m; i++) {
        Ap[i  ] += A[i  ]*p[i  ] + A[i+2*m]*p[i+m];
        Ap[i+m] += A[i+m]*p[i+m] + A[i+2*m]*p[i  ];
    }
}


double
sumsq_be(int dm[], double a[], double b[], double s[], double u[])
{
    double w00, w01, w02, w10, w11, w20;
    double ss = 0.0;
    int i, j;

    w00 = s[2]*( 6*s[0]*s[0]*s[0]*s[0] + 6*s[1]*s[1]*s[1]*s[1] +
                 8*s[0]*s[0]*s[1]*s[1]) + s[4];
    w01 = s[2]*(-4*s[0]*s[0]*s[1]*s[1] - 4*s[1]*s[1]*s[1]*s[1]);
    w02 = s[2]*(   s[1]*s[1]*s[1]*s[1]);
    w10 = s[2]*(-4*s[0]*s[0]*s[0]*s[0] - 4*s[0]*s[0]*s[1]*s[1]);
    w11 = s[2]*( 2*s[0]*s[0]*s[1]*s[1]);
    w20 = s[2]*(   s[0]*s[0]*s[0]*s[0]);

    for (j = 0; j < dm[1]; j++) {
        double *pux, *puy, *pbx, *pby, *paxx, *paxy, *payy;
        int jm2, jm1, jp1, jp2, im2, im1, ip1, ip2;

        pux  = u + dm[0]*j;
        puy  = u + dm[0]*(j + dm[1]);
        pbx  = b + dm[0]*j;
        pby  = b + dm[0]*(j + dm[1]);
        paxx = a + dm[0]*j;
        payy = a + dm[0]*(j + dm[1]);
        paxy = a + dm[0]*(j + dm[1]*2);

        jm1 = (BOUNDY(j-1,dm[1]) - j)*dm[0];
        jp1 = (BOUNDY(j+1,dm[1]) - j)*dm[0];
        jm2 = (BOUNDY(j-2,dm[1]) - j)*dm[0];
        jp2 = (BOUNDY(j+2,dm[1]) - j)*dm[0];

        for (i = 0; i < dm[0]; i++) {
            double *px = &pux[i], *py = &puy[i];
            double tmp;

            im1 = BOUNDX(i-1,dm[0]) - i;
            ip1 = BOUNDX(i+1,dm[0]) - i;
            im2 = BOUNDX(i-2,dm[0]) - i;
            ip2 = BOUNDX(i+2,dm[0]) - i;

            tmp = (w00 + paxx[i])*px[0] + paxy[i]*py[0] +
                  w01*(px[    jm1] + px[    jp1]) +
                  w02*(px[    jm2] + px[    jp2]) +
                  w10*(px[im1    ] + px[ip1    ]) +
                  w11*(px[im1+jm1] + px[ip1+jm1] + px[im1+jp1] + px[ip1+jp1]) +
                  w20*(px[im2    ] + px[ip2    ]) -
                  pbx[i];
            ss += tmp*tmp;

            tmp = (w00 + payy[i])*py[0] + paxy[i]*px[0] +
                  w01*(py[    jm1] + py[    jp1]) +
                  w02*(py[    jm2] + py[    jp2]) +
                  w10*(py[im1    ] + py[ip1    ]) +
                  w11*(py[im1+jm1] + py[ip1+jm1] + py[im1+jp1] + py[ip1+jp1]) +
                  w20*(py[im2    ] + py[ip2    ]) -
                  pby[i];
            ss += tmp*tmp;
        }
    }
    return(ss);
}


void
LtLf_be(int dm[], double f[], double s[], double g[])
{
    int i, j, jm2, jm1, jp1, jp2, im2, im1, ip1, ip2;
    double *pgx, *pgy, *pfx, *pfy;
    double w00, w01, w02, w10, w11, w20;
    /*
        syms s1 s2
        K1 = [0 -s1 0; 0 2*s1 0; 0 -s1 0];
        K2 = [0 0 0; -s2 2*s2 -s2; 0 0 0];
        K  = K1+K2
        % L  = conv2(K,K)
        L  = sym(zeros(5))
        for i=1:3,
            for j=1:3,
                L(i-1+1:i+1+1,j-1+1:j+1+1) = L(i-1+1:i+1+1,j-1+1:j+1+1) +
                                             K(i,j)*K;
            end;
        end;
    */

    w00 = s[2]*( 6*s[0]*s[0]*s[0]*s[0] + 6*s[1]*s[1]*s[1]*s[1] +
                 8*s[0]*s[0]*s[1]*s[1]) + s[4];
    w01 = s[2]*(-4*s[0]*s[0]*s[1]*s[1] - 4*s[1]*s[1]*s[1]*s[1]);
    w02 = s[2]*(   s[1]*s[1]*s[1]*s[1]);
    w10 = s[2]*(-4*s[0]*s[0]*s[0]*s[0] - 4*s[0]*s[0]*s[1]*s[1]);
    w11 = s[2]*( 2*s[0]*s[0]*s[1]*s[1]);
    w20 = s[2]*(   s[0]*s[0]*s[0]*s[0]);

#   pragma omp parallel for collapse(2) private(i)
    for (j = 0; j < dm[1]; j++) {
        pgx = g + dm[0]*j;
        pgy = g + dm[0]*(j + dm[1]);
        pfx = f + dm[0]*j;
        pfy = f + dm[0]*(j + dm[1]);

        jm1 = (BOUNDY(j-1,dm[1]) - j)*dm[0];
        jp1 = (BOUNDY(j+1,dm[1]) - j)*dm[0];
        jm2 = (BOUNDY(j-2,dm[1]) - j)*dm[0];
        jp2 = (BOUNDY(j+2,dm[1]) - j)*dm[0];

        for (i = 0; i < dm[0]; i++) {
            double *px = &pfx[i], *py = &pfy[i];

            im1 = BOUNDX(i-1,dm[0]) - i;
            ip1 = BOUNDX(i+1,dm[0]) - i;
            im2 = BOUNDX(i-2,dm[0]) - i;
            ip2 = BOUNDX(i+2,dm[0]) - i;

            pgx[i] = w00* px[0] +
                     w01*(px[    jm1] + px[    jp1]) +
                     w02*(px[    jm2] + px[    jp2]) +
                     w10*(px[im1    ] + px[ip1    ]) +
                     w11*(px[im1+jm1] + px[ip1+jm1] +
                          px[im1+jp1] + px[ip1+jp1]) +
                     w20*(px[im2    ] + px[ip2    ]);
            pgy[i] = w00* py[0] +
                     w01*(py[    jm1] + py[    jp1]) +
                     w02*(py[    jm2] + py[    jp2]) +
                     w10*(py[im1    ] + py[ip1    ]) +
                     w11*(py[im1+jm1] + py[ip1+jm1] +
                          py[im1+jp1] + py[ip1+jp1]) +
                     w20*(py[im2    ] + py[ip2    ]);
        }
    }
}


void
relax_be(int dm[], double a[], double b[], double s[], int nit, double u[])
{
    int it;
    double w00, w01, w02, w10, w11, w20;
    double lam;

    w00 = s[2]*( 6*s[0]*s[0]*s[0]*s[0] + 6*s[1]*s[1]*s[1]*s[1] +
                 8*s[0]*s[0]*s[1]*s[1]) + s[4];
    w01 = s[2]*(-4*s[0]*s[0]*s[1]*s[1] - 4*s[1]*s[1]*s[1]*s[1]);
    w02 = s[2]*(   s[1]*s[1]*s[1]*s[1]);
    w10 = s[2]*(-4*s[0]*s[0]*s[0]*s[0] - 4*s[0]*s[0]*s[1]*s[1]);
    w11 = s[2]*( 2*s[0]*s[0]*s[1]*s[1]);
    w20 = s[2]*(   s[0]*s[0]*s[0]*s[0]);
    lam = 2*(w20 + w02) - 2*(w10 + w01) + 4*w11 - w00;
    if (lam < 0.0) lam = 0.0;

#ifdef VERBOSE
    for(it=0; it< 10-(int)ceil(1.44269504088896*log((double)dm[0])); it++)
        printf("  ");
    printf("%dx%d (%g): ", dm[0],dm[1],lam);
#endif

    for (it = 0; it < 4*nit; it++) {
        int j, jstart,jend,jskip;
        int i, istart,iend,iskip;

#ifdef VERBOSE
        printf(" %g", sumsq_be(dm, a, b, s, u)); 
#endif

        if ((it/2)%2) {
            jstart = 0;
            jend   = dm[1];
            jskip  = 1;
        } else {
            jstart = dm[1]-1;
            jend   = -1;
            jskip  = -1;
        }
        if (it%2) {
            istart = 0;
            iend   = dm[0];
            iskip  = 1;
        } else {
            istart = dm[0]-1;
            iend   = -1;
            iskip  = -1;
        }

        for (j = jstart; j != jend; j+=jskip) {
            double *pux, *puy, *pbx, *pby, *paxx, *paxy, *payy;
            int jm2, jm1, jp1, jp2, im2, im1, ip1, ip2;

            pux  = u + dm[0]*j;
            puy  = u + dm[0]*(j + dm[1]);
            pbx  = b + dm[0]*j;
            pby  = b + dm[0]*(j + dm[1]);
            paxx = a + dm[0]*j;
            payy = a + dm[0]*(j + dm[1]);
            paxy = a + dm[0]*(j + dm[1]*2);

            jm1 = (BOUNDY(j-1,dm[1]) - j)*dm[0];
            jp1 = (BOUNDY(j+1,dm[1]) - j)*dm[0];
            jm2 = (BOUNDY(j-2,dm[1]) - j)*dm[0];
            jp2 = (BOUNDY(j+2,dm[1]) - j)*dm[0];

            for (i = istart; i != iend; i+=iskip) {
                double sux, suy, axx, ayy, axy, idt;
                double *px = &pux[i], *py = &puy[i];

                im1 = BOUNDX(i-1,dm[0]) - i;
                ip1 = BOUNDX(i+1,dm[0]) - i;
                im2 = BOUNDX(i-2,dm[0]) - i;
                ip2 = BOUNDX(i+2,dm[0]) - i;

                sux = pbx[i] - ((w00 + paxx[i])*px[0] + paxy[i]*py[0] +
                                w01*(px[    jm1] + px[    jp1]) +
                                w02*(px[    jm2] + px[    jp2]) +
                                w10*(px[im1    ] + px[ip1    ]) +
                                w11*(px[im1+jm1] + px[ip1+jm1] +
                                     px[im1+jp1] + px[ip1+jp1]) +
                                w20*(px[im2    ] + px[ip2    ]));
                suy = pby[i] - (paxy[i]*px[0] + (w00 + payy[i])*py[0] +
                                w01*(py[    jm1] + py[    jp1]) +
                                w02*(py[    jm2] + py[    jp2]) +
                                w10*(py[im1    ] + py[ip1    ]) +
                                w11*(py[im1+jm1] + py[ip1+jm1] +
                                     py[im1+jp1] + py[ip1+jp1]) +
                                w20*(py[im2    ] + py[ip2    ]));

                /*
                   syms axx ayy axy sux suy
                   A = [axx axy; axy ayy];
                   su = [sux ; suy]
                   inv(A)*su
                */
                axx  = paxx[i] + w00 + lam;
                ayy  = payy[i] + w00 + lam;
                axy  = paxy[i];
                idt  = 1.0/(axx*ayy - axy*axy);

                *px += idt*( ayy*sux - axy*suy);
                *py += idt*(-axy*sux + axx*suy);
            }
        }
    }
#ifdef VERBOSE
    printf(" %g\n", sumsq_be(dm, a, b, s, u)); 
#endif
}


void
Atimesp_be(int dm[], double A[], double param[], double p[], double Ap[])
{
    int i, m = dm[0]*dm[1];
    LtLf_be(dm, p, param, Ap);

    for (i = 0; i < m; i++) {
        Ap[i  ] = Ap[i  ] + A[i  ]*p[i  ] + A[i+2*m]*p[i+m];
        Ap[i+m] = Ap[i+m] + A[i+m]*p[i+m] + A[i+2*m]*p[i  ];
    }
}


void
solve22(double a[], double b[], double t,  double u[])
{
    double dt;
    double a0 = a[0]+t, a1 = a[1]+t;

    dt  = a0*a1*S - a[2]*a[2];
    if (dt < 1e-12*a0*a1)
        dt = 1e-12*a0*a1;
    u[0] = (   a1*b[0] - a[2]*b[1])/dt;
    u[1] = (-a[2]*b[0] +   a0*b[1])/dt;
}


double
dotprod(int m, double a[], double b[])
{
    int i;
    double dp = 0.0;

    for (i = 0; i < m; i++)
        dp += a[i]*b[i];

    return(dp);
}


void
addscaled(int m, double a[], double b[], double s)
{
    int i;

    for(i = 0; i < m; i++)
        a[i] += s*b[i];
}


double
norm(int m, double a[])
{
    int i;
    double dp = 0.0;

    for (i = 0; i < m; i++)
        dp += a[i]*a[i];

    return(sqrt(dp));
}


/*
% Solve A*x = b by the conjugate gradient method
% See Gilbert, Moler & Schreiber (1991)
% Sparse Matrices in Matlab: Design and Implementation
% SIAM Journal on Matrix Analysis and Applications
% http://citeseer.ist.psu.edu/gilbert91sparse.html

if nargin<3, tol = 1e-4; end;
if nargin<4, nit = 1000; end;

x    = zeros(size(b));
r    = b;
rtr  = r'*r;
p    = zeros(size(b));
beta = 0;
it   = 0;
while norm(r) > tol*norm(b),
    p      = r + beta*p;
    Ap     = A*p;
    alpha  = rtr/(p'*Ap);
    x      = x + alpha*p;
    r      = r - alpha*Ap;
    rtrold = rtr;
    rtr    = r'*r;
    beta   = rtr/rtrold;

    it = it+1;
    if it>nit, break; end;
end;
*/
void cgs2(int dm[], double A[], double b[], int rtype, double param[],
          double tol, int nit, double x[], double r[], double p[], double Ap[])
{
    int i, m = dm[0]*dm[1]*2, it;
    double rtr, nb, rtrold, alpha, beta;
    void (*Atimesp)();

    /* printf("\n **** %dx%d ****\n",dm[0],dm[1]); */
    if (rtype == 1)
        Atimesp = Atimesp_me;
    else
        Atimesp = Atimesp_be;

    nb      = tol*norm(m,b);

    if (0) {
        /* Assuming starting estimates of zeros */
        /* x    = zeros(size(b)); */
        for (i = 0; i < m; i++)
            x[i] = 0.0;

        /* r    = b; */
        for (i = 0; i < m; i++)
            r[i] = b[i];
    } else {
        /* Assume starting estimates are passed as arguments */
        /* r    = b-A*x; */
        Atimesp(dm, A, param, x, Ap);
        for (i = 0; i < m; i++)
            r[i] = b[i]-Ap[i];
    }

    /* rtr  = r'*r; */
    rtr     = dotprod(m, r, r);

    /* p    = zeros(size(b)); */
    for (i = 0; i < m; i++)
        p[i] = 0.0;

    /* beta = 0; */
    beta    = 0.0;

    /* for it=1:nit, */
    for (it = 0; it < nit; it++) {
        /* if norm(r) < tol*norm(b), break; end; */
        if (norm(m,r) < nb)
            break;

        /* p      = r + beta*p; */
        for (i = 0; i < m; i++)
            p[i]  = r[i] + beta*p[i];

        /* Ap     = A*p; */
        Atimesp(dm, A, param, p, Ap);

        /* alpha  = rtr/(p'*Ap); */
        alpha     = rtr/dotprod(m, p, Ap);

        /* x      = x + alpha*p; */
        addscaled(m, x, p, alpha);

        /* r      = r - alpha*Ap; */
        addscaled(m, r, Ap, -alpha);

        /* rtrold = rtr; */
        rtrold = rtr;

        /* rtr    = r'*r; */
        rtr       = dotprod(m, r, r);

        /* beta   = rtr/rtrold; */
        beta      = rtr/rtrold;

        printf("%d\t%g\t%g  %g %g\n",it, norm(m,r), nb/tol, alpha, beta);
    /* end; */
    }
    /* printf("Done after %d iterations (%g, %g).\n",it, norm(m,r), nb); */
}

/*******************************************************/


double
wt2(double x)
{
        x = fabs(x);
        if (x < 0.5)
                return(0.75 - x*x);
        if (x < 1.5) {
                x = 1.5 - x;
                return(0.5*x*x);
        }
        return(0.0);
}

void
resize(int na[], double *a, int nc[], double *c, double *b)
{
    int i, j, o,om,op;
    double loc, s, w, wm, wp;
    double *ap, *bp, *cp;
    /* a - na[0]*na[1]
     * c - nc[0]*nc[1]
     * b - na[0]*nc[1]
     */

    s = (double) na[1] / (double) nc[1];
    for (j = 0; j < nc[1]; j++) {
        loc = (j + 0.5)*s - 0.5;
        o   = floor(loc + 0.5);
        w   = wt2( o   -loc);
        wp  = wt2((o+1)-loc);
        wm  = wt2((o-1)-loc);
        for (i = 0, bp = b+j*na[0]; i < na[0]; i++, bp++) {
            *bp = wm*a[bound(i,o-1,na)] +
                  w *a[bound(i,o  ,na)] +
                  wp*a[bound(i,o+1,na)];
        }
    }
    s = (double) na[0] / (double) nc[0];
    for (i = 0; i < nc[0]; i++) {
        loc = (i + 0.5)*s - 0.5;
        o   = floor(loc + 0.5);
        om  = BOUNDX(o-1,na[0]);
        op  = BOUNDX(o+1,na[0]);
        w   = wt2( o   -loc);
        wp  = wt2((o+1)-loc);
        wm  = wt2((o-1)-loc);
        for (j = 0, bp = b, cp = c+i; j < nc[1]; j++, bp+=na[0], cp+=nc[0]) {
            *cp = wm*bp[om] + w*bp[o] + wp*bp[op];
        }
    }
}


void
rescale(int n, double *a, double s)
{
    int i;
    for (i = 0; i < n; i++)
        a[i] *= s;
}


void
restrict2(int n, int na[], double *a, int nc[], double *c, double *b)
{
    int i;
    for (i = 0; i < n; i++) {
        resize(na, a + i*na[0]*na[1], nc, c + i*nc[0]*nc[1], b);
     /* rescale(nc[0]*nc[1], c+i*nc[0]*nc[1], 4.0); */
    }
}


void
prolong(int n, int na[], double *a, int nc[], double *c, double *b)
{
    int i;
    for (i = 0; i < n; i++)
        resize(na, a + i*na[0]*na[1], nc, c + i*nc[0]*nc[1], b);
}


void
zeros(int n, double *a)
{
    int i;
    for (i = 0; i < n; i++)
        a[i] = 0.0;
}


void
copy(int n, double *a, double *b)
{
    int i;
    for (i = 0; i < n; i++)
        b[i] = a[i];
}


void
addto(int n, double *a, double *b)
{
    int i;
    for (i = 0; i < n; i++)
        a[i] += b[i];
}


int
fmg2_scratchsize(int n0[])
{
    int    n[32][2], m[32], bs, j;
    bs = 0;
    n[0][0] = n0[0];
    n[0][1] = n0[1];

    for (j = 1; j < 16; j++) {
        n[j][0] = ceil(n[j-1][0]/2.0);
        n[j][1] = ceil(n[j-1][1]/2.0);
        m[j]    = n[j][0]*n[j][1];
        bs += m[j];
        if (n[j][0] < 2 && n[j][1] < 2)
            break;
    }
    return((2*n0[0]*n0[1] + n[0][0]*n[1][1] + 9*bs));
}

/* Full Multigrid solver.  See Numerical Recipes (2nd ed) for more info */
void
fmg2(int n0[], double *a0, double *b0, int rtype, double param0[], int c,
     int nit, double *u0, double *scratch)
{
    int i, j, ng, bs;
    int    n[32][2], m[32];
    double *bo[32], *a[32], *b[32], *u[32], *res, *rbuf, param[32][5];
    void (*relax)(), (*Atimesp)();

    if (rtype == 0) {
        relax   = relax_le;
        Atimesp = Atimesp_le;
    } else if (rtype == 1) {
        relax   = relax_me;
        Atimesp = Atimesp_me;
    } else {
        relax   = relax_be;
        Atimesp = Atimesp_be;
    }

    bo[0]   = b0;
    b[0]    = b0;
    u[0]    = u0;
    a[0]    = a0;
    n[0][0] = n0[0];
    n[0][1] = n0[1];
    m[0]    = n0[0]*n0[1];
    param[0][0] = param0[0];
    param[0][1] = param0[1];
    param[0][2] = param0[2];
    param[0][3] = param0[3];
    param[0][4] = param0[4];

    ng = 1;
    bs = 0;
    for (j = 1; j < 16; j++) {
        n[j][0] = ceil(n[j-1][0]/2.0);
        n[j][1] = ceil(n[j-1][1]/2.0);
        m[j]    = n[j][0]*n[j][1];
        ng ++;
        bs += m[j];
        if (n[j][0] < 2 && n[j][1] < 2)
            break;
    }

    res    = scratch;
    rbuf   = scratch + 2*m[0];
    bo[1]  = scratch + 2*m[0] + n[0][0]*n[1][1];
    b[1]   = scratch + 2*m[0] + n[0][0]*n[1][1] + 2*bs;
    u[1]   = scratch + 2*m[0] + n[0][0]*n[1][1] + 4*bs;
    a[1]   = scratch + 2*m[0] + n[0][0]*n[1][1] + 6*bs;

    for (j = 2; j < ng; j++) {
        bo[j] = bo[j-1] + 2*m[j-1];
        b[j]  =  b[j-1] + 2*m[j-1];
        u[j]  =  u[j-1] + 2*m[j-1];
        a[j]  =  a[j-1] + 3*m[j-1];
    }

    for (j = 1; j < ng; j++) {
        param[j][0] = param0[0]*(double)n[j][0]/n0[0];
        param[j][1] = param0[1]*(double)n[j][1]/n0[1];
        param[j][2] = param[0][2];
        param[j][3] = param[0][3];
        param[j][4] = param[0][4];

        restrict2(2, n[j-1], bo[j-1], n[j], bo[j], rbuf);
        restrict2(3, n[j-1], a[j-1],  n[j], a[j],  rbuf);
    }

    solve22(a[ng-1], bo[ng-1], param0[4], u[ng-1]);

    for (j = ng-2; j >= 0; j--) {
        int jc;

        prolong(2, n[j+1], u[j+1], n[j], u[j], rbuf);
        if (j > 0) copy(2*m[j], bo[j], b[j]);
        for (jc = 0; jc < c; jc++) {
            int jj;
            for (jj = j; jj < ng-1; jj++) {
                relax(n[jj], a[jj], b[jj], param[jj], nit, u[jj]);
                Atimesp(n[jj], a[jj], param[jj], u[jj], res);
                for (i = 0; i < 2*m[jj]; i++)
                    res[i] = b[jj][i] - res[i];

                restrict2(2, n[jj], res, n[jj+1], b[jj+1], rbuf);
                zeros(2*m[jj+1], u[jj+1]);
            }
            solve22(a[ng-1], b[ng-1], param0[4], u[ng-1]);
            for (jj = ng-2; jj >= j; jj--) {
                prolong(2, n[jj+1], u[jj+1], n[jj], res, rbuf);
                addto(2*m[jj], u[jj], res);
                relax(n[jj], a[jj], b[jj], param[jj], nit, u[jj]);
            }
        }
    }
}

