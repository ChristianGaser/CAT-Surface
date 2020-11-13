/* $Id$ 
   largely modified version of diffeo2d.c from John Ashburner 
   (c) John Ashburner (2007) */
   
#include <math.h>
#include <stdio.h>
#include "dartel.h"
#include "CAT_Surf.h"
#include "CAT_Interpolate.h"
#ifdef _OPENMP
#include <omp.h>
#endif

extern double floor();


void
composition(int dm[], double *A, double *B, double *C)
{
    double *Ax, *Ay, *Bx, *By, *Cx, *Cy;
    int i, m = dm[0], n = dm[1], mm = m*n;

    Ax = A;
    Ay = &A[mm];
    Bx = B;
    By = &B[mm];
    Cx = C;
    Cy = &C[mm];

    for (i = 0; i < mm; i++) {
        double x, y;
        double k11,k12,k21,k22;
        double dx1, dx2, dy1, dy2;
        int ix, iy;

        x     = Ax[i]-1.0;
        y     = Ay[i]-1.0;
        ix    = (int)floor(x); dx1=x-ix; dx2=1.0-dx1;
        iy    = (int)floor(y); dy1=y-iy; dy2=1.0-dy1;

        k22   = Bx[bound(ix,  iy,  dm)] - 1.0;
        k12   = Bx[bound(ix+1,iy,  dm)] - 1.0;
        k21   = Bx[bound(ix,  iy+1,dm)] - 1.0;
        k11   = Bx[bound(ix+1,iy+1,dm)] - 1.0;
        k12   = k12 - floor((k12-k22)/m+0.5)*m;
        k21   = k21 - floor((k21-k22)/m+0.5)*m;
        k11   = k11 - floor((k11-k22)/m+0.5)*m;
        Cx[i] = (k22*dx2 + k12*dx1)*dy2 + (k21*dx2 + k11*dx1)*dy1 + 1.0;

        k22   = By[bound(ix,  iy,  dm)] - 1.0;
        k12   = By[bound(ix+1,iy,  dm)] - 1.0;
        k21   = By[bound(ix,  iy+1,dm)] - 1.0;
        k11   = By[bound(ix+1,iy+1,dm)] - 1.0;
        k12   = k12 - floor((k12-k22)/n+0.5)*n;
        k21   = k21 - floor((k21-k22)/n+0.5)*n;
        k11   = k11 - floor((k11-k22)/n+0.5)*n;
        Cy[i] = (k22*dx2 + k12*dx1)*dy2 + (k21*dx2 + k11*dx1)*dy1 + 1.0;
    }
}


void
composition_jacobian(int dm[],
                     double *A, double *JA, double *B, double *JB,
                     double *C, double *JC)
{
    double *Ax, *Ay, *JA00, *JA01, *JA10, *JA11;
    double *Bx, *By, *JB00, *JB01, *JB10, *JB11, jb[2][2];
    double *Cx, *Cy, *JC00, *JC01, *JC10, *JC11;
    int i, m = dm[0], n = dm[1], mm = m*n;
    Ax   =   A;
    Ay   =  &A[mm];
    JA00 = &JA[mm*0];
    JA01 = &JA[mm*1];
    JA10 = &JA[mm*2];
    JA11 = &JA[mm*3];

    Bx   =   B;
    By   =  &B[mm];
    JB00 = &JB[mm*0];
    JB01 = &JB[mm*1];
    JB10 = &JB[mm*2];
    JB11 = &JB[mm*3];

    Cx   =   C;
    Cy   =  &C[mm];
    JC00 = &JC[mm*0];
    JC01 = &JC[mm*1];
    JC10 = &JC[mm*2];
    JC11 = &JC[mm*3];

    for (i = 0; i < mm; i++) {
        double x, y;
        double k11,k12,k21,k22;
        double dx1, dx2, dy1, dy2;
        int ix, iy, o11,o12,o21,o22;

        x    = Ax[i]-1.0;
        y    = Ay[i]-1.0;
        ix   = (int)floor(x); dx1=x-ix; dx2=1.0-dx1;
        iy   = (int)floor(y); dy1=y-iy; dy2=1.0-dy1;

        o22   = bound(ix,  iy,  dm);
        o12   = bound(ix+1,iy,  dm);
        o21   = bound(ix,  iy+1,dm);
        o11   = bound(ix+1,iy+1,dm);

        k22  = Bx[o22]-1.0;
        k12  = Bx[o12]-1.0;
        k21  = Bx[o21]-1.0;
        k11  = Bx[o11]-1.0;
        k12  = k12-floor((k12-k22)/m+0.5)*m;
        k21  = k21-floor((k21-k22)/m+0.5)*m;
        k11  = k11-floor((k11-k22)/m+0.5)*m;
        Cx[i]= (k22*dx2 + k12*dx1)*dy2 + (k21*dx2 + k11*dx1)*dy1 + 1.0;

        k22  = By[o22]-1.0;
        k12  = By[o12]-1.0;
        k21  = By[o21]-1.0;
        k11  = By[o11]-1.0;
        k12  = k12-floor((k12-k22)/n+0.5)*n;
        k21  = k21-floor((k21-k22)/n+0.5)*n;
        k11  = k11-floor((k11-k22)/n+0.5)*n;
        Cy[i]= (k22*dx2 + k12*dx1)*dy2 + (k21*dx2 + k11*dx1)*dy1 + 1.0;

        k22  = JB00[o22];
        k12  = JB00[o12];
        k21  = JB00[o21];
        k11  = JB00[o11];
        jb[0][0] =  (k22*dx2 + k12*dx1)*dy2 + (k21*dx2 + k11*dx1)*dy1;

        k22  = JB01[o22];
        k12  = JB01[o12];
        k21  = JB01[o21];
        k11  = JB01[o11];
        jb[0][1] =  (k22*dx2 + k12*dx1)*dy2 + (k21*dx2 + k11*dx1)*dy1;

        k22  = JB10[o22];
        k12  = JB10[o12];
        k21  = JB10[o21];
        k11  = JB10[o11];
        jb[1][0] =  (k22*dx2 + k12*dx1)*dy2 + (k21*dx2 + k11*dx1)*dy1;

        k22  = JB11[o22];
        k12  = JB11[o12];
        k21  = JB11[o21];
        k11  = JB11[o11];
        jb[1][1] =  (k22*dx2 + k12*dx1)*dy2 + (k21*dx2 + k11*dx1)*dy1;

        JC00[i] = jb[0][0]*JA00[i] + jb[1][0]*JA01[i];
        JC01[i] = jb[0][1]*JA00[i] + jb[1][1]*JA01[i];
        JC10[i] = jb[0][0]*JA10[i] + jb[1][0]*JA11[i];
        JC11[i] = jb[0][1]*JA10[i] + jb[1][1]*JA11[i];
    }
}


void
composition_detjac(int dm[],
                     double *A, double *dA, double *B, double *dB,
                     double *C, double *dC)
{
    double *Ax, *Ay;
    double *Bx, *By, jb;
    double *Cx, *Cy;
    int i, m = dm[0], n = dm[1], mm = m*n;
    Ax   =   A;
    Ay   =  &A[mm];

    Bx   =   B;
    By   =  &B[mm];

    Cx   =   C;
    Cy   =  &C[mm];

    for (i = 0; i < mm; i++) {
        double x, y;
        double k11,k12,k21,k22;
        double dx1, dx2, dy1, dy2;
        int ix, iy, o11,o12,o21,o22;

        x    = Ax[i]-1.0;
        y    = Ay[i]-1.0;
        ix   = (int)floor(x); dx1=x-ix; dx2=1.0-dx1;
        iy   = (int)floor(y); dy1=y-iy; dy2=1.0-dy1;
        
        o22   = bound(ix,  iy,  dm);
        o12   = bound(ix+1,iy,  dm);
        o21   = bound(ix,  iy+1,dm);
        o11   = bound(ix+1,iy+1,dm);

        k22  = Bx[o22]-1.0;
        k12  = Bx[o12]-1.0;
        k21  = Bx[o21]-1.0;
        k11  = Bx[o11]-1.0;
        k12  = k12-floor((k12-k22)/m+0.5)*m;
        k21  = k21-floor((k21-k22)/m+0.5)*m;
        k11  = k11-floor((k11-k22)/m+0.5)*m;
        Cx[i]= (k22*dx2 + k12*dx1)*dy2 + (k21*dx2 + k11*dx1)*dy1 + 1.0;

        k22  = By[o22]-1.0;
        k12  = By[o12]-1.0;
        k21  = By[o21]-1.0;
        k11  = By[o11]-1.0;
        k12  = k12-floor((k12-k22)/n+0.5)*n;
        k21  = k21-floor((k21-k22)/n+0.5)*n;
        k11  = k11-floor((k11-k22)/n+0.5)*n;
        Cy[i]= (k22*dx2 + k12*dx1)*dy2 + (k21*dx2 + k11*dx1)*dy1 + 1.0;

        k22  = dB[o22];
        k12  = dB[o12];
        k21  = dB[o21];
        k11  = dB[o11];
        jb   =  (k22*dx2 + k12*dx1)*dy2 + (k21*dx2 + k11*dx1)*dy1;

        dC[i] = jb*dA[i];
    }
}


double
samp(int dm[], double f[], double x, double y)
{
    int ix, iy, o11, o12, o21, o22;
    double dx1, dx2, dy1, dy2;

    ix   = (int)floor(x); dx1=x-ix; dx2=1.0-dx1;
    iy   = (int)floor(y); dy1=y-iy; dy2=1.0-dy1;

    o22   = bound(ix,  iy,  dm);
    o12   = bound(ix+1,iy,  dm);
    o21   = bound(ix,  iy+1,dm);
    o11   = bound(ix+1,iy+1,dm);

    return((f[o22]*dx2 + f[o12]*dx1)*dy2 + (f[o21]*dx2 + f[o11]*dx1)*dy1);
}


void
expdef(int dm[], int k, double v[], double t0[], double t1[],
       double J0[], double J1[])
{
    double *optr;
    double td;
    int m = dm[0]*dm[1];
    int i, j;

    optr = t0;

    td = 1;
    for (i = 0; i < k; i++)
        td = td*2;
    td = 1.0/td;

    if (J0 != (double *) 0) {
        for (j = 0; j < dm[1]; j++) {
            for (i = 0; i < dm[0]; i++) {
                int o   = i + dm[0]*j;
                t0[o  ] = (i+1) + v[o  ]*td;
                t0[o+m] = (j+1) + v[o+m]*td;
                J0[o    ] = (v[bound(i+1,j  ,dm)  ] -
                             v[bound(i-1,j  ,dm)  ])*td/2 + 1.0;
                J0[o+  m] = (v[bound(i+1,j  ,dm)+m] -
                             v[bound(i-1,j  ,dm)+m])*td/2;
                J0[o+2*m] = (v[bound(i  ,j+1,dm)  ] -
                             v[bound(i  ,j-1,dm)  ])*td/2;
                J0[o+3*m] = (v[bound(i  ,j+1,dm)+m] -
                             v[bound(i  ,j-1,dm)+m])*td/2 + 1.0;
            }
        }
        for (i = 0; i < k; i++) {
            double *tmpp;
            composition_jacobian(dm, t0, J0, t0, J0, t1, J1);
            tmpp = t0; t0   = t1; t1   = tmpp;
            tmpp = J0; J0   = J1; J1   = tmpp;
        }
    } else {
        for (j = 0; j < dm[1]; j++) {
            for (i = 0; i < dm[0]; i++) {
                int o   = i + dm[0]*j;
                t0[o  ] = (i+1) + v[o  ]*td;
                t0[o+m] = (j+1) + v[o+m]*td;
            }
        }
        for (i = 0; i < k; i++) {
            double *tmpp;
            composition(dm, t0, t0, t1);
            tmpp = t0; t0   = t1; t1   = tmpp;
        }
    }

    if (optr != t0) {
        for (i = 0; i < 2*m; i++)
            t1[i] = t0[i];

        if (J0 != (double *) 0)
            for (i = 0; i < 4*m; i++)
                J1[i] = J0[i];
    }
}


void
expdefdet(int dm[], int k, double v[], double t0[], double t1[],
          double J0[], double J1[])
{
    double *optr;
    double td;
    int m = dm[0]*dm[1];
    int i, j;

    optr = t0;

    td = 1;
    for (i = 0; i < k; i++)
        td = td*2;
    td = 1.0/td;

    if (J0 != (double *) 0) {
        for (j = 0; j < dm[1]; j++) {
            for (i = 0; i < dm[0]; i++) {
                double j00, j01, j10, j11;
                int o   = i + dm[0]*j;
                t0[o  ] = (i+1) + v[o  ]*td;
                t0[o+m] = (j+1) + v[o+m]*td;

                j00 = (v[bound(i+1,j  ,dm)  ] -
                       v[bound(i-1,j  ,dm)  ])*td/2 + 1.0;
                j10 = (v[bound(i+1,j  ,dm)+m] -
                       v[bound(i-1,j  ,dm)+m])*td/2;
                j01 = (v[bound(i  ,j+1,dm)  ] -
                       v[bound(i  ,j-1,dm)  ])*td/2;
                j11 = (v[bound(i  ,j+1,dm)+m] -
                       v[bound(i  ,j-1,dm)+m])*td/2 + 1.0;

                J0[o]   = j00*j11 - j10*j01;
            }
        }
        for (i = 0; i < k; i++) {
            double *tmpp;
            composition_detjac(dm, t0, J0, t0, J0, t1, J1);
            tmpp = t0; t0   = t1; t1   = tmpp;
            tmpp = J0; J0   = J1; J1   = tmpp;
        }
    } else {
        for (j = 0; j < dm[1]; j++) {
            for (i = 0; i < dm[0]; i++) {
                int o   = i + dm[0]*j;
                t0[o  ] = (i+1) + v[o  ]*td;
                t0[o+m] = (j+1) + v[o+m]*td;
            }
        }
        for (i = 0; i < k; i++) {
            double *tmpp;
            composition(dm, t0, t0, t1);
            tmpp = t0; t0   = t1; t1   = tmpp;
        }
    }

    if (optr != t0) {
        for (i = 0; i < 2*m; i++)
            t1[i] = t0[i];

        if (J0 != (double *) 0)
            for (i = 0; i < 4*m; i++)
                J1[i] = J0[i];
    }
}


int
pow2(int k)
{
    int j0, td = 1;

    for (j0 = 0; j0 < k; j0++)
        td = td*2;
    return(td);
}


/*
 * J0 := J0*inv(I+diag(v0)*sc)
 */
void
jac_div_smalldef(int dm[], double sc, double v0[], double J0[])
{
    int j0, j1;
    int m = dm[0]*dm[1];
    double sc2 = sc/2.0;
    double *v1 = v0+m;

    for (j1 = 0; j1 < dm[1]; j1++) {
        for (j0 = 0; j0 < dm[0]; j0++) {
            int o, om1, op1;
            double j00, j01, j10, j11;
            double t00, t01, t10, t11;
            double idt;

            om1 = bound(j0-1,j1,dm);
            op1 = bound(j0+1,j1,dm);
            j00 = (v0[op1] - v0[om1])*sc2 + 1.0;
            j01 = (v1[op1] - v1[om1])*sc2;

            om1 = bound(j0,j1-1,dm);
            op1 = bound(j0,j1+1,dm);
            j10 = (v0[op1] - v0[om1])*sc2;
            j11 = (v1[op1] - v1[om1])*sc2 + 1.0;

            /*
            syms j00 j01 j10 j11
            syms t00 t01 t10 t11
            J1 = [j00 j01; j10 j11];
            J0 = [t00 t01; t10 t11];
            inv(J1)
            J1*J0
            */
            idt = 1.0/(j00*j11-j10*j01);
            t00 =  idt*j11;
            t01 = -idt*j10;
            t10 = -idt*j01;
            t11 =  idt*j00;


            o   = j0+dm[0]*j1;
            j00 = J0[o  ]; j01 = J0[o+m*2];
            j10 = J0[o+m]; j11 = J0[o+m*3];
            J0[o    ] = j00*t00 + j01*t10;
            J0[o+m  ] = j10*t00 + j11*t10;
            J0[o+m*2] = j00*t01 + j01*t11;
            J0[o+m*3] = j10*t01 + j11*t11;
        }
    }
}


double
initialise_objfun(int dm[], double f[], double g[], double t0[], double J0[],
                  double dj[], double b[], double A[])
{
    int j, m = dm[0]*dm[1];
    double ssl = 0.0, dt = 1.0;

    for (j = 0; j < m; j++) {
        double x, y;
        int    ix, iy;
        double k11, k12, k21, k22;
        int    o11, o12, o21, o22;
        double dx0, dx1, dx2, dy0, dy1, dy2;
        double d, dx, dy;

        x    = t0[j  ]-1.0;
        y    = t0[j+m]-1.0;
        ix   = (int)floor(x); dx1=x-ix; dx2=1.0-dx1;
        iy   = (int)floor(y); dy1=y-iy; dy2=1.0-dy1;

        o22   = bound(ix,  iy,  dm);
        o12   = bound(ix+1,iy,  dm);
        o21   = bound(ix,  iy+1,dm);
        o11   = bound(ix+1,iy+1,dm);

        k22  = f[o22];
        k12  = f[o12];
        k21  = f[o21];
        k11  = f[o11];

        d    = ((k11*dx1 + k21*dx2)*dy1 + (k12*dx1 + k22*dx2)*dy2) - g[j];
        dx0  = ((k11     - k21    )*dy1 + (k12     - k22    )*dy2);
        dy0  = ((k11*dx1 + k21*dx2)     - (k12*dx1 + k22*dx2)    );

        dx   = J0[j    ]*dx0 + J0[j+  m]*dy0;
        dy   = J0[j+2*m]*dx0 + J0[j+3*m]*dy0;

        if (dj != (double *)0)
            dt = dj[j];

        A[j    ] = dx*dx*dt;
        A[j+  m] = dy*dy*dt;
        A[j+2*m] = dx*dy*dt;

        b[j  ]   = dx*d*dt;
        b[j+m]   = dy*d*dt;

        ssl += d*d*dt;
    }
    return(0.5*ssl);
}

double
initialise_objfun_mn(int dm[], double f[], double g[], double t0[],
                     double J0[], double jd[], double b[], double A[])
{
    int j, m = dm[0]*dm[1];
    double ssl = 0.0;
 
    for (j = 0; j < m; j++) {
        double x, y;
        int    ix, iy, k;
        double dx1, dx2, dy1, dy2;
        double k22, k12, k21, k11;
        int    o11, o12, o21, o22;
        double dx0, dy0;
        double dx[128], dy[128], Ya[128], T[128], sT = 1.0, sY;
        double ta11, ta22, ta12;
        double tb1,  tb2, tss;
        double sk22 = 1.0, sk12 = 1.0, sk21 = 1.0, sk11 = 1.0;
        
        x    = t0[j    ]-1.0;
        y    = t0[j+m  ]-1.0;
        ix   = (int)floor(x); dx1=x-ix; dx2=1.0-dx1;
        iy   = (int)floor(y); dy1=y-iy; dy2=1.0-dy1;

        o22   = bound(ix,  iy,  dm);
        o12   = bound(ix+1,iy,  dm);
        o21   = bound(ix,  iy+1,dm);
        o11   = bound(ix+1,iy+1,dm);

        sY   = 0.0;
        
        sk22 = sk12 = sk21 = sk11 = 1.0;
        
        for (k = 0; k < dm[2]; k++) {
            int km = k*m;

            T[k]  = g[j + km];
            sT   -= T[k];
            k22  = f[o22 + km]; sk22 -= k22; k22 = LOG(k22);
            k12  = f[o12 + km]; sk12 -= k12; k12 = LOG(k12);
            k21  = f[o21 + km]; sk21 -= k21; k21 = LOG(k21);
            k11  = f[o11 + km]; sk11 -= k11; k11 = LOG(k11);

            Ya[k]  = exp((k22*dx2 + k12*dx1)*dy2 + (k21*dx2 + k11*dx1)*dy1);
            sY   += Ya[k];
            
            dx0   = ((k22     - k12    )*dy2 + (k21     - k11    )*dy1);
            dy0   = ((k22*dx2 + k12*dx1)     - (k21*dx2 + k11*dx1)    );

            dx[k] = -(J0[j    ]*dx0 + J0[j+  m]*dy0);
            dy[k] = -(J0[j+2*m]*dx0 + J0[j+3*m]*dy0);
        }
        
        k    = dm[2];
        T[k] = sT;

        k22 = LOG(sk22);
        k12 = LOG(sk21);
        k21 = LOG(sk12);
        k11 = LOG(sk11);

        Ya[k]  = exp((k22*dx2 + k12*dx1)*dy2 + (k21*dx2 + k11*dx1)*dy1);
        sY   += Ya[k];
            
        dx0   = ((k22     - k12    )*dy2 + (k21     - k11    )*dy1);
        dy0   = ((k22*dx2 + k12*dx1)     - (k21*dx2 + k11*dx1)    );

        dx[k] = -(J0[j    ]*dx0 + J0[j+  m]*dy0);
        dy[k] = -(J0[j+2*m]*dx0 + J0[j+3*m]*dy0);

        ta11 = ta22 = ta12 = 0.0;
        tb1  = tb2  = 0.0;
        tss  = 0.0;
        for (k = 0; k <= dm[2]; k++) {
            double wt;
            int k1;

            Ya[k] /= sY;
            tss  += log(Ya[k])*T[k];
            tb1  += (Ya[k] - T[k])*dx[k];
            tb2  += (Ya[k] - T[k])*dy[k];

            for (k1 = 0; k1 < k; k1++) {
                wt    =  -Ya[k]*Ya[k1];
                ta11 += wt* dx[k]*dx[k1]*2;
                ta22 += wt* dy[k]*dy[k1]*2;
                ta12 += wt*(dx[k]*dy[k1] + dx[k1]*dy[k]);
            }
            wt    = Ya[k]*(1.0 - Ya[k]);
            ta11 += wt*dx[k]*dx[k];
            ta22 += wt*dy[k]*dy[k];
            ta12 += wt*dx[k]*dy[k];
        }
        if (jd != (double *) 0) {
            double dt = jd[j];
            if (dt<0.0) dt = 0.0;
            A[j    ]  = ta11*dt;
            A[j+m  ]  = ta22*dt;
            A[j+m*2]  = ta12*dt;
            b[j    ]  = tb1*dt;
            b[j+m  ]  = tb2*dt;
            ssl      -= tss*dt;
        } else {
            A[j    ] = ta11;
            A[j+m  ] = ta22;
            A[j+m*2] = ta12;
            b[j    ] = tb1;
            b[j+m  ] = tb2;
            ssl     -= tss;
        }
    }
    return(ssl);
}


/*
 * t0 = Id + v0*sc
 * J0 = Id + I+diag(v0)*sc
 */
void
smalldef_jac(int dm[], double sc, double v0[], double t0[], double J0[])
{
    int i, j;
    int m = dm[0]*dm[1];
    double sc2 = sc/2.0;
    double *v1 = v0+m;

    for (j = 0; j < dm[1]; j++) {
        for (i = 0; i < dm[0]; i++) {
            int o   = i + dm[0]*j;
            int om1, op1;
            t0[o  ] = (i+1) + v0[o]*sc;
            t0[o+m] = (j+1) + v1[o]*sc;

            om1 = bound(i-1,j,dm);
            op1 = bound(i+1,j,dm);
            J0[o    ] = (v0[op1] - v0[om1])*sc2 + 1.0;
            J0[o+  m] = (v1[op1] - v1[om1])*sc2;

            om1 = bound(i,j-1,dm);
            op1 = bound(i,j+1,dm);
            J0[o+2*m] = (v0[op1] - v0[om1])*sc2;
            J0[o+3*m] = (v1[op1] - v1[om1])*sc2 + 1.0;
        }
    }
}


void
squaring(int dm[], int k, int save_transf, double b[], double A[],
         double t0[], double t1[], double J0[], double J1[])
{
    int i, j, m = dm[0]*dm[1];
    double *ptr = t0;

    for (i = 0; i < k; i++) {
        double *buf1, *buf2;
        buf1 = t1; /* Re-use some memory */
        buf2 = J1;

        for (j = 0; j < m; j++) {
            double tmp00, tmp01, tmp11, tmp1, tmp2, tmp3, tmp4;
            double x, y;
            double j11, j21, j12, j22, dt;

            x   = t0[j  ]-1.0;
            y   = t0[j+m]-1.0;
            j11 = J0[j  ]; j12 = J0[j+2*m];
            j21 = J0[j+m]; j22 = J0[j+3*m];
            dt  = j11*j22 - j12*j21;

            tmp1      = samp(dm,b  ,x,y);
            tmp2      = samp(dm,b+m,x,y);

            buf1[j  ] = dt*(tmp1*j11 + tmp2*j21);
            buf1[j+m] = dt*(tmp1*j12 + tmp2*j22);

            tmp00 = samp(dm,A    ,x,y);
            tmp11 = samp(dm,A+  m,x,y);
            tmp01 = samp(dm,A+2*m,x,y);

            tmp1  = tmp00*j11+tmp01*j21;
            tmp2  = tmp01*j11+tmp11*j21;
            tmp3  = tmp00*j12+tmp01*j22;
            tmp4  = tmp01*j12+tmp11*j22;

            buf2[j    ] = dt*(tmp1*j11+tmp2*j21);
            buf2[j+  m] = dt*(tmp3*j12+tmp4*j22);
            buf2[j+2*m] = dt*(tmp1*j12+tmp2*j22);
        }
        for (j = 0; j < 2*m; j++) b[j] += buf1[j];
        for (j = 0; j < 3*m; j++) A[j] += buf2[j];

        if (save_transf || (i < k-1)) {
            double *tmpp;
            composition_jacobian(dm, t0, J0, t0, J0, t1, J1);
            tmpp = t0; t0   = t1; t1   = tmpp;
            tmpp = J0; J0   = J1; J1   = tmpp;
        }
    }
    if (save_transf && ptr != t0) {
        for (j = 0; j < m*2; j++) t1[j] = t0[j];
        for (j = 0; j < m*4; j++) J1[j] = J0[j];
    }
}


void
unwrap(int dm[], double f[])
{
    int i, j;

    f[0] = f[0] - floor((f[0] - 0)/dm[0] + 0.5)*dm[0];
    for (j = 0; j < dm[1]; j++) {
        double *pt = &f[j*dm[0]];
        if (j == 0) {
            pt[0] = pt[0] - floor(pt[0]/dm[0] + 0.5)*dm[0];
        } else {
            pt[0] = pt[0] - floor((pt[0] - pt[-dm[0]])/dm[0] + 0.5)*dm[0];
        }
        for (i = 1; i < dm[0]; i++) {
            pt[i] = pt[i] - floor((pt[i] - pt[i-1])/dm[0] + 0.5)*dm[0];
        }
    }
    for (i = 0; i < dm[0]; i++) {
        double *pt = &f[i + dm[0]*dm[1]];
        if (i == 0) {
            pt[0] = pt[0] - floor(pt[0]/dm[1] + 0.5)*dm[1];
        } else {
            pt[0] = pt[0] - floor((pt[0] - pt[-1])/dm[1] + 0.5)*dm[1];
        }
        for (j = dm[0]; j < dm[1]*dm[0]; j+=dm[0]) {
            pt[j] = pt[j] - floor((pt[j] - pt[j-dm[0]])/dm[1] + 0.5)*dm[1];
        }
    }
}


int
dartel_scratchsize(int dm[], int code)
{
    int m1, m2;
    int m = dm[0]*dm[1];

    m1 = 15*m;
    if (code) m1 += 5*m;

    m2 = 5*m+fmg2_scratchsize(dm);
    if (m1 > m2)
        return(m1);
    else
        return(m2);
}

void
dartel(struct dartel_prm prm, int dm[], double v[], double g[], double f[],
       double dj[], double ov[], double ll[], double *buf)
{
    double *sbuf;
    double *b, *A, *b1, *A1;
    double *t0, *t1, *J0, *J1;
    double sc;
    double ssl, ssp;
    double normb, wphi;
    static double param[5] = {1.0,1.0,1.0,0.0,0.0};
    int i, j, k, m = dm[0]*dm[1];

    /*
        Allocate memory.
          0  1  2  3  4  5  6  7  8  9 10 11 12 13 14
        [ A  A  A  t  t  J  J  J  J  t  t  J  J  J  J] for computing derivatives
        [ A  A  A s1 s2 s3 s4 s5 s6 s7 s8] for CGS solver
    */
    b    = ov;
    A    = buf;
    t0   = buf +  3*m;
    J0   = buf +  5*m;
    t1   = buf +  9*m;
    J1   = buf + 11*m;
    A1   = buf + 15*m;
    b1   = buf + 18*m;
    sbuf = buf +  3*m;

    sc = 1.0/pow2(prm.k);

    expdef(dm, prm.k, v, t0, t1, J0, J1);
    
    jac_div_smalldef(dm, sc, v, J0);

    if (prm.code == 2)
        ssl = initialise_objfun_mn(dm, f, g, t0, J0, dj, b, A);
    else
        ssl = initialise_objfun(dm, f, g, t0, J0, dj, b, A);

    smalldef_jac(dm, -sc, v, t0, J0);

    squaring(dm, prm.k, prm.code==1, b, A, t0, t1, J0, J1);

    if (prm.code == 1) {
        jac_div_smalldef(dm, -sc, v, J0);
        ssl += initialise_objfun(dm, g, f, t0, J0, (double *)0, b1, A1);
        smalldef_jac(dm, sc, v, t0, J0);
        squaring(dm, prm.k, 0, b1, A1, t0, t1, J0, J1);
        for (j = 0; j < m*2; j++) b[j] -= b1[j];
        for (j = 0; j < m*3; j++) A[j] += A1[j];
    }

    param[2] = prm.rparam[2];
    param[3] = prm.rparam[3];
    param[4] = prm.rparam[4];

    if (prm.rtype == 0)
        LtLf_le(dm, v, param, t1);
    else if (prm.rtype == 1)
        LtLf_me(dm, v, param, t1);
    else
        LtLf_be(dm, v, param, t1);

    ssp = 0.0;
    for (j = 0; j < 2*m; j++) {
        b[j] = b[j]*sc + t1[j];
        ssp += t1[j]*v[j];
    }
    normb = norm(2*m,b);

    for (j = 0; j < 3*m; j++) A[j] *= sc;

    /* Solve equations for Levenberg-Marquardt update:
       v = v - inv(H + L'*L + R)*(d + L'*L*v)
          v: velocity or flow field
          H: matrix of second derivatives
          L: regularisation (L'*L is the inverse of the prior covariance)
          R: Levenberg-Marquardt regularisation
          d: vector of first derivatives
     */
     
    if (prm.lmreg>0.0) param[4] = param[4] + prm.lmreg;

    fmg2(dm, A, b, prm.rtype, param, prm.cycles, prm.its, sbuf, sbuf+2*m);

    /* update and locally weight velocity field: pole regions are weighted less 
       due to distortions from sphere to cartesian grid 
       See also Cheng et al. (2020) Cortical surface registration using unsupervised learning, Neuroimage */
    k = 0;
    for (j = 0; j < dm[1]; j++) {
        wphi = sin(_PI*(j+1)/(dm[1]));
/* weighting currently not used */
        wphi = 1.0;
        for (i = 0; i < dm[0]; i++) {
          ov[k  ] = v[k  ] - wphi*sbuf[k  ];
          ov[k+m] = v[k+m] - wphi*sbuf[k+m];
          k++;
        }
    }

    ll[0] = ssl;
    ll[1] = ssp*0.5;
    ll[2] = normb;
}

