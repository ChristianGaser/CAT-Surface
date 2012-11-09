/* %Id$ */
/* (c) John Ashburner (2007) */

#include <math.h>
#include <stdio.h>
#include "dartel.h"
#include "CAT_Surf.h"
#include "CAT_Interpolate.h"

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

int
pow2(int k)
{
    int j0, td = 1;

    for (j0 = 0; j0 < k; j0++)
        td = td*2;
    return(td);
}

void
LtLf_me_poly(polygons_struct *sphere, struct dartel_poly *dpoly,
             double f[], double s[], double g[])
{
    int i, mm = sphere->n_points;
    double w00,w01,w10;

printf("LtLf_me_poly\n");
    w00 = s[2] * (2*s[0]*s[0] + 2*s[1]*s[1]) + s[4];
    w01 = s[2] * (-s[1]*s[1]);
    w10 = s[2] * (-s[0]*s[0]);

    for (i = 0; i < mm; i++) {
        double kxn, kxp, kyn, kyp;

        kxn   = interp_point_unit_sphere(sphere, f, dpoly->ntheta[i]);
        kxp   = interp_point_unit_sphere(sphere, f, dpoly->ptheta[i]);
        kyn   = interp_point_unit_sphere(sphere, f, dpoly->nphi[i]);
        kyp   = interp_point_unit_sphere(sphere, f, dpoly->pphi[i]);
        g[i] = w00*f[i] + w01*(kyn + kyp) + w10*(kxn + kxp);

        kxn   = interp_point_unit_sphere(sphere, f+mm, dpoly->ntheta[i]);
        kxp   = interp_point_unit_sphere(sphere, f+mm, dpoly->ptheta[i]);
        kyn   = interp_point_unit_sphere(sphere, f+mm, dpoly->nphi[i]);
        kyp   = interp_point_unit_sphere(sphere, f+mm, dpoly->pphi[i]);
        g[i+mm] = w00*f[i+mm] + w01*(kyn + kyp) + w10*(kxn + kxp);
    }
}

void
Atimesp_me_poly(polygons_struct *sphere, struct dartel_poly *dpoly,
                double A[], double s[], double p[], double Ap[])
{
    int i, mm = sphere->n_points;

printf("Atimesp_me_poly\n");
    LtLf_me_poly(sphere, dpoly, p, s, Ap);
output_values_any_format("Ap0.txt", sphere->n_points, Ap, 1);

    for(i = 0; i < mm; i++) {
        Ap[i   ] += A[i   ]*p[i   ] + A[i+2*mm]*p[i+mm];
        Ap[i+mm] += A[i+mm]*p[i+mm] + A[i+2*mm]*p[i   ];
    }
}


void
cgs2_poly(polygons_struct *sphere, struct dartel_poly *dpoly, double A[],
          double b[], double param[], double tol, int nit,
          double x[], double r[], double p[], double Ap[])
{
    int i, m = sphere->n_points*2, it;
    double rtr, nb, rtrold, alpha, beta;

printf("cgs2_poly\n");
    nb      = tol*norm(m,b);

//for (i = 0; i < m; i++) x[i] = 0.0;

    if (0) {
        /* Assuming starting estimates of zeros */
        /* x    = zeros(size(b)); */
        for (i = 0; i < m; i++)
            x[i] = 0.0;

        /* r    = b; */
        for(i=0; i<m;i++)
            r[i] = b[i];
    } else {
        /* Assume starting estimates are passed as arguments */
        /* r    = b-A*x; */
        Atimesp_me_poly(sphere, dpoly, A, param, x, Ap);
        for (i = 0; i < m; i++)
            r[i] = b[i] - Ap[i];
    }

output_values_any_format("x.txt", sphere->n_points, x, 1);
output_values_any_format("Ap.txt", sphere->n_points, Ap, 1);
output_values_any_format("r.txt", sphere->n_points, r, 1);

    rtr     = dotprod(m, r, r);

    for (i = 0; i < m; i++)
        p[i] = 0.0;

    /* beta = 0; */
    beta    = 0.0;

    /* for it=1:nit, */
    for (it = 0; it < nit; it++) {
        /* if norm(r) < tol*norm(b), break; end; */
        if (norm(m, r) < nb)
            break;

        /* p      = r + beta*p; */
        for (i = 0; i < m; i++)
            p[i]  = r[i] + beta*p[i];

        /* Ap     = A*p; */
        Atimesp_me_poly(sphere, dpoly, A, param, p, Ap);

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


void
composition_poly(polygons_struct *sphere, double *A, double *B, double *C)
{
    double *Ax, *Ay, *Bx, *By, *Cx, *Cy;
    int i, mm = sphere->n_points;

printf("composition_poly\n");
    Ax = A;
    Ay = &A[mm];
    Bx = B;
    By = &B[mm];
    Cx = C;
    Cy = &C[mm];

    for(i=0; i<mm; i++) {
        double x, y;

        x     = Ax[i];
        y     = Ay[i];

        Cx[i]   = interp_uv_unit_sphere(sphere, Bx, x, y);
        Cy[i]   = interp_uv_unit_sphere(sphere, By, x, y);
    }
}


void
composition_jacobian_poly(polygons_struct *sphere,
                          double *A, double *JA, double *B, double *JB,
                          double *C, double *JC)
{
    double *Ax, *Ay, *JA00, *JA01, *JA10, *JA11;
    double *Bx, *By, *JB00, *JB01, *JB10, *JB11, jb[2][2];
    double *Cx, *Cy, *JC00, *JC01, *JC10, *JC11;
    int i, mm = sphere->n_points;

printf("composition_jacobian_poly: no Bx-correction!!!!!!!!!!\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n");
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

    for(i=0; i<mm; i++) {
        double x, y;

        x    = Ax[i];
        y    = Ay[i];

        Cx[i]    = interp_uv_unit_sphere(sphere, Bx, x, y);
        Cy[i]    = interp_uv_unit_sphere(sphere, By, x, y);
        jb[0][0] = interp_uv_unit_sphere(sphere, JB00, x, y);
        jb[0][1] = interp_uv_unit_sphere(sphere, JB01, x, y);
        jb[1][0] = interp_uv_unit_sphere(sphere, JB10, x, y);
        jb[1][1] = interp_uv_unit_sphere(sphere, JB11, x, y);

        JC00[i] = jb[0][0]*JA00[i] + jb[1][0]*JA01[i];
        JC01[i] = jb[0][1]*JA00[i] + jb[1][1]*JA01[i];
        JC10[i] = jb[0][0]*JA10[i] + jb[1][0]*JA11[i];
        JC11[i] = jb[0][1]*JA10[i] + jb[1][1]*JA11[i];
    }
}


void
composition_detjac_poly(polygons_struct *sphere, double *A, double *dA,
                        double *B, double *dB, double *C, double *dC)
{
    double *Ax, *Ay;
    double *Bx, *By, jb;
    double *Cx, *Cy;
    int i, mm = sphere->n_points;

printf("composition_detjac_poly\n");
    Ax   =   A;
    Ay   =  &A[mm];

    Bx   =   B;
    By   =  &B[mm];

    Cx   =   C;
    Cy   =  &C[mm];

    for (i = 0; i < mm; i++) {
        double x, y;

        x    = Ax[i];
        y    = Ay[i];

        Cx[i] = interp_uv_unit_sphere(sphere, Bx, x, y);
        Cy[i] = interp_uv_unit_sphere(sphere, By, x, y);
        jb    = interp_uv_unit_sphere(sphere, dB, x, y);
        dC[i] = jb*dA[i];
    }
}

double
samp_poly(polygons_struct *sphere, double f[], double x, double y)
{
    return(interp_uv_unit_sphere(sphere, f, x, y));
}


void
expdef_poly(polygons_struct *sphere, struct dartel_poly *dpoly, int k,
            double flow[], double t0[], double t1[], double J0[], double J1[])
{
    double *optr;
    double td;
    int m = dpoly->n_points;
    int i, j;

    optr = t0;

    td = 1;
    for (i = 0; i < k; i++)
        td = td*2;
    td = 1.0/td;

    if (J0 != (double *)0) {
printf("expdef_poly: J= != 0\n");
        for (i = 0; i < m; i++) {
            Point pdx, ndx, pdy, ndy;
            double u, v, px, nx, py, ny;

            t0[i  ] = dpoly->u[i] + flow[i  ]*td;
            t0[i+m] = dpoly->v[i] + flow[i+m]*td;

            px   = interp_point_unit_sphere(sphere, flow, dpoly->ptheta[i]);
            nx   = interp_point_unit_sphere(sphere, flow, dpoly->ntheta[i]);
            py   = interp_point_unit_sphere(sphere, flow, dpoly->pphi[i]);
            ny   = interp_point_unit_sphere(sphere, flow, dpoly->nphi[i]);
            J0[i    ] = (px - nx)/(2*UTHETA)*td/2 + 1.0;
            J0[i+2*m] = (py - ny)/(2*VPHI)*td/2;

            px   = interp_point_unit_sphere(sphere, flow + m, dpoly->ptheta[i]);
            nx   = interp_point_unit_sphere(sphere, flow + m, dpoly->ntheta[i]);
            py   = interp_point_unit_sphere(sphere, flow + m, dpoly->pphi[i]);
            ny   = interp_point_unit_sphere(sphere, flow + m, dpoly->nphi[i]);
            J0[i+  m] = (px - nx)/(2*UTHETA)*td/2;
            J0[i+3*m] = (py - ny)/(2*VPHI)*td/2 + 1.0;
        }
        for (i = 0; i < k; i++) {
            double *tmpp;
            composition_jacobian_poly(sphere, t0, J0, t0, J0, t1, J1);
            tmpp = t0; t0   = t1; t1   = tmpp;
            tmpp = J0; J0   = J1; J1   = tmpp;
        }
    } else {
printf("expdef_poly\n");
        for (i = 0; i < m; i++) {
            t0[i  ] = dpoly->u[i] + flow[i  ]*td;
            t0[i+m] = dpoly->v[i] + flow[i+m]*td;
        }
        for (i = 0; i < k; i++) {
            double *tmpp;
            composition_poly(sphere, t0, t0, t1);
            tmpp = t0; t0   = t1; t1   = tmpp;
        }
    }

    if (optr != t0) {
        for (i = 0; i < 2*m; i++)
            t1[i] = t0[i];

        if (J0 != (double *)0) {
            for (i = 0; i < 4*m; i++)
                J1[i] = J0[i];
        }
    }
}

void
expdefdet_poly(polygons_struct *sphere, struct dartel_poly *dpoly, int k,
               double flow[], double t0[], double t1[], double J0[],
               double J1[])
{
    double *optr;
    double td;
    int m = dpoly->n_points;
    int i, j;

printf("expdefdet_poly\n");
    optr = t0;

    td = 1;
    for(i=0; i<k; i++)
        td = td*2;
    td = 1.0/td;

    for(i = 0; i < m; i++) {
        double j00,j01,j10,j11;
        Point pdx, ndx, pdy, ndy;
        double u, v, px, nx, py, ny;

        t0[i  ] = dpoly->u[i] + flow[i  ]*td;
        t0[i+m] = dpoly->v[i] + flow[i+m]*td;

        px   = interp_point_unit_sphere(sphere, flow, dpoly->ptheta[i]);
        nx   = interp_point_unit_sphere(sphere, flow, dpoly->ntheta[i]);
        py   = interp_point_unit_sphere(sphere, flow, dpoly->pphi[i]);
        ny   = interp_point_unit_sphere(sphere, flow, dpoly->nphi[i]);
        j00 = (px - nx)/(2*UTHETA)*td/2 + 1.0;
        j01 = (py - ny)/(2*VPHI)*td/2;

        px   = interp_point_unit_sphere(sphere, flow + m, dpoly->ptheta[i]);
        nx   = interp_point_unit_sphere(sphere, flow + m, dpoly->ntheta[i]);
        py   = interp_point_unit_sphere(sphere, flow + m, dpoly->pphi[i]);
        ny   = interp_point_unit_sphere(sphere, flow + m, dpoly->nphi[i]);
        j10 = (px - nx)/(2*UTHETA)*td/2;
        j11 = (py - ny)/(2*VPHI)*td/2 + 1.0;

        J0[i]   = j00*j11 - j10*j01;
    }
    
    for(i=0; i<k; i++) {
        double *tmpp;
        composition_detjac_poly(sphere, t0, J0, t0, J0, t1, J1);
        tmpp = t0; t0   = t1; t1   = tmpp;
        tmpp = J0; J0   = J1; J1   = tmpp;
    }

    if (optr != t0) {
        for(i=0; i<2*m; i++)
            t1[i] = t0[i];

        for(i=0; i<4*m; i++)
            J1[i] = J0[i];
    }
}

/* J0 := J0*inv(I+diag(v0)*sc) */
void
jac_div_smalldef_poly(polygons_struct *sphere, struct dartel_poly *dpoly,
                      double sc, double v0[], double J0[])
{
    int i;
    int m = sphere->n_points;
    double sc2 = sc/2.0;
    double *v1 = v0+m;

printf("jac_div_smalldef_poly\n");
    for(i = 0; i < m; i++) {
        double j00,j01, j10,j11;
        double t00,t01, t10,t11;
        Point om1, op1;
        double u, v, idt, km1, kp1;

        km1 = interp_point_unit_sphere(sphere, v0, dpoly->ntheta[i]);
        kp1 = interp_point_unit_sphere(sphere, v0, dpoly->ptheta[i]);
        j00 = (kp1 - km1)*sc2 + 1.0;
        km1 = interp_point_unit_sphere(sphere, v1, dpoly->ntheta[i]);
        kp1 = interp_point_unit_sphere(sphere, v1, dpoly->ptheta[i]);
        j01 = (kp1 - km1)*sc2;

        km1 = interp_point_unit_sphere(sphere, v0, dpoly->nphi[i]);
        kp1 = interp_point_unit_sphere(sphere, v0, dpoly->pphi[i]);
        j10 = (kp1 - km1)*sc2;
        km1 = interp_point_unit_sphere(sphere, v1, dpoly->nphi[i]);
        kp1 = interp_point_unit_sphere(sphere, v1, dpoly->pphi[i]);
        j11 = (kp1 - km1)*sc2 + 1.0;

        idt = 1.0/(j00*j11 - j10*j01);
        t00 =  idt*j11;
        t01 = -idt*j10;
        t10 = -idt*j01;
        t11 =  idt*j00;

        j00 = J0[i  ]; j01 = J0[i+m*2];
        j10 = J0[i+m]; j11 = J0[i+m*3];
        J0[i    ] = j00*t00 + j01*t10;
        J0[i+m  ] = j10*t00 + j11*t10;
        J0[i+m*2] = j00*t01 + j01*t11;
        J0[i+m*3] = j10*t01 + j11*t11;
    }
}

double
init_objfun_poly(polygons_struct *sphere, struct dartel_poly *dpoly,
                 double f[], double g[], double t0[], double J0[],
                 double dj[], double b[], double A[])
{
    int j, m = dpoly->n_points;
    double ssl = 0.0, dt = 1.0;

printf("init_objfun_poly\n");
    for (j = 0; j < m; j++) {
        double x, y;
        int    ix, iy;
        double k11, k12, k21, k22;
        Point o11, o12, o21, o22;
        double dx0, dx1, dx2, dy0, dy1, dy2;
        double d, dx, dy;

        x    = t0[j  ];
        y    = t0[j+m];

        k22   = interp_uv_unit_sphere(sphere, f, x, y);
        k12   = interp_uv_unit_sphere(sphere, f, x+UTHETA, y);
        k21   = interp_uv_unit_sphere(sphere, f, x, y+VPHI);

        d    = k22 - g[j];
        dx0  = (k12 - k22)/(UTHETA);
        dy0  = (k21 - k22)/(VPHI);

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

/*
 * t0 = Id + v0*sc
 * J0 = Id + I+diag(v0)*sc
 */
void
smalldef_jac_poly(polygons_struct *sphere, struct dartel_poly *dpoly, double sc,
                  double v0[], double t0[], double J0[])
{
    int i, m = dpoly->n_points;
    double sc2 = sc/2.0;
    double *v1 = v0+m;

printf("smalldef_jac_poly\n");
    for(i = 0; i < m; i++) {
        Point om1, op1;
        double u, v, km1, kp1;

        //point_to_uv(&sphere->points[i], &u, &v);

        t0[i  ] = dpoly->u[i] + v0[i]*sc;
        t0[i+m] = dpoly->v[i] + v1[i]*sc;

        km1 = interp_point_unit_sphere(sphere, v0, dpoly->ntheta[i]);
        kp1 = interp_point_unit_sphere(sphere, v0, dpoly->ptheta[i]);
        J0[i    ] = (kp1 - km1)*sc2 + 1.0;
        km1 = interp_point_unit_sphere(sphere, v1, dpoly->ntheta[i]);
        kp1 = interp_point_unit_sphere(sphere, v1, dpoly->ptheta[i]);
        J0[i+  m] = (kp1 - km1)*sc2;

        km1 = interp_point_unit_sphere(sphere, v0, dpoly->nphi[i]);
        kp1 = interp_point_unit_sphere(sphere, v0, dpoly->pphi[i]);
        J0[i+2*m] = (kp1 - km1)*sc2;
        km1 = interp_point_unit_sphere(sphere, v1, dpoly->nphi[i]);
        kp1 = interp_point_unit_sphere(sphere, v1, dpoly->pphi[i]);
        J0[i+3*m] = (kp1 - km1)*sc2 + 1.0;
    }
}


void
squaring_poly(polygons_struct *sphere, struct dartel_poly *dpoly, int k,
              int save_transf, double b[], double A[], double t0[], double t1[],
              double J0[], double J1[])
{
    int i, j, m = dpoly->n_points;
    double *ptr = t0;

printf("squaring_poly\n");
    for(i=0; i<k; i++) {
        double *buf1, *buf2;
        buf1 = t1; /* Re-use some memory */
        buf2 = J1;

        for(j=0; j<m; j++) {
            double tmp00, tmp01, tmp11, tmp1, tmp2, tmp3, tmp4;
            double x, y;
            double j11, j21, j12, j22, dt;

            x   = t0[j  ];
            y   = t0[j+m];
            j11 = J0[j  ]; j12 = J0[j+2*m];
            j21 = J0[j+m]; j22 = J0[j+3*m];
            dt  = j11*j22 - j12*j21;

         /* if (dt < 1e-9) dt = 1e-9;
            if (dt > 1e9 ) dt = 1e9; */

            tmp1      = samp_poly(sphere, b  , x, y);
            tmp2      = samp_poly(sphere, b+m, x, y);

            buf1[j  ] = dt*(tmp1*j11+tmp2*j21);
            buf1[j+m] = dt*(tmp1*j12+tmp2*j22);

            tmp00 = samp_poly(sphere, A      , x, y);
            tmp11 = samp_poly(sphere, A +   m, x, y);
            tmp01 = samp_poly(sphere, A + 2*m, x, y);

            tmp1  = tmp00*j11 + tmp01*j21;
            tmp2  = tmp01*j11 + tmp11*j21;
            tmp3  = tmp00*j12 + tmp01*j22;
            tmp4  = tmp01*j12 + tmp11*j22;

         /* if (dt < 1e-9) dt = 1e-9; */
         /* if (dt > 1e9 ) dt = 1e9; */

            buf2[j    ] = dt*(tmp1*j11 + tmp2*j21);
            buf2[j+  m] = dt*(tmp3*j12 + tmp4*j22);
            buf2[j+2*m] = dt*(tmp1*j12 + tmp2*j22);
        }
        for(j=0; j<2*m; j++) b[j] += buf1[j];
        for(j=0; j<3*m; j++) A[j] += buf2[j];

        if (save_transf || (i<k-1)) {
            double *tmpp;
            composition_jacobian_poly(sphere, t0, J0, t0, J0, t1, J1);
            tmpp = t0; t0   = t1; t1   = tmpp;
            tmpp = J0; J0   = J1; J1   = tmpp;
        }
    }
    if (save_transf && ptr!=t0) {
        for(j=0; j<m*2; j++) t1[j] = t0[j];
        for(j=0; j<m*4; j++) J1[j] = J0[j];
    }
}


void
init_dartel_poly(polygons_struct *sphere, struct dartel_poly *dpoly)
{
        int i;
        double x, y, z, xo, yo, zo;

        dpoly->n_points = sphere->n_points;
        dpoly->u = (double *) malloc(sizeof(double) * dpoly->n_points);
        dpoly->v = (double *) malloc(sizeof(double) * dpoly->n_points);
        dpoly->ntheta = (Point *) malloc(sizeof(Point) * dpoly->n_points);
        dpoly->ptheta = (Point *) malloc(sizeof(Point) * dpoly->n_points);
        dpoly->nphi = (Point *) malloc(sizeof(Point) * dpoly->n_points);
        dpoly->pphi = (Point *) malloc(sizeof(Point) * dpoly->n_points);

        for (i = 0; i < sphere->n_points; i++) {
                xo = Point_x(sphere->points[i]);
                yo = Point_y(sphere->points[i]);
                zo = Point_z(sphere->points[i]);

                point_to_uv(&sphere->points[i], &dpoly->u[i], &dpoly->v[i]);

               /* +THETA */
                x = xo * cos(THETA)
                  + zo * sin(THETA);
                y = xo * -sin(THETA) * sin(THETA)
                  + yo * cos(THETA)
                  + zo * cos(THETA) * sin(THETA);
                z = xo * -sin(THETA) * cos(THETA)
                  + yo * -sin(THETA)
                  + zo * cos(THETA) * cos(THETA);
                fill_Point(dpoly->ptheta[i], x, y, z);

               /* -THETA */
                x = xo * cos(-THETA)
                  + zo * sin(-THETA);
                y = xo * -sin(-THETA) * sin(-THETA)
                  + yo * cos(-THETA)
                  + zo * cos(-THETA) * sin(-THETA);
                z = xo * -sin(-THETA) * cos(-THETA)
                  + yo * -sin(-THETA)
                  + zo * cos(-THETA) * cos(-THETA);
                fill_Point(dpoly->ntheta[i], x, y, z);

                /* +PHI */
                x = xo * cos(PHI)
                  + yo * sin(PHI);
                y = xo * -sin(PHI)
                  + yo * cos(PHI);
                z = zo;
                fill_Point(dpoly->pphi[i], x, y, z);

                /* -PHI */
                x = xo * cos(-PHI)
                  + yo * sin(-PHI);
                y = xo * -sin(-PHI)
                  + yo * cos(-PHI);
                fill_Point(dpoly->nphi[i], x, y, z);
        }
}

int
cgs2_scratchsize(polygons_struct *sphere)
{
    int    n[32], m[32], bs, j;
    bs = 0;
    n[0] = sphere->n_points;

    for(j=1; j<16; j++) {
        n[j] = ceil(n[j-1]/2.0);
        m[j] = n[j]*n[j];
        bs += m[j];
        if (n[j] < 2)
            break;
    }
    return((2*n[0]*n[0] + n[0]*n[1] + 9*bs));
}


int
dartel_scratchsize_poly(polygons_struct *sphere, int code)
{
        int m1, m2;
        int m = sphere->n_points;

        m1 = 15*m;
        if (code) m1 += 5*m;
        m2 = 5*m + cgs2_scratchsize(sphere);

        if (m1 > m2)
                return(m1);
        else
                return(m2);
}

void
dartel_poly2(polygons_struct *sphere, struct dartel_poly *dpoly,
            struct dartel_prm prm, double v[], double g[], double f[],
            double dj[], double ov[], double ll[], double *buf)
{
    double *sbuf;
    double *b, *A, *b1, *A1;
    double *t0, *t1, *J0, *J1;
    double sc;
    double ssl, ssp;
    double normb;
    int j, m = sphere->n_points;

    /*
     *  Allocate memory.
     *    0  1  2  3  4  5  6  7  8  9 10 11 12 13 14
     *  [ A  A  A  t  t  J  J  J  J  t  t  J  J  J  J] for computing derivatives
     *  [ A  A  A s1 s2 s3 s4 s5 s6 s7 s8] for CGS solver
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

output_values_any_format("J00.txt", sphere->n_points, J0, 1);
output_values_any_format("t00.txt", sphere->n_points, t0, 1);
output_values_any_format("v.txt", sphere->n_points, v, 1);
    expdef_poly(sphere, dpoly, prm.k, v, t0, t1, J0, J1);

    printf("expdef: v[0] = %f, v[m/2] = %f, v[m] = %f\n", v[0], v[m/2], v[m]);
    printf("expdef: t0[0] = %f, t0[m/2] = %f, t0[m] = %f\n", t0[0], t0[m/2], t0[m]);
    printf("expdef: t1[0] = %f, t1[m/2] = %f,t1[m] = %f\n", t1[0], t1[m/2], t1[m]);
    printf("expdef: J0[0] = %f, J0[m/2] = %f,J0[m] = %f\n", J0[0], J0[m/2], J0[m]);
    printf("expdef: J1[0] = %f, J1[m/2] = %f,J1[m] = %f\n", J1[0], J1[m/2], J1[m]);

output_values_any_format("J01.txt", sphere->n_points, J0, 1);
output_values_any_format("t01x.txt", sphere->n_points, t0, 1);
output_values_any_format("t01y.txt", sphere->n_points, t0+m, 1);
    jac_div_smalldef_poly(sphere, dpoly, sc, v, J0);
output_values_any_format("J02.txt", sphere->n_points, J0, 1);
    
    ssl = init_objfun_poly(sphere, dpoly, f, g, t0, J0, dj, b, A);
output_values_any_format("A0.txt", sphere->n_points, A, 1);
output_values_any_format("b0x.txt", sphere->n_points, b, 1);
output_values_any_format("b0y.txt", sphere->n_points, b+m, 1);
    
    smalldef_jac_poly(sphere, dpoly, -sc, v, t0, J0);
    
    squaring_poly(sphere, dpoly, prm.k, prm.code==1, b, A, t0, t1, J0, J1);
output_values_any_format("A1.txt", sphere->n_points, A, 1);
output_values_any_format("b1x.txt", sphere->n_points, b, 1);
output_values_any_format("b1y.txt", sphere->n_points, b+m, 1);
output_values_any_format("t02.txt", sphere->n_points, t0, 1);
output_values_any_format("t10.txt", sphere->n_points, t1, 1);
    
    jac_div_smalldef_poly(sphere, dpoly, -sc, v, J0);
    
    ssl += init_objfun_poly(sphere, dpoly, g, f, t0, J0, (double *)0, b1, A1);
    
    smalldef_jac_poly(sphere, dpoly, sc, v, t0, J0);
    
    squaring_poly(sphere, dpoly, prm.k, 0, b1, A1, t0, t1, J0, J1);

    for(j = 0; j < m*2; j++) b[j] -= b1[j];
    for(j = 0; j < m*3; j++) A[j] += A1[j];

    LtLf_me_poly(sphere, dpoly, v, prm.rparam, t1);

    ssp = 0.0;
    for(j = 0; j < 2*m; j++) {
        b[j] = b[j]*sc + t1[j];
        ssp += t1[j]*v[j];
    }
    normb = norm(2*m, b);

    for(j = 0; j < 3*m; j++) A[j] *= sc;
    for(j = 0; j < 2*m; j++) A[j] += prm.lmreg;

output_values_any_format("b2.txt", sphere->n_points, b, 1);
output_values_any_format("A2.txt", sphere->n_points, A, 1);
output_values_any_format("v0.txt", sphere->n_points, v, 1);
    /* Solve equations for Levenberg-Marquardt update:
     * v = v - inv(H + L'*L + R)*(d + L'*L*v)
     *     v: velocity or flow field
     *     H: matrix of second derivatives
     *     L: regularisation (L'*L is the inverse of the prior covariance)
     *     R: Levenberg-Marquardt regularisation
     *     d: vector of first derivatives
     */
output_values_any_format("J03.txt", sphere->n_points, J0, 1);
output_values_any_format("J11.txt", sphere->n_points, J1, 1);
output_values_any_format("sbuf0.txt", sphere->n_points, sbuf, 1);
    cgs2_poly(sphere, dpoly, A, b, prm.rparam, 1e-6, 100, sbuf, sbuf+2*m,
              sbuf+4*m, sbuf+6*m);
    //fmg2(sphere, A, b, prm.rtype, prm.rparam, prm.cycles, prm.its, sbuf, sbuf+2*m);

output_values_any_format("sbuf.txt", sphere->n_points, sbuf, 1);
output_values_any_format("sbuf2.txt", sphere->n_points, sbuf+2*m, 1);
output_values_any_format("sbuf4.txt", sphere->n_points, sbuf+4*m, 1);
output_values_any_format("sbuf6.txt", sphere->n_points, sbuf+6*m, 1);
    for(j = 0; j < 2*m; j++) ov[j] = v[j] - sbuf[j];
     ll[0] = ssl;
    ll[1] = ssp*0.5;
    ll[2] = normb;
}

