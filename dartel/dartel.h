/* (c) John Ashburner (2007) */

#include <bicpl.h>

#define LOG(x) (((x)>0) ? log(x+0.001): -6.9078)
#define WRAP(i,m) (((i)>=0) ? (i)%(m) : ((m)+(i)%(m))%m)

#define S 1.0000000001

#define THETA 0.0175 /* theta step.. 1 degree */
#define PHI 0.0175 /* phi step 1 degree */
#define UTHETA 0.0028
#define VPHI 0.0028

struct dartel_prm {
  int rtype;         /* regularization type: 0 - linear elastic energy; */
                     /* 1 - membrane energy; 2 - bending energy */
  double rparam[5];  /* regularization parameters: width height mu lambda id */
  double lmreg;      /* LM regularization */
  int cycles;        /* # of cycles for full multi grid (FMG) */
  int its;           /* # of relaxation iterations in each multigrid cycle */
  int k;             /* time steps for solving the PDE */
  int code;          /* objective function: 0 - sum of squares; */
                     /* 1 - symmetric sum of squares */
};

struct dartel_poly {
        int n_points;
        int *n_neighbours;
        int **neighbours;

        double *u;
        double *v;
        Point *ntheta;
        Point *ptheta;
        Point *nphi;
        Point *pphi;

        double **du;
        double **dv;
        Vector *uv;
};


/* functions in diffeosphere.c */
void composition(int dm[], double *A, double *B, double *C);
void composition_jacobian(int dm[], double *A, double *JA, double *B,
                          double *JB, double *C, double *JC);
void composition_detjac(int dm[], double *A, double *dA, double *B, double *dB,
                        double *C, double *dC);

double samp(int dm[], double f[], double x, double y);

void expdef(int dm[], int k, double v[], double t0[], double t1[],
            double J0[], double J1[]);
void expdefdet(int dm[], int k, double v[], double t0[], double t1[],
            double J0[], double J1[]);

int pow2(int k);

void jac_div_smalldef(int dm[], double sc, double v0[], double J0[]);

double initialise_objfun2(int dm[], double f[], double g[], double t0[],
                          double J0[], double dj[], double b[], double A[]);
double initialise_objfun(int dm[], double f[], double g[], double t0[],
                         double J0[], double dj[], double b[], double A[]);
double initialise_objfun_mn(int dm[], double f[], double g[], double t0[],
                            double J0[], double jd[], double b[], double A[]);

void smalldef_jac(int dm[], double sc, double v0[], double t0[], double J0[]);

void squaring(int dm[], int k, int save_transf, double b[], double A[],
              double t0[], double t1[], double J0[], double J1[]);

void unwrap(int dm[], double f[]);

int dartel_scratchsize(int dm[], int code);
void dartel(struct dartel_prm prm, int dm[], double v[], double g[],
            double f[], double dj[], double ov[], double ll[], double *buf);
void init_dartel_poly(polygons_struct *, struct dartel_poly *);


/* functions in optimizersphere.c */
static int neumann(int i, int m);

double sumsq_le(int dm[], double a[], double b[], double s[], double u[]);
void LtLf_le(int dm[], double f[], double s[], double g[]);
void relax_le(int [], double a[], double b[], double s[], int nit, double u[]);
void Atimesp_le(int dm[], double A[], double s[], double p[], double Ap[]);

double sumsq_me(int dm[], double a[], double b[], double s[], double u[]);
void LtLf_me(int dm[], double f[], double s[], double g[]);
void relax_me(int [], double a[], double b[], double s[], int nit, double u[]);
void Atimesp_me(int dm[], double A[], double s[], double p[], double Ap[]);

double sumsq_be(int dm[], double a[], double b[], double s[], double u[]);
void LtLf_be(int dm[], double f[], double s[], double g[]);
void relax_be(int [], double a[], double b[], double s[], int nit, double u[]);
void Atimesp_be(int [], double A[], double param[], double p[], double Ap[]);

void solve22(double a[], double b[], double t,  double u[]);
double dotprod(int m, double a[], double b[]);
void addscaled(int m, double a[], double b[], double s);
double norm(int m, double a[]);

void cgs2(int dm[], double A[], double b[], int rtype, double param[],
          double tol, int nit, double x[], double r[], double p[], double Ap[]);

double wt2(double x);

void resize(int na[], double *a, int nc[], double *c, double *b);
void rescale(int n, double *a, double s);
void restrict2(int n, int na[], double *a, int nc[], double *c, double *b);
void prolong(int n, int na[], double *a, int nc[], double *c, double *b);

static void zeros(int n, double *a);
static void copy(int n, double *a, double *b);
static void addto(int n, double *a, double *b);

int fmg2_scratchsize(int n0[]);
void fmg2(int n0[], double *a0, double *b0, int rtype, double param0[], int c,
          int nit, double *u0, double *scratch);
