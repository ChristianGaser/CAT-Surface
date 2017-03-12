/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 *
*/

#include <bicpl.h>
#include <fftw3.h>

#include "s2kit10/makeweights.h"
#include "s2kit10/FST_semi_fly.h"

#define BINTREE_FACTOR   0.5

int read_SPHxyz(char *, int, double **, double **, double **, double **,
                double **, double **);
int write_SPHxyz(char *, int, double *, double *, double *, double *,
                 double *, double *);

void sample_sphere_from_sph(double *, double *, double *, polygons_struct *,
                            int, polygons_struct *, int);

void replaceSPH(int, int, double *, double *);

void shape_description(int, double *, double *, double *, double *, double *, 
                       double *, double *);
void butterworth_filter(int, int, double *, double *);
void bandpass_bandwidth(int, int, int, double *, double *);
void limit_bandwidth(int, int, double *, double *);

void get_sph_coeffs_of_realdata(double *, int, int, double *, double *);
void get_realdata_from_sph_coeffs(double *, int, int, double *, double *);

void get_equally_sampled_coords_of_polygon(polygons_struct *,
                                           polygons_struct *, int, double [],
                                           double [], double []);
void get_equally_sampled_coords_holes(polygons_struct *, polygons_struct *,
                                      int *, int, int *, int, double [],
                                      double [], double []);

object_struct ** create_equally_sampled_unit_sphere(int, int);

double * laplace2d(double *, unsigned char *, int, int, double);
