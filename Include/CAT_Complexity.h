/* Rachel Yotter - rachel.yotter@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 *
 */

#ifndef _CAT_COMPLEXITY_H_
#define _CAT_COMPLEXITY_H_

#define DATAFORMAT 1 /* 1 = double data, 0 = complex data */
#define SPH_ITERS 10
#define BW 1024
#define FWHM 30.0

double slope(double *, double *, int);
double get_globalfd(double *, double *, int);
void get_localfd(polygons_struct *, double *, double **, int, double *, int);
int min_triangles_update(int *, int *, int *);
double fractal_dimension(polygons_struct *, polygons_struct *, int, char *,
                         int, int);
void get_smoothed_values(polygons_struct *, double *, double);
double fractal_dimension_sph(polygons_struct *, polygons_struct *, char *, int,
                             polygons_struct *, int, int);

#endif
