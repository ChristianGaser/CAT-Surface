/* Rachel Yotter - rachel.yotter@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id: CAT_Complexity.h 139 2009-10-02 15:39:28Z raytrace $
 *
 */

#define DATAFORMAT 1 /* 1 = real data, 0 = complex data */

#define SPH_ITERS 10
#define BW 1536
#define FWHM 30.0

double slope(double *, double *, int);
double get_globalfd(double *, double *, int);
void get_localfd(polygons_struct *, double *, double **, int, double *, int);

object_struct ** create_resampling_sphere(double, int);
int min_triangles_update(int *, int *, int *, int *);
object_struct ** resample_surface(polygons_struct *, polygons_struct *,
                                  int, double *, double *);
double fractal_dimension(polygons_struct *, polygons_struct *, int, char *,
                         int, int);

void get_smoothed_values(polygons_struct *, double *, double);
double fractal_dimension_sph(polygons_struct *, polygons_struct *, char *, int,
                             polygons_struct *, int, int);
