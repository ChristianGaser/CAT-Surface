/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 *
 */

Vector projectToPlane(Vector projected, Vector basis[2]);
Vector projection(Vector vector, Vector normal);

void leastSquares(const int num, Vector dc[], Vector dn[], Real *k1, Real *k2);
void compute_points_centroid_and_normal_cg(polygons_struct *, int, int, int [],
                                           Point *, Vector *, Real *, int,
                                           Real *);
void get_polygon_vertex_curvatures_cg(polygons_struct *, int [], int *[],
                                      Real, int, Real []);

void get_smoothed_curvatures(polygons_struct *, double *,
                             double, int);
void compute_local_sharpness(polygons_struct *, int [], int * [], double *);
void calc_convexity(polygons_struct *, int [], int *[], double, double *);
