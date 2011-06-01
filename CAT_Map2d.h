/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 *
*/

typedef struct {
	long x;
	long y;
} Header;

typedef struct {
	double x;
	double y;
} Vector2D;

void point_to_uv(Point *, double *, double *);
void uv_to_point(double, double, Point *);
void wrap_sheet2d(double *data, int *dm, int wrap, double *wdata);
void unwrap_sheet2d(double *wdata, int *wdm, int wrap, double *data);

void map_smoothed_curvature_to_sphere(polygons_struct *, polygons_struct *,
                                      double *, double *, double, int *, int);
void map_sheet2d_to_sphere(double *, double *, polygons_struct *, int, int *);
void smooth_sheet2d(double *, int[] , double);
