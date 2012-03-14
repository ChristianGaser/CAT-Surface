/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 *
 */

void get_all_polygon_point_neighbours(polygons_struct *, int *[], int **[]);
void heatkernel_blur_points(int, Point [], double [], int, int *, int, Real,
                            Point *, double *);
void smooth_heatkernel(polygons_struct *, double *, double);

