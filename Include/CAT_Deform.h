/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry, University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 *
 */

#ifndef _CAT_DEFORM_H_
#define _CAT_DEFORM_H_

void smooth_displacement_field(double (*displacement_field)[3], polygons_struct *polygons, 
                      int *n_neighbours, int **neighbours, int iterations, 
                      double sigma);
void surf_deform(polygons_struct *polygons, float *input, nifti_image *nii_ptr, 
                      double w[3], double sigma, float lim, int it, int remove_intersections, 
                      int verbose);
void surf_deform_dual(polygons_struct *polygons1, polygons_struct *polygons2, 
                      polygons_struct *polygons_orig, float *input, nifti_image *nii_ptr, 
                      double w[3], double sigma, float lim1, float lim2, double *target_distance,
                      int it, int verbose);

#endif
