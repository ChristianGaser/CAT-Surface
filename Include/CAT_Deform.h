/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry, University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 *
 */

void smooth_displacement_field(double (*displacement_field)[3], polygons_struct *polygons, 
                               int *n_neighbours, int **neighbours, int iterations, double sigma);
void surf_deform(polygons_struct *polygons, float *input, nifti_image *nii_ptr, 
                 double w[3], double sigma, float lim, int it, int verbose);                               