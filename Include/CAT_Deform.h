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

/**
 * \brief Refine pial and white surfaces toward image intensity edges using
 *        normal-ray edge search (FreeSurfer-inspired approach).
 *
 * After surf_deform_dual has positioned the surfaces near their targets,
 * this function performs additional refinement by searching along each
 * vertex normal for the steepest intensity gradient (the actual tissue
 * boundary), then nudging vertices toward that edge location.
 *
 * \param polygons1        (in/out) pial surface mesh
 * \param polygons2        (in/out) white surface mesh
 * \param input            (in)     intensity volume data
 * \param nii_ptr          (in)     NIfTI header for dimensions and transforms
 * \param lim1             (in)     intensity threshold for pial surface
 * \param lim2             (in)     intensity threshold for white surface
 * \param target_distance  (in)     desired pial-white spacing per vertex
 * \param it               (in)     number of refinement iterations
 * \param verbose          (in)     if nonzero, print progress
 */
void surf_deform_gradient_dual(polygons_struct *polygons1, polygons_struct *polygons2,
                               float *input, nifti_image *nii_ptr,
                               float lim1, float lim2,
                               double *target_distance, int it, int verbose);

#endif
