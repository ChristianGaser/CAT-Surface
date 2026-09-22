/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry, University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 *
 */

#ifndef _CAT_DEFORM_H_
#define _CAT_DEFORM_H_

#include <bicpl.h>
#include "nifti1_io.h"

/**
 * \brief Smooth a displacement field using neighborhood averaging with exponential decay.
 *
 * \param displacement_field (in/out) double[n_points][3]; per-vertex displacement vectors
 * \param polygons            (in)    surface mesh
 * \param n_neighbours        (in)    int[n_points]; vertex degree (number of neighbors)
 * \param neighbours          (in)    int*[n_points]; neighbor indices per vertex
 * \param iterations          (in)    number of smoothing passes
 * \param sigma               (in)    exponential decay parameter of attenuation
 */
void smooth_displacement_field(double (*displacement_field)[3],
                               polygons_struct *polygons, int *n_neighbours,
                               int **neighbours, int iterations, double sigma);
/**
 * \brief Deform a surface mesh toward intensity gradients from an external volume.
 *
 * \param polygons            (in/out) surface mesh (modified in-place)
 * \param input               (in)     float[nvoxels]; intensity volume data
 * \param nii_ptr             (in)     NIfTI header with volume dimensions and voxel size
 * \param w                   (in)     double[3]; weight factors {smoothness, gradient_strength, ?}
 * \param sigma               (in)     Gaussian smoothing parameter for displacement field
 * \param lim                 (in)     intensity threshold controlling deformation magnitude
 * \param it                  (in)     number of deformation iterations
 * \param remove_selfintersect (in)    boolean; enable self-intersection removal
 * \param verbose             (in)     boolean; print iteration progress if true
 */
void surf_deform(polygons_struct *polygons, float *input, nifti_image *nii_ptr,
                 double w[3], double sigma, float lim, int it,
                 int remove_selfintersect, int verbose);
/**
 * \brief Simultaneously deform two surfaces (e.g., white+pial) toward image gradients.
 *
 * \param polygons1         (in/out) first deformed surface (e.g., pial matter)
 * \param polygons2         (in/out) second deformed surface (e.g., white surface)
 * \param polygons_orig     (in)     reference surface (starting template)
 * \param input             (in)     float[nvoxels]; intensity volume data
 * \param nii_ptr           (in)     NIfTI header with volume dimensions and voxel size
 * \param w                 (in)     double[3]; weight factors for deformation forces
 * \param sigma             (in)     Gaussian smoothing parameter for displacement field
 * \param lim1              (in)     intensity threshold for first surface deformation
 * \param lim2              (in)     intensity threshold for second surface deformation
 * \param target_distance   (in/out) double[n_points]; desired separation between surfaces
 * \param it                (in)     number of deformation iterations
 * \param verbose           (in)     boolean; print iteration progress if true
 */
void surf_deform_dual(polygons_struct *polygons1, polygons_struct *polygons2,
                      polygons_struct *polygons_orig, float *input,
                      nifti_image *nii_ptr, double w[3], double sigma,
                      float lim1, float lim2, double *target_distance, int it,
                      int verbose);

#endif
