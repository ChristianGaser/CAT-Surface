/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 *
 */

#ifndef _CAT_NIFTILIB_H_
#define _CAT_NIFTILIB_H_

#include "nifti1_io.h"
#include "nifti1_local.h"

#include <float.h>

#ifndef MAX
#define MAX(A,B) ((A) > (B) ? (A) : (B))
#endif

#ifndef MIN
#define MIN(A,B) ((A) < (B) ? (A) : (B))
#endif

/**
 * \brief Compare image dimensions and voxel sizes of two NIfTI images.
 *
 * \param nii_ptr  (in) first NIfTI image
 * \param nii_ptr2 (in) second NIfTI image to compare
 * \return 1 if dimensions match (within tolerance), 0 if mismatched
 */
int equal_image_dimensions(nifti_image *nii_ptr, nifti_image *nii_ptr2);
/**
 * \brief Initialize NIfTI image header structure to reasonable defaults.
 *
 * \param nii_ptr (in/out) NIfTI image structure to initialize
 */
void init_nifti_header(nifti_image *nii_ptr);
/**
 * \brief Write double-precision volume data to NIfTI file with type conversion.
 *
 * \param output_filename (in) path to output .nii/.nii.gz file
 * \param image           (in) linear array of voxel values (size dim[0]*dim[1]*dim[2])
 * \param data_type       (in) target NIfTI datatype (DT_UINT8, DT_INT16, DT_FLOAT32, etc.)
 * \param slope           (in) optional scaling slope; 0.0 for auto-scaling
 * \param dim            (in) array of 3 image dimensions [nx, ny, nz]
 * \param vox            (in) array of 3 voxel spacings [dx, dy, dz]
 * \param in_ptr         (in) optional template NIfTI header; NULL to use defaults
 * \return 1 on success; 0 on error (invalid datatype, no extension, write failed)
 */
int write_nifti_double(const char *output_filename, double image[],
                       int data_type, double slope, int dim[], double vox[],
                       nifti_image *in_ptr);
/**
 * \brief Write single-precision volume data to NIfTI file with type conversion.
 *
 * \param output_filename (in) path to output .nii/.nii.gz file
 * \param image           (in) linear array of float voxel values
 * \param data_type       (in) target NIfTI datatype
 * \param slope           (in) optional scaling slope; 0.0 for auto-scaling
 * \param dim            (in) array of 3 image dimensions [nx, ny, nz]
 * \param vox            (in) array of 3 voxel spacings [dx, dy, dz]
 * \param in_ptr         (in) optional template NIfTI header; NULL to use defaults
 * \return 1 on success; 0 on error (invalid datatype, no extension, write failed)
 */
int write_nifti_float(const char *output_filename, float image[], int data_type,
                      double slope, int dim[], double vox[],
                      nifti_image *in_ptr);
/**
 * \brief Read NIfTI image file and load data into double-precision array.
 *
 * \param input_filename (in)  path to NIfTI file
 * \param image          (out) pointer to allocated double array (size nx*ny*nz)
 * \param read_data      (in)  if non-zero, read voxel data; if 0, read header only
 * \return pointer to NIfTI image structure on success; NULL on error
 */
nifti_image *read_nifti_double(const char *input_filename, double *image[],
                               int read_data);
/**
 * \brief Read NIfTI image file and load data into single-precision array.
 *
 * \param input_filename (in)  path to NIfTI file
 * \param image          (out) pointer to allocated float array (size nx*ny*nz)
 * \param read_data      (in)  if non-zero, read voxel data; if 0, read header only
 * \return pointer to NIfTI image structure on success; NULL on error
 */
nifti_image *read_nifti_float(const char *input_filename, float *image[],
                              int read_data);

#endif
