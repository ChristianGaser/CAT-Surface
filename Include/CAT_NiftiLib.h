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
 * \brief Public API for equal_image_dimensions.
 *
 * This function is part of the CAT-Surface public library interface and is used by command-line tools.
 *
 * \param nii_ptr (in/out) Parameter of equal_image_dimensions.
 * \param nii_ptr2 (in/out) Parameter of equal_image_dimensions.
 * \return Return value of equal_image_dimensions.
 */
int equal_image_dimensions(nifti_image *nii_ptr, nifti_image *nii_ptr2);
void init_nifti_header(nifti_image *nii_ptr);
/**
 * \brief Public API for write_nifti_double.
 *
 * This function is part of the CAT-Surface public library interface and is used by command-line tools.
 *
 * \param output_filename (in/out) Parameter of write_nifti_double.
 * \param image (in/out) Parameter of write_nifti_double.
 * \param data_type (in/out) Parameter of write_nifti_double.
 * \param slope (in/out) Parameter of write_nifti_double.
 * \param dim (in/out) Parameter of write_nifti_double.
 * \param vox (in/out) Parameter of write_nifti_double.
 * \param in_ptr (in/out) Parameter of write_nifti_double.
 * \return Return value of write_nifti_double.
 */
int write_nifti_double( const char *output_filename, double image[], int data_type, double slope, int dim[], double vox[], nifti_image *in_ptr);
/**
 * \brief Public API for write_nifti_float.
 *
 * This function is part of the CAT-Surface public library interface and is used by command-line tools.
 *
 * \param output_filename (in/out) Parameter of write_nifti_float.
 * \param image (in/out) Parameter of write_nifti_float.
 * \param data_type (in/out) Parameter of write_nifti_float.
 * \param slope (in/out) Parameter of write_nifti_float.
 * \param dim (in/out) Parameter of write_nifti_float.
 * \param vox (in/out) Parameter of write_nifti_float.
 * \param in_ptr (in/out) Parameter of write_nifti_float.
 * \return Return value of write_nifti_float.
 */
int write_nifti_float( const char *output_filename, float image[], int data_type, double slope, int dim[], double vox[], nifti_image *in_ptr);
nifti_image *read_nifti_double( const char *input_filename, double *image[], int read_data);
nifti_image *read_nifti_float( const char *input_filename, float *image[], int read_data);

#endif
