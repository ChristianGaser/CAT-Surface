/*
 * Christian Gaser
 * $Id: niilib.h 191 2012-07-31 15:41:36Z gaser $ 
 *
 */

#ifndef _NIILIB_H_
#define _NIILIB_H_

#include <stdlib.h>

#include "nifti/nifti1_io.h"
#include "nifti/nifti1_local.h"

#include <float.h>

#ifndef MAX
#define MAX(A,B) ((A) > (B) ? (A) : (B))
#endif

#ifndef MIN
#define MIN(A,B) ((A) < (B) ? (A) : (B))
#endif

int
equal_image_dimensions(nifti_image *nii_ptr, nifti_image *nii_ptr2);

int
write_nifti_double( const char *output_filename, double image[], int data_type, double slope, int dim[], double vox[], nifti_image *in_ptr);

int
write_nifti_float( const char *output_filename, float image[], int data_type, double slope, int dim[], double vox[], nifti_image *in_ptr);

nifti_image
*read_nifti_double( const char *input_filename, double *image[], int read_data);

nifti_image
*read_nifti_float( const char *input_filename, float *image[], int read_data);

#endif
