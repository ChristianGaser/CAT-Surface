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

#include <stdlib.h>

#include "3rdparty/nifti/nifti1_io.h"
#include "3rdparty/nifti/nifti1_local.h"

#include <float.h>

#ifndef MAX
#define MAX(A,B) ((A) > (B) ? (A) : (B))
#endif

#ifndef MIN
#define MIN(A,B) ((A) < (B) ? (A) : (B))
#endif

int equal_image_dimensions(nifti_image *nii_ptr, nifti_image *nii_ptr2);
void init_nifti_header(nifti_image *nii_ptr);
int write_nifti_double( const char *output_filename, double image[], int data_type, double slope, int dim[], double vox[], nifti_image *in_ptr);
int write_nifti_float( const char *output_filename, float image[], int data_type, double slope, int dim[], double vox[], nifti_image *in_ptr);
nifti_image *read_nifti_double( const char *input_filename, double *image[], int read_data);
nifti_image *read_nifti_float( const char *input_filename, float *image[], int read_data);

#endif
