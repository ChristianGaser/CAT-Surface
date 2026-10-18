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
 * \brief Make sto_xyz the voxel-to-world matrix FSL and nibabel use.
 *
 * libCAT maps voxels to world coordinates with sto_xyz, which nifti1_io sets
 * only from an sform. This uses the sform when sform_code > 0, otherwise the
 * qform, and updates sto_ijk. sform_code is left unchanged. The readers
 * read_nifti_float() and read_nifti_double() apply it; call it yourself on a
 * header from nifti_image_read().
 *
 * \param nii_ptr (in/out) header whose sto_xyz and sto_ijk are set
 */
void set_world_matrix(nifti_image *nii_ptr);
/**
 * \brief Matrix from voxel indices to FSL's scaled voxel coordinates.
 *
 * FSL (flirt, mcflirt, applyxfm4D) states transforms in scaled voxel
 * coordinates: voxel index times voxel size, with the x axis reversed when the
 * voxel-to-world matrix (sto_xyz) has a positive determinant, i.e. for images
 * stored neurologically. FSL coordinates therefore always run radiologically.
 *
 * \param nii_ptr (in)  image header (sto_xyz, nx, dx, dy, dz)
 * \param M       (out) 4x4 row-major matrix, FSL coordinate = M * voxel index
 */
void vox2fsl_matrix(const nifti_image *nii_ptr, double M[16]);
/**
 * \brief Convert a world-space transform between two images to FSL's convention.
 *
 * ras maps world coordinates of the source (input) image to world coordinates
 * of the reference, the direction of a flirt matrix. The result maps FSL
 * coordinates of the source to FSL coordinates of the reference, as flirt
 * -init and -applyxfm read it:
 *   fsl = F_ref * V_ref^-1 * ras * V_src * F_src^-1
 * with V the voxel-to-world (sto_xyz) and F the vox2fsl_matrix() of an image.
 *
 * \param ras (in)  4x4 row-major world transform, source -> reference
 * \param src (in)  header of the source image
 * \param ref (in)  header of the reference image
 * \param fsl (out) 4x4 row-major FSL transform, source -> reference
 */
void ras_to_fsl_matrix(const double ras[16], const nifti_image *src,
                       const nifti_image *ref, double fsl[16]);
/**
 * \brief Convert an FSL transform between two images to world space.
 *
 * Inverse of ras_to_fsl_matrix().
 *
 * \param fsl (in)  4x4 row-major FSL transform, source -> reference
 * \param src (in)  header of the source image
 * \param ref (in)  header of the reference image
 * \param ras (out) 4x4 row-major world transform, source -> reference
 */
void fsl_to_ras_matrix(const double fsl[16], const nifti_image *src,
                       const nifti_image *ref, double ras[16]);
/**
 * \brief Write a 4x4 matrix as FSL writes it (flirt -omat, mcflirt -mats).
 *
 * Four lines of four "%f " values and no comment lines, which FSL's matrix
 * reader would not skip.
 *
 * \param filename (in) output file
 * \param m        (in) 4x4 row-major matrix
 * \return 0 on success, -1 if the file cannot be written
 */
int write_fsl_matrix(const char *filename, const double m[16]);
/**
 * \brief Initialize NIfTI image header structure to reasonable defaults.
 *
 * \param nii_ptr (in/out) NIfTI image structure to initialize
 */
void init_nifti_header(nifti_image *nii_ptr);
/**
 * \brief Write double-precision volume data to NIfTI file with type conversion.
 *
 * The file always holds a single 3D volume, also when the template header comes
 * from a 4D file.
 *
 * \param output_filename (in) path to output .nii/.nii.gz file
 * \param image           (in) linear array of voxel values (size dim[0]*dim[1]*dim[2])
 * \param data_type       (in) target NIfTI datatype (DT_UINT8, DT_INT16, DT_FLOAT32, etc.)
 * \param slope           (in) optional scaling slope; 0.0 for auto-scaling
 * \param dim            (in) array of 3 image dimensions [nx, ny, nz]
 * \param vox            (in) array of 3 voxel spacings [dx, dy, dz]
 * \param in_ptr         (in/out) optional template NIfTI header; NULL to use defaults
 * \return 1 on success; 0 on error (invalid datatype, no extension, write failed)
 */
int write_nifti_double(const char *output_filename, double image[],
                       int data_type, double slope, int dim[], double vox[],
                       nifti_image *in_ptr);
/**
 * \brief Write single-precision volume data to NIfTI file with type conversion.
 *
 * The file always holds a single 3D volume, also when the template header comes
 * from a 4D file.
 *
 * \param output_filename (in) path to output .nii/.nii.gz file
 * \param image           (in) linear array of float voxel values
 * \param data_type       (in) target NIfTI datatype
 * \param slope           (in) optional scaling slope; 0.0 for auto-scaling
 * \param dim            (in) array of 3 image dimensions [nx, ny, nz]
 * \param vox            (in) array of 3 voxel spacings [dx, dy, dz]
 * \param in_ptr         (in/out) optional template NIfTI header; NULL to use defaults
 * \return 1 on success; 0 on error (invalid datatype, no extension, write failed)
 */
int write_nifti_float(const char *output_filename, float image[], int data_type,
                      double slope, int dim[], double vox[],
                      nifti_image *in_ptr);
/**
 * \brief Write a single-precision 4D time series to a NIfTI file.
 *
 * Like write_nifti_float(), but the file holds nt volumes, stored one after the
 * other in image as read_nifti_float() returns them (volume v at v*nx*ny*nz).
 * The time step of a 4D template header is kept, otherwise it is set to 1.
 * Float32 output is written from image without a copy.
 *
 * \param output_filename (in) path to output .nii/.nii.gz file
 * \param image           (in) nt volumes of dim[0]*dim[1]*dim[2] float values
 * \param data_type       (in) target NIfTI datatype
 * \param slope           (in) optional scaling slope; 0.0 for auto-scaling
 * \param dim             (in) array of 3 image dimensions [nx, ny, nz]
 * \param nt              (in) number of volumes (at least 1; 1 writes a 3D file)
 * \param vox             (in) array of 3 voxel spacings [dx, dy, dz]
 * \param in_ptr          (in/out) optional template NIfTI header; NULL to use defaults
 * \return 1 on success; 0 on error (invalid datatype or nt, no extension, write failed)
 */
int write_nifti_float_4d(const char *output_filename, float image[], int data_type,
                         double slope, int dim[], int nt, double vox[],
                         nifti_image *in_ptr);
/**
 * \brief Read NIfTI image file and load data into double-precision array.
 *
 * A 4D file is returned as its nt volumes one after the other: volume v starts at
 * image[v * nvox], where nvox of the returned header is the voxel count of ONE
 * volume (nx*ny*nz) and nt the number of volumes. Non-finite values of float
 * images are read as 0 (the NIfTI library replaces NaN and Inf). sto_xyz of the
 * returned header is the voxel-to-world matrix, also for a qform-only file
 * (see set_world_matrix()).
 *
 * \param input_filename (in)  path to NIfTI file
 * \param image          (out) pointer to allocated double array (size nx*ny*nz*nt)
 * \param read_data      (in)  if non-zero, keep the raw NIfTI data in the returned
 *                             header; if 0, release it (image is returned either way)
 * \return pointer to NIfTI image structure on success; NULL on error
 */
nifti_image *read_nifti_double(const char *input_filename, double *image[],
                               int read_data);
/**
 * \brief Read NIfTI image file and load data into single-precision array.
 *
 * A 4D file is returned as its nt volumes one after the other, as by
 * read_nifti_double(): volume v starts at image[v * nvox], where nvox of the
 * returned header is the voxel count of ONE volume. sto_xyz is the
 * voxel-to-world matrix, also for a qform-only file (see set_world_matrix()).
 *
 * \param input_filename (in)  path to NIfTI file
 * \param image          (out) pointer to allocated float array (size nx*ny*nz*nt)
 * \param read_data      (in)  if non-zero, keep the raw NIfTI data in the returned
 *                             header; if 0, release it (image is returned either way)
 * \return pointer to NIfTI image structure on success; NULL on error
 */
nifti_image *read_nifti_float(const char *input_filename, float *image[],
                              int read_data);

#endif
