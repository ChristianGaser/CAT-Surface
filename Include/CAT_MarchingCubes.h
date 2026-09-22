/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 *
 */

#ifndef _CAT_MARCHINGCUBES_H_
#define _CAT_MARCHINGCUBES_H_

#include <bicpl/marching.h>
#include "CAT_NiftiLib.h"
#include "CAT_PpmSulci.h"

#define CHUNK_SIZE    1000000

/** Dark-sheet response a defect must run along before it is cut, not filled. */
#define CAT_TOPO_PRECUT_THRESH 0.3
#define IDX(x, y, z, nx, ny) ((z) * (nx) * (ny) + (y) * (nx) + (x))

/**
 * \brief Correct topological defects in binary volume using Euler characteristic.
 *
 * \param volume        (in/out) floating-point volume (probability or binary)
 * \param thresh        (in)     voxels >= thresh are foreground
 * \param dims          (in)     [nx, ny, nz] volume dimensions
 * \param conn_arr      (in)     two connectivity values, e.g. {18, 26}
 */
void correct_topology(float *volume, float *vol_changed, float thresh,
                      int dims[3], int conn_arr[2]);

/**
 * \brief Extract brain surface mesh with advanced preprocessing and topology correction.
 *
 * \param input_float       (in)  input 3D probabilistic tissue segmentation
 * \param nii_ptr           (in)  NIfTI image header with voxel dimensions and affine
 * \param label             (in)  optional tissue label mask (NULL to skip)
 * \param min_threshold     (in)  isosurface threshold value (typically 0.5 for probabilities)
 * \param pre_fwhm          (in)  Gaussian smoothing FWHM in mm (0 to skip)
 * \param iter_laplacian    (in)  number of Laplacian smoothing iterations
 * \param dist_morph        (in)  distance offset for morphological expansion (mm)
 * \param n_median_filter   (in)  iterations of median filtering to apply
 * \param n_iter            (in)  total outer loop iterations
 * \param strength_gyri_mask (in) weighting factor for gyral preservation masking (0-1)
 * \param sulci_opts        (in)  buried-sulcus correction on the PPM, or NULL to skip
 *                                    it. A buried sulcus is a valley in the PPM whose floor
 *                                    never drops below the isovalue, so the two banks fuse
 *                                    when the isosurface is extracted. No intensity image is
 *                                    needed: the PPM carries the geometry itself, and a
 *                                    Hessian sheetness filter finds the valley. The field is
 *                                    used three times -- to push those floors below the
 *                                    isovalue, to damp the gyral boost above (which would
 *                                    otherwise lift a sulcal floor back over it), and to
 *                                    orient the median filter so it cannot close what was
 *                                    just opened. Note sulci_opts->sheet_strength: the raw
 *                                    response on real data sits well below the thresholds,
 *                                    so a gain of 1 leaves the whole correction inert.
 * \param verbose           (in)  1 to print progress, 0 for silent
 * \return Allocated object_struct containing pial surface polygons; caller must free
 */
object_struct *apply_marching_cubes(float *input_float, nifti_image *nii_ptr,
                                    float *label, double min_threshold,
                                    double pre_fwhm, int iter_laplacian,
                                    double dist_morph, int n_median_filter,
                                    int n_iter, double strength_gyri_mask,
                                    const CAT_PpmSulciOpts *sulci_opts,
                                    double topo_sheet, int verbose);

/**
 * \brief Fast surface mesh extraction with minimal preprocessing.
 *
 * \param input_float     (in)  input 3D volume (typically already binary/thresholded)
 * \param nii_ptr         (in)  NIfTI image header with dimensions and affine transform
 * \param min_threshold   (in)  isosurface threshold value
 * \param iter_laplacian  (in)  number of Laplacian smoothing iterations post-extraction
 * \param verbose         (in)  1 to print progress, 0 for silent operation
 * \return Allocated object_struct containing surface polygons; caller must free
 */
object_struct *apply_marching_cubes_fast(float *input_float,
                                         nifti_image *nii_ptr,
                                         double min_threshold,
                                         int iter_laplacian, int verbose);

/**
 * \brief Extract polygonal surface mesh from volumetric data using marching cubes.
 *
 * \param vol             (in)  input 3D volume as linearized float array
 * \param sizes           (in)  array [nx, ny, nz] dimensions of volume
 * \param min_label       (in)  minimum label value for extraction (-1 to ignore)
 * \param max_label       (in)  maximum label value for extraction (-1 to ignore)
 * \param nii_mat         (in)  NIfTI affine 4x4 matrix (voxel to world coordinates)
 * \param method          (in)  marching cubes algorithm variant
 * \param binary_flag     (in)  treat volume as binary (0/1) if TRUE
 * \param min_threshold   (in)  minimum intensity threshold for surface
 * \param max_threshold   (in)  maximum intensity threshold for surface
 * \param valid_low       (in)  lowest valid data value
 * \param valid_high      (in)  highest valid data value
 * \param polygons        (out) output mesh structure; allocated and populated by function
 * \param verbose         (in)  1 to print progress, 0 for silent operation
 */
void extract_isosurface(float *vol, int sizes[3], double min_label,
                        double max_label, mat44 nii_mat,
                        Marching_cubes_methods method, BOOLEAN binary_flag,
                        double min_threshold, double max_threshold,
                        double valid_low, double valid_high,
                        polygons_struct *polygons, int verbose);

#endif
