/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 *
*/

#ifndef _CAT_SPH_H_
#define _CAT_SPH_H_

#include <bicpl.h>
#include <fftw3.h>

#include "makeweights.h"
#include "FST_semi_fly.h"

#define BINTREE_FACTOR   0.5

/**
 * \brief Read spherical harmonics coefficients from file (x, y, z components).
 *
 * \param file       (in) filename containing SPH coefficients in binary format
 * \param bandwidth  (in) spherical harmonics bandwidth (determines coefficient count)
 * \param rcx,rcy,rcz (out) real SH coefficients for x, y, z components
 * \param icx,icy,icz (out) imaginary SH coefficients for x, y, z components
 * \return 0 on success, 1 on file I/O error
 */
int read_SPHxyz(char *file, int bandwidth, double **rcx, double **rcy,
                double **rcz, double **icx, double **icy, double **icz);
/**
 * \brief Write spherical harmonics coefficients to file (x, y, z components).
 *
 * \param file       (in) output filename for SPH coefficients
 * \param bandwidth  (in) spherical harmonics bandwidth
 * \param rcx,rcy,rcz (in) real SH coefficients for x, y, z components
 * \param icx,icy,icz (in) imaginary SH coefficients for x, y, z components
 * \return file close status (0 on success)
 */
int write_SPHxyz(char *file, int bandwidth, double *rcx, double *rcy,
                 double *rcz, double *icx, double *icy, double *icz);
/**
 * \brief Sample 3D shape from spherical harmonics coefficients onto meshed sphere.
 *
 * \param rdatax,rdatay,rdataz (in) real-valued SH data for x, y, z coordinates
 * \param sphere               (out) output mesh with coordinates from SH reconstruction
 * \param n_triangles          (in)  number of triangles in sphere (for vertex count)
 * \param reparam              (in)  reference sphere topology/remapping structure
 * \param bandwidth            (in)  spherical harmonics bandwidth
 */
void sample_sphere_from_sph(double *rdatax, double *rdatay, double *rdataz,
                            polygons_struct *sphere, int n_triangles,
                            polygons_struct *reparam, int bandwidth);
/**
 * \brief Replace spherical harmonics coefficients with bandwidth-limited version.
 *
 * \param bandwidth         (in) total spherical harmonics bandwidth
 * \param bandwidth_limited (in) maximum degree to copy from filter
 * \param coeffs            (in/out) coefficient array to update
 * \param coeffs_filter     (in) source coefficients for limited bandwidth
 */
void replaceSPH(int bandwidth, int bandwidth_limited, double *coeffs,
                double *coeffs_filter);
/**
 * \brief Compute energy descriptor of 3D shape from SH coefficients per degree.
 *
 * \param bandwidth (in) spherical harmonics bandwidth
 * \param rcx,rcy,rcz (in) real SH coefficients for x, y, z
 * \param icx,icy,icz (in) imaginary SH coefficients for x, y, z
 * \param shape_desc (out) energy per degree l (bandwidth elements)
 */
void shape_description(int bandwidth, double *rcx, double *rcy, double *rcz,
                       double *icx, double *icy, double *icz,
                       double *shape_desc);
/**
 * \brief Apply Butterworth low-pass filter to spherical harmonics coefficients.
 *
 * \param bandwidth         (in) total spherical harmonics bandwidth
 * \param bandwidth_limited (in) cutoff frequency for filter
 * \param coeffs_old        (in) input coefficients
 * \param coeffs_new        (out) filtered coefficients
 */
void butterworth_filter(int bandwidth, int bandwidth_limited,
                        double *coeffs_old, double *coeffs_new);
/**
 * \brief Apply bandpass boxcar filter to spherical harmonics coefficients.
 *
 * \param bandwidth   (in) total spherical harmonics bandwidth
 * \param bw_lo       (in) minimum degree to keep
 * \param bw_hi       (in) maximum degree to keep
 * \param coeffs_old  (in) input coefficients
 * \param coeffs_new  (out) bandpass-filtered coefficients
 */
void bandpass_bandwidth(int bandwidth, int bw_lo, int bw_hi, double *coeffs_old,
                        double *coeffs_new);
/**
 * \brief Apply low-pass boxcar filter to spherical harmonics coefficients.
 *
 * \param bandwidth         (in) total spherical harmonics bandwidth
 * \param bandwidth_limited (in) maximum degree to keep
 * \param coeffs_old        (in) input coefficients
 * \param coeffs_new        (out) bandwidth-limited coefficients
 */
void limit_bandwidth(int bandwidth, int bandwidth_limited, double *coeffs_old,
                     double *coeffs_new);
/**
 * \brief Compute spherical harmonics coefficients from real-valued data using FFT.
 *
 * \param rdata     (in)  real-valued spherical data (bandwidth^2 points)
 * \param bandwidth (in)  spherical harmonics bandwidth (determines grid resolution)
 * \param dataformat (in) format identifier (0=real, 1=complex)
 * \param rc        (out) real parts of SH coefficients
 * \param ic        (out) imaginary parts of SH coefficients
 */
void get_sph_coeffs_of_realdata(double *rdata, int bandwidth, int dataformat,
                                double *rc, double *ic);
/**
 * \brief Reconstruct real-valued data from spherical harmonics coefficients using IFFT.
 *
 * \param rdata     (out) reconstructed real-valued spherical data
 * \param bandwidth (in)  spherical harmonics bandwidth
 * \param dataformat (in) format identifier (0=real, 1=complex)
 * \param rc        (in)  real parts of SH coefficients
 * \param ic        (in)  imaginary parts of SH coefficients
 */
void get_realdata_from_sph_coeffs(double *rdata, int bandwidth, int dataformat,
                                  double *rc, double *ic);
/**
 * \brief Extract equally sampled 2D coordinates from 3D polygonal mesh for SH analysis.
 *
 * \param polygons (in)  source 3D polygonal mesh
 * \param sphere   (in)  reference sphere for coordinate mapping
 * \param bandwidth (in) grid resolution (determines 2D sampling density)
 * \param xcoord,ycoord,zcoord (out) uniform 2D gridded coordinates from mesh
 */
void get_equally_sampled_coords_of_polygon(polygons_struct *polygons,
                                           polygons_struct *sphere,
                                           int bandwidth, double xcoord[],
                                           double ycoord[], double zcoord[]);
/**
 * \brief Maps 3D mesh vertices onto uniformly sampled 2D grid with explicit handling of defined
 *        topological defects. Defect-aware coordinate remapping prevents discontinuities at holes/handles.
 *        Produces regular grid coordinates suitable for SH analysis even with mesh singularities.
 *
 * \param polygons     (in) source 3D polygonal mesh
 * \param defects      (in) per-vertex defect labels
 * \param n_defects    (in) number of defects to separate
 * \param holes        (in) classification of defects as holes (1) or handles (2)
 * \param bandwidth    (in) grid resolution (determines point density)
 * \param xcoord,ycoord,zcoord (out) uniform 2D sampled coordinates
 * \param force        (in) flag for forcing specific handling mode
 */
void get_equally_sampled_coords_holes(polygons_struct *polygons,
                                      polygons_struct *sphere, int *defects,
                                      int n_defects, int *holes, int bandwidth,
                                      double xcoord[], double ycoord[],
                                      double zcoord[], int force);
/**
 * \brief Unit sphere sampled on a regular (theta, phi) grid.
 *
 * \param n_theta (in) number of samples in longitude
 * \param n_phi   (in) number of samples in latitude, poles included
 * \return new object array containing a single POLYGONS object
 */
object_struct ** create_equally_sampled_unit_sphere(int n_theta, int n_phi);
/**
 * \brief Apply 2D Laplacian filter with continuous boundary handling.
 *
 * \param im   (in)  input floating-point 2D image
 * \param msk  (in)  binary mask (1=valid, 0=invalid/masked)
 * \param dimx,dimy (in) image dimensions
 * \param TH   (in)  threshold or tolerance parameter
 * \return Laplacian array (dynamically allocated, same size as input)
 */
double * laplace2d(double *im, unsigned char *msk, int dimx, int dimy,
                   double TH);
/**
 * \brief Compute 2D gradient magnitude with continuous boundary handling.
 *
 * \param im        (in)  input floating-point 2D image
 * \param dimx,dimy (in)  image dimensions (width x height)
 * \return gradient magnitude array (dynamically allocated, same dimensions as input)
 */
double * gradient_magnitude(double *im, int dimx, int dimy);
/**
 * \brief Threshold 2D image to binary (0/1) mask using fixed threshold value.
 *
 * \param im        (in)  input floating-point 2D image
 * \param dimx,dimy (in)  image dimensions (width x height)
 * \param threshold (in)  threshold value for binarization
 * \return binary image as unsigned char array (dynamically allocated)
 */
unsigned char * threshold_image(double *im, int dimx, int dimy,
                                double threshold);
/**
 * \brief Convert linear array index to 2D coordinates for spherical grid.
 *
 * \param i   (in)  linear array index (0-based)
 * \param x   (out) x coordinate (0 to sx-1)
 * \param y   (out) y coordinate (0 to sy-1)
 * \param sxy (in)  total array size (width * height)
 * \param sy  (in)  grid height (number of columns)
 */
void ind2sub2D(int i, int *x, int *y, int sxy, int sy);

#endif
