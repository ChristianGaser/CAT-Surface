/* Rachel Yotter - rachel.yotter@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 *
 */

#ifndef _CAT_COMPLEXITY_H_
#define _CAT_COMPLEXITY_H_

#include <bicpl.h>

#define DATAFORMAT 1 /* 1 = double data, 0 = complex data */
#define SPH_ITERS 10
#define BW 1024
#define FWHM 30.0

/**
 * \brief Compute average discrete slope between successive samples.
 *
 * \param x    (in)  x-axis samples
 * \param y    (in)  y-axis samples
 * \param len  (in)  number of samples
 * \return Average slope value
 */
double slope(double *x, double *y, int len);
/**
 * \brief Compute global fractal dimension from log-log slope.
 *
 * \param x    (in)  scale values (e.g., bandwidth or resolution)
 * \param y    (in)  corresponding area ratios
 * \param len  (in)  number of samples
 * \return Global fractal dimension estimate
 */
double get_globalfd(double *x, double *y, int len);
/**
 * \brief Compute per-vertex local fractal dimension values.
 *
 * \param polygons   (in)  surface mesh
 * \param x          (in)  scale values (bandwidths or resolutions)
 * \param areas      (in)  area ratios per scale and polygon
 * \param x_len      (in)  number of scales
 * \param fd         (out) per-vertex fractal dimension values
 * \param smoothflag (in)  non-zero to apply smoothing
 */
void get_localfd(polygons_struct *polygons, double *x, double **areas,
                 int x_len, double *fd, int smoothflag);
/**
 * \brief Select next triangle count for multi-resolution resampling.
 *
 * \param base6   (in/out) base triangle count for 6-subdivision
 * \param base8   (in/out) base triangle count for 8-subdivision
 * \param base20  (in/out) base triangle count for 20-subdivision
 * \return Selected triangle count for this iteration
 */
int min_triangles_update(int *base6, int *base8, int *base20);
/**
 * \brief Compute fractal dimension by multi-resolution surface resampling.
 *
 * \param surface    (in)  input surface mesh
 * \param sphere     (in)  spherical parameterization of the surface
 * \param maxiters   (in)  number of resampling iterations
 * \param file       (in)  output filename for local FD values
 * \param smoothflag (in)  non-zero to smooth local FD values
 * \param debugflag  (in)  non-zero to write debug outputs
 * \return Global fractal dimension estimate
 */
double fractal_dimension(polygons_struct *surface, polygons_struct *sphere,
                         int maxiters, char *file, int smoothflag,
                         int debugflag);
/**
 * \brief Smooth per-vertex values using heat kernel smoothing.
 *
 * \param polygons  (in)  surface mesh
 * \param values    (in/out) per-vertex values to smooth
 * \param fwhm      (in)  full-width at half-maximum in mm
 */
void get_smoothed_values(polygons_struct *polygons, double *values,
                         double fwhm);
/**
 * \brief Compute fractal dimension using spherical harmonics bandwidth reduction.
 *
 * \param surface     (in)  input surface mesh
 * \param sphere      (in)  spherical parameterization
 * \param file        (in)  output filename for local FD values
 * \param n_triangles (in)  triangle count for resampled surfaces
 * \param reparam     (in)  reparameterized sphere for local FD computation
 * \param smoothflag  (in)  non-zero to smooth local FD values
 * \param debugflag   (in)  non-zero to write debug outputs
 * \return Global fractal dimension estimate
 */
double fractal_dimension_sph(polygons_struct *surface, polygons_struct *sphere,
                             char *file, int n_triangles,
                             polygons_struct *reparam, int smoothflag,
                             int debugflag);

#endif
