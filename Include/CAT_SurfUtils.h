/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 *
 */

#ifndef _CAT_SURFUTILS_H_
#define _CAT_SURFUTILS_H_

#include <bicpl.h>

/*
 * Lightweight surface-utility declarations.
 *
 * This header intentionally contains only small, broadly used helpers that are
 * implemented in Lib/CAT_Surf.c, but are needed by multiple other library
 * modules. It exists to avoid pulling the full umbrella header CAT_Surf.h
 * (and its large include/dependency tree) into modules that only need one or
 * two simple utilities.
 */

#ifndef BINTREE_FACTOR
#define BINTREE_FACTOR 0.5
#endif

#ifndef MAX_NEIGHBOURS
#define MAX_NEIGHBOURS 2000
#endif

/**
 * \brief Apply mixed boundary conditions to grid indices.
 *
 * \param i  input x index (can be negative/out of bounds)
 * \param j  input y index (can be negative/out of bounds)
 * \param dm two-element array with lattice extents \c {nx, ny}
 * \return flattened index in \c [0, nx*ny)
 */
int bound(int i, int j, int dm[]);
/**
 * \brief Axis-aligned bounding box of a mesh.
 *
 * \param polygons (in)  mesh
 * \param bounds   (out) {xmin, xmax, ymin, ymax, zmin, zmax}
 */
void get_bounds(polygons_struct *polygons, double bounds[6]);
/**
 * \brief Set or scale a geometric quantity.
 *
 * \param p (Point *)
 * \param newLength (double)
 */
void set_vector_length(Point *p, double newLength);

/**
 * \brief Translate or align a mesh.
 *
 * \param polygons (polygons_struct *)
 */
void translate_to_center_of_mass(polygons_struct *polygons);
/**
 * \brief Compute or return a derived quantity from the mesh.
 *
 * \param polygons (polygons_struct *)
 * \param area_values (double *)
 * \return See function description for return value semantics.
 */
double get_area_of_polygons(polygons_struct *polygons, double *area_values);
/**
 * \brief Resample per-vertex area values to the sphere and normalize to equal area.
 *
 * \param polygons        source mesh.
 * \param sphere          target spherical mesh (same topology).
 * \param area_values     output array (length \c n_points).
 * \return total surface area of the resampled spherical mesh.
 */
double get_area_of_points_normalized_to_sphere(polygons_struct *polygons,
                                               polygons_struct *sphere,
                                               double *area_values);
/**
 * \brief Multi-stage pipeline to convert a mesh into (increasingly smoothed/inflated) spherical form.
 *
 * \param stop_at stage index (1..5) controlling how far to proceed.
 * \param verbose print stage info and iteration scaling for large meshes.
 */
void surf_to_sphere(polygons_struct *polygons, int stop_at, int verbose);

#endif
