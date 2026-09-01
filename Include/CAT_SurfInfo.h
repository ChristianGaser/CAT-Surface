/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 *
 */

#ifndef _CAT_SURFINFO_H_
#define _CAT_SURFINFO_H_

#include <stdio.h>
#include <bicpl.h>

#include "CAT_SurfaceIO.h"

/** Sentinel stored in the self-intersection fields when the test was skipped. */
#define CAT_SURFINFO_SKIPPED (-1)

/**
 * \brief Geometric and topological summary of a surface mesh.
 *
 * Filled by compute_surf_info() and rendered by print_surf_info() or
 * print_surf_info_tabular().  All lengths are in the units of the mesh
 * coordinates (mm for the meshes produced by this package).
 */
typedef struct {
    /* counts */
    int n_points;               /**< number of vertices                      */
    int n_polygons;             /**< number of faces (polygons)              */
    int n_triangles;            /**< faces with exactly three vertices       */
    int n_nontriangular;        /**< faces with more than three vertices     */
    int n_edges;                /**< number of distinct undirected edges     */
    int n_boundary_edges;       /**< edges used by exactly one polygon       */
    int n_nonmanifold_edges;    /**< edges used by more than two polygons    */
    int n_unreferenced_points;  /**< vertices not used by any polygon        */
    int n_degenerate_polygons;  /**< polygons of (numerically) zero area     */

    /* topology */
    int euler;                  /**< Euler characteristic chi = V - E + F    */
    int n_components;           /**< connected components of the mesh graph  */
    int closed;                 /**< 1 if closed and manifold, else 0        */
    int genus;                  /**< (2*components - chi)/2; valid if closed */

    /* geometry */
    double surface_area;        /**< total surface area                      */
    double volume;              /**< signed enclosed volume (closed meshes)  */
    double area_min;            /**< smallest polygon area                   */
    double area_max;            /**< largest polygon area                    */
    double area_mean;           /**< mean polygon area                       */
    double edge_min;            /**< shortest edge length                    */
    double edge_max;            /**< longest edge length                     */
    double edge_mean;           /**< mean edge length                        */
    double bounds[6];           /**< xmin,xmax,ymin,ymax,zmin,zmax           */
    double centroid[3];         /**< area-weighted centroid                  */

    /* self-intersections; CAT_SURFINFO_SKIPPED when the test was not run */
    int n_self_intersections;   /**< connected self-intersecting regions     */
    int n_triangle_hits;        /**< triangle-triangle intersection hits     */
    int n_intersecting_polygons;/**< polygons taking part in an intersection */
} surf_info_struct;

/**
 * \brief Collect geometric and topological information about a mesh.
 *
 * Computes counts, edge topology, Euler characteristic, connected components,
 * surface area, enclosed volume, bounding box and (optionally) the number of
 * self-intersections.  The mesh is not modified.
 *
 * n_triangle_hits is the raw number of hits reported by find_selfintersections()
 * and is what CAT_SurfSelfIntersect prints as "All Triangle Intersections".  It
 * is a count of detections, not of distinct triangle pairs: a pair that spans
 * several cells of the octree is counted once per shared cell.  Use
 * n_self_intersections and n_intersecting_polygons, both of which are derived
 * from the per-polygon labels, when an exact number is needed.
 *
 * \param polygons            (in)  input mesh
 * \param info                (out) filled summary, see surf_info_struct
 * \param check_intersections (in)  if non-zero, run the self-intersection test
 * \param verbose             (in)  if non-zero, print progress information
 * \return void (no return value).
 */
void compute_surf_info(polygons_struct *polygons, surf_info_struct *info,
                       int check_intersections, int verbose);

/**
 * \brief Print a human-readable report of a surf_info_struct.
 *
 * \param info (in) summary filled by compute_surf_info()
 * \param name (in) name shown in the header, may be NULL
 * \param fp   (in) output stream
 * \return void (no return value).
 */
void print_surf_info(const surf_info_struct *info, const char *name, FILE *fp);

/**
 * \brief Print a surf_info_struct as tab-separated key/value lines.
 *
 * Machine-readable counterpart of print_surf_info(), one "key<TAB>value" pair
 * per line.  Fields that were not computed are omitted.
 *
 * \param info (in) summary filled by compute_surf_info()
 * \param name (in) value of the "file" key, may be NULL
 * \param fp   (in) output stream
 * \return void (no return value).
 */
void print_surf_info_tabular(const surf_info_struct *info, const char *name,
                             FILE *fp);

/**
 * \brief Print the DataArrays a GIFTI file carries besides its mesh.
 *
 * A .gii holds the mesh in two DataArrays but may embed per-vertex data next
 * to it -- thickness, curvature, labels -- which the mesh itself does not
 * reveal.  Lists every array with its intent, shape, encoding and value range.
 *
 * \param arrays   (in) descriptions from input_gifti_darrays()
 * \param n_arrays (in) number of entries in arrays
 * \param fp       (in) output stream
 * \return void (no return value).
 */
void print_gifti_darrays(const gifti_darray_info *arrays, int n_arrays,
                         FILE *fp);

/**
 * \brief Print the GIFTI DataArrays as tab-separated key/value lines.
 *
 * \param arrays   (in) descriptions from input_gifti_darrays()
 * \param n_arrays (in) number of entries in arrays
 * \param fp       (in) output stream
 * \return void (no return value).
 */
void print_gifti_darrays_tabular(const gifti_darray_info *arrays, int n_arrays,
                                 FILE *fp);

#endif
