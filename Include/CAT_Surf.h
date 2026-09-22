/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 *
 */

#ifndef _CAT_SURF_H_
#define _CAT_SURF_H_

#include <stdbool.h>
#include <float.h>
#include "bicpl.h"
#include "CAT_Vol.h"
#include "CAT_Math.h"
#include "CAT_Map.h"
#include "CAT_Smooth.h"
#include "CAT_Resample.h"
#include "CAT_Intersect.h"
#include "CAT_Curvature.h"
#include "CAT_Defect.h"
#include "CAT_Deform.h"
#include "CAT_SurfaceIO.h"
#include "CAT_SurfUtils.h"

#define BW 1024
#define FWHM 30.0
#define DATAFORMAT 1 /* 1 = double data, 0 = complex data */
#define _PI 3.14159265358979323846264338327510

/* edge incidence */
typedef struct { int a,b,c, idx; } keytri_t;
typedef struct {
    int a,b;      /* a<b undirected edge key */
    int tri;      /* triangle index */
    double area;  /* triangle area for ranking */
} edge_occ;

/**
 * \brief Estimate local surface ratio within a spherical neighborhood.
 *
 * \param radius    (in)  neighborhood radius in voxels; if <0, auto-estimate
 * \param polygons  (in)  input surface mesh
 * \param normalize (in)  if non-zero, normalize radius by total surface area
 * \return newly allocated per-vertex surface ratio array
 */
double * get_surface_ratio(double radius, polygons_struct *polygons,
                           int normalize);
/**
 * \brief Compute per-vertex areas using equal polygon distribution.
 *
 * \param polygons      (in)  input mesh
 * \param vertex_areas  (out) per-vertex areas (length n_points)
 * \return total surface area of the mesh
 */
double get_vertex_areas(polygons_struct *polygons, double *vertex_areas);
/**
 * \brief Compute per-vertex area (FreeSurfer-compatible barycentric method).
 *
 * \param polygons    (in)  input mesh
 * \param area_values (out) per-vertex areas (length n_points)
 * \return total surface area of the mesh
 */
double get_area_of_points(polygons_struct *polygons, double *area_values);
/**
 * \brief Compute per-vertex radius values from the origin.
 *
 * \param polygons (in)  input mesh
 * \param radius   (out) per-vertex radii (length n_points)
 */
void get_radius_of_points(polygons_struct *polygons, double *radius);
/**
 * \brief Per-vertex local statistics on a surface (fixed 1-ring).
 *
 * \param polygons  (in)     mesh
 * \param input     (in/out) polygons->n_points per-vertex values
 * \param mask      (in)     optional mask of n_points entries; 0 = skip (NULL: all)
 * \param stat_func (in)     F_MEAN, F_MEDIAN, F_STD, F_MIN or F_MAX
 * \param iters     (in)     number of iterations (>= 1)
 */
void localstat_surface_double(polygons_struct *polygons, double *input,
                              unsigned char *mask, int stat_func, int iters);
/**
 * \brief Mixed boundary condition index mapping for a 2‑D lattice.
 *
 * \param polygons (polygons_struct *)
 * \param target (polygons_struct *)
 */
void correct_bounds_to_target(polygons_struct *polygons,
                              polygons_struct *target);
/**
 * \brief Mixed boundary condition index mapping for a 2‑D lattice.
 *
 * \param polygons (polygons_struct *)
 * \param target (polygons_struct *)
 */
void correct_bounds_to_target_with_scaling(polygons_struct *polygons,
                                           polygons_struct *target);
/**
 * \brief Compute the radius of a sphere enclosing the mesh.
 *
 * \param polygons (in) input mesh
 * \return radius (max vertex distance)
 */
double get_sphere_radius(polygons_struct *polygons);
/**
 * \brief Compute vertex-based Hausdorff distance between meshes.
 *
 * \param p       (in)  first mesh
 * \param p2      (in)  second mesh
 * \param hd      (out) per-vertex distances from p to p2 (length n_points)
 * \param verbose (in)  if non-zero, print summary
 * \return Hausdorff distance (max of forward distances)
 */
double compute_point_hausdorff(polygons_struct *p, polygons_struct *p2,
                               double *hd, int verbose);
/**
 * \brief Compute exact point-wise Hausdorff distance between meshes.
 *
 * \param p       (in)  first mesh
 * \param p2      (in)  second mesh (same vertex count)
 * \param hd      (out) per-vertex distances (length n_points)
 * \param verbose (in)  if non-zero, print summary
 * \return Hausdorff distance (max of hd)
 */
double compute_exact_hausdorff(polygons_struct *p, polygons_struct *p2,
                               double *hd, int verbose);
/**
 * \brief Compute closest-point distances from one mesh to another.
 *
 * \param p       (in)  source mesh
 * \param p2      (in)  target mesh
 * \param hd      (out) per-vertex distances (length n_points)
 * \param verbose (in)  if non-zero, print summary
 * \return mean of hd
 */
double compute_point_distance(polygons_struct *p, polygons_struct *p2,
                              double *hd, int verbose);
/**
 * \brief Compute mean symmetric closest-point distance between two meshes.
 *
 * \param p       (in)  first mesh
 * \param p2      (in)  second mesh (same vertex count)
 * \param dist    (out) per-vertex mean distances (length n_points)
 * \param verbose (in)  if non-zero, print summary
 * \return mean symmetric closest-point distance
 */
double compute_point_distance_mean(polygons_struct *p, polygons_struct *p2,
                                   double *dist, int verbose);
/**
 * \brief Compute Euler characteristic \f$\chi = F + V - E\f$ for a mesh.
 *
 * \param polygons (in) input mesh
 * \param verbose  (in) if non-zero, print duplicate edge info
 * \return Euler characteristic
 */
int euler_characteristic(polygons_struct *polygons, int verbose);
/**
 * \brief Project an ellipsoid-like mesh onto a sphere with a target area.
 *
 * \param polygons           (in/out) mesh to project
 * \param desiredSurfaceArea (in)     desired surface area
 */
void convert_ellipsoid_to_sphere_with_surface_area(polygons_struct *polygons,
                                                   double desiredSurfaceArea);
/**
 * \brief Linear (umbrella) smoothing with optional edge-only passes.
 *
 * \param polygons                   (in/out) mesh to smooth
 * \param strength                   (in)     in (0,1]; larger moves more toward neighbor average
 * \param iters                      (in)     number of iterations
 * \param smoothEdgesEveryXIters     (in)     smooth only on these iterations (0 disables)
 * \param smoothOnlyTheseNodes       (in)     optional mask (length n_points) for selective smoothing
 * \param projectToSphereEveryXIters (in)     project to current sphere radius every X iterations (0 disables)
 */
void linear_smoothing(polygons_struct *polygons, double strength, int iters,
                      int smoothEdgesEveryXIters, int *smoothOnlyTheseNodes,
                      int projectToSphereEveryXIters);
/**
 * \brief Areal smoothing with weights based on local triangle areas.
 *
 * \param polygons                   (in/out) mesh to smooth
 * \param strength                   (in)     smoothing strength
 * \param iters                      (in)     number of iterations
 * \param smoothEdgesEveryXIters     (in)     smooth only on these iterations (0 disables)
 * \param smoothOnlyTheseNodes       (in)     optional mask (length n_points) for selective smoothing
 * \param projectToSphereEveryXIters (in)     project to current sphere radius every X iterations (0 disables)
 */
void areal_smoothing(polygons_struct *polygons, double strength, int iters,
                     int smoothEdgesEveryXIters, int *smoothOnlyTheseNodes,
                     int projectToSphereEveryXIters);
/**
 * \brief Distance-weighted smoothing using Manhattan neighbor distances.
 *
 * \param polygons                   (in/out) mesh to smooth
 * \param strength                   (in)     smoothing strength
 * \param iters                      (in)     number of iterations
 * \param smoothEdgesEveryXIters     (in)     smooth only on these iterations (0 disables)
 * \param smoothOnlyTheseNodes       (in)     optional mask (length n_points) for selective smoothing
 * \param projectToSphereEveryXIters (in)     project to current sphere radius every X iterations (0 disables)
 */
void distance_smoothing(polygons_struct *polygons, double strength, int iters,
                        int smoothEdgesEveryXIters, int *smoothOnlyTheseNodes,
                        int projectToSphereEveryXIters);
/**
 * \brief Inflate a surface while smoothing highly distorted regions.
 *
 * \param polygonsIn           (in/out) mesh to inflate/smooth
 * \param n_smoothingCycles    (in)     number of inflation cycles
 * \param regSmoothStrength    (in)     regular smoothing strength
 * \param regSmoothIters       (in)     regular smoothing iterations
 * \param inflationFactorIn    (in)     inflation factor (>1 inflates)
 * \param compStretchThresh    (in)     distortion threshold for targeted smoothing
 * \param fingerSmoothStrength (in)     strength of targeted smoothing
 * \param fingerSmoothIters    (in)     iterations of targeted smoothing
 */
void inflate_surface_and_smooth_fingers(polygons_struct *polygonsIn,
                                        const int n_smoothingCycles,
                                        const double regSmoothStrength,
                                        const int regSmoothIters,
                                        const double inflationFactorIn,
                                        const double compStretchThresh,
                                        const double fingerSmoothStrength,
                                        const int fingerSmoothIters);                                        
/**
 * \brief Create a new pial/white surface from a central surface without modifying the input.
 *
 * \param polygons         (in) central surface, left unchanged
 * \param thickness_values (in) per-vertex thickness values
 * \param extents          (in) per-vertex displacement multipliers (0.5: pial,
 *                                  -0.5: white)
 * \param check_intersects (in) if non-zero, remove near self-intersections
 * \param sigma            (in) smoothing sigma for the displacement field
 * \param iterations       (in) smoothing iterations
 * \param verbose          (in) verbosity flag
 * \return new object array containing a single POLYGONS object
 */
object_struct ** central_to_new_pial(polygons_struct *polygons,
                                     double *thickness_values, double *extents,
                                     int check_intersects, double sigma,
                                     int iterations, int verbose);
/**
 * \brief Displace a central surface along normals to estimate pial/white surfaces.
 *
 * \param polygons         (in/out) mesh to deform
 * \param thickness_values (in)     per-vertex thickness values
 * \param extents          (in)     per-vertex displacement multipliers
 * \param check_intersects (in)     if non-zero, remove near self-intersections
 * \param sigma            (in)     smoothing sigma for displacement field
 * \param iterations       (in)     smoothing iterations
 * \param verbose          (in)     verbosity flag
 */
void central_to_pial(polygons_struct *polygons, double *thickness_values,
                     double *extents, int check_intersects, double sigma,
                     int iterations, int verbose);
/**
 * \brief Compute per-vertex areas for a thickness-shifted surface.
 *
 * \param polygons         (in)  input mesh
 * \param area             (out) per-vertex areas (length n_points)
 * \param thickness_values (in)  per-vertex thickness values
 * \param extent           (in)  displacement extent
 * \return total surface area of the displaced surface
 */
double get_area_of_points_central_to_pial(polygons_struct *polygons,
                                          double *area,
                                          double *thickness_values,
                                          double extent);
/**
 * \brief Correct folded mesh regions using curvature-driven displacement.
 *
 * \param polygons           (in/out) mesh to correct
 * \param polygons_reference (in)     reference mesh (can be NULL if volume is used)
 * \param vol                (in)     volume data (can be NULL if reference mesh is used)
 * \param nii_ptr            (in)     NIfTI header for sampling (required if vol is used)
 * \param isovalue           (in)     target isovalue in volume
 * \return EXIT_SUCCESS on success, EXIT_FAILURE on failure
 */
int correct_mesh_folding(polygons_struct *polygons,
                         polygons_struct *polygons_reference, float *vol,
                         nifti_image *nii_ptr, double isovalue);
/**
 * \brief Reduce mesh complexity using Quadric Error Metrics (QEM).
 *
 * \param polygons       (in/out) mesh to simplify
 * \param target_faces   (in)     desired triangle count (<=0 uses half)
 * \param aggressiveness (in)     simplifier aggressiveness (larger => stronger)
 * \param preserve_sharp (in)     if non-zero, preserve sharp features
 * \param verbose        (in)     if non-zero, print progress
 * \return 0 on success, -1 on failure
 */
int reduce_mesh_quadrics(polygons_struct *polygons, int target_faces,
                         double aggressiveness, int preserve_sharp,
                         int verbose);

#endif
