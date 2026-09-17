/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "CAT_SurfPialWhite.h"
#include "CAT_Surf.h"
#include "CAT_Vol.h"
#include "CAT_Smooth.h"
#include "CAT_Deform.h"
#include "CAT_Curvature.h"
#include "CAT_SurfLaplacian.h"
#include "CAT_Intersect.h"

/* Tissue class thresholds */
#ifndef CGM
#define CGM 1.5
#endif
#ifndef GWM
#define GWM 2.5
#endif

/**
 * \brief Initialize pial/white estimation options with defaults.
 *
 * \param opts (out) options structure to initialize
 * \return void
 */
void CAT_PialWhiteOptionsInit(CAT_PialWhiteOptions *opts)
{
    if (!opts)
        return;
    opts->w1 = 0.05;
    opts->w2 = 0.05;
    opts->w3 = 0.1;
    opts->sigma = 0.2;
    opts->iterations = 200;
    opts->method = 2;
    opts->pial_profile = 1;
    CAT_PialProfileOptionsInit(&opts->profile);
    opts->profile.isovalue = CGM;
    opts->remove_intersect = 0;
    opts->verbose = 0;
}

/**
 * \brief Estimate pial and white surfaces from a central surface.
 *
 * Creates initial pial/white surfaces from thickness values, applies
 * curvature-guided smoothing to the pial surface and deforms the white
 * surface driven by tissue labels and gradients.  The pial surface is then
 * placed by profile search (opts->pial_profile), which reaches the CSF/GM
 * boundary at gyral crowns and stops at the valley bottom of glued sulci,
 * and is then smoothed with 2 iterations of HC Laplacian smoothing to remove
 * the per-vertex noise of the placement.  With pial_profile = 0 it is
 * deformed jointly with the white surface.
 *
 * \param central          (in)  central surface mesh
 * \param thickness_values (in)  per-vertex thickness values
 * \param labels           (in)  tissue label volume
 * \param nii_ptr          (in)  NIfTI header for coordinate transforms
 * \param pial_out         (out) pial surface mesh
 * \param white_out        (out) white surface mesh
 * \param opts             (in)  algorithm options
 * \return 0 on success, non-zero on error
 */
int CAT_SurfEstimatePialWhite(
    polygons_struct *central,
    const double *thickness_values,
    float *labels,
    nifti_image *nii_ptr,
    polygons_struct *pial_out,
    polygons_struct *white_out,
    const CAT_PialWhiteOptions *opts)
{
    int p;
    int n_points;
    double *extents = NULL;
    double *weight = NULL;
    int *n_neighbours = NULL;
    int **neighbours = NULL;
    object_struct **objects_out;
    polygons_struct *polygons_pial = NULL;
    polygons_struct *polygons_white = NULL;
    polygons_struct *polygons_smoothed = NULL;
    double weights[3];
    /* Target offsets of the balloon deformation: the pial one only applies to
     * the legacy pial path.  The white offset compensates the inward bias of
     * the deformation; with the ADE start 0.1 gives the lowest white error
     * (0.2 put the surface 0.10 mm inside the GM/WM boundary). */
    double shifting[2] = {-0.25, 0.1};

    if (!central || !thickness_values || !labels || !nii_ptr ||
        !pial_out || !white_out || !opts)
        return -1;

    n_points = central->n_points;

    /* ------ ADE start surfaces: both (method 1) or white only (method 2) ------ */
    if (opts->method > 0)
    {
        /* method 2 starts the pial surface from the thickness, so its
         * streamlines would be traced for nothing */
        int rc = surf_ade_pial_white(
            central, labels, nii_ptr,
            CGM, GWM,
            thickness_values,
            opts->method == 1 ? pial_out : NULL, white_out, opts->verbose);
        if (rc != 0)
            return rc;

        /* Continue with the standard deformation pipeline using
         * the ADE surfaces as starting surfaces. */
        if (opts->method == 1)
            polygons_pial = pial_out;
        polygons_white = white_out;
    }

    /* ------ Deformation method (method == 0, default) ------ */

    /* Allocate working arrays */
    extents = (double *)malloc(sizeof(double) * n_points);
    weight = (double *)malloc(sizeof(double) * n_points);
    polygons_smoothed = (polygons_struct *)malloc(sizeof(polygons_struct));

    if (!extents || !weight || !polygons_smoothed)
    {
        if (extents)
            free(extents);
        if (weight)
            free(weight);
        if (polygons_smoothed)
            free(polygons_smoothed);
        return -2;
    }

    /* Get neighbours for curvature computation */
    get_all_polygon_point_neighbours(central, &n_neighbours, &neighbours);

    /* Compute curvature-based weights (negative mean curvature -> smoothing) */
    get_polygon_vertex_curvatures_cg(central, n_neighbours, neighbours, 3.0, 0.0, weight);
    for (p = 0; p < n_points; p++)
    {
        weight[p] = fmin(0.0, weight[p]);   /* Only negative curvatures */
        weight[p] = fmax(-90.0, weight[p]); /* Clip at -90 */
        weight[p] /= -90.0;                 /* Normalize to [0..1] */
    }

    /* Initial estimate of pial surface */
    if (opts->method == 0 || opts->method == 2)
    {
        for (p = 0; p < n_points; p++)
            extents[p] = 0.5;
        objects_out = central_to_new_pial(central, (double *)thickness_values, extents,
                                          1, 0.5 * opts->sigma, 5, opts->verbose);
        polygons_pial = get_polygons_ptr(objects_out[0]);
    }

    /* Smooth pial surface based on local curvature */
    copy_polygons(polygons_pial, polygons_smoothed);
    smooth_heatkernel(polygons_smoothed, NULL, 5.0);

    /* Blend original and smoothed using curvature weights */
    for (p = 0; p < n_points; p++)
    {
        Point_x(polygons_pial->points[p]) = weight[p] * Point_x(polygons_smoothed->points[p]) +
                                            (1.0 - weight[p]) * Point_x(polygons_pial->points[p]);
        Point_y(polygons_pial->points[p]) = weight[p] * Point_y(polygons_smoothed->points[p]) +
                                            (1.0 - weight[p]) * Point_y(polygons_pial->points[p]);
        Point_z(polygons_pial->points[p]) = weight[p] * Point_z(polygons_smoothed->points[p]) +
                                            (1.0 - weight[p]) * Point_z(polygons_pial->points[p]);
    }

    /* Initial estimate of white surface */
    if (opts->method == 0)
    {
        for (p = 0; p < n_points; p++)
            extents[p] = -0.5;
        objects_out = central_to_new_pial(central, (double *)thickness_values, extents,
                                          0, 0.5 * opts->sigma, 5, opts->verbose);
        polygons_white = get_polygons_ptr(objects_out[0]);
    }

    if (!polygons_pial || !polygons_white)
    {
        free(extents);
        free(weight);
        free(polygons_smoothed);
        return -3;
    }
    
    /* Deformation: white surface always, pial surface only in legacy mode */
    weights[0] = opts->w1;
    weights[1] = opts->w2;
    weights[2] = opts->w3;
    surf_deform_dual(opts->pial_profile ? NULL : polygons_pial, polygons_white,
                     central, labels, nii_ptr,
                     weights, opts->sigma, CGM + shifting[0], GWM + shifting[1],
                     (double *)thickness_values, opts->iterations, opts->verbose);

    /* Profile-based pial placement.  It starts from the thickness-based
     * estimate and does not need the balloon-force deformation above. */
    if (opts->pial_profile)
    {
        CAT_PialProfileOptions profile = opts->profile;
        profile.verbose = opts->verbose;
        if (CAT_SurfDeformPialProfile(polygons_pial, labels, nii_ptr, &profile) != 0)
        {
            free(extents);
            free(weight);
            free(polygons_smoothed);
            return -4;
        }

        /* The placement is per vertex and leaves the surface much noisier
         * than the central surface.  Two iterations of HC Laplacian smoothing
         * bring the roughness back to that of the central surface while the
         * label at the vertices stays unchanged on average.  The white
         * surface is already smoothed at the end of surf_deform_dual(). */
        smooth_laplacian(polygons_pial, 2, 0.1, 0.5);
    }

    /* Copy results to output.
     * In method 1, polygons_pial/polygons_white may already be pial_out/white_out. */
    if (polygons_pial != pial_out)
        copy_polygons(polygons_pial, pial_out);
    if (polygons_white != white_out)
        copy_polygons(polygons_white, white_out);

    /* Optionally repair self-intersections of both surfaces.  The pial surface
     * is pushed outwards into tight sulci, and the two sides of a thin gyral
     * blade can be driven through each other on either surface -- inside the
     * blade the white target is never reached.  Local smoothing cannot
     * separate such crossed sheets, so defects that survive it retreat
     * towards the central surface both surfaces started from.  The repair is
     * topology preserving, so the vertex correspondence between central, pial
     * and white surfaces - and with it the per-vertex thickness - is kept. */
    if (opts->remove_intersect)
    {
        if (opts->verbose)
            fprintf(stdout, "Remove self-intersections of pial surface\n");
        remove_intersections_ref(pial_out, central->points, 10, 50, opts->verbose);

        if (opts->verbose)
            fprintf(stdout, "Remove self-intersections of white surface\n");
        remove_intersections_ref(white_out, central->points, 10, 50, opts->verbose);
    }

    /* Cleanup */
    free(extents);
    free(weight);
    free(polygons_smoothed);
//    if (n_neighbours && neighbours)
//        delete_polygon_point_neighbours(central, n_neighbours, neighbours, NULL, NULL);

    return 0;
}
