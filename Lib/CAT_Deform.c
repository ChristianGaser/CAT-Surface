/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry, University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 *
 */
#include <bicpl.h>

#include "CAT_Math.h"
#include "CAT_NiftiLib.h"
#include "CAT_Smooth.h"
#include "CAT_Vol.h"
#include "CAT_Intersect.h"
#include "CAT_MeshClean.h"
#include "CAT_Curvature.h"

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/**
 * \brief Smooths a 3D displacement field with Jacobian- and curvature-based blending.
 *
 * This function performs iterative smoothing of a displacement field associated with a surface mesh.
 * Blending weights increase smoothing in regions with small Jacobian determinants (non-diffeomorphic
 * risk) and optionally in regions with a selected curvature sign (e.g., sulci vs gyri), controlled by
 * a global curvature blending weight.
 *
 * \param displacement_field  2D array [n_points][3] of displacements per vertex (in/out).
 * \param polygons            Pointer to the surface mesh (polygons_struct).
 * \param n_neighbours        Array with number of neighbors per vertex.
 * \param neighbours          2D ragged array of neighbor indices per vertex.
 * \param iterations          Number of smoothing iterations to perform.
 * \param sigma               Smoothing decay factor applied to the neighbor average.
 * \param min_det             Minimum acceptable Jacobian determinant for diffeomorphism.
 * \param blend_strength      Logistic sharpness for Jacobian blending (higher = sharper).
 * \param curvature           Optional curvature array (size = polygons->n_points). Pass NULL to disable.
 * \param curvature_sign      Curvature sign selection: -1 = use only negative, +1 = only positive,
 *                            0 = absolute (both). Values outside {-1,0,1} are clamped to 0.
 * \param curvature_weight    Global scaling for curvature-based blending contribution. Positive
 *                            values increase smoothing in selected curvature regions; negative
 *                            values decrease smoothing there. Magnitudes >1 are clipped by
 *                            the final [0..1] clamp of the blend.
 */
void smooth_displacement_field_blended(double (*displacement_field)[3],
                                       polygons_struct *polygons,
                                       int *n_neighbours, int **neighbours,
                                       int iterations, double sigma,
                                       double min_det, double blend_strength,
                                       const double *curvature, int curvature_sign,
                                       double curvature_weight)
{
    int it, v, j, k, pidx;
    double (*new_disp)[3] = malloc(sizeof(double[3]) * polygons->n_points);
    const double curv_clip = 90.0; // Clip magnitude before normalization

    // Normalize curvature_sign to {-1,0,1}
    if (curvature_sign > 0)
        curvature_sign = 1;
    else if (curvature_sign < 0)
        curvature_sign = -1;
    else
        curvature_sign = 0;

    for (it = 0; it < iterations; it++)
    {
        for (v = 0; v < polygons->n_points; v++)
        {
            double J[3][3], detJ;

            // Estimate local Jacobian using 3 neighboring displacements
            for (k = 0; k < 3; k++)
            {
                J[k][0] = displacement_field[v][k] - displacement_field[neighbours[v][0]][k];
                J[k][1] = displacement_field[v][k] - displacement_field[neighbours[v][1]][k];
                J[k][2] = displacement_field[v][k] - displacement_field[neighbours[v][2]][k];
            }

            detJ = J[0][0] * (J[1][1] * J[2][2] - J[1][2] * J[2][1]) -
                   J[0][1] * (J[1][0] * J[2][2] - J[1][2] * J[2][0]) +
                   J[0][2] * (J[1][0] * J[2][1] - J[1][1] * J[2][0]);

            // Jacobian-based blend: 1 = smooth fully, 0 = keep original
            double blend_jac = 1.0 / (1.0 + exp(blend_strength * (detJ - min_det)));

            // Curvature-based blend in [0..1]
            double blend_curv = 0.0;
            if (curvature && curvature_weight != 0.0)
            {
                double c = curvature[v];
                if (curvature_sign < 0)
                {
                    // Only negative curvature: clip to [-curv_clip, 0] and map to [0..1]
                    c = fmin(0.0, c);
                    c = fmax(-curv_clip, c);
                    blend_curv = (-c) / curv_clip;
                }
                else if (curvature_sign > 0)
                {
                    // Only positive curvature: clip to [0, curv_clip] and map to [0..1]
                    c = fmax(0.0, c);
                    c = fmin(curv_clip, c);
                    blend_curv = c / curv_clip;
                }
                else
                {
                    // Absolute curvature: clip to [0, curv_clip] and map to [0..1]
                    c = fabs(c);
                    c = fmin(curv_clip, c);
                    blend_curv = c / curv_clip;
                }
            }

            // Combine blends; clamp to [0..1]
            double blend = blend_jac + curvature_weight * blend_curv;
            if (blend < 0.0)
                blend = 0.0;
            if (blend > 1.0)
                blend = 1.0;

            double smoothed[3] = {0.0, 0.0, 0.0};
            int count = 0;

            for (j = 0; j < n_neighbours[v]; j++)
            {
                pidx = neighbours[v][j];
                for (k = 0; k < 3; k++)
                {
                    smoothed[k] += displacement_field[pidx][k];
                }
                count++;
            }

            for (k = 0; k < 3; k++)
            {
                smoothed[k] = (smoothed[k] / count) * exp(-sigma);
                // Blend between original and smoothed displacement
                new_disp[v][k] = (1.0 - blend) * displacement_field[v][k] + blend * smoothed[k];
            }
        }

        // Commit updates
        memcpy(displacement_field[0], new_disp[0], sizeof(double[3]) * polygons->n_points);
    }

    free(new_disp);
}

/**
 * \brief Smooth a displacement field using neighborhood averaging with exponential decay.
 *
 * Iteratively blurs a per-vertex displacement field by averaging over vertex neighborhoods,
 * then attenuating by exponential factor exp(-sigma). Used to regularize surface deformations
 * and prevent erratic vertex motion during iterative mesh deformation.
 *
 * Algorithm:
 *  1. For each vertex: accumulate displacements from all neighbors
 *  2. Average and apply exponential attenuation exp(-sigma)
 *  3. Repeat for specified number of iterations
 *
 * \param displacement_field (in/out) double[n_points][3]; per-vertex displacement vectors
 * \param polygons            (in)    surface mesh
 * \param n_neighbours        (in)    int[n_points]; vertex degree (number of neighbors)
 * \param neighbours          (in)    int*[n_points]; neighbor indices per vertex
 * \param iterations          (in)    number of smoothing passes
 * \param sigma               (in)    exponential decay parameter of attenuation
 */
void smooth_displacement_field(double (*displacement_field)[3], polygons_struct *polygons,
                               int *n_neighbours, int **neighbours, int iterations, double sigma)
{
    int v, j, k, it, pidx;
    double (*new_displacement)[3] = malloc(sizeof(double[3]) * polygons->n_points);

    for (it = 0; it < iterations; it++)
    {
        for (v = 0; v < polygons->n_points; v++)
        {
            double smoothed[3] = {0.0, 0.0, 0.0};
            int count = 0;

            for (j = 0; j < n_neighbours[v]; j++)
            {
                pidx = neighbours[v][j];
                for (k = 0; k < 3; k++)
                    smoothed[k] += displacement_field[pidx][k];
                count++;
            }

            for (k = 0; k < 3; k++)
                new_displacement[v][k] = (smoothed[k] / count) * exp(-sigma);
        }

        // Update displacement field
        memcpy(displacement_field[0], new_displacement[0], sizeof(double[3]) * polygons->n_points);
    }

    free(new_displacement);
}

/**
 * \brief Deform a surface mesh toward intensity gradients from an external volume.
 *
 * Iteratively moves mesh vertices using active contour principles: balancing internal
 * smoothness constraints with external image gradient forces. Optionally checks and
 * removes self-intersections using a spatial voxel grid to prevent mesh degeneracy.
 *
 * Algorithm:
 *  1. Compute per-voxel intensity gradients from input volume
 *  2. For each iteration:
 *  3.   Compute local smoothness forces from vertex neighborhoods
 *  4.   Accumulate external forces from image gradients
 *  5.   Move vertices along combined force direction
 *  6.   Check and correct self-intersections via grid-based collision detection
 *
 * \param polygons            (in/out) surface mesh (modified in-place)
 * \param input               (in)     float[nvoxels]; intensity volume data
 * \param nii_ptr             (in)     NIfTI header with volume dimensions and voxel size
 * \param w                   (in)     double[3]; weight factors {smoothness, gradient_strength, ?}
 * \param sigma               (in)     Gaussian smoothing parameter for displacement field
 * \param lim                 (in)     intensity threshold controlling deformation magnitude
 * \param it                  (in)     number of deformation iterations
 * \param remove_intersections (in)    boolean; enable self-intersection removal
 * \param verbose             (in)     boolean; print iteration progress if true
 */
void surf_deform(polygons_struct *polygons, float *input, nifti_image *nii_ptr,
                 double w[3], double sigma, float lim, int it, int remove_intersections, int verbose)
{
    int i, j, k, v, dims[3], nvox, pidx;
    int *n_neighbours, **neighbours;
    float *gradient_x, *gradient_y, *gradient_z;
    double vx[3], s;
    Point points[MAX_POINTS_PER_POLYGON];
    polygons_struct *polygons_orig;
    object_struct *orig_object;

    orig_object = create_object(POLYGONS);
    polygons_orig = get_polygons_ptr(orig_object);
    copy_polygons(polygons, polygons_orig);

    // Extract image dimensions and voxel size
    dims[0] = nii_ptr->nx;
    dims[1] = nii_ptr->ny;
    dims[2] = nii_ptr->nz;
    nvox = dims[0] * dims[1] * dims[2];

    vx[0] = nii_ptr->dx;
    vx[1] = nii_ptr->dy;
    vx[2] = nii_ptr->dz;

    // Allocate memory for gradient images
    gradient_x = (float *)malloc(sizeof(float) * nvox);
    gradient_y = (float *)malloc(sizeof(float) * nvox);
    gradient_z = (float *)malloc(sizeof(float) * nvox);

    if (!gradient_x || !gradient_y || !gradient_z)
    {
        fprintf(stderr, "Memory allocation error\n");
        exit(EXIT_FAILURE);
    }

    // Compute gradient of the input volume
    gradient3D(input, NULL, gradient_x, gradient_y, gradient_z, dims, vx);

    // Compute surface normals and neighbors
    compute_polygon_normals(polygons);
    create_polygon_point_neighbours(polygons, TRUE, &n_neighbours, &neighbours, NULL, NULL);

    // Allocate displacement field
    double (*displacement_field)[3] = malloc(sizeof(double[3]) * polygons->n_points);
    if (!displacement_field)
    {
        fprintf(stderr, "Memory allocation error\n");
        exit(EXIT_FAILURE);
    }

    // Iterative deformation process
    int counter = 0;
    double s_prev = FLT_MAX;
    for (i = 0; i < it; i++)
    {
        s = 0.0;

        for (v = 0; v < polygons->n_points; v++)
        {
            // Compute centroid of neighboring vertices for smoothing
            double c[3] = {0.0, 0.0, 0.0};
            for (j = 0; j < n_neighbours[v]; j++)
            {
                pidx = neighbours[v][j];
                c[0] += Point_x(polygons->points[pidx]);
                c[1] += Point_y(polygons->points[pidx]);
                c[2] += Point_z(polygons->points[pidx]);
            }
            for (k = 0; k < 3; k++)
            {
                c[k] /= (float)n_neighbours[v];
            }

            // Get current vertex position
            double p[3] = {Point_x(polygons->points[v]),
                           Point_y(polygons->points[v]),
                           Point_z(polygons->points[v])};

            // Get vertex normal
            double n[3] = {Point_x(polygons->normals[v]),
                           Point_y(polygons->normals[v]),
                           Point_z(polygons->normals[v])};

            // Compute external force based on image gradient
            float di = isoval(input, p[0], p[1], p[2], dims, nii_ptr) - lim;
            float fx = isoval(gradient_x, p[0], p[1], p[2], dims, nii_ptr);
            float fy = isoval(gradient_y, p[0], p[1], p[2], dims, nii_ptr);
            float fz = isoval(gradient_z, p[0], p[1], p[2], dims, nii_ptr);
            float f3 = ((di / 1.0));
            float f2 = fmax(-1.0, fmin(1.0, fx * n[0] + fy * n[1] + fz * n[2]));

            // Dynamic boosting (optional: limit max)
            float boost = 1.0 + tanh(fabs(di));
            float w3_scaled = fmin(w[2] * boost, 5.0 * w[2]);

            // Compute vertex displacement and store in the displacement field
            for (k = 0; k < 3; k++)
            {
                displacement_field[v][k] = w[0] * (c[k] - p[k]) +
                                           ((w[1] * f2 + w3_scaled) * f3) * n[k];
            }

            s += di * di;
        }

        // Stop if no further minimization is seen and s increases instead
        if (s > s_prev)
            counter++;
        if (counter > 10)
            break;

        // Apply smoothing to the displacement field
        smooth_displacement_field(displacement_field, polygons, n_neighbours, neighbours, 5, sigma);

        // Apply the final displacement to vertices
        for (v = 0; v < polygons->n_points; v++)
        {
            Point_x(polygons->points[v]) += displacement_field[v][0];
            Point_y(polygons->points[v]) += displacement_field[v][1];
            Point_z(polygons->points[v]) += displacement_field[v][2];
        }

        int n_self_hits = 0;
        int *flags = find_near_self_intersections(polygons, 0.75, &n_self_hits);
        for (v = 0; v < polygons->n_points; v++)
        {
            if (flags[v])
            {
                Point_x(polygons->points[v]) -= displacement_field[v][0];
                Point_y(polygons->points[v]) -= displacement_field[v][1];
                Point_z(polygons->points[v]) -= displacement_field[v][2];
            }
        }

        // Update normals for next iteration
        compute_polygon_normals(polygons);
        if (verbose)
        {
            fprintf(stdout, "\rMesh: deform: iter %03d | Error: %6.4f", i + 1,
                    sqrt(s / polygons->n_points));
            fflush(stdout); // Force output update
        }
        s_prev = s;
    }
    if (verbose)
        fprintf(stdout, "\n");

    for (v = 0; v < polygons->n_points; v++)
    {
        displacement_field[v][0] = Point_x(polygons->points[v]) - Point_x(polygons_orig->points[v]);
        displacement_field[v][1] = Point_y(polygons->points[v]) - Point_y(polygons_orig->points[v]);
        displacement_field[v][2] = Point_z(polygons->points[v]) - Point_z(polygons_orig->points[v]);
    }

    // Get the squared sum of displacement
    double *displacement = malloc(sizeof(double) * polygons->n_points);

    double prctile[2] = {95.0, 95.0};

    /* Get percentile for x-displacement */
    for (v = 0; v < polygons->n_points; v++)
        displacement[v] = displacement_field[v][0];
    double threshold_x[2];
    get_prctile(displacement, polygons->n_points, threshold_x, prctile, 1, DT_FLOAT64);

    /* Get percentile for y-displacement */
    for (v = 0; v < polygons->n_points; v++)
        displacement[v] = displacement_field[v][1];
    double threshold_y[2];
    get_prctile(displacement, polygons->n_points, threshold_y, prctile, 1, DT_FLOAT64);

    /* Get percentile for z-displacement */
    for (v = 0; v < polygons->n_points; v++)
        displacement[v] = displacement_field[v][2];
    double threshold_z[2];
    get_prctile(displacement, polygons->n_points, threshold_z, prctile, 1, DT_FLOAT64);

    for (v = 0; v < polygons->n_points; v++)
    {
        displacement[v] = SQR(displacement_field[v][0]) +
                          SQR(displacement_field[v][1]) +
                          SQR(displacement_field[v][2]);
    }

    /* Get percentile for squared sum of displacements */
    double threshold[2];
    get_prctile(displacement, polygons->n_points, threshold, prctile, 1, DT_FLOAT64);

    /* If squared sum of displacements is exceeding 95% percentile we limit all
      displacements to 95% percentile to prevent that these outliers cause
      self-intersections */
    for (v = 0; v < polygons->n_points; v++)
    {
        if (displacement[v] > threshold[1])
        {
            displacement_field[v][0] = threshold_x[1];
            displacement_field[v][1] = threshold_y[1];
            displacement_field[v][2] = threshold_z[1];
        }
    }

    // Apply very slight smoothing to the displacement field to smooth replacements
    smooth_displacement_field(displacement_field, polygons, n_neighbours, neighbours, 5, 0.05 * sigma);

    // Apply the final displacement to vertices
    for (v = 0; v < polygons->n_points; v++)
    {
        Point_x(polygons_orig->points[v]) += displacement_field[v][0];
        Point_y(polygons_orig->points[v]) += displacement_field[v][1];
        Point_z(polygons_orig->points[v]) += displacement_field[v][2];
    }
    copy_polygons(polygons_orig, polygons);

    // Use MeshFix approach to remove self-intersections
    if (remove_intersections)
    {
        if (verbose)
            fprintf(stdout, "\n");
        CAT_MeshCleanOptions opts;
        CAT_MeshCleanOptionsInit(&opts);
        int result = CAT_SurfMeshClean(polygons, &opts);
    }

    // Free allocated memory
    free(gradient_x);
    free(gradient_y);
    free(gradient_z);
    free(displacement);
    free(displacement_field);
    delete_polygon_point_neighbours(polygons, n_neighbours, neighbours, NULL, NULL);
}
/**
 * \brief Jointly deforms two homologous surface meshes (pial and white) with image forces,
 *        distance constraint, and curvature-/Jacobian-blended smoothing.
 *
 * This routine performs a coupled deformation of two meshes that share topology and vertex
 * correspondence (e.g., pial and white matter surfaces). For each vertex, the total
 * displacement combines:
 *  - Internal smoothing (neighbor centroid attraction)
 *  - External image forces (gradient alignment and balloon force based on intensity offset)
 *  - A distance constraint to maintain target inter-surface spacing
 *
 * Displacement fields are iteratively smoothed using a blended scheme that increases smoothing
 * where local Jacobians suggest non-diffeomorphic risk, and modulates smoothing by curvature:
 *  - polygons1: less smoothing in positive curvature (gyri) regions
 *  - polygons2: less smoothing in negative curvature (sulci) regions
 * Curvature is computed once from the provided baseline mesh and reused each iteration.
 *
 * Parameters
 *  - polygons1, polygons2: Input/output meshes to be jointly deformed (same topology).
 *  - polygons_orig: Baseline mesh used to initialize both surfaces and compute curvature.
 *                   Must be non-NULL and share topology with polygons1/2.
 *  - input: Pointer to the scalar 3D image used to derive external forces.
 *  - nii_ptr: NIfTI image header (provides dimensions and voxel size).
 *  - w[4]: Weights for force terms:
 *      w[0]: internal smoothing weight
 *      w[1]: gradient alignment weight (edge attraction along normal)
 *      w[2]: balloon force weight (proportional to intensity difference)
 *  - sigma: Smoothing decay factor for displacement field averaging.
 *  - lim1, lim2: Isovalue targets for pial and white surfaces (controls balloon term).
 *  - target_distance: Array of desired distances between corresponding vertices on polygons1/2.
 *  - it: Number of deformation iterations.
 *  - verbose: If non-zero, prints per-iteration error metrics.
 *
 * Behavior and safeguards
 *  - External forces are limited to avoid instabilities (e.g., clamped normal alignment).
 *  - Balloon term is adaptively boosted when far from the isosurface, with an upper cap.
 *  - Blended smoothing uses a logistic of the Jacobian determinant and curvature-dependent
 *    weights (positive curvature → less smoothing for polygons1; negative for polygons2).
 *  - Self-intersections are detected and locally counteracted each iteration.
 *
 * Note
 *  This function expects polygons_orig to be provided as a shared baseline. Curvature is
 *  taken from this geometry once (fixed across iterations) because polygons1 and polygons2
 *  are similarly folded, which reduces compute and improves stability.
 */
/**
 * \brief Simultaneously deform two surfaces (e.g., white+pial) toward image gradients.
 *
 * Extends surf_deform to maintain consistent spacing between two surfaces (white matter
 * and pial surface) while both deform toward intensity features. Both surfaces start
 * from the same reference surface and deform with potentially different intensity thresholds,
 * while maintaining a target distance between them.
 *
 * Algorithm:
 *  1. Copy reference surface to both output surfaces
 *  2. For each iteration: deform both surfaces toward image gradients
 *  3. Optionally constrain distance between surfaces to match cortical thickness
 *  4. Check and correct self-intersections in both surfaces
 *  5. Update target distance array based on achieved surface separation
 *
 * \param polygons1         (in/out) first deformed surface (e.g., pial matter)
 * \param polygons2         (in/out) second deformed surface (e.g., white surface)
 * \param polygons_orig     (in)     reference surface (starting template)
 * \param input             (in)     float[nvoxels]; intensity volume data
 * \param nii_ptr           (in)     NIfTI header with volume dimensions and voxel size
 * \param w                 (in)     double[3]; weight factors for deformation forces
 * \param sigma             (in)     Gaussian smoothing parameter for displacement field
 * \param lim1              (in)     intensity threshold for first surface deformation
 * \param lim2              (in)     intensity threshold for second surface deformation
 * \param target_distance   (in/out) double[n_points]; desired separation between surfaces
 * \param it                (in)     number of deformation iterations
 * \param verbose           (in)     boolean; print iteration progress if true
 */
void surf_deform_dual(polygons_struct *polygons1, polygons_struct *polygons2,
                      polygons_struct *polygons_orig, float *input, nifti_image *nii_ptr,
                      double w[3], double sigma, float lim1, float lim2,
                      double *target_distance, int it, int verbose)
{
    int i, j, k, v, dims[3], nvox, pidx, n_self_hits;
    int *n_neighbours, **neighbours;
    double vx[3], w2[3], scale_field;
    int have1 = (polygons1 != NULL);
    int have2 = (polygons2 != NULL);

    /* Need at least one surface */
    if (!have1 && !have2)
        return;

    /* We need other weightings for white surface */
    w2[0] = 0.3 * w[0];
    w2[1] = 0.3 * w[1];
    w2[2] = 3.0 * w[2];

    /* The "active" surface is the one used for topology (neighbours etc.) */
    polygons_struct *active = have1 ? polygons1 : polygons2;
    int n_points = active->n_points;

    // Create backup copies of original meshes
    object_struct *orig_object1 = NULL;
    object_struct *orig_object2 = NULL;
    polygons_struct *polygons1_orig = NULL;
    polygons_struct *polygons2_orig = NULL;

    if (have1)
    {
        orig_object1 = create_object(POLYGONS);
        polygons1_orig = get_polygons_ptr(orig_object1);
        copy_polygons(polygons_orig, polygons1_orig);
    }
    if (have2)
    {
        orig_object2 = create_object(POLYGONS);
        polygons2_orig = get_polygons_ptr(orig_object2);
        copy_polygons(polygons_orig, polygons2_orig);
    }

    // Extract image dimensions and voxel size
    dims[0] = nii_ptr->nx;
    dims[1] = nii_ptr->ny;
    dims[2] = nii_ptr->nz;
    nvox = dims[0] * dims[1] * dims[2];

    vx[0] = nii_ptr->dx;
    vx[1] = nii_ptr->dy;
    vx[2] = nii_ptr->dz;

    // Gradient volume
    float *gradient_x = (float *)malloc(sizeof(float) * nvox);
    float *gradient_y = (float *)malloc(sizeof(float) * nvox);
    float *gradient_z = (float *)malloc(sizeof(float) * nvox);
    int *flags = (int *)calloc(n_points, sizeof(int));

    if (!gradient_x || !gradient_y || !gradient_z || !flags)
    {
        fprintf(stderr, "Memory allocation error\n");
        exit(EXIT_FAILURE);
    }

    // Compute gradient of the input volume
    gradient3D(input, NULL, gradient_x, gradient_y, gradient_z, dims, vx);

    // Compute surface normals and neighbors
    if (have1)
        compute_polygon_normals(polygons1);
    if (have2)
        compute_polygon_normals(polygons2);

    // Neighbours from active surface (all share the same topology)
    create_polygon_point_neighbours(active, TRUE, &n_neighbours, &neighbours, NULL, NULL);

    // Allocate displacement fields
    double (*displacement_field1)[3] = NULL;
    double (*displacement_field2)[3] = NULL;

    if (have1)
    {
        displacement_field1 = malloc(sizeof(double[3]) * n_points);
        if (!displacement_field1)
        {
            fprintf(stderr, "Memory allocation error\n");
            exit(EXIT_FAILURE);
        }
    }
    if (have2)
    {
        displacement_field2 = malloc(sizeof(double[3]) * n_points);
        if (!displacement_field2)
        {
            fprintf(stderr, "Memory allocation error\n");
            exit(EXIT_FAILURE);
        }
    }

    // Curvature array computed once from the initial/reference mesh.
    double *curv = (double *)malloc(sizeof(double) * n_points);
    if (!curv)
    {
        fprintf(stderr, "Memory allocation error (curvature)\n");
        exit(EXIT_FAILURE);
    }
    get_polygon_vertex_curvatures_cg(polygons_orig, n_neighbours, neighbours, 3.0, 0, curv);

    // Iterative deformation process
    double s1_prev = FLT_MAX, s2_prev = FLT_MAX;
    int counter1 = 0, counter2 = 0;
    double curv_weight_less = -0.9; // tuneable: negative reduces smoothing

    for (i = 0; i < it; i++)
    {
        double s1 = 0.0, s2 = 0.0;

        // Process surfaces
        for (v = 0; v < n_points; v++)
        {
            // --- Surface 1 (pial) ---
            if (have1)
            {
                double c1[3] = {0.0, 0.0, 0.0};
                for (j = 0; j < n_neighbours[v]; j++)
                {
                    pidx = neighbours[v][j];
                    c1[0] += Point_x(polygons1->points[pidx]);
                    c1[1] += Point_y(polygons1->points[pidx]);
                    c1[2] += Point_z(polygons1->points[pidx]);
                }
                for (k = 0; k < 3; k++)
                    c1[k] /= n_neighbours[v];

                double p1[3] = {Point_x(polygons1->points[v]),
                                Point_y(polygons1->points[v]),
                                Point_z(polygons1->points[v])};
                double n1[3] = {Point_x(polygons1->normals[v]),
                                Point_y(polygons1->normals[v]),
                                Point_z(polygons1->normals[v])};

                float di1 = isoval(input, p1[0], p1[1], p1[2], dims, nii_ptr) - lim1;
                float fx1 = isoval(gradient_x, p1[0], p1[1], p1[2], dims, nii_ptr);
                float fy1 = isoval(gradient_y, p1[0], p1[1], p1[2], dims, nii_ptr);
                float fz1 = isoval(gradient_z, p1[0], p1[1], p1[2], dims, nii_ptr);
                float f2_1 = fmax(-1.0, fmin(1.0, fx1 * n1[0] + fy1 * n1[1] + fz1 * n1[2]));
                float boost1 = 0.5 + tanh(fabs(di1));
                float w3_scaled1 = fmin(w[2] * boost1, 5.0 * w[2]);

                if (counter1 == 0)
                {
                    for (k = 0; k < 3; k++)
                        displacement_field1[v][k] = 2.0 * w[0] * (c1[k] - p1[k]) +
                                                    ((w[1] * f2_1 + w3_scaled1) * di1) *
                                                        n1[k];
                }
                else
                {
                    for (k = 0; k < 3; k++)
                        displacement_field1[v][k] = 0.0;
                }
                s1 += di1 * di1;
            }

            // --- Surface 2 (white) ---
            if (have2)
            {
                double c2[3] = {0.0, 0.0, 0.0};
                for (j = 0; j < n_neighbours[v]; j++)
                {
                    pidx = neighbours[v][j];
                    c2[0] += Point_x(polygons2->points[pidx]);
                    c2[1] += Point_y(polygons2->points[pidx]);
                    c2[2] += Point_z(polygons2->points[pidx]);
                }
                for (k = 0; k < 3; k++)
                    c2[k] /= n_neighbours[v];

                double p2[3] = {Point_x(polygons2->points[v]),
                                Point_y(polygons2->points[v]),
                                Point_z(polygons2->points[v])};
                double n2[3] = {Point_x(polygons2->normals[v]),
                                Point_y(polygons2->normals[v]),
                                Point_z(polygons2->normals[v])};

                float di2 = isoval(input, p2[0], p2[1], p2[2], dims, nii_ptr) - lim2;
                float fx2 = isoval(gradient_x, p2[0], p2[1], p2[2], dims, nii_ptr);
                float fy2 = isoval(gradient_y, p2[0], p2[1], p2[2], dims, nii_ptr);
                float fz2 = isoval(gradient_z, p2[0], p2[1], p2[2], dims, nii_ptr);
                float f2_2 = fmax(-1.0, fmin(1.0, fx2 * n2[0] + fy2 * n2[1] + fz2 * n2[2]));
                float boost2 = 0.5 + tanh(fabs(di2));
                float w3_scaled2 = fmin(w2[2] * boost2, 5.0 * w2[2]);

                if (counter2 == 0)
                {
                    for (k = 0; k < 3; k++)
                        displacement_field2[v][k] = w2[0] * (c2[k] - p2[k]) +
                                                    ((w2[1] * f2_2 + w3_scaled2) * di2) *
                                                        n2[k];
                }
                else
                {
                    for (k = 0; k < 3; k++)
                        displacement_field2[v][k] = 0.0;
                }
                s2 += di2 * di2;
            }
        }

        if (have1 && s1 > s1_prev)
            counter1++;
        if (have2 && s2 > s2_prev)
            counter2++;
        // Stop only after 5 consecutive non-improving iterations on both surfaces
        if ((!have1 || counter1 > 1) && (!have2 || counter2 > 1))
            break;

        // Curvature-aware smoothing
        if (have1)
            smooth_displacement_field_blended(displacement_field1, polygons1,
                                              n_neighbours, neighbours, 5, sigma, 0.1, 10,
                                              curv, +1, curv_weight_less);
        if (have2)
            smooth_displacement_field_blended(displacement_field2, polygons2,
                                              n_neighbours, neighbours, 5, sigma, 0.1, 10,
                                              curv, -1, curv_weight_less);

        // Apply displacement
        if (have1)
        {
            for (v = 0; v < n_points; v++)
            {
                if ((i>0) && (flags[v]>0)) {
                    scale_field = 0.1*fmax(0.0, 10.0-(float)flags[v]);
                } else
                    scale_field = 1.0;
                Point_x(polygons1->points[v]) += scale_field*displacement_field1[v][0];
                Point_y(polygons1->points[v]) += scale_field*displacement_field1[v][1];
                Point_z(polygons1->points[v]) += scale_field*displacement_field1[v][2];
            }
        }
        if (have2)
        {
            for (v = 0; v < n_points; v++)
            {
                Point_x(polygons2->points[v]) += displacement_field2[v][0];
                Point_y(polygons2->points[v]) += displacement_field2[v][1];
                Point_z(polygons2->points[v]) += displacement_field2[v][2];
            }
        }

        // Minimize self-intersections
        if (have1)
        {
            n_self_hits = 0;
            int *flags1 = find_near_self_intersections(polygons1, 0.75, &n_self_hits);
            for (v = 0; v < n_points; v++)
            {
                if (flags1[v])
                {
                    flags[v]++;
                    Point_x(polygons1->points[v]) -= 1.5*displacement_field1[v][0];
                    Point_y(polygons1->points[v]) -= 1.5*displacement_field1[v][1];
                    Point_z(polygons1->points[v]) -= 1.5*displacement_field1[v][2];
                } else flags[v] = 0;
            }
            free(flags1);
        }

        if (have2)
        {
            n_self_hits = 0;
            int *flags2 = find_near_self_intersections(polygons2, 0.75, &n_self_hits);
            for (v = 0; v < n_points; v++)
            {
                if (flags2[v])
                {
                    flags[v]++;
                    Point_x(polygons2->points[v]) -= displacement_field2[v][0];
                    Point_y(polygons2->points[v]) -= displacement_field2[v][1];
                    Point_z(polygons2->points[v]) -= displacement_field2[v][2];
                } else if (!have1) flags[v] = 0;
            }
            free(flags2);
        }

        if (have1)
            compute_polygon_normals(polygons1);
        if (have2)
            compute_polygon_normals(polygons2);

        if (verbose)
        {
            if (have1 && have2)
                fprintf(stdout, "\rMesh: deform: iter %03d | Errors: %6.4f/%6.4f", i + 1,
                        sqrt(s1 / n_points), sqrt(s2 / n_points));
            else if (have1)
                fprintf(stdout, "\rMesh: deform: iter %03d | Error pial: %6.4f", i + 1,
                        sqrt(s1 / n_points));
            else
                fprintf(stdout, "\rMesh: deform: iter %03d | Error white: %6.4f", i + 1,
                        sqrt(s2 / n_points));
            fflush(stdout);
        }

        if (have1)
            s1_prev = s1;
        if (have2)
            s2_prev = s2;
    }
    if (verbose)
        fprintf(stdout, "\n");

    // Post-processing: compute total displacement, smooth, and re-apply
    if (have1)
    {
        for (v = 0; v < n_points; v++)
        {
            displacement_field1[v][0] = Point_x(polygons1->points[v]) -
                                        Point_x(polygons1_orig->points[v]);
            displacement_field1[v][1] = Point_y(polygons1->points[v]) -
                                        Point_y(polygons1_orig->points[v]);
            displacement_field1[v][2] = Point_z(polygons1->points[v]) -
                                        Point_z(polygons1_orig->points[v]);
        }
        smooth_displacement_field_blended(displacement_field1, polygons1, n_neighbours,
                                          neighbours, 5, 0.1 * sigma, 0.1, 10,
                                          curv, +1, curv_weight_less / 2.0);
        for (v = 0; v < n_points; v++)
        {
            Point_x(polygons1_orig->points[v]) += displacement_field1[v][0];
            Point_y(polygons1_orig->points[v]) += displacement_field1[v][1];
            Point_z(polygons1_orig->points[v]) += displacement_field1[v][2];
        }
        copy_polygons(polygons1_orig, polygons1);
    }

    if (have2)
    {
        for (v = 0; v < n_points; v++)
        {
            displacement_field2[v][0] = Point_x(polygons2->points[v]) -
                                        Point_x(polygons2_orig->points[v]);
            displacement_field2[v][1] = Point_y(polygons2->points[v]) -
                                        Point_y(polygons2_orig->points[v]);
            displacement_field2[v][2] = Point_z(polygons2->points[v]) -
                                        Point_z(polygons2_orig->points[v]);
        }
        smooth_displacement_field_blended(displacement_field2, polygons2, n_neighbours,
                                          neighbours, 5, 0.1 * sigma, 0.1, 10,
                                          curv, -1, curv_weight_less / 2.0);
        for (v = 0; v < n_points; v++)
        {
            Point_x(polygons2_orig->points[v]) += displacement_field2[v][0];
            Point_y(polygons2_orig->points[v]) += displacement_field2[v][1];
            Point_z(polygons2_orig->points[v]) += displacement_field2[v][2];
        }
        copy_polygons(polygons2_orig, polygons2);
    }

    if (verbose)
        fprintf(stdout, "\n");
    if (have1)
        remove_near_intersections(polygons1, 0.75, verbose);
    if (have2)
        remove_near_intersections(polygons2, 0.75, verbose);

    // Final Laplacian smoothing
    if (have1)
        smooth_laplacian(polygons1, 10, 0.1, 0.5);
    if (have2)
        smooth_laplacian(polygons2, 10, 0.1, 0.5);

    // Free allocated memory
    free(gradient_x);
    free(gradient_y);
    free(gradient_z);
    free(flags);
    if (displacement_field1)
        free(displacement_field1);
    if (displacement_field2)
        free(displacement_field2);
    free(curv);
    delete_polygon_point_neighbours(active, n_neighbours, neighbours, NULL, NULL);
}

/**
 * \brief Refine pial and white surfaces toward image intensity edges using
 *        normal-ray edge search (FreeSurfer-inspired approach).
 *
 * After surf_deform_dual has positioned the surfaces near their targets,
 * this function performs additional refinement by searching along each
 * vertex normal for the steepest intensity gradient (the actual tissue
 * boundary), then nudging vertices toward that edge location.
 *
 * Unlike the balloon + gradient force approach of surf_deform_dual which
 * tends to over-smooth, this function directly locates edges by profiling
 * intensity along the normal ray, similar to FreeSurfer's
 * mrisComputeTargetLocationTerm.
 *
 * Each iteration:
 *  1. For each vertex, sample intensity at multiple offsets along the normal
 *  2. Compute finite-difference gradient magnitude at each sample point
 *  3. Select the offset with maximum gradient magnitude as the edge location
 *  4. Verify the edge is consistent with the expected intensity threshold
 *  5. Compute displacement toward the edge with tangential smoothing
 *  6. Apply a thickness constraint to prevent pial-white collapse
 *  7. Check and revert self-intersecting vertices
 *
 * \param polygons1        (in/out) pial surface mesh
 * \param polygons2        (in/out) white surface mesh
 * \param input            (in)     intensity volume data
 * \param nii_ptr          (in)     NIfTI header for dimensions and transforms
 * \param lim1             (in)     intensity threshold for pial surface
 * \param lim2             (in)     intensity threshold for white surface
 * \param target_distance  (in)     desired pial-white spacing per vertex
 * \param it               (in)     number of refinement iterations
 * \param verbose          (in)     if nonzero, print progress
 */
void surf_deform_gradient_dual(polygons_struct *polygons1, polygons_struct *polygons2,
                               float *input, nifti_image *nii_ptr,
                               float lim1, float lim2,
                               double *target_distance, int it, int verbose)
{
    int i, j, k, v, dims[3], pidx, n_self_hits;
    int *n_neighbours, **neighbours;

    int n_points = polygons1->n_points;

    dims[0] = nii_ptr->nx;
    dims[1] = nii_ptr->ny;
    dims[2] = nii_ptr->nz;

    /* Setup neighbors (shared topology) */
    compute_polygon_normals(polygons1);
    compute_polygon_normals(polygons2);
    create_polygon_point_neighbours(polygons1, TRUE, &n_neighbours, &neighbours, NULL, NULL);

    /* Displacement fields */
    double (*disp1)[3] = malloc(sizeof(double[3]) * n_points);
    double (*disp2)[3] = malloc(sizeof(double[3]) * n_points);
    if (!disp1 || !disp2)
    {
        fprintf(stderr, "Memory allocation error\n");
        exit(EXIT_FAILURE);
    }

    /* Normal-ray search parameters */
    const double search_dist = 1.0; /* search ±1mm along normal */
    const int n_samples = 11;       /* sample points along ray */
    const double sample_step = 2.0 * search_dist / (n_samples - 1);
    const double step_frac = 0.5;          /* fraction of edge distance to apply */
    const double max_disp = 0.5;           /* max displacement per iteration (mm) */
    const double smooth_weight = 0.05;     /* tangential smoothing weight */
    const double min_thickness_frac = 0.5; /* minimum thickness fraction */
    const double min_grad_mag = 0.05;      /* min gradient to accept as edge */

    double err1_prev = FLT_MAX, err2_prev = FLT_MAX;
    int counter1 = 0, counter2 = 0;

    for (i = 0; i < it; i++)
    {
        double err1 = 0.0, err2 = 0.0;

        for (v = 0; v < n_points; v++)
        {
            /* ---- PIAL SURFACE (polygons1) ---- */
            double p1[3] = {Point_x(polygons1->points[v]),
                            Point_y(polygons1->points[v]),
                            Point_z(polygons1->points[v])};
            double n1[3] = {Point_x(polygons1->normals[v]),
                            Point_y(polygons1->normals[v]),
                            Point_z(polygons1->normals[v])};

            double len = sqrt(n1[0] * n1[0] + n1[1] * n1[1] + n1[2] * n1[2]);
            if (len > 1e-10)
            {
                n1[0] /= len;
                n1[1] /= len;
                n1[2] /= len;
            }

            /* Sample intensity along normal ray and find closest edge with
             * correct transition direction.
             * For both pial and white surfaces the outward normal points away
             * from the brain, so at the true tissue boundary intensity must
             * decrease along the outward normal (GM->CSF or WM->GM).
             * We therefore only accept sample points where the signed
             * directional derivative dI/dt < 0 (negative gradient along
             * outward normal).  Among all valid edges we pick the one
             * closest to the current vertex position (smallest |t|) to
             * avoid locking onto the opposite sulcal wall. */
            double best_t1 = 0.0, best_grad1 = 0.0;
            {
                double intensities[11]; /* n_samples */
                int si;
                for (si = 0; si < n_samples; si++)
                {
                    double t = -search_dist + si * sample_step;
                    intensities[si] = isoval(input,
                                             p1[0] + t * n1[0],
                                             p1[1] + t * n1[1],
                                             p1[2] + t * n1[2], dims, nii_ptr);
                }
                /* Finite-difference signed gradient along ray (dI/dt) */
                double closest_abs_t = search_dist + 1.0; /* > any valid |t| */
                for (si = 1; si < n_samples - 1; si++)
                {
                    double grad_signed = (intensities[si + 1] - intensities[si - 1]) /
                                         (2.0 * sample_step);
                    double grad_abs = fabs(grad_signed);
                    double t = -search_dist + si * sample_step;
                    /* Accept only edges where intensity decreases outward
                     * (grad_signed < 0) and gradient is strong enough */
                    if (grad_signed < 0.0 && grad_abs > min_grad_mag &&
                        fabs(t) < closest_abs_t)
                    {
                        closest_abs_t = fabs(t);
                        best_grad1 = grad_abs;
                        best_t1 = t;
                    }
                }
            }

            /* Compute displacement toward edge */
            double d1 = 0.0;
            if (counter1 == 0 && best_grad1 > min_grad_mag)
                d1 = fmax(-max_disp, fmin(max_disp, step_frac * best_t1));

            /* Tangential smoothing: centroid attraction onto tangent plane */
            double c1[3] = {0.0, 0.0, 0.0};
            for (j = 0; j < n_neighbours[v]; j++)
            {
                pidx = neighbours[v][j];
                c1[0] += Point_x(polygons1->points[pidx]);
                c1[1] += Point_y(polygons1->points[pidx]);
                c1[2] += Point_z(polygons1->points[pidx]);
            }
            for (k = 0; k < 3; k++)
                c1[k] /= n_neighbours[v];

            double tc1[3] = {c1[0] - p1[0], c1[1] - p1[1], c1[2] - p1[2]};
            double dot1 = tc1[0] * n1[0] + tc1[1] * n1[1] + tc1[2] * n1[2];
            for (k = 0; k < 3; k++)
                tc1[k] -= dot1 * n1[k];

            for (k = 0; k < 3; k++)
                disp1[v][k] = d1 * n1[k] + smooth_weight * tc1[k];

            double di1 = isoval(input, p1[0], p1[1], p1[2], dims, nii_ptr) - lim1;
            err1 += di1 * di1;

            /* ---- WHITE SURFACE (polygons2) ---- */
            double p2[3] = {Point_x(polygons2->points[v]),
                            Point_y(polygons2->points[v]),
                            Point_z(polygons2->points[v])};
            double n2[3] = {Point_x(polygons2->normals[v]),
                            Point_y(polygons2->normals[v]),
                            Point_z(polygons2->normals[v])};

            len = sqrt(n2[0] * n2[0] + n2[1] * n2[1] + n2[2] * n2[2]);
            if (len > 1e-10)
            {
                n2[0] /= len;
                n2[1] /= len;
                n2[2] /= len;
            }

            double best_t2 = 0.0, best_grad2 = 0.0;
            {
                double intensities[11];
                int si;
                for (si = 0; si < n_samples; si++)
                {
                    double t = -search_dist + si * sample_step;
                    intensities[si] = isoval(input,
                                             p2[0] + t * n2[0],
                                             p2[1] + t * n2[1],
                                             p2[2] + t * n2[2], dims, nii_ptr);
                }
                /* Same direction-aware search: closest edge with negative
                 * directional derivative (intensity decreasing outward) */
                double closest_abs_t = search_dist + 1.0;
                for (si = 1; si < n_samples - 1; si++)
                {
                    double grad_signed = (intensities[si + 1] - intensities[si - 1]) /
                                         (2.0 * sample_step);
                    double grad_abs = fabs(grad_signed);
                    double t = -search_dist + si * sample_step;
                    if (grad_signed < 0.0 && grad_abs > min_grad_mag &&
                        fabs(t) < closest_abs_t)
                    {
                        closest_abs_t = fabs(t);
                        best_grad2 = grad_abs;
                        best_t2 = t;
                    }
                }
            }

            double d2 = 0.0;
            if (counter2 == 0 && best_grad2 > min_grad_mag)
                d2 = fmax(-max_disp, fmin(max_disp, step_frac * best_t2));

            double c2[3] = {0.0, 0.0, 0.0};
            for (j = 0; j < n_neighbours[v]; j++)
            {
                pidx = neighbours[v][j];
                c2[0] += Point_x(polygons2->points[pidx]);
                c2[1] += Point_y(polygons2->points[pidx]);
                c2[2] += Point_z(polygons2->points[pidx]);
            }
            for (k = 0; k < 3; k++)
                c2[k] /= n_neighbours[v];

            double tc2[3] = {c2[0] - p2[0], c2[1] - p2[1], c2[2] - p2[2]};
            double dot2 = tc2[0] * n2[0] + tc2[1] * n2[1] + tc2[2] * n2[2];
            for (k = 0; k < 3; k++)
                tc2[k] -= dot2 * n2[k];

            for (k = 0; k < 3; k++)
                disp2[v][k] = d2 * n2[k] + smooth_weight * tc2[k];

            double di2 = isoval(input, p2[0], p2[1], p2[2], dims, nii_ptr) - lim2;
            err2 += di2 * di2;

            /* ---- THICKNESS CONSTRAINT ---- */
            double p1_new[3], p2_new[3];
            for (k = 0; k < 3; k++)
            {
                p1_new[k] = p1[k] + disp1[v][k];
                p2_new[k] = p2[k] + disp2[v][k];
            }
            double ddx = p1_new[0] - p2_new[0];
            double ddy = p1_new[1] - p2_new[1];
            double ddz = p1_new[2] - p2_new[2];
            double new_dist = sqrt(ddx * ddx + ddy * ddy + ddz * ddz);

            if (new_dist < min_thickness_frac * target_distance[v])
            {
                for (k = 0; k < 3; k++)
                {
                    disp1[v][k] *= 0.5;
                    disp2[v][k] *= 0.5;
                }
            }
        }

        /* Early termination if both surfaces stopped improving */
        if (err1 > err1_prev)
            counter1++;
        if (err2 > err2_prev)
            counter2++;
        if (counter1 > 3 && counter2 > 3)
            break;

        /* Apply displacements directly — no global displacement smoothing
         * since tangential smoothing is already built into each vertex's
         * displacement and the edge-seeking signal should not be blurred. */
        for (v = 0; v < n_points; v++)
        {
            Point_x(polygons1->points[v]) += disp1[v][0];
            Point_y(polygons1->points[v]) += disp1[v][1];
            Point_z(polygons1->points[v]) += disp1[v][2];
        }
        for (v = 0; v < n_points; v++)
        {
            Point_x(polygons2->points[v]) += disp2[v][0];
            Point_y(polygons2->points[v]) += disp2[v][1];
            Point_z(polygons2->points[v]) += disp2[v][2];
        }

        /* Self-intersection check and revert */
        n_self_hits = 0;
        int *flags1 = find_near_self_intersections(polygons1, 0.75, &n_self_hits);
        for (v = 0; v < n_points; v++)
        {
            if (flags1[v])
            {
                Point_x(polygons1->points[v]) -= disp1[v][0];
                Point_y(polygons1->points[v]) -= disp1[v][1];
                Point_z(polygons1->points[v]) -= disp1[v][2];
            }
        }
        free(flags1);

        int *flags2 = find_near_self_intersections(polygons2, 0.75, &n_self_hits);
        for (v = 0; v < n_points; v++)
        {
            if (flags2[v])
            {
                Point_x(polygons2->points[v]) -= disp2[v][0];
                Point_y(polygons2->points[v]) -= disp2[v][1];
                Point_z(polygons2->points[v]) -= disp2[v][2];
            }
        }
        free(flags2);

        /* Update normals for next iteration */
        compute_polygon_normals(polygons1);
        compute_polygon_normals(polygons2);

        if (verbose)
        {
            fprintf(stdout, "\rMesh: gradient refine: iter %03d | Errors: %6.4f/%6.4f",
                    i + 1, sqrt(err1 / n_points), sqrt(err2 / n_points));
            fflush(stdout);
        }

        err1_prev = err1;
        err2_prev = err2;
    }
    if (verbose)
        fprintf(stdout, "\n");

    /* Final light intersection removal and Laplacian smoothing */
    remove_near_intersections(polygons1, 0.75, verbose);
    remove_near_intersections(polygons2, 0.75, verbose);
    smooth_laplacian(polygons1, 5, 0.1, 0.5);
    smooth_laplacian(polygons2, 5, 0.1, 0.5);

    /* Cleanup */
    free(disp1);
    free(disp2);
    delete_polygon_point_neighbours(polygons1, n_neighbours, neighbours, NULL, NULL);
}
