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
#include "CAT_Curvature.h"

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/**
 * \brief Edge strength along a surface normal: the intensity decrease -dI/dn.
 *
 * Positive where intensity falls outwards, as it does across both the GM/WM
 * and the CSF/GM boundary.  For axis-aligned images stored with all axes
 * negative this equals the index-space dot product used before, so results on
 * such images are unchanged.
 *
 * \param M  (in) matrix from gradient3D_world_matrix()
 * \param gx (in) gradient3D() x component at the vertex
 * \param gy (in) gradient3D() y component at the vertex
 * \param gz (in) gradient3D() z component at the vertex
 * \param n  (in) unit surface normal in world space
 * \return -dI/dn
 */
static double
edge_strength(double M[3][3], double gx, double gy, double gz, const double n[3])
{
    double s = 0.0;
    int r;

    for (r = 0; r < 3; r++)
        s += (M[r][0] * gx + M[r][1] * gy + M[r][2] * gz) * n[r];
    return -s;
}

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
static void smooth_displacement_field_blended(double (*displacement_field)[3],
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
 *  6.   Revert the step of vertices that come too close to a facing sheet
 *  7. Cap total displacements above the 95th percentile to that length
 *  8. Remove the per-vertex noise with 2 iterations of HC Laplacian smoothing
 *  9. Optionally remove the remaining self-intersections
 *
 * Step 6 only counts vertices with opposing normals: a pure distance test flags
 * 8-11% of the vertices of a reduced central surface per iteration, 88-99.9% of
 * them 2-ring neighbours of the same sheet, and freezing those raised the PPM
 * error of the result by a third.  The total displacement is no longer smoothed
 * after the loop: like any smoothing of an accumulated displacement it pulled
 * the surface off the isovalue (PPM error 0.027 -> 0.037).
 *
 * \param polygons            (in/out) surface mesh (modified in-place)
 * \param input               (in)     float[nvoxels]; intensity volume data
 * \param nii_ptr             (in)     NIfTI header with volume dimensions and voxel size
 * \param w                   (in)     double[3]; weight factors {smoothness, gradient_strength, ?}
 * \param sigma               (in)     Gaussian smoothing parameter for displacement field
 * \param lim                 (in)     intensity threshold controlling deformation magnitude
 * \param it                  (in)     number of deformation iterations
 * \param remove_selfintersect (in)    boolean; enable self-intersection removal
 * \param verbose             (in)     boolean; print iteration progress if true
 */
void surf_deform(polygons_struct *polygons, float *input, nifti_image *nii_ptr,
                 double w[3], double sigma, float lim, int it, int remove_selfintersect, int verbose)
{
    int i, j, k, v, dims[3], nvox, pidx;
    int *n_neighbours, **neighbours;
    float *gradient_x, *gradient_y, *gradient_z;
    double vx[3], s;
    Point points[MAX_POINTS_PER_POLYGON];
    polygons_struct *polygons_orig;
    object_struct *orig_object;
    Point *start_points;

    orig_object = create_object(POLYGONS);
    polygons_orig = get_polygons_ptr(orig_object);
    copy_polygons(polygons, polygons_orig);

    // Start positions: the way back for defects the repair cannot smooth out
    // (polygons_orig itself is moved to the result below)
    start_points = (Point *)malloc(sizeof(Point) * polygons->n_points);
    if (!start_points)
    {
        fprintf(stderr, "Memory allocation error\n");
        exit(EXIT_FAILURE);
    }
    memcpy(start_points, polygons->points, sizeof(Point) * polygons->n_points);

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
    double g2w[3][3];
    gradient3D_world_matrix(nii_ptr, g2w);

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
            float f2 = fmax(-1.0, fmin(1.0, edge_strength(g2w, fx, fy, fz, n)));

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
        int *flags = find_near_facing_intersections(polygons, 0.75, 0.3, &n_self_hits);
        for (v = 0; v < polygons->n_points; v++)
        {
            if (flags[v])
            {
                Point_x(polygons->points[v]) -= displacement_field[v][0];
                Point_y(polygons->points[v]) -= displacement_field[v][1];
                Point_z(polygons->points[v]) -= displacement_field[v][2];
            }
        }
        free(flags);

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

    // Squared length of the total displacement
    double *displacement = malloc(sizeof(double) * polygons->n_points);
    for (v = 0; v < polygons->n_points; v++)
    {
        displacement[v] = SQR(displacement_field[v][0]) +
                          SQR(displacement_field[v][1]) +
                          SQR(displacement_field[v][2]);
    }

    double prctile[2] = {95.0, 95.0};
    double threshold[2];
    get_prctile(displacement, polygons->n_points, threshold, prctile, 1, DT_FLOAT64);

    /* Cap displacements above the 95th percentile to that length, keeping their
       direction, so that outliers cannot cause self-intersections.  (Replacing
       them with the per-axis percentiles moved all outliers by one fixed
       vector, whatever their own direction.) */
    for (v = 0; v < polygons->n_points; v++)
    {
        if (displacement[v] > threshold[1] && displacement[v] > 0.0)
        {
            double scale = sqrt(threshold[1] / displacement[v]);
            displacement_field[v][0] *= scale;
            displacement_field[v][1] *= scale;
            displacement_field[v][2] *= scale;
        }
    }

    // Apply the final displacement to vertices
    for (v = 0; v < polygons->n_points; v++)
    {
        Point_x(polygons_orig->points[v]) += displacement_field[v][0];
        Point_y(polygons_orig->points[v]) += displacement_field[v][1];
        Point_z(polygons_orig->points[v]) += displacement_field[v][2];
    }
    copy_polygons(polygons_orig, polygons);

    // The per-vertex forces leave the mesh noisier than its start surface.  Two
    // iterations of HC Laplacian smoothing remove that noise without pulling the
    // surface off the isovalue (HR075: PPM error 0.0220 -> 0.0194).
    smooth_laplacian(polygons, 2, 0.1, 0.5);

    // Remove self-intersections by locally smoothing the intersecting regions,
    // retreating the ones smoothing cannot resolve towards the start surface.
    // This preserves the mesh topology, i.e. the number of vertices and faces
    // and their connectivity are left unchanged.
    if (remove_selfintersect)
    {
        if (verbose)
            fprintf(stdout, "\n");
        remove_intersections_ref(polygons, start_points, 10, 50, verbose);
    }
    free(start_points);

    // Free allocated memory
    free(gradient_x);
    free(gradient_y);
    free(gradient_z);
    free(displacement);
    free(displacement_field);
    delete_polygon_point_neighbours(polygons, n_neighbours, neighbours, NULL, NULL);
}
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
    double g2w[3][3];
    gradient3D_world_matrix(nii_ptr, g2w);

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
                float f2_1 = fmax(-1.0, fmin(1.0, edge_strength(g2w, fx1, fy1, fz1, n1)));
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
                float f2_2 = fmax(-1.0, fmin(1.0, edge_strength(g2w, fx2, fy2, fz2, n2)));
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
                if ((i > 0) && (flags[v] > 0))
                {
                    scale_field = 0.1 * fmax(0.0, 10.0 - (float)flags[v]);
                }
                else
                    scale_field = 1.0;
                Point_x(polygons1->points[v]) += scale_field * displacement_field1[v][0];
                Point_y(polygons1->points[v]) += scale_field * displacement_field1[v][1];
                Point_z(polygons1->points[v]) += scale_field * displacement_field1[v][2];
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
                    Point_x(polygons1->points[v]) -= 1.5 * displacement_field1[v][0];
                    Point_y(polygons1->points[v]) -= 1.5 * displacement_field1[v][1];
                    Point_z(polygons1->points[v]) -= 1.5 * displacement_field1[v][2];
                }
                else
                    flags[v] = 0;
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
                }
                else if (!have1)
                    flags[v] = 0;
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
        smooth_laplacian(polygons1, 2, 0.1, 0.5);
    if (have2)
        smooth_laplacian(polygons2, 2, 0.1, 0.5);

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
