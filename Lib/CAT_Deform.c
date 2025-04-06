/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry, University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 *
 */


#include "CAT_Surf.h"
#include "CAT_NiftiLib.h"
#include "CAT_Intersect.h"

/**
 * @brief Smooths the displacement field to avoid abrupt vertex movements.
 *
 * The displacement field is iteratively updated to ensure smooth transitions between neighboring vertices.
 *
 * @param displacement_field 3D displacement field for each vertex.
 * @param polygons The 3D surface mesh.
 * @param n_neighbours Array storing the number of neighbors for each vertex.
 * @param neighbours 2D array storing the indices of each vertex’s neighbors.
 * @param iterations Number of smoothing iterations.
 * @param sigma Smoothing factor.
 */
void smooth_displacement_field(double (*displacement_field)[3], polygons_struct *polygons, 
                               int *n_neighbours, int **neighbours, int iterations, double sigma) 
{
    int v, j, k, it, pidx;
    double (*new_displacement)[3] = malloc(sizeof(double[3]) * polygons->n_points);

    for (it = 0; it < iterations; it++) {
        for (v = 0; v < polygons->n_points; v++) {
            double smoothed[3] = {0.0, 0.0, 0.0};
            int count = 0;

            for (j = 0; j < n_neighbours[v]; j++) {
                pidx = neighbours[v][j];
                for (k = 0; k < 3; k++)
                    smoothed[k] += displacement_field[pidx][k];
                count++;
            }

            for (k = 0; k < 3; k++)
                new_displacement[v][k] = (smoothed[k] / count) * exp(-sigma);
        }

        // Update displacement field
        for (v = 0; v < polygons->n_points; v++) {
            for (k = 0; k < 3; k++) {
                displacement_field[v][k] = new_displacement[v][k];
            }
        }
    }

    free(new_displacement);
}

/**
 * @brief Deforms a 3D surface mesh using an external force field derived from an input volume.
   This approach deforms a 3D surface mesh by computing forces from a reference image.
   It moves each vertex based on a balance of:
    - Internal forces → Maintain smoothness.
    - External forces → Derived from image intensity gradients.
    - Self-intersection prevention → Ensures a valid mesh.
   This method is related to Active Contour Models (Snakes) and Level Set Methods, 
   but explicitly operates on a mesh representation.
 
  This function iteratively moves each vertex of the input `polygons` mesh according to:
  1. **Internal force**: Keeps the mesh smooth by averaging neighboring vertices.
  2. **External force**: Derived from the intensity and gradient of an image.
  3. **Constraints**: Prevents self-intersections using a spatial grid.
 
  @param polygons Pointer to the surface mesh (polygons_struct).
  @param input Pointer to the intensity volume (3D image).
  @param nii_ptr Pointer to the NIfTI image structure (contains volume metadata).
  @param w Weight factors for internal smoothing and external forces (array of 3 floats).
  @param lim Intensity threshold that controls the deformation limit.
  @param it Number of deformation iterations.
 */
void surf_deform(polygons_struct *polygons, float *input, nifti_image *nii_ptr, 
                 double w[3], double sigma, float lim, int it, int del_intersections, int verbose)
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

    if (!gradient_x || !gradient_y || !gradient_z) {
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
    if (!displacement_field) {
        fprintf(stderr, "Memory allocation error\n");
        exit(EXIT_FAILURE);
    }

    // Iterative deformation process
    int counter = 0;
    double s_prev = FLT_MAX;
    for (i = 0; i < it; i++) {
        s = 0.0;

        for (v = 0; v < polygons->n_points; v++) {
            // Compute centroid of neighboring vertices for smoothing
            double c[3] = {0.0, 0.0, 0.0};
            for (j = 0; j < n_neighbours[v]; j++) {
                pidx = neighbours[v][j];
                c[0] += Point_x(polygons->points[pidx]);
                c[1] += Point_y(polygons->points[pidx]);
                c[2] += Point_z(polygons->points[pidx]);
            }
            for (k = 0; k < 3; k++) {
                c[k] /= (float)n_neighbours[v];
            }

            // Get current vertex position
            double p[3] = { Point_x(polygons->points[v]), 
                            Point_y(polygons->points[v]), 
                            Point_z(polygons->points[v]) };

            // Get vertex normal
            double n[3] = { Point_x(polygons->normals[v]), 
                            Point_y(polygons->normals[v]), 
                            Point_z(polygons->normals[v]) };

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
            for (k = 0; k < 3; k++) {
                displacement_field[v][k] = w[0] * (c[k] - p[k]) + 
                                          ((w[1] * f2 + w3_scaled) * f3) * n[k];
            }

            s += di * di;
        }

        // Stop if no further minimization is seen and s increases instead
        if (s > s_prev) counter++;
        if (counter > 10) break;

        // Apply smoothing to the displacement field
        smooth_displacement_field(displacement_field, polygons, n_neighbours, neighbours, 5, sigma);

        // Apply the final displacement to vertices
        for (v = 0; v < polygons->n_points; v++) {
            Point_x(polygons->points[v]) += displacement_field[v][0];
            Point_y(polygons->points[v]) += displacement_field[v][1];
            Point_z(polygons->points[v]) += displacement_field[v][2];
        }

        int n_self_hits = 0;
        int *flags = find_near_self_intersections(polygons, 0.75, &n_self_hits);
        for (v = 0; v < polygons->n_points; v++) {
            if (flags[v]) {
                Point_x(polygons->points[v]) -= displacement_field[v][0];
                Point_y(polygons->points[v]) -= displacement_field[v][1];
                Point_z(polygons->points[v]) -= displacement_field[v][2];
            }
        }

        // Update normals for next iteration
        compute_polygon_normals(polygons);
        if (verbose) {
            fprintf(stdout, "\rMesh: deform: iter %03d | Error: %6.4f", i+1, 
                            sqrt(s / polygons->n_points));
            fflush(stdout);  // Force output update
        }
        s_prev = s;
    }
    if (verbose) fprintf(stdout, "\n");

    for (v = 0; v < polygons->n_points; v++) {
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

    for (v = 0; v < polygons->n_points; v++) {
        displacement[v] = SQR(displacement_field[v][0]) + 
                          SQR(displacement_field[v][1]) + 
                          SQR(displacement_field[v][2]);
    }
    
    /* Get percentile for squared sum of displacements */
    double threshold[2];
    get_prctile(displacement, polygons->n_points, threshold, prctile, 1, DT_FLOAT64);

    /* If squaed sum of displacements is exceeding 95% percentile we limit all
      displacements to 95% percentile to prevent that these outliers cause 
      self-intersections */
    for (v = 0; v < polygons->n_points; v++) {
        if (displacement[v] > threshold[1]) {
            displacement_field[v][0] = threshold_x[1];
            displacement_field[v][1] = threshold_y[1];
            displacement_field[v][2] = threshold_z[1];
        }
    }

    // Apply very slight smoothing to the displacement field to smooth replacements
    smooth_displacement_field(displacement_field, polygons, n_neighbours, neighbours, 5, 0.05*sigma);

    // Apply the final displacement to vertices
    for (v = 0; v < polygons->n_points; v++) {
        Point_x(polygons_orig->points[v]) += displacement_field[v][0];
        Point_y(polygons_orig->points[v]) += displacement_field[v][1];
        Point_z(polygons_orig->points[v]) += displacement_field[v][2];
    }
    copy_polygons(polygons_orig, polygons);

    if (del_intersections) {
        if (verbose) fprintf(stdout,"\n");
        remove_near_intersections(polygons, 0.75, verbose);
    }

    // Free allocated memory
    free(gradient_x);
    free(gradient_y);
    free(gradient_z);
    free(displacement);
    free(displacement_field);
    delete_polygon_point_neighbours(polygons, n_neighbours, neighbours, NULL, NULL);
}

void surf_deform_dual(polygons_struct *polygons1, polygons_struct *polygons2, 
                      float *input, nifti_image *nii_ptr, 
                      double w[4], double sigma, float lim1, float lim2, 
                      double *target_distance, int it, int verbose) 
{
    int i, j, k, v, dims[3], nvox, pidx, n_self_hits;
    int *n_neighbours, **neighbours;
    float *gradient_x, *gradient_y, *gradient_z;
    double vx[3];
    
    // Create backup copies of original meshes
    polygons_struct *polygons1_orig, *polygons2_orig;
    object_struct *orig_object1 = create_object(POLYGONS);
    object_struct *orig_object2 = create_object(POLYGONS);
    polygons1_orig = get_polygons_ptr(orig_object1);
    polygons2_orig = get_polygons_ptr(orig_object2);
    copy_polygons(polygons1, polygons1_orig);
    copy_polygons(polygons2, polygons2_orig);

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

    if (!gradient_x || !gradient_y || !gradient_z) {
        fprintf(stderr, "Memory allocation error\n");
        exit(EXIT_FAILURE);
    }

    // Compute gradient of the input volume
    gradient3D(input, NULL, gradient_x, gradient_y, gradient_z, dims, vx);

    // Compute surface normals and neighbors
    compute_polygon_normals(polygons1);
    compute_polygon_normals(polygons2);

    // Both surfaces have the same topology, thus neighbours can be estimated from 1st
    create_polygon_point_neighbours(polygons1, TRUE, &n_neighbours, &neighbours, NULL, NULL);

    // Allocate displacement fields
    double (*displacement_field1)[3] = malloc(sizeof(double[3]) * polygons1->n_points);
    double (*displacement_field2)[3] = malloc(sizeof(double[3]) * polygons2->n_points);

    if (!displacement_field1 || !displacement_field2) {
        fprintf(stderr, "Memory allocation error\n");
        exit(EXIT_FAILURE);
    }

    // Iterative deformation process
    int counter1 = 0, counter2 = 0;
    double s1_prev = FLT_MAX, s2_prev = FLT_MAX;
    for (i = 0; i < it; i++) {
        double s1 = 0.0, s2 = 0.0;

        // Process both surfaces
        for (v = 0; v < polygons1->n_points; v++) {
            // Compute centroids for smoothing
            double c1[3] = {0.0, 0.0, 0.0};
            double c2[3] = {0.0, 0.0, 0.0};

            for (j = 0; j < n_neighbours[v]; j++) {
                pidx = neighbours[v][j];
                c1[0] += Point_x(polygons1->points[pidx]);
                c1[1] += Point_y(polygons1->points[pidx]);
                c1[2] += Point_z(polygons1->points[pidx]);
                c2[0] += Point_x(polygons2->points[pidx]);
                c2[1] += Point_y(polygons2->points[pidx]);
                c2[2] += Point_z(polygons2->points[pidx]);
            }
            for (k = 0; k < 3; k++) {
                c1[k] /= (float)n_neighbours[v];
                c2[k] /= (float)n_neighbours[v];
            }

            // Get current positions
            double p1[3] = { Point_x(polygons1->points[v]), 
                             Point_y(polygons1->points[v]), 
                             Point_z(polygons1->points[v]) };
            double p2[3] = { Point_x(polygons2->points[v]), 
                             Point_y(polygons2->points[v]), 
                             Point_z(polygons2->points[v]) };

            // Get normals
            double n1[3] = { Point_x(polygons1->normals[v]), 
                             Point_y(polygons1->normals[v]), 
                             Point_z(polygons1->normals[v]) };
            double n2[3] = { Point_x(polygons2->normals[v]), 
                             Point_y(polygons2->normals[v]), 
                             Point_z(polygons2->normals[v]) };

            // Compute external forces
            float di1 = isoval(input, p1[0], p1[1], p1[2], dims, nii_ptr) - lim1;
            float di2 = isoval(input, p2[0], p2[1], p2[2], dims, nii_ptr) - lim2;

            float fx1 = isoval(gradient_x, p1[0], p1[1], p1[2], dims, nii_ptr);
            float fy1 = isoval(gradient_y, p1[0], p1[1], p1[2], dims, nii_ptr);
            float fz1 = isoval(gradient_z, p1[0], p1[1], p1[2], dims, nii_ptr);

            float fx2 = isoval(gradient_x, p2[0], p2[1], p2[2], dims, nii_ptr);
            float fy2 = isoval(gradient_y, p2[0], p2[1], p2[2], dims, nii_ptr);
            float fz2 = isoval(gradient_z, p2[0], p2[1], p2[2], dims, nii_ptr);

            float f2_1 = fmax(-1.0, fmin(1.0, fx1 * n1[0] + fy1 * n1[1] + fz1 * n1[2]));
            float f2_2 = fmax(-1.0, fmin(1.0, fx2 * n2[0] + fy2 * n2[1] + fz2 * n2[2]));

            float f3_1 = di1;
            float f3_2 = di2;

            // Dynamic boosting (optional: limit max)
            float boost1 = 0.5 + tanh(fabs(di1));
            float boost2 = 0.5 + tanh(fabs(di2));
            float w3_scaled1 = fmin(w[2] * boost1, 5.0 * w[2]);
            float w3_scaled2 = fmin(w[2] * boost2, 5.0 * w[2]);

            // Compute Euclidean distance between corresponding points
            double current_distance = sqrt(
                (p2[0] - p1[0]) * (p2[0] - p1[0]) +
                (p2[1] - p1[1]) * (p2[1] - p1[1]) +
                (p2[2] - p1[2]) * (p2[2] - p1[2])
            );

            // Compute distance constraint force
            double distance_force[3];
            double distance_factor = (current_distance - target_distance[v]) / target_distance[v];
            for (k = 0; k < 3; k++) {
                if ((counter1 > 4) || (counter2 > 4))
                    distance_force[k] = 0.0;
                else
                    distance_force[k] = (p2[k] - p1[k]) * distance_factor;
            }

            // Compute final displacement
            for (k = 0; k < 3; k++) {
                displacement_field1[v][k] = w[0] * (c1[k] - p1[k]) + 
                                          ((w[1] * f2_1 + w3_scaled1) * f3_1) * 
                                           n1[k] + w[3] * distance_force[k];
                displacement_field2[v][k] = w[0] * (c2[k] - p2[k]) + 
                                          ((w[1] * f2_2 + w3_scaled2) * f3_2) * 
                                           n2[k] + w[3] * distance_force[k];
            }

            s1 += di1 * di1;
            s2 += di2 * di2;
        }

        if (s1 > s1_prev) counter1++;
        if (s2 > s2_prev) counter2++;

        // Stop if no further minimization is seen and s1 and s2 increase instead
        //if ((counter1 > 5) && (counter2 > 5)) break;

        // Apply smoothing to the (pial) displacement field
        smooth_displacement_field(displacement_field1, polygons1, n_neighbours, neighbours, 5, sigma);

        // Apply the final displacement to vertices
        for (v = 0; v < polygons1->n_points; v++) {
            Point_x(polygons1->points[v]) += displacement_field1[v][0];
            Point_y(polygons1->points[v]) += displacement_field1[v][1];
            Point_z(polygons1->points[v]) += displacement_field1[v][2];
        }
        
        n_self_hits = 0;
        int *flags1 = find_near_self_intersections(polygons1, 0.75, &n_self_hits);
        for (v = 0; v < polygons1->n_points; v++) {
            if (flags1[v]) {
                Point_x(polygons1->points[v]) -= displacement_field1[v][0];
                Point_y(polygons1->points[v]) -= displacement_field1[v][1];
                Point_z(polygons1->points[v]) -= displacement_field1[v][2];
            }
        }
        compute_polygon_normals(polygons1);

        // Apply smoothing to the (white) displacement field
        smooth_displacement_field(displacement_field2, polygons2, n_neighbours, neighbours, 5, sigma);

        // Apply the final displacement to vertices
        for (v = 0; v < polygons2->n_points; v++) {
            Point_x(polygons2->points[v]) += displacement_field2[v][0];
            Point_y(polygons2->points[v]) += displacement_field2[v][1];
            Point_z(polygons2->points[v]) += displacement_field2[v][2];
        }

        int *flags2 = find_near_self_intersections(polygons2, 0.75, &n_self_hits);
        for (v = 0; v < polygons2->n_points; v++) {
            if (flags2[v]) {
                Point_x(polygons2->points[v]) -= displacement_field2[v][0];
                Point_y(polygons2->points[v]) -= displacement_field2[v][1];
                Point_z(polygons2->points[v]) -= displacement_field2[v][2];
            }
        }
        compute_polygon_normals(polygons2);

        if (verbose) {
            fprintf(stdout, "\rMesh: deform: iter %03d | Errors: %6.4f/%6.4f", i+1, 
                      sqrt(s1 / polygons1->n_points),sqrt(s2 / polygons1->n_points));
            fflush(stdout);
        }
        s1_prev = s1;
        s2_prev = s2;
    }
    if (verbose) fprintf(stdout, "\n");

    for (v = 0; v < polygons1->n_points; v++) {
        displacement_field1[v][0] = Point_x(polygons1->points[v]) - Point_x(polygons1_orig->points[v]);
        displacement_field1[v][1] = Point_y(polygons1->points[v]) - Point_y(polygons1_orig->points[v]);
        displacement_field1[v][2] = Point_z(polygons1->points[v]) - Point_z(polygons1_orig->points[v]);
        displacement_field2[v][0] = Point_x(polygons2->points[v]) - Point_x(polygons2_orig->points[v]);
        displacement_field2[v][1] = Point_y(polygons2->points[v]) - Point_y(polygons2_orig->points[v]);
        displacement_field2[v][2] = Point_z(polygons2->points[v]) - Point_z(polygons2_orig->points[v]);
    }

    double *displacement1 = malloc(sizeof(double) * polygons1->n_points);
    double *displacement2 = malloc(sizeof(double) * polygons2->n_points);

    double prctile[2] = {99.9, 99.9};
    
    /* Get percentile for x-displacement */
if (0) {
      for (v = 0; v < polygons1->n_points; v++)
        displacement1[v] = displacement_field1[v][0];
    double threshold_x1[2];
    get_prctile(displacement1, polygons1->n_points, threshold_x1, prctile, 1, DT_FLOAT64);
    for (v = 0; v < polygons2->n_points; v++)
        displacement2[v] = displacement_field2[v][0];
    double threshold_x2[2];
    get_prctile(displacement2, polygons2->n_points, threshold_x2, prctile, 1, DT_FLOAT64);
    
    /* Get percentile for y-displacement */
    for (v = 0; v < polygons1->n_points; v++)
        displacement1[v] = displacement_field1[v][1];
    double threshold_y1[2];
    get_prctile(displacement1, polygons1->n_points, threshold_y1, prctile, 1, DT_FLOAT64);
    for (v = 0; v < polygons2->n_points; v++)
        displacement2[v] = displacement_field2[v][1];
    double threshold_y2[2];
    get_prctile(displacement2, polygons2->n_points, threshold_y2, prctile, 1, DT_FLOAT64);

    /* Get percentile for z-displacement */
    for (v = 0; v < polygons1->n_points; v++)
        displacement1[v] = displacement_field1[v][2];
    double threshold_z1[2];
    get_prctile(displacement1, polygons1->n_points, threshold_z1, prctile, 1, DT_FLOAT64);
    for (v = 0; v < polygons2->n_points; v++)
        displacement2[v] = displacement_field2[v][2];
    double threshold_z2[2];
    get_prctile(displacement2, polygons2->n_points, threshold_z2, prctile, 1, DT_FLOAT64);

    // Get the squared sum of displacement
    for (v = 0; v < polygons1->n_points; v++) {
        displacement1[v] = SQR(displacement_field1[v][0]) + 
                           SQR(displacement_field1[v][1]) + 
                           SQR(displacement_field1[v][2]);
    }
    
    /* Get percentile for squared sum of displacements */
    double threshold1[2];
    get_prctile(displacement1, polygons1->n_points, threshold1, prctile, 1, DT_FLOAT64);

    // Get the squared sum of displacement
    for (v = 0; v < polygons2->n_points; v++) {
        displacement2[v] = SQR(displacement_field2[v][0]) + 
                           SQR(displacement_field2[v][1]) + 
                           SQR(displacement_field2[v][2]);
    }
    
    /* Get percentile for squared sum of displacements */
    double threshold2[2];
    get_prctile(displacement2, polygons2->n_points, threshold2, prctile, 1, DT_FLOAT64);


    /* If squaed sum of displacements is exceeding 95% percentile we limit all
      displacements to 95% percentile to prevent that these outliers cause 
      self-intersections */
    for (v = 0; v < polygons1->n_points; v++) {
        if (displacement1[v] > threshold1[1]) {
            displacement_field1[v][0] = threshold_x1[1];
            displacement_field1[v][1] = threshold_y1[1];
            displacement_field1[v][2] = threshold_z1[1];
        }
    }

    for (v = 0; v < polygons2->n_points; v++) {
        if (displacement2[v] > threshold2[1]) {
            displacement_field2[v][0] = threshold_x2[1];
            displacement_field2[v][1] = threshold_y2[1];
            displacement_field2[v][2] = threshold_z2[1];
        }
    }
}
    // Apply very slight smoothing to the displacement field to smooth replacements
    // Estimate weight w.r.t. curvature and esnure minimum of min_weight
    smooth_displacement_field(displacement_field1, polygons1, n_neighbours, neighbours, 5, 0.25*sigma);
    smooth_displacement_field(displacement_field2, polygons2, n_neighbours, neighbours, 5, 0.25*sigma);

    // Apply the final displacement to vertices
    for (v = 0; v < polygons1->n_points; v++) {
        Point_x(polygons1_orig->points[v]) += displacement_field1[v][0];
        Point_y(polygons1_orig->points[v]) += displacement_field1[v][1];
        Point_z(polygons1_orig->points[v]) += displacement_field1[v][2];
    }
    copy_polygons(polygons1_orig, polygons1);

    for (v = 0; v < polygons2->n_points; v++) {
        Point_x(polygons2_orig->points[v]) += displacement_field2[v][0];
        Point_y(polygons2_orig->points[v]) += displacement_field2[v][1];
        Point_z(polygons2_orig->points[v]) += displacement_field2[v][2];
    }
    copy_polygons(polygons2_orig, polygons2);

    if (verbose) fprintf(stdout,"\n");
    remove_near_intersections(polygons1, 0.75, verbose);
    remove_near_intersections(polygons2, 0.75, verbose);

    // Free allocated memory
    free(gradient_x);
    free(gradient_y);
    free(gradient_z);
    free(displacement1);
    free(displacement2);
    free(displacement_field1);
    free(displacement_field2);
    delete_polygon_point_neighbours(polygons1, n_neighbours, neighbours, NULL, NULL);
}

