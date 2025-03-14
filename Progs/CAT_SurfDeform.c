/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 *
 */

#include "CAT_SurfaceIO.h"
#include "CAT_Surf.h"
#include "CAT_Intersect.h"
#include "CAT_NiftiLib.h"

void
usage(char *executable)
{
    fprintf(stderr, "%s volume_file surface_file output_surface_file [w[0] w[1] w[2] threshold iterations]\n", executable);
}

/**
 * @brief Deforms a 3D surface mesh using an external force field derived from an input volume.
   This approach deforms a 3D surface mesh by computing forces from a reference image.
   It moves each vertex based on a balance of:
    - Internal forces → Maintain smoothness.
    - External forces → Derived from image intensity gradients.
    - Self-intersection prevention → Ensures a valid mesh.
   This method is related to Active Contour Models (Snakes) and Level Set Methods, but explicitly operates on a mesh representation.
 
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
void surf_deform_orig(polygons_struct *polygons, float *input, nifti_image *nii_ptr, 
                 float w[3], float lim, int it, int check_intersections, int verbose)
{
    int i, j, k, v, dims[3], nvox, pidx, counter;
    int *n_neighbours, **neighbours;
    int *defects, *polydefects, n_intersects;
    float *gradient_x, *gradient_y, *gradient_z;
    float d[3], n[3], p[3];
    double vx[3];
    Point points[MAX_POINTS_PER_POLYGON];
    polygons_struct  *polygons_deform;

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

    // Check for memory allocation failure
    if (!gradient_x || !gradient_y || !gradient_z) {
        fprintf(stderr, "Memory allocation error\n");
        exit(EXIT_FAILURE);
    }

    // Compute gradient of the input volume
    gradient3D(input, NULL, gradient_x, gradient_y, gradient_z, dims, vx);

    // Compute surface normals
    compute_polygon_normals(polygons);
        
    // Create neighbor structure for each vertex
    create_polygon_point_neighbours(polygons, TRUE, &n_neighbours, &neighbours, NULL, NULL);

    // Check for selfintersections
    if (check_intersections) {
        defects = (int *) malloc(sizeof(int) * polygons->n_points);
        polydefects = (int *) malloc(sizeof(int) * polygons->n_items);

        /* Apply deformationq and find self intersections to limit search for it
           to these areas to save time */

        polygons_deform = (polygons_struct *) malloc(sizeof(polygons_struct));
        copy_polygons(polygons, polygons_deform);
        
        surf_deform_orig(polygons_deform, input, nii_ptr, w, lim, it, 
                    0, 0);
        n_intersects = find_selfintersections_masked(polygons_deform, defects, polydefects);
//            n_intersects = find_selfintersections(polygons_deform, defects, polydefects, 1);
        
        // Indicate that these defects can be skipped to save time
        for (v = 0; v < polygons->n_items; v++)
            if (polydefects[v] <= 0) polydefects[v] = -1;
                  
        free(polygons_deform);
    }

    double s = 0.0;
    for (i = 0; i < it; i++) {
        s = 0.0;
        if (check_intersections && (i % check_intersections) == 0)
            n_intersects = find_selfintersections(polygons, defects, polydefects, 0);

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
            p[0] = Point_x(polygons->points[v]);
            p[1] = Point_y(polygons->points[v]);
            p[2] = Point_z(polygons->points[v]);

            // Get vertex normal
            n[0] = Point_x(polygons->normals[v]);
            n[1] = Point_y(polygons->normals[v]);
            n[2] = Point_z(polygons->normals[v]);

            // Compute external force based on image gradient
            float di = isoval(input,      p[0], p[1], p[2], dims, nii_ptr) - lim;
            float fx = isoval(gradient_x, p[0], p[1], p[2], dims, nii_ptr);
            float fy = isoval(gradient_y, p[0], p[1], p[2], dims, nii_ptr);
            float fz = isoval(gradient_z, p[0], p[1], p[2], dims, nii_ptr);            
            float f3 = ((di / 10.0));
            float f2 = fmax(-1.0, fmin(1.0, fx * n[0] + fy * n[1] + fz * n[2]));
        
            // Compute vertex displacement
            for (k = 0; k < 3; k++) {
                d[k] = w[0] * (c[k] - p[k]) + ((w[1] * f2 + w[2]) * f3) * n[k];
            }

            // Move vertex while preventing intersections
            double new_p[3] = {p[0] + d[0], p[1] + d[1], p[2] + d[2]};
            
            if (check_intersections) {
                if (n_intersects == 0 || defects[v] <= 0)
                    fill_Point(polygons->points[v], new_p[0], new_p[1], new_p[2]);
                
            } else  fill_Point(polygons->points[v], new_p[0], new_p[1], new_p[2]);

            s += di * di;
        }

        // Update normals for next iteration
        compute_polygon_normals(polygons);
        if (verbose) fprintf(stdout, "Mesh: deform: iter %03d %g\n", i, sqrt(s / ((float)polygons->n_points)));
    }
    
    // Free allocated memory
    free(gradient_x);
    free(gradient_y);
    free(gradient_z);

    if (check_intersections) {
        free(defects);
        free(polydefects);
    }
        
    delete_polygon_point_neighbours(polygons, n_neighbours, neighbours, NULL, NULL);
}

void enforce_diffeomorphism(double (*displacement_field)[3], polygons_struct *polygons, 
                            int *n_neighbours, int **neighbours, double min_det) {
    int v, k, max_iterations = 4; // Prevent infinite loops
    double J[3][3], detJ, scale_factor;

    for (v = 0; v < polygons->n_points; v++) {
        int iter = 0;  // Track iterations per vertex

        do {
            // Compute the Jacobian matrix for this vertex
            for (k = 0; k < 3; k++) {
                J[k][0] = displacement_field[v][k] - displacement_field[neighbours[v][0]][k];
                J[k][1] = displacement_field[v][k] - displacement_field[neighbours[v][1]][k];
                J[k][2] = displacement_field[v][k] - displacement_field[neighbours[v][2]][k];
            }

            // Compute determinant of Jacobian
            detJ = J[0][0] * (J[1][1] * J[2][2] - J[1][2] * J[2][1]) -
                   J[0][1] * (J[1][0] * J[2][2] - J[1][2] * J[2][0]) +
                   J[0][2] * (J[1][0] * J[2][1] - J[1][1] * J[2][0]);

            // If determinant is below threshold, reduce displacement
            if (detJ < min_det) {
                scale_factor = 0.9; // Reduce displacement by 10%
                for (k = 0; k < 3; k++) {
                    displacement_field[v][k] *= scale_factor;
                }
            }

            iter++;
        } while (detJ < min_det && iter < max_iterations);
    }
}

void smooth_displacement_field(double (*displacement_field)[3], polygons_struct *polygons, 
                               int *n_neighbours, int **neighbours, int iterations, double sigma) {
    int v, j, k, pidx;
    double new_displacement[polygons->n_points][3];

    for (int it = 0; it < iterations; it++) {
        for (v = 0; v < polygons->n_points; v++) {
            double smoothed[3] = {0.0, 0.0, 0.0};
            int count = 0;

            for (j = 0; j < n_neighbours[v]; j++) {
                pidx = neighbours[v][j];
                for (k = 0; k < 3; k++) {
                    smoothed[k] += displacement_field[pidx][k];
                }
                count++;
            }

            for (k = 0; k < 3; k++) {
                new_displacement[v][k] = (smoothed[k] / count) * exp(-sigma);
            }
        }

        for (v = 0; v < polygons->n_points; v++) {
            for (k = 0; k < 3; k++) {
                displacement_field[v][k] = new_displacement[v][k];
            }
        }
    }
}

void surf_deform(polygons_struct *polygons, float *input, nifti_image *nii_ptr, 
                 float w[3], float lim, int it, int verbose)
{
    int i, j, k, v, dims[3], nvox, pidx;
    int *n_neighbours, **neighbours;
    float *gradient_x, *gradient_y, *gradient_z;
    double vx[3], s;
    Point points[MAX_POINTS_PER_POLYGON];

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
            float f3 = ((di / 10.0));
            float f2 = fmax(-1.0, fmin(1.0, fx * n[0] + fy * n[1] + fz * n[2]));

            // Compute vertex displacement and store in the displacement field
            for (k = 0; k < 3; k++) {
                displacement_field[v][k] = w[0] * (c[k] - p[k]) + ((w[1] * f2 + w[2]) * f3) * n[k];
            }

            s += di * di;
        }

        // Apply smoothing to the displacement field
        smooth_displacement_field(displacement_field, polygons, n_neighbours, neighbours, 5, 0.1);

        // Enforce diffeomorphic constraints
        enforce_diffeomorphism(displacement_field, polygons, n_neighbours, neighbours, 0.0001);

        // Apply the final displacement to vertices
        for (v = 0; v < polygons->n_points; v++) {
            Point_x(polygons->points[v]) += displacement_field[v][0];
            Point_y(polygons->points[v]) += displacement_field[v][1];
            Point_z(polygons->points[v]) += displacement_field[v][2];
        }

        // Update normals for next iteration
        compute_polygon_normals(polygons);
        if (verbose) fprintf(stdout, "Mesh: deform: iter %03d %g\n", i, sqrt(s / polygons->n_points));
    }

    // Free allocated memory
    free(gradient_x);
    free(gradient_y);
    free(gradient_z);
    free(displacement_field);
    delete_polygon_point_neighbours(polygons, n_neighbours, neighbours, NULL, NULL);
}


int
main(int argc, char *argv[])
{
    char *volume_file = NULL;
    char *input_file = NULL, *output_surface_file = NULL;
    float *input, w[3], lim;
    int n_objects, it, check_intersections = 10;
    File_formats file_format;
    object_struct **object_list;
    polygons_struct *polygons;
    nifti_image *nii_ptr;

    initialize_argument_processing(argc, argv);

    get_string_argument( NULL, &volume_file);

    /* read first image to get image parameters */
    nii_ptr = read_nifti_float(volume_file, &input, 0);
    if(nii_ptr == NULL) {
        fprintf(stderr,"Error reading %s.\n", volume_file);
        return(EXIT_FAILURE);
    }

    get_string_argument( NULL, &input_file);
    get_string_argument( NULL, &output_surface_file);
    
    get_real_argument(0.01, &w[0]);
    get_real_argument(0.1, &w[1]);
    get_real_argument(0.2, &w[2]);
    get_real_argument(0.5, &lim);
    get_int_argument(100, &it);

    if (input_graphics_any_format(input_file, &file_format,
              &n_objects, &object_list) == ERROR ||
              n_objects != 1 || object_list[0]->object_type != POLYGONS) {
        fprintf(stderr, "File must contain 1 polygons struct.\n");
        exit(EXIT_FAILURE);
    }

    polygons = get_polygons_ptr(object_list[0]);    

    surf_deform(polygons, input, nii_ptr, w, lim, it, 1);
    //remove_intersections(polygons, 1);
    
    if (output_graphics_any_format(output_surface_file, ASCII_FORMAT,
                     n_objects, object_list, NULL) == ERROR)
        exit(EXIT_FAILURE);

    return(EXIT_SUCCESS);
}
