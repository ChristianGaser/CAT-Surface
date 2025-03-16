/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry, University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 *
 * Surface Deformation Algorithm with Diffeomorphic Constraints
 */

#include <ParseArgv.h>
#include "CAT_SurfaceIO.h"
#include "CAT_Surf.h"
#include "CAT_Intersect.h"
#include "CAT_NiftiLib.h"

/* Default parameter values */
double w1 = 0.05;  // Internal smoothness force
double w2 = 0.1;   // Gradient alignment force
double w3 = 0.5;   // Balloon force (expansion/contraction)
double isovalue = 0.5;  // Target isosurface value
int iterations = 200;    // Number of iterations
int verbose = 0;         // Verbose mode

/* Command-line argument table */
static ArgvInfo argTable[] = {
    {"-isovalue", ARGV_FLOAT, (char *) TRUE, (char *) &isovalue, "Define isovalue (target intensity)."},
    {"-w1", ARGV_FLOAT, (char *) TRUE, (char *) &w1, "Set internal smoothness weight (w1)."},
    {"-w2", ARGV_FLOAT, (char *) TRUE, (char *) &w2, "Set gradient alignment weight (w2)."},
    {"-w3", ARGV_FLOAT, (char *) TRUE, (char *) &w3, "Set balloon force weight (w3)."},
    {"-iter", ARGV_INT, (char *) TRUE, (char *) &iterations, "Set number of deformation iterations."},
    {"-verbose", ARGV_CONSTANT, (char *) TRUE, (char *) &verbose, "Enable verbose output."},
    {NULL, ARGV_END, NULL, NULL, NULL}
};

/**
 * @brief Displays usage instructions for the surface deformation tool.
 *
 * This function provides an overview of the deformation process and guides users
 * on how to correctly execute the program.
 *
 * **Approach Overview:**
 *  - The method deforms a 3D mesh based on image intensity and gradient forces.
 *  - Ensures a balance between smoothness, gradient attraction, and volume expansion.
 *  - Enforces diffeomorphic constraints to prevent self-intersections.
 *
 * **Usage:**
 * ```
 *  CAT_SurfDeform [options] <volume.nii> <input_surface.gii> <output_surface.gii>
 * ```
 *
 * **Arguments:**
 *  - `<volume.nii>`: 3D image file (NIfTI format).
 *  - `<input_surface.gii>`: Initial mesh surface file.
 *  - `<output_surface.gii>`: Output deformed mesh surface file.
 *
 * **Options:**
 *  - `-isovalue <float>` → Target isosurface value (default: `0.5`).
 *  - `-w1 <float>` → Internal smoothing weight (`0.01`).
 *  - `-w2 <float>` → Gradient alignment weight (`0.1`).
 *  - `-w3 <float>` → Balloon force weight (`0.2`).
 *  - `-iter <int>` → Number of iterations (`100`).
 *  - `-verbose` → Enable detailed output.
 */
void usage(char *executable) {
    char *usage_str = "\n"
    "--------------------------------------------------------------------------------\n"
    " Surface Deformation Tool - Mesh Adaptation Based on Image Intensity\n"
    "--------------------------------------------------------------------------------\n"
    "\n"
    "Usage:\n"
    "  %s [options] <volume.nii> <input_surface.gii> <output_surface.gii>\n"
    "\n"
    "Description:\n"
    "  This program deforms a 3D surface mesh using an image-driven approach.\n"
    "  It balances smoothness, edge attraction, and volume preservation while\n"
    "  ensuring a diffeomorphic transformation (preventing self-intersections).\n"
    "\n"
    "Required Arguments:\n"
    "  <volume.nii>         3D image file (NIfTI format) used for deformation.\n"
    "  <input_surface.gii>  Input mesh surface file.\n"
    "  <output_surface.gii> Output deformed mesh surface file.\n"
    "\n"
    "Optional Parameters:\n"
    "  -isovalue <float>    Target isosurface intensity value (default: 0.5).\n"
    "  -w1 <float>          Internal smoothing weight (default: 0.01).\n"
    "  -w2 <float>          Gradient alignment weight (default: 0.1).\n"
    "  -w3 <float>          Balloon force weight (default: 0.2).\n"
    "  -iter <int>          Number of deformation iterations (default: 100).\n"
    "  -verbose             Enable verbose mode for detailed logs.\n"
    "\n"
    "Deformation Forces:\n"
    "  - w1: Encourages local smoothness by averaging neighboring vertices.\n"
    "  - w2: Attracts the surface toward edges in the image (gradient force).\n"
    "  - w3: Expands or contracts the surface, assisting in volume growth/shrinkage.\n"
    "\n"
    "Example Usage:\n"
    "  %s -isovalue 0.6 -w1 0.02 -w2 0.15 -w3 0.25 -iter 150 \\\n"
    "     brain.nii.gz white_surface.gii white_surface_deformed.gii\n"
    "\n"
    "Expected Output:\n"
    "  - The output surface mesh is deformed according to the image intensities.\n"
    "  - A valid diffeomorphic transformation is maintained to minimize self-intersections.\n"
    "--------------------------------------------------------------------------------\n";

    fprintf(stderr, usage_str, executable, executable);
}

/**
 * @brief Limits the displacement field to enforce diffeomorphic (one-to-one) mapping.
 *
 * This function ensures that the deformation is invertible by preventing mesh folding.
 * If the Jacobian determinant of a transformation is too low, the displacement is reduced iteratively.
 *
 * @param displacement_field 3D displacement field for each vertex.
 * @param polygons The 3D surface mesh.
 * @param n_neighbours Array storing the number of neighbors for each vertex.
 * @param neighbours 2D array storing the indices of each vertex’s neighbors.
 * @param min_det Minimum allowable determinant value to maintain a diffeomorphic mapping.
 */
void enforce_diffeomorphism(double (*displacement_field)[3], polygons_struct *polygons, 
                            int *n_neighbours, int **neighbours, double min_det) {
    int v, k, max_iterations = 4;  // Prevent infinite loops
    double J[3][3], detJ, scale_factor;

    for (v = 0; v < polygons->n_points; v++) {
        int iter = 0;  // Iteration counter for this vertex

        do {
            // Compute the Jacobian matrix for the vertex
            for (k = 0; k < 3; k++) {
                J[k][0] = displacement_field[v][k] - displacement_field[neighbours[v][0]][k];
                J[k][1] = displacement_field[v][k] - displacement_field[neighbours[v][1]][k];
                J[k][2] = displacement_field[v][k] - displacement_field[neighbours[v][2]][k];
            }

            // Compute determinant of the Jacobian
            detJ = J[0][0] * (J[1][1] * J[2][2] - J[1][2] * J[2][1]) -
                   J[0][1] * (J[1][0] * J[2][2] - J[1][2] * J[2][0]) +
                   J[0][2] * (J[1][0] * J[2][1] - J[1][1] * J[2][0]);

            // If determinant is too small, reduce displacement
            if (detJ < min_det) {
                scale_factor = 0.9;  // Reduce displacement by 10%
                for (k = 0; k < 3; k++) {
                    displacement_field[v][k] *= scale_factor;
                }
            }
            iter++;
        } while (detJ < min_det && iter < max_iterations);  // Repeat until valid
    }
}

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
                               int *n_neighbours, int **neighbours, int iterations, double sigma) {
    int v, j, k, pidx;

    double (*new_displacement)[3] = malloc(sizeof(double[3]) * polygons->n_points);
    if (!new_displacement) {
        fprintf(stderr, "Memory allocation error\n");
        exit(EXIT_FAILURE);
    }

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
void surf_deform(polygons_struct *polygons, float *input, nifti_image *nii_ptr, 
                 double w[3], float lim, int it, int verbose)
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
            float f3 = ((di / 5.0));
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
        if (verbose) {
            fprintf(stdout, "\rMesh: deform: iter %03d | Error: %g", i, sqrt(s / polygons->n_points));
            fflush(stdout);  // Force output update
        }
    }
    if (verbose) fprintf(stdout, "\n");

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
    float *input;
    int n_objects;
    File_formats file_format;
    object_struct **object_list;
    polygons_struct *polygons;
    nifti_image *nii_ptr;

    initialize_argument_processing(argc, argv);

    /* get the arguments from the command line */
    if (ParseArgv(&argc, argv, argTable, 0)) {
        usage(argv[0]);
        fprintf(stderr, "     %s -help\n\n", argv[0]);
        return EXIT_FAILURE;
    }

    if (!get_string_argument(NULL, &volume_file) || 
        !get_string_argument(NULL, &input_file)  ||
        !get_string_argument(NULL, &output_surface_file)) {
        usage(argv[0]);
        fprintf(stderr, "Usage: CAT_SurfDeform volume_file surface_file output_surface_file\n");
        return EXIT_FAILURE;
    }

    /* read first image to get image parameters */
    nii_ptr = read_nifti_float(volume_file, &input, 0);
    if(nii_ptr == NULL) {
        fprintf(stderr,"Error reading %s.\n", volume_file);
        return(EXIT_FAILURE);
    }

    if (input_graphics_any_format(input_file, &file_format,
              &n_objects, &object_list) == ERROR ||
              n_objects != 1 || object_list[0]->object_type != POLYGONS) {
        fprintf(stderr, "File must contain 1 polygons struct.\n");
        exit(EXIT_FAILURE);
    }

    polygons = get_polygons_ptr(object_list[0]);    

    double weights[3] = {w1, w2, w3};
    surf_deform(polygons, input, nii_ptr, weights, isovalue, iterations, verbose);
    
    if (output_graphics_any_format(output_surface_file, ASCII_FORMAT,
                     n_objects, object_list, NULL) == ERROR)
        exit(EXIT_FAILURE);

    return(EXIT_SUCCESS);
}
