/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Surface reconstruction pipeline: Estimate pial and white matter surfaces
 * from a central surface using cortical thickness and volume information.
 *
 * Copyright Christian Gaser, University of Jena.
 */

#include <bicpl.h>
#include <ParseArgv.h>

#include "CAT_Surf.h"
#include "CAT_Vol.h"
#include "CAT_SurfaceIO.h"
#include "CAT_Smooth.h"
#include "CAT_Deform.h"
#include "CAT_Intersect.h"
#include "CAT_Curvature.h"

// ---------------------------------------------
// Default arguments for the surface deformation
// ---------------------------------------------
double w1 = 0.1;  // Internal smoothness force
double w2 = 0.2;  // Gradient alignment force
double w3 = 0.4;  // Balloon force (based on intensity deviation)
double w4 = 0.0;  // Connection force (between pial and white)
double sigma = 0.3; 
int verbose = 0; 
int iterations = 50; 

// Argument table for command-line parsing
static ArgvInfo argTable[] = {
    {"-w1", ARGV_FLOAT, (char *) TRUE, (char *) &w1, "Set internal smoothness weight (w1)."},
    {"-w2", ARGV_FLOAT, (char *) TRUE, (char *) &w2, "Set gradient alignment weight (w2)."},
    {"-w3", ARGV_FLOAT, (char *) TRUE, (char *) &w3, "Set balloon force weight (w3)."},
    {"-w4", ARGV_FLOAT, (char *) TRUE, (char *) &w4, "Set connection force weight (w4)."},
    {"-sigma", ARGV_FLOAT, (char *) TRUE, (char *) &sigma, "Define sigma for smoothing the displacement field."},
    {"-iter", ARGV_INT, (char *) TRUE, (char *) &iterations, "Set number of deformation iterations."},
    {"-verbose", ARGV_CONSTANT, (char *) TRUE, (char *) &verbose, "Enable verbose output."},
    {NULL, ARGV_END, NULL, NULL, NULL}
};

// ---------------------------------------------
// Print usage/help text
// ---------------------------------------------
void usage(char *executable) {
    fprintf(stderr,
        "\nUsage: %s surface_file thickness_file label_file output_pial_file output_white_file\n"
        "\n"
        "Estimate pial and white matter surfaces from a central surface using:\n"
        "- Cortical thickness values.\n"
        "- A label image that encodes tissue classes.\n\n"
        "This tool performs the following steps:\n"
        "1. Estimate preliminary pial and white surfaces using thickness.\n"
        "2. Smooth pial surface with curvature-guided correction.\n"
        "3. Perform joint deformation of both surfaces using the image intensity and gradient field,\n"
        "   while preserving topology and maintaining surface distance (cortical thickness).\n\n"
        "Key deformation forces:\n"
        "  -w1  Internal smoothness term (e.g. 0.1).\n"
        "  -w2  Gradient alignment force (edges attraction).\n"
        "  -w3  Balloon force, based on isovalue distance.\n"
        "  -w4  Connection force to maintain consistent pial/white spacing.\n"
        "  -sigma  Controls displacement smoothing.\n"
        "  -iter   Number of iterations (e.g. 50).\n"
        "  -verbose  Show iteration-wise output.\n\n",
        executable
    );
}

// ---------------------------------------------
// Main entry point
// ---------------------------------------------
int main(int argc, char *argv[]) {
    int p, n_objects, n_values;
    char *src_file, *pial_file, *white_file, *values_file, *label_file;
    float *labels;
    double *thickness_values, *extents;
    double *weight;
    File_formats format;
    Status status;
    nifti_image *nii_ptr;
    polygons_struct *polygons;
    object_struct **object_list, **objects_out;

    // Parse command-line arguments
    if (ParseArgv(&argc, argv, argTable, 0)) {
        usage(argv[0]);
        exit(EXIT_FAILURE);
    }

    initialize_argument_processing(argc, argv);

    // Input filenames
    if (!get_string_argument(NULL, &src_file) ||
        !get_string_argument(NULL, &values_file) ||
        !get_string_argument(NULL, &label_file) ||
        !get_string_argument(NULL, &pial_file) ||
        !get_string_argument(NULL, &white_file)) {
        usage(argv[0]);
        exit(EXIT_FAILURE);
    }

    // Load surface and check
    if (input_graphics_any_format(src_file, &format, &n_objects, &object_list) != OK)
        exit(EXIT_FAILURE);

    if (n_objects > 1) {
        fprintf(stderr, "Error: Only one surface object allowed.\n");
        exit(EXIT_FAILURE);
    }

    polygons = get_polygons_ptr(object_list[0]);

    // Load thickness values
    if (input_values_any_format(values_file, &n_values, &thickness_values) != OK)
        exit(EXIT_FAILURE);

    if (polygons->n_points != n_values) {
        fprintf(stderr, "Error: Number of surface vertices does not match number of thickness values.\n");
        exit(EXIT_FAILURE);
    }

    // Load NIfTI label volume
    nii_ptr = read_nifti_float(label_file, &labels, 0);
    if (!nii_ptr) {
        fprintf(stderr, "Error reading label volume: %s.\n", label_file);
        return EXIT_FAILURE;
    }

    // Allocate memory for surfaces and weights
    object_struct *object_pial = create_object(POLYGONS);
    object_struct *object_white = create_object(POLYGONS);
    object_struct *object_smoothed = create_object(POLYGONS);

    polygons_struct *polygons_pial;
    polygons_struct *polygons_white;
    polygons_struct *polygons_smoothed;

    int *n_neighbours;
    int **neighbours;

    get_all_polygon_point_neighbours(polygons, &n_neighbours, &neighbours);

    weight = malloc(sizeof(double) * polygons->n_points);

    // Use negative mean curvature to drive smoothing
    get_polygon_vertex_curvatures_cg(polygons, n_neighbours, neighbours, 3.0, 0.0, weight);
    for (p = 0; p < polygons->n_points; p++) {
        weight[p] = fmin(0.0, weight[p]);   // Only negative curvatures
        weight[p] = fmax(-90.0, weight[p]); // Clip at -90
        weight[p] /= -90.0;                 // Normalize to [0..1]
    }

    // Initial estimate of pial surface
    extents = malloc(sizeof(double) * polygons->n_points);
    for (p = 0; p < polygons->n_points; p++) extents[p] = 0.75;

    objects_out = central_to_new_pial(polygons, thickness_values, extents, 1, 0.5*sigma, 5, verbose);
    object_pial = objects_out[0];
    polygons_pial = get_polygons_ptr(object_pial);

    // Smooth pial surface based on local curvature
    polygons_smoothed = get_polygons_ptr(object_smoothed);
    copy_polygons(polygons_pial, polygons_smoothed);
    smooth_heatkernel(polygons_smoothed, NULL, 5.0);

    // Blend original and smoothed surface using curvature-based weights
    for (p = 0; p < polygons->n_points; p++) {
        Point_x(polygons_pial->points[p]) = weight[p]*Point_x(polygons_smoothed->points[p]) + 
                                            (1.0-weight[p])*Point_x(polygons_pial->points[p]);
        Point_y(polygons_pial->points[p]) = weight[p]*Point_y(polygons_smoothed->points[p]) + 
                                            (1.0-weight[p])*Point_y(polygons_pial->points[p]);
        Point_z(polygons_pial->points[p]) = weight[p]*Point_z(polygons_smoothed->points[p]) + 
                                            (1.0-weight[p])*Point_z(polygons_pial->points[p]);
    }

    // Initial estimate of white surface
    for (p = 0; p < polygons->n_points; p++) extents[p] = -0.5;
    objects_out = central_to_new_pial(polygons, thickness_values, extents, 0, 0.5*sigma, 5, verbose);
    object_white = objects_out[0];
    polygons_white = get_polygons_ptr(object_white);

    // Perform final dual-surface deformation using slightly deviating isovalues which
    // are optimized w.r.t. used parameters
    double weights[4] = {w1, w2, w3, w4};
    double shifting[2] = {-0.1, 0.1};
    surf_deform_dual(polygons_pial, polygons_white, polygons, labels, nii_ptr,
                     weights, sigma, CGM+shifting[0], GWM+shifting[1], thickness_values,
                     iterations, verbose);

    // Save output surfaces
    if (output_graphics_any_format(pial_file, format, 1, &object_pial, NULL) != OK ||
        output_graphics_any_format(white_file, format, 1, &object_white, NULL) != OK) {
        fprintf(stderr, "Error writing output surfaces.\n");
        exit(EXIT_FAILURE);
    }

    // Cleanup
    free(extents);
    free(object_pial);
    free(object_white);

    return OK;
}
