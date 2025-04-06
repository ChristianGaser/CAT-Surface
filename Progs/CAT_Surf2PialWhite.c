/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 *
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

/* argument defaults */
double w1 = 0.1;  // Internal smoothness force
double w2 = 0.2;  // Gradient alignment force
double w3 = 0.4; // Balloon force (expansion/contraction)
double w4 = 0.0; // Connection force
double sigma = 0.3; 
int verbose = 0; 
int iterations = 50; 

/* the argument table */
static ArgvInfo argTable[] = {
    {"-w1", ARGV_FLOAT, (char *) TRUE, (char *) &w1, "Set internal smoothness weight (w1)."},
    {"-w2", ARGV_FLOAT, (char *) TRUE, (char *) &w2, "Set gradient alignment weight (w2)."},
    {"-w3", ARGV_FLOAT, (char *) TRUE, (char *) &w3, "Set balloon force weight (w3)."},
    {"-w4", ARGV_FLOAT, (char *) TRUE, (char *) &w4, "Set conncetion force weight (w4)."},
    {"-sigma", ARGV_FLOAT, (char *) TRUE, (char *) &sigma, "Define sigma for smoothing the displacement field."},
    {"-iter", ARGV_INT, (char *) TRUE, (char *) &iterations, "Set number of deformation iterations."},
    {"-verbose", ARGV_CONSTANT, (char *) TRUE, (char *) &verbose, "Enable verbose output."},
   {NULL, ARGV_END, NULL, NULL, NULL}
};

void
usage(char *executable)
{
    static char *usage_str = "\n\
Usage: %s  surface_file thickness_file label_file output_pial_file output_white_file\n\
Estimate pial or white surface from central surface using cortical thickness values. In order to estimate the pial surface an extent of 0.5 (default) should be used, while an extent of -0.5 results in the estimation of the white matter surface.\n\
The equi-volume model optionally allows to correct the position of the surface around gyri and sulci. The area of the inner (white) and outer (pial) surface is used for this correction.\n\
Furthermore, you can weight the extent of equi-volume correction which is helpful to correct the initial central surface in CAT12 in highly folded areas with high mean curvature.\n\n";

    fprintf(stderr, usage_str, executable);
}

int
main(int argc, char *argv[])
{
    int p, n_objects, n_values;
    int *n_neighbours, **neighbours;
    double *thickness_values, *extents;
    double value, surface_area, pos;
    float *labels;
    Status status;
    char *src_file, *pial_file, *white_file, *values_file, *label_file;
    File_formats format;
    nifti_image *nii_ptr;
    polygons_struct *polygons;
    object_struct **object_list, **objects_out;

    /* get the arguments from the command line */
    if (ParseArgv(&argc, argv, argTable, 0)) {
        usage(argv[0]);
        fprintf(stderr, "   %s -help\n\n", argv[0]);
        exit(EXIT_FAILURE);
    }

    initialize_argument_processing(argc, argv);
    
    if (!get_string_argument( NULL, &src_file) ||
      !get_string_argument( NULL, &values_file) ||
      !get_string_argument( NULL, &label_file) ||
      !get_string_argument( NULL, &pial_file) ||
      !get_string_argument( NULL, &white_file)) {
        usage(argv[0]);
        fprintf(stderr, "   %s -help\n\n", argv[0]);
        exit(EXIT_FAILURE);
    }

    if (input_graphics_any_format(src_file, &format, &n_objects,
                    &object_list) != OK)
        exit(EXIT_FAILURE);

    if (input_values_any_format(values_file, &n_values, &thickness_values) != OK)
        exit(EXIT_FAILURE);
    
    if (n_objects > 1) {
        fprintf(stderr,"Only one object allowed.\n");
        exit(EXIT_FAILURE);
    }

    polygons = get_polygons_ptr(object_list[0]);

    if (polygons->n_points != n_values) {
        fprintf(stderr,"Number of points differs from number of values.\n");
        exit(EXIT_FAILURE);
    }

    nii_ptr = read_nifti_float(label_file, &labels, 0);
    if (!nii_ptr) {
        fprintf(stderr, "Error reading %s.\n", label_file);
        return (EXIT_FAILURE);
    }

    polygons_struct *polygons_pial, *polygons_white, *polygons_smoothed;
    object_struct *object_pial = create_object(POLYGONS);
    object_struct *object_white = create_object(POLYGONS);
    object_struct *object_smoothed = create_object(POLYGONS);

    // Use neg. mean curvature in sulci to weight smoothing
    double *weight = (double *) malloc(sizeof(double) * polygons->n_points);
    get_all_polygon_point_neighbours(polygons, &n_neighbours, &neighbours);
    get_polygon_vertex_curvatures_cg(polygons, n_neighbours, neighbours,
                     3.0, 0.0, weight);
    
    // Estimate weight by using neg. mean curvature and scale it to 0..1
    for (p = 0; p < polygons->n_points; p++) {
        weight[p] = fmin(0.0, weight[p]);
        weight[p] = fmax(-80.0, weight[p]);
        weight[p] /= -80.0;
    }

    extents = (double *) malloc(sizeof(double) * polygons->n_points);

    /* We have to use a much larger extent (theoretically a value of 0.5 would be correct)
       since we smooth displacements to minimize self-intersections */
    for (p = 0; p < polygons->n_points; p++) extents[p] = 1.0;
    objects_out = central_to_new_pial(polygons, thickness_values, extents, 1, 0.5*sigma, 5, verbose);
    object_pial = objects_out[0];
    polygons_pial = get_polygons_ptr(object_pial);

    // Created smoothed pial surface
    polygons_smoothed = get_polygons_ptr(object_smoothed);
    copy_polygons(polygons_pial, polygons_smoothed);
    smooth_heatkernel(polygons_smoothed, NULL, 5.0);
    
    // Weight between pial and smoothe surface w.r.t. curvature
    for (p = 0; p < polygons->n_points; p++) {
        Point_x(polygons_pial->points[p]) = weight[p]*Point_x(polygons_smoothed->points[p]) + 
                                            (1.0-weight[p])*Point_x(polygons_pial->points[p]);
        Point_y(polygons_pial->points[p]) = weight[p]*Point_y(polygons_smoothed->points[p]) +
                                            (1.0-weight[p])*Point_y(polygons_pial->points[p]);
        Point_z(polygons_pial->points[p]) = weight[p]*Point_z(polygons_smoothed->points[p]) +
                                            (1.0-weight[p])*Point_z(polygons_pial->points[p]);
    }

    // Get the white matter surface again with a much larger extent 
    for (p = 0; p < polygons->n_points; p++) extents[p] = -1.0;
    objects_out = central_to_new_pial(polygons, thickness_values, extents, 0, 0.5*sigma, 5, verbose);
    object_white = objects_out[0];
    polygons_white = get_polygons_ptr(object_white);

    // Build weighting matrix and deform pial and white matter surface to isovalues
    double weights[4] = {w1, w2, w3, w4};
    surf_deform_dual(polygons_pial, polygons_white, labels, nii_ptr, 
                      weights, sigma, 1.4, 2.6, thickness_values, iterations, verbose);

    if(output_graphics_any_format(pial_file, format, 1, &object_pial, NULL) != OK)
        exit(EXIT_FAILURE);

    if(output_graphics_any_format(white_file, format, 1, &object_white, NULL) != OK)
        exit(EXIT_FAILURE);
            
    free(extents);
    free(object_pial);
    free(object_white);
    free(weight);
    return(status != OK);
}
