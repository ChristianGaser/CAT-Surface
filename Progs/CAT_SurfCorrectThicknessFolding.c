/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 *
 */

#include <bicpl.h>
#include <float.h>
#include <ParseArgv.h>

#include "CAT_SurfaceIO.h"
#include "CAT_Math.h"
#include "CAT_Curvature.h"
#include "CAT_Smooth.h"

double max_dist = 6.0;           /* maximal thickness */

/* the argument table */
static ArgvInfo argTable[] = {
  {"-max", ARGV_FLOAT, (char *) TRUE, (char *) &max_dist, "Define maximum thickness, where all values exceeding that will be cut."},
  { NULL, ARGV_END, NULL, NULL, NULL }
};

void
usage(char *executable)
{
    char *usage_str = "\n\
Usage: %s  surface_file [options] surface_file thickness_file output_thickness_file\n\n\
   This program corrects cortical thickness that is influenced by folding using\n\
   the ideas from that paper:\n\
   Nagehan Demirci, Timothy S. Coalson, Maria A. Holland, David C. Van Essen, Matthew F. Glasser\n\
   Compensating Cortical Thickness for Cortical Folding-Related Variation\n\
   https://doi.org/10.1101/2025.05.03.651968.\n";

    fprintf(stderr, usage_str, executable);
}
  
int
main(int argc, char *argv[])
{
    double *thickness, *curvatures;
    double **G, **invG, *beta;
    int i, j, k, n_objects, n_vals;
    int *n_neighbours, **neighbours;
    Status status;
    char *surface_file, *output_thickness_file, *thickness_file;
    File_formats format;
    object_struct **objects;
    polygons_struct *polygons;

    /* Call ParseArgv */
    if (ParseArgv(&argc, argv, argTable, 0)) {
        usage(argv[0]);
        fprintf(stderr, "       %s -help\n\n", argv[0]);
        exit(EXIT_FAILURE);
    }

    initialize_argument_processing(argc, argv);

    if (!get_string_argument( NULL, &surface_file) ||
        !get_string_argument( NULL, &thickness_file) ||
        !get_string_argument( NULL, &output_thickness_file)) {
          usage(argv[0]);
          exit(EXIT_FAILURE);
    }

    if (input_graphics_any_format(surface_file, &format, &n_objects,
                    &objects) != OK)
        exit(EXIT_FAILURE);

    if (input_values_any_format(thickness_file, &n_vals, &thickness) != OK)
        exit(EXIT_FAILURE);

    polygons = get_polygons_ptr(objects[0]);

    // remove mean from thickness
    double mean_thickness = get_mean_double(thickness, n_vals, 0);
    for (i = 0; i < n_vals; i++) thickness[i] -= mean_thickness;

    get_all_polygon_point_neighbours(polygons, &n_neighbours, &neighbours);

    curvatures = (double *) malloc(sizeof(double) * n_vals);

    // add gaussian curvature, curvedness, shape index, mean curvature to G
    int curvtype[4] = {1,2,3,4};
    int n_curvtypes = sizeof(curvtype) / sizeof(curvtype[0]);
    int n_beta = 1 + 2*(n_curvtypes);

    ALLOC2D(G, n_vals, n_beta);
    ALLOC2D(invG, n_beta, n_vals);
    beta = malloc(sizeof(double) * n_beta);

    // Create design matrix G by using linear and squared terms for 4 different foldings
    for (j = 0; j < n_curvtypes; j++) {
        get_polygon_vertex_curvatures_cg(polygons, n_neighbours, neighbours,
                         0.0, curvtype[j], curvatures);
        
        // Smooth foldings with FWHM of 3mm
        smooth_heatkernel(polygons, curvatures, 3.0);
        
        // Add  linear term
        normalize_double(curvatures, n_vals);
        for (i = 0; i < n_vals; i++) G[i][2*j] = curvatures[i];

        // Add squared term
        for (i = 0; i < n_vals; i++) curvatures[i] *= curvatures[i];
        normalize_double(curvatures, n_vals);
        for (i = 0; i < n_vals; i++) G[i][2*j+1] = curvatures[i];
    }
    // also add constant as last column
    for (i = 0; i < n_vals; i++) G[i][n_beta-1] = 1.0;

    /* Compute pseudo inverse from design matrix */        
    (void) pinv(n_vals, n_beta, G, invG);

    /* Get betas */
    for (i = 0; i < n_beta; i++) {
        beta[i] = 0.0;
        for (j = 0; j < n_vals; j++)
            beta[i] += invG[i][j] * thickness[j];        
    }

    // Correct thickness by removing effects due to folding
    for (i = 0; i < n_vals; i++) {
        double proj = 0.0;
        for (j = 0; j < n_beta; j++) proj += G[i][j] * beta[j];

        thickness[i] -= proj;        
    }

    // Add mean again to thickness
    for (i = 0; i < n_vals; i++) thickness[i] += mean_thickness;

    clip_data(thickness, n_vals, 0, max_dist, DT_FLOAT64); 
    
    output_values_any_format(output_thickness_file, n_vals,
                 thickness, TYPE_DOUBLE);

    delete_object_list(n_objects, objects);
    free(curvatures);
    free(beta);
    FREE2D(G);
    FREE2D(invG);

    return(OK);

}
