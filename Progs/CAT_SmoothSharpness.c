/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 *
 */

/*
 * Heat kernel smoothing is based on matlab code from Moo K. Chung:
 * Chung, M.K., Robbins,S., Dalton, K.M., Davidson, R.J., Evans, A.C. (2005) 
 * Cortical thickness analysis in autism via heat kernel smoothing. NeuroImage. 
 * http://www.stat.wisc.edu/~mchung/papers/ni_heatkernel.pdf
 */

#include <bicpl.h>

#include "CAT_Smooth.h"
#include "CAT_SurfaceIO.h"
#include "CAT_Curvature.h"

#ifndef MIN
#define MIN(A,B) ((A) > (B) ? (B) : (A))
#endif

#ifndef MAX
#define MAX(A,B) ((A) > (B) ? (A) : (B))
#endif

void
usage(char *executable)
{
    char *usage_str = "\n\
Usage: %s surface_file output_surface_file fwhm\n\n\
   Diffusion smoothing of surface points w.r.t. sharpness.\n\n";

    fprintf(stderr, usage_str, executable);
}

int
main(int argc, char *argv[])
{
    char       *input_file, *output_surface_file;
    int        n_objects;
    int        *n_neighbours, **neighbours;
    int        p, k;
    File_formats   format;
    object_struct  **object_list;
    polygons_struct  *polygons, *smoothed_polygons;
    float      fwhm;
    double       *sharpness, min, max;

    initialize_argument_processing(argc, argv);

    if (!get_string_argument(NULL, &input_file) ||
      !get_string_argument(NULL, &output_surface_file)) {
        usage(argv[0]);
        exit(EXIT_FAILURE);
    }

    get_real_argument(25.0, &fwhm);

    if (input_graphics_any_format(input_file, &format, &n_objects,
                    &object_list) != OK ||
      n_objects != 1 || get_object_type(object_list[0]) != POLYGONS) {
        fprintf(stderr, "Error reading %s.\n", input_file);
        exit(EXIT_FAILURE);
    }

    polygons = get_polygons_ptr(object_list[0]);

    sharpness = (double *)malloc(sizeof(double)*polygons->n_points);
      
    get_all_polygon_point_neighbours(polygons, &n_neighbours, &neighbours);

    compute_local_sharpness(polygons, n_neighbours, neighbours, sharpness);
    
    /* use squared inverted sharpness as weighting */
    min = 0.0; max = 0.0;
    for (p = 0; p < polygons->n_points; p++) {
        /* use srq to locally change weighting */
        sharpness[p] = pow(sharpness[p],2.0);
        min = MIN(sharpness[p], min);
        max = MAX(sharpness[p], max);
    }

    /* scale weighting to 0..1 */
    for (p = 0; p < polygons->n_points; p++)
        sharpness[p] = ((sharpness[p] - min)/(max - min));

    smoothed_polygons = get_polygons_ptr(create_object(POLYGONS));
    copy_polygons(polygons, smoothed_polygons);
    smooth_heatkernel(smoothed_polygons, NULL, fwhm);
    
    /* weighted averaging between unsmoothed and smoothed surface */
    for (p = 0; p < polygons->n_points; p++) {
        for (k = 0; k < 3; k++) 
            Point_coord(polygons->points[p], k) = sharpness[p]*Point_coord(smoothed_polygons->points[p], k) + ((1-sharpness[p])*Point_coord(polygons->points[p], k));
    }

    compute_polygon_normals(polygons);

    if(output_graphics_any_format(output_surface_file, format, 1, 
            object_list, NULL) != OK)
        exit(EXIT_FAILURE);
    
    free(sharpness);

    return(EXIT_SUCCESS);
}
