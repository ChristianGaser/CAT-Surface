/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 *
 *   References:
 *
 *  M. Boucher, S. Whitesides, and A. Evans, "Depth potential function for 
 *    folding pattern representation, registration and analysis", Medical
 *    Image Analysis, Volume 13, Issue 2, pp 203-214, April 2009.
 *
 */

#include <bicpl.h>

#include "CAT_Smooth.h"
#include "CAT_Surf.h"
#include "CAT_DepthPotential.h"
#include "CAT_SurfaceIO.h"
  
int
main(int argc, char *argv[])
{
    char         *object_file, *output_surface_file;
    File_formats     format;
    int          poly, n_objects;
    object_struct    **objects;
    polygons_struct    *polygons;
    float        alpha;
    double         *measure;

    initialize_argument_processing(argc, argv);

    if (!get_string_argument(NULL, &object_file) ||
      !get_string_argument(NULL, &output_surface_file)) {
        fprintf(stderr, "Usage: %s  surface_file output_values_file alpha\n",
            argv[0]);
        exit(EXIT_FAILURE);
    }

    if (input_graphics_any_format(object_file, &format,
                    &n_objects, &objects) != OK) {
        fprintf(stderr, "Error reading %s.\n", object_file);
        exit(EXIT_FAILURE);
    }

    if (n_objects != 1 || get_object_type(objects[0]) != POLYGONS) {
        printf("File must contain 1 polygons object.\n");
        exit(EXIT_FAILURE);
    }

    get_real_argument(0.0015, &alpha);

    polygons = get_polygons_ptr(objects[0]);

    measure = compute_depth_potential( polygons, alpha); 
  
    output_values_any_format(output_surface_file, polygons->n_points, measure,
                 TYPE_DOUBLE);
    
    free(measure);
    delete_object_list(n_objects, objects);

    return(EXIT_SUCCESS);
}
