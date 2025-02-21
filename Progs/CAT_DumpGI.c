/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 *
 */

/* Heat kernel smoothing is based on matlab code from Moo K. Chung:
 *  Chung, M.K., Robbins,S., Dalton, K.M., Davidson, R.J., Evans, A.C. (2004) 
 *  Cortical thickness analysis in autism via heat kernel smoothing.
 * NeuroImage, submitted. 
 *  http://www.stat.wisc.edu/~mchung/papers/ni_heatkernel.pdf
 */

#include <bicpl.h>

#include "CAT_Smooth.h"
#include "CAT_Curvature.h"
#include "CAT_SurfaceIO.h"

void
usage(char *executable)
{
    char *usage_str = "\n\
Usage: %s  surface_file output_surface_file fwhm [fwhm_surf] [curvtype]\n\n\
   Calculate gyrification index using the ratio of local surface area\n\
   and local inflated surface area. Local surface area can be approximated\n\
   by use of different curve types (default: mean curvature averaged over\n\
   3mm):\n\
   curvtype:  0 - mean curvature (averaged over 3mm, in degrees)\n\
        2 - curvedness\n\
        3 - shape index\n\
        4 - mean curvature (in radians)\n\n";

    fprintf(stderr, usage_str, executable);
}

int
main(int argc, char *argv[])
{
    char         *object_file, *output_surface_file;
    File_formats     format;
    int          i, n_objects, curvtype;
    int          *n_neighbours, **neighbours;
    object_struct    **objects;
    polygons_struct    *polygons;
    double           fwhm, fwhm_surf, *curvatures, *curvs_inflated;
    double           *GI, distance;

    initialize_argument_processing(argc, argv);

    fprintf(stderr,"Experimental function that is not yet working.\n");

    if (!get_string_argument(NULL, &object_file) ||
      !get_string_argument(NULL, &output_surface_file)) {
        usage(argv[0]);
        exit(EXIT_FAILURE);
    }

    get_real_argument(30.0, &fwhm);
    get_real_argument(100.0, &fwhm_surf);
    get_int_argument(0, &curvtype);
  
    if (input_graphics_any_format(object_file, &format,
                    &n_objects, &objects) != OK)
        exit(EXIT_FAILURE);

    if (n_objects != 1 || get_object_type(objects[0]) != POLYGONS) {
        printf("File must contain 1 polygons object.\n");
        exit(EXIT_FAILURE);
    }

    polygons = get_polygons_ptr(objects[0]);

    compute_polygon_normals(polygons);

    ALLOC(curvatures, polygons->n_points);
    ALLOC(curvs_inflated, polygons->n_points);
    ALLOC(GI, polygons->n_points);

    get_all_polygon_point_neighbours(polygons, &n_neighbours, &neighbours);
  
    // get curvature values
    if (curvtype == 0)
        distance = 3.0;
    else
        distance = 0.0;
    get_polygon_vertex_curvatures_cg(polygons, n_neighbours, neighbours,
                     distance, curvtype, curvatures);

    /* inflate surface by smoothing with FWHM of 150mm */
    smooth_heatkernel(polygons, NULL, fwhm_surf);

    compute_polygon_normals(polygons);

    /* get curvature values of inflated surface */
    get_polygon_vertex_curvatures_cg(polygons, n_neighbours, neighbours,
                     distance, curvtype, curvs_inflated);

    /* calculate ratio of absolute values */
    for (i = 0; i < polygons->n_points; i++) {
        if (curvs_inflated[i] != 0)
            GI[i] = fabs(curvatures[i]/curvs_inflated[i]);
        else
            GI[i] = 0.0;
    }
  
    smooth_heatkernel(polygons, GI, fwhm);
    
    output_values_any_format(output_surface_file, polygons->n_points,
                 GI, TYPE_DOUBLE);

    delete_polygon_point_neighbours(polygons, n_neighbours,
                    neighbours, NULL, NULL);
    delete_object_list(n_objects, objects);
    FREE(curvatures);
    FREE(curvs_inflated);
    FREE(GI);

    return(EXIT_SUCCESS);
}
