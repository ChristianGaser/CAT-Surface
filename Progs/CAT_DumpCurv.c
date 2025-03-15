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
 *  Chung, M.K., Robbins,S., Dalton, K.M., Davidson, R.J., Evans, A.C. (2004) 
 *  Cortical thickness analysis in autism via heat kernel smoothing.
 *  NeuroImage, submitted. 
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
Usage: %s  surface_file output_values_file [curvtype] [fwhm] [use_abs_vals] [-1|1]\n\n\
   Calculate different curvature parameters (default: mean curvature\n\
   averaged over 3mm) from a given obj file.  If smoothing filter [fwhm]\n\
   is defined (in FWHM) a diffusion heat kernel will be applied. Optionally\n\
   the absolute value can ba calculated if use_abs_values is defined and\n\
   only positive or negative values of the values can be saved using the\n\
   last option (-1 or 1).\n\
   curvtype:  0 - mean curvature (averaged over 3mm, in degrees)\n\
        1 - gaussian curvature\n\
        2 - curvedness\n\
        3 - shape index\n\
        4 - mean curvature (in radians)\n\
        5 - sulcal depth like estimator\n\
         >5 - depth potential with alpha = 1/curvtype (recommended value curvtype=650)\n\n";

    fprintf(stderr, usage_str, executable);
}

int
main(int argc, char *argv[])
{
    char *object_file, *output_surface_file;
    File_formats format;
    int i, n_objects, curvtype, sign;
    int *n_neighbours, **neighbours, use_abs_values;
    object_struct **objects;
    polygons_struct *polygons;
    double *curvatures, distance;
    float fwhm;
    BOOLEAN smoothing;

    initialize_argument_processing(argc, argv);

    if (!get_string_argument(NULL, &object_file) ||
      !get_string_argument(NULL, &output_surface_file)) {
        usage(argv[0]);
        exit(EXIT_FAILURE);
    }

    get_int_argument(0, &curvtype);
    get_real_argument(0.0, &fwhm);
    get_int_argument(0, &use_abs_values);
    get_int_argument(0, &sign);
  
    if (fwhm > 0)
        smoothing = 1;
    else smoothing = 0;

    if (input_graphics_any_format(object_file, &format, &n_objects,
                    &objects ) != OK)
        exit(EXIT_FAILURE);

    if (n_objects != 1 || get_object_type(objects[0]) != POLYGONS) {
        printf("File must contain 1 polygons object.\n");
        exit(EXIT_FAILURE);
    }

    polygons = get_polygons_ptr(objects[0]);

    curvatures = (double *) malloc(sizeof(double) * polygons->n_points);
  
    get_all_polygon_point_neighbours(polygons, &n_neighbours, &neighbours);

    if (curvtype == 0)
        distance = 3.0;
    else distance = 0.0;

    get_polygon_vertex_curvatures_cg(polygons, n_neighbours, neighbours,
                     distance, curvtype, curvatures);

    /* limit range to values between -1..1 for all curvtypes > 0 */
    /*  (don't ask me where the large values come from...) */
    if ((curvtype > 0) && (curvtype < 5)) { 
        for (i = 0; i < polygons->n_points; i++) {
            if (curvatures[i] < -1)
                curvatures[i] = -1;
            if (curvatures[i] >  1)
                curvatures[i] = 1;
        }
    }
  
    /* use absolute value */
    if (use_abs_values) { 
        for (i = 0; i < polygons->n_points; i++)
            curvatures[i] = fabs(curvatures[i]);
    }

    /* use positive or negative values only if option is used */
    if (sign > 0) {
        for (i = 0; i < polygons->n_points; i++) {
            if (curvatures[i] < 0)
                curvatures[i] = 0;
        }
    }

    if (sign < 0) {
        for (i = 0; i < polygons->n_points; i++) {
            if (curvatures[i] > 0)
                curvatures[i] = 0;
        }
    }

    /* and smooth curvatures */
    if (smoothing)
        smooth_heatkernel(polygons, curvatures, fwhm);
  
    output_values_any_format(output_surface_file, polygons->n_points,
                 curvatures, TYPE_DOUBLE);

    delete_polygon_point_neighbours(polygons, n_neighbours,
                    neighbours, NULL, NULL);
    delete_object_list(n_objects, objects);
    free(curvatures);
  
    return(EXIT_SUCCESS);
}
