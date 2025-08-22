/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 *
 */

#include <bicpl.h>

#include "CAT_Vol.h"
#include "CAT_Smooth.h"
#include "CAT_Curvature.h"
#include "CAT_SurfaceIO.h"

void
usage(char *executable)
{
    char *usage_str = "\n\
Usage: %s  surface_file output_values_file [curvtype] [fwhm] [use_abs_vals]\n\n\
   Calculate different curvature parameters (default: mean curvature\n\
   averaged over 3mm) from a given surface mesh file. If smoothing filter [fwhm]\n\
   is defined (in FWHM) a diffusion heat kernel will be applied. Optionally\n\
   the absolute value can ba calculated if use_abs_values is defined.\n\
   curvtype:\n\
        0 - mean curvature (averaged over 3mm, in degrees) (k1+k2)/2\n\
        1 - gaussian curvature k1*k2\n\
        2 - curvedness sqrt(0.5*(k1*k1+k2*k2))\n\
        3 - shape index atan((k1+k2)/(k2-k1))\n\
        4 - mean curvature (in radians) (k1+k2)/2\n\
        5 - sulcal depth like estimator\n\
        6 - bending energy k1*k1 + k2*k2\n\
        7 - sharpness (k1 - k2)^2\n\
        8 - folding index |k1|*(|k1| - |k2|)\n\
        9 - minimum curvature k2\n\
       10 - maximum curvature k1\n\
      >11 - depth potential with alpha = 1/curvtype (recommended value curvtype=650)\n\n";

    fprintf(stderr, usage_str, executable);
}

int
main(int argc, char *argv[])
{
    char *object_file, *output_surface_file;
    File_formats format;
    int i, n_objects, curvtype;
    int *n_neighbours, **neighbours, use_abs_values;
    object_struct **objects;
    polygons_struct *polygons;
    double *curvatures, distance, fwhm;
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

    /* limit range to values between -1..1 for some normalized curvtypes */
    if (((curvtype > 0) && (curvtype < 4)) || (curvtype == 10)) { 
        clip_data(curvatures, polygons->n_points, -1, 1, DT_FLOAT64);
    }
  
    /* use absolute value */
    if (use_abs_values) { 
        for (i = 0; i < polygons->n_points; i++)
            curvatures[i] = fabs(curvatures[i]);
    }

    /* and optionally smooth curvatures */
    if (smoothing)
        smooth_heatkernel(polygons, curvatures, fwhm);
  
    output_values_any_format(output_surface_file, polygons->n_points,
                 curvatures, TYPE_DOUBLE);

    delete_object_list(n_objects, objects);
    free(curvatures);
  
    return(EXIT_SUCCESS);
}
