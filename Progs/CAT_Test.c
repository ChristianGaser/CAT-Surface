/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 *
 */


#include <bicpl.h>
#include "CAT_SurfaceIO.h"
#include "CAT_Deform.h"
#include "CAT_Curvature.h"
#include "CAT_Smooth.h"
#include "CAT_Intersect.h"
#include "CAT_Defect.h"

void
usage(char *executable)
{
    char *usage_str = "\n\
Usage: %s surface_file output_file fwhm [values_file] [mask_file]\n\n\
   Diffusion smoothing of values or surface points using\n\
   heat kernel. If values are defined then values will be\n\
   smoothed, otherwise only surface points.\n\n";

    fprintf(stderr, usage_str, executable);
}

int
main(int argc, char *argv[])
{
    char       *input_file, *values_file;
    int        i, n_objects, n_values;
    int        *n_neighbours, **neighbours;
    File_formats   format;
    object_struct  **object_list;
    polygons_struct  *polygons;

    initialize_argument_processing(argc, argv);

    if (!get_string_argument(NULL, &input_file) ||
      !get_string_argument(NULL, &values_file)) {
        usage(argv[0]);
        exit(EXIT_FAILURE);
    }

    if (input_graphics_any_format(input_file, &format, &n_objects,
                    &object_list) != OK ||
      n_objects != 1 || get_object_type(object_list[0]) != POLYGONS) {
        fprintf(stderr, "Error reading %s.\n", input_file);
        exit(EXIT_FAILURE);
    }

    polygons = get_polygons_ptr(object_list[0]);

    compute_polygon_normals(polygons);
    int n_self_hits;
    int *defects = find_near_self_intersections(polygons, 0.95, &n_self_hits);

    output_values_any_format(values_file, polygons->n_points, defects, TYPE_INTEGER);



    free(defects);

    return(EXIT_SUCCESS);
}
