/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 *
 */

#include <bicpl.h>

#include "CAT_Surf.h"
#include "CAT_SurfaceIO.h"
#include "CAT_Resample.h"

void
usage(char *executable)
{
        char *usage_str = "\n\
Usage: %s  surface_file sphere_file target_sphere_file resampled_output_surface_file [input_values output_values]\n\
Resamples a spherical inflated surface to an external defined sphere.\n\
\n\n";

        fprintf(stderr, usage_str, executable);
}

int
main(int argc, char *argv[])
{
        char             *surface_file, *sphere_file, *output_surface_file, *target_sphere_file;
        char             *input_values_file, *output_values_file;
        File_formats     format;
        int              n_objects;
        int              i, j, k;
        int              n_points, n_values;
        object_struct    **objects, **objects_src_sphere, **objects_target_sphere;
        polygons_struct  *polygons, *polygons_sphere, *target_sphere;
        double           *input_values = NULL, *output_values = NULL;
        BOOLEAN          values_specified;

        initialize_argument_processing(argc, argv);

        if (!get_string_argument(NULL, &surface_file    ) ||
            !get_string_argument(NULL, &sphere_file     ) ||
            !get_string_argument(NULL, &target_sphere_file     ) ||
            !get_string_argument(NULL, &output_surface_file)) {
                usage(argv[0]);
                exit(EXIT_FAILURE);
        }

        if (input_graphics_any_format(surface_file, &format,
                                      &n_objects, &objects) != OK ||
            n_objects != 1 || get_object_type(objects[0]) != POLYGONS) {
                fprintf(stderr, "File %s must contain 1 polygons object.\n",
                        surface_file);
                exit(EXIT_FAILURE);
        }

        if (strcmp(sphere_file,"NULL")) {
                if (input_graphics_any_format(sphere_file, &format,
                    &n_objects,
                    &objects_src_sphere) != OK ||
                    n_objects != 1 ||
                    get_object_type(objects_src_sphere[0]) != POLYGONS ) {
                        fprintf(stderr, "File %s must contain 1 polygons object.\n",
                        sphere_file);
                        exit(EXIT_FAILURE);
                }
                polygons_sphere = get_polygons_ptr(objects_src_sphere[0]);
        } else  polygons_sphere = NULL;

        if (input_graphics_any_format(target_sphere_file, &format,
                                      &n_objects,
                                      &objects_target_sphere) != OK ||
            n_objects != 1 ||
            get_object_type(objects_target_sphere[0]) != POLYGONS ) {
                fprintf(stderr, "File %s must contain 1 polygons object.\n",
                        target_sphere_file);
                exit(EXIT_FAILURE);
        }

        polygons        = get_polygons_ptr(objects[0]);
        target_sphere   = get_polygons_ptr(objects_target_sphere[0]);

        values_specified = get_string_argument(NULL, &input_values_file) &&
                           get_string_argument(NULL, &output_values_file);

        if (values_specified) {
                input_values  = (double *) malloc(sizeof(double) * polygons->n_points);
                output_values = (double *) malloc(sizeof(double) * target_sphere->n_points);

                if (input_values_any_format(input_values_file, &n_values, &input_values) != OK) {
                        fprintf(stderr, "Cannot read values in %s.\n", input_values_file);
    	               exit(EXIT_FAILURE);
                }
        }
        
        objects_target_sphere = resample_surface_to_target_sphere(polygons, polygons_sphere, 
                target_sphere, input_values, output_values);

        if(output_graphics_any_format(output_surface_file, format, 1,
                                   objects_target_sphere, NULL) != OK)
                    exit(EXIT_FAILURE);
    
        if (values_specified) {
                output_values_any_format(output_values_file,
                                         target_sphere->n_points,
                                         output_values, TYPE_DOUBLE);
                FREE(input_values);
                FREE(output_values);
        }
        
        return(EXIT_SUCCESS);
}
