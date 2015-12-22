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

#ifdef _MSC_VER
  static const unsigned long __nan[2] = {0xffffffff, 0x7fffffff};
  #define FNAN (*(const float *) __nan)
#else
  #define FNAN 0.0f/0.0f
#endif

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
        char             *input_file, *output_surface_file, *values_file, *mask_file;
        int              i, n_objects, n_values, n_mask;
        int              *n_neighbours, **neighbours;
        File_formats     format;
        object_struct    **object_list;
        polygons_struct  *polygons;
        double           fwhm, *values, *mask;
        BOOLEAN          values_present, mask_present;

        initialize_argument_processing(argc, argv);

        if (!get_string_argument(NULL, &input_file) ||
            !get_string_argument(NULL, &output_surface_file) ||
            !get_real_argument(0.0, &fwhm)) {
                usage(argv[0]);
                exit(EXIT_FAILURE);
        }

        values_present = get_string_argument(NULL, &values_file);
        mask_present   = get_string_argument(NULL, &mask_file);

        if (input_graphics_any_format(input_file, &format, &n_objects,
                                      &object_list) != OK ||
            n_objects != 1 || get_object_type(object_list[0]) != POLYGONS) {
                fprintf(stderr, "Error reading %s.\n", input_file);
                exit(EXIT_FAILURE);
        }

        polygons = get_polygons_ptr(object_list[0]);

        if (values_present) {
                ALLOC(values, polygons->n_points);

                if (input_values_any_format(values_file, &n_values, &values) != OK) {
                        fprintf(stderr, "Cannot read values in %s.\n", values_file);
    	                exit(EXIT_FAILURE);
                }
                if (n_values != polygons->n_points) {
                        fprintf(stderr, "Number of values in %s differs from surface.\n", values_file);
    	                exit(EXIT_FAILURE);
                }
                
                /* if mask exists set all values to NaN where mask is zero */
                if (mask_present) {
                        ALLOC(mask, polygons->n_points);

                        if (input_values_any_format(mask_file, &n_mask, &mask) != OK) {
                                fprintf(stderr, "Cannot read values in %s.\n", mask_file);
    	                        exit(EXIT_FAILURE);
                        }
                        if (n_mask != polygons->n_points) {
                                fprintf(stderr, "Number of mask values in %s differs from surface.\n", mask_file);
    	                        exit(EXIT_FAILURE);
                        }
                        for (i=0; i<n_mask; i++) 
                                if(mask[i]==0) values[i] = FNAN;
                        FREE(mask);
                }

        } else values = NULL;

        smooth_heatkernel(polygons, values, fwhm);

        if (values_present) {
                output_values_any_format(output_surface_file, polygons->n_points,
                                         values, TYPE_DOUBLE);
                FREE(values);
        } else {

                compute_polygon_normals(polygons);

                if(output_graphics_any_format(output_surface_file, format, 1, 
                                object_list, NULL) != OK)
                        exit(EXIT_FAILURE);
        }

        return(EXIT_SUCCESS);
}
