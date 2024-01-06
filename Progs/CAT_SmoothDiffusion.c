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
#include <ParseArgv.h>

#include "CAT_Smooth.h"
#include "CAT_Curvature.h"
#include "CAT_SurfaceIO.h"
#include "CAT_Surf.h"

/* argument defaults */
double fwhm   = 12.0;
int  correct_mesh = 0;
char *values_file = NULL;
char *mask_file   = NULL;

/* the argument table */
static ArgvInfo argTable[] = {
  {"-fwhm", ARGV_FLOAT, (char *) TRUE, (char *) &fwhm,
     "Filter size in FWHM."},
  {"-values", ARGV_STRING, (char *) 1, (char *) &values_file, 
     "Optionally define values that will be smoothed instead of the surface"},
  {"-mask", ARGV_STRING, (char *) 1, (char *) &mask_file, 
     "Optional mask to limit smoothing of values within mask only and not smooth across the border"},
  {"-correct_mesh", ARGV_CONSTANT, (char *) TRUE, (char *) &correct_mesh,
     "Correct mesh in folded areas to compensate for the averaging effect in gyri and sulci.\n\
      Do not use this correction for smoothing sizes > 3 mm because it only works reliably for smaller smoothing."},
   {NULL, ARGV_END, NULL, NULL, NULL}
};

void
usage(char *executable)
{
        char *usage_str = "\n\
Usage: %s [options] surface_file output_file\n\n\
     Diffusion smoothing of values or surface points using heat kernel.\n\
     If values are defined with the values option then values will be\n\
     smoothed, otherwise only surface points.\n\n";
        fprintf(stderr, usage_str, executable);
}

int
main(int argc, char *argv[])
{
        char             *input_file, *output_surface_file;
        int              i, p, n_objects, n_values, n_mask;
        File_formats     format;
        object_struct    **object_list;
        polygons_struct  *polygons, *polygons_orig;
        double           distance, *values, *mask;
        BOOLEAN          values_present, mask_present;

        /* get the arguments from the command line */
        if (ParseArgv(&argc, argv, argTable, 0)) {
                usage(argv[0]);
                fprintf(stderr, "     %s -help\n\n", argv[0]);
                exit(EXIT_FAILURE);
        }

        initialize_argument_processing(argc, argv);

        if (!get_string_argument(NULL, &input_file) ||
            !get_string_argument(NULL, &output_surface_file)) {
                usage(argv[0]);
                fprintf(stderr, "     %s -help\n\n", argv[0]);
                exit(EXIT_FAILURE);
        }

        values_present = (values_file != NULL);
        mask_present   = (mask_file   != NULL);

        /* Don't use mesh correction for smoothing of values */
        if (values_present && correct_mesh > 0.0) {
                correct_mesh = 0.0;
                fprintf(stderr, "Correction for mesh folding will be disabled because \n\
                    it only works for mesh smoothing but not for smoothing of values.\n");
        }

        if (input_graphics_any_format(input_file, &format, &n_objects,
                                      &object_list) != OK ||
            n_objects != 1 || get_object_type(object_list[0]) != POLYGONS) {
                fprintf(stderr, "Error reading %s.\n", input_file);
                exit(EXIT_FAILURE);
        }

        polygons = get_polygons_ptr(object_list[0]);

        if (values_present) {

                if (input_values_any_format(values_file, &n_values, &values) != OK) {
                        fprintf(stderr, "Cannot read values in %s.\n", values_file);
                      exit(EXIT_FAILURE);
                }
                if (n_values != polygons->n_points) {
                        fprintf(stderr, "Number of values (n=%d) in %s differs from surface (n=%d).\n", n_values, values_file, polygons->n_points);
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

        /* save original mesh for comparison */
        if (correct_mesh) {
                polygons_orig = (polygons_struct *) malloc(sizeof(polygons_struct));
                copy_polygons(polygons, polygons_orig);
        }

        if (fwhm > 0.0)
                smooth_heatkernel(polygons, values, fwhm);

        /* Correct mesh in folded areas to compensate for the averaging effect in gyri and sulci.
           We use a folding measure (i.e. mean curvature averaged) to estimate the compensation. 
           The amount of compensation is automatically estimated using the difference between 
           the original and the smoothed surface. */
        if (correct_mesh) {
                correct_mesh_folding(polygons, polygons_orig, NULL, NULL, 0);
                free(polygons_orig);
        }

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
