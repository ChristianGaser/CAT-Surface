/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 *
 */

#include <bicpl.h>
#include <ParseArgv.h>

#include "CAT_Surf.h"
#include "CAT_SurfaceIO.h"
#include "CAT_Resample.h"

int label_interpolation = 0;

static ArgvInfo argTable[] = {
  {"-label", ARGV_CONSTANT, (char *) TRUE, (char *) &label_interpolation,
     "Use label (categorical) interpolation (i.e. for labeled atlas data). This option is automatically used for annot (label) files."},
   {NULL, ARGV_END, NULL, NULL, NULL}
};

void
usage(char *executable)
{
        char *usage_str = "\n\
Usage: %s  surface_file|NULL sphere_file|NULL target_sphere_file resampled_output_surface_file|NULL [input_values|input_annot output_values|output_annot]\n\
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
        int              i;
        int              n_points, n_values, n_table, *out_annot, *in_annot;
        object_struct    **objects, **objects_src_sphere, **objects_target_sphere;
        polygons_struct  *polygons, *polygons_sphere, *target_sphere;
        double           *input_values = NULL, *output_values = NULL;
        BOOLEAN          values_specified;
        ATABLE           *atable;

        initialize_argument_processing(argc, argv);

        if (ParseArgv(&argc, argv, argTable, 0) ||
            !get_string_argument(NULL, &surface_file) ||
            !get_string_argument(NULL, &sphere_file) ||
            !get_string_argument(NULL, &target_sphere_file) ||
            !get_string_argument(NULL, &output_surface_file)) {
                usage(argv[0]);
                fprintf(stderr, "     %s -help\n\n", argv[0]);
                exit(EXIT_FAILURE);
        }

        if (strcmp(surface_file,"NULL")) {
                if (input_graphics_any_format(surface_file, &format,
                                              &n_objects, &objects) != OK ||
                    n_objects != 1 || get_object_type(objects[0]) != POLYGONS) {
                        fprintf(stderr, "File %s must contain 1 polygons object.\n",
                                surface_file);
                        exit(EXIT_FAILURE);
                }
        }

        if ((!strcmp(sphere_file,"NULL")) && (!strcmp(surface_file,"NULL"))) {
          fprintf(stderr,"%s\n%s\n",sphere_file,surface_file);
                fprintf(stderr, "You have to define either surface or sphere.\n");
                exit(EXIT_FAILURE);
        }
        
        if ((strcmp(output_surface_file,"NULL")) && (!strcmp(surface_file,"NULL"))) {
                fprintf(stderr, "You have to define surface if resampled mesh should be written.\n");
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

        if (strcmp(surface_file,"NULL")) {
                polygons = get_polygons_ptr(objects[0]);
        } else  polygons = NULL;
        
        target_sphere   = get_polygons_ptr(objects_target_sphere[0]);

        values_specified = get_string_argument(NULL, &input_values_file) &&
                           get_string_argument(NULL, &output_values_file);

        if (values_specified) {                
                input_values  = (double *) malloc(sizeof(double) * polygons_sphere->n_points);
                output_values = (double *) malloc(sizeof(double) * target_sphere->n_points);

                if (filename_extension_matches(input_values_file, "annot")) {
                        /* this is not yet working under Windows */
                        fprintf(stderr, "Resampling of annotation files not yet working.\n");
                        exit(EXIT_FAILURE);

                        label_interpolation = 1;
                        read_annotation_table(input_values_file, &n_values, &in_annot, &n_table, &atable);
                        for (i = 0 ; i < n_values ; i++) 
                                input_values[i] = (double)in_annot[i];
                } else {
                    if (input_values_any_format(input_values_file, &n_values, &input_values) != OK) {
                            fprintf(stderr, "Cannot read values in %s.\n", input_values_file);
                           exit(EXIT_FAILURE);
                    }
                }
        }
        
        objects_target_sphere = resample_surface_to_target_sphere(polygons, polygons_sphere, 
                target_sphere, input_values, output_values, label_interpolation);

        /* skip writing resampled surface file if output name is NULL */
        if (strcmp(output_surface_file, "NULL" )) {
                if(output_graphics_any_format(output_surface_file, format, 1,
                                   objects_target_sphere, NULL) != OK)
                    exit(EXIT_FAILURE);
        }
    
        if (values_specified) {
                if ((filename_extension_matches(output_values_file, "annot")) || (label_interpolation)) {
                        out_annot  = (int *) malloc(sizeof(int) * target_sphere->n_points);
                        for (i = 0 ; i < target_sphere->n_points ; i++) 
                                out_annot[i] = (int)round(output_values[i]);
                }
                
                if (filename_extension_matches(output_values_file, "annot")) {
                        write_annotation_table(output_values_file, target_sphere->n_points, out_annot, n_table, atable);
                } else {
                        if (label_interpolation)
                                output_values_any_format(output_values_file,
                                         target_sphere->n_points,
                                         out_annot, TYPE_INTEGER);
                        else
                                output_values_any_format(output_values_file,
                                         target_sphere->n_points,
                                         output_values, TYPE_DOUBLE);
                }
                        
                if ((filename_extension_matches(output_values_file, "annot")) || (label_interpolation))
                        FREE(out_annot);
                FREE(input_values);
                FREE(output_values);
        }
        
        return(EXIT_SUCCESS);
}
