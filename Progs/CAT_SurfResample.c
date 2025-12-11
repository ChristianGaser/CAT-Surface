/**
 * Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 *
 */

/**
 * @file CAT_SurfResample.c
 * 
 * This program resamples a spherical surface to a target sphere, optionally processing 
 * associated value files.
 * 
 * Usage:
 *    CAT_SurfResample surface_file|NULL sphere_file|NULL target_sphere_file 
 *              resampled_output_surface_file|NULL [input_values|input_annot output_values|output_annot|NULL]
 * 
 * The program allows for various inputs including the surface file, sphere file, target sphere file,
 * and the output surface file. It can also handle value files (input_values and output_values) or 
 * annotation files (input_annot and output_annot). The program supports label (categorical) interpolation
 * for handling labeled atlas data or annotation files.
 * 
 * Key Functionalities:
 * - Read and validate input files (surface, sphere, target sphere, value/annotation files).
 * - Resample the surface to the target sphere.
 * - Handle label interpolation for annotation files.
 * - Output the resampled surface and/or values to specified files.
 * 
 * 
 */

#include <bicpl.h>
#include <ParseArgv.h>

#include "CAT_Surf.h"
#include "CAT_SurfaceIO.h"
#include "CAT_Resample.h"

int label_interpolation = 0;
int areal_interpolation = 0;

static ArgvInfo argTable[] = {
  {"-label", ARGV_CONSTANT, (char *) TRUE, (char *) &label_interpolation,
    "Use label (categorical) interpolation (i.e. for labeled atlas data). This option is automatically used for annot (label) files."},
  {"-areal", ARGV_CONSTANT, (char *) TRUE, (char *) &areal_interpolation,
    "Use areal (sum-preserving) interpolation for scalar values."},
   {NULL, ARGV_END, NULL, NULL, NULL}
};

/**
 * Prints usage information for the program.
 * 
 * @param executable The name of the executable.
 */
void usage(char *executable)
{
    char *usage_str = 
    "\nUsage: %s surface_file|NULL sphere_file|NULL target_sphere_file resampled_output_surface_file|NULL [input_values|input_annot output_values|output_annot|NULL]\n\n\
    Resamples a spherical surface to an externally defined sphere. The program allows for various input and output combinations, handling both surface data and associated value or annotation files.\n\n\
    Arguments:\n\
      surface_file: Name of the input surface file or 'NULL' if not applicable.\n\
      sphere_file: Name of the input sphere file or 'NULL' if not applicable.\n\
      target_sphere_file: Name of the target sphere file.\n\
      resampled_output_surface_file: Name for the output resampled surface file or 'NULL' to skip saving the resampled surface.\n\
      input_values|input_annot: Optional name to a file containing input values or annotations.\n\
      output_values|output_annot: Optional name for outputting resampled values or annotations. 'NULL' to add values to the resampled surface.\n\n\
    Options:\n\
    -label: Enables label (categorical) interpolation, typically used for labeled atlas data or annot (label) files. This option is automatically used when working with annotation files.\n\n\
    -areal: Use sum-preserving areal interpolation for scalar values.\n\n\
    Example:\n\
    %s surface.gii sphere.gii target_sphere.gii output_surface.gii input_values.txt output_values.txt\n\n\
    Note:\n\
      At least one of surface_file or sphere_file must be defined.\n\
      If using annotation files, output values must be saved in a separate file.\n\
      The resampled values can be added directly to the resampled surface file if the output name for values is omitted.\n";

    fprintf(stderr, usage_str, executable, executable);
}

int main(int argc, char *argv[])
{
    char *surface_file, *sphere_file, *output_surface_file, *target_sphere_file;
    char *input_values_file, *output_values_file;
    File_formats format;
    int i, n_objects;
    int n_points, n_values, n_table, *out_annot, *in_annot;
    object_struct **objects, **objects_src_sphere, **objects_target_sphere;
    polygons_struct *polygons, *polygons_sphere, *target_sphere;
    double *input_values = NULL, *output_values = NULL;
    BOOLEAN values_defined, output_values_defined;
    ATABLE *atable;

    // Initialize argument processing and parse command line arguments
    initialize_argument_processing(argc, argv);

    // Parse arguments and handle errors if arguments are missing or invalid
    if (ParseArgv(&argc, argv, argTable, 0) ||
        !get_string_argument(NULL, &surface_file) ||
        !get_string_argument(NULL, &sphere_file) ||
        !get_string_argument(NULL, &target_sphere_file) ||
        !get_string_argument(NULL, &output_surface_file)) {
        usage(argv[0]);
        fprintf(stderr, "   %s -help\n\n", argv[0]);
        exit(EXIT_FAILURE);
    }

    // Determine if input and output value files are defined
    values_defined = get_string_argument(NULL, &input_values_file);
    output_values_defined = get_string_argument(NULL, &output_values_file);

    // Ensure at least one of the surface or sphere files is defined
    if ((!strcmp(sphere_file, "NULL")) && (!strcmp(surface_file, "NULL"))) {
        fprintf(stderr, "You have to define either surface or sphere.\n");
        exit(EXIT_FAILURE);
    }

    // Ensure output file is defined for annotation resampling
    if ((filename_extension_matches(input_values_file, "annot")) && !output_values_defined) {
        fprintf(stderr, "You have to define output for resampling of annot files.\n");
        exit(EXIT_FAILURE);
    }

    // Ensure at least one output (surface or values) is defined
    if ((!strcmp(output_surface_file, "NULL")) && values_defined && !output_values_defined) {
        fprintf(stderr, "You have to define either output surface or output values.\n");
        exit(EXIT_FAILURE);
    }

    // Ensure surface file is defined if resampled mesh is to be written
    if ((strcmp(output_surface_file, "NULL")) && (!strcmp(surface_file, "NULL"))) {
        fprintf(stderr, "You have to define surface if resampled mesh should be written.\n");
        exit(EXIT_FAILURE);
    }

    // Read input surface file, handling errors for invalid input
    if (strcmp(surface_file, "NULL")) {
        if (input_graphics_any_format(surface_file, &format,
                        &n_objects, &objects) != OK ||
            n_objects != 1 || get_object_type(objects[0]) != POLYGONS) {
            fprintf(stderr, "File %s must contain 1 polygons object.\n",
                    surface_file);
            exit(EXIT_FAILURE);
        }
    }

    // Read sphere file, handling errors for invalid input
    if (strcmp(sphere_file, "NULL")) {
        if (input_graphics_any_format(sphere_file, &format,
              &n_objects, &objects_src_sphere) != OK ||
              n_objects != 1 ||
              get_object_type(objects_src_sphere[0]) != POLYGONS) {
            fprintf(stderr, "File %s must contain 1 polygons object.\n",
                    sphere_file);
            exit(EXIT_FAILURE);
        }
        polygons_sphere = get_polygons_ptr(objects_src_sphere[0]);
    } else  polygons_sphere = NULL;

    // Read target sphere file, handling errors for invalid input
    if (input_graphics_any_format(target_sphere_file, &format,
                    &n_objects, &objects_target_sphere) != OK ||
        n_objects != 1 || get_object_type(objects_target_sphere[0]) != POLYGONS) {
        fprintf(stderr, "File %s must contain 1 polygons object.\n",
                target_sphere_file);
        exit(EXIT_FAILURE);
    }

    // Extract polygons from surface and target sphere
    if (strcmp(surface_file, "NULL")) {
        polygons = get_polygons_ptr(objects[0]);
    } else  polygons = NULL;
    target_sphere = get_polygons_ptr(objects_target_sphere[0]);

    // Process input values or annotation file
    if (values_defined) {
        // Allocate memory for input and output values
        input_values = (double *) malloc(sizeof(double) * polygons_sphere->n_points);
        output_values = (double *) malloc(sizeof(double) * target_sphere->n_points);

        // Handle annotation files, with a note on current platform limitation
        if (filename_extension_matches(input_values_file, "annot")) {
            /* this is not yet working under Windows */
            //fprintf(stderr, "Resampling of annotation files not yet working.\n");
            //exit(EXIT_FAILURE);

            label_interpolation = 1;
            read_annotation_table(input_values_file, &n_values, &in_annot, &n_table, &atable);
            for (i = 0 ; i < n_values ; i++)
                input_values[i] = (double)in_annot[i];
        } else {
            // Handle standard value files
            if (input_values_any_format(input_values_file, &n_values, &input_values) != OK) {
                fprintf(stderr, "Cannot read values in %s.\n", input_values_file);
                 exit(EXIT_FAILURE);
            }
        }
    }
    
    // Resample the surface to the target sphere
    objects_target_sphere = resample_surface_to_target_sphere(polygons, polygons_sphere, 
        target_sphere, input_values, output_values, label_interpolation, areal_interpolation);

    // Add resampled values to the surface or save surface without values
    if (values_defined && !output_values_defined) {
        if(output_graphics_any_format(output_surface_file, format, 1,
                   objects_target_sphere, output_values) != OK)
          exit(EXIT_FAILURE);
    } else {
        if (strcmp(output_surface_file, "NULL" )) {
            if(output_graphics_any_format(output_surface_file, format, 1,
                       objects_target_sphere, NULL) != OK)
              exit(EXIT_FAILURE);
        }
    }

    // Save resampled values as a separate file, handling different data types
    if (values_defined) {                
        if (output_values_defined) {
            // Convert output values to integers for annotation files and label interpolation
            if ((filename_extension_matches(output_values_file, "annot")) || (label_interpolation)) {
                out_annot = (int *) malloc(sizeof(int) * target_sphere->n_points);
                for (i = 0; i < target_sphere->n_points; i++) 
                    out_annot[i] = (int)round(output_values[i]);
            }
            
            // Write output values to annotation or value files
            if (filename_extension_matches(output_values_file, "annot")) {
                write_annotation_table(output_values_file, target_sphere->n_points, out_annot, n_table, atable);
            } else {
                if (label_interpolation)
                    output_values_any_format(output_values_file,
                         target_sphere->n_points, out_annot, TYPE_INTEGER);
                else
                    output_values_any_format(output_values_file,
                         target_sphere->n_points, output_values, TYPE_DOUBLE);
            }
            // Free allocated memory for output annotations if used
            if ((filename_extension_matches(output_values_file, "annot")) || (label_interpolation))
                FREE(out_annot);
        }
            
        // Free allocated memory for input and output values
        FREE(input_values);
        FREE(output_values);
    }   
    
    return(EXIT_SUCCESS);
}
