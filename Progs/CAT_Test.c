/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 *
 */

#include "CAT_SurfaceIO.h"

void
usage(char *executable)
{
        char *usage_str = "\n\
Usage: %s  sphere_file target_sphere_file annot_file input_values_file\n\
Text.\n\
\n\n";

        fprintf(stderr, usage_str, executable);
}

int  main(int   argc, char  *argv[])
{
        char             *sphere_file, *target_sphere_file, *annot_filename, *input_values_file;
        File_formats     format;
        int              i, *array, n_table, n_array;
        int              n_objects, n_values;
        char             **struct_names;
        double           *input_values;
        object_struct    **objects_src_sphere, **objects_target_sphere;
        polygons_struct  *polygons_sphere, *target_sphere;
        ATABLE           *atable;
    
        initialize_argument_processing( argc, argv );
    
        if (!get_string_argument(NULL, &sphere_file) ||
            !get_string_argument(NULL, &target_sphere_file) ||
            !get_string_argument(NULL, &annot_filename) ||
            !get_string_argument(NULL, &input_values_file)) {
                usage(argv[0]);
                exit(EXIT_FAILURE);
        }

        if (input_graphics_any_format(sphere_file, &format,
                                      &n_objects, &objects_src_sphere) != OK ||
            n_objects != 1 || get_object_type(objects_src_sphere[0]) != POLYGONS) {
                fprintf(stderr, "File %s must contain 1 polygons object.\n",
                        sphere_file);
                exit(EXIT_FAILURE);
        }

        if (input_graphics_any_format(target_sphere_file, &format,
                                      &n_objects,
                                      &objects_target_sphere) != OK ||
            n_objects != 1 || get_object_type(objects_target_sphere[0]) != POLYGONS ) {
                fprintf(stderr, "File %s must contain 1 polygons object.\n",
                        target_sphere_file);
                exit(EXIT_FAILURE);
        }

        target_sphere   = get_polygons_ptr(objects_target_sphere[0]);
        polygons_sphere = get_polygons_ptr(objects_src_sphere[0]);

        input_values  = (double *) malloc(sizeof(double) * polygons_sphere->n_points);

        if (input_values_any_format(input_values_file, &n_values, &input_values) != OK) {
                fprintf(stderr, "Cannot read values in %s.\n", input_values_file);
               exit(EXIT_FAILURE);
        }

        read_annotation_table(annot_filename, &n_array, &array, &n_table, &atable);
        write_annotation_table("test.annot", n_array, array, n_table, atable);

        if(target_sphere->n_points != n_array) {
                    fprintf(stderr,"Annotation file and target sphere have different dimensions.\n");
                    exit(EXIT_FAILURE);
        }

        for (i = 0 ; i < n_table ; i++) 
            printf("%d %s\n",atable[i].annotation,atable[i].name);
    
    /*    for (i = 0 ; i < n ; i++) 
    printf("%d %d\n",i,array[i]);
      */  
        return( 0 );
}
