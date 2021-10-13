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
Usage: %s  annot_file annot_out_file\n\
Text.\n\
\n\n";

        fprintf(stderr, usage_str, executable);
}

int  main(int   argc, char  *argv[])
{
        char             *sphere_file, *target_sphere_file, *annot_filename, *annot_out_file;
        File_formats     format;
        int              i, *array, n_labels, n_array;
        int              n_objects, n_values;
        char             **struct_names;
        double           *input_values;
        object_struct    **objects_src_sphere, **objects_target_sphere;
        polygons_struct  *polygons_sphere, *target_sphere;
        ATABLE           *atable;
    
        initialize_argument_processing( argc, argv );
    
        if (!get_string_argument(NULL, &annot_filename) ||
            !get_string_argument(NULL, &annot_out_file)) {
                usage(argv[0]);
                exit(EXIT_FAILURE);
        }

        read_annotation_table(annot_filename, &n_array, &array, &n_labels, &atable);

        printf("%d %d\n",n_labels, n_array);
        for (i = 0 ; i < n_labels ; i++) 
            printf("%d (%d %d %d) %s\n",atable[i].annotation,atable[i].r,atable[i].g,atable[i].b,atable[i].name);

        write_annotation_table(annot_out_file, n_array, array, n_labels, atable);
    
        return( 0 );
}
