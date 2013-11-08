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

void find_conformal_map(polygons_struct *polygons);

void
usage(char *executable)
{
        char *usage_str = "\n\
Usage: %s input.obj output.obj\n\n\
     Generate the Laplace-Beltrami conformal map of the mesh specified in infile.\n\n\
     Results are saved in outfile.\n";

        fprintf(stderr, usage_str, executable);
}

int
main(int argc, char** argv)
{
        char *input_file, *output_file;
        object_struct **objects;
        int n_objects;
        File_formats format;
        polygons_struct *polygons;
    
        initialize_argument_processing(argc, argv);
        if (!get_string_argument(NULL, &input_file) ||
            !get_string_argument(NULL, &output_file)) {
                usage(argv[0]);
                exit(EXIT_FAILURE);
        }

        if (input_graphics_any_format(input_file, &format,
                                      &n_objects, &objects) != OK) {
                printf("Error reading input file\n");
                exit(EXIT_FAILURE);
        }

        if (n_objects != 1 || get_object_type(objects[0]) != POLYGONS) {
                printf("Input file must contain one polygon object.\n");
                exit(EXIT_FAILURE);
        }

        polygons = get_polygons_ptr(objects[0]);

        find_conformal_map(polygons);

        if(output_graphics_any_format(output_file, format, 1, 
                        objects, NULL) != OK)
                    exit(EXIT_FAILURE);
        exit(EXIT_SUCCESS);
}
