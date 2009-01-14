/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 *
 */

#include <volume_io/internal_volume_io.h>
#include <bicpl.h>

void
usage(char *executable)
{
        char *usage_str = "\n\
Usage: %s  input.ext output.ext\n\n\
    The surface format will be recognized by the file extension:\n\
       .obj - BIC object\n\
       .off - Geomview OOGL\n\
    or otherwise Freesurfer\n\
         \n\n";

        fprintf(stderr, usage_str, executable);
}

int
main(int argc, char *argv[])
{
        char                 *input_file, *output_file;
        File_formats         format;
        int                  status, n_objects;
        object_struct        **objects;
        polygons_struct      *polygons;
        signed char          *done_flags;

        initialize_argument_processing(argc, argv);

        if (!get_string_argument(NULL, &input_file) ||
            !get_string_argument(NULL, &output_file)) {
                usage(argv[0]);
                return(1);
        }

        if (input_graphics_any_format(input_file, &format,
                                      &n_objects, &objects) != OK)
                return(1);

        if (n_objects != 1 || get_object_type(objects[0]) != POLYGONS) {
                printf("File must contain 1 polygons object.\n");
                return(1);
        }

        if (output_graphics_any_format(output_file, format,
                                       n_objects, objects) != OK)
                return(1);

        delete_object_list(n_objects, objects);

        return(0);
}
