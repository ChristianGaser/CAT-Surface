/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id: CAT_AddValuesToSurf.c 261 2012-04-10 10:24:11Z gaser $
 *
 */

#include <bicpl.h>
#include <float.h>

#include "CAT_SurfaceIO.h"


void
usage(char *executable)
{
        char *usage_str = "\n\
Usage: %s  src.obj values_file dest.gii\n\n";

        fprintf(stderr, usage_str, executable);
}
    
int
main(int argc, char *argv[])
{
        double               value, *values;
        Status               status;
        char                 *src_file, *dest_file, *values_file;
        int                  i, p, n_objects, n_pts, n_values;
        Point                *pts;
        File_formats         format;
        object_struct        **object_list;

        initialize_argument_processing(argc, argv);

        if (!get_string_argument( NULL, &src_file) ||
            !get_string_argument( NULL, &values_file) ||
            !get_string_argument( NULL, &dest_file)) {
                usage(argv[0]);
                exit(EXIT_FAILURE);
        }

        if (filename_extension_matches(dest_file, "gii") != 1) {
                fprintf(stderr,"Only gifti output format allowed (use .gii as extension).\n");
                exit(EXIT_FAILURE);
        }
        
        if (input_graphics_any_format(src_file, &format, &n_objects,
                                      &object_list) != OK)
                exit(EXIT_FAILURE);

        if (input_values_any_format(values_file, &n_values, &values) != OK)
                exit(EXIT_FAILURE);
      
        if (n_objects > 1) {
                fprintf(stderr,"Only one object allowed.\n");
                exit(EXIT_FAILURE);
        }

        n_pts = get_object_points(object_list[1], &pts);
        if (n_pts != n_values) {
                fprintf(stderr,"Number of points differs from number of values.\n");
                exit(EXIT_FAILURE);
        }

        status = output_graphics_any_format(dest_file, format,
                                            n_objects, object_list, values);

        return(status != OK);
}
