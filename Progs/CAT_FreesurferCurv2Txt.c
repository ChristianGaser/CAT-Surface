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

void
usage(char *executable)
{
    char *usage_str = "\n\
Usage: %s  freesurfer_curv_file output_txt_file\n\n";

    fprintf(stderr, usage_str, executable);
}

int
main(int argc, char *argv[])
{
    char         *input_file, *output_surface_file;
    File_formats     format;
    FILE         *fp;
    int          n_values, i;
    double           *values;

    initialize_argument_processing(argc, argv);

    if (!get_string_argument(NULL, &input_file) ||
      !get_string_argument(NULL, &output_surface_file)) {
        usage(argv[0]);
        exit(EXIT_FAILURE);
    }

    if (input_freesurfer_curv( input_file, &n_values, &values) != OK) {
        printf("Error while reading %s.\n", input_file);
        exit(EXIT_FAILURE);
    }

    fp = fopen(output_surface_file, "w");
    if (!fp) exit(EXIT_FAILURE);

    for (i = 0; i < n_values; i++) {
        if ((fprintf(fp, " %g", values[i] ) <= 0) ||
          (fprintf(fp, "\n" ) <= 0))
            break;
    }

    fclose(fp);    
    return(EXIT_SUCCESS);
}
