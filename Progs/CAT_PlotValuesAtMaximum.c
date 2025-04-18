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
Usage: %s values_ref.txt output_values.txt values1.txt [values2.txt .. valuesn.txt]\n\n\
     Plot values at maximum of reference for each file.\n\n";

    fprintf(stderr, usage_str, executable);
}

int
main(int argc, char *argv[])
{
    char         *values_file, *input_file, *output_surface_file;
    FILE         *infp, *outfp;
    File_formats     format;
    int          n_files, n_values, i, j, max_index;
    double           output_value, *values, max_value;

    initialize_argument_processing(argc, argv);

    if (!get_string_argument(NULL, &values_file) ||
      !get_string_argument(NULL, &output_surface_file)) {
        usage(argv[0]);
        exit(EXIT_FAILURE);
    }

    if (input_values_any_format(values_file, &n_values, &values) != OK) {
        fprintf(stderr, "Cannot read values in %s.\n", values_file);
        exit(EXIT_FAILURE);
    }

    n_files = argc - 3;
    
    // find index of maximum
    max_value = -1e10;
    for (i = 0; i < n_values; i++) {
        if (values[i] > max_value) {
            max_value = values[i];
            max_index = i;
        }
    }
  
    fprintf(stderr, "Maximum value of %3.2f found at line %d.\n",
        max_value, max_index);
  
    outfp = fopen(output_surface_file, "w");
    if (!outfp) exit(EXIT_FAILURE);

    for (i = 0; i < n_files; i++) {
        get_string_argument(NULL, &input_file);
        if (input_values_any_format(input_file, &n_values, &values) != OK) {
            fprintf(stderr, "Cannot read values in %s.\n", values_file);
            exit(EXIT_FAILURE);
        }
        if ((fprintf(outfp, " %g", values[max_index] ) <= 0) ||
            (fprintf(outfp, "\n" ) <= 0))
            exit(EXIT_FAILURE);

    }

    fclose(outfp);
    return(EXIT_SUCCESS);
}
