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

#include "CAT_SurfaceIO.h"

int
main(int argc, char *argv[])
{
        FILE             *fp;
        char             *file, *output_file;
        int              counter, i, j, n_values, prev_n_values, n_plus;
        int              n_minus, n_infiles;
        Real             *values, *result;
        BOOLEAN          plus = 0, minus = 0;
        char             **infiles;

        initialize_argument_processing(argc, argv);

        if (!get_string_argument(NULL, &output_file) || (argc < 3)) {
                printf("Usage: %s  output.txt file1 +|- file2 +|- file3 ..\n",
                      argv[0]);
                exit(EXIT_FAILURE);
        }

        counter = 0;
        n_plus = 1;
        n_minus = 0;
        n_infiles = 0;

        infiles = &argv[1];

        while (get_string_argument(NULL, &file)) {
                n_infiles++;
                if (equal_strings(file, "+"))
                        n_plus++;
                else if (equal_strings(file, "-"))
                        n_minus++;
        }

        for (i = 1; i <= n_infiles; i++) {
                if (equal_strings(infiles[i], "+" )) {
                        plus  = 1;
                        minus = 0;
                } else if (equal_strings(infiles[i], "-" )) {
                        plus  = 0;
                        minus = 1;
                } else {
                        if (input_values_any_format(infiles[i], &n_values,
                                                 &values ) != OK) {
                                fprintf(stderr, "Error reading file %s\n",
                                        infiles[i]);
                                exit(EXIT_FAILURE);
                        }
                        if (counter == 0) {
                                ALLOC(result, n_values);

                                for (j=0; j < n_values; j++) {
                                        if (plus)
                                                result[j] = values[j]/n_plus;
                                        else if (minus)
                                                result[j] = - values[j]/n_minus;
                                        else
                                                result[j] = values[j]/n_plus;
                                }

                        } else {
                                if (prev_n_values != n_values) {
                                        fprintf(stderr, "Wrong number of values in %s: %d instead of %d\n", file, n_values, prev_n_values);
                                        exit(EXIT_FAILURE);
                                }
                                for (j = 0; j < n_values; j++) {
                                        if (plus)
                                                result[j] += values[j]/n_plus;
                                        else if (minus)
                                                result[j] -= values[j]/n_minus;
                                        else
                                                result[j] += values[j]/n_plus;
                                }
                        }
                        prev_n_values = n_values;
                        counter++;
                }
        }

        output_values_any_format(output_file, n_values, result, TYPE_REAL);

        FREE(result);
        return(EXIT_SUCCESS);
}
