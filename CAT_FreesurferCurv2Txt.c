/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 */

#include <volume_io/internal_volume_io.h>
#include <bicpl.h>

void
usage(char *executable)
{
        char *usage_str = "\n\
Usage: %s  freesurfer_curv txt_file\n\n";

        fprintf(stderr, usage_str, executable);
}

int
main(int argc, char *argv[])
{
        char                 *input_file, *output_file;
        File_formats         format;
        FILE                 *fp;
        int                  n_values, i;
        Real                 *values;

        initialize_argument_processing(argc, argv);

        if (!get_string_argument(NULL, &input_file) ||
            !get_string_argument(NULL, &output_file)) {
                usage(argv[0]);
                return(1);
        }

        if (input_freesurfer_curv( input_file, &n_values, &values) != OK) {
                printf("Error while reading %s.\n", input_file);
                return(1);
        }

        if (open_file(output_file, WRITE_FILE, ASCII_FORMAT, &fp) != OK)
                return(1);

        for (i = 0; i < n_values; i++) {
                if (output_real(fp, values[i]) != OK ||
                    output_newline(fp) != OK)
                        break;
        }

        close_file(fp);    
        return(0);
}
