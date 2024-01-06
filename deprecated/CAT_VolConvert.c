/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * most of the code is modified from
 * caret/caret_brain_set/BrainModelSurface.cxx.
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 *
 */

#include <bicpl.h>
#include <float.h>

#include "CAT_NiftiIO.h"

void
usage(char *executable)
{
        char *usage_str = "\n\
Usage: %s input_volume output_volume_file\n\n\
     Converts minc and nifti data.\n";

        fprintf(stderr, usage_str, executable);
}

int
main(int argc, char *argv[])
{
        char             *input_volume_file, *output_volume_file;
        File_formats     format;
        Volume           volume;        
    
        initialize_argument_processing(argc, argv);

        if (!get_string_argument(NULL, &input_volume_file) ||
            !get_string_argument(NULL, &output_volume_file)) {
                usage(argv[0]);
                exit(EXIT_FAILURE);
        }
    
        if (input_volume_all(input_volume_file, 3, File_order_dimension_names,
                            NC_UNSPECIFIED, FALSE, 0.0, 0.0,
                            TRUE, &volume, NULL) != OK)
                exit(EXIT_FAILURE);

        output_volume_all(output_volume_file, NC_FLOAT, 0, 0.0, 100000.0, volume, "CAT_VolConvert\n", NULL); 

        delete_volume( volume );

        return(EXIT_SUCCESS);
}
