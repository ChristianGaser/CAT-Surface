/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 *
 */

#include <bicpl.h>

#include "CAT_Smooth.h"
#include "CAT_SurfaceIO.h"
#include "CAT_Separate.h"

void
usage(char *executable) {
    char *usage_str = "\n\
Usage: %s  surface_file output_surface_file index \n\n\
   Separates polygons into its disjoint parts. Use a value of -1 for index to write the largest component.\n\n";

    fprintf(stderr, usage_str, executable);
}

int
main(int argc, char *argv[])
{
    char       *input_filename, *output_filename;
    char       out_file[512];
    int        n_objects, n_out;
    int        i, desired_index;
    File_formats   format;
    object_struct  **object_list, **object, *object2;
    polygons_struct  *polygons;

    initialize_argument_processing(argc, argv);

    if( !get_string_argument( NULL, &input_filename ) ||
      !get_string_argument( NULL, &output_filename ) )
    {
        usage(argv[0]);
        exit(EXIT_FAILURE);
    }

    get_int_argument(-1, &desired_index);

    if (input_graphics_any_format(input_filename, &format, &n_objects,
                    &object_list) != OK ||
      n_objects < 1 || get_object_type(object_list[0]) != POLYGONS) {
        fprintf(stderr, "File must have a polygons structure.\n");
        exit(EXIT_FAILURE);
    }

    polygons = get_polygons_ptr(object_list[0]);
    check_polygons_neighbours_computed(polygons);

    n_out = separate_polygons( polygons, desired_index, &object );

    if( n_out > 1) fprintf(stderr,"Extract largest of %d components.\n",n_out);
  
  
/*
Don't know whether this makes trouble or not, but sometimes this tool crashes for no reasons

    object2  = create_object( POLYGONS );
    triangulate_polygons( get_polygons_ptr(object[0]), get_polygons_ptr(object2) );  
    (void) output_graphics_any_format( output_filename, format, 1, &object2, NULL);
*/
    (void) output_graphics_any_format( output_filename, format, 1, &object[0], NULL);

    return(EXIT_SUCCESS);
}