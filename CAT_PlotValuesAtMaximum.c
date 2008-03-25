/* Christian Gaser - christian.gaser@uni-jena.de                             */
/* Department of Psychiatry                                                  */
/* University of Jena                                                        */
/*                                                                           */
/* Copyright Christian Gaser, University of Jena.                            */

#include  <volume_io/internal_volume_io.h>
#include  <bicpl.h>

private  void  usage(
    STRING   executable )
{
    STRING   usage_str = "\n\
Usage: %s values_ref.txt output_values.txt values1.txt [values2.txt .. valuesn.txt]\n\n\
     Plot values at maximum of reference for each file.\n\n";

    print_error( usage_str, executable );
}

int  main(
    int   argc,
    char  *argv[] )
{
    STRING               values_filename, input_filename, output_filename;
    FILE                 *file, *out_file;
    File_formats         format;
    int                  n_files, n_values, i, j, max_index;
    Real		         output_value, *values, max_value;

    initialize_argument_processing( argc, argv );

    if( !get_string_argument( NULL, &values_filename ) ||
        !get_string_argument( NULL, &output_filename ))
    {
        usage( argv[0] );
        return( 1 );
    }

    if( input_texture_values( values_filename, &n_values, &values ) != OK ) {
        print_error( "Cannot read values.\n" );
    	return( 1 );
    }

    n_files = argc - 3;
        
    // find index of maximum
    max_value = -1e10;
    for_less( i, 0, n_values )
    {
		if(values[i] > max_value) {
			max_value = values[i];
			max_index = i;
		}
    }
    
    fprintf(stderr,"Maximum value of %3.2f found at line %d.\n",max_value,max_index);
    
    if( open_file( output_filename, WRITE_FILE, ASCII_FORMAT, &out_file ) != OK)
    	return( 1 );

    for_less(i, 0, n_files)
    {
        get_string_argument( NULL, &input_filename );
        if((file = fopen(input_filename, "r")) == 0) {
            fprintf(stderr, "Couldn't open file %s.\n", input_filename);
            return(1);
        }
        for_less(j, 0, max_index+1) {
            if( input_real( file, &values[j] ) != OK )
                return( 1 );
        }
        if( output_real( out_file, values[max_index] ) != OK ||
        	output_newline( out_file ) != OK )
                return( 1 );
        (void) close_file( file );
    }
    (void) close_file( out_file );

    return( 0 );
}
