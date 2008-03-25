/* Christian Gaser - christian.gaser@uni-jena.de                             */
/* Department of Psychiatry                                                  */
/* University of Jena                                                        */
/*                                                                           */
/* Copyright Christian Gaser, University of Jena.                            */

#include  <volume_io/internal_volume_io.h>
#include  <bicpl.h>

int  main(
    int    argc,
    char   *argv[] )
{
    FILE             *file;
    STRING           filename, output_filename;
    int              counter, i, j, n_values, prev_n_values, n_plus;
    int              n_minus, n_infiles;
    Real             *values, *result;
    BOOLEAN          plus = 0, minus = 0;
    char             **infiles;

    initialize_argument_processing( argc, argv );

    if( !get_string_argument( NULL, &output_filename ) || (argc < 3)) {
        print( "Usage: %s  output.txt file1 +|- file2 +|- file3 ..\n",
               argv[0] );
        return( 1 );
    }

    counter = 0;
    n_plus = 1;
    n_minus = 0;
    n_infiles = 0;
    
    infiles = &argv[1];

    while (get_string_argument( NULL, &filename)) {
        n_infiles++;
        if( equal_strings( filename, "+" ) ) n_plus++;
        else if( equal_strings( filename, "-" ) ) n_minus++;
    } /* while */

    for (i=1; i <= (n_infiles); i++) {
        if( equal_strings( infiles[i], "+" )) {
			plus  = 1;
			minus = 0;
        } else if( equal_strings( infiles[i], "-" )) {
			plus  = 0;
			minus = 1;
	} else {
		if( input_texture_values( infiles[i], &n_values, &values ) != OK ) {
                print_error( "Error reading file %s\n", infiles[i]);
                return( 1 );
		}
	    if( counter == 0) {
	        ALLOC(result, n_values);
	        for (j=0; j < n_values; j++)
	            if( plus ) result[j] = values[j]/n_plus;
		    else if( minus) result[j] = - values[j]/n_minus;
		    else result[j] = values[j]/n_plus;
	    } else {
	        if( prev_n_values != n_values) {
	            print_error("Wrong number of values in %s: %d instead of %d\n",filename, n_values, prev_n_values);
				return( 1 );
			}
	        for (j=0; j < n_values; j++)
	            if( plus ) result[j] += values[j]/n_plus;
				else if( minus) result[j] -= values[j]/n_minus;
				else result[j] += values[j]/n_plus;
            }	
	    prev_n_values = n_values;
	    counter++;
		}
    }

    if( open_file( output_filename, WRITE_FILE, ASCII_FORMAT, &file ) != OK )
        return( 1 );
    for (i=0; i < n_values; i++) {
        if( output_real( file, result[i] ) != OK ||
            output_newline( file ) != OK )
            break;
    }
    (void) close_file( file );	    
    
    FREE(result);
}
