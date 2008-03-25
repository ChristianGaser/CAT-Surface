/* Christian Gaser - christian.gaser@uni-jena.de                             */
/* Department of Psychiatry                                                  */
/* University of Jena                                                        */
/*                                                                           */
/* Copyright Christian Gaser, University of Jena.                            */

#include  <volume_io/internal_volume_io.h>
#include  <bicpl.h>
#include  <CAT_SurfaceIO.h>

private  void  usage(
    STRING   executable )
{
    STRING   usage_str = "\n\
Usage: %s  freesurfer_curv txt_file\n\n";

    print_error( usage_str, executable );
}

int  main(
    int   argc,
    char  *argv[] )
{
    STRING               input_filename, output_filename;
    File_formats         format;
    FILE                 *file;
    int                  n_values, i;
    Real				 *values;

    initialize_argument_processing( argc, argv );

    if( !get_string_argument( NULL, &input_filename ) ||
        !get_string_argument( NULL, &output_filename ))
    {
        usage( argv[0] );
        return( 1 );
    }

    if( input_freesurfer_curv( input_filename , &n_values, &values ) != OK )
    {
        print( "Error while reading %s.\n", input_filename );
        return( 1 );
    }

	if( open_file( output_filename, WRITE_FILE, ASCII_FORMAT, &file ) != OK )
		return( 1 );

	for (i=0; i<n_values; i++) {
		if( output_real( file, values[i] ) != OK || output_newline( file ) != OK )
			break;		    
    }
    
	(void) close_file( file );    

    return( 0 );
}
