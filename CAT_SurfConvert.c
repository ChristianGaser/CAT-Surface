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
Usage: %s  input.ext output.ext\n\
\n\
The surface format will be recognized by the file extension:\n\
   .obj - BIC object\n\
   .off - Geomview OOGL\n\
    or otherwise Freesurfer\n\
         \n\n";

    print_error( usage_str, executable );
}

int  main(
    int   argc,
    char  *argv[] )
{
    STRING               input_filename, output_filename;
    File_formats         format;
    int                  status, n_objects;
    object_struct        **objects;
    polygons_struct      *polygons;
    Smallest_int	     *done_flags;

    initialize_argument_processing( argc, argv );

    if( !get_string_argument( NULL, &input_filename ) ||
        !get_string_argument( NULL, &output_filename ))
    {
        usage( argv[0] );
        return( 1 );
    }

    if( input_graphics_any_format( input_filename,
                             &format, &n_objects, &objects ) != OK )
        return( 1 );

    if( n_objects != 1 || get_object_type( objects[0] ) != POLYGONS )
    {
        print( "File must contain 1 polygons object.\n" );
        return( 1 );
    }
    
    if( output_graphics_any_format( output_filename, format,
                                   n_objects, objects ) != OK)
        return(1);

    delete_object_list( n_objects, objects );
    
    return( 0 );
}
