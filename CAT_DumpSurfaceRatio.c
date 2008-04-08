/* Christian Gaser - christian.gaser@uni-jena.de                             */
/* Department of Psychiatry                                                  */
/* University of Jena                                                        */
/*                                                                           */
/* Copyright Christian Gaser, University of Jena.                            */

#include  <volume_io/internal_volume_io.h>
#include  <bicpl.h>
#include  "CAT_Surf.h"
    
private  void  usage(
    STRING   executable )
{
    static  STRING  usage_str = "\n\
Usage: %s  surface.obj output_values.txt\n\
Computes surface ratio based on the method of Toro et al. 2008.\n\
\n\n";

    print_error( usage_str, executable );
}


int  main(
    int   argc,
    char  *argv[] )
{
    STRING               object_filename, output_filename;
    FILE                 *file;
    File_formats         format;
    int                  n_objects, point_index;
    object_struct        **objects;
    polygons_struct      *polygons;
    float                *area_values;
    Real                 radius;

    initialize_argument_processing( argc, argv );

    if( !get_string_argument( NULL, &object_filename ) ||
        !get_string_argument( NULL, &output_filename ) )
    {
        usage( argv[0] );
        return( 1 );
    }

    (void) get_real_argument( 20.0, &radius );

    if( input_graphics_any_format( object_filename,
                             &format, &n_objects, &objects ) != OK )
        return( 1 );

	if( open_file( output_filename, WRITE_FILE, ASCII_FORMAT, &file ) != OK )

    if( n_objects != 1 || get_object_type( objects[0] ) != POLYGONS )
    {
        print( "File must contain 1 polygons object.\n" );
        return( 1 );
    }

    polygons = get_polygons_ptr( objects[0] );
    
    ALLOC( area_values, polygons->n_points );
    area_values = get_surface_ratio( radius, polygons );
    
    for_less( point_index, 0, polygons->n_points )
    {
        if( output_real( file, area_values[point_index] ) != OK || output_newline( file ) != OK )
            break;		    
    }
	(void) close_file( file );    

    delete_object_list( n_objects, objects );

    return( 0 );
}
