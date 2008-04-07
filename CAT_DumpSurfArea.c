/* Christian Gaser - christian.gaser@uni-jena.de                             */
/* Department of Psychiatry                                                  */
/* University of Jena                                                        */
/*                                                                           */
/* Copyright Christian Gaser, University of Jena.                            */

#include  <volume_io/internal_volume_io.h>
#include  <bicpl.h>
#include  <CAT_Blur2d.h>
#include  <CAT_Surf.h>
    
int  main(
    int   argc,
    char  *argv[] )
{
    STRING               object_filename, output_filename;
    FILE                 *file;
    File_formats         format;
    int                  poly, n_objects, point_index, vertex_index, size;
    object_struct        **objects;
    polygons_struct      *polygons;
    float                 poly_size, area, surface_area;
    float                 *area_values;
    Point                points[MAX_POINTS_PER_POLYGON];
    BOOLEAN              all_values;
    Smallest_int         *point_done;

    initialize_argument_processing( argc, argv );

    if( !get_string_argument( NULL, &object_filename ))
    {
        print_error(
              "Usage: %s  object_file [output_file]\n",
              argv[0] );
        return( 1 );
    }

    if( input_graphics_any_format( object_filename,
                             &format, &n_objects, &objects ) != OK )
        return( 1 );

    if( get_string_argument( NULL, &output_filename ) )
    {
	all_values = TRUE;
	if( open_file( output_filename, WRITE_FILE, ASCII_FORMAT, &file ) != OK )
		return( 1 );
    }
    else all_values = FALSE;

    if( n_objects != 1 || get_object_type( objects[0] ) != POLYGONS )
    {
        print( "File must contain 1 polygons object.\n" );
        return( 1 );
    }

    polygons = get_polygons_ptr( objects[0] );
    
    ALLOC( area_values, polygons->n_points );
    surface_area = get_area_of_points( polygons, area_values );
    
    if( all_values )
    {
        for_less( point_index, 0, polygons->n_points )
        {
            if( output_real( file, area_values[point_index] ) != OK || output_newline( file ) != OK )
                break;		    
        }
	(void) close_file( file );    
    }

    print( "Total surface area: %g\n", surface_area );

    delete_object_list( n_objects, objects );

    return( 0 );
}
