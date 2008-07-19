/* Christian Gaser - christian.gaser@uni-jena.de                             */
/* Department of Psychiatry                                                  */
/* University of Jena                                                        */
/*                                                                           */
/* Copyright Christian Gaser, University of Jena.                            */

#include  <volume_io/internal_volume_io.h>
#include  <bicpl.h>
#include  <CAT_SheetIO.h>
#include  <CAT_Map2d.h>
#include  <CAT_Blur2d.h>

#define  BINTREE_FACTOR   0.5
    
private  void  usage(
    STRING   executable )
{
    STRING  usage_str = "\n\
Usage: %s surface.obj  input.pgm output.txt [0|1]\n\n\
     Maps an PGM-image to a surface. The output is saved as texture values.\n\n\
     Optionally linear interpolation can be skipped with a 0 as 4th argument\n\n";

    print_error( usage_str, executable );
}

int  main(
    int   argc,
    char  *argv[] )
{
    STRING               surface_filename, output_filename, input_filename;
    FILE		         *outputfile;
    File_formats         format;
    polygons_struct      *polygons;
    int                  i, n_objects;
    int                  degree, poly, size_map[2], ind, n_done, interpolate;
    double		         *image;
    object_struct        **objects;
    double               *values;

    /*--- get the arguments from the command line */

    initialize_argument_processing( argc, argv );

    if( !get_string_argument( NULL, &surface_filename ) ||
        !get_string_argument( NULL, &input_filename ) ||
        !get_string_argument( NULL, &output_filename ) )
    {
        usage( argv[0] );
        return( 1 );
    }

	get_int_argument(1, &interpolate);
    /*--- input the surface */
    if( input_graphics_any_format( surface_filename, &format, &n_objects, &objects )
         != OK )
        return( 1 );

    /*--- check that the surface file contains a polyhedron */
    if( n_objects != 1 || get_object_type( objects[0] ) != POLYGONS )
    {
        print( "Surface file must contain 1 polygons object.\n" );
        return( 1 );
    }

    /*--- get a pointer to the surface */
    polygons = get_polygons_ptr( objects[0] );
    
    values = (double *)malloc(sizeof(double)*polygons->n_points);

    if ((image = read_pgm( input_filename, &size_map[0], &size_map[1] )) == NULL )
        return( 1 );
        
    map_sheet2d_to_sphere(image, values, polygons, interpolate, size_map);

    if( open_file( output_filename, WRITE_FILE, ASCII_FORMAT, &outputfile ) != OK )
        return( 1 );
    
    for ( i=0; i<polygons->n_points; i++ ) {
		if( output_real( outputfile, values[i] ) != OK || output_newline( outputfile ) != OK )
            break;
    }
    (void) close_file( outputfile );

    delete_object_list( n_objects, objects );

    free( image );
    free( values );

    return( 0 );
}
