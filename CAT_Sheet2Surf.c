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
    polygons_struct      *polygons, unit_sphere;
    int                  i, n_objects, nx, ny, x, y, point_index;
    int                  degree, poly, size, ind, n_done, interpolate;
    int                  *n_neighbours, **neighbours;
    double		         *image;
    object_struct        **objects;
    BOOLEAN              use_volume;
    Point                unit_point, on_sphere_point, centre;
    Point                poly_points[1000], centroid;
    double               u, v, *values;
    Vector               normal;
    progress_struct      progress;
    double               tmp_x, tmp_y;

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
    
    create_polygon_point_neighbours( polygons, TRUE, &n_neighbours,
                                     &neighbours, FALSE, NULL );

    /*--- create a unit sphere with same number of triangles as skin surface */
    fill_Point( centre, 0.0, 0.0, 0.0 );

    create_tetrahedral_sphere( &centre, 1.0, 1.0, 1.0,
                               polygons->n_items, &unit_sphere );

    create_polygons_bintree( &unit_sphere,
                             ROUND( (double) unit_sphere.n_items *
                                    BINTREE_FACTOR ) );

    ALLOC(values, polygons->n_points);

    if ((image = read_pgm( input_filename, &nx, &ny )) == NULL )
        return( 1 );

    initialize_progress_report( &progress, FALSE, polygons->n_points, "Mapping" );

    if( open_file( output_filename, WRITE_FILE, ASCII_FORMAT, &outputfile ) != OK )
        return( 1 );

	tmp_x = (double)nx - 1.0;
	tmp_y = (double)ny - 1.0;

    for ( i=0; i<polygons->n_points; i++ ) {

        point_to_uv(&unit_sphere.points[i], &u, &v);
		x = (int) round(u*tmp_x - 1.0);
		y = (int) round(v*tmp_y - 1.0);
		if (interpolate) {
			double xp = u*tmp_x - 1.0 - x;
			double yp = v*tmp_y - 1.0 - y;
			double xm = 1.0 - xp;
			double ym = 1.0 - yp;
			double H00 = image[x  +    y*nx];
			double H01 = image[x  +(y+1)*nx];
			double H10 = image[x+1+    y*nx];
			double H11 = image[x+1+(y+1)*nx];

			values[i] = (ym * ( xm * H00 + xp * H10) + 
		 		yp * ( xm * H01 + xp * H11));
		} else values[i] = image[x + y*nx];
		
		if( output_real( outputfile, values[i] ) != OK || output_newline( outputfile ) != OK )
            break;

        update_progress_report( &progress, i + 1 );
    }
    
    (void) close_file( outputfile );

    terminate_progress_report( &progress );

    delete_polygons( &unit_sphere );
    delete_object_list( n_objects, objects );

    FREE( image );
    FREE( values );

    return( 0 );
}
