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
#include  <float.h>

#define  BINTREE_FACTOR   0.5
    
private  void  usage(
    STRING   executable )
{
    STRING  usage_str = "\n\
Usage: %s surface.obj|none values_file|none output.pgm [nx ny] [curvature_distance]\n\n\
     Maps a surface to a flat sheet image. If values_file is not specified\n\
     then mean curvature is used as color. In this case you have to define surface_file.\n\
     Sheet image will be saved as PNG-file with dimensions nx x ny (default 256 x 256).\n\n";

    print_error( usage_str, executable );
}

int  main(
    int   argc,
    char  *argv[] )
{
    STRING               surface_filename, output_filename, values_filename;
    File_formats         format;
    polygons_struct      *polygons, unit_sphere;
    int                  i, n_objects, nx, ny, x, y, point_index;
    int                  poly, size, ind, n_items;
    int			         n_values, value_index;
    int                  *n_neighbours, **neighbours;
    Real                 value, u, v, curvature_distance;
    Real                 weights[1000], *values, *data, mn, mx;
    object_struct        **objects, *object;
    Point                unit_point, on_sphere_point, centre;
    Point                poly_points[1000];
    Vector               normal;
    progress_struct      progress;
    BOOLEAN              values_specified;

    /*--- get the arguments from the command line */

    initialize_argument_processing( argc, argv );

    if( !get_string_argument( NULL, &surface_filename ) ||
        !get_string_argument( NULL, &values_filename ) ||
        !get_string_argument( NULL, &output_filename ) )
    {
        usage( argv[0] );
        return( 1 );
    }

    get_int_argument( 256, &nx );
    get_int_argument( 256, &ny );
    get_real_argument( 3.0, &curvature_distance );

    if( equal_strings( surface_filename, "none" ) && 
        equal_strings( values_filename, "none" ) ) {
            print( "You have to either define surface_file or values_file.\n" );
            return( 1 );
    }

    if( equal_strings( values_filename, "none" ) ) {
        values_specified = FALSE;
    } else {
        if( input_texture_values( values_filename, &n_values, &values ) != OK )
            return( 1 );
        values_specified = TRUE;
    }

    if( equal_strings( surface_filename, "none" ) ) {
        /* if no surface_file is given then create tetra */
        object = create_object( POLYGONS );
        polygons = get_polygons_ptr(object);
        fill_Point( centre, 0.0, 0.0, 0.0 );
        create_tetrahedral_sphere( &centre, 1, 1, 1, 2*(n_values-2), polygons );
    } else {
        if( values_specified ) {
            if( input_texture_values( values_filename, &n_values, &values ) != OK )
                return( 1 );
        }
        if( input_graphics_any_format( surface_filename, &format, &n_objects, &objects ) != OK )
            return( 1 );
        /*--- check that the surface file contains a polyhedron */
        if( n_objects != 1 || get_object_type( objects[0] ) != POLYGONS ) {
            print( "Surface file must contain 1 polygons object.\n" );
            return( 1 );
        }
        /*--- get a pointer to the surface */
        polygons = get_polygons_ptr( objects[0] );
    }
    
    if( ! values_specified) {
        create_polygon_point_neighbours( polygons, TRUE, &n_neighbours,
                                     &neighbours, FALSE, NULL );
        ALLOC( values, polygons->n_points );
        get_polygon_vertex_curvatures( polygons, n_neighbours, neighbours,
                                   curvature_distance, 0.0, values );
    }
    
    /*--- create a unit sphere with same number of triangles as skin surface */
    fill_Point( centre, 0.0, 0.0, 0.0 );

    create_tetrahedral_sphere( &centre, 1.0, 1.0, 1.0,
                               polygons->n_items, &unit_sphere );

    create_polygons_bintree( &unit_sphere,
                             ROUND( (Real) unit_sphere.n_items *
                                    BINTREE_FACTOR ) );
    
    if( nx < 1 ) nx = 1;
    if( ny < 1 ) ny = 1;

    ALLOC(data, nx*ny);
    
    initialize_progress_report( &progress, FALSE, nx, "Mapping" );

    for_less( x, 0, nx ) {
        for_less( y, 0, ny ) {
        
            u = ((Real) x)/(Real) (nx - 1);
            v = ((Real) y)/(Real) (ny - 1);
  
            uv_to_point(u, v, &unit_point );
            
            poly = find_closest_polygon_point( &unit_point, &unit_sphere,
                                               &on_sphere_point );
            
            size = get_polygon_points( &unit_sphere, poly, poly_points );

            get_polygon_interpolation_weights( &on_sphere_point, size,
                                                   poly_points, weights );

            value = 0.0;
            for_less( i, 0, size ) {
                ind = unit_sphere.indices[
                         POINT_INDEX(unit_sphere.end_indices,poly,i)];
                value += weights[i] * values[ind];
            }
            value_index = x + (nx*y);
            data[value_index] = value;  
        }

        update_progress_report( &progress, x + 1 );
    }

	// scale data to uint8 range
	mn = FLT_MAX; mx = -FLT_MAX;
	for( i=0; i<nx*ny; i++) {
		if(data[i] > mx) mx = data[i];
		if(data[i] < mn) mn = data[i];
	}
	
	for( i=0; i<nx*ny; i++) 
		data[i] = 255.0*(data[i]-mn)/(mx - mn);
		
    terminate_progress_report( &progress );

    if( write_pgm(output_filename, data, nx, ny) != 0 )
        return( 1 );
    
    delete_polygons( &unit_sphere );
    delete_object_list( n_objects, objects );
    FREE( data );
    
    if( ! values_specified) FREE( values );
    
    return( 0 );
}
