/* Christian Gaser - christian.gaser@uni-jena.de                             */
/* Department of Psychiatry                                                  */
/* University of Jena                                                        */
/*                                                                           */
/* Copyright Christian Gaser, University of Jena.                            */

#include  <volume_io/internal_volume_io.h>
#include  <bicpl.h>
#include "Cat_Surf.h"

private  void  usage(
    STRING   executable )
{
    static  STRING  usage_str = "\n\
Usage: %s  surface.obj sphere.obj output.obj n_points [input_values.txt output_values.txt]\n\
Resamples a spherical inflated surface to a sphere.\n\
\n\n";

    print_error( usage_str, executable );
}

int  main(
    int   argc,
    char  *argv[] )
{
    STRING               surface_filename, sphere_filename, output_filename;
    STRING               input_values_filename, output_values_filename;
    File_formats         format;
    FILE                 *file;
    int                  n_objects, n_objects_src_sphere, n_triangles, point, i, j, k;
    int                  poly, n_points, n_intersections;
    int                  *n_neighbours, **neighbours;
    Point                centre, point_on_src_sphere, poly_points[MAX_POINTS_PER_POLYGON];
    Point				 poly_points_src[MAX_POINTS_PER_POLYGON], *new_points;
    Point                scaled_point;
    object_struct        **objects, **objects_src_sphere, *objects_dest_sphere;
    polygons_struct      *polygons_fiducial, *polygons_src_sphere, *polygons_dest_sphere;
    Real                 *input_values, *output_values, dist;
	Vector 				 normal;
    BOOLEAN				 values_specified;
    Real                 weights[MAX_POINTS_PER_POLYGON];

    initialize_argument_processing( argc, argv );

    if( !get_string_argument( NULL, &surface_filename ) ||
        !get_string_argument( NULL, &sphere_filename ) ||
        !get_string_argument( NULL, &output_filename ) )
    {
        usage( argv[0] );
        return( 1 );
    }

    (void) get_int_argument( 81920, &n_triangles );

	// check tetrahedral topology
	// best areal distribution of triangles is achieved for 20 edges
	int n_polygons;
	n_polygons = n_triangles;
	
    while( n_polygons != 20 && n_polygons > 8 && n_polygons % 4 == 0 )
         n_polygons /= 4;

	if( n_polygons != 20 ) {
		fprintf(stderr,"Warning: Number of triangles %d is not recommend because\n", n_triangles);
		fprintf(stderr,"tetrahedral topology is not optimal.\n");
		fprintf(stderr,"Please try 20*(4*x) triangles (e.g. 81920).\n");
	}
	
    if( input_graphics_any_format( surface_filename,
                             &format, &n_objects, &objects ) != OK ||
        n_objects != 1 || get_object_type(objects[0]) != POLYGONS ) {
        print_error( "File %s must contain 1 polygons object.\n",
                     surface_filename );
        return( 1 );
    }

    if( input_graphics_any_format( sphere_filename,
                             &format, &n_objects_src_sphere, &objects_src_sphere ) != OK ||
        n_objects_src_sphere != 1 || get_object_type(objects_src_sphere[0]) != POLYGONS ) {
        print_error( "File %s must contain 1 polygons object.\n",
                     sphere_filename );
        return( 1 );
    }

    polygons_fiducial = get_polygons_ptr( objects[0] );
    polygons_src_sphere = get_polygons_ptr( objects_src_sphere[0] );

    values_specified = get_string_argument( NULL, &input_values_filename ) &&
                       get_string_argument( NULL, &output_values_filename );

    if( values_specified )
    {
        ALLOC( input_values, polygons_fiducial->n_points );
        ALLOC( output_values, polygons_src_sphere->n_points );

        if( open_file( input_values_filename, READ_FILE, ASCII_FORMAT, &file ) != OK )
            return( 1 );

        for_less( i, 0, polygons_fiducial->n_points )
        {
            if( input_real( file, &input_values[i] ) != OK )
                return( 1 );
        }
        (void) close_file( file );

    }

    // Determine radius for the output sphere 
    // sphere is not always perfectly spherical, thus use average radius
	float sphereRadius = 0.0, r;
	for( i=0; i<polygons_src_sphere->n_points; i++) {
    	r = 0.0;
    	for( k=0; k<3; k++) 
    		r += Point_coord(polygons_src_sphere->points[i],k)*Point_coord(polygons_src_sphere->points[i],k);
    	sphereRadius += sqrt(r);
	}
	sphereRadius /= polygons_src_sphere->n_points;
		
    // Determine centre of sphere based on bounds of input (to correct for shiftings)
    float bounds[6];
    get_bounds(polygons_src_sphere, bounds);
    fill_Point( centre, bounds[0]+bounds[1], bounds[2]+bounds[3], bounds[4]+bounds[5] );
    
    objects_dest_sphere = create_object( POLYGONS );
    polygons_dest_sphere = get_polygons_ptr(objects_dest_sphere);
    
    // make radius slightly smaller to get sure that the inner side of handles will be found as
    // nearest point on the surface
    sphereRadius *= 0.975;
    create_tetrahedral_sphere( &centre, sphereRadius, sphereRadius, sphereRadius,
                               n_triangles, polygons_dest_sphere );

    create_polygons_bintree( polygons_src_sphere,
                             ROUND( (Real) polygons_src_sphere->n_items * 0.5 ) );

	ALLOC(new_points, polygons_dest_sphere->n_points);

    for_less( i, 0, polygons_dest_sphere->n_points ) {
        poly =  find_closest_polygon_point( &polygons_dest_sphere->points[i], polygons_src_sphere,
                                              &point_on_src_sphere );
		
		n_points = get_polygon_points( polygons_src_sphere, poly, poly_points_src );
        get_polygon_interpolation_weights( &point_on_src_sphere, n_points, poly_points_src,
                                           weights );

        if( get_polygon_points( polygons_fiducial, poly, poly_points ) != n_points )
            handle_internal_error( "map_point_between_polygons" );

        fill_Point( new_points[i], 0.0, 0.0, 0.0 );
        if( values_specified )
            output_values[i] = 0.0;

        for_less( k, 0, n_points )
        {
            SCALE_POINT( scaled_point, poly_points[k], weights[k] );
            ADD_POINTS( new_points[i], new_points[i], scaled_point );
            if( values_specified )
                output_values[i] += weights[k] * input_values[polygons_fiducial->indices[
                            POINT_INDEX(polygons_fiducial->end_indices,poly,k)]];
        }
   }

    create_polygon_point_neighbours( polygons_dest_sphere, TRUE, &n_neighbours,
                                     &neighbours, NULL, NULL );

    for_less( i, 0, polygons_dest_sphere->n_points ) {
			polygons_dest_sphere->points[i] = new_points[i];
	}
		
    compute_polygon_normals( polygons_dest_sphere );

    output_graphics_any_format( output_filename, format, 1, &objects_dest_sphere );
    
    if( values_specified )
    {
        if( open_file( output_values_filename, WRITE_FILE, ASCII_FORMAT, &file ) != OK)
            return( 1 );

        for_less( i, 0, polygons_dest_sphere->n_points )
        {
            if( output_real( file, output_values[i] ) != OK ||
                output_newline( file ) != OK )
                return( 1 );
        }

        (void) close_file( file );
        FREE( input_values );
        FREE( output_values );
    }
        
	FREE(new_points);
    return( 0 );
}
