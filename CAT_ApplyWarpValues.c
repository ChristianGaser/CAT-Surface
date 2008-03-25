/* Christian Gaser - christian.gaser@uni-jena.de                             */
/* Department of Psychiatry                                                  */
/* University of Jena                                                        */
/*                                                                           */
/* Copyright Christian Gaser, University of Jena.                            */

#include  <volume_io/internal_volume_io.h>
#include  <bicpl.h>
#include  <CAT_Map2d.h>

#define  BINTREE_FACTOR   0.5
        
private  void  usage(
    STRING   executable )
{
    STRING  usage_str = "\n\
Usage: %s input.txt ShiftField output.txt\n\n\
     Applies deformations of warping to surface values.\n\n";

    print_error( usage_str, executable );
}

int  main(
    int   argc,
    char  *argv[] )
{
    STRING               output_filename, vector_filename, values_filename;
    FILE		         *outputfile, *inputfile;
    File_formats         format;
    polygons_struct      *polygons, unit_sphere;
    int                  i, x, y, point_index;
    int                  degree, poly, size, ind, n_done;
    int			         n_values;
    Vector2D		     *image;
    BOOLEAN              use_volume;
    Point                unit_point, on_sphere_point, centre;
    Point                poly_points[1000], centroid;
    Real                 u, v, *values, *input_values, **sheet;
    Real                 weights[1000], value;
    object_struct        *object;    
    Vector               normal;
    progress_struct      progress;
    Header               raw_header;
    Real                 tmp_x, tmp_y;
    
    /*--- get the arguments from the command line */

    initialize_argument_processing( argc, argv );

    if( !get_string_argument( NULL, &values_filename ) ||
        !get_string_argument( NULL, &vector_filename ) ||
        !get_string_argument( NULL, &output_filename ))
    {
        usage( argv[0] );
        return( 1 );
    }

    if( input_texture_values( values_filename, &n_values, &input_values ) != OK )
        return( 1 );

    object = create_object( POLYGONS );
    polygons = get_polygons_ptr(object);
    fill_Point( centre, 0.0, 0.0, 0.0 );
    create_tetrahedral_sphere( &centre, 1, 1, 1, 2*(n_values-2), polygons );

    if((inputfile = fopen(vector_filename, "rb")) == NULL)
    {
        fprintf(stderr, "Error: Couldn't read file %s.\n", vector_filename);
        return(1);
    }

    fread(&raw_header, 1, sizeof(raw_header), inputfile);
    ALLOC(image, raw_header.x*raw_header.y*2);
    fread(image, raw_header.x*raw_header.y*2, sizeof(float), inputfile);
    fclose(inputfile);

    /* don't ask me why, but sometimes there are artifacts in the vectorfile indicated
    by very large/small values. We set these values to zero */
    for_less( i, 0, raw_header.x*raw_header.y )
    {
        if(( image[i].x > 30) || ( image[i].x < -30))
            image[i].x = 0;
        if(( image[i].y > 30) || ( image[i].y < -30))
            image[i].y = 0;
    }

    /*--- create a unit sphere with same number of triangles as skin surface */
    fill_Point( centre, 0.0, 0.0, 0.0 );

    create_tetrahedral_sphere( &centre, 1.0, 1.0, 1.0,
                               polygons->n_items, &unit_sphere );

    create_polygons_bintree( &unit_sphere,
                             ROUND( (Real) unit_sphere.n_items *
                                    BINTREE_FACTOR ) );

    ALLOC(values, polygons->n_points);
    
    initialize_progress_report( &progress, FALSE, raw_header.x, "Mapping to sheet" );

    if( open_file( output_filename, WRITE_FILE, ASCII_FORMAT, &outputfile ) != OK ) return( 1 );

    ALLOC2D(sheet, raw_header.x,raw_header.y);

	tmp_x = raw_header.x - 1;
	tmp_y = raw_header.y - 1;

    /*--- map to flat sheet image */
    for ( x=0; x<raw_header.x; x++ ) {
        for ( y=0; y<raw_header.y; y++ ) {
        
            u = ((Real) x)/(Real) tmp_x;
            v = ((Real) y)/(Real) tmp_y;
  
            uv_to_point(u, v, &unit_point );
            poly = find_closest_polygon_point( &unit_point, &unit_sphere,
                                               &on_sphere_point );
            
            size = get_polygon_points( &unit_sphere, poly, poly_points );

            get_polygon_interpolation_weights( &on_sphere_point, size,
                                                   poly_points, weights );

            value = 0.0;
            for ( i=0; i<size; i++ ) {
                ind = unit_sphere.indices[
                         POINT_INDEX(unit_sphere.end_indices,poly,i)];
                value += weights[i] * input_values[ind];
            }
            sheet[x][y] = value;
        }
        update_progress_report( &progress, x + 1 );
    }

	terminate_progress_report( &progress );

    initialize_progress_report( &progress, FALSE, polygons->n_points, "Mapping to sphere" );
	
	/*--- remap to sphere */
    for ( i=0; i<polygons->n_points; i++ ) {

        point_to_uv(&unit_sphere.points[i], &u, &v);
        
        ind = ROUND(u*tmp_x) + raw_header.x*ROUND(v*tmp_y);
        
        x = (int) ROUND(u*tmp_x - image[ind].x);
        y = (int) ROUND(v*tmp_y - image[ind].y);
        if(x > tmp_x) x = x - raw_header.x; 
        if(y > tmp_y) y = y - raw_header.y; 
        if(x < 0) x = raw_header.x + x; 
        if(y < 0) y = raw_header.y + y; 
	
        values[i] = sheet[x][y];
		if( output_real( outputfile, values[i] ) != OK || output_newline( outputfile ) != OK )
            break;
 
        update_progress_report( &progress, i + 1 );
    }

    terminate_progress_report( &progress );

    (void) close_file( outputfile );

    delete_polygons( &unit_sphere );

    FREE( image );
    FREE( values );
 
    return( 0 );
}

