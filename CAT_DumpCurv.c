/* Christian Gaser - christian.gaser@uni-jena.de                             */
/* Department of Psychiatry                                                  */
/* University of Jena                                                        */
/*                                                                           */
/* Copyright Christian Gaser, University of Jena.                            */

/* Heat kernel smoothing is based on matlab code from Moo K. Chung:
    Chung, M.K., Robbins,S., Dalton, K.M., Davidson, R.J., Evans, A.C. (2004) 
    Cortical thickness analysis in autism via heat kernel smoothing. NeuroImage, submitted. 
    http://www.stat.wisc.edu/~mchung/papers/ni_heatkernel.pdf */

#include  <volume_io/internal_volume_io.h>
#include  <bicpl.h>
#include  <CAT_Blur2d.h>
#include  <CAT_Curvature.h>

private  void  usage(
    STRING   executable )
{
    STRING   usage_str = "\n\
Usage: %s  object_file output_file [curvtype] [fwhm] [use_abs_values] [-1|1]\n\n\
     Calculate different curvature parameters (default: mean curvature averaged over 3mm)\n\
     from a given obj file.If smoothing filter [fwhm] is defined (in FWHM) a diffusion heat\n\
     kernel will be applied. Optionally the absolute value can ba calculated if\n\
     use_abs_values is defined and only positive or negative values of the values\n\
     can be saved using the last option (-1 or 1).\n\
     curvtype:  0 - mean curvature (averaged over 3mm, in degrees)\n\
                1 - gaussian curvature\n\
                2 - curvedness\n\
                3 - shape index\n\
                4 - mean curvature (in radians)\n\n";

    print_error( usage_str, executable );
}

int  main(
    int   argc,
    char  *argv[] )
{
    STRING               object_filename, output_filename;
    FILE                 *file;
    File_formats         format;
    int                  i, j, n_iter, n_objects, curvtype, sign;
    int                  *n_neighbours, **neighbours, use_abs_values;
    object_struct        **objects;
    polygons_struct      *polygons;
    Smallest_int	     *done_flags;
    Real                 fwhm, sigma, *curvatures, value, distance;
    BOOLEAN		         smoothing;
    progress_struct	     progress;

    initialize_argument_processing( argc, argv );

    if( !get_string_argument( NULL, &object_filename ) ||
        !get_string_argument( NULL, &output_filename ) )
    {
        usage( argv[0] );
        return( 1 );
    }

    (void) get_int_argument( 0, &curvtype );
    (void) get_real_argument( 0.0, &fwhm );
    (void) get_int_argument( 0, &use_abs_values );
    (void) get_int_argument( 0, &sign );
    
    if( fwhm > 0 ) smoothing = 1; else smoothing = 0;
    if( input_graphics_any_format( object_filename,
                             &format, &n_objects, &objects ) != OK )
        return( 1 );

    if( n_objects != 1 || get_object_type( objects[0] ) != POLYGONS )
    {
        print( "File must contain 1 polygons object.\n" );
        return( 1 );
    }

    polygons = get_polygons_ptr( objects[0] );

    ALLOC( curvatures, polygons->n_points );
    
    get_all_polygon_point_neighbours( polygons, &n_neighbours, &neighbours);

    if( curvtype == 0 ) distance = 3.0; else distance = 0.0;
    get_polygon_vertex_curvatures_cg( polygons, n_neighbours, neighbours,
                                   distance, curvtype, curvatures );

    /* limit range to values between -1..1 for all curvtypes > 0 
    (don't ask me where the large values come from...) */
    if( curvtype > 0 ) 
        for_less( i, 0, polygons->n_points ) {
            if(curvatures[i] < -1) curvatures[i] = -1;
            if(curvatures[i] >  1) curvatures[i] = 1;
        }
    
    /* use absolute value */
    if( use_abs_values ) 
        for_less( i, 0, polygons->n_points )
            curvatures[i] = fabs(curvatures[i]);

    /* use positive or negative values only if option is used */
    if(sign > 0)
        for_less( i, 0, polygons->n_points )
            if(curvatures[i] < 0) curvatures[i] = 0;

    if(sign < 0)
        for_less( i, 0, polygons->n_points )
            if(curvatures[i] > 0) curvatures[i] = 0;

    /* and smooth curvatures */
    if( smoothing )
    {
        /* calculate n_iter for sigma = 1.0 */
        n_iter = ROUND(fwhm/2.35482 * fwhm/2.35482);
        if( n_iter == 0 )
            n_iter = 1;
        sigma = 1.0;
        ALLOC( done_flags, polygons->n_points );
        initialize_progress_report( &progress, FALSE, n_iter*polygons->n_points,
                                "Blurring" );
				
        for_less( j, 0, n_iter )
        {
            for_less( i, 0, polygons->n_points )
            {
                heatkernel_blur_points( polygons->n_points, polygons->points, curvatures,
                              n_neighbours[i], neighbours[i],i, sigma, NULL, &value );
	            curvatures[i] = value;
                update_progress_report( &progress, j*polygons->n_points + i + 1 );
            }
        }
        terminate_progress_report( &progress );
        FREE(done_flags);

    }
    
    if( open_file( output_filename, WRITE_FILE, ASCII_FORMAT, &file ) != OK )
        return( 1 );

    for_less( i, 0, polygons->n_points )
    {
        if( output_real( file, curvatures[i] ) != OK ||
            output_newline( file ) != OK )
            break;
    }

    (void) close_file( file );

    delete_object_list( n_objects, objects );
    FREE(curvatures);
    
    return( 0 );
}
