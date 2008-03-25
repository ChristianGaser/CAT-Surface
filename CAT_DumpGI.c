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
Usage: %s  object_file output_file fwhm [fwhm_surf] [curvtype]\n\n\
     Calculate gyrification index using the ratio of local surface area and local inflated\n\
     surface area. Local surface area can be approximated by use of different curve types\n\
     (default: mean curvature averaged over 3mm):\n\
     curvtype:  0 - mean curvature (averaged over 3mm, in degrees)\n\
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
    int                  i, j, n_iter, n_objects, curvtype;
    int                  *n_neighbours, **neighbours;
    object_struct        **objects;
    polygons_struct      *polygons;
    Point                *smooth_points, point;
    Smallest_int	     *done_flags;
    Real                 fwhm, fwhm_surf, *curvatures, *curvatures_inflated;
    Real                 *GI, value, distance, sigma;
    progress_struct	     progress;

    initialize_argument_processing( argc, argv );

    if( !get_string_argument( NULL, &object_filename ) ||
        !get_string_argument( NULL, &output_filename ) )
    {
        usage( argv[0] );
        return( 1 );
    }

    (void) get_real_argument( 30.0, &fwhm );
    (void) get_real_argument( 100.0, &fwhm_surf );
    (void) get_int_argument( 0, &curvtype );
    
    if( input_graphics_file( object_filename,
                             &format, &n_objects, &objects ) != OK )
        return( 1 );

    if( n_objects != 1 || get_object_type( objects[0] ) != POLYGONS )
    {
        print( "File must contain 1 polygons object.\n" );
        return( 1 );
    }

    polygons = get_polygons_ptr( objects[0] );

    compute_polygon_normals( polygons );

    ALLOC( curvatures, polygons->n_points );
    ALLOC( curvatures_inflated, polygons->n_points );
    ALLOC( GI, polygons->n_points );
    ALLOC( smooth_points, polygons->n_points );
    ALLOC( done_flags, polygons->n_points );
    
    get_all_polygon_point_neighbours( polygons, &n_neighbours, &neighbours);
    
    // get curvature values
    if( curvtype == 0 ) distance = 3.0; else distance = 0.0;
    get_polygon_vertex_curvatures_cg( polygons, n_neighbours, neighbours,
                                   distance, curvtype, curvatures );

    // inflate surface by smoothing with FWHM of 150mm
    
    // calculate n_iter with regard to sigma
    sigma = 8.0;
    n_iter = ROUND(fwhm_surf/2.35482 * fwhm_surf/2.35482/sigma);

    for_less( i, 0, polygons->n_points )
        done_flags[i] = FALSE;

    initialize_progress_report( &progress, FALSE, n_iter*polygons->n_points,
                                "Blurring surface" );

	/* diffusion smoothing using heat kernel */			
    for_less( j, 0, n_iter )
    {
        for_less( i, 0, polygons->n_points )
        {
            heatkernel_blur_points( polygons->n_points, polygons->points, NULL,
                n_neighbours[i], neighbours[i], i, sigma, &point, &value );
            smooth_points[i] = point;

            update_progress_report( &progress, j*polygons->n_points + i + 1 );
        }
        for_less( i, 0, polygons->n_points )
        {
            polygons->points[i] = smooth_points[i];
        }
    }
    terminate_progress_report( &progress );

    polygons->points = smooth_points;
    compute_polygon_normals( polygons );

    // get curvature values of inflated surface
    get_polygon_vertex_curvatures_cg( polygons, n_neighbours, neighbours,
                                   distance, curvtype, curvatures_inflated );

    // calculate ratio of absolute values
    for_less( i, 0, polygons->n_points )
        if (curvatures_inflated[i] != 0)
            GI[i] = ABS(curvatures[i])/ABS(curvatures_inflated[i]);
        else
            GI[i] = 0.0;
    
    /* and smooth GI values */
    /* calculate n_iter for sigma = 1.0 */
    n_iter = ROUND(fwhm/2.35482 * fwhm/2.35482);
    sigma = 1.0;

    for_less( i, 0, polygons->n_points )
        done_flags[i] = FALSE;

    initialize_progress_report( &progress, FALSE, n_iter*polygons->n_points,
                                "Blurring values" );
				
    for_less( j, 0, n_iter )
    {
        for_less( i, 0, polygons->n_points )
        {
            heatkernel_blur_points( polygons->n_points, polygons->points, GI,
                          n_neighbours[i], neighbours[i],
                          i, sigma, NULL, &value );
            GI[i] = value;
            update_progress_report( &progress, j*polygons->n_points + i + 1 );
        }
    }
    terminate_progress_report( &progress );

    
    if( open_file( output_filename, WRITE_FILE, ASCII_FORMAT, &file ) != OK )
        return( 1 );

    for_less( i, 0, polygons->n_points )
    {
        if( output_real( file, GI[i] ) != OK ||
            output_newline( file ) != OK )
            break;
    }

    (void) close_file( file );

    delete_object_list( n_objects, objects );
    FREE(curvatures);
    FREE(curvatures_inflated);
    FREE(GI);
    FREE(done_flags);
    
    return( 0 );
}
