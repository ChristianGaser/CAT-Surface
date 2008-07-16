/* Christian Gaser - christian.gaser@uni-jena.de                             */
/* Department of Psychiatry                                                  */
/* University of Jena                                                        */
/*                                                                           */
/* most of the code is modified from caret/caret_brain_set/BrainModelSurface.cxx  */
/*                                                                           */
/* Copyright Christian Gaser, University of Jena.                            */

#include  <volume_io/internal_volume_io.h>
#include  <bicpl.h>
#include  <CAT_Surf.h>
#include <float.h>

private  void  usage(
    STRING   executable )
{
    STRING  usage_str = "\n\
Usage: %s surface.obj sphere.obj [stop_at]\n\n\
     Maps a surface to a sphere using the caret inflating approach.\n\n\
     The inflating can be limited using stop_at (default 5), where\n\
       1 - Low smooth\n\
       2 - Inflating\n\
       3 - Very inflating\n\
       4 - High smoothing\n\
       5 - Ellipsoid\n";

    print_error( usage_str, executable );
}

int  main(
    int    argc,
    char   *argv[] )
{
    STRING           input_filename, output_filename;
    int              n_objects, i, stop_at;
    FILE             *file;
    File_formats     format;
    object_struct    **object_list;
    polygons_struct  *polygons;
    BOOLEAN          enableFingerSmoothing = 1;
    int              fingerSmoothingIterations;
    
    initialize_argument_processing( argc, argv );

    if( !get_string_argument( NULL, &input_filename ) ||
        !get_string_argument( NULL, &output_filename ) )
    {
        usage( argv[0] );
        return( 1 );
    }

    get_int_argument(5, &stop_at);
    
    if( input_graphics_any_format( input_filename, &format, &n_objects,
                             &object_list ) != OK || n_objects != 1 ||
        get_object_type(object_list[0]) != POLYGONS )
    {
        print_error( "Error reading %s.\n", input_filename );
        return( 1 );
    }

    polygons = get_polygons_ptr( object_list[0] );


    if( euler_characteristic( polygons ) !=2) {
        print_error( "Euler characteristic of %s must be 2.\n", input_filename );
    }
     
    float fiducialSurfaceArea = get_polygons_surface_area( polygons );

    // low smooth
	fprintf(stderr,"%20s\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b","Low smoothing...    ");
    inflate_surface_and_smooth_fingers(polygons,
                                                   1,     // number of cycles
                                                   0.2,   // regular smoothing strength
                                                   50,    // regular smoothing iterations
                                                   1.0,   // inflation factor
                                                   3.0,   // finger compress/stretch threshold
                                                   1.0,   // finger smoothing strength
                                                   0);     // finger smoothing iterations
                                                   
    if (stop_at > 1) {
      // inflated                                                      
      fingerSmoothingIterations = 0;
      if (enableFingerSmoothing)
        fingerSmoothingIterations = 30;
	  fprintf(stderr,"%20s\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b","Inflating...        ");
      inflate_surface_and_smooth_fingers(polygons,
                                                   2,     // number of cycles
                                                   1.0,   // regular smoothing strength
                                                   30,    // regular smoothing iterations
                                                   1.4,   // inflation factor
                                                   3.0,   // finger compress/stretch threshold
                                                   1.0,   // finger smoothing strength
                                                   fingerSmoothingIterations);     // finger smoothing iterations
    }                                             
    
    if (stop_at > 2) {
      // very inflated
      fprintf(stderr,"%20s\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b","Very inflating...   ");
      inflate_surface_and_smooth_fingers(polygons,
                                                      4,     // number of cycles
                                                      1.0,   // regular smoothing strength
                                                      30,    // regular smoothing iterations
                                                      1.1,   // inflation factor
                                                      3.0,   // finger compress/stretch threshold
                                                      1.0,   // finger smoothing strength
                                                      0);     // finger smoothing iterations
    }
    
    if (stop_at > 3) {
      // high smooth
      fingerSmoothingIterations = 0;
      if (enableFingerSmoothing) 
    	fingerSmoothingIterations = 60;
	  fprintf(stderr,"%20s\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b","High smoothing...   ");
      inflate_surface_and_smooth_fingers(polygons,
                                                      6,     // number of cycles
                                                      1.0,   // regular smoothing strength
                                                      60,    // regular smoothing iterations
                                                      1.6,   // inflation factor
                                                      3.0,   // finger compress/stretch threshold
                                                      1.0,   // finger smoothing strength
                                                      fingerSmoothingIterations);     // finger smoothing iterations

    }

    if (stop_at > 4) {
      // ellipsoid
      fprintf(stderr,"%20s\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b","Ellipsoid...        ");
      inflate_surface_and_smooth_fingers(polygons,
                                                      6,     // number of cycles
                                                      1.0,   // regular smoothing strength
                                                      50,    // regular smoothing iterations
                                                      1.4,   // inflation factor
                                                      4.0,   // finger compress/stretch threshold
                                                      1.0,   // finger smoothing strength
                                                      fingerSmoothingIterations);     // finger smoothing iterations
      convert_ellipsoid_to_sphere_with_surface_area( polygons, fiducialSurfaceArea);
    }
    
	fprintf(stderr,"Done                \n");

    compute_polygon_normals( polygons );

    output_graphics_any_format( output_filename, format, 1, object_list );

    return( 0 );
}
