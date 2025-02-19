/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 *
 */

#include "CAT_Refine.h"
#include "CAT_SurfaceIO.h"

private  void  usage(
    STRING   executable )
{
    STRING  usage_str = "\n\
Usage: %s  surface_file output_surface_file  max_length [weight_curvatures]\n\
  Refine mesh so that max_length is the largest distance between two points. If curvature\n\
  weighting is used regions in sulci and gyri with large absolute mean curvature values\n\
  will result in a smaller max_length to allow better surface deforming.\n\n";

    fprintf(stderr,"%s\n %s\n",usage_str, executable);
}

int  main(
    int    argc,
    char   *argv[] )
{
    STRING             input_filename, output_filename;
    int                n_objects, n_done, i;
    File_formats       format;
    object_struct      **object_list;
    polygons_struct    *polygons, new_polygons;
    Point              *length_points;
    float              max_length, weight_curvatures;

    initialize_argument_processing( argc, argv );

    if( !get_string_argument( NULL, &input_filename ) ||
        !get_string_argument( NULL, &output_filename ) ||
        !get_real_argument( 0.0, &max_length ) )
    {
        usage( argv[0] );
        return( 1 );
    }

    get_real_argument( 0.0, &weight_curvatures );

    if( input_graphics_any_format( input_filename, &format, &n_objects,
                             &object_list ) != OK ||
        n_objects < 1 || get_object_type(object_list[0]) != POLYGONS )
    {
        fprintf(stderr, "Error in input file.\n");
        return( 1 );
    }

    polygons = get_polygons_ptr( object_list[0] );

    SET_ARRAY_SIZE( length_points, 0, polygons->n_points, DEFAULT_CHUNK_SIZE );
    for_less( i, 0, polygons->n_points )
        length_points[i] = polygons->points[i];

    do
    {
        n_done = refine_mesh( &length_points, polygons, max_length,
                              &new_polygons, weight_curvatures );

        delete_polygons( polygons );
        *polygons = new_polygons;
    }
    while( n_done > 0 );

    fprintf(stdout,"Resampled into %d polygons.\n", polygons->n_items );

    (void) output_graphics_any_format( output_filename, format, n_objects,
                                 object_list, NULL);

    return( 0 );
}