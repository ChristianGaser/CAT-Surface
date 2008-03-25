/* Christian Gaser - christian.gaser@uni-jena.de                             */
/* Department of Psychiatry                                                  */
/* University of Jena                                                        */
/*                                                                           */
/* Copyright Christian Gaser, University of Jena.                            */

#include  <volume_io/internal_volume_io.h>
#include  <bicpl.h>

private  void  usage(
    STRING   executable )
{
    STRING   usage_str = "\n\
Usage: %s object_file x y z values1.txt [values2.txt .. valuesn.txt]\n\n\
     Plot values at coordinate x y z for each file. The object file is used to\n\
     link the coordinates at the surface model to the values.\n\n";

    print_error( usage_str, executable );
}

int  main(
    int   argc,
    char  *argv[] )
{
    STRING               object_filename, input_filename;
    FILE                 *file;
    File_formats         format;
    int                  poly, n_objects, point_index, min_index, n_files, i, j;
    object_struct        **objects;
    polygons_struct      *polygons;
    Real		         x, y, z, dist, min_dist, value;
    Point		         point;
    char                 line[256];

    initialize_argument_processing( argc, argv );

    if( !get_string_argument( NULL, &object_filename ))
    {
        usage( argv[0] );
        return( 1 );
    }

    if( input_graphics_file( object_filename,
                             &format, &n_objects, &objects ) != OK )
        return( 1 );

    if( !get_real_argument( 0.0, &x ) ||
        !get_real_argument( 0.0, &y ) ||
        !get_real_argument( 0.0, &z ) )
    {
        usage( argv[0] );
        return( 1 );
    }

    n_files = argc - 5;
        
    if( n_objects != 1 || get_object_type( objects[0] ) != POLYGONS )
    {
        print( "File must contain 1 polygons object.\n" );
        return( 1 );
    }

    fill_Point( point, x, y, z );
    
    polygons = get_polygons_ptr( objects[0] );
    
    min_dist = 1e15;
    for_less( point_index, 0, polygons->n_points )
    {
        dist = distance_between_points( &polygons->points[point_index], &point );
        if(dist < min_dist)
        {
            min_dist = dist;
            min_index = point_index;
        }
    }

    print("Closest distance of %3.2f found at index %d\n",min_dist,min_index);
    
    for_less(i, 0, n_files)
    {
        get_string_argument( NULL, &input_filename );
        if((file = fopen(input_filename, "r")) == 0) {
            fprintf(stderr, "Couldn't open file %s.\n", input_filename);
            return(0);
        }
        for_less(j, 0, min_index+1)
            fgets(line, 256, file);
        print("%s",line);
        (void) close_file( file );
    }


    delete_object_list( n_objects, objects );

    return( 0 );
}
