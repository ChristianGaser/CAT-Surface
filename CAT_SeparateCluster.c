#include <volume_io/internal_volume_io.h>
#include <bicpl.h>
#include <CAT_Blur2d.h>

private  int   separate_polygons(
    polygons_struct    *polygons,
    int                desired_index,
    Real			   *values);

private  void  usage(
    STRING   executable )
{
    STRING  usage_str = "\n\
Usage: %s  input.obj input.txt output_prefix [which] \n\
\n\
     Separates polygons into its disjoint parts.\n\n";

    print_error( usage_str, executable );
}

int  main(
    int    argc,
    char   *argv[] )
{
    STRING           input_filename, output_prefix, values_filename;
    char             out_filename[EXTREMELY_LARGE_STRING_SIZE];
    int              n_objects, n_out;
    int              i, desired_index, n_values;
    Real			 *values;
    File_formats     format;
    object_struct    **object_list;
    polygons_struct  *polygons;

    initialize_argument_processing( argc, argv );

    if( !get_string_argument( NULL, &input_filename ) ||
        !get_string_argument( NULL, &values_filename ) ||
        !get_string_argument( NULL, &output_prefix ) )
    {
        usage( argv[0] );
        return( 1 );
    }

    (void) get_int_argument( -1, &desired_index );

    if( input_graphics_any_format( input_filename, &format, &n_objects,
                             &object_list ) != OK || n_objects < 1 ||
        get_object_type( object_list[0] ) != POLYGONS )
    {
        print_error( "File must have a polygons structure.\n" );
        return( 1 );
    }

    if( input_texture_values( values_filename, &n_values, &values ) != OK ) {
        print_error( "Cannot read values.\n" );
    	return( 1 );
    }

    polygons = get_polygons_ptr( object_list[0] );

    check_polygons_neighbours_computed( polygons );

    n_out = separate_polygons( polygons, desired_index, values);
fprintf(stderr,"%d\n",n_out);

    for_less( i, 0, n_out )
    {
/*        if( n_out == 1 )
            (void) sprintf( out_filename, "%s", output_prefix );
        else
            (void) sprintf( out_filename, "%s_%d.obj", output_prefix, i+1 );

        (void) output_graphics_any_format( out_filename, format, 1, &out[i] );
*/
		fprintf(stderr,"%g\n",i);
    }

    return( 0 );
}

private  int   make_connected_components(
    polygons_struct    *polygons,
    int       point_classes[],
    Real			   *values,
    int                n_in_class[] )
{
    int                point, edge, size;
    int                neigh, j;
    int                *n_neighbours, **neighbours;
    int                n_components, not_done;
    QUEUE_STRUCT(int)  queue;

    n_components = 0;

    not_done = 16383;

    get_all_polygon_point_neighbours( polygons, &n_neighbours, &neighbours);

    for_less( point, 0, polygons->n_points )
        point_classes[point] = not_done;

    for_less( point, 0, polygons->n_points )
    {
        if( point_classes[point] != not_done )
            continue;

        if( n_components == 16383 )
        {
            ++n_components;
            break;
        }

        point_classes[point] = n_components;
        n_in_class[n_components] = 1;

            
		if (n_neighbours[point] > 0) { 
            for_less( j, 0, n_neighbours[point] ) {
		        if( point_classes[j] == not_done) 
		        {
		         	if (values[neighbours[point][j]]>0.0 && values[point]>0.0)
	        		{
		         		point_classes[j] = n_components;
        	     		++n_in_class[n_components];
        			} else         ++n_components;
        		}

	        }
		}
			
		fprintf(stderr,"%d: %d\n",point, point_classes[point]);
    }

    return( n_components );
}

private  int   separate_polygons(
    polygons_struct    *polygons,
    int                desired_index,
    Real			   *values)
{
    int                ind, p_ind, point, vertex, size, i, j, tmp;
    int                point_index, *new_point_ids, n_objects, comp, c;
    int                biggest;
    int       *point_classes;
    int                n_components, *n_in_class, *ordered;

    ALLOC( point_classes, polygons->n_points );
    ALLOC( n_in_class, 16384 );
    ALLOC( ordered, 16384 );

    n_components = make_connected_components( polygons, point_classes, values,
                                              n_in_class );

    for_less( i, 0, n_components )
        ordered[i] = i;

    for_less( i, 0, n_components-1 )
    {
        biggest = i;
        for_less( j, i+1, n_components )
        {
            if( n_in_class[ordered[j]] > n_in_class[ordered[biggest]] )
                biggest = j;
        }

        tmp = ordered[i];
        ordered[i] = ordered[biggest];
        ordered[biggest] = tmp;
    }

    FREE( point_classes );
    FREE( n_in_class );
    FREE( ordered );

    return( n_objects );
}
