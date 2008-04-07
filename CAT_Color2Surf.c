/* Christian Gaser - christian.gaser@uni-jena.de                             */
/* Department of Psychiatry                                                  */
/* University of Jena                                                        */
/*                                                                           */
/* Copyright Christian Gaser, University of Jena.                            */

#include  <volume_io/internal_volume_io.h>
#include  <bicpl.h>
#include <float.h>

private  void  usage(
    STRING   executable )
{
    STRING  usage_str = "\n\
Usage: %s  src.obj values_file dest.obj\n\
           gray|hot|spectral|blue|green|red|user [user.map] [low high]\n\
           [under|BG] [over|BG] [opacity] [replace||add]\n\n";

    print_error( usage_str, executable );
}

#define  GRAY_STRING       "gray"
#define  HOT_STRING        "hot"
#define  SPECTRAL_STRING   "spectral"
#define  RED_STRING        "red"
#define  GREEN_STRING      "green"
#define  BLUE_STRING       "blue"
#define  USER_STRING       "user"

typedef  enum { REPLACE_COLOUR, ADD_COLOUR }
              Composite_methods;
    
int  main(
    int    argc,
    char   *argv[] )
{
    Real                 value, *values, min_range, max_range;
    Status               status;
    STRING               src_filename, dest_filename, values_filename;
    STRING               under_colour_name, over_colour_name;
    STRING               user_def_filename;
    int                  i, p, n_objects, n_points, n_values, value_index, value_index2;
    Point                *points;
    File_formats         format;
    object_struct        **object_list;
    Colour               *colours, under_colour, over_colour, prev_colour, col;
    STRING               default_over;
    Colour_coding_types  coding_type;
    colour_coding_struct colour_coding;
    Colour_flags         *colour_flag_ptr;
    STRING               coding_type_string;
    Real                 low, high, r, g, b, a, opacity;
    BOOLEAN              per_vertex;
    Composite_methods    composite_method;
    STRING               composite_method_name;

    initialize_argument_processing( argc, argv );

    if( !get_string_argument( NULL, &src_filename ) ||
        !get_string_argument( NULL, &values_filename ) ||
        !get_string_argument( NULL, &dest_filename ) ||
        !get_string_argument( NULL, &coding_type_string ) )
    {
        usage( argv[0] );
        return( 1 );
    }

    default_over = "WHITE";

    if( equal_strings( coding_type_string, GRAY_STRING ) )
        coding_type = GRAY_SCALE;
    else if( equal_strings( coding_type_string, HOT_STRING ) )
        coding_type = HOT_METAL;
    else if( equal_strings( coding_type_string, SPECTRAL_STRING ) )
        coding_type = SPECTRAL;
    else if( equal_strings( coding_type_string, RED_STRING ) )
        coding_type = RED_COLOUR_MAP;
    else if( equal_strings( coding_type_string, GREEN_STRING ) )
        coding_type = GREEN_COLOUR_MAP;
    else if( equal_strings( coding_type_string, BLUE_STRING ) )
        coding_type = BLUE_COLOUR_MAP;
    else if( equal_strings( coding_type_string, USER_STRING ) )
    {
        coding_type = USER_DEFINED_COLOUR_MAP;
        if( !get_string_argument( NULL, &user_def_filename ) )
        {
            usage( argv[0] );
            print_error( "Error in user def colour coding argument.\n");
            return( 1 );
        }
    }
    else
    {
        coding_type = SINGLE_COLOUR_SCALE;
        default_over = coding_type_string;
    }

    get_real_argument( 0.0, &low );
    get_real_argument( 0.0, &high );

    (void) get_string_argument( "BG", &under_colour_name );
    (void) get_string_argument( default_over, &over_colour_name );
    (void) get_real_argument( 1.0, &opacity );
    (void) get_string_argument( "replace", &composite_method_name );

    if( equal_strings( composite_method_name, "replace" ) )
        composite_method = REPLACE_COLOUR;
    else if( equal_strings( composite_method_name, "add" ) )
        composite_method = ADD_COLOUR;
    else
    {
        usage( argv[0] );
        return( 1 );
    }
    
    under_colour = convert_string_to_colour( under_colour_name );
    over_colour = convert_string_to_colour( over_colour_name );

    initialize_colour_coding( &colour_coding, coding_type,
                              under_colour, over_colour, low, high );

    if( coding_type == USER_DEFINED_COLOUR_MAP )
    {
        if( input_user_defined_colour_coding( &colour_coding,
                                              user_def_filename ) != OK)                {
            print_error( "Error in user defined colour map: %s\n",
                          user_def_filename );
            return( 1 );
        }
    }

    if( input_graphics_any_format( src_filename, &format, &n_objects,
                             &object_list ) != OK )
        return( 1 );

    if( input_texture_values( values_filename, &n_values, &values ) != OK )
        return( 1 );
      
    value_index = 0;

    for_less( i, 0, n_objects )
    {
        n_points = get_object_points( object_list[i], &points );

        colour_flag_ptr = get_object_colours( object_list[i], &colours );

        if( *colour_flag_ptr == PER_VERTEX_COLOURS )
            per_vertex = TRUE;
        else
        {
            per_vertex = FALSE;
            prev_colour = colours[0];
            REALLOC( colours, n_points );
            set_object_colours( object_list[i], colours );
        }

        *colour_flag_ptr = PER_VERTEX_COLOURS;

        min_range = FLT_MAX;
        max_range = FLT_MIN;
        
        value_index2 = 0;
        for_less( p, 0, n_points )
        {
            value = values[value_index2];
            ++value_index2;
            if( value != 0.0 && value < min_range )
                min_range = value;
            if( value > max_range )
                max_range = value;
        }

        print( "Value range (excluding zero): %g %g\n", min_range, max_range );
        
        if(( low == 0.0 ) && ( high == 0.0)) {
            colour_coding.min_value = min_range;
            colour_coding.max_value = max_range;
        }

        for_less( p, 0, n_points )
        {
            if( value_index >= n_values )
            {
                print_error( "Insufficient number of values in file.\n" );
                return( 1 );
            }

            value = values[value_index];
            ++value_index;

            col = get_colour_code( &colour_coding, value );

            if( per_vertex )
                prev_colour = colours[p];

            r = get_Colour_r_0_1( col );
            g = get_Colour_g_0_1( col );
            b = get_Colour_b_0_1( col );
            a = opacity * get_Colour_a_0_1( col );

            switch( composite_method )
            {
            case REPLACE_COLOUR:
                /* replace values outside range with background if specified */
                if((equal_strings( under_colour_name, "BG" ) && value < low) || (equal_strings( over_colour_name, "BG" ) && value > high))
                {
                    r = get_Colour_r_0_1( prev_colour );
                    g = get_Colour_g_0_1( prev_colour );
                    b = get_Colour_b_0_1( prev_colour );
                    a = 1.0;
                } else
                {
                    col = make_rgba_Colour_0_1( r, g, b, a);
                    COMPOSITE_COLOURS( col, col, prev_colour );
                    r = get_Colour_r_0_1( col );
                    g = get_Colour_g_0_1( col );
                    b = get_Colour_b_0_1( col );
                    a = 1.0;
                }
                break;

            case ADD_COLOUR:
                r = r + get_Colour_r_0_1( prev_colour );
                g = g + get_Colour_g_0_1( prev_colour );
                b = b + get_Colour_b_0_1( prev_colour );
                a = a + get_Colour_a_0_1( prev_colour );
                break;
            }

            colours[p] = make_rgba_Colour_0_1( r, g, b, a );
        }
    }

    status = output_graphics_any_format( dest_filename, format,
                                   n_objects, object_list );

    return( status != OK );
}
