/* Christian Gaser - christian.gaser@uni-jena.de                             */
/* Department of Psychiatry                                                  */
/* University of Jena                                                        */
/*                                                                           */
/* Copyright Christian Gaser, University of Jena.                            */

#include  <volume_io/internal_volume_io.h>
#include  <bicpl.h>
#include <float.h>

#include  "ParseArgv.h"

#define  GET_grid_POINT( result, grid_origin, normal, length ) \
    { \
        (result)[X] = RPoint_x(grid_origin) + (length) * Vector_x(normal); \
        (result)[Y] = RPoint_y(grid_origin) + (length) * Vector_y(normal); \
        (result)[Z] = RPoint_z(grid_origin) + (length) * Vector_z(normal); \
    }


#define F_AVERAGE 0
#define F_RANGE   1
#define F_MAX     2
#define F_MIN     3
#define F_EXP     4
#define F_SUM     5

#define LOG05     -0.69314718
#define MAX_N_ARRAY 250

/* argument defaults */
int     degrees_continuity = -1; /* interpolation - default: nearest neighbour */
int     map_func = F_MAX;       /* default mapping function: maximum value */
Real    grid_res = 1.0;         /* resolution of grid along normals in mm */
Real    grid_length = 10.0;      /* length of grid along normals */
Real    grid_origin = 0.0;      /* origin of grid along normals */
Real    range[2] = {FLT_MAX, FLT_MAX};
Real    exp_half = FLT_MAX;     /* */

/* the argument table */
ArgvInfo argTable[] = {
   {"-res", ARGV_FLOAT, (char *) 1, (char *) &grid_res,
       "Resolution of grid along normals [mm]."},
   {"-origin", ARGV_FLOAT, (char *) 1, (char *) &grid_origin,
       "Origin (start point) of grid along normals [mm]. Give negative values for \n\t\torigin outside the surface."},
   {"-length", ARGV_FLOAT, (char *) 1, (char *) &grid_length,
       "Lenght of grid along normals [mm]."},
   {NULL, ARGV_HELP, (char *) NULL, (char *) NULL, 
       "Interpolation options:"},
  { "-linear", ARGV_CONSTANT, (char *) 0, 
    (char *) &degrees_continuity,
    "Use linear interpolation." },
  { "-nearest_neighbour", ARGV_CONSTANT, (char *) -1, 
    (char *) &degrees_continuity,
    "Use nearest neighbour interpolation (Default)." },
  { "-cubic", ARGV_CONSTANT, (char *) 2,
    (char *) &degrees_continuity,
        "Use cubic interpolation." },
   {NULL, ARGV_HELP, (char *) NULL, (char *) NULL, 
       "Mapping function options:"},
  { "-average", ARGV_CONSTANT, (char *) 0, 
    (char *) &map_func,
    "Use average for mapping along normals." },
  { "-range", ARGV_FLOAT, (char *) 2, 
    (char *) range,
    "Count number of values in range for mapping along normals. If any value is out of range \n\t\tvalues will be counted only until this point" },
  { "-max", ARGV_CONSTANT, (char *) 2, 
    (char *) &map_func,
    "Use maximum value for mapping along normals (Default)." },
  { "-min", ARGV_CONSTANT, (char *) 3, 
    (char *) &map_func,
    "Use minimum value for mapping along normals." },
  { "-exp", ARGV_FLOAT, (char *) 1, 
    (char *) &exp_half,
    "Use exponential average of values for mapping along normals. The argument defines the \n\t\tdistance in mm where values are decayed to 50% (recommended value is 10mm)." },
  { "-sum", ARGV_CONSTANT, (char *) 5, 
    (char *) &map_func,
    "Use sum of values for mapping along normals." },
  
  { NULL, ARGV_END, NULL, NULL, NULL }
};


private Real evaluate_function(
    Real    val_array[],
    int     n_val,
    int     map_func,
    Real    exp_array[])
{
    int     i, in_range;
    Real    result;
    
    switch( map_func )
    {
    case F_AVERAGE:
        result = 0.0;
        for_less( i, 0, n_val )
            result += val_array[i]; 
        result /= (Real) n_val;
        break;
    case F_RANGE:
        /* count only if values are in range until any value is out of range */
        in_range = 0;
        result = 0.0;
        for_less( i, 0, n_val )
        {
            /* stop counting if values are leaving range */
            if((in_range == 1) && ((val_array[i] < range[0]) || (val_array[i] > range[1])))
                break;
            /* are values for the first time in range? */
            if((val_array[i] >= range[0]) && (val_array[i] <= range[1]))
                in_range = 1;
            /* count values in range */
            if( in_range )
                result += grid_res;
        }
        break;
    case F_MAX:
        result = -FLT_MAX;
        for_less( i, 0, n_val )
        {
            if( val_array[i] > result)
                result = val_array[i];
        } 
        break;
    case F_MIN:
        result = FLT_MAX;
        for_less( i, 0, n_val )
        {
            if( val_array[i] < result)
                result = val_array[i];
        } 
        break;
    case F_EXP:
        /* exponential average */
        result = 0.0;
        for_less( i, 0, n_val )
            result += val_array[i]*exp_array[i];
        break;
    case F_SUM:
        result = 0.0;
        for_less( i, 0, n_val )
            result += val_array[i]; 
        break;
    }
    return(result);
}

int  main(
    int   argc,
    char  *argv[] )
{
    STRING               input_volume_filename, object_filename;
    STRING               output_filename;
    File_formats         format;
    Volume               volume;
    int                  i, j, n_values, n_objects;
    object_struct        **objects;
    polygons_struct      *polygons;
    Real                 value, voxel[N_DIMENSIONS], val_array[MAX_N_ARRAY], length_array[MAX_N_ARRAY];
    Real                 exp_sum, exp_array[MAX_N_ARRAY];
    FILE                 *file;
    Vector               normal;

    /* Call ParseArgv */
    if ( ParseArgv( &argc, argv, argTable, 0 ) || ( argc != 4 ) ) {
      (void) fprintf( stderr,"\nUsage: %s [options] <volume.mnc> <object.obj> <output.txt>\n", argv[0] );
      (void) fprintf( stderr,"\nMap data from a minc volume to a surface.\n" );
      return( 1 );
    }

    initialize_argument_processing( argc, argv );

    if( !get_string_argument( NULL, &input_volume_filename ) ||
        !get_string_argument( NULL, &object_filename ) ||
        !get_string_argument( NULL, &output_filename ) )
    {
        print_error(
           "Usage: %s  volume.mnc  object.obj  output_file.txt\n", argv[0]);
        return( 1 );
    }
    
    /* calculate number of values needed to subdivide grid along normals */ 
    n_values = (grid_length/grid_res) + 1;
    
    /* check minimum number of values */
    if( n_values < 2 )
    {
        print_error("Resolution of grid is too low.\n");
        return( 1 );
    }
    
    /* .. and also check maximum number of values */
    if( n_values > MAX_N_ARRAY )
    {
        print_error("Resolution of grid is too high.\n");
        return( 1 );
    }

    /* if exp is given use exponential mapping function */
    if( exp_half != FLT_MAX )
        map_func = F_EXP;

    /* if range is given use range mapping function */
    if(( range[0] != FLT_MAX ) && ( range[1] != FLT_MAX))
        map_func = F_RANGE;
        
    /* check range values */
    if( range[0] > range[1])
    {
        print_error("First value of range is larger than second one.\n");
        return( 1 );
    }

    /* initialize values for lengths (starting with origin) */
    fprintf(stderr,"Calculate values along:\n");
    for_less( j, 0, n_values)
    {
        length_array[j] = grid_origin + ((Real)j/(Real)(n_values-1)*grid_length);
        fprintf(stderr,"%3.2f ",length_array[j]);
    }
    fprintf(stderr,"\n");

    /* calculate exponential decay if exp function is defined */
    if( exp_half != FLT_MAX )
    {
        exp_sum = 0.0;
        for_less( j, 0, n_values)
        {
            exp_array[j] = exp(LOG05/exp_half*length_array[j]);
            exp_sum += exp_array[j];
        }
        /* scale sum of exponential function to 1 */
        for_less( j, 0, n_values)
            exp_array[j] /= exp_sum;
    }
    
    if( input_volume( input_volume_filename, 3, File_order_dimension_names,
                      NC_UNSPECIFIED, FALSE, 0.0, 0.0,
                      TRUE, &volume, NULL ) != OK )
        return( 1 );

    if( input_graphics_file( object_filename,
                             &format, &n_objects, &objects ) != OK )
        return( 1 );

    if( n_objects != 1 )
        print( "Warning, more than one object in file: %s\n", object_filename );

    polygons = get_polygons_ptr( objects[0] );

    compute_polygon_normals( polygons );

    if( open_file( output_filename, WRITE_FILE, ASCII_FORMAT, &file ) != OK )
        return( 1 );

    for_less( i, 0, polygons->n_points )
    {
        /* look only for inward normals */
        SCALE_VECTOR( normal, polygons->normals[i], -1.0 );
        for_less( j, 0, n_values)
        {
            /* get point from origin in normal direction */
            GET_grid_POINT(voxel, polygons->points[i], normal, length_array[j]);
            evaluate_volume_in_world( volume,
                                  voxel[X], voxel[Y], voxel[Z],
                                  degrees_continuity, 
                                  FALSE, 0.0,
                                  &value,
                                  NULL, NULL, NULL,
                                  NULL, NULL, NULL, NULL, NULL, NULL );
            val_array[j] = value;
        }
        /* evaluate function */
        value = evaluate_function(val_array, n_values, map_func, exp_array);
        if( output_real( file, value ) != OK || output_newline( file ) != OK )
            return( 1 );
    }

    (void) close_file( file );

    delete_object_list( n_objects, objects );

    return( 0 );
}
