
#include  <volume_io/internal_volume_io.h>
#include  <bicpl.h>
#include "CAT_NiftiIO.h"

int  main(
    int   argc,
    char  *argv[] )
{
    STRING     input_filename;
    double       x, y, z, value, voxel[3];
    Volume     volume;
    BOOLEAN   interpolating_dimensions[MAX_DIMENSIONS];
    int d,axis;

    initialize_argument_processing( argc, argv );

    if( !get_string_argument( "", &input_filename ) )
    {
        return( 1 );
    }
 
    if( input_volume_all( input_filename, 3, XYZ_dimension_names,
                      NC_UNSPECIFIED, FALSE, 0.0, 0.0,
                      TRUE, &volume, (minc_input_options *) NULL ) != OK )
        return( 1 );

    for_less( d, 0, N_DIMENSIONS )
    {
        axis = volume->spatial_axes[d];
        if( axis < 0 )
        {
            print_error(
                  "evaluate_volume_in_world(): must have 3 spatial axes.\n" );
            return(1);
        }

        interpolating_dimensions[axis] = TRUE;
    }

    while( get_real_argument( 0.0, &x ) &&
           get_real_argument( 0.0, &y ) &&
           get_real_argument( 0.0, &z ) )
    {
        evaluate_volume_in_world( volume, x, y, z, -1, FALSE, 0.0,
                                  &value,
                                  NULL, NULL, NULL,
                                  NULL, NULL, NULL, NULL, NULL, NULL );
                        convert_world_to_voxel( volume, x,y,z, voxel );
                        printf("%f %f %f %f %f %f\t%f\n",x,y,z,voxel[X],voxel[Y],voxel[Z],value);

        print( "%g\n", value );

voxel[X] = x;
voxel[Y] = y;
voxel[Z] = z;
        evaluate_volume( volume, voxel, interpolating_dimensions, 0, FALSE, 0.0,
                                  &value,
                                  NULL, NULL );

        print( "%g\n", value );
    }

    return( 0 );
}
