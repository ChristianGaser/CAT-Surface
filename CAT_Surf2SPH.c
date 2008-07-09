/* Christian Gaser - christian.gaser@uni-jena.de                             */
/* Department of Psychiatry                                                  */
/* University of Jena                                                        */
/*                                                                           */
/* Copyright Christian Gaser, University of Jena.                            */

#include  <volume_io/internal_volume_io.h>
#include  <bicpl.h>
#include  <fftw3.h>
#include  "CAT_SheetIO.h"
#include  "CAT_Map2d.h"
#include  "CAT_Blur2d.h"
#include  "CAT_SPH.h"

private  void  usage(
    STRING   executable )
{
    static  STRING  usage_str = "\n\
Usage: %s  surface.obj sphere.obj SPH.dat bandwidth\n\
Extract spherical harmonic coefficients of a surface using 2*bandwidth x 2*bandwidth sample points.\n\
\n\n";

    print_error( usage_str, executable );
}

int  main(
    int   argc,
    char  *argv[] )
{
    STRING               surface_filename, sphere_filename, SPH_filename;
    File_formats         format;
    polygons_struct      *polygons, *polygons_sphere;
    int                  i, j, k, n_objectspoint_index;
    int                  n_items, bandwidth, bandwidth2;
    double               *rcoeffsx, *icoeffsx, *rcoeffsy, *icoeffsy, *rcoeffsz, *icoeffsz;
    double               *rdatax, *rdatay, *rdataz;
    object_struct        **objects;
    int                  n_objects, n_dims = 3;
    int                  dataformat; /* dataformat =0 -> samples are complex, =1 -> samples real */

    /*--- get the arguments from the command line */

    initialize_argument_processing( argc, argv );

    if( !get_string_argument( NULL, &surface_filename ) ||
        !get_string_argument( NULL, &sphere_filename ) ||
        !get_string_argument( NULL, &SPH_filename ) )
    {
        usage( argv[0] );
        return( 1 );
    }

    get_int_argument( 256, &bandwidth );
    bandwidth2 = 2*bandwidth;
    
    if( input_graphics_any_format( surface_filename, &format, &n_objects, &objects ) != OK )
            return( 1 );
    /* check that the surface file contains a polyhedron */
    if( n_objects != 1 || get_object_type( objects[0] ) != POLYGONS ) {
            print( "Surface file must contain 1 polygons object.\n" );
            return( 1 );
    }
    /* get a pointer to the surface */
    polygons = get_polygons_ptr( objects[0] );
    
    if( input_graphics_any_format( sphere_filename, &format, &n_objects, &objects ) != OK )
            return( 1 );
    /* check that the surface file contains a polyhedron */
    if( n_objects != 1 || get_object_type( objects[0] ) != POLYGONS ) {
            print( "Surface file must contain 1 polygons object.\n" );
            return( 1 );
    }
    /* get a pointer to the surface */
    polygons_sphere = get_polygons_ptr( objects[0] );

    /* check that surface and sphere are same size */
    if( polygons->n_items != polygons_sphere->n_items ) {
            print( "Surface and sphere must have same size.\n" );
            return( 1 );
    }
    
    rdatax   = (double *) malloc(sizeof(double)*bandwidth2*bandwidth2);
    rdatay   = (double *) malloc(sizeof(double)*bandwidth2*bandwidth2);
    rdataz   = (double *) malloc(sizeof(double)*bandwidth2*bandwidth2);
    rcoeffsx = (double *) malloc(sizeof(double)*bandwidth*bandwidth);
    rcoeffsy = (double *) malloc(sizeof(double)*bandwidth*bandwidth);
    rcoeffsz = (double *) malloc(sizeof(double)*bandwidth*bandwidth);
    icoeffsx = (double *) malloc(sizeof(double)*bandwidth*bandwidth);
    icoeffsy = (double *) malloc(sizeof(double)*bandwidth*bandwidth);
    icoeffsz = (double *) malloc(sizeof(double)*bandwidth*bandwidth);
    		
    fprintf(stderr,"%30s\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b","Sample spherical coordinates..");
    get_equally_sampled_coords_of_polygon(polygons, polygons_sphere, bandwidth, 
        rdatax, rdatay, rdataz);

    // dataformat indicates real data
    dataformat = 1;
    fprintf(stderr,"%30s\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b","Forward SPH transform...      ");
    get_sph_coeffs_of_realdata(rdatax, bandwidth, dataformat, rcoeffsx, icoeffsx);
    get_sph_coeffs_of_realdata(rdatay, bandwidth, dataformat, rcoeffsy, icoeffsy);
    get_sph_coeffs_of_realdata(rdataz, bandwidth, dataformat, rcoeffsz, icoeffsz);

    if (write_SPHxyz(SPH_filename, bandwidth, rcoeffsx, rcoeffsy, rcoeffsz,
        icoeffsx, icoeffsy, icoeffsz) != OK)
          return( 1 );

    /* clean up */        
    free( rcoeffsx );
    free( rcoeffsy );
    free( rcoeffsz );
    free( icoeffsx );
    free( icoeffsy );
    free( icoeffsz );
    free( rdatax );
    free( rdatay );
    free( rdataz );

	fprintf(stderr,"%30s\n","Done                          ");

    return( 0 );
}
