/* Christian Gaser - christian.gaser@uni-jena.de                             */
/* Department of Psychiatry                                                  */
/* University of Jena                                                        */
/*                                                                           */
/* Copyright Christian Gaser, University of Jena.                            */

#include  <volume_io/internal_volume_io.h>
#include  <bicpl.h>
#include  <ParseArgv.h>
#include  "CAT_SheetIO.h"
#include  "CAT_Map2d.h"
#include  "CAT_Blur2d.h"
#include  "CAT_SPH.h"

/* argument defaults */
int   bandwidth = 256;
int   bandwidth_limited = 64;
int   n_triangles = 81920;

/* the argument table */
ArgvInfo argTable[] = {
  { "-bw", ARGV_INT, (char *) 1, 
    (char *) &bandwidth,
    "Bandwidth of coefficients for spherical harmonic expansion." },
  { "-lim", ARGV_INT, (char *) 1, 
    (char *) &bandwidth_limited,
    "Limit bandwidth of spherical harmonic expansion." },
  { "-n", ARGV_INT, (char *) 1, 
    (char *) &n_triangles,
    "Number of triangles for sampled surface." },
  { NULL, ARGV_END, NULL, NULL, NULL }
};

private  void  usage(
    STRING   executable )
{
    static  STRING  usage_str = "\n\
Usage: %s  surface.obj sphere.obj output.obj\n\
Use spherical harmonic coefficients to create an equally sampled surface.\n\
\n\n";

    print_error( usage_str, executable );
}

int  main(
    int   argc,
    char  *argv[] )
{
    STRING               surface_filename, sphere_filename, output_filename;
    File_formats         format;
    int                  dataformat; /* dataformat =0 -> samples are complex, =1 -> samples real */
    int                  i, bandwidth2, n_dims, n_objects;
    double               *rcoeffsx, *icoeffsx, *rcoeffsy, *icoeffsy, *rcoeffsz, *icoeffsz;
    double               *rcoeffsx2, *icoeffsx2, *rcoeffsy2, *icoeffsy2, *rcoeffsz2, *icoeffsz2;
    double               *rdatax, *rdatay, *rdataz;
    polygons_struct      *polygons, *polygons_sphere, *polygons_output;
    object_struct        **objects, *objects_output;

    /* Call ParseArgv */
    if ( ParseArgv( &argc, argv, argTable, 0 ) || ( argc != 4 ) ) {
        usage( argv[0] );
        fprintf(stderr, "       %s -help\n\n", argv[0]);
        return( 1 );
    }

    initialize_argument_processing( argc, argv );

    if( !get_string_argument( NULL, &surface_filename ) ||
        !get_string_argument( NULL, &sphere_filename ) ||
        !get_string_argument( NULL, &output_filename ) )
    {
        usage( argv[0] );
        return( 1 );
    }
     
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
    
    bandwidth2 = 2*bandwidth;

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

    if ( bandwidth_limited > 0 && bandwidth_limited < bandwidth) {
        fprintf(stderr,"%30s\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b","Limit bandwidth...            ");
        limit_bandwidth(bandwidth, bandwidth_limited, rcoeffsx, rcoeffsx);
        limit_bandwidth(bandwidth, bandwidth_limited, rcoeffsy, rcoeffsy);
        limit_bandwidth(bandwidth, bandwidth_limited, rcoeffsz, rcoeffsz);
        limit_bandwidth(bandwidth, bandwidth_limited, icoeffsx, icoeffsx);
        limit_bandwidth(bandwidth, bandwidth_limited, icoeffsy, icoeffsy);
        limit_bandwidth(bandwidth, bandwidth_limited, icoeffsz, icoeffsz);
    }
    
    fprintf(stderr,"%30s\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b","Inverse SPH transform...      ");
    get_realdata_from_sph_coeffs(rdatax, bandwidth, dataformat, rcoeffsx, icoeffsx);
    get_realdata_from_sph_coeffs(rdatay, bandwidth, dataformat, rcoeffsy, icoeffsy);
    get_realdata_from_sph_coeffs(rdataz, bandwidth, dataformat, rcoeffsz, icoeffsz);

    objects_output = create_object( POLYGONS );
    polygons_output = get_polygons_ptr(objects_output);

    fprintf(stderr,"%30s\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b","Resample surface...           ");
    sample_sphere_from_sph(rdatax, rdatay, rdataz, polygons_output, n_triangles, bandwidth);

    if( output_graphics_any_format( output_filename, ASCII_FORMAT, 1,
                &objects_output ) != OK )
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

    return(0);    
}