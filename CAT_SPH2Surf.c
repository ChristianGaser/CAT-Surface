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
int   bandwidth_limited = 0;
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
Usage: %s  SPH.dat output.obj [bandwidth_limited n_points]\n\
Use spherical harmonic coefficients to create a surface.\n\
\n\n";

    print_error( usage_str, executable );
}

int  main(
    int   argc,
    char  *argv[] )
{
    FILE                 *fp;
    char                 line[256];
    STRING               surface_filename, SPH_filename;
    int                  dataformat; /* dataformat =0 -> samples are complex, =1 -> samples real */
    int                  i, bandwidth2, n_dims;
    double               *rcoeffsx, *icoeffsx, *rcoeffsy, *icoeffsy, *rcoeffsz, *icoeffsz;
    double               *rcoeffsx2, *icoeffsx2, *rcoeffsy2, *icoeffsy2, *rcoeffsz2, *icoeffsz2;
    double               *rdatax, *rdatay, *rdataz;
    polygons_struct      *polygons_sphere;
    object_struct        *objects_sphere;

    /* Call ParseArgv */
    if ( ParseArgv( &argc, argv, argTable, 0 ) || ( argc != 3 ) ) {
        usage( argv[0] );
        fprintf(stderr, "       %s -help\n\n", argv[0]);
        return( 1 );
    }

    initialize_argument_processing( argc, argv );

    if( !get_string_argument( NULL, &SPH_filename ) ||
        !get_string_argument( NULL, &surface_filename ) )
    {
        usage( argv[0] );
        return( 1 );
    }
     
    // read coefficients
    if((fp = fopen(SPH_filename, "rb")) == NULL) {
        fprintf(stderr, "Error opening file %s.\n", SPH_filename);
        return(1);
    }

    fgets(line, 256, fp);
    if (strncmp(line,"SPH", 3)) {
        fprintf(stderr, "Wrong file format.\n");
        return(1);
    }
    fgets(line, 256, fp);
    sscanf(line, "%d %d %d", &bandwidth, &n_dims, &dataformat);

    if (dataformat !=1) {
        fprintf(stderr, "Wrong dataformat: Data should be real and not complex.\n");
        return(1);
    }

    if (n_dims !=3) {
        fprintf(stderr, "Data dimension should be 3.\n");
        return(1);
    }

    int size = (bandwidth) * (bandwidth);
    rcoeffsx = (double *) malloc(sizeof(double )*size);
    rcoeffsy = (double *) malloc(sizeof(double )*size);
    rcoeffsz = (double *) malloc(sizeof(double )*size);
    icoeffsx = (double *) malloc(sizeof(double )*size);
    icoeffsy = (double *) malloc(sizeof(double )*size);
    icoeffsz = (double *) malloc(sizeof(double )*size);

    if (fread(rcoeffsx, sizeof(double), size, fp) != size) {
        fprintf(stderr, "Error reading data.\n");
        return(1);
    }
    if (fread(icoeffsx, sizeof(double), size, fp) != size) {
        fprintf(stderr, "Error reading data.\n");
        return(1);
    }
    if (fread(rcoeffsy, sizeof(double), size, fp) != size) {
        fprintf(stderr, "Error reading data.\n");
        return(1);
    }
    if (fread(icoeffsy, sizeof(double), size, fp) != size) {
        fprintf(stderr, "Error reading data.\n");
        return(1);
    }
    if (fread(rcoeffsz, sizeof(double), size, fp) != size) {
        fprintf(stderr, "Error reading data.\n");
        return(1);
    }
    if (fread(icoeffsz, sizeof(double), size, fp) != size) {
        fprintf(stderr, "Error reading data.\n");
        return(1);
    }
        
    fclose(fp);

    bandwidth2 = 2*bandwidth;

    rdatax   = (double *) malloc(sizeof(double)*bandwidth2*bandwidth2);
    rdatay   = (double *) malloc(sizeof(double)*bandwidth2*bandwidth2);
    rdataz   = (double *) malloc(sizeof(double)*bandwidth2*bandwidth2);

    if ( bandwidth_limited > 0 && bandwidth_limited < bandwidth) {
        fprintf(stderr,"%30s\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b","Limit bandwidth...            ");
        limit_bandwidth(bandwidth, bandwidth_limited, rcoeffsx, rcoeffsx);
        limit_bandwidth(bandwidth, bandwidth_limited, rcoeffsy, rcoeffsy);
        limit_bandwidth(bandwidth, bandwidth_limited, rcoeffsz, rcoeffsz);
        limit_bandwidth(bandwidth, bandwidth_limited, icoeffsx, icoeffsx);
        limit_bandwidth(bandwidth, bandwidth_limited, icoeffsy, icoeffsy);
        limit_bandwidth(bandwidth, bandwidth_limited, icoeffsz, icoeffsz);
    }
    
    dataformat = 1;
    fprintf(stderr,"%30s\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b","Inverse SPH transform...      ");
    get_realdata_from_sph_coeffs(rdatax, bandwidth, dataformat, rcoeffsx, icoeffsx);
    get_realdata_from_sph_coeffs(rdatay, bandwidth, dataformat, rcoeffsy, icoeffsy);
    get_realdata_from_sph_coeffs(rdataz, bandwidth, dataformat, rcoeffsz, icoeffsz);

    objects_sphere = create_object( POLYGONS );
    polygons_sphere = get_polygons_ptr(objects_sphere);
    fprintf(stderr,"%30s\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b","Resample surface...           ");
    sample_sphere_from_sph(rdatax, rdatay, rdataz, polygons_sphere, n_triangles, bandwidth);

    if( output_graphics_any_format( surface_filename, ASCII_FORMAT, 1,
                &objects_sphere ) != OK )
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