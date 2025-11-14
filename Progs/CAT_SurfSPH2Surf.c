/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 *
 */

#include <bicpl.h>
#include <ParseArgv.h>

#include "CAT_SurfaceIO.h"
#include "CAT_Map.h"
#include "CAT_Smooth.h"
#include "CAT_SPH.h"
#include "CAT_Surf.h"
#include "CAT_SafeAlloc.h"

/* argument defaults */
int bandwidth = 256;
int bandwidth_limited = 0;
int n_triangles = 81920;
int correct_bb = 0;

/* the argument table */
ArgvInfo argTable[] = {
  { "-bw", ARGV_INT, (char *) 1, 
  (char *) &bandwidth,
  "Bandwidth of coefficients for spherical harmonic expansion." },
  { "-lim", ARGV_INT, (char *) 1, 
  (char *) &bandwidth_limited,
  "Limit bandwidth of spherical harmonic expansion." },
  { "-correct", ARGV_CONSTANT, (char *) TRUE, 
  (char *) &correct_bb,
  "Correct bounding box according to unlimited spherical harmonic expansion. This accounts for the smaller bounding box after applying bandwidth limitation." },
  { "-n", ARGV_INT, (char *) 1, 
  (char *) &n_triangles,
  "Number of triangles for sampled surface." },
  { NULL, ARGV_END, NULL, NULL, NULL }
};

void
usage(char *executable)
{
    static char *usage_str = "\n\
Usage: %s  SPH.txt output_surface_file [-bw bandwidth -lim bandwidth_limited -correct -n n_points]\n\
Use spherical harmonic coefficients to create a surface.\n\n\n";

    fprintf(stderr, usage_str, executable);
}

int
main(int argc, char *argv[])
{
    FILE         *fp;
    char         line[256];
    char         *surface_file, *SPH_file;
    int          dataformat; /* = 0 -> samples are complex,
                      * = 1 -> samples are real */
    int          i, bandwidth2, n_dims;
    double         *rcx, *icx, *rcy, *icy, *rcz, *icz;
    double         *rdatax, *rdatay, *rdataz;
    float        rx,ix,ry,iy,rz,iz;
    polygons_struct    *polygons_bw, *polygons_full;
    object_struct    *objects_polygon_bw, *objects_polygon_full;
   
    /* Call ParseArgv */
    if (ParseArgv(&argc, argv, argTable, 0) || argc != 3) {
        usage(argv[0]);
        fprintf(stderr, "   %s -help\n\n", argv[0]);
        exit(EXIT_FAILURE);
    }

    initialize_argument_processing(argc, argv);

    if (!get_string_argument(NULL, &SPH_file) ||
      !get_string_argument(NULL, &surface_file)) {
        usage(argv[0]);
        exit(EXIT_FAILURE);
    }
   
    // read coefficients
    fp = SAFE_FOPEN(SPH_file, "rb");

    fgets(line, 256, fp);
    if (strncmp(line,"SPH", 3)) {
        fprintf(stderr, "Wrong file format.\n");
        exit(EXIT_FAILURE);
    }
    fgets(line, 256, fp);
    sscanf(line, "%d %d %d", &bandwidth, &n_dims, &dataformat);

    if (dataformat != 1) {
        fprintf(stderr, "Wrong dataformat: Data should be real and not complex.\n");
        exit(EXIT_FAILURE);
    }

    if (n_dims != 3) {
        fprintf(stderr, "Data dimension should be 3.\n");
        exit(EXIT_FAILURE);
    }

    int size = bandwidth * bandwidth;
    rcx = (double *) malloc(sizeof(double) * size);
    rcy = (double *) malloc(sizeof(double) * size);
    rcz = (double *) malloc(sizeof(double) * size);
    icx = (double *) malloc(sizeof(double) * size);
    icy = (double *) malloc(sizeof(double) * size);
    icz = (double *) malloc(sizeof(double) * size);

    for (i = 0; i < size; i++) {
        fgets(line, 256, fp);
        sscanf(line, "%f %f %f %f %f %f", &rx, &ix, &ry, &iy, &rz, &iz);
        rcx[i] = rx; icx[i] = ix; rcy[i] = ry; icy[i] = iy; rcz[i] = rz; icz[i] = iz; 
    }

    fclose(fp);

    dataformat = 1;
    bandwidth2 = 2*bandwidth;

    rdatax   = (double *) malloc(sizeof(double) * bandwidth2*bandwidth2);
    rdatay   = (double *) malloc(sizeof(double) * bandwidth2*bandwidth2);
    rdataz   = (double *) malloc(sizeof(double) * bandwidth2*bandwidth2);
    
    if (correct_bb) {
        printf("%30s\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b","Inverse SPH transform full...    ");
        get_realdata_from_sph_coeffs(rdatax, bandwidth, dataformat, rcx, icx);
        get_realdata_from_sph_coeffs(rdatay, bandwidth, dataformat, rcy, icy);
        get_realdata_from_sph_coeffs(rdataz, bandwidth, dataformat, rcz, icz);
    
        objects_polygon_full = create_object( POLYGONS );
        polygons_full = get_polygons_ptr(objects_polygon_full);
        printf("%30s\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b","Resample full surface...       ");
        sample_sphere_from_sph(rdatax, rdatay, rdataz, polygons_full,
                     n_triangles, NULL, bandwidth);
    
    }

    if (bandwidth_limited > 0 && bandwidth_limited < bandwidth) {
        printf("%30s\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b","Limit bandwidth...      ");
        limit_bandwidth(bandwidth, bandwidth_limited, rcx, rcx);
        limit_bandwidth(bandwidth, bandwidth_limited, rcy, rcy);
        limit_bandwidth(bandwidth, bandwidth_limited, rcz, rcz);
        limit_bandwidth(bandwidth, bandwidth_limited, icx, icx);
        limit_bandwidth(bandwidth, bandwidth_limited, icy, icy);
        limit_bandwidth(bandwidth, bandwidth_limited, icz, icz);
    }
  
    printf("%30s\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b","Inverse SPH transform...    ");
    get_realdata_from_sph_coeffs(rdatax, bandwidth, dataformat, rcx, icx);
    get_realdata_from_sph_coeffs(rdatay, bandwidth, dataformat, rcy, icy);
    get_realdata_from_sph_coeffs(rdataz, bandwidth, dataformat, rcz, icz);

    objects_polygon_bw = create_object( POLYGONS );
    polygons_bw = get_polygons_ptr(objects_polygon_bw);
    printf("%30s\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b","Resample surface...       ");
    sample_sphere_from_sph(rdatax, rdatay, rdataz, polygons_bw,
                 n_triangles, NULL, bandwidth);

    if (correct_bb) {
        printf("%30s\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b","Correct bounding box...       ");
        correct_bounds_to_target_with_scaling(polygons_bw, polygons_full);
    }
    
    if (output_graphics_any_format(surface_file, ASCII_FORMAT, 1,
                     &objects_polygon_bw, NULL) != OK)
        exit(EXIT_FAILURE);

    free(rcx);
    free(rcy);
    free(rcz);
    free(icx);
    free(icy);
    free(icz);
    free(rdatax);
    free(rdatay);
    free(rdataz);
  
    printf("%30s\n","Done              ");

    return(EXIT_SUCCESS);  
}
