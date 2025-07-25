/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 *
 */

#include <ParseArgv.h>

#include "CAT_SurfaceIO.h"
#include "CAT_Map.h"
#include "CAT_Smooth.h"
#include "CAT_SPH.h"

/* argument defaults */
int bandwidth = 256;
int bandwidth_limited = 128;
int n_triangles = 81920;

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

void
usage(char *executable)
{
    static char *usage_str = "\n\
Usage: %s surface_file sphere_file output_surface_file\n\
Use spherical harmonic coefficients to create an equally sampled surface.\n\
\n\n";

     fprintf(stderr, usage_str, executable);
}

int
main(int argc, char *argv[])
{
    char         *surface_file, *sphere_file, *output_surface_file;
    File_formats     format;
    int          dataformat; /* = 0 -> samples are complex,
                      * = 1 -> samples are real */
    int          i, j, bandwidth2, n_objects;
    double         *rcx, *icx, *rcy, *icy, *rcz, *icz;
    double         *rdatax, *rdatay, *rdataz;
    double         r, value;
    polygons_struct    *polygons, *sphere, *polygons_output;
    object_struct    **objects, *objects_output;

    /* Call ParseArgv */
    if (ParseArgv(&argc, argv, argTable, 0) || argc != 4) {
        usage(argv[0]);
        fprintf(stderr, "   %s -help\n\n", argv[0]);
        exit(EXIT_FAILURE);
    }

    initialize_argument_processing(argc, argv);

    if (!get_string_argument(NULL, &surface_file) ||
      !get_string_argument(NULL, &sphere_file) ||
      !get_string_argument(NULL, &output_surface_file)) {
        usage(argv[0]);
        exit(EXIT_FAILURE);
    }
   
    if (input_graphics_any_format(surface_file, &format,
                    &n_objects, &objects) != OK)
        exit(EXIT_FAILURE);

    /* check that the surface file contains a polyhedron */
    if (n_objects != 1 || get_object_type(objects[0]) != POLYGONS) {
        printf("Surface file must contain 1 polygons object.\n");
        exit(EXIT_FAILURE);
    }
    /* get a pointer to the surface */
    polygons = get_polygons_ptr(objects[0]);
  
    if (input_graphics_any_format(sphere_file, &format,
                    &n_objects, &objects) != OK)
        exit(EXIT_FAILURE);

    /* check that the surface file contains a polyhedron */
    if (n_objects != 1 || get_object_type(objects[0]) != POLYGONS) {
        printf("Surface file must contain 1 polygons object.\n");
        exit(EXIT_FAILURE);
    }
    /* get a pointer to the surface */
    sphere = get_polygons_ptr(objects[0]);

    /* check that surface and sphere are same size */
    if (polygons->n_items != sphere->n_items) {
        printf("Surface and sphere must have same size.\n");
        exit(EXIT_FAILURE);
    }
  
    bandwidth2 = 2*bandwidth;

    rdatax = (double *) malloc(sizeof(double) * bandwidth2 * bandwidth2);
    rdatay = (double *) malloc(sizeof(double) * bandwidth2 * bandwidth2);
    rdataz = (double *) malloc(sizeof(double) * bandwidth2 * bandwidth2);
    rcx    = (double *) malloc(sizeof(double) * bandwidth * bandwidth);
    rcy    = (double *) malloc(sizeof(double) * bandwidth * bandwidth);
    rcz    = (double *) malloc(sizeof(double) * bandwidth * bandwidth);
    icx    = (double *) malloc(sizeof(double) * bandwidth * bandwidth);
    icy    = (double *) malloc(sizeof(double) * bandwidth * bandwidth);
    icz    = (double *) malloc(sizeof(double) * bandwidth * bandwidth);

    printf("Samp SPH coords..");
    get_equally_sampled_coords_of_polygon(polygons, sphere, bandwidth,
                        rdatax, rdatay, rdataz);

    /* dataformat indicates real data */
    dataformat = 1;
    printf("Fwd SPH xform..");
    get_sph_coeffs_of_realdata(rdatax, bandwidth, dataformat, rcx, icx);
    get_sph_coeffs_of_realdata(rdatay, bandwidth, dataformat, rcy, icy);
    get_sph_coeffs_of_realdata(rdataz, bandwidth, dataformat, rcz, icz);

    if (bandwidth_limited > 0 && bandwidth_limited < bandwidth) {
        printf("Limit BW..");
        butterworth_filter(bandwidth, bandwidth_limited, rcx, rcx);
        butterworth_filter(bandwidth, bandwidth_limited, rcy, rcy);
        butterworth_filter(bandwidth, bandwidth_limited, rcz, rcz);
        butterworth_filter(bandwidth, bandwidth_limited, icx, icx);
        butterworth_filter(bandwidth, bandwidth_limited, icy, icy);
        butterworth_filter(bandwidth, bandwidth_limited, icz, icz);
    }
  
    printf("Inv SPH xform..");
    get_realdata_from_sph_coeffs(rdatax, bandwidth, dataformat, rcx, icx);
    get_realdata_from_sph_coeffs(rdatay, bandwidth, dataformat, rcy, icy);
    get_realdata_from_sph_coeffs(rdataz, bandwidth, dataformat, rcz, icz);

    objects_output = create_object(POLYGONS);
    polygons_output = get_polygons_ptr(objects_output);

    printf("Resamp surf..");
    sample_sphere_from_sph(rdatax, rdatay, rdataz,
                 polygons_output, n_triangles, NULL, bandwidth);

    if (output_graphics_any_format(output_surface_file, ASCII_FORMAT, 1,
                     &objects_output, NULL) != OK)
        exit(EXIT_FAILURE);

    /* clean up */
    free(rcx);
    free(rcy);
    free(rcz);
    free(icx);
    free(icy);
    free(icz);
    free(rdatax);
    free(rdatay);
    free(rdataz);
  
    printf("Done\n");

    return(EXIT_SUCCESS);  
}
