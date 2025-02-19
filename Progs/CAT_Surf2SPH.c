/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 *
 */

#include <bicpl.h>
#include <fftw3.h>

#include "CAT_Map.h"
#include "CAT_Smooth.h"
#include "CAT_SPH.h"
#include "CAT_SurfaceIO.h"

void
usage(char *executable)
{
    static char *usage_str = "\n\
Usage: %s  surface_file sphere_file SPH.txt bandwidth\n\
    Extract spherical harmonic coefficients of a surface using\n\
    2*bandwidth x 2*bandwidth sample points.\n\
\n\n";

    fprintf(stderr, usage_str, executable);
}

int
main(int argc, char *argv[])
{
    char         *surface_file, *sphere_file, *SPH_file;
    File_formats     format;
    polygons_struct    *polygons, *polygons_sphere;
    int          i, j, k;
    int          n_items, bandwidth, bandwidth2;
    double         *rcx, *icx, *rcy, *icy, *rcz, *icz;
    double         *rdatax, *rdatay, *rdataz;
    object_struct    **objects;
    int          n_objects, n_dims = 3;
    int          dataformat; /* = 0 -> samples are complex, */
                     /* = 1 -> samples are real */

    initialize_argument_processing(argc, argv);

    if (!get_string_argument(NULL, &surface_file) ||
      !get_string_argument(NULL, &sphere_file) ||
      !get_string_argument(NULL, &SPH_file)) {
        usage(argv[0]);
        exit(EXIT_FAILURE);
    }

    get_int_argument(256, &bandwidth);
    bandwidth2 = 2 * bandwidth;

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
    polygons_sphere = get_polygons_ptr(objects[0]);

    /* check that surface and sphere are same size */
    if (polygons->n_items != polygons_sphere->n_items) {
        printf("Surface and sphere must have same size.\n");
        exit(EXIT_FAILURE);
    }

    rdatax   = (double *) malloc(sizeof(double) * bandwidth2*bandwidth2);
    rdatay   = (double *) malloc(sizeof(double) * bandwidth2*bandwidth2);
    rdataz   = (double *) malloc(sizeof(double) * bandwidth2*bandwidth2);
    rcx    = (double *) malloc(sizeof(double) * bandwidth*bandwidth);
    rcy    = (double *) malloc(sizeof(double) * bandwidth*bandwidth);
    rcz    = (double *) malloc(sizeof(double) * bandwidth*bandwidth);
    icx    = (double *) malloc(sizeof(double) * bandwidth*bandwidth);
    icy    = (double *) malloc(sizeof(double) * bandwidth*bandwidth);
    icz    = (double *) malloc(sizeof(double) * bandwidth*bandwidth);

    printf("%30s\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b", "Sample spherical coordinates..");
    get_equally_sampled_coords_of_polygon(polygons, polygons_sphere,
                        bandwidth, 
                        rdatax, rdatay, rdataz);

    /* dataformat indicates real data */
    dataformat = 1;
    printf("%30s\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b", "Forward SPH transform...    ");
    get_sph_coeffs_of_realdata(rdatax, bandwidth, dataformat, rcx, icx);
    get_sph_coeffs_of_realdata(rdatay, bandwidth, dataformat, rcy, icy);
    get_sph_coeffs_of_realdata(rdataz, bandwidth, dataformat, rcz, icz);

    if (write_SPHxyz(SPH_file, bandwidth, rcx, rcy, rcz,
             icx, icy, icz) != OK)
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

    printf("%30s\n", "Done              ");

    return(EXIT_SUCCESS);
}
