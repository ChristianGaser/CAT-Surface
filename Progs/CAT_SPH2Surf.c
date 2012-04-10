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

/* argument defaults */
int bandwidth = 256;
int bandwidth_limited = 0;
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
Usage: %s  SPH.dat output.obj [bandwidth_limited n_points]\n\
Use spherical harmonic coefficients to create a surface.\n\n\n";

        fprintf(stderr, usage_str, executable);
}

int
main(int argc, char *argv[])
{
        FILE                 *fp;
        char                 line[256];
        char                 *surface_file, *SPH_file;
        int                  dataformat; /* = 0 -> samples are complex,
                                          * = 1 -> samples are real */
        int                  i, bandwidth2, n_dims;
        double               *rcx, *icx, *rcy, *icy, *rcz, *icz;
        double               *rdatax, *rdatay, *rdataz;
        polygons_struct      *polygons_sphere;
        object_struct        *objects_sphere;
   
        /* Call ParseArgv */
        if (ParseArgv(&argc, argv, argTable, 0) || argc != 3) {
                usage(argv[0]);
                fprintf(stderr, "       %s -help\n\n", argv[0]);
                exit(EXIT_FAILURE);
        }

        initialize_argument_processing(argc, argv);

        if (!get_string_argument(NULL, &SPH_file) ||
            !get_string_argument(NULL, &surface_file)) {
                usage(argv[0]);
                exit(EXIT_FAILURE);
        }
     
        // read coefficients
        if ((fp = fopen(SPH_file, "rb")) == NULL) {
                fprintf(stderr, "Error opening file %s.\n", SPH_file);
                exit(EXIT_FAILURE);
        }

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

        if (fread(rcx, sizeof(double), size, fp) != size) {
                fprintf(stderr, "Error reading data.\n");
                exit(EXIT_FAILURE);
        }
        if (fread(icx, sizeof(double), size, fp) != size) {
                fprintf(stderr, "Error reading data.\n");
                exit(EXIT_FAILURE);
        }
        if (fread(rcy, sizeof(double), size, fp) != size) {
                fprintf(stderr, "Error reading data.\n");
                exit(EXIT_FAILURE);
        }
        if (fread(icy, sizeof(double), size, fp) != size) {
                fprintf(stderr, "Error reading data.\n");
                exit(EXIT_FAILURE);
        }
        if (fread(rcz, sizeof(double), size, fp) != size) {
                fprintf(stderr, "Error reading data.\n");
                exit(EXIT_FAILURE);
        }
        if (fread(icz, sizeof(double), size, fp) != size) {
                fprintf(stderr, "Error reading data.\n");
                exit(EXIT_FAILURE);
        }

        fclose(fp);

        bandwidth2 = 2*bandwidth;

        rdatax   = (double *) malloc(sizeof(double) * bandwidth2*bandwidth2);
        rdatay   = (double *) malloc(sizeof(double) * bandwidth2*bandwidth2);
        rdataz   = (double *) malloc(sizeof(double) * bandwidth2*bandwidth2);

        if (bandwidth_limited > 0 && bandwidth_limited < bandwidth) {
                fprintf(stderr,"%30s\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b","Limit bandwidth...            ");
                limit_bandwidth(bandwidth, bandwidth_limited, rcx, rcx);
                limit_bandwidth(bandwidth, bandwidth_limited, rcy, rcy);
                limit_bandwidth(bandwidth, bandwidth_limited, rcz, rcz);
                limit_bandwidth(bandwidth, bandwidth_limited, icx, icx);
                limit_bandwidth(bandwidth, bandwidth_limited, icy, icy);
                limit_bandwidth(bandwidth, bandwidth_limited, icz, icz);
        }
    
        dataformat = 1;
        fprintf(stderr,"%30s\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b","Inverse SPH transform...      ");
        get_realdata_from_sph_coeffs(rdatax, bandwidth, dataformat, rcx, icx);
        get_realdata_from_sph_coeffs(rdatay, bandwidth, dataformat, rcy, icy);
        get_realdata_from_sph_coeffs(rdataz, bandwidth, dataformat, rcz, icz);

        objects_sphere = create_object( POLYGONS );
        polygons_sphere = get_polygons_ptr(objects_sphere);
        fprintf(stderr,"%30s\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b","Resample surface...           ");
        sample_sphere_from_sph(rdatax, rdatay, rdataz, polygons_sphere,
                               n_triangles, NULL, bandwidth);

        if (output_graphics_any_format(surface_file, ASCII_FORMAT, 1,
                                       &objects_sphere) != OK)
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
    
        fprintf(stderr,"%30s\n","Done                          ");

        return(EXIT_SUCCESS);    
}
