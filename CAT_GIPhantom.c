/* Rachel Yotter - rachel.yotter@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 */

#include <volume_io/internal_volume_io.h>
#include <bicpl.h>
#include <bicpl/deform.h>
#include <ParseArgv.h>
#include <math.h>

#include "CAT_Map2d.h"
#include "CAT_Blur2d.h"
#include "CAT_SurfaceIO.h"

/* argument defaults */
int ntheta = 4;
int gtheta = 0;
int nphi = 0;
int gphi = 0;
double amplitude = 0.05;
double radius = 65.0;
BOOLEAN gamp = 0;
int n_triangles = 327680;

#define PINF  1.7976931348623157e+308 /* for doubles */
#define NINF -1.7976931348623157e+308 /* for doubles */

/* the argument table */
ArgvInfo argTable[] = {
  { "-ntheta", ARGV_INT, (char *) 1, 
    (char *) &ntheta,
    "Number of oscillations (or rings) wrt theta." },
  { "-gtheta", ARGV_INT, (char *) 1, 
    (char *) &gtheta,
    "Size of theta gradient (frequencies from ntheta to ntheta*gtheta." },
  { "-nphi", ARGV_INT, (char *) 1, 
    (char *) &nphi,
    "Number of oscillations (or rings) wrt phi." },
  { "-gphi", ARGV_INT, (char *) 1, 
    (char *) &gphi,
    "Size of phi gradient (frequencies from nphi to nphi*gphi." },
  { "-amplitude", ARGV_FLOAT, (char *) 1, 
    (char *) &amplitude,
    "Amplitude of oscillations, between 0 and 1." },
  { "-gamp", ARGV_CONSTANT, (char *) 1, 
    (char *) &gamp,
    "Flag to use amplitude gradient along y-axis." },
  { "-radius", ARGV_FLOAT, (char *) 1, 
    (char *) &radius,
    "Radius of output surface." },
  { "-ntriangles", ARGV_INT, (char *) 1, 
    (char *) &n_triangles,
    "Number of triangles for output surface." },
  { NULL, ARGV_END, NULL, NULL, NULL }
};

void
usage(char *executable)
{
        static char *usage_str = "\n\
Usage: %s output.obj\n\
Generate a GI phantom with oscillations along theta and/or phi.\n\n\n";

       fprintf(stderr, usage_str, executable);
}

int
main(int argc, char *argv[])
{
        char                 *output_file;
        File_formats         format;
        int                  n_objects, n_polygons, p;
        polygons_struct      *surface;
        object_struct        **objects;
        int                  *n_neighbours, **neighbours;
        Point                centre;
        double               a, theta, phi, x, y, z, gx, gz;

        /* Call ParseArgv */
        if (ParseArgv(&argc, argv, argTable, 0) || argc != 2) {
                usage(argv[0]);
                fprintf(stderr, "       %s -help\n\n", argv[0]);
                exit(EXIT_FAILURE);
        }

        initialize_argument_processing(argc, argv);

        if (!get_string_argument(NULL, &output_file)) {
                usage(argv[0]);
                exit(EXIT_FAILURE);
        }

        if (ntheta < 0) {
                fprintf(stderr,"Error: theta oscillations must be > 0.\n");
                exit(EXIT_FAILURE);
        }
        if (nphi < 0) {
                fprintf(stderr,"Error: phi oscillations must be > 0.\n");
                exit(EXIT_FAILURE);
        }
        if (amplitude < 0 || amplitude > 1) {
                fprintf(stderr,"Error: amplitude must be in the range (0,1].\n");
                exit(EXIT_FAILURE);
        }

        /* check tetrahedral topology */
        /* best areal distribution is achieved for 20 edges */
        n_polygons = n_triangles;

        while (n_polygons % 4 == 0)
                n_polygons /= 4;

        if (n_polygons != 5) {
                fprintf(stderr,"Warning: Number of triangles %d is",
                               n_triangles);
                fprintf(stderr," not recommend because\ntetrahedral");
                fprintf(stderr," topology is not optimal.\nPlease try");
                fprintf(stderr," 20*(4*x) triangles (e.g. 81920).\n");
        }

        objects = (object_struct **) malloc(sizeof(object_struct *));
        *objects = create_object(POLYGONS);
        surface = get_polygons_ptr(*objects);
        initialize_polygons(surface, WHITE, NULL);
   
        fill_Point(centre, 0.0, 0.0, 0.0);
        create_tetrahedral_sphere(&centre, 1.0, 1.0, 1.0, n_triangles, surface);
        compute_polygon_normals(surface);

        for (p = 0; p < surface->n_points; p++) {
                x = Point_x(surface->normals[p]);
                y = Point_y(surface->normals[p]);
                z = Point_z(surface->points[p]);

                theta = acos(x);
                phi = acos(z);

                if (gamp) {
                        if (y != 0)
                                a = amplitude *
                                    (asin(y / sqrt(x*x + y*y))/M_PI + 0.5);
                        else
                                a = amplitude * 0.5;
                } else {
                        a = amplitude; /* default value */
                }

                /* set the gradients */
                gx = 1 + gtheta * (theta/M_PI);
                gz = 1 + gphi * (phi/M_PI);

                a = a * ((cos(theta * gx * 2 * ntheta) +
                          cos(phi * gz * 2 * nphi))/2) + 1;
                set_vector_length(&surface->points[p], radius*a);
        }

        compute_polygon_normals(surface);

        if (output_graphics_any_format(output_file, ASCII_FORMAT, 1,
                                       objects) != OK)
                exit(EXIT_FAILURE);

        /* clean up */
        delete_object_list(1, objects);

        return(EXIT_SUCCESS);    
}
