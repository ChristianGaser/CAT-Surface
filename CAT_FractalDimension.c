/* Rachel Yotter - rachel.yotter@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id: CAT_FractalDimension.c 140 2009-11-09 12:54:40Z raytrace $
 *
 */

/* Calculate the fractal dimension. */

#include <volume_io/internal_volume_io.h>
#include <volume_io/geometry.h>
#include <bicpl.h>
#include <ParseArgv.h>

#include "CAT_SheetIO.h"
#include "CAT_Map2d.h"
#include "CAT_Blur2d.h"
#include "CAT_Surf.h"
#include "CAT_SPH.h"
#include "CAT_Intersect.h"
#include "CAT_SurfaceIO.h"
#include "CAT_Defect.h"
#include "CAT_Complexity.h"

/* argument defaults */
int maxiters = 30; /* the number of FD's to calculate */
BOOLEAN sph = 1; /* 1=sph, 0=straight resampling */
int n_triangles = 327680; /* # of triangles */
//int n_triangles = 81920; /* # of triangles */
int bw = 1024; /* starting bandwidth, only for SPH */

/* the argument table */
ArgvInfo argTable[] = {
  { "-iters", ARGV_INT, (char *) 1, 
    (char *) &maxiters,
    "Number of iterations." },
  { "-bw", ARGV_INT, (char *) 1, 
    (char *) &bw,
    "Coefficient bandwidth for spherical harmonic expansion (SPH only)." },
  { "-n", ARGV_INT, (char *) 1,
    (char *) &n_triangles,
    "Number of triangles for sampled surface." },
  { "-sph", ARGV_CONSTANT, (char *) TRUE,
    (char *) &sph,
    "Use SPH for FD calculation [default]." },
  { "-nosph", ARGV_CONSTANT, (char *) FALSE,
    (char *) &sph,
    "Use direct resampling for FD calculation." },
  { NULL, ARGV_END, NULL, NULL, NULL }
};

void
usage(char *executable)
{
        static char *usage_str = "\n\
Usage: %s surface.obj sphere.obj out.txt\n\
Calculate the fractal dimension.  Output are surface areas.\n\n\n";

       fprintf(stderr, usage_str, executable);
}

int
main(int argc, char *argv[])
{
        char                 *surface_file, *sphere_file, *output_file;
        File_formats         format;
        int                  n_objects;
        polygons_struct      *surface, *sphere;
        object_struct        **surf_objects, **sphere_objects;
        double               fd;

        /* Call ParseArgv */
        if (ParseArgv(&argc, argv, argTable, 0) || argc != 4) {
                usage(argv[0]);
                fprintf(stderr, "       %s -help\n\n", argv[0]);
                exit(EXIT_FAILURE);
        }

        initialize_argument_processing(argc, argv);

        if (!get_string_argument(NULL, &surface_file) ||
            !get_string_argument(NULL, &sphere_file) ||
            !get_string_argument(NULL, &output_file)) {
                usage(argv[0]);
                exit(EXIT_FAILURE);
        }

        if (input_graphics_any_format(surface_file, &format,
                                      &n_objects, &surf_objects) != OK)
                exit(EXIT_FAILURE);

        /* check that the surface file contains a polyhedron */
        if (n_objects != 1 || get_object_type(surf_objects[0]) != POLYGONS) {
                fprintf(stderr,"Surface file must contain 1 polygons object.\n");
                exit(EXIT_FAILURE);
        }

        /* get a pointer to the surface */
        surface = get_polygons_ptr(surf_objects[0]);
    
        if (input_graphics_any_format(sphere_file, &format,
                                      &n_objects, &sphere_objects) != OK)
                exit(EXIT_FAILURE);

        /* check that the surface file contains a polyhedron */
        if (n_objects != 1 || get_object_type(sphere_objects[0]) != POLYGONS) {
                fprintf(stderr,"Surface file must contain 1 polygons object.\n");
                exit(EXIT_FAILURE);
        }

        /* get a pointer to the surface */
        sphere = get_polygons_ptr(sphere_objects[0]);

        /* check that surface and sphere are same size */
        if (surface->n_items != sphere->n_items) {
                fprintf(stderr,"Surface and sphere must have same size.\n");
                exit(EXIT_FAILURE);
        }

        if (sph) {
                fd = fractal_dimension_sph(surface, sphere, output_file,
                                           n_triangles);
        } else {
                fd = fractal_dimension(surface, sphere, maxiters, output_file);
        }

        printf("fd = %f\n", fd);

        /* clean up */
        delete_object_list(1, surf_objects);
        delete_object_list(1, sphere_objects);

        return(EXIT_SUCCESS);    
}
