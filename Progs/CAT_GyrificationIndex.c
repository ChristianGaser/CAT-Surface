/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id: CAT_GyrificationsIndex.c.c 261 2012-04-10 10:24:11Z gaser $
 *
 */

/* Calculate the gyrification index using a convex hull. */

#include <bicpl.h>
#include <ParseArgv.h>

#include "CAT_Map.h"
#include "CAT_Smooth.h"
#include "CAT_Surf.h"
#include "CAT_SPH.h"
#include "CAT_ConvexHull.h"
#include "CAT_SurfaceIO.h"

/* argument defaults */
int n_triangles = 327680; /* # of triangles */
BOOLEAN smooth = 1; /* smooth the values */
char *reparam_file = NULL;

/* the argument table */
ArgvInfo argTable[] = {
  {"-sphere", ARGV_STRING, (char *) 1, (char *) &reparam_file,
     "Sphere object for reparameterization."},
  { "-n", ARGV_INT, (char *) 1,
    (char *) &n_triangles,
    "Number of triangles for sampled surface." },
  { "-nosmooth", ARGV_CONSTANT, (char *) FALSE,
    (char *) &smooth,
    "Turn off smoothing of FD values." },
  { NULL, ARGV_END, NULL, NULL, NULL }
};

void
usage(char *executable)
{
        static char *usage_str = "\n\
Usage: %s surface.obj sphere.obj out.txt\n\
Calculate the gyrification index (GI) using a convex hull.  Output are local GI values.\n\n\n";

       fprintf(stderr, usage_str, executable);
}

int
main(int argc, char *argv[])
{
        char                 *surface_file, *sphere_file, *output_file;
        File_formats         format;
        int                  n_objects, i;
        polygons_struct      *surface, *sphere, *reparam;
        object_struct        **objects, **sphere_objects, **reparam_objects;
        double               gi;

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
                                      &n_objects, &objects) != OK)
                exit(EXIT_FAILURE);

        /* check that the surface file contains a polyhedron */
        if (n_objects != 1 || get_object_type(objects[0]) != POLYGONS) {
                fprintf(stderr,"Surface file must contain 1 polygons object.\n");
                exit(EXIT_FAILURE);
        }

        /* get a pointer to the surface */
        surface = get_polygons_ptr(objects[0]);
    
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
        if (reparam_file != NULL) {
                if (input_graphics_any_format(reparam_file, &format, &n_objects,
                                              &reparam_objects) != OK)
                        exit(EXIT_FAILURE);
                reparam = get_polygons_ptr(*reparam_objects);
                n_triangles = reparam->n_items;
        } else {
                reparam_objects = (object_struct **)
                                  malloc(sizeof(object_struct *));
                *reparam_objects = create_object(POLYGONS);
                reparam = get_polygons_ptr(*reparam_objects);
                copy_polygons(sphere, reparam);
                compute_polygon_normals(reparam);
        }

        /* center and scale the reparameterizing sphere */
        translate_to_center_of_mass(reparam);
        for (i = 0; i < reparam->n_points; i++)
                set_vector_length(&reparam->points[i], 1.0);

        calculate_convex_hull(surface);

//        gi = gyrification_index_sph(surface, sphere, output_file,
//                                           n_triangles, reparam, smooth);

        printf("global GI: %f\n", gi);

        /* clean up */
        delete_object_list(1, objects);
        delete_object_list(1, sphere_objects);

        return(EXIT_SUCCESS);
}
