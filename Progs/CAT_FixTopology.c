/* Rachel Yotter - rachel.yotter@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 */

#include "CAT_FixTopology.h"
#include "CAT_NiftiIO.h"
#include "CAT_SurfaceIO.h"
#include "CAT_Defect.h"

/* the argument table */
ArgvInfo argTable[] = {
  {"-deform", ARGV_CONSTANT, (char *) TRUE, 
     (char *) &do_surface_deform,
     "Deform corrected surface to its uncorrected version in aeras where topology correction was made."},
  { "-bw", ARGV_INT, (char *) 1, 
    (char *) &bw,
    "Bandwidth of coefficients for spherical harmonic expansion." },
  { "-lim", ARGV_INT, (char *) 1, 
    (char *) &lim,
    "Limit bandwidth of spherical harmonic expansion." },
  { "-n", ARGV_INT, (char *) 1, 
    (char *) &n_triangles,
    "Number of triangles for sampled surface." },
  {"-refine_length", ARGV_FLOAT, (char *) TRUE, 
    (char *) &max_refine_length,
     "Maximal length of vertex side after refinement (use negative values for no refinement)."},
  {"-sphere", ARGV_STRING, (char *) 1, (char *) &reparam_file,
     "Sphere object for reparameterization."},
  { "-fill", ARGV_CONSTANT, (char *) TRUE, 
    (char *) &holes,
     "Force filling topology errors."},
  { "-cut", ARGV_CONSTANT, (char *) TRUE, 
    (char *) &handles,
     "Force cutting topology errors."},
  { "-laplace", ARGV_FLOAT, (char *) 1, 
    (char *) &laplace_thresh,
    "Apply 2D laplace filter to minimize irregularities in coordinate maps for SPH. Use 0 for no filtering." },
  { NULL, ARGV_END, NULL, NULL, NULL }
};

void
usage(char *executable)
{
        static char *usage_str = "\n\
Usage: %s surface_file sphere_file output_surface_file\n\
Correct the topology of a brain surface.\n\n\n";

       fprintf(stderr, usage_str, executable);
}


int
main(int argc, char *argv[])
{
        char                 *surface_file, *sphere_file, *output_surface_file;
        File_formats         format;
        int                  n_objects, force;
        polygons_struct      *surface, *sphere;
        object_struct        **surf_objects, **sphere_objects, **objects;

        /* Call ParseArgv */
        if (ParseArgv(&argc, argv, argTable, 0) || argc != 4) {
                usage(argv[0]);
                fprintf(stderr, "       %s -help\n\n", argv[0]);
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
                                      &n_objects, &surf_objects) != OK)
                exit(EXIT_FAILURE);

        /* check that the surface file contains a polyhedron */
        if (n_objects != 1 || get_object_type(surf_objects[0]) != POLYGONS) {
                fprintf(stderr,"Surface file must contain 1 polygons object.\n");
                exit(EXIT_FAILURE);
        }

        /* check that the surface file contains a polyhedron */
        if (handles== 1 && holes==1) {
                fprintf(stderr,"Use either -holes or -handles as option.\n");
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

        /* force to either to assume holes, handles or estimate automated */
        force = 0;
        if (holes)   force = HOLE;
        if (handles) force = HANDLE;
        
        objects = fix_topology_sph(surface, sphere, n_triangles, bw, lim, reparam_file, max_refine_length, 
            do_surface_deform, force, laplace_thresh);

        if (output_graphics_any_format(output_surface_file, ASCII_FORMAT, 1,
                                       objects, NULL) != OK)
                exit(EXIT_FAILURE);

        /* clean up */

        delete_object_list(1, surf_objects);
        delete_object_list(1, sphere_objects);
        delete_object_list(1, objects);

        return(EXIT_SUCCESS);    
}
