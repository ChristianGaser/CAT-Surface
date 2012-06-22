/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 *
 */

#include    <bicpl.h>
#include    <float.h>

#include    "CAT_Surf.h"
#include    "CAT_SPH.h"
    
void
usage(char *executable)
{
        static  char  *usage_str = "\n\
Usage: %s  source.obj source_sphere.obj target.obj target_sphere.obj warped.obj [cutoff]\n\
        Use spherical harmonic coefficients to warp a surface to a given template.\n\
\n\n";

        fprintf(stderr, usage_str, executable);
}
int
main(int argc, char *argv[])
{
        char             *source_file, *source_sphere_file, *warped_sphere_file;
        char             *target_file, *target_sphere_file, *warped_file;
        File_formats     format;
        object_struct    **object_list;
        polygons_struct  *polygons, *polygonsIn;
        Point            *smooth_pts;
        Real             fwhm;
        double           *sulc_depth, *faces, *vertices;

        initialize_argument_processing(argc, argv);

        if (!get_string_argument(NULL, &source_file) ||
            !get_string_argument(NULL, &source_sphere_file) ||
            !get_string_argument(NULL, &target_file) ||
            !get_string_argument(NULL, &target_sphere_file) ||    
            !get_string_argument(NULL, &warped_file) ||    
            !get_string_argument(NULL, &warped_sphere_file)) {
                usage(argv[0]);
                exit(EXIT_FAILURE);
        }

        get_int_argument(32, &cutoff);
        bandwidth = 256;
        
        bandwidth2  = bandwidth*2;
        n_triangles = 81920;

        if (input_graphics_any_format(target_file, &format,
                                      &n_objects, &objects) != OK)
                exit(EXIT_FAILURE);

        /* check that the surface file contains a polyhedron */
        if (n_objects != 1 || get_object_type(objects[0]) != POLYGONS) {
                print("Surface file must contain 1 polygons object.\n");
                exit(EXIT_FAILURE);
        }

        polygons = get_polygons_ptr(object_list[0]);
        polygonsIn = get_polygons_ptr(create_object(POLYGONS));
        copy_polygons(polygons, polygonsIn);

        if (input_graphics_any_format(target_sphere_file, &format,
                                      &n_objects, &objects) != OK)
                exit(EXIT_FAILURE);
        
        /* check that the surface file contains a polyhedron */
        if (n_objects != 1 || get_object_type(objects[0]) != POLYGONS) {
                print("Surface file must contain 1 polygons object.\n");
                exit(EXIT_FAILURE);
        }

 //       get_all_polygon_point_neighbours(polygons, &n_neighbours, &neighbours);

fwhm = 5;
        smooth_heatkernel(polygons, NULL, fwhm);

        compute_polygon_normals(polygons);

        for (i = 0; i < polygons->n_points; i++) {
                sulc_depth[i] = ((Point_x(polygons->points[i]) - Point_x(polygonsIn->points[i]))*Point_x(polygons->normals[i]) +
                                 (Point_y(polygons->points[i]) - Point_y(polygonsIn->points[i]))*Point_y(polygons->normals[i]) +
                                 (Point_z(polygons->points[i]) - Point_z(polygonsIn->points[i]))*Point_z(polygons->normals[i]));
        }

        for (i = 0; i < warped_sphere->n_points; i++) 
                set_vector_length(&warped_sphere->points[i], 1.0);

        if (output_graphics_any_format(warped_file, ASCII_FORMAT,
                                       1, &object) != OK)
                exit(EXIT_FAILURE);

        if (output_graphics_any_format(warped_sphere_file, ASCII_FORMAT,
                                       1, &object2) != OK)
                exit(EXIT_FAILURE);

        delete_object_list(n_objects, objects);
               
        return(EXIT_SUCCESS);
}
