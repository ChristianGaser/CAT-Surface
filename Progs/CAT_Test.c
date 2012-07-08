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
        polygons_struct  *source, *source_sphere, *warped_sphere;
        polygons_struct  *target, *target_sphere, *warped;
        int              x, y, i, j, l, m;
        int              xy_size, n_objects, bandwidth, bandwidth2, cutoff, n_points;
        int              size_map[2], n_triangles, n_polygons, poly;
        double           u, v, u2, v2, xp, yp, xm, ym, weight, *ux, *vy;
        object_struct    **objects, *object, *object2;
        int              dataformat;
        Point            sphere_point, new_point, closest;    

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

        /* get a pointer to the surface */
        target = get_polygons_ptr(objects[0]);

        if (input_graphics_any_format(target_sphere_file, &format,
                                      &n_objects, &objects) != OK)
                exit(EXIT_FAILURE);
        
        /* check that the surface file contains a polyhedron */
        if (n_objects != 1 || get_object_type(objects[0]) != POLYGONS) {
                print("Surface file must contain 1 polygons object.\n");
                exit(EXIT_FAILURE);
        }

        /* get a pointer to the surface */
        target_sphere = get_polygons_ptr(objects[0]);

        if (input_graphics_any_format(source_file, &format,
                                      &n_objects, &objects) != OK)
                exit(EXIT_FAILURE);
                
        /* check that the surface file contains a polyhedron */
        if (n_objects != 1 || get_object_type(objects[0]) != POLYGONS) {
                print("Surface file must contain 1 polygons object.\n");
                exit(EXIT_FAILURE);
        }
        
        /* get a pointer to the surface */
        source = get_polygons_ptr(objects[0]);
    
        if (input_graphics_any_format(source_sphere_file, &format,
                                      &n_objects, &objects) != OK)
                exit(EXIT_FAILURE);
        
        /* check that the surface file contains a polyhedron */
        if (n_objects != 1 || get_object_type(objects[0]) != POLYGONS) {
                print("Surface file must contain 1 polygons object.\n");
                exit(EXIT_FAILURE);
        }
        
        /* get a pointer to the surface */
        source_sphere = get_polygons_ptr(objects[0]);

        /* create sphere */
        object = create_object(POLYGONS);
        warped = get_polygons_ptr(object);
        copy_polygons(source,warped);
        
        object2 = create_object(POLYGONS);
        warped_sphere = get_polygons_ptr(object2);
        copy_polygons(source,warped_sphere);

        surf_to_sphere(warped, 3);
        surf_to_sphere(warped_sphere, 3);

        translate_to_center_of_mass(warped);
        translate_to_center_of_mass(warped_sphere);
        
        /* set radius to 1 */
        for (i = 0; i < source_sphere->n_points; i++) 
                set_vector_length(&source_sphere->points[i], 1.0);
        for (i = 0; i < target_sphere->n_points; i++) 
                set_vector_length(&target_sphere->points[i], 1.0);

        create_polygons_bintree(warped, round((double) warped->n_items *
                                              BINTREE_FACTOR));
        create_polygons_bintree(target, round((double) target->n_items *
                                              BINTREE_FACTOR));

        for (i = 0; i < source->n_points; i++) {

                poly = find_closest_polygon_point(&warped->points[i], warped_sphere, &closest);
//                n_points = get_polygon_points(target_sphere, poly, &sphere_point);
                sphere_point = target_sphere->points[target_sphere->indices[POINT_INDEX(target_sphere->end_indices,poly,0)]];
                
                fprintf(stderr,"%d/%d: %3.4f %3.4f %3.4f\n",i,source->n_points,Point_x(source_sphere->points[i])-Point_x(sphere_point),Point_y(source_sphere->points[i])-Point_y(sphere_point),Point_z(source_sphere->points[i])-Point_z(sphere_point));
                warped_sphere->points[i] = sphere_point;
                
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
