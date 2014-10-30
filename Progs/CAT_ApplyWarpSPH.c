/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 *
 */

#include <bicpl.h>

#include "CAT_Surf.h"
#include "CAT_SurfaceIO.h"

void
usage(char *executable)
{
        char *usage_str = "\n\
Usage: %s  sphere_file1 sphere_file2 input_values_file output_values_file [n_triangles]\n\
Resamples a spherical inflated surface to a sphere.\n\
\n\n";

        fprintf(stderr, usage_str, executable);
}

int
main(int argc, char *argv[])
{
        char             *surface_file, *sphere_file;
        char             *input_values_file, *output_values_file;
        File_formats     format;
        FILE             *fp;
        int              n_objects, n_objects_src_sphere, point;
        int              n_triangles, n_polygons;
        int              i, j, k;
        int              poly, n_points, n_intersections, n_values;
        int              *n_neighbours, **neighbours;
        Point            centre, point_on_src_sphere;
        Point            poly_points[MAX_POINTS_PER_POLYGON];
        Point            poly_points_src[MAX_POINTS_PER_POLYGON];
        object_struct    **objects, **objects_src_sphere, *objects_dest_sphere;
        polygons_struct  *polygons, *poly_src_sphere, *poly_dest_sphere;
        double             *input_values, *output_values, *output_values2, dist;
        double             weights[MAX_POINTS_PER_POLYGON];
        double           sphereRadius, r;

        initialize_argument_processing(argc, argv);

        if (!get_string_argument(NULL, &surface_file) ||
            !get_string_argument(NULL, &sphere_file) ||
            !get_string_argument(NULL, &input_values_file) ||
            !get_string_argument(NULL, &output_values_file)) {
                usage(argv[0]);
                exit(EXIT_FAILURE);
        }

        get_int_argument(81920, &n_triangles);

        /*
         * Check tetrahedral topology. Best areal distribution of triangles
         * is achieved for 20 edges
         */
        n_polygons = n_triangles;
        while (n_polygons != 20 && n_polygons > 8 && n_polygons % 4 == 0)
                n_polygons /= 4;

        if (n_polygons != 20) {
                fprintf(stderr, "Warning: Number of triangles %d", n_triangles);
                fprintf(stderr," is not recommended because\ntetrahedral ");
                fprintf(stderr,"topology is not optimal.\n");
                fprintf(stderr,"Please try 20*(4*x) triangles (e.g. 81920).\n");
        }
	
        if (input_graphics_any_format(surface_file, &format,
                                      &n_objects, &objects) != OK ||
            n_objects != 1 || get_object_type(objects[0]) != POLYGONS) {
                fprintf(stderr, "File %s must contain 1 polygons object.\n",
                        surface_file);
                exit(EXIT_FAILURE);
        }

        if (input_graphics_any_format(sphere_file, &format,
                                      &n_objects_src_sphere,
                                      &objects_src_sphere) != OK ||
            n_objects_src_sphere != 1 ||
            get_object_type(objects_src_sphere[0]) != POLYGONS ) {
                fprintf(stderr, "File %s must contain 1 polygons object.\n",
                        sphere_file);
                exit(EXIT_FAILURE);
        }

        polygons = get_polygons_ptr(objects[0]);
        poly_src_sphere = get_polygons_ptr(objects_src_sphere[0]);

        ALLOC(input_values, polygons->n_points);
        ALLOC(output_values, poly_src_sphere->n_points);
        ALLOC(output_values2, poly_src_sphere->n_points);

        if (input_values_any_format(input_values_file, &n_values, &input_values) != OK) {
                fprintf(stderr, "Cannot read values in %s.\n", input_values_file);
    	        exit(EXIT_FAILURE);
        }

        /*
         * Determine radius for the output sphere.  The sphere is not always
         * perfectly spherical, thus use average radius
         */
        sphereRadius = 0.0;
        for (i = 0; i < poly_src_sphere->n_points; i++) {
                r = 0.0;
                for (k = 0; k < 3; k++) 
                        r += Point_coord(poly_src_sphere->points[i], k) *
                             Point_coord(poly_src_sphere->points[i], k);
                sphereRadius += sqrt(r);
        }
        sphereRadius /= poly_src_sphere->n_points;
    
        /*
         * Make radius slightly smaller to get sure that the
         * inner side of handles will be found as nearest point on the surface
         */
        for (i = 0; i < polygons->n_points; i++) 
                set_vector_length(&polygons->points[i],sphereRadius);

        sphereRadius *= 0.975;
        objects_dest_sphere = create_object(POLYGONS);
        poly_dest_sphere = get_polygons_ptr(objects_dest_sphere);

        create_tetrahedral_sphere(&centre, sphereRadius, sphereRadius,
                                  sphereRadius, n_triangles, poly_dest_sphere);

        create_polygons_bintree(poly_src_sphere,
                                ROUND((Real) poly_src_sphere->n_items * 0.5));

        create_polygons_bintree(polygons,
                                ROUND((Real) polygons->n_items * 0.5));

        for (i = 0; i < polygons->n_points; i++) {
                poly = find_closest_polygon_point(&poly_dest_sphere->points[i],
                                                  poly_src_sphere,
                                                  &point_on_src_sphere);
		
                n_points = get_polygon_points(poly_src_sphere, poly,
                                              poly_points_src);
                get_polygon_interpolation_weights(&point_on_src_sphere,
                                                  n_points, poly_points_src,
                                                  weights);

                if (get_polygon_points(polygons, poly, poly_points) != n_points)
                        handle_internal_error("map_point_between_polygons");

                output_values[i] = 0.0;

                for (k = 0; k < n_points; k++) {
                        output_values[i] += weights[k] * input_values[polygons->indices[POINT_INDEX(polygons->end_indices,poly,k)]];
                }

        }

        for (i = 0; i < polygons->n_points; i++) {
                poly = find_closest_polygon_point(&poly_dest_sphere->points[i],
                                                  polygons,
                                                  &point_on_src_sphere);
		
                n_points = get_polygon_points(polygons, poly,
                                              poly_points_src);
                get_polygon_interpolation_weights(&point_on_src_sphere,
                                                  n_points, poly_points_src,
                                                  weights);

                if (get_polygon_points(polygons, poly, poly_points) != n_points)
                        handle_internal_error("map_point_between_polygons");

                output_values2[i] = 0.0;

                for (k = 0; k < n_points; k++) {
                        output_values2[i] += weights[k] * output_values[polygons->indices[POINT_INDEX(polygons->end_indices,poly,k)]];
                }
        }
		    
        output_values_any_format(output_values_file, polygons->n_points,
                                 output_values2, TYPE_DOUBLE);
        FREE(input_values);
        FREE(output_values);
        FREE(output_values2);
        
        return(EXIT_SUCCESS);
}
