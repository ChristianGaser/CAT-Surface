/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 *
 */

#include <volume_io/internal_volume_io.h>
#include <bicpl.h>

#include "CAT_Surf.h"
#include "CAT_SurfaceIO.h"

void
usage(char *executable)
{
        char *usage_str = "\n\
Usage: %s  src.obj src_sphere.obj target_sphere.obj resampled_output.obj [input_values.txt output_values.txt]\n\
Resamples a spherical inflated surface to an external defined sphere.\n\
\n\n";

        fprintf(stderr, usage_str, executable);
}

int
main(int argc, char *argv[])
{
        char             *surface_file, *sphere_file, *output_file, *target_sphere_file;
        char             *input_values_file, *output_values_file;
        File_formats     format;
        int              n_objects, point;
        int              i, j, k;
        int              poly, n_points, n_intersections, n_values;
        int              *n_neighbours, **neighbours;
        Point            center, point_on_src_sphere, scaled_point;
        Point            poly_points[MAX_POINTS_PER_POLYGON];
        Point            poly_points_src[MAX_POINTS_PER_POLYGON];
        Point            *new_points;
        object_struct    **objects, **objects_src_sphere, **objects_target_sphere;
        polygons_struct  *polygons, *poly_src_sphere, *poly_target_sphere;
        double             *input_values, *output_values, dist;
        BOOLEAN          values_specified;
        double             weights[MAX_POINTS_PER_POLYGON];
        double           bounds_dest[6], bounds_src[6];

        initialize_argument_processing(argc, argv);

        if (!get_string_argument(NULL, &surface_file    ) ||
            !get_string_argument(NULL, &sphere_file     ) ||
            !get_string_argument(NULL, &target_sphere_file     ) ||
            !get_string_argument(NULL, &output_file)) {
                usage(argv[0]);
                exit(EXIT_FAILURE);
        }

        if (input_graphics_any_format(surface_file, &format,
                                      &n_objects, &objects) != OK ||
            n_objects != 1 || get_object_type(objects[0]) != POLYGONS) {
                fprintf(stderr, "File %s must contain 1 polygons object.\n",
                        surface_file);
                exit(EXIT_FAILURE);
        }

        if (input_graphics_any_format(sphere_file, &format,
                                      &n_objects,
                                      &objects_src_sphere) != OK ||
            n_objects != 1 ||
            get_object_type(objects_src_sphere[0]) != POLYGONS ) {
                fprintf(stderr, "File %s must contain 1 polygons object.\n",
                        sphere_file);
                exit(EXIT_FAILURE);
        }

        if (input_graphics_any_format(target_sphere_file, &format,
                                      &n_objects,
                                      &objects_target_sphere) != OK ||
            n_objects != 1 ||
            get_object_type(objects_target_sphere[0]) != POLYGONS ) {
                fprintf(stderr, "File %s must contain 1 polygons object.\n",
                        target_sphere_file);
                exit(EXIT_FAILURE);
        }

        polygons = get_polygons_ptr(objects[0]);
        poly_src_sphere = get_polygons_ptr(objects_src_sphere[0]);
        poly_target_sphere = get_polygons_ptr(objects_target_sphere[0]);

        values_specified = get_string_argument(NULL, &input_values_file) &&
                           get_string_argument(NULL, &output_values_file);

        if (values_specified) {
                ALLOC(input_values, polygons->n_points);
                ALLOC(output_values, poly_target_sphere->n_points);

                if (input_values_any_format(input_values_file, &n_values, &input_values) != OK) {
                        fprintf(stderr, "Cannot read values in %s.\n", input_values_file);
    	               exit(EXIT_FAILURE);
                }

        }
        
        /*
         * set spheres to same radius
         */

        for (i = 0; i < poly_src_sphere->n_points; i++) 
                set_vector_length(&poly_src_sphere->points[i], 100.0);

        for (i = 0; i < poly_target_sphere->n_points; i++) 
                set_vector_length(&poly_target_sphere->points[i], 100.0);

        /* Calc. sphere center based on bounds of input (correct for shifts) */    
        get_bounds(poly_target_sphere, bounds_dest);
        get_bounds(poly_src_sphere, bounds_src);
        fill_Point(center, bounds_src[0]-bounds_dest[0],
                           bounds_src[2]-bounds_dest[2], bounds_src[4]-bounds_dest[4]);

        for (i = 0; i < poly_src_sphere->n_points; i++) {
                for (j = 0; j < 3; j++)
                        Point_coord(poly_src_sphere->points[i], j) -= Point_coord(center, j);
        }

        create_polygons_bintree(poly_src_sphere,
                                ROUND((Real) poly_src_sphere->n_items * 0.5));

        ALLOC(new_points, poly_target_sphere->n_points);

        for (i = 0; i < poly_target_sphere->n_points; i++) {
                poly = find_closest_polygon_point(&poly_target_sphere->points[i],
                                                  poly_src_sphere,
                                                  &point_on_src_sphere);
		
                n_points = get_polygon_points(poly_src_sphere, poly,
                                              poly_points_src);
                get_polygon_interpolation_weights(&point_on_src_sphere,
                                                  n_points, poly_points_src,
                                                  weights);

                if (get_polygon_points(polygons, poly, poly_points) != n_points)
                        handle_internal_error("map_point_between_polygons");

                fill_Point(new_points[i], 0.0, 0.0, 0.0);
                if (values_specified)
                        output_values[i] = 0.0;

                for (k = 0; k < n_points; k++) {
                        SCALE_POINT(scaled_point, poly_points[k], weights[k]);
                        ADD_POINTS(new_points[i], new_points[i], scaled_point);
                        if (values_specified)
                                output_values[i] += weights[k] *
input_values[polygons->indices[POINT_INDEX(polygons->end_indices,poly,k)]];
                }
       }

        create_polygon_point_neighbours(poly_target_sphere, TRUE, &n_neighbours,
                                        &neighbours, NULL, NULL);

        for (i = 0; i < poly_target_sphere->n_points; i++) {
                poly_target_sphere->points[i] = new_points[i];
        }
		
        compute_polygon_normals(poly_target_sphere);

        if(output_graphics_any_format(output_file, format, 1,
                                   objects_target_sphere) != OK)
                    exit(EXIT_FAILURE);
    
        if (values_specified) {
                output_values_any_format(output_values_file,
                                         poly_target_sphere->n_points,
                                         output_values, TYPE_DOUBLE);
                FREE(input_values);
                FREE(output_values);
        }
        
        delete_the_bintree(&poly_src_sphere->bintree);
        FREE(new_points);
        return(EXIT_SUCCESS);
}
