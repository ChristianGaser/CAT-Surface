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
Usage: %s  surface_file sphere_file output_surface_file n_triangles [input_values_file output_values_file]\n\
Resamples a spherical inflated surface to a sphere.\n\
\n\n";

    fprintf(stderr, usage_str, executable);
}

int
main(int argc, char *argv[])
{
    char       *surface_file, *sphere_file, *output_surface_file;
    char       *input_values_file, *output_values_file;
    File_formats   format;
    int        n_objects, point;
    int        n_triangles, n_polygons;
    int        i, j, k;
    int        poly, n_points, n_intersections, n_values;
    int        *n_neighbours, **neighbours;
    Point      center, point_on_src_sphere, scaled_point;
    Point      poly_points[MAX_POINTS_PER_POLYGON];
    Point      poly_points_src[MAX_POINTS_PER_POLYGON];
    Point      *new_points;
    object_struct  **objects, **objects_src_sphere, *objects_target_sphere;
    polygons_struct  *polygons, *polygons_sphere, *target_sphere;
    double       *input_values, *output_values, dist;
    BOOLEAN      values_specified;
    double       weights[MAX_POINTS_PER_POLYGON];
    double       sphereRadius, r, bounds[6];

    initialize_argument_processing(argc, argv);

    if (!get_string_argument(NULL, &surface_file) ||
      !get_string_argument(NULL, &sphere_file) ||
      !get_string_argument(NULL, &output_surface_file)) {
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
        printf("Warning: Number of triangles %d", n_triangles);
        printf(" is not recommended because\ntetrahedral ");
        printf("topology is not optimal.\n");
        printf("Please try 20*(4*x) triangles (e.g. 81920).\n");
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

    polygons = get_polygons_ptr(objects[0]);
    polygons_sphere = get_polygons_ptr(objects_src_sphere[0]);

    values_specified = (get_string_argument(NULL, &input_values_file) &&
              get_string_argument(NULL, &output_values_file));

    if (values_specified) {
        ALLOC(input_values, polygons->n_points);

        if (input_values_any_format(input_values_file, &n_values, &input_values) != OK) {
            fprintf(stderr, "Cannot read values in %s.\n", input_values_file);
           exit(EXIT_FAILURE);
        }

    }

    /*
     * Determine radius for the output sphere.  The sphere is not always
     * perfectly spherical, thus use average radius
     */
    sphereRadius = 0.0;
    for (i = 0; i < polygons_sphere->n_points; i++) {
        r = 0.0;
        for (k = 0; k < 3; k++) 
            r += Point_coord(polygons_sphere->points[i], k) *
               Point_coord(polygons_sphere->points[i], k);
        sphereRadius += sqrt(r);
    }
    sphereRadius /= polygons_sphere->n_points;

    /* Calc. sphere center based on bounds of input (correct for shifts) */
    get_bounds(polygons_sphere, bounds);
    fill_Point(center, bounds[0]+bounds[1],
               bounds[2]+bounds[3], bounds[4]+bounds[5]);
  
    objects_target_sphere = create_object(POLYGONS);
    target_sphere = get_polygons_ptr(objects_target_sphere);
  
    /*
     * Make radius slightly smaller to get sure that the
     * inner side of handles will be found as nearest point on the surface
     */
    sphereRadius *= 0.975;
    create_tetrahedral_sphere(&center, sphereRadius, sphereRadius,
                  sphereRadius, n_triangles, target_sphere);

    if (values_specified)
        ALLOC(output_values, target_sphere->n_points);

    create_polygons_bintree(polygons_sphere,
                ROUND((Real) polygons_sphere->n_items * 0.5));

    ALLOC(new_points, target_sphere->n_points);

    for (i = 0; i < target_sphere->n_points; i++) {
        poly = find_closest_polygon_point(&target_sphere->points[i],
                          polygons_sphere,
                          &point_on_src_sphere); 
        n_points = get_polygon_points(polygons_sphere, poly,
                        poly_points_src);
        get_polygon_interpolation_weights(&point_on_src_sphere,
                          n_points, poly_points_src,
                          weights);

        if (get_polygon_points(polygons, poly, poly_points) != n_points)
            fprintf(stderr,"map_point_between_polygons\n");

        fill_Point(new_points[i], 0.0, 0.0, 0.0);
        if (values_specified)
            output_values[i] = 0.0;

        for (k = 0; k < n_points; k++) {
            SCALE_POINT(scaled_point, poly_points[k], (double)weights[k]);
            ADD_POINTS(new_points[i], new_points[i], scaled_point);
            if (values_specified)
                output_values[i] += weights[k] * 
                    input_values[polygons->indices[POINT_INDEX(polygons->end_indices,poly,k)]];
        }
     }

    create_polygon_point_neighbours(target_sphere, TRUE, &n_neighbours,
                    &neighbours, NULL, NULL);

    for (i = 0; i < target_sphere->n_points; i++) {
        target_sphere->points[i] = new_points[i];
    }

    compute_polygon_normals(target_sphere);

    if(output_graphics_any_format(output_surface_file, format, 1,
                   &objects_target_sphere, NULL) != OK)
          exit(EXIT_FAILURE);
  
    if (values_specified) {
        output_values_any_format(output_values_file,
                     target_sphere->n_points,
                     output_values, TYPE_DOUBLE);
        FREE(input_values);
        FREE(output_values);
    }
    
    delete_the_bintree(&polygons_sphere->bintree);
    FREE(new_points);
    return(EXIT_SUCCESS);
}
