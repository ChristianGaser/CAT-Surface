/* Rachel Yotter - rachel.yotter@uni-jena.de
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
#include "CAT_Complexity.h"

object_struct **
resample_tetrahedron(polygons_struct *surface, polygons_struct *sphere,
                     int n_triangles, double *invals, double *outvals)
{
        int i, j, t, poly, n_points;
        Point point, scaled_point, centre;
        Point *new_points, poly_points[MAX_POINTS_PER_POLYGON];
        object_struct **objects, *scaled_objects;
        polygons_struct *platonic_solid, *scaled_sphere;
        double radius, r;
        double weights[MAX_POINTS_PER_POLYGON];

        /* Check tetrahedral topology. Best areal distribution of triangles
         * is achieved for 20 edges
         */ /*
        t = n_triangles;
        while (t != 20 && t > 8 && t % 4 == 0)
                t /= 4;

        if (t != 20) {
                fprintf(stderr, "Warning: Number of triangles %d", n_triangles);
                fprintf(stderr," is not recommended because\ntetrahedral ");
                fprintf(stderr,"topology is not optimal.\n");
                fprintf(stderr,"Please try 20*(4*x) triangles (e.g. 81920).\n");
        } */

        scaled_objects = create_object(POLYGONS);
        scaled_sphere = get_polygons_ptr(scaled_objects);
        copy_polygons(sphere, scaled_sphere);

        translate_to_center_of_mass(scaled_sphere);
        for (i = 0; i < scaled_sphere->n_points; i++)
                set_vector_length(&scaled_sphere->points[i], 100.0);

        objects = (object_struct **) malloc(sizeof(object_struct *));
        *objects = create_object(POLYGONS);
        platonic_solid = get_polygons_ptr(*objects);
        fill_Point(centre, 0.0, 0.0, 0.0);
        create_tetrahedral_sphere(&centre, 100.0, 100.0, 100.0, n_triangles,
                                  platonic_solid);

        create_polygons_bintree(sphere, ROUND((Real) sphere->n_items * 0.5));

        new_points = (Point *) malloc(sizeof(Point) * platonic_solid->n_points);
        if (invals != NULL) {
                outvals = (double *) malloc(sizeof(double) *
                                            platonic_solid->n_points);
        }

        for (i = 0; i < platonic_solid->n_points; i++) {
                poly = find_closest_polygon_point(&platonic_solid->points[i],
                                                  sphere, &point);
		
                n_points = get_polygon_points(sphere, poly, poly_points);
                get_polygon_interpolation_weights(&point, n_points, poly_points,
                                                  weights);

                if (get_polygon_points(surface, poly, poly_points) != n_points)
                        handle_internal_error("map_point_between_polygons");

                fill_Point(new_points[i], 0.0, 0.0, 0.0);
                if (invals != NULL)
                        outvals[i] = 0.0;

                for (j = 0; j < n_points; j++) {
                        SCALE_POINT(scaled_point, poly_points[j], weights[j]);
                        ADD_POINTS(new_points[i], new_points[i], scaled_point);
                        if (invals != NULL) {
                                outvals[i] += weights[j] *
                                               invals[surface->indices[
                                     POINT_INDEX(surface->end_indices,poly,j)]];
                        }
                }
        }

        free(platonic_solid->points);
        platonic_solid->points = new_points;

        compute_polygon_normals(platonic_solid);

        return(objects);
}

void
resample_noscale(polygons_struct *source, polygons_struct *target,
               double *invals, double *outvals)
{
        int i, j, t, poly, n_points;
        Point point, centre;
        Point *new_points, poly_points[MAX_POINTS_PER_POLYGON];
        double radius, r;
        double weights[MAX_POINTS_PER_POLYGON];

        if (source->bintree == NULL) {
                create_polygons_bintree(source,
                                        ROUND((Real) source->n_items * 0.5));
        }

        for (i = 0; i < target->n_points; i++) {
                poly = find_closest_polygon_point(&target->points[i],
                                                  source, &point);
		
                n_points = get_polygon_points(source, poly, poly_points);
                get_polygon_interpolation_weights(&point, n_points, poly_points,
                                                  weights);

                outvals[i] = 0.0;
                for (j = 0; j < n_points; j++) {
                        outvals[i] += weights[j] * invals[source->indices[
                                                POINT_INDEX(source->end_indices,
                                                            poly,j)]];
                }
        }
        delete_the_bintree(&source->bintree);
}

/* resample values from source sphere onto target sphere */
void
resample_values(polygons_struct *source, polygons_struct *target,
                double *invals, double *outvals)
{
        int i;
        object_struct **source_objects, **target_objects;
        polygons_struct *scaled_source, *scaled_target;

        source_objects = (object_struct **) malloc(sizeof(object_struct *));
        *source_objects = create_object(POLYGONS);
        scaled_source = get_polygons_ptr(*source_objects);
        copy_polygons(source, scaled_source);
        translate_to_center_of_mass(scaled_source);
        for (i = 0; i < scaled_source->n_points; i++)
                set_vector_length(&scaled_source->points[i], 100.0);

        target_objects = (object_struct **) malloc(sizeof(object_struct *));
        *target_objects = create_object(POLYGONS);
        scaled_target = get_polygons_ptr(*target_objects);
        copy_polygons(target, scaled_target);
        translate_to_center_of_mass(scaled_target);
        for (i = 0; i < scaled_target->n_points; i++)
                set_vector_length(&scaled_target->points[i], 100.0);

        resample_noscale(scaled_source, scaled_target, invals, outvals);

        delete_object_list(1, source_objects);
        delete_object_list(1, target_objects);
}

/* resample back into the original object space... polygons has tetrahedral
 * topology as its sphere. */
object_struct **
resample_surface_sphere(polygons_struct *polygons, polygons_struct *sphere)
{
        int i, j, t, poly, n_points;
        Point point, scaled_point, centre;
        Point *new_points, poly_points[MAX_POINTS_PER_POLYGON];
        object_struct **objects, **scaled_objects;
        polygons_struct *polygons_sphere, *scaled_sphere;
        double weights[MAX_POINTS_PER_POLYGON];

        objects = (object_struct **) malloc(sizeof(object_struct *));
        *objects = create_object(POLYGONS);
        polygons_sphere = get_polygons_ptr(*objects);
        fill_Point(centre, 0.0, 0.0, 0.0);
        create_tetrahedral_sphere(&centre, 100.0, 100.0, 100.0,
                                  polygons->n_items, polygons_sphere);

        /* scale the re-parameterizing sphere... also the return object. */
        scaled_objects = (object_struct **) malloc(sizeof(object_struct *));
        *scaled_objects = create_object(POLYGONS);
        scaled_sphere = get_polygons_ptr(*scaled_objects);
        copy_polygons(sphere, scaled_sphere);

        translate_to_center_of_mass(scaled_sphere);
        for (i = 0; i < scaled_sphere->n_points; i++)
                set_vector_length(&scaled_sphere->points[i], 100.0);

        create_polygons_bintree(scaled_sphere,
                                ROUND((Real) scaled_sphere->n_items * 0.5));

        new_points = (Point *) malloc(sizeof(Point) * scaled_sphere->n_points);

        for (i = 0; i < scaled_sphere->n_points; i++) {
                poly = find_closest_polygon_point(&scaled_sphere->points[i],
                                                  polygons_sphere, &point);
		
                n_points = get_polygon_points(polygons_sphere, poly,
                                              poly_points);
                get_polygon_interpolation_weights(&point, n_points, poly_points,
                                                  weights);

                get_polygon_points(polygons, poly, poly_points);

                fill_Point(new_points[i], 0.0, 0.0, 0.0);

                for (j = 0; j < n_points; j++) {
                        SCALE_POINT(scaled_point, poly_points[j], weights[j]);
                        ADD_POINTS(new_points[i], new_points[i], scaled_point);
                }
        }

        free(scaled_sphere->points);
        scaled_sphere->points = new_points;

        compute_polygon_normals(scaled_sphere);

        delete_object_list(1, objects);

        return(scaled_objects);
}


object_struct **
resample_surface(polygons_struct *surface, polygons_struct *sphere,
                 int n_triangles, double *invals, double *outvals)
{
        int i, j, t, poly, n_points;
        Point point, scaled_point, centre;
        Point *new_points, poly_points[MAX_POINTS_PER_POLYGON];
        object_struct **objects, *scaled_objects;
        polygons_struct *platonic_solid, *scaled_sphere;
        double radius, r;
        double weights[MAX_POINTS_PER_POLYGON];

        /* Check tetrahedral topology. Best areal distribution of triangles
         * is achieved for 20 edges
         */ /*
        t = n_triangles;
        while (t != 20 && t > 8 && t % 4 == 0)
                t /= 4;

        if (t != 20) {
                fprintf(stderr, "Warning: Number of triangles %d", n_triangles);
                fprintf(stderr," is not recommended because\ntetrahedral ");
                fprintf(stderr,"topology is not optimal.\n");
                fprintf(stderr,"Please try 20*(4*x) triangles (e.g. 81920).\n");
        } */

        scaled_objects = create_object(POLYGONS);
        scaled_sphere = get_polygons_ptr(scaled_objects);
        copy_polygons(sphere, scaled_sphere);

        translate_to_center_of_mass(scaled_sphere);
        for (i = 0; i < scaled_sphere->n_points; i++)
                set_vector_length(&scaled_sphere->points[i], 100.0);

        objects = (object_struct **) malloc(sizeof(object_struct *));
        *objects = create_object(POLYGONS);
        platonic_solid = get_polygons_ptr(*objects);
        fill_Point(centre, 0.0, 0.0, 0.0);
        create_tetrahedral_sphere(&centre, 100.0, 100.0, 100.0, n_triangles,
                                  platonic_solid);

        create_polygons_bintree(sphere, ROUND((Real) sphere->n_items * 0.5));

        new_points = (Point *) malloc(sizeof(Point) * platonic_solid->n_points);
        if (invals != NULL) {
                outvals = (double *) malloc(sizeof(double) *
                                            platonic_solid->n_points);
        }

        for (i = 0; i < platonic_solid->n_points; i++) {
                poly = find_closest_polygon_point(&platonic_solid->points[i],
                                                  sphere, &point);
		
                n_points = get_polygon_points(sphere, poly, poly_points);
                get_polygon_interpolation_weights(&point, n_points, poly_points,
                                                  weights);

                if (get_polygon_points(surface, poly, poly_points) != n_points)
                        handle_internal_error("map_point_between_polygons");

                fill_Point(new_points[i], 0.0, 0.0, 0.0);
                if (invals != NULL)
                        outvals[i] = 0.0;

                for (j = 0; j < n_points; j++) {
                        SCALE_POINT(scaled_point, poly_points[j], weights[j]);
                        ADD_POINTS(new_points[i], new_points[i], scaled_point);
                        if (invals != NULL) {
                                outvals[i] += weights[j] *
                                               invals[surface->indices[
                                     POINT_INDEX(surface->end_indices,poly,j)]];
                        }
                }
        }

        free(platonic_solid->points);
        platonic_solid->points = new_points;

        compute_polygon_normals(platonic_solid);

        return(objects);
}
