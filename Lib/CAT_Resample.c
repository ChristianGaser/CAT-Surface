/* Rachel Yotter - rachel.yotter@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 *
 */

#include "CAT_Resample.h"

/* correct shifting and scaling of source sphere w.r.t. target sphere */
void
correct_shift_scale_sphere(polygons_struct *source_sphere, polygons_struct *target_sphere,
                polygons_struct **out_source_sphere, polygons_struct **out_target_sphere)
{
        int i;
        
        object_struct **source_sphere_objects, **target_sphere_objects;
        polygons_struct *scaled_source_sphere, *scaled_target_sphere;

        source_sphere_objects = (object_struct **) malloc(sizeof(object_struct *));
        *source_sphere_objects = create_object(POLYGONS);
        scaled_source_sphere = get_polygons_ptr(*source_sphere_objects);
        copy_polygons(source_sphere, scaled_source_sphere);
        for (i = 0; i < scaled_source_sphere->n_points; i++)
                set_vector_length(&scaled_source_sphere->points[i], 100.0);

        target_sphere_objects = (object_struct **) malloc(sizeof(object_struct *));
        *target_sphere_objects = create_object(POLYGONS);
        scaled_target_sphere = get_polygons_ptr(*target_sphere_objects);
        copy_polygons(target_sphere, scaled_target_sphere);
        for (i = 0; i < scaled_target_sphere->n_points; i++)
                set_vector_length(&scaled_target_sphere->points[i], 100.0);
                
        correct_bounds_to_target(scaled_source_sphere, scaled_target_sphere);

        *out_source_sphere = scaled_source_sphere;
        *out_target_sphere = scaled_target_sphere;

//        delete_object_list(1, source_sphere_objects);
//        delete_object_list(1, target_sphere_objects);

}

void
resample_values_sphere_noscale(polygons_struct *source_sphere, polygons_struct *target_sphere,
               double *invals, double *outvals)
{
        int    i, j, t, poly, n_points;
        Point  point;
        Point  poly_points[MAX_POINTS_PER_POLYGON];
        double weights[MAX_POINTS_PER_POLYGON];

        if (source_sphere->bintree == NULL) {
                create_polygons_bintree(source_sphere,
                                        ROUND((double) source_sphere->n_items * 0.5));
        }

        for (i = 0; i < target_sphere->n_points; i++) {
                poly = find_closest_polygon_point(&target_sphere->points[i],
                                                  source_sphere, &point);
    
                n_points = get_polygon_points(source_sphere, poly, poly_points);
                get_polygon_interpolation_weights(&point, n_points, poly_points,
                                                  weights);

                outvals[i] = 0.0;
                for (j = 0; j < n_points; j++) {
                        outvals[i] += weights[j] * invals[source_sphere->indices[
                                                POINT_INDEX(source_sphere->end_indices,
                                                            poly,j)]];
                }
        }
        delete_the_bintree(&source_sphere->bintree);
}

/* resample values from source sphere onto target sphere */
void
resample_values_sphere(polygons_struct *source_sphere, polygons_struct *target_sphere,
                double *invals, double *outvals, int scale_and_shift)
{
        int i;
        polygons_struct *scaled_source_sphere, *scaled_target_sphere;
        
        if (scale_and_shift) { /* correct shifts and vector length */
                correct_shift_scale_sphere(source_sphere, target_sphere, &scaled_source_sphere,
                        &scaled_target_sphere);

                resample_values_sphere_noscale(scaled_source_sphere, scaled_target_sphere, 
                        invals, outvals);
                
        } else  resample_values_sphere_noscale(source_sphere, target_sphere, invals, outvals);
}

/* resample surface and values to space defined by target sphere */
object_struct **
resample_surface_to_target_sphere(polygons_struct *polygons, polygons_struct *polygons_sphere, polygons_struct *target_sphere, 
                                  double *input_values, double *output_values, int label_interpolation)
{
        int             j, k, poly, n_points;
        int             label_values[65536], n_labels;
        int             *n_neighbours, **neighbours, *histo;
        signed long     i, min_val = 2147483647, max_val = -2147483647, longval;
        Point           point, scaled_point, center;
        Point           *new_points, poly_points[MAX_POINTS_PER_POLYGON];
        object_struct   **objects, **scaled_objects;
        polygons_struct *scaled_target_sphere, *scaled_polygons_sphere;
        double          max_prob, val, weights[MAX_POINTS_PER_POLYGON];

        objects  = (object_struct **) malloc(sizeof(object_struct *));
        *objects = create_object(POLYGONS);
        scaled_polygons_sphere = get_polygons_ptr(*objects);
        
        if (label_interpolation) {
                /* get maximum and minimum (!=0) value for label interpolation */
                for (i = 0; i < polygons_sphere->n_points; i++) {
                        longval = (signed long)input_values[i];
                        if ((min_val > longval) && (longval != 0))
                                min_val = longval;
                        if (max_val < longval)
                                max_val = longval;
                }
                /* build histogram of label values which is necessary for large and sparse label entries */
                histo = (int *) malloc(sizeof(int) * (max_val + 1 - min_val));
                memset(histo, 0, sizeof(int) * (max_val + 1 - min_val));
                for (i = 0; i < polygons_sphere->n_points; i++) {
                        longval = (signed long)input_values[i];
                        if (longval == 0) continue;
                        histo[longval - min_val]++;
                }
                /* collect label values */
                n_labels = 0;
                for (i = 0; i < (max_val + 1 - min_val); i++) {
                        if (histo[i] == 0) continue;
                        label_values[n_labels] = i + min_val;
                        n_labels++;
                }
        }
        /* if no source sphere is defined a tetrahedral topology of the surface is assumed
           where the corresponding sphere can be simply estimated by its tetrahedral topology */
        if (polygons_sphere == NULL) {
                fill_Point(center, 0.0, 0.0, 0.0);
                create_tetrahedral_sphere(&center, 100.0, 100.0, 100.0,
                                  polygons->n_items, scaled_polygons_sphere);
        } else {
                copy_polygons(polygons_sphere, scaled_polygons_sphere);
                for (i = 0; i < scaled_polygons_sphere->n_points; i++)
                        set_vector_length(&scaled_polygons_sphere->points[i], 100.0);
        }
        
        /* scale the re-parameterizing sphere... also the returned object. */
        scaled_objects = (object_struct **) malloc(sizeof(object_struct *));
        *scaled_objects = create_object(POLYGONS);
        scaled_target_sphere = get_polygons_ptr(*scaled_objects);
        copy_polygons(target_sphere, scaled_target_sphere);

        for (i = 0; i < scaled_target_sphere->n_points; i++)
                set_vector_length(&scaled_target_sphere->points[i], 100.0);

        /* correct sphere center based on bounds of target */
        correct_bounds_to_target(scaled_polygons_sphere, scaled_target_sphere);    

        create_polygons_bintree(scaled_polygons_sphere,
                                ROUND((double) scaled_polygons_sphere->n_items * 0.5));

        new_points = (Point *) malloc(sizeof(Point) * scaled_target_sphere->n_points);

        for (i = 0; i < scaled_target_sphere->n_points; i++) {
                poly = find_closest_polygon_point(&scaled_target_sphere->points[i],
                                                  scaled_polygons_sphere, &point);
                                                  
                n_points = get_polygon_points(scaled_polygons_sphere, poly, poly_points);
                get_polygon_interpolation_weights(&point, n_points, poly_points, weights);

                if (polygons != NULL) 
                        if (get_polygon_points(polygons, poly, poly_points) != n_points)
                                fprintf(stderr,"map_point_between_polygons\n");

                fill_Point(new_points[i], 0.0, 0.0, 0.0);
                
                if (input_values != NULL) output_values[i] = 0.0;
                

                /* resample mesh if defined */
                if (polygons != NULL) {
                        for (j = 0; j < n_points; j++) {
                                SCALE_POINT(scaled_point, poly_points[j], weights[j]);
                                ADD_POINTS(new_points[i], new_points[i], scaled_point);
                        }
                }
 
                /* apply interpolation for each ROI label separately */
                if (label_interpolation) {

                        if (input_values != NULL) {
                                max_prob = 0.0;
                                for (k = 0; k < n_labels; k++) {
                                        val = 0.0;
                                        for (j = 0; j < n_points; j++) {
                                                if (label_values[k] == (int)input_values[scaled_polygons_sphere->indices[
                                                        POINT_INDEX(scaled_polygons_sphere->end_indices,poly,j)]])
                                                        val += weights[j];
                                        }
                                        if (max_prob < val) {
                                                max_prob = val;
                                                output_values[i] = (double)label_values[k];
                                        }
                                }
                        }
                } else {
                            
                        if (input_values != NULL) output_values[i] = 0.0;        
                        for (j = 0; j < n_points; j++) {
                                if (input_values != NULL)
                                        output_values[i] += weights[j] * input_values[scaled_polygons_sphere->indices[
                                                        POINT_INDEX(scaled_polygons_sphere->end_indices,poly,j)]];
                        }
                }
        }

        if (polygons != NULL) {
                free(scaled_target_sphere->points);
                scaled_target_sphere->points = new_points;
        }

        if (label_interpolation) free(histo);
                
        compute_polygon_normals(scaled_target_sphere);

        delete_the_bintree(&scaled_polygons_sphere->bintree);
        delete_object_list(1, objects);

        return(scaled_objects);
}


/* resample surface and values to space defined by tetrahedral sphere */
object_struct **
resample_surface(polygons_struct *surface, polygons_struct *sphere,
                 int n_triangles, double *invals, double *outvals)
{
        int             i, j, poly, n_points;
        Point           point, scaled_point, center;
        Point           *new_points, poly_points[MAX_POINTS_PER_POLYGON];
        object_struct   **objects, *scaled_objects;
        polygons_struct *output_surface, *scaled_sphere;
        double          weights[MAX_POINTS_PER_POLYGON];

        scaled_objects = create_object(POLYGONS);
        scaled_sphere  = get_polygons_ptr(scaled_objects);
        copy_polygons(sphere, scaled_sphere);

        translate_to_center_of_mass(scaled_sphere);
        for (i = 0; i < scaled_sphere->n_points; i++)
                set_vector_length(&scaled_sphere->points[i], 100.0);

        objects  = (object_struct **) malloc(sizeof(object_struct *));
        *objects = create_object(POLYGONS);
        output_surface = get_polygons_ptr(*objects);
        fill_Point(center, 0.0, 0.0, 0.0);
        create_tetrahedral_sphere(&center, 100.0, 100.0, 100.0, n_triangles,
                                  output_surface);

        create_polygons_bintree(sphere, ROUND((double) sphere->n_items * 0.5));

        new_points = (Point *) malloc(sizeof(Point) * output_surface->n_points);
        
        if (invals != NULL) 
                outvals = (double *) malloc(sizeof(double) * output_surface->n_points);

        for (i = 0; i < output_surface->n_points; i++) {
                poly = find_closest_polygon_point(&output_surface->points[i],
                                                  sphere, &point);
    
                n_points = get_polygon_points(sphere, poly, poly_points);
                get_polygon_interpolation_weights(&point, n_points, poly_points,
                                                  weights);

                if (get_polygon_points(surface, poly, poly_points) != n_points)
                        fprintf(stderr,"map_point_between_polygons\n");

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

        free(output_surface->points);
        output_surface->points = new_points;

        compute_polygon_normals(output_surface);

        delete_the_bintree(&sphere->bintree);

        return(objects);
}
