/* Rachel Yotter - rachel.yotter@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 *
 */

#include "CAT_Resample.h"

/* Multithreading stuff*/
#if defined(_WIN32) || defined(_WIN64)
    #define THREAD_HANDLE HANDLE
    #define THREAD_RETURN DWORD WINAPI
#else
    #include <pthread.h>
    #define THREAD_HANDLE pthread_t
    #define THREAD_RETURN void *
#endif

// Define the struct to pass arguments to threads
typedef struct {
    int start_idx;
    int end_idx;
    polygons_struct *source_sphere;
    polygons_struct *target_sphere;
    double *invals;
    double *outvals;
} ThreadArgs_process_target_points;

// Define the struct to pass arguments to threads
typedef struct {
    int start_idx;
    int end_idx;
    polygons_struct *polygons;
    polygons_struct *poly_src_sphere;
    polygons_struct *resampled_source;
    double *input_values;
    double *output_values;
    Point *new_points;
} ThreadArgs;

// Thread function for parallel execution
#if !defined(_WIN32) && !defined(_WIN64)
THREAD_RETURN process_resample_surface(void *args) {
    ThreadArgs *thread_args = (ThreadArgs *)args;
    int i, k, poly, n_points;
    Point point_on_src_sphere, scaled_point;
    Point poly_points[MAX_POINTS_PER_POLYGON];
    Point poly_points_src[MAX_POINTS_PER_POLYGON];
    float weights[MAX_POINTS_PER_POLYGON];

    for (i = thread_args->start_idx; i < thread_args->end_idx; i++) {
        poly = find_closest_polygon_point(&thread_args->resampled_source->points[i],
                                          thread_args->poly_src_sphere,
                                          &point_on_src_sphere);

        n_points = get_polygon_points(thread_args->poly_src_sphere, poly, poly_points_src);
        get_polygon_interpolation_weights(&point_on_src_sphere, n_points, poly_points_src, weights);

        if (get_polygon_points(thread_args->polygons, poly, poly_points) != n_points) {
            fprintf(stderr, "map_point_between_polygons\n");
        }

        fill_Point(thread_args->new_points[i], 0.0, 0.0, 0.0);
        if (thread_args->input_values != NULL) {
            thread_args->output_values[i] = 0.0;
        }

        for (k = 0; k < n_points; k++) {
            SCALE_POINT(scaled_point, poly_points[k], (double)weights[k]);
            ADD_POINTS(thread_args->new_points[i], thread_args->new_points[i], scaled_point);
            if (thread_args->input_values != NULL) {
                thread_args->output_values[i] += weights[k] *
                    thread_args->input_values[thread_args->polygons->indices[
                    POINT_INDEX(thread_args->polygons->end_indices, poly, k)]];
            }
        }
    }
    return NULL;
}
#endif

void resample_spherical_surface(polygons_struct *polygons,
                                polygons_struct *poly_src_sphere,
                                polygons_struct *resampled_source,
                                double *input_values, double *output_values,
                                int n_triangles) {
    int i, t, num_threads = 8; // Adjust based on hardware
#if !defined(_WIN32) && !defined(_WIN64)
    THREAD_HANDLE threads[num_threads];
    ThreadArgs thread_args[num_threads];
#endif

    double sphereRadius = 0.0, r, bounds[6];

    // Compute average sphere radius
    for (i = 0; i < poly_src_sphere->n_points; i++) {
        r = 0.0;
        for (int k = 0; k < 3; k++)
            r += Point_coord(poly_src_sphere->points[i], k) * Point_coord(poly_src_sphere->points[i], k);
        sphereRadius += sqrt(r);
    }
    sphereRadius /= poly_src_sphere->n_points;

    // Calculate sphere center
    get_bounds(poly_src_sphere, bounds);
    Point centre;
    fill_Point(centre, bounds[0] + bounds[1], bounds[2] + bounds[3], bounds[4] + bounds[5]);

    // Slightly reduce radius
    sphereRadius *= 0.975;
    create_tetrahedral_sphere(&centre, sphereRadius, sphereRadius, sphereRadius, n_triangles, resampled_source);

    create_polygons_bintree(poly_src_sphere, ROUND((float) poly_src_sphere->n_items * 0.5));

    Point *new_points;
    ALLOC(new_points, resampled_source->n_points);
    if (input_values != NULL)
        ALLOC(output_values, resampled_source->n_points);

    int total_points = resampled_source->n_points;
    int chunk_size = total_points / num_threads;
    int remainder = total_points % num_threads;

#if defined(_WIN32) || defined(_WIN64)
    // **Sequential Execution for Windows**
    for (i = 0; i < resampled_source->n_points; i++) {
        int poly, n_points;
        Point point_on_src_sphere, scaled_point;
        Point poly_points[MAX_POINTS_PER_POLYGON];
        Point poly_points_src[MAX_POINTS_PER_POLYGON];
        float weights[MAX_POINTS_PER_POLYGON];

        poly = find_closest_polygon_point(&resampled_source->points[i], poly_src_sphere, &point_on_src_sphere);

        n_points = get_polygon_points(poly_src_sphere, poly, poly_points_src);
        get_polygon_interpolation_weights(&point_on_src_sphere, n_points, poly_points_src, weights);

        if (get_polygon_points(polygons, poly, poly_points) != n_points) {
            fprintf(stderr, "map_point_between_polygons\n");
        }

        fill_Point(new_points[i], 0.0, 0.0, 0.0);
        if (input_values != NULL)
            output_values[i] = 0.0;

        for (int k = 0; k < n_points; k++) {
            SCALE_POINT(scaled_point, poly_points[k], (double)weights[k]);
            ADD_POINTS(new_points[i], new_points[i], scaled_point);
            if (input_values != NULL)
                output_values[i] += weights[k] * input_values[polygons->indices[
                                      POINT_INDEX(polygons->end_indices, poly, k)]];
        }
    }
#else
    // **Parallel Execution for Non-Windows Systems**
    for (t = 0; t < num_threads; t++) {
        thread_args[t].start_idx = t * chunk_size;
        thread_args[t].end_idx = (t == num_threads - 1) ? (t + 1) * chunk_size + remainder : (t + 1) * chunk_size;
        thread_args[t].polygons = polygons;
        thread_args[t].poly_src_sphere = poly_src_sphere;
        thread_args[t].resampled_source = resampled_source;
        thread_args[t].input_values = input_values;
        thread_args[t].output_values = output_values;
        thread_args[t].new_points = new_points;

        if (pthread_create(&threads[t], NULL, process_resample_surface, &thread_args[t]) != 0) {
            perror("pthread_create");
            exit(EXIT_FAILURE);
        }
    }

    // Join all threads
    for (t = 0; t < num_threads; t++) {
        if (pthread_join(threads[t], NULL) != 0) {
            perror("pthread_join");
            exit(EXIT_FAILURE);
        }
    }
#endif

    // Copy new points to the resampled surface
    for (i = 0; i < resampled_source->n_points; i++) {
        resampled_source->points[i] = new_points[i];
    }

    compute_polygon_normals(resampled_source);
    free(new_points);
    delete_the_bintree(&poly_src_sphere->bintree);
}

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

// Thread function (used only on non-Windows systems)
#if !defined(_WIN32) && !defined(_WIN64)
THREAD_RETURN process_target_points(void *args) {
    ThreadArgs_process_target_points *thread_args = (ThreadArgs_process_target_points *)args;

    int start = thread_args->start_idx;
    int end = thread_args->end_idx;
    polygons_struct *source_sphere = thread_args->source_sphere;
    polygons_struct *target_sphere = thread_args->target_sphere;
    double *invals = thread_args->invals;
    double *outvals = thread_args->outvals;

    int i, j, poly, n_points;
    Point point;
    Point poly_points[MAX_POINTS_PER_POLYGON];
    float weights[MAX_POINTS_PER_POLYGON];

    for (i = start; i < end; i++) {
        poly = find_closest_polygon_point(&target_sphere->points[i], source_sphere, &point);

        n_points = get_polygon_points(source_sphere, poly, poly_points);
        get_polygon_interpolation_weights(&point, n_points, poly_points, weights);

        outvals[i] = 0.0;
        for (j = 0; j < n_points; j++) {
            outvals[i] += (double)weights[j] * invals[source_sphere->indices[
                                   POINT_INDEX(source_sphere->end_indices, poly, j)]];
        }
    }

    return NULL;
}
#endif

// Parallelized or sequential version of the function
void resample_values_sphere_noscale(polygons_struct *source_sphere, 
                                    polygons_struct *target_sphere, 
                                    double *invals, double *outvals) {
    int i, j, t, num_threads = 8; // Number of threads (adjust based on system hardware)
#if !defined(_WIN32) && !defined(_WIN64)
    THREAD_HANDLE threads[num_threads];
    ThreadArgs_process_target_points thread_args[num_threads];
#endif

    // Create bintree if it doesn't already exist
    if (source_sphere->bintree == NULL) {
        create_polygons_bintree(source_sphere,
                                ROUND((double)source_sphere->n_items * 0.5));
    }

    // Divide the workload among threads
    int total_points = target_sphere->n_points;
    int chunk_size = total_points / num_threads;
    int remainder = total_points % num_threads;

#if defined(_WIN32) || defined(_WIN64)
    // **Sequential Execution for Windows**
    for (i = 0; i < total_points; i++) {
        int poly, n_points;
        Point point;
        Point poly_points[MAX_POINTS_PER_POLYGON];
        float weights[MAX_POINTS_PER_POLYGON];

        poly = find_closest_polygon_point(&target_sphere->points[i], source_sphere, &point);

        n_points = get_polygon_points(source_sphere, poly, poly_points);
        get_polygon_interpolation_weights(&point, n_points, poly_points, weights);

        outvals[i] = 0.0;
        for (j = 0; j < n_points; j++) {
            outvals[i] += (double)weights[j] * invals[source_sphere->indices[
                                   POINT_INDEX(source_sphere->end_indices, poly, j)]];
        }

    }
#else
    // **Parallel Execution for Non-Windows Systems**
    for (t = 0; t < num_threads; t++) {
        thread_args[t].start_idx = t * chunk_size;
        thread_args[t].end_idx = (t == num_threads - 1) ? (t + 1) * chunk_size + remainder
                                                        : (t + 1) * chunk_size;
        thread_args[t].source_sphere = source_sphere;
        thread_args[t].target_sphere = target_sphere;
        thread_args[t].invals = invals;
        thread_args[t].outvals = outvals;

        if (pthread_create(&threads[t], NULL, process_target_points, &thread_args[t]) != 0) {
            perror("pthread_create");
            exit(EXIT_FAILURE);
        }
    }

    // Join all threads
    for (t = 0; t < num_threads; t++) {
        if (pthread_join(threads[t], NULL) != 0) {
            perror("pthread_join");
            exit(EXIT_FAILURE);
        }
    }
#endif

    // Delete the bintree after processing
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
    int       j, k, poly, n_points;
    int       label_values[65536], n_labels;
    int       *n_neighbours, **neighbours, *histo;
    signed long   i, min_val = 2147483647, max_val = -2147483647, longval;
    Point     point, scaled_point, center;
    Point     *new_points, poly_points[MAX_POINTS_PER_POLYGON];
    object_struct **objects, **scaled_objects;
    polygons_struct *scaled_target_sphere, *scaled_polygons_sphere;
    double      max_prob, val;
    float     weights[MAX_POINTS_PER_POLYGON];

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
                SCALE_POINT(scaled_point, poly_points[j], (double)weights[j]);
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

void
resample_spherical_surface_orig(polygons_struct *polygons,
               polygons_struct *poly_src_sphere,
               polygons_struct *resampled_source,
               double *input_values, double *output_values,
               int n_triangles)
{
    int    i, k, poly, n_points;
    int    *n_neighbours, **neighbours;
    Point  centre, point_on_src_sphere, scaled_point;
    Point  poly_points[MAX_POINTS_PER_POLYGON];
    Point  poly_points_src[MAX_POINTS_PER_POLYGON];
    Point  *new_points;
    float  weights[MAX_POINTS_PER_POLYGON];
    double sphereRadius, r, bounds[6];

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

    /* Calc. sphere center based on bounds of input (correct for shifts) */
    get_bounds(poly_src_sphere, bounds);
    fill_Point(centre, bounds[0]+bounds[1],
               bounds[2]+bounds[3], bounds[4]+bounds[5]);
  
    /*
     * Make radius slightly smaller to get sure that the
     * inner side of handles will be found as nearest point on the surface
     */
    sphereRadius *= 0.975;
    create_tetrahedral_sphere(&centre, sphereRadius, sphereRadius,
                  sphereRadius, n_triangles, resampled_source);

    create_polygons_bintree(poly_src_sphere,
                ROUND((float) poly_src_sphere->n_items * 0.5));

    ALLOC(new_points, resampled_source->n_points);
    if (input_values != NULL)
        ALLOC(output_values, resampled_source->n_points);

    for (i = 0; i < resampled_source->n_points; i++) {
        poly = find_closest_polygon_point(&resampled_source->points[i],
                          poly_src_sphere,
                          &point_on_src_sphere);
  
        n_points = get_polygon_points(poly_src_sphere, poly,
                        poly_points_src);
        get_polygon_interpolation_weights(&point_on_src_sphere,
                          n_points, poly_points_src,
                          weights);

        if (get_polygon_points(polygons, poly, poly_points) != n_points)
            fprintf(stderr,"map_point_between_polygons\n");

        fill_Point(new_points[i], 0.0, 0.0, 0.0);
        if (input_values != NULL)
            output_values[i] = 0.0;

        for (k = 0; k < n_points; k++) {
            SCALE_POINT(scaled_point, poly_points[k], (double)weights[k]);
            ADD_POINTS(new_points[i], new_points[i], scaled_point);
            if (input_values != NULL)
                output_values[i] += weights[k] *
                  input_values[polygons->indices[
                  POINT_INDEX(polygons->end_indices,poly,k)]];
        }
     }

    create_polygon_point_neighbours(resampled_source, TRUE, &n_neighbours,
                    &neighbours, NULL, NULL);

    for (i = 0; i < resampled_source->n_points; i++) {
        resampled_source->points[i] = new_points[i];
    }
  
    compute_polygon_normals(resampled_source);
    free(new_points);
    delete_the_bintree(&poly_src_sphere->bintree);
}

/* resample surface and values to space defined by tetrahedral sphere */
object_struct **
resample_surface(polygons_struct *surface, polygons_struct *sphere,
                 int n_triangles, double *invals, double *outvals)
{
    int       i, j, poly, n_points;
    Point     point, scaled_point, center;
    Point     *new_points, poly_points[MAX_POINTS_PER_POLYGON];
    object_struct **objects, *scaled_objects;
    polygons_struct *output_surface, *scaled_sphere;
    float      weights[MAX_POINTS_PER_POLYGON];

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
            SCALE_POINT(scaled_point, poly_points[j], (double)weights[j]);
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
