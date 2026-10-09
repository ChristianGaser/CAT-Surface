/* Rachel Yotter - rachel.yotter@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 *
 */

#include "CAT_Resample.h"
#include "CAT_SafeAlloc.h"

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
typedef struct
{
    int start_idx;
    int end_idx;
    polygons_struct *source_sphere;
    polygons_struct *target_sphere;
    double *invals;
    double *outvals;
} ThreadArgs_process_target_points;

// Define the struct to pass arguments to threads
typedef struct
{
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
static THREAD_RETURN process_resample_surface(void *args)
{
    ThreadArgs *thread_args = (ThreadArgs *)args;
    int i, k, poly, n_points;
    Point point_on_src_sphere, scaled_point;
    Point poly_points[MAX_POINTS_PER_POLYGON];
    Point poly_points_src[MAX_POINTS_PER_POLYGON];
    double weights[MAX_POINTS_PER_POLYGON];

    for (i = thread_args->start_idx; i < thread_args->end_idx; i++)
    {
        poly = find_closest_polygon_point(&thread_args->resampled_source->points[i],
                                          thread_args->poly_src_sphere,
                                          &point_on_src_sphere);

        n_points = get_polygon_points(thread_args->poly_src_sphere, poly, poly_points_src);
        get_polygon_interpolation_weights(&point_on_src_sphere, n_points, poly_points_src, weights);

        if (get_polygon_points(thread_args->polygons, poly, poly_points) != n_points)
        {
            fprintf(stderr, "map_point_between_polygons\n");
        }

        fill_Point(thread_args->new_points[i], 0.0, 0.0, 0.0);
        if (thread_args->input_values != NULL)
        {
            thread_args->output_values[i] = 0.0;
        }

        for (k = 0; k < n_points; k++)
        {
            SCALE_POINT(scaled_point, poly_points[k], (double)weights[k]);
            ADD_POINTS(thread_args->new_points[i], thread_args->new_points[i], scaled_point);
            if (thread_args->input_values != NULL)
            {
                thread_args->output_values[i] += weights[k] *
                                                 thread_args->input_values[thread_args->polygons->indices[POINT_INDEX(thread_args->polygons->end_indices, poly, k)]];
            }
        }
    }
    return NULL;
}
#endif

/**
 * \brief Resample surface data onto a tetrahedral sphere mesh with parallel processing.
 *
 * Maps surface data from a source polygonal mesh onto a tetrahedral sphere using
 * barycentric interpolation within matching surface polygons. Uses multithreading
 * on non-Windows systems for efficient batch processing of large meshes.
 *
 * Algorithm:
 *  1. Compute average radius of source sphere mesh
 *  2. Create tetrahedral sphere output mesh (0.975 * source radius)
 *  3. Build spatial index (bintree) for source sphere
 *  4. For each output vertex: find closest source polygon, interpolate surface position
 *  5. If input values provided, interpolate values via barycentric weights
 *  6. Normalize polygon normals for output mesh
 *
 * \param polygons         (in)  source surface mesh with data to map
 * \param poly_src_sphere  (in)  source spherical mesh aligned with polygons
 * \param resampled_source (out) output tetrahedral sphere mesh (modified in-place)
 * \param input_values     (in)  double[n_points]; surface data values (can be NULL)
 * \param output_values    (out) double[n_points]; interpolated values on output mesh
 * \param n_triangles      (in)  desired number of triangles in tetrahedral output sphere
 */
void resample_spherical_surface(polygons_struct *polygons,
                                polygons_struct *poly_src_sphere,
                                polygons_struct *resampled_source,
                                double *input_values, double *output_values,
                                int n_triangles)
{
    int i, t, k, num_threads = 8; // Adjust based on hardware
#if !defined(_WIN32) && !defined(_WIN64)
    THREAD_HANDLE threads[num_threads];
    ThreadArgs thread_args[num_threads];
#endif

    double sphereRadius = 0.0, r, bounds[6];

    // Compute average sphere radius
    for (i = 0; i < poly_src_sphere->n_points; i++)
    {
        r = 0.0;
        for (k = 0; k < 3; k++)
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

    create_polygons_bintree(poly_src_sphere, ROUND((float)poly_src_sphere->n_items * 0.5));

    Point *new_points;
    ALLOC(new_points, resampled_source->n_points);
    if (input_values != NULL)
        ALLOC(output_values, resampled_source->n_points);

    int total_points = resampled_source->n_points;
    int chunk_size = total_points / num_threads;
    int remainder = total_points % num_threads;

#if defined(_WIN32) || defined(_WIN64)
    // **Sequential Execution for Windows**
    for (i = 0; i < resampled_source->n_points; i++)
    {
        int poly, n_points;
        Point point_on_src_sphere, scaled_point;
        Point poly_points[MAX_POINTS_PER_POLYGON];
        Point poly_points_src[MAX_POINTS_PER_POLYGON];
        double weights[MAX_POINTS_PER_POLYGON];

        poly = find_closest_polygon_point(&resampled_source->points[i], poly_src_sphere, &point_on_src_sphere);

        n_points = get_polygon_points(poly_src_sphere, poly, poly_points_src);
        get_polygon_interpolation_weights(&point_on_src_sphere, n_points, poly_points_src, weights);

        if (get_polygon_points(polygons, poly, poly_points) != n_points)
        {
            fprintf(stderr, "map_point_between_polygons\n");
        }

        fill_Point(new_points[i], 0.0, 0.0, 0.0);
        if (input_values != NULL)
            output_values[i] = 0.0;

        for (k = 0; k < n_points; k++)
        {
            SCALE_POINT(scaled_point, poly_points[k], (double)weights[k]);
            ADD_POINTS(new_points[i], new_points[i], scaled_point);
            if (input_values != NULL)
                output_values[i] += weights[k] * input_values[polygons->indices[POINT_INDEX(polygons->end_indices, poly, k)]];
        }
    }
#else
    // **Parallel Execution for Non-Windows Systems**
    for (t = 0; t < num_threads; t++)
    {
        thread_args[t].start_idx = t * chunk_size;
        thread_args[t].end_idx = (t == num_threads - 1) ? (t + 1) * chunk_size + remainder : (t + 1) * chunk_size;
        thread_args[t].polygons = polygons;
        thread_args[t].poly_src_sphere = poly_src_sphere;
        thread_args[t].resampled_source = resampled_source;
        thread_args[t].input_values = input_values;
        thread_args[t].output_values = output_values;
        thread_args[t].new_points = new_points;

        if (pthread_create(&threads[t], NULL, process_resample_surface, &thread_args[t]) != 0)
        {
            perror("pthread_create");
            exit(EXIT_FAILURE);
        }
    }

    // Join all threads
    for (t = 0; t < num_threads; t++)
    {
        if (pthread_join(threads[t], NULL) != 0)
        {
            perror("pthread_join");
            exit(EXIT_FAILURE);
        }
    }
#endif

    // Copy new points to the resampled surface
    for (i = 0; i < resampled_source->n_points; i++)
    {
        resampled_source->points[i] = new_points[i];
    }

    compute_polygon_normals(resampled_source);
    free(new_points);
    delete_the_bintree(&poly_src_sphere->bintree);
}

/**
 * \brief Normalize source and target spheres to consistent scale and center.
 *
 * Corrects for differences in sphere size and position between source and target
 * by scaling both meshes to unit radius (radius 100.0 in this implementation) and
 * aligning their centers using bounds-based registration.
 *
 * Algorithm:
 *  1. Copy source and target sphere meshes
 *  2. Scale both copies to unit radius (100.0)
 *  3. Align centers using boundary-based correction
 *  4. Return scaled copies as output pointers
 *
 * \param source_sphere     (in)  source spherical mesh (unmodified)
 * \param target_sphere     (in)  target spherical mesh (unmodified)
 * \param out_source_sphere (out) scaled and centered source sphere copy
 * \param out_target_sphere (out) scaled and centered target sphere copy
 */
void correct_shift_scale_sphere(polygons_struct *source_sphere, polygons_struct *target_sphere,
                                polygons_struct **out_source_sphere, polygons_struct **out_target_sphere)
{
    int i;

    object_struct **source_sphere_objects, **target_sphere_objects;
    polygons_struct *scaled_source_sphere, *scaled_target_sphere;

    source_sphere_objects = (object_struct **)malloc(sizeof(object_struct *));
    *source_sphere_objects = create_object(POLYGONS);
    scaled_source_sphere = get_polygons_ptr(*source_sphere_objects);
    copy_polygons(source_sphere, scaled_source_sphere);
    for (i = 0; i < scaled_source_sphere->n_points; i++)
        set_vector_length(&scaled_source_sphere->points[i], 100.0);

    target_sphere_objects = (object_struct **)malloc(sizeof(object_struct *));
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

/**
 * \brief Sum-preserving (areal) resampling of per-vertex values between spheres.
 *
 * Extensive measures such as the local surface area have to keep their sum
 * when they are moved to another mesh.  Every target vertex is located in the
 * source sphere (forward mapping), and h[s] sums the barycentric weights with
 * which the target vertices reference source vertex s:
 *  - h[s] >= 1: s is divided among the target vertices that reference it,
 *    each receiving weight / h[s] of its value.
 *  - h[s] <  1: the target mesh is locally coarser than the source mesh.  The
 *    referencing target vertices receive weight times the value, and the
 *    remaining (1 - h[s]) of the value goes to the vertices of the target
 *    triangle that contains s (reverse mapping, barycentric weights).
 * Each source vertex therefore contributes exactly its value, and the total is
 * preserved also for source vertices that no target vertex reaches.  Without
 * the reverse mapping these vertices were lost, which for a 32k target mesh
 * removed 5-25% of the native surface area, the more the larger the brain.
 * This is the barycentric counterpart of the forward and reverse
 * nearest-neighbour mapping with Jacobian correction in FreeSurfer
 * (mri_surf2surf --mapmethod nnfr --jac).
 *
 * Both spheres have to be aligned and of the same radius.  Bintrees that
 * already exist are used and kept.
 *
 * \param source_sphere (in)  source spherical mesh
 * \param target_sphere (in)  target spherical mesh
 * \param invals        (in)  double[source n_points]; source vertex values
 * \param outvals       (out) double[target n_points]; resampled values
 */
static void
resample_values_areal(polygons_struct *source_sphere, polygons_struct *target_sphere,
                      double *invals, double *outvals)
{
    int i, j, poly, n_points, idx;
    int own_src_tree = 0, own_trg_tree = 0;
    Point point, poly_points[MAX_POINTS_PER_POLYGON];
    double weights[MAX_POINTS_PER_POLYGON];
    double *hits, rest;

    if (source_sphere->bintree == NULL)
    {
        create_polygons_bintree(source_sphere,
                                ROUND((double)source_sphere->n_items * 0.5));
        own_src_tree = 1;
    }

    hits = SAFE_MALLOC(double, source_sphere->n_points);
    memset(hits, 0, sizeof(double) * source_sphere->n_points);

    /* forward mapping: weights with which the target vertices reference
       each source vertex */
    for (i = 0; i < target_sphere->n_points; i++)
    {
        poly = find_closest_polygon_point(&target_sphere->points[i], source_sphere, &point);
        n_points = get_polygon_points(source_sphere, poly, poly_points);
        get_polygon_interpolation_weights(&point, n_points, poly_points, weights);

        for (j = 0; j < n_points; j++)
            hits[source_sphere->indices[POINT_INDEX(source_sphere->end_indices, poly, j)]] += weights[j];
    }

    for (i = 0; i < target_sphere->n_points; i++)
    {
        poly = find_closest_polygon_point(&target_sphere->points[i], source_sphere, &point);
        n_points = get_polygon_points(source_sphere, poly, poly_points);
        get_polygon_interpolation_weights(&point, n_points, poly_points, weights);

        outvals[i] = 0.0;
        for (j = 0; j < n_points; j++)
        {
            idx = source_sphere->indices[POINT_INDEX(source_sphere->end_indices, poly, j)];
            if (hits[idx] >= 1.0)
                outvals[i] += weights[j] * invals[idx] / hits[idx];
            else
                outvals[i] += weights[j] * invals[idx];
        }
    }

    /* reverse mapping: the part of a source vertex that the target vertices
       do not reference goes to the target triangle that contains it */
    if (target_sphere->bintree == NULL)
    {
        create_polygons_bintree(target_sphere,
                                ROUND((double)target_sphere->n_items * 0.5));
        own_trg_tree = 1;
    }

    for (i = 0; i < source_sphere->n_points; i++)
    {
        if (hits[i] >= 1.0)
            continue;

        rest = (1.0 - hits[i]) * invals[i];
        poly = find_closest_polygon_point(&source_sphere->points[i], target_sphere, &point);
        n_points = get_polygon_points(target_sphere, poly, poly_points);
        get_polygon_interpolation_weights(&point, n_points, poly_points, weights);

        for (j = 0; j < n_points; j++)
            outvals[target_sphere->indices[POINT_INDEX(target_sphere->end_indices, poly, j)]] += weights[j] * rest;
    }

    if (own_trg_tree)
        delete_the_bintree(&target_sphere->bintree);
    if (own_src_tree)
        delete_the_bintree(&source_sphere->bintree);
    FREE(hits);
}

// Thread function (used only on non-Windows systems)
#if !defined(_WIN32) && !defined(_WIN64)
static THREAD_RETURN process_target_points(void *args)
{
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
    double weights[MAX_POINTS_PER_POLYGON];

    for (i = start; i < end; i++)
    {
        poly = find_closest_polygon_point(&target_sphere->points[i], source_sphere, &point);

        n_points = get_polygon_points(source_sphere, poly, poly_points);
        get_polygon_interpolation_weights(&point, n_points, poly_points, weights);

        outvals[i] = 0.0;
        for (j = 0; j < n_points; j++)
        {
            outvals[i] += (double)weights[j] * invals[source_sphere->indices[POINT_INDEX(source_sphere->end_indices, poly, j)]];
        }
    }

    return NULL;
}
#endif

/**
 * \brief Interpolate scalar values from source sphere vertices to target sphere vertices.
 *
 * Maps scalar values (e.g., cortical thickness) from source sphere mesh to target
 * sphere mesh using barycentric interpolation within closest polygons. Assumes both
 * spheres are already aligned (use correct_shift_scale_sphere if alignment needed).
 * Supports optional sum-preserving areal interpolation (see resample_values_areal).
 *
 * Algorithm:
 *  1. Build spatial index (bintree) for source sphere if needed
 *  2. For each target vertex: find closest polygon in source sphere
 *  3. Compute barycentric coordinates within that polygon
 *  4. Interpolate value as weighted sum of vertex values
 *  5. If areal_interpolation: distribute the source values with forward and
 *     reverse mapping so that their sum is preserved
 *
 * \param source_sphere       (in)  source spherical mesh with data
 * \param target_sphere       (in)  target spherical mesh (aligned)
 * \param invals              (in)  double[n_points]; source vertex values
 * \param outvals             (out) double[n_points]; interpolated target vertex values
 * \param areal_interpolation (in)  0/1; apply sum-preserving areal interpolation if 1
 */
void resample_values_sphere_noscale(polygons_struct *source_sphere,
                                    polygons_struct *target_sphere,
                                    double *invals, double *outvals,
                                    int areal_interpolation)
{
    int t, num_threads = 8; // Number of threads (adjust based on system hardware)
#if !defined(_WIN32) && !defined(_WIN64)
    THREAD_HANDLE threads[num_threads];
    ThreadArgs_process_target_points thread_args[num_threads];
#endif

    // Create bintree if it doesn't already exist
    if (source_sphere->bintree == NULL)
    {
        create_polygons_bintree(source_sphere,
                                ROUND((double)source_sphere->n_items * 0.5));
    }

    if (areal_interpolation && invals != NULL)
    {
        resample_values_areal(source_sphere, target_sphere, invals, outvals);
        delete_the_bintree(&source_sphere->bintree);
        return;
    }

    // Divide the workload among threads
    int total_points = target_sphere->n_points;
    int chunk_size = total_points / num_threads;
    int remainder = total_points % num_threads;

#if defined(_WIN32) || defined(_WIN64)
    // **Sequential Execution for Windows**
    int i, j;
    for (i = 0; i < total_points; i++)
    {
        int poly, n_points;
        Point point;
        Point poly_points[MAX_POINTS_PER_POLYGON];
        double weights[MAX_POINTS_PER_POLYGON];

        poly = find_closest_polygon_point(&target_sphere->points[i], source_sphere, &point);

        n_points = get_polygon_points(source_sphere, poly, poly_points);
        get_polygon_interpolation_weights(&point, n_points, poly_points, weights);

        outvals[i] = 0.0;
        for (j = 0; j < n_points; j++)
        {
            outvals[i] += (double)weights[j] * invals[source_sphere->indices[POINT_INDEX(source_sphere->end_indices, poly, j)]];
        }
    }
#else
    // **Parallel Execution for Non-Windows Systems**
    for (t = 0; t < num_threads; t++)
    {
        thread_args[t].start_idx = t * chunk_size;
        thread_args[t].end_idx = (t == num_threads - 1) ? (t + 1) * chunk_size + remainder
                                                        : (t + 1) * chunk_size;
        thread_args[t].source_sphere = source_sphere;
        thread_args[t].target_sphere = target_sphere;
        thread_args[t].invals = invals;
        thread_args[t].outvals = outvals;

        if (pthread_create(&threads[t], NULL, process_target_points, &thread_args[t]) != 0)
        {
            perror("pthread_create");
            exit(EXIT_FAILURE);
        }
    }

    // Join all threads
    for (t = 0; t < num_threads; t++)
    {
        if (pthread_join(threads[t], NULL) != 0)
        {
            perror("pthread_join");
            exit(EXIT_FAILURE);
        }
    }
#endif

    // Delete the bintree after processing
    delete_the_bintree(&source_sphere->bintree);
}

/**
 * \brief Interpolate scalar values from source to target sphere with optional alignment.
 *
 * High-level wrapper for value resampling that optionally normalizes sphere scale
 * and position before interpolation. Useful when source and target spheres may have
 * different radii or centers.
 *
 * Algorithm:
 *  1. If scale_and_shift: correct sphere alignment using correct_shift_scale_sphere
 *  2. Call resample_values_sphere_noscale for actual interpolation
 *  3. If not scaling: directly call resample_values_sphere_noscale
 *
 * \param source_sphere        (in)  source spherical mesh
 * \param target_sphere        (in)  target spherical mesh
 * \param invals               (in)  double[n_points]; source vertex values
 * \param outvals              (out) double[n_points]; interpolated target values
 * \param scale_and_shift      (in)  0/1; normalize sphere alignment if 1
 * \param areal_interpolation  (in)  0/1; use sum-preserving areal interpolation if 1
 */
void resample_values_sphere(polygons_struct *source_sphere, polygons_struct *target_sphere,
                            double *invals, double *outvals, int scale_and_shift, int areal_interpolation)
{
    int i;
    polygons_struct *scaled_source_sphere, *scaled_target_sphere;

    if (scale_and_shift)
    { /* correct shifts and vector length */
        correct_shift_scale_sphere(source_sphere, target_sphere, &scaled_source_sphere,
                                   &scaled_target_sphere);

        resample_values_sphere_noscale(scaled_source_sphere, scaled_target_sphere,
                                       invals, outvals, areal_interpolation);
    }
    else
        resample_values_sphere_noscale(source_sphere, target_sphere, invals, outvals, areal_interpolation);
}

/**
 * \brief Resample surface mesh and its data to a target spherical coordinate system.
 *
 * Maps a surface mesh and optional scalar data onto a target spherical mesh, using
 * per-triangle label interpolation when multiple labels exist within a triangle.
 * Handles both explicit sphere meshes and implicit tetrahedral sphere definitions.
 *
 * Algorithm:
 *  1. Create or copy tetrahedral sphere matching target geometry
 *  2. Build spatial index for source sphere
 *  3. For each target vertex: find closest source polygon
 *  4. Interpolate vertex position via barycentric weighting
 *  5. If values provided: interpolate with optional label-aware refinement
 *  6. Free spatial index and return resampled surface
 *
 * \param polygons            (in)  source surface mesh with data to map
 * \param polygons_sphere     (in)  source sphere (can be NULL for implicit tetrahedral)
 * \param target_sphere       (in)  target spherical coordinate system
 * \param input_values        (in)  double[n_points]; scalar data on source surface
 * \param output_values       (out) double[n_points]; interpolated values on output surface
 * \param label_interpolation (in)  0/1; use per-triangle label refinement if 1
 * \param areal_interpolation (in)  0/1; sum-preserving interpolation of values
 *                                  (e.g. surface area) if 1, see resample_values_areal
 * \return                          object_struct** containing resampled surface mesh
 */
object_struct **
resample_surface_to_target_sphere(polygons_struct *polygons, polygons_struct *polygons_sphere, polygons_struct *target_sphere,
                                  double *input_values, double *output_values, int label_interpolation,
                                  int areal_interpolation)
{
    int j, k, poly, n_points, use_areal;
    int *n_neighbours, **neighbours;
    signed long i;
    Point point, scaled_point, center;
    Point *new_points, poly_points[MAX_POINTS_PER_POLYGON];
    object_struct **objects, **scaled_objects;
    polygons_struct *scaled_target_sphere, *scaled_polygons_sphere;
    double max_prob, val;
    double weights[MAX_POINTS_PER_POLYGON];

    objects = SAFE_MALLOC(object_struct *, 1);
    *objects = create_object(POLYGONS);
    scaled_polygons_sphere = get_polygons_ptr(*objects);

    /* Per-triangle label interpolation: collect local non-zero labels only */
    /* if no source sphere is defined a tetrahedral topology of the surface is assumed
       where the corresponding sphere can be simply estimated by its tetrahedral topology */
    if (polygons_sphere == NULL)
    {
        fill_Point(center, 0.0, 0.0, 0.0);
        create_tetrahedral_sphere(&center, 100.0, 100.0, 100.0,
                                  polygons->n_items, scaled_polygons_sphere);
    }
    else
    {
        copy_polygons(polygons_sphere, scaled_polygons_sphere);
        for (i = 0; i < scaled_polygons_sphere->n_points; i++)
            set_vector_length(&scaled_polygons_sphere->points[i], 100.0);
    }

    /* scale the re-parameterizing sphere... also the returned object. */
    scaled_objects = SAFE_MALLOC(object_struct *, 1);
    *scaled_objects = create_object(POLYGONS);
    scaled_target_sphere = get_polygons_ptr(*scaled_objects);
    copy_polygons(target_sphere, scaled_target_sphere);

    for (i = 0; i < scaled_target_sphere->n_points; i++)
        set_vector_length(&scaled_target_sphere->points[i], 100.0);

    /* correct sphere center based on bounds of target */
    correct_bounds_to_target(scaled_polygons_sphere, scaled_target_sphere);

    create_polygons_bintree(scaled_polygons_sphere,
                            ROUND((double)scaled_polygons_sphere->n_items * 0.5));

    new_points = SAFE_MALLOC(Point, scaled_target_sphere->n_points);

    use_areal = areal_interpolation && input_values != NULL && !label_interpolation;

    for (i = 0; i < scaled_target_sphere->n_points; i++)
    {
        poly = find_closest_polygon_point(&scaled_target_sphere->points[i],
                                          scaled_polygons_sphere, &point);

        n_points = get_polygon_points(scaled_polygons_sphere, poly, poly_points);
        get_polygon_interpolation_weights(&point, n_points, poly_points, weights);

        if (polygons != NULL)
            if (get_polygon_points(polygons, poly, poly_points) != n_points)
                fprintf(stderr, "map_point_between_polygons\n");

        fill_Point(new_points[i], 0.0, 0.0, 0.0);

        if (input_values != NULL)
            output_values[i] = 0.0;

        /* resample mesh if defined */
        if (polygons != NULL)
        {
            for (j = 0; j < n_points; j++)
            {
                SCALE_POINT(scaled_point, poly_points[j], (double)weights[j]);
                ADD_POINTS(new_points[i], new_points[i], scaled_point);
            }
        }

        /* apply interpolation for each ROI label separately (local per-triangle) */
        if (label_interpolation)
        {
            if (input_values != NULL)
            {
                int cand_labels[MAX_POINTS_PER_POLYGON];
                double cand_probs[MAX_POINTS_PER_POLYGON];
                int n_cand = 0;

                /* accumulate weights per unique non-zero label from the triangle */
                for (j = 0; j < n_points; j++)
                {
                    int idx = scaled_polygons_sphere->indices[POINT_INDEX(scaled_polygons_sphere->end_indices, poly, j)];
                    int lbl = (int)input_values[idx];
                    if (lbl == 0)
                        continue; /* exclude background from candidates */

                    int found = -1;
                    for (k = 0; k < n_cand; k++)
                    {
                        if (cand_labels[k] == lbl)
                        {
                            found = k;
                            break;
                        }
                    }
                    if (found >= 0)
                    {
                        cand_probs[found] += weights[j];
                    }
                    else
                    {
                        cand_labels[n_cand] = lbl;
                        cand_probs[n_cand] = weights[j];
                        n_cand++;
                    }
                }

                /* default remains 0.0 when no non-zero labels present */
                if (n_cand > 0)
                {
                    /* sort candidates ascending by label value to preserve tie-breaking */
                    for (k = 0; k < n_cand - 1; k++)
                    {
                        int m;
                        for (m = k + 1; m < n_cand; m++)
                        {
                            if (cand_labels[m] < cand_labels[k])
                            {
                                int tl = cand_labels[k];
                                double tp = cand_probs[k];
                                cand_labels[k] = cand_labels[m];
                                cand_probs[k] = cand_probs[m];
                                cand_labels[m] = tl;
                                cand_probs[m] = tp;
                            }
                        }
                    }

                    /* choose label with max weighted probability; strict '>' keeps smallest on ties */
                    max_prob = 0.0;
                    for (k = 0; k < n_cand; k++)
                    {
                        val = cand_probs[k];
                        if (max_prob < val)
                        {
                            max_prob = val;
                            output_values[i] = (double)cand_labels[k];
                        }
                    }
                }
            }
        }
        else if (input_values != NULL && !use_areal)
        {
            /* areal values are resampled separately below */
            output_values[i] = 0.0;
            for (j = 0; j < n_points; j++)
            {
                int idx = scaled_polygons_sphere->indices[POINT_INDEX(scaled_polygons_sphere->end_indices, poly, j)];
                output_values[i] += weights[j] * input_values[idx];
            }
        }
    }

    /* sum-preserving resampling of values; this needs the target sphere, so it
       has to be done before its points are replaced by the resampled surface */
    if (use_areal)
        resample_values_areal(scaled_polygons_sphere, scaled_target_sphere,
                              input_values, output_values);

    if (polygons != NULL)
    {
        free(scaled_target_sphere->points);
        scaled_target_sphere->points = new_points;
    }

    /* no global histogram used; nothing to free for labels */

    compute_polygon_normals(scaled_target_sphere);

    delete_the_bintree(&scaled_polygons_sphere->bintree);
    delete_object_list(1, objects);

    return (scaled_objects);
}

/**
 * \brief Resample surface mesh to a tetrahedral sphere using barycentric interpolation.
 *
 * Maps a surface mesh and optional scalar data onto a new tetrahedral spherical
 * mesh using barycentric coordinate interpolation. Creates output sphere with
 * specified triangulation and scales to match source sphere radius.
 *
 * Algorithm:
 *  1. Copy source sphere and scale to unit radius (100.0)
 *  2. Create new tetrahedral sphere output mesh
 *  3. Build spatial index (bintree) for source sphere
 *  4. For each output vertex: find closest source polygon
 *  5. Compute barycentric weights and interpolate vertex position
 *  6. If values provided: interpolate values using same weights
 *
 * \param surface     (in)  source surface mesh with data
 * \param sphere      (in)  source or reference spherical mesh
 * \param n_triangles (in)  number of triangles for output tetrahedral sphere
 * \param invals      (in)  double[n_points]; surface data values (can be NULL)
 * \param outvals     (out) double[n_points]; interpolated values on output sphere
 * \return                  object_struct** containing resampled surface and sphere
 */
object_struct **
resample_surface(polygons_struct *surface, polygons_struct *sphere,
                 int n_triangles, double *invals, double *outvals)
{
    int i, j, poly, n_points;
    Point point, scaled_point, center;
    Point *new_points, poly_points[MAX_POINTS_PER_POLYGON];
    object_struct **objects, *scaled_objects;
    polygons_struct *output_surface, *scaled_sphere;
    double weights[MAX_POINTS_PER_POLYGON];

    scaled_objects = create_object(POLYGONS);
    scaled_sphere = get_polygons_ptr(scaled_objects);
    copy_polygons(sphere, scaled_sphere);

    translate_to_center_of_mass(scaled_sphere);
    for (i = 0; i < scaled_sphere->n_points; i++)
        set_vector_length(&scaled_sphere->points[i], 100.0);

    objects = (object_struct **)malloc(sizeof(object_struct *));
    *objects = create_object(POLYGONS);
    output_surface = get_polygons_ptr(*objects);
    fill_Point(center, 0.0, 0.0, 0.0);
    create_tetrahedral_sphere(&center, 100.0, 100.0, 100.0, n_triangles,
                              output_surface);

    create_polygons_bintree(sphere, ROUND((double)sphere->n_items * 0.5));

    new_points = (Point *)malloc(sizeof(Point) * output_surface->n_points);

    if (invals != NULL)
        ALLOC(outvals, output_surface->n_points);

    for (i = 0; i < output_surface->n_points; i++)
    {
        poly = find_closest_polygon_point(&output_surface->points[i], sphere, &point);
        n_points = get_polygon_points(sphere, poly, poly_points);
        get_polygon_interpolation_weights(&point, n_points, poly_points, weights);

        fill_Point(new_points[i], 0.0, 0.0, 0.0);
        if (invals != NULL)
            outvals[i] = 0.0;

        for (j = 0; j < n_points; j++)
        {
            SCALE_POINT(scaled_point, poly_points[j], (double)weights[j]);
            ADD_POINTS(new_points[i], new_points[i], scaled_point);
            if (invals != NULL)
                outvals[i] += weights[j] * invals[surface->indices[POINT_INDEX(surface->end_indices, poly, j)]];
        }
    }

    for (i = 0; i < output_surface->n_points; i++)
        output_surface->points[i] = new_points[i];

    compute_polygon_normals(output_surface);
    free(new_points);
    delete_the_bintree(&sphere->bintree);
    return (objects);
}
