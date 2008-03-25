/* Christian Gaser - christian.gaser@uni-jena.de                             */
/* Department of Psychiatry                                                  */
/* University of Jena                                                        */
/*                                                                           */
/* Copyright Christian Gaser, University of Jena.                            */

#include  <volume_io/internal_volume_io.h>
#include  <CAT_Deform.h>

private  void  perturb_points(
    polygons_struct              *polygons,
    Point                        new_points[],
    Real                         fractional_step,
    Real                         max_step,
    Real                         max_search_distance,
    int                          degrees_continuity,
    deform_data_struct           *deform_data,
    boundary_definition_struct   *boundary_def,
    deformation_model_struct     *deformation_model,
    Real                         movement_threshold,
    float                        prev_movements[],
    deform_stats                 *stats );
private  Real  one_iteration_polygons(
    polygons_struct   *polygons,
    deform_struct     *deform_parms,
    int               iteration );
private  void    check_polygons_shape_integrity(
    polygons_struct *polygons,
    Point           new_points[] );

public  void  deform_polygons(
    polygons_struct   *polygons,
    deform_struct     *deform_parms )
{
    int                iteration, countdown, countdown2;
    Real               avg_error, prev_avg_error, rate;

    iteration = 0;
    prev_avg_error = 1e10;
    countdown = 0;
    countdown2 = 0;
    
    print("\n");
    do
    {
        ++iteration;

        avg_error = one_iteration_polygons( polygons, deform_parms, iteration );
        rate = (prev_avg_error - avg_error)/(prev_avg_error + avg_error);
		if (rate < deform_parms->stop_threshold)
			{countdown = countdown + 1;}
		else {countdown = 0;}
        if (prev_avg_error < avg_error)
			{countdown2 = countdown2 + 1;}
		else {countdown2 = 0;}
        prev_avg_error = avg_error;
    }
    while( (countdown < 4 || countdown2 < 4) &&
           countdown < 10 &&
           iteration < deform_parms->max_iterations );
}

void  deform_polygons_one_iteration(
    polygons_struct   *polygons,
    deform_struct     *deform_parms,
    int               iteration )
{
    (void) one_iteration_polygons( polygons, deform_parms, iteration );
}

private  Real  one_iteration_polygons(
    polygons_struct   *polygons,
    deform_struct     *deform_parms,
    int               iteration )
{
    int                i;
    Point              *new_points, *tmp;
    deform_stats       stats;
    

    if( polygons->n_points <= 0 )
        return( 0.0 );

    if( !check_correct_deformation_polygons( polygons,
                                        &deform_parms->deformation_model ) )
        return( 0.0 );

    if( deform_parms->n_movements_alloced != polygons->n_points )
    {
        if( deform_parms->n_movements_alloced > 0 )
            FREE( deform_parms->prev_movements );

        deform_parms->n_movements_alloced = polygons->n_points;

        ALLOC( deform_parms->prev_movements, polygons->n_points );

        for_less( i, 0, polygons->n_points )
            deform_parms->prev_movements[i] =
                      (float) deform_parms->movement_threshold + 1.0f;
    }

    ALLOC( new_points, polygons->n_points );

    check_polygons_neighbours_computed( polygons );

    initialize_deform_stats( &stats );

    /* --- every 50 iterations do all points */

    if( iteration % 50 == 0 )
    {
        for_less( i, 0, polygons->n_points )
        {
            deform_parms->prev_movements[i] =
                         (float) deform_parms->movement_threshold + 1.0f;
        }
    }

    perturb_points( polygons, new_points,
                    deform_parms->fractional_step,
                    deform_parms->max_step,
                    deform_parms->max_search_distance,
                    deform_parms->degrees_continuity,
                    &deform_parms->deform_data,
                    &deform_parms->boundary_definition,
                    &deform_parms->deformation_model,
                    deform_parms->movement_threshold,
                    deform_parms->prev_movements,
                    &stats );

    print( "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b");
    print( "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b");
    print("Iteration %d:\t",iteration );
    print("avg: %3.4f\tmax: %3.4f\n",stats.average/polygons->n_points,stats.maximum);
//    print_deform_stats( &stats, polygons->n_points );
    if( iteration % 20 == 0 )
        (void) flush_file( stdout );

    check_polygons_shape_integrity( polygons, new_points );

    tmp = polygons->points;
    polygons->points = new_points;
    new_points = tmp;

    FREE( new_points );

    if( deform_parms->movement_threshold < 0.0 )
    {
        deform_parms->movement_threshold = stats.average / 10.0;
        if( deform_parms->movement_threshold > 0.1 )
            deform_parms->movement_threshold = 0.1;
    }

    return( stats.average );
}

#define  MAX_NEIGHBOURS  2000

private  void  get_polygon_equilibrium_point(
    polygons_struct              *polygons,
    int                          poly,
    int                          vertex_index,
    Real                         curvature_factors[],
    Real                         max_search_distance,
    int                          degrees_continuity,
    deform_data_struct           *deform_data,
    boundary_definition_struct   *boundary_def,
    deformation_model_struct     *deformation_model,
    Point                        *equilibrium_point )
{
    int              point_index;
    int              n_neighbours, neighbours[MAX_NEIGHBOURS];
    BOOLEAN          interior_flag, found;
    Real             curvature_factor, model_distance, boundary_distance;
    Real             base_length;
    Point            centroid, model_point, search_origin;
    Vector           normal, pos_model_dir, neg_model_dir;

    point_index = polygons->indices[
                          POINT_INDEX(polygons->end_indices,poly,vertex_index)];

    n_neighbours = get_subsampled_neighbours_of_point( deformation_model,
                                 polygons, poly,
                                 vertex_index, neighbours, MAX_NEIGHBOURS,
                                 &interior_flag );

    compute_points_centroid_and_normal( polygons, point_index,
                              n_neighbours, neighbours, &centroid,
                              &normal, &base_length, &curvature_factor );

    get_model_point( deformation_model, polygons->points,
                     point_index, n_neighbours, neighbours, curvature_factors,
                     &centroid, &normal, base_length, &model_point );

    compute_model_dirs( &centroid, &normal, base_length, &model_point,
                        &model_distance, &search_origin,
                        &pos_model_dir, &neg_model_dir );

    found = find_boundary_in_direction( deform_data->volume,
                                        deform_data->label_volume,
                                        NULL, NULL, NULL,
                                        model_distance, &search_origin,
                                        &pos_model_dir, &neg_model_dir,
                                        max_search_distance,
                                        max_search_distance,
                                        degrees_continuity,
                                        boundary_def,
                                        &boundary_distance );

    compute_equilibrium_point( point_index, found, boundary_distance,
                               base_length, model_distance,
                               &pos_model_dir, &neg_model_dir,
                               &centroid, deformation_model, equilibrium_point);
}

private  BOOLEAN  counter_clockwise_neighbours(
    Point        *centroid,
    Vector       *normal,
    Point        points[],
    int          n_neighbours,
    int          neighbours[],
    Smallest_int point_error[] )
{
    Vector    to_neighbour, prev_to_neighbour, up, offset;
    Real      len;
    int       i;
    BOOLEAN   counter_clockwise;

    counter_clockwise = TRUE;

    fill_Vector( to_neighbour, 0.0, 0.0, 0.0 );
    for_less( i, 0, n_neighbours + 1 )
    {
        prev_to_neighbour = to_neighbour;
        SUB_VECTORS( to_neighbour, points[neighbours[i%n_neighbours]],
                     *centroid );
        len = DOT_VECTORS( to_neighbour, *normal );
        SCALE_VECTOR( offset, *normal, len );
        SUB_VECTORS( to_neighbour, to_neighbour, offset );

        if( i != 0 )
        {
            CROSS_VECTORS( up, prev_to_neighbour, to_neighbour );
            if( DOT_VECTORS( up, *normal ) < 0.0 )
            {
                ++point_error[neighbours[i%n_neighbours]];
                ++point_error[neighbours[(i-1+n_neighbours)%n_neighbours]];
                counter_clockwise = FALSE;
            }
        }
    }

    return( counter_clockwise );
}

private  void  check_polygons_shape_integrity(
    polygons_struct              *polygons,
    Point                        new_points[] )
{
    Smallest_int     *point_done;
    int              vertex_index, point_index, poly, size;
    Point            *centroids;
    Vector           normal;
    progress_struct  progress;
    Real             base_length, curvature_factor;
    int              n_neighbours, neighbours[MAX_NEIGHBOURS];
    BOOLEAN          interior_flag;
    Smallest_int     *point_error;
#ifdef  DEBUG
    int              n_errors, n_bad_points;
#endif

    ALLOC( point_done, polygons->n_points );
    ALLOC( point_error, polygons->n_points );
    ALLOC( centroids, polygons->n_points );

    for_less( point_index, 0, polygons->n_points )
    {
        point_done[point_index] = FALSE;
        point_error[point_index] = 0;
    }

    initialize_progress_report( &progress, TRUE, polygons->n_items,
                                "Checking Integrity" );

    for_less( poly, 0, polygons->n_items )
    {
        size = GET_OBJECT_SIZE( *polygons, poly );

        for_less( vertex_index, 0, size )
        {
            point_index = polygons->indices[
                      POINT_INDEX(polygons->end_indices,poly,vertex_index)];

            if( !point_done[point_index] )
            {
                point_done[point_index] = TRUE;

                compute_polygon_point_centroid( polygons, poly,
                                vertex_index, point_index,
                                &centroids[point_index], &normal, &base_length,
                                &curvature_factor );
                n_neighbours = get_neighbours_of_point( polygons,
                                 poly, vertex_index, neighbours, MAX_NEIGHBOURS,
                                 &interior_flag );

#ifdef CHECK_CLOCKWISE_NEIGHBOURS
                if( !counter_clockwise_neighbours( &centroids[point_index],
                            &normal,
                            new_points, n_neighbours, neighbours, point_error ))
                {
                    point_error[point_index] = TRUE;
                    ++n_errors;
                }
#else
                (void) counter_clockwise_neighbours( &centroids[point_index],
                            &normal, new_points, n_neighbours,
                            neighbours, point_error );
#endif
            }
        }

        update_progress_report( &progress, poly+1 );
    }

    terminate_progress_report( &progress );

#ifdef DEBUG
    n_errors = 0;
    n_bad_points = 0;
#endif
    for_less( point_index, 0, polygons->n_points )
    {
        if( point_error[point_index] > 0 )
        {
#ifdef DEBUG
            ++n_errors;
            n_bad_points += point_error[point_index];
            if( n_errors < 10 )
                print( " %d", point_index );
#endif
            new_points[point_index] = centroids[point_index];
        }
    }

#ifdef DEBUG
    if( n_errors > 0 )
        print( ": Shape errors %d/%d\n", n_errors, n_bad_points );
#endif

    FREE( point_error );
    FREE( centroids );
    FREE( point_done );
}

private  void  perturb_points(
    polygons_struct              *polygons,
    Point                        new_points[],
    Real                         fractional_step,
    Real                         max_step,
    Real                         max_search_distance,
    int                          degrees_continuity,
    deform_data_struct           *deform_data,
    boundary_definition_struct   *boundary_def,
    deformation_model_struct     *deformation_model,
    Real                         movement_threshold,
    float                        prev_movements[],
    deform_stats                 *stats )
{
    Real             *curvature_factors;
    Smallest_int     *point_done;
    int              vertex_index, point_index, poly, size, n1, n2;
    Point            centroid;
    Vector           normal;
    progress_struct  progress;
    Point            equilibrium_point;
    Real             dist_to_equil, base_length;
    float            *movements;

    ALLOC( curvature_factors, polygons->n_points );
    ALLOC( point_done, polygons->n_points );

    if( deformation_model_includes_average( deformation_model ) )
    {
        for_less( point_index, 0, polygons->n_points )
            point_done[point_index] = FALSE;

        initialize_progress_report( &progress, TRUE, polygons->n_items,
                                    "Computing Curvatures" );

        for_less( poly, 0, polygons->n_items )
        {
            size = GET_OBJECT_SIZE( *polygons, poly );

            for_less( vertex_index, 0, size )
            {
                point_index = polygons->indices[
                          POINT_INDEX(polygons->end_indices,poly,vertex_index)];

                if( !point_done[point_index] )
                {
                    point_done[point_index] = TRUE;

                    compute_polygon_point_centroid( polygons, poly,
                                      vertex_index, point_index,
                                      &centroid, &normal, &base_length,
                                      &curvature_factors[point_index] );
                }
            }

            update_progress_report( &progress, poly+1 );
        }

        terminate_progress_report( &progress );
    }

    ALLOC( movements, polygons->n_points );

    for_less( point_index, 0, polygons->n_points )
    {
        movements[point_index] = 0.0f;
        point_done[point_index] = FALSE;
        new_points[point_index] = polygons->points[point_index];
    }

    for_less( poly, 0, polygons->n_items )
    {
        size = GET_OBJECT_SIZE( *polygons, poly );

        for_less( vertex_index, 0, size )
        {
            point_index = polygons->indices[
                          POINT_INDEX(polygons->end_indices,poly,vertex_index)];

            if( point_done[point_index] )
                continue;

            n1 = polygons->indices[
                       POINT_INDEX(
                       polygons->end_indices,poly,(vertex_index+1)%size )];
            n2 = polygons->indices[
                       POINT_INDEX(
                       polygons->end_indices,poly,(vertex_index-1+size)%size )];

            if( (Real) prev_movements[point_index] < movement_threshold &&
                (Real) prev_movements[n1] < movement_threshold &&
                (Real) prev_movements[n2] < movement_threshold )
                continue;

            point_done[point_index] = TRUE;

            get_polygon_equilibrium_point( polygons, poly, vertex_index,
                                           curvature_factors,
                                           max_search_distance,
                                           degrees_continuity,
                                           deform_data, boundary_def,
                                           deformation_model,
                                           &equilibrium_point );

            dist_to_equil = deform_point( point_index, polygons->points,
                                      &equilibrium_point,
                                      fractional_step, max_step,
                                      deformation_model->position_constrained,
                                      deformation_model->max_position_offset,
                                      deformation_model->original_positions,
                                      &new_points[point_index] );

            movements[point_index] = (float) distance_between_points(
                                         &polygons->points[point_index],
                                         &new_points[point_index] );

            record_error_in_deform_stats( stats, dist_to_equil );
        }
    }

    for_less( point_index, 0, polygons->n_points )
    {
        if( !point_done[point_index] )
            record_error_in_deform_stats( stats, 0.0 );
        prev_movements[point_index] = movements[point_index];
    }

    FREE( movements );
    FREE( point_done );
    FREE( curvature_factors );
}
