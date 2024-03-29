/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 *
 */

#include "CAT_DeformPolygons.h"
#include "CAT_Surf.h"
#include "CAT_Intersect.h"

/* #define DEBUG 1 */

void
deform_polygons(polygons_struct *polygons, deform_struct *deform_parms)
{
        int                iter, countdown, countdown2;
        double               avg_err, prev_avg_err, rate;

        iter = 0;
        prev_avg_err = 1e10;
        countdown = 0;
        countdown2 = 0;
    
#ifdef DEBUG
        printf("\n");
#endif
        do {
                iter++;

                avg_err = one_iter_polygons(polygons, deform_parms, iter);
                rate = (prev_avg_err - avg_err) / (prev_avg_err + avg_err);

                if (rate < deform_parms->stop_threshold) {
                        countdown = countdown + 1;
                } else countdown = 0;

                if (prev_avg_err < avg_err) {
                        countdown2 = countdown2 + 1;
                } else countdown2 = 0;

                prev_avg_err = avg_err;
        } while ((countdown < 4 || countdown2 < 4) && countdown < 10 &&
                 iter < deform_parms->max_iterations);
}

void
deform_polygons_one_iter(polygons_struct *polygons,
                         deform_struct *deform_parms, int iter)
{
        one_iter_polygons(polygons, deform_parms, iter);
}

Real
one_iter_polygons(polygons_struct *polygons, deform_struct *deform_parms,
                  int iter)
{
        int                i;
        Point              *new_pts, *tmp;
        deform_stats       stats;
    

        if (polygons->n_points <= 0)
                return(0.0);

        if (!check_correct_deformation_polygons(polygons,
                                              &deform_parms->deformation_model))
                return(0.0);

        if (deform_parms->n_movements_alloced != polygons->n_points) {
                if (deform_parms->n_movements_alloced > 0)
                        FREE(deform_parms->prev_movements);

                deform_parms->n_movements_alloced = polygons->n_points;

                ALLOC(deform_parms->prev_movements, polygons->n_points);

                for (i = 0; i < polygons->n_points; i++)
                        deform_parms->prev_movements[i] =
                               (double) deform_parms->movement_threshold + 1.0f;
        }

        ALLOC(new_pts, polygons->n_points);

        check_polygons_neighbours_computed(polygons);

        initialize_deform_stats(&stats);

        /* every 50 iterations do all points */
        if (iter % 50 == 0) {
                for (i = 0; i < polygons->n_points; i++) {
                        deform_parms->prev_movements[i] =
                               (double) deform_parms->movement_threshold + 1.0f;
                }
        }

        perturb_points(polygons, new_pts,
                       deform_parms->fractional_step,
                       deform_parms->max_step,
                       deform_parms->max_search_distance,
                       deform_parms->degrees_continuity,
                       &deform_parms->deform_data,
                       &deform_parms->boundary_definition,
                       &deform_parms->deformation_model,
                       deform_parms->movement_threshold,
                       deform_parms->prev_movements,
                       &stats);

#ifdef DEBUG
        printf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b");
        printf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b");
        printf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b");
        printf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b");
        printf("Iteration %d:\t",iter);
        printf("avg: %3.4f\tmax: %3.4f",
               stats.average/polygons->n_points, stats.maximum);

        /* print_deform_stats(&stats, polygons->n_points); */
        if (iter % 20 == 0)
                flush_file(stdout);
#endif

        check_polygons_shape_integrity(polygons, new_pts);

        tmp = polygons->points;
        polygons->points = new_pts;
        new_pts = tmp;

        FREE(new_pts);

        if (deform_parms->movement_threshold < 0.0) {
                deform_parms->movement_threshold = stats.average / 10.0;
                if (deform_parms->movement_threshold > 0.1)
                        deform_parms->movement_threshold = 0.1;
        }

        return(stats.average);
}

void
get_polygon_equilibrium_point(polygons_struct *polygons, int poly,
                              int vertidx, double curv_factors[],
                              double max_search_dist, int degrees_continuity,
                              deform_data_struct *deform_data,
                              boundary_definition_struct *bound_def,
                              deformation_model_struct *deform_model,
                              Point *equil_pt)
{
        int              ptidx;
        int              n_nb, neighbours[MAX_NEIGHBOURS];
        BOOLEAN          interior_flag, found;
        double             curv_factor, model_dist, bound_dist;
        double             base_length;
        Point            centroid, model_point, search_origin;
        Vector           normal, pos_model_dir, neg_model_dir;

        ptidx = polygons->indices[
                          POINT_INDEX(polygons->end_indices, poly, vertidx)];

        n_nb = get_subsampled_neighbours_of_point(deform_model, polygons, poly,
                                                  vertidx, neighbours,
                                                  MAX_NEIGHBOURS,
                                                  &interior_flag);

        compute_points_centroid_and_normal(polygons, ptidx, n_nb, neighbours,
                                           &centroid, &normal,
                                           &base_length, &curv_factor);

        get_model_point(deform_model, polygons->points, ptidx, n_nb,
                        neighbours, curv_factors,
                        &centroid, &normal, base_length, &model_point);

        compute_model_dirs(&centroid, &normal, base_length, &model_point,
                           &model_dist, &search_origin,
                           &pos_model_dir, &neg_model_dir);

        found = find_boundary_in_direction(deform_data->volume,
                                           deform_data->label_volume,
                                           NULL, NULL, NULL,
                                           model_dist, &search_origin,
                                           &pos_model_dir, &neg_model_dir,
                                           max_search_dist, max_search_dist,
                                           degrees_continuity,
                                           bound_def, &bound_dist);

        compute_equilibrium_point(ptidx, found, bound_dist, base_length,
                                  model_dist, &pos_model_dir, &neg_model_dir,
                                  &centroid, deform_model, equil_pt);
}

BOOLEAN
ccw_neighbours(Point *centroid, Vector *normal, Point points[],
               int n_nb, int neighbours[], signed char point_error[])
{
        Vector    to_nb, prev_to_nb, up, offset;
        double      len;
        int       i;
        BOOLEAN   ccw;

        ccw = TRUE;
        fill_Vector(to_nb, 0.0, 0.0, 0.0);

        for (i = 0; i < n_nb + 1; i++) {
                prev_to_nb = to_nb;
                SUB_VECTORS(to_nb, points[neighbours[i % n_nb]], *centroid);
                len = DOT_VECTORS(to_nb, *normal);
                SCALE_VECTOR(offset, *normal, len);
                SUB_VECTORS(to_nb, to_nb, offset);

                if (i != 0) {
                        CROSS_VECTORS(up, prev_to_nb, to_nb);
                        if (DOT_VECTORS(up, *normal) < 0.0) {
                                ++point_error[neighbours[i % n_nb]];
                                ++point_error[neighbours[(i-1 + n_nb) % n_nb]];
                                ccw = FALSE;
                        }
                }
        }

        return(ccw);
}

void
check_polygons_shape_integrity(polygons_struct *polygons, Point new_points[])
{
        signed char      *point_done;
        int              vertidx, ptidx, poly, size;
        Point            *centroids;
        Vector           normal;
        progress_struct  progress;
        double           base_length, curv_factor;
        int              n_nb, neighbours[MAX_NEIGHBOURS];
        BOOLEAN          interior_flag;
        signed char      *point_error;
#ifdef  DEBUG
        int              n_errors, n_bad_points;
#endif

        ALLOC(point_done, polygons->n_points);
        ALLOC(point_error, polygons->n_points);
        ALLOC(centroids, polygons->n_points);

        for (ptidx = 0; ptidx <  polygons->n_points; ptidx++) {
                point_done[ptidx] = FALSE;
                point_error[ptidx] = 0;
        }

        initialize_progress_report(&progress, TRUE, polygons->n_items,
                                   "Checking Integrity");

        for (poly = 0; poly < polygons->n_items; poly++) {
                size = GET_OBJECT_SIZE(*polygons, poly);

                for (vertidx = 0; vertidx <  size; vertidx++) {
                        ptidx = polygons->indices[
                          POINT_INDEX(polygons->end_indices, poly, vertidx)];

                        if (!point_done[ptidx]) {
                                point_done[ptidx] = TRUE;

                                compute_polygon_point_centroid(polygons, poly,
                                                              vertidx, ptidx,
                                                              &centroids[ptidx],
                                                              &normal,
                                                              &base_length,
                                                              &curv_factor);
                                n_nb = get_neighbours_of_point(polygons, poly,
                                                               vertidx,
                                                               neighbours,
                                                               MAX_NEIGHBOURS,
                                                               &interior_flag);

#ifdef CHECK_CLOCKWISE_NEIGHBOURS
                                if (!ccw_neighbours(&centroids[ptidx], &normal,
                                                    new_points, n_nb,
                                                    neighbours, point_error)) {
                                        point_error[ptidx] = TRUE;
                                        ++n_errors;
                                }
#else
                                ccw_neighbours(&centroids[ptidx], &normal,
                                               new_points, n_nb,
                                               neighbours, point_error);
#endif
                        }
                }

                update_progress_report(&progress, poly+1);
        }

        terminate_progress_report(&progress);

#ifdef DEBUG
        n_errors = 0;
        n_bad_points = 0;
#endif
        for (ptidx = 0; ptidx < polygons->n_points; ptidx++) {
                if (point_error[ptidx] > 0) {
#ifdef DEBUG
                        ++n_errors;
                        n_bad_points += point_error[ptidx];
                        if (n_errors < 10)
                                printf(" %d", ptidx);
#endif
                        new_points[ptidx] = centroids[ptidx];
                }
        }

#ifdef DEBUG
        if (n_errors > 0)
                printf(": Shape errors %d/%d\n", n_errors, n_bad_points);
#endif

        FREE(point_error);
        FREE(centroids);
        FREE(point_done);
}

void
perturb_points(polygons_struct *polygons, Point new_points[],
               double fractional_step, double max_step, double max_search_dist,
               int degrees_continuity, deform_data_struct *deform_data,
               boundary_definition_struct *bound_def,
               deformation_model_struct *deform_model,
               double movement_thresh, float prev_movements[],
               deform_stats *stats)
{
        double             *curv_factors;
        signed char      *point_done;
        int              vertidx, ptidx, poly, size, n1, n2;
        Point            centroid;
        Vector           normal;
        progress_struct  progress;
        Point            equil_pt;
        double             dist_to_equil, base_length;
        double           *movements;

        ALLOC(curv_factors, polygons->n_points);
        ALLOC(point_done, polygons->n_points);

        if (deformation_model_includes_average(deform_model)) {
                for (ptidx = 0; ptidx < polygons->n_points; ptidx++)
                        point_done[ptidx] = FALSE;

                initialize_progress_report(&progress, TRUE, polygons->n_items,
                                           "Computing Curvatures");

                for (poly = 0; poly <  polygons->n_items; poly++) {
                        size = GET_OBJECT_SIZE(*polygons, poly);

                        for (vertidx = 0; vertidx < size; vertidx++) {
                                ptidx = polygons->indices[POINT_INDEX(
                                                          polygons->end_indices,
                                                          poly, vertidx)];

                                if (!point_done[ptidx]) {
                                        point_done[ptidx] = TRUE;

                                        compute_polygon_point_centroid(polygons,
                                                         poly, vertidx, ptidx,
                                                         &centroid, &normal,
                                                         &base_length,
                                                         &curv_factors[ptidx]);
                                }
                        }

                        update_progress_report(&progress, poly+1);
                }

                terminate_progress_report(&progress);
        }

        ALLOC(movements, polygons->n_points);

        for (ptidx = 0; ptidx < polygons->n_points; ptidx++) {
                movements[ptidx] = 0.0f;
                point_done[ptidx] = FALSE;
                new_points[ptidx] = polygons->points[ptidx];
        }

        for (poly = 0; poly <  polygons->n_items; poly++) {
                size = GET_OBJECT_SIZE(*polygons, poly);

                for (vertidx = 0; vertidx < size; vertidx++) {
                        ptidx = polygons->indices[POINT_INDEX(
                                                 polygons->end_indices,
                                                 poly, vertidx)];

                        if (point_done[ptidx])
                                continue;

                        n1 = polygons->indices[POINT_INDEX(
                                               polygons->end_indices,
                                               poly, (vertidx+1) % size)];
                        n2 = polygons->indices[POINT_INDEX(
                                               polygons->end_indices,
                                               poly, (vertidx-1+size)%size)];

                        if ((Real) prev_movements[ptidx] < movement_thresh &&
                            (Real) prev_movements[n1] < movement_thresh &&
                            (Real) prev_movements[n2] < movement_thresh)
                                continue;

                        point_done[ptidx] = TRUE;

                        get_polygon_equilibrium_point(polygons, poly,
                                                      vertidx, curv_factors,
                                                      max_search_dist,
                                                      degrees_continuity,
                                                      deform_data, bound_def,
                                                      deform_model, &equil_pt);

                        dist_to_equil = deform_point(ptidx, polygons->points,
                                             &equil_pt,
                                             fractional_step, max_step,
                                             deform_model->position_constrained,
                                             deform_model->max_position_offset,
                                             deform_model->original_positions,
                                             &new_points[ptidx]);

                        movements[ptidx] = (double) distance_between_points(
                                                       &polygons->points[ptidx],
                                                       &new_points[ptidx]);

                        record_error_in_deform_stats(stats, dist_to_equil);
                }
        }

        for (ptidx = 0; ptidx <  polygons->n_points; ptidx++) {
                if (!point_done[ptidx])
                        record_error_in_deform_stats(stats, 0.0);
                prev_movements[ptidx] = movements[ptidx];
        }

        FREE(movements);
        FREE(point_done);
        FREE(curv_factors);
}

/* deform polygons and check for self intersections every check_every_iteration and optionally 
   force that no self intersections exist anymore */
void
deform_polygons_check_selfintersection(polygons_struct *polygons, deform_struct *deform_parms,
                       int check_every_iteration, int force_no_selfintersections)
{
        int                iter, countdown, countdown2, p, n_intersects, counter;
        int                *defects, *polydefects;
        int                *n_neighbours, **neighbours;
        double             avg_err, prev_avg_err, rate;
        Point              *pts;

        pts         = (Point *) malloc(sizeof(Point) * polygons->n_points);
        defects     = (int *) malloc(sizeof(int) * polygons->n_points);
        polydefects = (int *) malloc(sizeof(int) * polygons->n_items);

        create_polygon_point_neighbours(polygons, TRUE, &n_neighbours,
                                        &neighbours, NULL, NULL);
                
        for (p = 0; p < polygons->n_points; p++)
                pts[p]  = polygons->points[p];
        
        iter = 0;
        prev_avg_err = 1e10;
        countdown = 0;
        countdown2 = 0;
        counter = 0;
    
#ifdef DEBUG
        printf("\n");
#endif
        do {
                iter++;

                avg_err = one_iter_polygons(polygons, deform_parms, iter);

                if ((check_every_iteration) && (((iter % check_every_iteration) == 0))) {
                
                        counter = 0;
                        n_intersects = find_selfintersections(polygons, defects, polydefects);
                        n_intersects = join_intersections(polygons, defects, polydefects,
                                              n_neighbours, neighbours);
                                              
                        do {
                                counter++;
                                
                                if (n_intersects > 0) {
                                        printf("%3d self intersections found that will be corrected.\n", n_intersects);
                
                                        n_intersects = smooth_selfintersections(polygons, defects, polydefects,
                                                     n_intersects, n_neighbours,
                                                     neighbours, 50);
        
                                }
                        } while (n_intersects > 0 && counter < 10);

                }
                
                rate = (prev_avg_err - avg_err) / (prev_avg_err + avg_err);

                if (rate < deform_parms->stop_threshold) {
                        countdown = countdown + 1;
                } else countdown = 0;

                if (prev_avg_err < avg_err) {
                        countdown2 = countdown2 + 1;
                } else countdown2 = 0;

                prev_avg_err = avg_err;

        } while ((countdown < 4 || countdown2 < 4) && countdown < 10 &&
                 iter < deform_parms->max_iterations);
                 
        if (iter < deform_parms->max_iterations)
                printf("Stopped after %d iterations\n", iter);

        if (force_no_selfintersections) {
                do {
                        counter++;
                        n_intersects = find_selfintersections(polygons, defects, polydefects);
                        n_intersects = join_intersections(polygons, defects, polydefects,
                                                  n_neighbours, neighbours);
                        
                        if (n_intersects > 0) {
                                printf("%3d self intersections found that will be corrected.\n", n_intersects);
        
                                n_intersects = smooth_selfintersections(polygons, defects, polydefects,
                                             n_intersects, n_neighbours,
                                             neighbours, 50);

                        }
                } while (n_intersects > 0 && counter < 10);
        }
        printf("\n\n");
        free(pts);

}

/* deform polygons and check for self intersections every check_every_iteration and optionally 
   force that no self intersections exist anymore */
void
deform_polygons_check_selfintersection_old(polygons_struct *polygons, deform_struct *deform_parms,
                       int check_every_iteration, int force_no_selfintersections)
{
        int                iter, countdown, countdown2, p, n_intersects, counter;
        int                *defects, *polydefects;
        int                *n_neighbours, **neighbours;
        double             avg_err, prev_avg_err, rate;
        Point              *pts;

        pts         = (Point *) malloc(sizeof(Point) * polygons->n_points);
        defects     = (int *) malloc(sizeof(int) * polygons->n_points);
        polydefects = (int *) malloc(sizeof(int) * polygons->n_items);

        create_polygon_point_neighbours(polygons, TRUE, &n_neighbours,
                                        &neighbours, NULL, NULL);
                
        for (p = 0; p < polygons->n_points; p++)
                pts[p]  = polygons->points[p];
        
        iter = 0;
        prev_avg_err = 1e10;
        countdown = 0;
        countdown2 = 0;
        counter = 0;
    
#ifdef DEBUG
        printf("\n");
#endif
        do {
                iter++;

                avg_err = one_iter_polygons(polygons, deform_parms, iter);

                if ((check_every_iteration > 0) && (((iter % check_every_iteration) == 0) || (iter == 1))) {
                
                        n_intersects = find_selfintersections(polygons, defects, polydefects);
                        n_intersects = join_intersections(polygons, defects, polydefects,
                                              n_neighbours, neighbours);
                                              
                        if (n_intersects > 0) {
                                printf("%d self intersections found at iteration %d\n", n_intersects, iter);
                                for (p = 0; p < polygons->n_points; p++)
                                        if (defects[p] > 0)
                                                polygons->points[p] = pts[p]; /* restore */
                        }
                        
                        /* save points for next check */
                        for (p = 0; p < polygons->n_points; p++)
                                pts[p] = polygons->points[p];
                }
                
                rate = (prev_avg_err - avg_err) / (prev_avg_err + avg_err);

                if (rate < deform_parms->stop_threshold) {
                        countdown = countdown + 1;
                } else countdown = 0;

                if (prev_avg_err < avg_err) {
                        countdown2 = countdown2 + 1;
                } else countdown2 = 0;

                prev_avg_err = avg_err;

        } while ((countdown < 4 || countdown2 < 4) && countdown < 10 &&
                 iter < deform_parms->max_iterations);
                 
        if (iter < deform_parms->max_iterations)
                printf("Stopped after %d iterations\n", iter);

        if (force_no_selfintersections) {
                do {
                        counter++;
                        n_intersects = find_selfintersections(polygons, defects, polydefects);
                        n_intersects = join_intersections(polygons, defects, polydefects,
                                                  n_neighbours, neighbours);
                        
                        if (n_intersects > 0) {
                                printf("%d self intersections found that will be corrected.\n", n_intersects);
        
                                n_intersects = smooth_selfintersections(polygons, defects, polydefects,
                                             n_intersects, n_neighbours,
                                             neighbours, 50);

                        }
                } while (n_intersects > 0 && counter < 10);
        }
        printf("\n\n");
        free(pts);

}
