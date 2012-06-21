/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 *
 */

#include <bicpl.h>

void
get_all_polygon_point_neighbours(polygons_struct *polygons,
                                 int *n_point_neighbours_ptr[],
                                 int **point_neighbours_ptr[] )
{
        int         p, poly, size, vertex;
        int         *n_neighbours, **neighbours;
        int         *total_neighbours, total_n_neighbours;
        BOOLEAN     interior;

        check_polygons_neighbours_computed(polygons);

        ALLOC(n_neighbours, polygons->n_points);
        ALLOC(neighbours, polygons->n_points);
        for (p = 0; p < polygons->n_points; p++)
                n_neighbours[p] = 0;

        total_n_neighbours = 0;
        for (poly = 0; poly < polygons->n_items; poly++) {
                size = GET_OBJECT_SIZE(*polygons, poly);
                for (vertex = 0; vertex < size; vertex++) {
                        p = polygons->indices[
                              POINT_INDEX(polygons->end_indices,poly,vertex)];
                        if (n_neighbours[p] > 0)
                                continue;

                        n_neighbours[p] = get_neighbours_of_point(polygons,
                                                    poly, vertex, NULL, 0,
                                                    &interior);
                        total_n_neighbours += n_neighbours[p];
                }
        }

        ALLOC(total_neighbours, total_n_neighbours);
        total_n_neighbours = 0;
        for (p = 0; p < polygons->n_points; p++) {
                neighbours[p] = &total_neighbours[total_n_neighbours];
                total_n_neighbours += n_neighbours[p];
        }

        for (p = 0; p < polygons->n_points; p++)
                n_neighbours[p] = 0;

        for (poly = 0; poly < polygons->n_items; poly++) {
                size = GET_OBJECT_SIZE(*polygons, poly);
                for (vertex = 0; vertex < size; vertex++) {
                        p = polygons->indices[
                              POINT_INDEX(polygons->end_indices,poly,vertex)];

                        if (n_neighbours[p] > 0)
                                continue;

                        n_neighbours[p] = get_neighbours_of_point(polygons,
                                                    poly, vertex, neighbours[p],
                                                    total_n_neighbours,
                                                    &interior);
                }
        }
    
        *n_point_neighbours_ptr = n_neighbours;
        *point_neighbours_ptr = neighbours;
}

Real
evaluate_heatkernel(double x, double sigma)
{
    return(exp(-x / (2 * sigma * sigma)));
}

void
heatkernel_blur_points(int n_polygon_pts, Point polygon_pts[],
                       double values[], int n_neighbours, int *neighbours,
                       int ptidx, double sigma, Point *smooth_point,
                       double *value)
{
        double   sum[3], weight, sum_weight, point_dist;
        int    i, c, n_pts, neigh;

        if (sigma <= 0.0)
                sigma = 1e-20;

        sum[0] = 0.0;
        sum[1] = 0.0;
        sum[2] = 0.0;
        sum_weight = 0.0;
    
        for (i = 0; i < n_neighbours+1; i++)
    	        *value = 1.0;

        for (i = 0; i < n_neighbours+1; i++) {
                if (i > 0) {
                        neigh = neighbours[i-1];
                        point_dist = distance_between_points(
                                        &polygon_pts[ptidx],
                                        &polygon_pts[neigh]);
                } else {
                        point_dist = 0.0;
                        neigh = ptidx;
                }
                weight = evaluate_heatkernel(point_dist, sigma);

                if (values != NULL) {
                        sum[0] += weight * values[neigh];
                } else {
                        for (c = 0; c < N_DIMENSIONS; c++)
                                sum[c] += weight * (Real)
                                          Point_coord(polygon_pts[neigh], c);
                }

                sum_weight += weight;
        }

        if (values != NULL) {
                *value = sum[0] / sum_weight;
        } else {
                for (c = 0; c < N_DIMENSIONS; c++)
                        Point_coord(*smooth_point, c) = (Point_coord_type)
                                                        (sum[c] / sum_weight);
        }
}

void
smooth_heatkernel(polygons_struct *polygons, double *values, double fwhm)
{
        double             sigma, value, *smooth_values;
        Point            point, *smooth_pts;
        int              n_iter, i, j;
        int              *n_neighbours, **neighbours;
        BOOLEAN          values_present;
        progress_struct  progress;
        
        get_all_polygon_point_neighbours(polygons, &n_neighbours, &neighbours);
           
        /* calculate n_iter with regard to fwhm */    
        sigma = 1.0;
        n_iter = ceil(fwhm*fwhm*3/(8*log(2)));

        initialize_progress_report(&progress, FALSE, n_iter*polygons->n_points,
                                   "Blurring");
        if (values != NULL) 
                smooth_values = (double *) malloc(sizeof(double) * polygons->n_points);
        else    smooth_pts = (Point *) malloc(sizeof(Point) * polygons->n_points);
        
        
        /* diffusion smoothing using heat kernel */
        for (j = 0; j < n_iter; j++) {
                for (i = 0; i < polygons->n_points; i++) {
                        heatkernel_blur_points(polygons->n_points,
                                               polygons->points, values,
                                               n_neighbours[i], neighbours[i],
                                               i, sigma, &point, &value);

                        if (values != NULL)
                                smooth_values[i] = value;
                        else    smooth_pts[i] = point;

                        update_progress_report(&progress,
                                               j*polygons->n_points + i + 1);
                }
                for (i = 0; i < polygons->n_points; i++) {
                        if (values != NULL)
                                values[i] = smooth_values[i];
                        else
                                polygons->points[i] = smooth_pts[i];
                }
        }

        terminate_progress_report(&progress);

        if (values != NULL) 
                free(smooth_values);
        else    free(smooth_pts);
        
}
