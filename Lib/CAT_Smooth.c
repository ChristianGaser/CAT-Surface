/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 *
 */

#include <bicpl.h>

#ifndef isnan
#define isnan(a) ((a)!=(a)) 
#endif

#ifdef _MSC_VER
  static const unsigned long __nan[2] = {0xffffffff, 0x7fffffff};
  #define FNAN (*(const float *) __nan)
#else
  #define FNAN 0.0f/0.0f
#endif

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

double
evaluate_heatkernel(double x, double sigma)
{
    return(exp(-(x * x) / (2.0 * sigma * sigma)));
}

void
heatkernel_blur_points(int n_polygon_pts, Point polygon_pts[],
                       double values[], int n_neighbours, int *neighbours,
                       int ptidx, double sigma, Point *smooth_point,
                       double *value)
{
        double   sum[3], weight, sum_weight, point_dist;
        int      i, c, n_pts, neigh;

        if (sigma <= 0.0)
                sigma = 1e-20;

        sum[0] = sum[1] = sum[2] = sum_weight = 0.0;
    
        for (i = 0; i < n_neighbours+1; i++)
                *value = 1.0;

        for (i = 0; i < n_neighbours+1; i++) {
                if (i > 0) {    /* neighbouring points */
                        neigh = neighbours[i-1];
                        point_dist = distance_between_points(
                                        &polygon_pts[ptidx],
                                        &polygon_pts[neigh]);
                } else {        /* center point */
                        point_dist = 0.0;
                        neigh = ptidx;
                }
                
                if (values != NULL) {
                        /* use Gaussian kernel for values */
                        weight = evaluate_heatkernel(point_dist, sigma);
                        
                        /* only consider values that are not NaN */
                        if (!isnan(values[neigh])) {
                                sum[0] += weight * values[neigh];
                                sum_weight += weight;
                        }
                } else {
                        /* this is rather based on empirically estimated values to be 
                           compatible to older versions */
                        if (i > 0) weight = 0.25;
                        else weight = 1.0;
                        
                        for (c = 0; c < N_DIMENSIONS; c++)
                                sum[c] += weight * (double)
                                          Point_coord(polygon_pts[neigh], c);
                        sum_weight += weight;
                }
        }
        
        if (values != NULL)
                *value = sum[0] / sum_weight;
        else {
                for (c = 0; c < N_DIMENSIONS; c++)
                        Point_coord(*smooth_point, c) = (Point_coord_type)(sum[c] / sum_weight);
        }
}

void
smooth_heatkernel(polygons_struct *polygons, double *values, double fwhm)
{
        double           sigma, value, point_dist, sum_dist, *smooth_values;
        Point            point, *smooth_pts;
        int              n_iter, i, j, size, p1, p2, n;
        int              *n_neighbours, **neighbours;
        BOOLEAN          values_present;
        progress_struct  progress;
        
        get_all_polygon_point_neighbours(polygons, &n_neighbours, &neighbours);
           
        /* estimate mean distance between the points to define optimal sigma */
        n = 0;
        sum_dist = 0.0;
        for (i = 0; i < polygons->n_items; i++) {
                size = GET_OBJECT_SIZE( *polygons, i );
    
                for (j = 0; j < size; j++) {
                        p1 = polygons->indices[POINT_INDEX(polygons->end_indices,i,j)];
                        p2 = polygons->indices[POINT_INDEX(polygons->end_indices,i,(j+1)%size)];

                        point_dist = distance_between_points(
                                                &polygons->points[p1],
                                                &polygons->points[p2]);
                        n++;
                        sum_dist += point_dist;
                }
        }
        
        /* use mean distance as sigma */
        sigma = sum_dist/(double)n;
        
        /* calculate n_iter with regard to fwhm   
           0.541011 equals to 3.0/(8.0*log(2.0))
           see SurfStatSmooth.m in surfstat from Keith Worsley */
        n_iter = ceil(fwhm*fwhm*0.541011/(sigma*sigma));
        if (n_iter < 1) n_iter = 1;

        /* recalibrate sigma to ensure integer numbers for iterations */
        sigma = (fwhm*0.7355345)/sqrt((double)n_iter);

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

                        if (values != NULL) {
                                if(!isnan(values[i]))
                                        smooth_values[i] = value;
                                else    smooth_values[i] = FNAN;
                        } else    smooth_pts[i] = point;

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
