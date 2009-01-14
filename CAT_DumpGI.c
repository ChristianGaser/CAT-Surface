/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 *
 */

/* Heat kernel smoothing is based on matlab code from Moo K. Chung:
 *  Chung, M.K., Robbins,S., Dalton, K.M., Davidson, R.J., Evans, A.C. (2004) 
 *  Cortical thickness analysis in autism via heat kernel smoothing.
 * NeuroImage, submitted. 
 *  http://www.stat.wisc.edu/~mchung/papers/ni_heatkernel.pdf
 */

#include <volume_io/internal_volume_io.h>
#include <bicpl.h>

#include "CAT_Blur2d.h"
#include "CAT_Curvature.h"

void
usage(char *executable)
{
        char *usage_str = "\n\
Usage: %s  object_file output_file fwhm [fwhm_surf] [curvtype]\n\n\
     Calculate gyrification index using the ratio of local surface area\n\
     and local inflated surface area. Local surface area can be approximated\n\
     by use of different curve types (default: mean curvature averaged over\n\
     3mm):\n\
     curvtype:  0 - mean curvature (averaged over 3mm, in degrees)\n\
                2 - curvedness\n\
                3 - shape index\n\
                4 - mean curvature (in radians)\n\n";

        fprintf(stderr, usage_str, executable);
}

int
main(int argc, char *argv[])
{
        char                 *object_file, *output_file;
        File_formats         format;
        int                  i, j, n_iter, n_objects, curvtype;
        int                  *n_neighbours, **neighbours;
        object_struct        **objects;
        polygons_struct      *polygons;
        Point                *smooth_points, point;
        signed char          *done_flags;
        Real                 fwhm, fwhm_surf, *curvatures, *curvs_inflated;
        Real                 *GI, value, distance, sigma;
        progress_struct      progress;

        initialize_argument_processing(argc, argv);

        if (!get_string_argument(NULL, &object_file) ||
            !get_string_argument(NULL, &output_file)) {
                usage(argv[0]);
                return(1);
        }

        get_real_argument(30.0, &fwhm);
        get_real_argument(100.0, &fwhm_surf);
        get_int_argument(0, &curvtype);
    
        if (input_graphics_any_format(object_file, &format,
                                      &n_objects, &objects) != OK)
                return(1);

        if (n_objects != 1 || get_object_type(objects[0]) != POLYGONS) {
                printf("File must contain 1 polygons object.\n");
                return(1);
        }

        polygons = get_polygons_ptr(objects[0]);

        compute_polygon_normals(polygons);

        ALLOC(curvatures, polygons->n_points);
        ALLOC(curvs_inflated, polygons->n_points);
        ALLOC(GI, polygons->n_points);
        ALLOC(smooth_points, polygons->n_points);
        ALLOC(done_flags, polygons->n_points);

        get_all_polygon_point_neighbours(polygons, &n_neighbours, &neighbours);
    
        // get curvature values
        if (curvtype == 0)
                distance = 3.0;
        else
                distance = 0.0;
        get_polygon_vertex_curvatures_cg(polygons, n_neighbours, neighbours,
                                         distance, curvtype, curvatures);

        /* inflate surface by smoothing with FWHM of 150mm */
    
        /* calculate n_iter with regard to sigma */
        sigma = 8.0;
        n_iter = ROUND(fwhm_surf/2.35482 * fwhm_surf/2.35482/sigma);

        for (i = 0; i < polygons->n_points; i++)
                done_flags[i] = FALSE;

        initialize_progress_report(&progress, FALSE, n_iter*polygons->n_points,
                                   "Blurring surface" );

        /* diffusion smoothing using heat kernel */                        
        for (j = 0; j < n_iter; j++) {
                for (i = 0; i < polygons->n_points; i++) {
                        heatkernel_blur_points(polygons->n_points,
                                               polygons->points, NULL,
                                               n_neighbours[i], neighbours[i],
                                               i, sigma, &point, &value);
                        smooth_points[i] = point;

                        update_progress_report(&progress,
                                               j*polygons->n_points + i + 1);
                }
                for (i = 0; i <  polygons->n_points; i++) {
                        polygons->points[i] = smooth_points[i];
                }
        }
        terminate_progress_report(&progress);

        polygons->points = smooth_points;
        compute_polygon_normals(polygons);

        /* get curvature values of inflated surface */
        get_polygon_vertex_curvatures_cg(polygons, n_neighbours, neighbours,
                                         distance, curvtype, curvs_inflated);

        /* calculate ratio of absolute values */
        for (i = 0; i < polygons->n_points; i++) {
                if (curvs_inflated[i] != 0)
                        GI[i] = fabs(curvatures[i]/curvs_inflated[i]);
                else
                        GI[i] = 0.0;
        }
    
        /* smooth GI values */
        /* calculate n_iter for sigma = 1.0 */
        n_iter = ROUND(fwhm/2.35482 * fwhm/2.35482);
        sigma = 1.0;

        for (i = 0; i < polygons->n_points; i++)
                done_flags[i] = FALSE;

        initialize_progress_report(&progress, FALSE, n_iter*polygons->n_points,
                                   "Blurring values");

        for (j = 0; j < n_iter; j++) {
                for (i = 0; i < polygons->n_points; i++) {
                        heatkernel_blur_points(polygons->n_points,
                                               polygons->points, GI,
                                               n_neighbours[i], neighbours[i],
                                               i, sigma, NULL, &value);
                        GI[i] = value;
                        update_progress_report(&progress,
                                               j*polygons->n_points + i + 1 );
                }
        }

        terminate_progress_report(&progress);

        output_values_any_format(output_file, polygons->n_points, GI);

        delete_object_list(n_objects, objects);
        FREE(curvatures);
        FREE(curvs_inflated);
        FREE(GI);
        FREE(done_flags);

        return(0);
}
