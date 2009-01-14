/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 *
 */

/*
 * Heat kernel smoothing is based on matlab code from Moo K. Chung:
 * Chung, M.K., Robbins,S., Dalton, K.M., Davidson, R.J., Evans, A.C. (2005) 
 * Cortical thickness analysis in autism via heat kernel smoothing. NeuroImage. 
 * http://www.stat.wisc.edu/~mchung/papers/ni_heatkernel.pdf
 */

#include <volume_io/internal_volume_io.h>
#include <bicpl.h>

#include "CAT_Blur2d.h"
#include "CAT_SurfaceIO.h"
#include "CAT_Curvature.h"

void
usage(char *executable)
{
        char *usage_str = "\n\
Usage: %s object_file output_file fwhm [gausscurv_threshold]\n\n\
     Diffusion smoothing of surface points w.r.t. gaussian curvature using\n\
     heat kernel.\n\n";

        fprintf(stderr, usage_str, executable);
}

int
main(int argc, char *argv[])
{
        char             *input_file, *output_file;
        int              n_objects, i, j, n_iter;
        int              *n_neighbours, **neighbours;
        File_formats     format;
        object_struct    **object_list;
        polygons_struct  *polygons;
        Point            *smooth_pts, point;
        Real             fwhm, value;
        progress_struct  progress;
        double           *gc_strength, gc_threshold, sigma;

        initialize_argument_processing(argc, argv);

        if (!get_string_argument(NULL, &input_file) ||
            !get_string_argument(NULL, &output_file)) {
                usage(argv[0]);
                return(1);
        }

        get_real_argument(25.0, &fwhm);
        get_real_argument(0.01, &gc_threshold);

        if (input_graphics_any_format(input_file, &format, &n_objects,
                                      &object_list) != OK ||
            n_objects != 1 || get_object_type(object_list[0]) != POLYGONS) {
                fprintf(stderr, "Error reading %s.\n", input_file);
                return(1);
        }

        polygons = get_polygons_ptr(object_list[0]);

        gc_strength = (double *)malloc(sizeof(double)*polygons->n_points);
        smooth_pts  = (Point *)malloc(sizeof(Point)*polygons->n_points);
            
        get_all_polygon_point_neighbours(polygons, &n_neighbours, &neighbours);

        get_polygon_vertex_curvatures_cg(polygons, n_neighbours, neighbours,
                                         0.0, 1, gc_strength);

        sigma = 5.0;
        /* use absolute gaussian gc_strength if values are above threshold otherwise don't smooth */
        for (i=0; i<polygons->n_points; i++) 
                gc_strength[i] = (fabs(gc_strength[i]) > gc_threshold) ? sigma*fabs(gc_strength[i]) : 0.0;       

        /* calculate n_iter in relation to sigma */
        n_iter = ROUND(fwhm/2.35482 * fwhm/2.35482/sigma);
        if (n_iter == 0)
                n_iter = 1;

        initialize_progress_report(&progress, FALSE, n_iter*polygons->n_points,
                                   "Blurring");

        /* diffusion smoothing using heat kernel */
        for (j = 0; j < n_iter; j++) {
                for (i = 0; i < polygons->n_points; i++) {
                        /* smooth only if strength is > 0 */
                        if (gc_strength[i] > 0.0) {
                                
                                heatkernel_blur_points(polygons->n_points,
                                               polygons->points, NULL,
                                               n_neighbours[i], neighbours[i],
                                               i, gc_strength[i], &point, &value);
                                smooth_pts[i] = point;
                        } else smooth_pts[i] = polygons->points[i];

                        update_progress_report(&progress, j*polygons->n_points + i + 1);
                }
                for (i = 0; i < polygons->n_points; i++)
                        polygons->points[i] = smooth_pts[i];
        }

        terminate_progress_report(&progress);

        polygons->points = smooth_pts;

        compute_polygon_normals(polygons);

        output_graphics_any_format(output_file, format, 1, object_list);
        
        free(smooth_pts);
        free(gc_strength);

        return(0);
}
