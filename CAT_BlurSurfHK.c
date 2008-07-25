/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
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

void
usage(char *executable)
{
        char *usage_str = "\n\
Usage: %s object_file output_file fwhm [values_file]\n\n\
     Diffusion smoothing of values or surface points using\n\
     heat kernel. If values are defined then values will be\n\
     smoothed, otherwise only surface points.\n\n";

        fprintf(stderr, usage_str, executable);
}

int
main(int argc, char *argv[])
{
        char             *input_file, *output_file, *values_file;
        int              n_objects, i, j, n_iter;
        int              *n_neighbours, **neighbours;
        FILE             *fp;
        File_formats     format;
        object_struct    **object_list;
        polygons_struct  *polygons;
        Point            *smooth_pts, point;
        Real             fwhm, *values, value, *smooth_values, sigma;
        BOOLEAN          values_present;
        progress_struct  progress;

        initialize_argument_processing(argc, argv);

        if (!get_string_argument(NULL, &input_file) ||
            !get_string_argument(NULL, &output_file) ||
            !get_real_argument(0.0, &fwhm)) {
                usage(argv[0]);
                return(1);
        }

        values_present = get_string_argument(NULL, &values_file);

        if (input_graphics_any_format(input_file, &format, &n_objects,
                                      &object_list) != OK ||
            n_objects != 1 || get_object_type(object_list[0]) != POLYGONS) {
                fprintf(stderr, "Error reading %s.\n", input_file);
                return(1);
        }

        polygons = get_polygons_ptr(object_list[0]);

        if (values_present) {
                ALLOC(values, polygons->n_points);
                ALLOC(smooth_values, polygons->n_points);

                if (open_file(values_file, READ_FILE, ASCII_FORMAT, &fp) != OK)
                        return(1);

                for (i = 0; i < polygons->n_points; i++) {
                        if (input_real(fp, &values[i]) != OK)
                                return(1);
                }

                close_file(fp);

                smooth_pts = NULL;
        } else {
                ALLOC(smooth_pts, polygons->n_points);
                values = NULL;
        }

        get_all_polygon_point_neighbours(polygons, &n_neighbours, &neighbours);

        /* select sigma according fwhm */
        if (fwhm > 50.0)
                sigma = 8.0;
        else if (fwhm > 30.0)
                sigma = 3.0;
        else if (fwhm > 20.0)
                sigma = 2.0;
        else sigma = 1.0;

        /* calculate n_iter in relation to sigma */
        n_iter = ROUND(fwhm/2.35482 * fwhm/2.35482/sigma);
        if (n_iter == 0)
                n_iter = 1;
    
        initialize_progress_report(&progress, FALSE, n_iter*polygons->n_points,
                                   "Blurring");

        /* diffusion smoothing using heat kernel */
        for (j = 0; j < n_iter; j++) {
                for (i = 0; i < polygons->n_points; i++) {
                        heatkernel_blur_points(polygons->n_points,
                                               polygons->points, values,
                                               n_neighbours[i], neighbours[i],
                                               i, sigma, &point, &value);
                        if (values_present)
                                smooth_values[i] = value;
                        else
                                smooth_pts[i] = point;

                        update_progress_report(&progress,
                                               j*polygons->n_points + i + 1);
                }
                for (i = 0; i < polygons->n_points; i++) {
                        if (values_present)
                                values[i] = smooth_values[i];
                        else
                                polygons->points[i] = smooth_pts[i];
                }
        }

        terminate_progress_report(&progress);

        if (values_present) {
                if (open_file(output_file, WRITE_FILE, ASCII_FORMAT, &fp) != OK)
                        return(1);

                for (i = 0; i < polygons->n_points; i++) {
                        if (output_real(fp, smooth_values[i]) != OK ||
                            output_newline(fp) != OK)
                                return(1);
                }

                close_file(fp);
                FREE(smooth_values);
                FREE(values);
        } else {
                polygons->points = smooth_pts;

                compute_polygon_normals(polygons);

                output_graphics_any_format(output_file, format, 1, object_list);
                FREE(smooth_pts);
        }

        return(0);
}
