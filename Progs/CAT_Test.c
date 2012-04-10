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

#include <bicpl.h>

#include "CAT_Smooth.h"
#include "CAT_SurfaceIO.h"

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
        char             *input_file, *output_file;
        int              n_objects, i;
        int              *n_neighbours, **neighbours;
        File_formats     format;
        object_struct    **object_list, *out_object;
        polygons_struct  *polygons, *polygonsIn;
        Point            *smooth_pts;
        Real             fwhm;
        double           *sulc_depth, *faces, *vertices;

        initialize_argument_processing(argc, argv);

        if (!get_string_argument(NULL, &input_file) ||
            !get_string_argument(NULL, &output_file)) {
                usage(argv[0]);
                exit(EXIT_FAILURE);
        }

        if (input_graphics_any_format(input_file, &format, &n_objects,
                                      &object_list) != OK ||
            n_objects != 1 || get_object_type(object_list[0]) != POLYGONS) {
                fprintf(stderr, "Error reading %s.\n", input_file);
                exit(EXIT_FAILURE);
        }

        polygons = get_polygons_ptr(object_list[0]);
        out_object = create_object(POLYGONS);
        polygonsIn = get_polygons_ptr(out_object);
        copy_polygons(polygons, polygonsIn);

        sulc_depth = (double *) malloc(sizeof(double) * polygons->n_points);

        get_all_polygon_point_neighbours(polygons, &n_neighbours, &neighbours);

        smooth_heatkernel(polygons, NULL, fwhm);

        for (i = 0; i < polygons->n_points; i++) {
                sulc_depth[i] = sqrt((Point_x(polygons->points[i]) - Point_x(polygonsIn->points[i]))*(Point_x(polygons->points[i]) - Point_x(polygonsIn->points[i])) +
                                (Point_y(polygons->points[i]) - Point_y(polygonsIn->points[i]))*(Point_y(polygons->points[i]) - Point_y(polygonsIn->points[i])) +
                                (Point_z(polygons->points[i]) - Point_z(polygonsIn->points[i]))*(Point_z(polygons->points[i]) - Point_z(polygonsIn->points[i])));
        }

        output_values_any_format(output_file, polygons->n_points,
                                 sulc_depth, TYPE_DOUBLE);


        return(EXIT_SUCCESS);
}
