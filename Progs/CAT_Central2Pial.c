/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id: CAT_Central2Pial.c 325 2014-10-30 15:39:48Z gaser $
 *
 */

#include <bicpl.h>

#include "CAT_Surf.h"
#include "CAT_SurfaceIO.h"
#include "CAT_Smooth.h"

void
usage(char *executable)
{
        static char *usage_str = "\n\
Usage: %s  surface_file values_file output_surface_file\n\
Estimate pial surface from central surface using cortical thickness.\n\n\n";

       fprintf(stderr, usage_str, executable);
}

int
main(int argc, char *argv[])
{
        double               value, *values;
        Status               status;
        char                 *src_file, *out_file, *values_file;
        int                  i, p, n_objects, n_pts, n_values;
        Point                *pts;
        File_formats         format;
        polygons_struct      *polygons;
        object_struct        **object_list;

        initialize_argument_processing(argc, argv);

        if (!get_string_argument( NULL, &src_file) ||
            !get_string_argument( NULL, &values_file) ||
            !get_string_argument( NULL, &out_file)) {
                usage(argv[0]);
                exit(EXIT_FAILURE);
        }
        
        if (input_graphics_any_format(src_file, &format, &n_objects,
                                      &object_list) != OK)
                exit(EXIT_FAILURE);

        if (input_values_any_format(values_file, &n_values, &values) != OK)
                exit(EXIT_FAILURE);
      
        if (n_objects > 1) {
                fprintf(stderr,"Only one object allowed.\n");
                exit(EXIT_FAILURE);
        }

        n_pts = get_object_points(object_list[0], &pts);
        if (n_pts != n_values) {
                fprintf(stderr,"Number of points differs from number of values.\n");
                exit(EXIT_FAILURE);
        }

        polygons = get_polygons_ptr(object_list[0]);
        compute_polygon_normals(polygons);

        for (p = 0; p < polygons->n_points; p++) {
                Point_x(polygons->points[p]) += 0.5*values[p]*Point_x(polygons->normals[p]);
                Point_y(polygons->points[p]) += 0.5*values[p]*Point_y(polygons->normals[p]);
                Point_z(polygons->points[p]) += 0.5*values[p]*Point_z(polygons->normals[p]);
        }
        
        /* smooth surface slightly */
        smooth_heatkernel(polygons, NULL, 1);

        compute_polygon_normals(polygons);

        if(output_graphics_any_format(out_file, format, 1, object_list, NULL) != OK)
                exit(EXIT_FAILURE);

        return(status != OK);
}
