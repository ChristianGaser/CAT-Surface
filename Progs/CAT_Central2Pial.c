/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 *
 */

#include <bicpl.h>

#include "CAT_Surf.h"
#include "CAT_SurfaceIO.h"
#include "CAT_Smooth.h"
#include "CAT_DeformPolygons.h"

void
usage(char *executable)
{
        static char *usage_str = "\n\
Usage: %s  surface_file values_file output_surface_file [-1]\n\
Estimate pial surface from central surface using cortical thickness. If you add -1 as last parameter the white matter surface will be estimated.\n\n\n";

       fprintf(stderr, usage_str, executable);
}

int
main(int argc, char *argv[])
{
        double               value, *values;
        Status               status;
        char                 *src_file, *out_file, *values_file;
        int                  i, p, n_steps, scale, n_objects, n_values, direction;
        File_formats         format;
        polygons_struct      *polygons, *tmp_polygon;
        Point                *new_pts;
        object_struct        **object_list, **objects;

        initialize_argument_processing(argc, argv);

        if (!get_string_argument( NULL, &src_file) ||
            !get_string_argument( NULL, &values_file) ||
            !get_string_argument( NULL, &out_file)) {
                usage(argv[0]);
                exit(EXIT_FAILURE);
        }
        get_int_argument(1, &direction);
        
        if (input_graphics_any_format(src_file, &format, &n_objects,
                                      &object_list) != OK)
                exit(EXIT_FAILURE);

        if (input_values_any_format(values_file, &n_values, &values) != OK)
                exit(EXIT_FAILURE);
      
        if (n_objects > 1) {
                fprintf(stderr,"Only one object allowed.\n");
                exit(EXIT_FAILURE);
        }

        polygons = get_polygons_ptr(object_list[0]);

        if (polygons->n_points != n_values) {
                fprintf(stderr,"Number of points differs from number of values.\n");
                exit(EXIT_FAILURE);
        }

        compute_polygon_normals(polygons);
        check_polygons_neighbours_computed(polygons);

        objects  = (object_struct **) malloc(sizeof(object_struct *));
        *objects = create_object(POLYGONS);
        tmp_polygon = get_polygons_ptr(*objects);
        
        /* use 10 steps to add thickness values to central surface and check in each step shape integrity */
        n_steps = 10;
        if (direction < 0) scale = -1*n_steps; else scale = n_steps; 
        for (i = 0; i < n_steps; i++) {
                copy_polygons(polygons, tmp_polygon);

                /* add half of thickness value in normal direction */ 
                for (p = 0; p < polygons->n_points; p++) {
                        Point_x(tmp_polygon->points[p]) += 0.5/scale*values[p]*Point_x(polygons->normals[p]);
                        Point_y(tmp_polygon->points[p]) += 0.5/scale*values[p]*Point_y(polygons->normals[p]);
                        Point_z(tmp_polygon->points[p]) += 0.5/scale*values[p]*Point_z(polygons->normals[p]);
                }
                /* get new points and check shape integrity */
                new_pts = tmp_polygon->points;
                check_polygons_shape_integrity(polygons, new_pts);
                polygons->points = new_pts;
        }

        /* smooth final surface */
        smooth_heatkernel(polygons, NULL, 1);

        compute_polygon_normals(polygons);

        if(output_graphics_any_format(out_file, format, 1, object_list, NULL) != OK)
                exit(EXIT_FAILURE);
                
        delete_object_list(1, objects);

        return(status != OK);
}
