/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 *
 */

#include  <volume_io/internal_volume_io.h>
#include  <bicpl.h>

#include "CAT_Curvature.h"
#include "CAT_Blur2d.h"
#include "CAT_Surf.h"
#include "CAT_SurfaceIO.h"
    
int
main(int argc, char *argv[])
{
        STRING               object_filename, object2_filename, output_filename;
        FILE                 *file;
        File_formats         format;
        int                  poly, n_objects, n_iters;
        object_struct        **objects;
        polygons_struct      *polygons;
        double               surface_area, *area_values;

        initialize_argument_processing(argc, argv);

        if(!get_string_argument(NULL, &object_filename) || (!get_string_argument(NULL, &output_filename))) {
                printf("Usage: %s  object_file output_file [n_iters]\n", argv[0]);
                return(1);
        }

        get_int_argument(4000, &n_iters);
    
        if(input_graphics_any_format(object_filename, &format, &n_objects, &objects) != OK)
                return(1);

        if(n_objects != 1 || get_object_type(objects[0]) != POLYGONS) {
                printf("File must contain 1 polygons object.\n");
                return(1);
        }

        polygons = get_polygons_ptr(objects[0]);

        area_values = (double *)malloc(sizeof(double)*polygons->n_points);
        surface_area = get_area_of_points(polygons, area_values);
        free(area_values);
            
        areal_smoothing(polygons, 1.0, n_iters, 1, NULL, 0);
                        
        convert_ellipsoid_to_sphere_with_surface_area(polygons, surface_area);

        compute_polygon_normals(polygons);
        output_graphics_any_format(output_filename, format, 1, objects);

        delete_object_list(n_objects, objects);

        return(0);
}
