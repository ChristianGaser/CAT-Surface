/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 */

#include <volume_io/internal_volume_io.h>
#include <bicpl.h>

#include "CAT_Blur2d.h"
#include "CAT_Surf.h"
    
int
main(int argc, char *argv[])
{
        char                 *object_file, *output_file;
        FILE                 *fp;
        File_formats         format;
        int                  poly, n_objects, ptidx, vertidx, size;
        object_struct        **objects;
        polygons_struct      *polygons;
        double               poly_size, area, surface_area;
        double               *area_values;
        Point                points[MAX_POINTS_PER_POLYGON];
        BOOLEAN              all_values;
        signed char          *point_done;

        initialize_argument_processing(argc, argv);

        if (!get_string_argument(NULL, &object_file)) {
                fprintf(stderr, "Usage: %s  object_file [output_file]\n",
                        argv[0]);
                return(1);
        }

        if (input_graphics_any_format(object_file, &format,
                                      &n_objects, &objects) != OK)
                return(1);

        if (get_string_argument(NULL, &output_file)) {
                all_values = TRUE;
                if (open_file(output_file, WRITE_FILE, ASCII_FORMAT, &fp) != OK)
                        return(1);
        } else all_values = FALSE;

        if (n_objects != 1 || get_object_type(objects[0]) != POLYGONS) {
                printf("File must contain 1 polygons object.\n");
                return(1);
        }

        polygons = get_polygons_ptr(objects[0]);

        ALLOC(area_values, polygons->n_points);
        surface_area = get_area_of_points(polygons, area_values);
    
        if (all_values) {
                for (ptidx = 0; ptidx < polygons->n_points; ptidx++) {
                        if (output_real(fp, area_values[ptidx]) != OK ||
                            output_newline(fp) != OK)
                                break;
                }
                close_file(fp);
        }

        printf("Total surface area: %g\n", surface_area);

        delete_object_list(n_objects, objects);

        return(0);
}
