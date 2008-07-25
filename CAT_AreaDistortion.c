/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 */

#include <volume_io/internal_volume_io.h>
#include <bicpl.h>
#include <ParseArgv.h>

#include "CAT_Surf.h"

typedef enum { RATIO, LOG_RATIO, PERCENTAGE }
             Distortion_method;

/* argument defaults */
Distortion_method method = RATIO; /* area disortion method - default: ratio */

/* the argument table */
ArgvInfo argTable[] = {
  { "-ratio", ARGV_CONSTANT, (char *) 0, 
    (char *) &method,
    "Use area ratio (Default): ratio = a2/a1." },
  { "-log", ARGV_CONSTANT, (char *) 1, 
    (char *) &method,
    "Use log of area ratio: ratio = log10(a2/a1)." },
  { "-percent", ARGV_CONSTANT, (char *) 2, 
    (char *) &method,
    "Use percent change of area: ratio = 100*(a2-a1)/a1." },
  
  { NULL, ARGV_END, NULL, NULL, NULL }
};


int
main(int argc, char *argv[])
{
       char              *object_file, *object2_file, *output_file;
       FILE              *fp;
       File_formats      format;
       int               poly, n_objects, ptidx, size;
       object_struct     **objects, **objects2;
       polygons_struct   *polygons, *polygons2;
       double            poly_size, area, surface_area, surface_area2;
       double            *area_values, *area_values2, ratio, value;
       Point             points[MAX_POINTS_PER_POLYGON];

        /* Call ParseArgv */
        if (ParseArgv(&argc, argv, argTable, 0) || (argc < 3)) {
                fprintf(stderr,"\nUsage: %s [options] object_file object_file2 output_file\n", argv[0]);
                fprintf( stderr,"\nCalculate area distortion between two surfaces.\n");
                fprintf(stderr, "       %s -help\n\n", argv[0]);
                return(1);
        }

        initialize_argument_processing(argc, argv);

        if (!get_string_argument(NULL, &object_file) ||
            !get_string_argument(NULL, &object2_file) ||
            !get_string_argument(NULL, &output_file)) {
                fprintf(stderr,
                        "Usage: %s  object_file object_file2 output_file\n",
                        argv[0]);
                return(1);
        }

        if (input_graphics_any_format(object_file, &format,
                                      &n_objects, &objects) != OK)
                return(1);

        if (n_objects != 1 || get_object_type(objects[0]) != POLYGONS) {
                printf("File must contain 1 polygons object.\n");
                return(1);
        }

        if (input_graphics_any_format(object2_file, &format,
                                      &n_objects, &objects2) != OK)
                return(1);

        if (n_objects != 1 || get_object_type(objects2[0]) != POLYGONS) {
                printf("File must contain 1 polygons object.\n");
                return(1);
        }

        if (open_file(output_file, WRITE_FILE, ASCII_FORMAT, &fp) != OK)
                return(1);

        polygons = get_polygons_ptr(objects[0]);
        polygons2 = get_polygons_ptr(objects2[0]);

        ALLOC(area_values, polygons->n_points);
        surface_area = get_area_of_points(polygons, area_values);

        ALLOC(area_values2, polygons2->n_points);
        surface_area2 = get_area_of_points(polygons2, area_values2);

        ratio = surface_area / surface_area2;

        for (ptidx = 0; ptidx < polygons->n_points; ptidx++) {
                switch (method) {
                case RATIO:
                        value = ratio * area_values2[ptidx] /
                                area_values[ptidx];
                        break;
                case LOG_RATIO:
                       value = log10(ratio * area_values2[ptidx] /
                                     area_values[ptidx]);
                       break;
                case PERCENTAGE:
                       value = 100 * (ratio * area_values2[ptidx] -
                                      area_values[ptidx]) /
                               area_values[ptidx];
                       break;
                }

                if (output_real(fp, value) != OK || output_newline(fp) != OK)
                        break;
        }

        close_file(fp);

        delete_object_list(n_objects, objects);
        delete_object_list(n_objects, objects2);

        return(0);
}

