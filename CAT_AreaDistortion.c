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


#define PINF  1.7976931348623157e+308 /* for doubles */
#define NINF -1.7976931348623157e+308 /* for doubles */

typedef enum { RATIO, LOG_RATIO, PERCENTAGE }
             Distortion_method;

/* argument defaults */
Distortion_method method = RATIO; /* area disortion method - default: ratio */
BOOLEAN PerPoly = 0; /* by default, generate area distortion per point */


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
  { "-polygon", ARGV_CONSTANT, (char *) 1, 
    (char *) &PerPoly,
    "Calculate area distortion on a per-polygon basis." },
  
  { NULL, ARGV_END, NULL, NULL, NULL }
};


int
main(int argc, char *argv[])
{
        char              *object_file, *object2_file, *output_file;
        FILE              *fp;
        File_formats      format;
        int               poly, n_objects, n_obj, i, size, vertidx;
        object_struct     **objects, **objects2;
        polygons_struct   *polygons, *polygons2;
        double            area, area2, surface_area, surface_area2;
        double            *area_values, *area_values2, ratio, value;
        double            distortion = 0;
        Point             pts[MAX_POINTS_PER_POLYGON];
        Point             pts2[MAX_POINTS_PER_POLYGON];
        signed char       *point_done;

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

        if (polygons->n_items != polygons2->n_items ||
            polygons->n_points != polygons2->n_points) {
                fprintf(stderr, "Input polygons don't match. Exiting.\n");
                return(1);
        }

        if (PerPoly) {
                ALLOC(area_values, polygons->n_items);
                ALLOC(area_values2, polygons2->n_items);
                n_obj = polygons->n_items;
        } else {
                ALLOC(area_values, polygons->n_points);
                ALLOC(area_values2, polygons2->n_points);
                ALLOC(point_done, polygons->n_points);
                for (i = 0; i < polygons->n_points; i++)
                        point_done[i] = FALSE;
                n_obj = polygons->n_points;
        }

        for (poly = 0; poly < polygons->n_items; poly++) {
                size = get_polygon_points(polygons, poly, pts);
                area = get_polygon_surface_area(size, pts);
                surface_area += area;

                size = get_polygon_points(polygons2, poly, pts2);
                area2 = get_polygon_surface_area(size, pts2);
                surface_area2 += area2;

                size = GET_OBJECT_SIZE(*polygons, poly);

                if (PerPoly) {
                        area_values[poly] = area;
                        area_values2[poly] = area2;
                } else {
                        for (vertidx = 0; vertidx < size; vertidx++) {
                                i = polygons->indices[
                                              POINT_INDEX(polygons->end_indices,
                                              poly, vertidx)];
                                if (!point_done[i]) {
                                        point_done[i] = TRUE;
                                        area_values[i] = area;
                                        area_values2[i] = area2;
                                }
                        }
                }
        }

        ratio = surface_area / surface_area2;

        for (i = 0; i < n_obj; i++) {
                switch (method) {
                case RATIO:
                        value = ratio * area_values2[i] / area_values[i];
                        break;
                case LOG_RATIO:
                       value = log10(ratio * area_values2[i] / area_values[i]);
                       break;
                case PERCENTAGE:
                       value = 100 *
                               (ratio * (area_values2[i] - area_values[i])) /
                               area_values[i];
                       break;
                }
                if (value < PINF && value > NINF)
                        distortion += fabs(value);

                if (output_real(fp, value) != OK || output_newline(fp) != OK)
                        break;
        }

        printf("Area distortion = %f\n", distortion / n_obj);

        close_file(fp);

        delete_object_list(n_objects, objects);
        delete_object_list(n_objects, objects2);

        FREE(area_values);
        FREE(area_values2);
        if (!PerPoly)
                FREE(point_done);

        return(0);
}
