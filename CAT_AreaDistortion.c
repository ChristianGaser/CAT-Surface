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
        double            area, area2;
        double            *ad_values, ratio, value;
        double            distortion = 0;
        Point             pts[MAX_POINTS_PER_POLYGON];
        Point             pts2[MAX_POINTS_PER_POLYGON];
        int               *n_polys;

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
                ALLOC(ad_values, polygons->n_items);
                n_obj = polygons->n_items;
        } else {
                ALLOC(ad_values, polygons->n_points);
                ALLOC(n_polys, polygons->n_points);
                for (i = 0; i < polygons->n_points; i++) {
                        n_polys[i] = 0;
                        ad_values[i] = 0;
                }
                n_obj = polygons->n_points;
        }

        ratio = get_polygons_surface_area(polygons) /
                get_polygons_surface_area(polygons2);

        for (poly = 0; poly < polygons->n_items; poly++) {
                size = get_polygon_points(polygons, poly, pts);
                area = get_polygon_surface_area(size, pts);

                size = get_polygon_points(polygons2, poly, pts2);
                area2 = get_polygon_surface_area(size, pts2);

                size = GET_OBJECT_SIZE(*polygons, poly);

                switch (method) {
                case RATIO:
                        value = ratio * area2 / area;
                        break;
                case LOG_RATIO:
                       value = log10(ratio * area2 / area);
                       break;
                case PERCENTAGE:
                       value = 100 * (ratio * (area2 - area)) / area;
                       break;
                }

                if (value < PINF && value > NINF)
                        distortion += fabs(value);

                if (PerPoly) {
                        ad_values[poly] = value;
                } else {
                        for (vertidx = 0; vertidx < size; vertidx++) {
                                i = polygons->indices[
                                              POINT_INDEX(polygons->end_indices,
                                              poly, vertidx)];
                                n_polys[i]++;
                                ad_values[i] += value;
                        }
                }
        }


        for (i = 0; i < n_obj; i++) {
                if (!PerPoly && n_polys[i] > 0) /* average distortion */
                        ad_values[i] /= n_polys[i];

                if (output_real(fp, ad_values[i]) != OK ||
                    output_newline(fp) != OK)
                        return(1);
        }

        printf("Area distortion = %f\n", distortion / polygons->n_items);

        close_file(fp);

        delete_object_list(n_objects, objects);
        delete_object_list(n_objects, objects2);

        FREE(ad_values);
        if (!PerPoly)
                FREE(n_polys);

        return(0);
}
