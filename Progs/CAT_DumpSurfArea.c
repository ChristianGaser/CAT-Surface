/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 *
 */

#include <bicpl.h>
#include <ParseArgv.h>

#include "CAT_Smooth.h"
#include "CAT_Surf.h"
#include "CAT_SurfaceIO.h"
    
/* argument defaults */

char *sphere_file = NULL;
int use_log = 0;

/* the argument table */
ArgvInfo argTable[] = {
  {"-sphere", ARGV_STRING, (char *) 1, (char *) &sphere_file,
     "If a sphere is given the local surface area is based on a re-parameterized tetrahedral sphere."},
  { "-log", ARGV_INT, (char *) 1,
    (char *) &use_log,
    "Obtain log10-scaled surface area." },
  { NULL, ARGV_END, NULL, NULL, NULL }
};


void
usage(char *executable)
{
        static char *usage_str = "\n\
Usage: %s surface_file output_values_file\n\
Calculate local surface area. If a sphere is given the local surface area is based on a re-parameterized tetrahedral sphere.\n\n\n";

       fprintf(stderr, usage_str, executable);
}

int
main(int argc, char *argv[])
{
        char                 *surface_file, *output_values_file;
        File_formats         format;
        int                  poly, n_objects, i;
        object_struct        **objects, **sphere_objects;
        polygons_struct      *polygons, *sphere;
        double               area, surface_area;
        double               *area_values;
        Point                points[MAX_POINTS_PER_POLYGON];
        signed char          *point_done;

        initialize_argument_processing(argc, argv);

        if (ParseArgv(&argc, argv, argTable, 0) || argc < 2) {
                usage(argv[0]);
                fprintf(stderr, "       %s -help\n\n", argv[0]);
                exit(EXIT_FAILURE);
        }

        if (!get_string_argument(NULL, &surface_file) ||
            !get_string_argument(NULL, &output_values_file)) {
                usage(argv[0]);
                fprintf(stderr, "       %s -help\n\n", argv[0]);
                exit(EXIT_FAILURE);
        }

        if (input_graphics_any_format(surface_file, &format,
                                      &n_objects, &objects) != OK)
                exit(EXIT_FAILURE);

        if (n_objects != 1 || get_object_type(objects[0]) != POLYGONS) {
                fprintf(stderr,"File must contain 1 polygons object.\n");
                exit(EXIT_FAILURE);
        }

        polygons = get_polygons_ptr(objects[0]);

        area_values = (double *) malloc(sizeof(double) * polygons->n_points);
        
       if (sphere_file != NULL) {
                if (input_graphics_any_format(sphere_file, &format, &n_objects,
                                              &sphere_objects) != OK)
                        exit(EXIT_FAILURE);
                sphere = get_polygons_ptr(*sphere_objects);
                
                if (polygons->n_items != sphere->n_items) {
                        fprintf(stderr,"Surface and sphere must have same size.\n");
                        exit(EXIT_FAILURE);
                }
                
                surface_area = get_area_of_points_normalized_to_sphere(polygons, sphere, area_values);
        } else  surface_area = get_area_of_points(polygons, area_values);

        if (use_log) 
                for(i=0; i<polygons->n_points; i++)
                        area_values[i] = log10(area_values[i]);

        output_values_any_format(output_values_file, polygons->n_points,
                                         area_values, TYPE_DOUBLE);
        
        if (use_log) 
                printf("Total log. surface area: %g\n", log10(surface_area));
        else
                printf("Total surface area: %g\n", surface_area);

        delete_object_list(n_objects, objects);
        free(area_values);

        return(EXIT_SUCCESS);
}
