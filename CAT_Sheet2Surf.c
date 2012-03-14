/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 *
 */

#include <volume_io/internal_volume_io.h>
#include <bicpl.h>

#include "CAT_Map2d.h"
#include "CAT_Blur2d.h"
#include "CAT_SurfaceIO.h"

#define BINTREE_FACTOR 0.5
    
void
usage(char *executable)
{
        char *usage_str = "\n\
Usage: %s surface.obj  input.pgm output.txt [0|1]\n\n\
     Maps an PGM-image to a surface. The output is saved as texture values.\n\n\
     Optionally linear interpolation can be skipped with a 0 as 4th argument\n\n";

        fprintf(stderr, usage_str, executable);
}

int
main(int argc, char *argv[])
{
        char                 *surface_file, *output_file, *input_file;
        File_formats         format;
        polygons_struct      *polygons;
        int                  i, n_objects;
        int                  degree, poly, size_map[2];
        int                  ind, n_done, interpolate;
        double               *image;
        object_struct        **objects;
        double               *values;

        initialize_argument_processing(argc, argv);

        if (!get_string_argument(NULL, &surface_file) ||
            !get_string_argument(NULL, &input_file) ||
            !get_string_argument(NULL, &output_file)) {
                usage(argv[0]);
                exit(EXIT_FAILURE);
        }

        get_int_argument(1, &interpolate);

        if (input_graphics_any_format(surface_file, &format,
                                      &n_objects, &objects) != OK)
                exit(EXIT_FAILURE);

        /* check that the surface file contains a polyhedron */
        if (n_objects != 1 || get_object_type(objects[0]) != POLYGONS) {
                printf("Surface file must contain 1 polygons object.\n");
                exit(EXIT_FAILURE);
        }

        /* get a pointer to the surface */
        polygons = get_polygons_ptr(objects[0]);
    
        values = (double *) malloc(sizeof(double) * polygons->n_points);

        if ((image = read_pgm(input_file, &size_map[0], &size_map[1])) == NULL)
                exit(EXIT_FAILURE);
        
        map_sheet2d_to_sphere(image, values, polygons, interpolate, size_map);

        output_values_any_format(output_file, polygons->n_points,
                                 values, TYPE_DOUBLE);

        delete_object_list(n_objects, objects);

        free(image);
        free(values);

        return(EXIT_SUCCESS);
}
