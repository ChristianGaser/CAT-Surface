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

#include "CAT_Smooth.h"

int
separate_polygons(polygons_struct *polygons, int desired_index, Real *values);

void
usage(char *executable) {
        char *usage_str = "\n\
Usage: %s  input.obj input.txt output_prefix [which] \n\n\
     Separates polygons into its disjoint parts.\n\n";

        fprintf(stderr, usage_str, executable);
}

int
make_connected_components(polygons_struct *polygons, int point_classes[],
                          Real *values, int n_in_class[])
{
        int                point, edge, size;
        int                neigh, j;
        int                *n_neighbours, **neighbours;
        int                n_parts = 0, not_done = 16383;
        QUEUE_STRUCT(int)  queue;

        get_all_polygon_point_neighbours(polygons, &n_neighbours, &neighbours);

        for (point = 0; point < polygons->n_points; point++)
                point_classes[point] = not_done;

        for (point = 0; point < polygons->n_points; point++) {
                if (point_classes[point] != not_done)
                        continue;

                if (n_parts == 16383) {
                        n_parts++;
                        break;
                }

                point_classes[point] = n_parts;
                n_in_class[n_parts] = 1;

                if (n_neighbours[point] > 0) { 
                        for (j = 0; j < n_neighbours[point]; j++) {
                                if (point_classes[j] == not_done) {
                                        if (values[neighbours[point][j]] >
                                             0.0 && values[point] > 0.0) {
                                                 point_classes[j] = n_parts;
                                                 n_in_class[n_parts]++;
                                        } else {
                                                n_parts++;
                                        }
                                }
                        }
                }
            
                fprintf(stderr, "%d: %d\n", point, point_classes[point]);
        }

        return(n_parts);
}

int
separate_polygons(polygons_struct *polygons, int desired_index, Real *values)
{
        int       ind, p_ind, point, vertex, size, i, j, tmp;
        int       *new_point_ids, n_objects, comp, c;
        int       biggest;
        int       *point_classes;
        int       n_parts, *n_in_class, *ordered;

        ALLOC(point_classes, polygons->n_points);
        ALLOC(n_in_class, 16384);
        ALLOC(ordered, 16384);

        n_parts = make_connected_components(polygons, point_classes, values,
                                            n_in_class);

        for (i = 0; i < n_parts; i++)
                ordered[i] = i;

        for (i = 0; i < n_parts-1; i++) {
                biggest = i;
                for (j = i+1; j < n_parts; j++) {
                        if (n_in_class[ordered[j]] >
                            n_in_class[ordered[biggest]])
                        biggest = j;
                }

                tmp = ordered[i];
                ordered[i] = ordered[biggest];
                ordered[biggest] = tmp;
        }

        FREE(point_classes);
        FREE(n_in_class);
        FREE(ordered);

        return(n_objects);
}

int
main(int argc, char *argv[])
{
        char             *input_file, *output_prefix, *values_file;
        char             out_file[512];
        int              n_objects, n_out;
        int              i, desired_index, n_values;
        Real             *values;
        File_formats     format;
        object_struct    **object_list;
        polygons_struct  *polygons;

        initialize_argument_processing(argc, argv);

        if (!get_string_argument(NULL, &input_file) ||
            !get_string_argument(NULL, &values_file) ||
            !get_string_argument(NULL, &output_prefix)) {
                usage(argv[0]);
                exit(EXIT_FAILURE);
        }

        get_int_argument(-1, &desired_index);

        if (input_graphics_any_format(input_file, &format, &n_objects,
                                      &object_list) != OK ||
            n_objects < 1 || get_object_type(object_list[0]) != POLYGONS) {
                fprintf(stderr, "File must have a polygons structure.\n");
                exit(EXIT_FAILURE);
        }

        if (input_values_any_format(values_file, &n_values, &values) != OK) {
                fprintf(stderr, "Cannot read values.\n");
                exit(EXIT_FAILURE);
        }

        polygons = get_polygons_ptr(object_list[0]);

        check_polygons_neighbours_computed(polygons);

        n_out = separate_polygons(polygons, desired_index, values);
        fprintf(stderr,"%d\n",n_out);

        for (i = 0; i < n_out; i++) {
                fprintf(stderr,"%d\n",i);
        }

        return(EXIT_SUCCESS);
}