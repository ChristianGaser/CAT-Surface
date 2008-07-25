/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 */

#include <volume_io/internal_volume_io.h>
#include <bicpl.h>

void
usage(char *executable)
{
        char *usage_str = "\n\
Usage: %s object_file x y z values1.txt [values2.txt .. valuesn.txt]\n\n\
     Plot values at coordinate x y z for each file. The object file is used to\n\
     link the coordinates at the surface model to the values.\n\n";

        fprintf(stderr, usage_str, executable);
}

int
main(int argc, char *argv[])
{
        char              *object_file, *input_file;
        FILE              *infp;
        File_formats      format;
        int               poly, n_objects, pidx, min_index, n_files;
        int               i, j;
        object_struct     **objects;
        polygons_struct   *polygons;
        Real              x, y, z, dist, min_dist, value;
        Point             point;
        char              line[256];

        initialize_argument_processing(argc, argv);

        if (!get_string_argument(NULL, &object_file)) {
                usage(argv[0]);
                return(1);
        }

        if (input_graphics_any_format(object_file, &format,
                                      &n_objects, &objects ) != OK)
                return(1);

        if (!get_real_argument(0.0, &x) ||
            !get_real_argument(0.0, &y) ||
            !get_real_argument(0.0, &z)) {
                usage(argv[0]);
                return(1);
        }

        n_files = argc - 5;
        
        if (n_objects != 1 || get_object_type(objects[0]) != POLYGONS) {
                printf("File must contain 1 polygons object.\n");
                return(1);
        }

        fill_Point(point, x, y, z);
    
        polygons = get_polygons_ptr(objects[0]);
    
        min_dist = 1e15;
        for (pidx = 0; pidx <  polygons->n_points; pidx++) {
                dist = distance_between_points(&polygons->points[pidx], &point);
                if (dist < min_dist) {
                        min_dist = dist;
                        min_index = pidx;
                }
        }

        printf("Closest distance of %3.2f found at index %d\n",
               min_dist, min_index);
    
        for (i = 0; i < n_files; i++) {
                get_string_argument(NULL, &input_file);
                if ((infp = fopen(input_file, "r")) == 0) {
                        fprintf(stderr, "Couldn't open file %s.\n",
                                input_file);
                        return(0);
                }
                for (j = 0; j <  min_index+1; j++)
                        fgets(line, 256, infp);
                printf("%s",line);
                close_file(infp);
        }

        delete_object_list(n_objects, objects);
        return(0);
}
