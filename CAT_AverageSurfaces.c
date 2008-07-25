/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 */

/* modified from average_objects.c to enable different surface formats */

#include <bicpl.h>

int
main(int argc, char *argv[])
{
        Status           status;
        char             *infile, *outfile;
        int              i, n_objects;
        File_formats     format;
        int              n_sets, n_pts, n_set_pts;
        object_struct    *out_object;
        object_struct    **object_list;
        Point            *pts, *set_pts;

        initialize_argument_processing(argc, argv);

        if (!get_string_argument(NULL, &outfile)) {
            fprintf(stderr, "Usage: %s output\n", argv[0]);
            fprintf(stderr, "         [input1] [input2] ...\n");
            return(1);
        }

        n_sets = 0;

        while (get_string_argument(NULL, &infile)) {
                if (input_graphics_any_format(infile, &format, &n_objects,
                                              &object_list) != OK) {
	                fprintf(stderr, "Could not read input %s\n", infile);
	                return(2);
	        }

                n_pts = get_object_points(object_list[0], &pts);

                if (n_sets == 0) {
                        out_object = object_list[0];
                        n_set_pts = n_pts;
                        set_pts = pts;
                } else {
                        if (n_pts != n_set_pts) {
                                printf("N points mismatch\n");
                                return(1);
                        }

                        for (i = 0; i < n_pts; i++) {
                                ADD_POINTS(set_pts[i], set_pts[i], pts[i]);
                        }
                }

                printf("%d:  %s\n", n_sets, infile);

                if (n_sets > 0)
                    delete_object_list(n_objects, object_list);

                ++n_sets;
        }

        if (n_sets == 0) {
	        fprintf(stderr, "No input?!  No output for you!!\n");
	        return(1);
        }

        for (i = 0; i < n_pts; i++) {
	        SCALE_POINT(set_pts[i], set_pts[i], (1.0/n_sets));
        }

        if (get_object_type(out_object) == POLYGONS)
                compute_polygon_normals(get_polygons_ptr(out_object));

        status = output_graphics_any_format(outfile, ASCII_FORMAT,
                                            1, &out_object);

        return(status != OK);
}
