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

#include "CAT_SurfaceIO.h"

/* argument defaults */
char *outfile  = NULL;
char *rms_file = NULL;

/* the argument table */
static ArgvInfo argTable[] = {
  {"-avg", ARGV_STRING, (char *) 1, (char *) &outfile, 
     "Define name for average file."},
  {"-rms", ARGV_STRING, (char *) 1, (char *) &rms_file, 
     "Define name for root mean square error file."},
   {NULL, ARGV_END, NULL, NULL, NULL}
};

void
usage(char *executable)
{
        char *usage_str = "\n\
Usage: %s -avg average_surface|-rms rms_values_file surface_file1 surface_file2 [... surface_fileN]\n\n\
     Calculate average and root mean square error of surfaces.\n\n";

        fprintf(stderr, usage_str, executable);
}

int
main(int argc, char *argv[])
{
        Status           status;
        char             *infile;
        int              i, n_objects;
        File_formats     format;
        int              n_sets, n_pts, n_avg_pts;
        object_struct    *out_object;
        object_struct    **object_list;
        Point            *pts, *avg_pts, *sqr_pts, sqr_pt;
        double           value, *values;

        /* Call ParseArgv */
        if (ParseArgv(&argc, argv, argTable, 0) ||
            (rms_file == NULL && outfile == NULL)) {
                fprintf(stderr, "Either average or rms file must be defined.\n");
                usage(argv[0]);
                fprintf(stderr, "     %s -help\n\n", argv[0]);
                exit(EXIT_FAILURE);
        }

        initialize_argument_processing(argc, argv);

        n_sets = 0;

        while (get_string_argument(NULL, &infile)) {
                if (input_graphics_any_format(infile, &format, &n_objects,
                                              &object_list) != OK) {
	                fprintf(stderr, "Could not read input %s\n", infile);
	                exit(EXIT_FAILURE);
	        }

                n_pts = get_object_points(object_list[0], &pts);

                if (n_sets == 0) {
                        out_object = object_list[0];
                        n_avg_pts = n_pts;
                        avg_pts = pts;
                        if (rms_file != NULL) {
                                sqr_pts = (Point *) malloc(sizeof(Point) * n_pts);
                                for (i = 0; i < n_pts; i++)
                                        POINT_EXP2(sqr_pts[i], pts[i], *, pts[i]);
                        }
                } else {
                        if (n_pts != n_avg_pts) {
                                printf("N points mismatch\n");
                                exit(EXIT_FAILURE);
                        }

                        for (i = 0; i < n_pts; i++) {
                                ADD_POINTS(avg_pts[i], avg_pts[i], pts[i]);
                                if (rms_file != NULL) {
                                        POINT_EXP2(sqr_pt, pts[i], *, pts[i]);
                                        ADD_POINTS(sqr_pts[i], sqr_pts[i], sqr_pt);
                                }
                        }
                }

                printf("%d:  %s\n", n_sets, infile);

                if (n_sets > 0)
                    delete_object_list(n_objects, object_list);

                ++n_sets;
        }

        if (n_sets == 0) {
                fprintf(stderr, "No input?!  No output for you!!\n");
                exit(EXIT_FAILURE);
        }

        for (i = 0; i < n_pts; i++) {
                SCALE_POINT(avg_pts[i], avg_pts[i], (1.0/n_sets));
        }

        if (outfile != NULL) {
                if (get_object_type(out_object) == POLYGONS)
                        compute_polygon_normals(get_polygons_ptr(out_object));

                status = output_graphics_any_format(outfile, ASCII_FORMAT,
                                            1, &out_object, NULL);
        }
        
        if (rms_file != NULL) {
                values = (double *) malloc(sizeof(double) * n_pts);
                for (i = 0; i < n_pts; i++) {
                        SCALE_POINT(sqr_pts[i], sqr_pts[i], (1.0/(n_sets-1.0)));
                        POINT_EXP2(avg_pts[i], avg_pts[i], *, avg_pts[i]);
                        SCALE_POINT(avg_pts[i], avg_pts[i], (n_sets/(n_sets-1.0)));
                        SUB_POINTS(sqr_pts[i], sqr_pts[i], avg_pts[i]);
                        value = sqrt(Point_x(sqr_pts[i]) + Point_y(sqr_pts[i]) + Point_z(sqr_pts[i]));
                        values[i] = value;
                }
                output_values_any_format(rms_file, n_pts, values, TYPE_DOUBLE);
                free(values);
                free(sqr_pts);
        }
        
  return(EXIT_SUCCESS);
}
