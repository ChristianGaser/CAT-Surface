/* Christian Gaser - christian.gaser@uni-jena.de
/* Department of Psychiatry
 * University of Jena
 *
 * most of the code is modified from
 * caret/caret_brain_set/BrainModelSurface.cxx.
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 *
 */

#include <volume_io/internal_volume_io.h>
#include <bicpl.h>
#include <float.h>

#include "CAT_Surf.h"

void
usage(char *executable)
{
        char *usage_str = "\n\
Usage: %s surface.obj sphere.obj [stop_at]\n\n\
     Maps a surface to a sphere using the caret inflating approach.\n\n\
     The inflating can be limited using stop_at (default 5), where\n\
       1 - Low smooth\n\
       2 - Inflating\n\
       3 - Very inflating\n\
       4 - High smoothing\n\
       5 - Ellipsoid\n";

        fprintf(stderr, usage_str, executable);
}

int
main(int argc, char *argv[])
{
        char             *input_file, *output_file;
        int              n_objects, i, stop_at;
        File_formats     format;
        object_struct    **object_list;
        polygons_struct  *polygons;
        BOOLEAN          enableFingerSmoothing = 1;
        int              fingerSmoothingIters;
        double           surfarea;
    
        initialize_argument_processing(argc, argv);

        if (!get_string_argument(NULL, &input_file) ||
            !get_string_argument(NULL, &output_file)) {
                usage(argv[0]);
                return(1);
        }

        get_int_argument(5, &stop_at);
    
        if (input_graphics_any_format(input_file, &format, &n_objects,
                                      &object_list) != OK || n_objects != 1 ||
            get_object_type(object_list[0]) != POLYGONS) {
                fprintf(stderr, "Error reading %s.\n", input_file);
                return(1);
        }

        polygons = get_polygons_ptr(object_list[0]);

        if (euler_characteristic(polygons) != 2) {
                fprintf(stderr, "Euler characteristic of %s must be 2.\n",
                            input_file);
        }
     
        surfarea = get_polygons_surface_area(polygons);

        /* low smooth */
	fprintf(stderr, "%20s\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b",
                "Low smoothing...    ");
        inflate_surface_and_smooth_fingers(polygons,
                              /* cycles */ 1,
          /* regular smoothing strength */ 0.2,
             /* regular smoothing iters */ 50,
                    /* inflation factor */ 1.0,
          /* finger comp/stretch thresh */ 3.0,
              /* finger smooth strength */ 1.0,
                 /* finger smooth iters */ 0);

        if (stop_at > 1) {
                /* inflated */
                fingerSmoothingIters = 0;
                if (enableFingerSmoothing)
                        fingerSmoothingIters = 30;
	        fprintf(stderr, "%20s\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b",
                        "Inflating...        ");
                inflate_surface_and_smooth_fingers(polygons,
                                      /* cycles */ 2,
                  /* regular smoothing strength */ 1.0,
                     /* regular smoothing iters */ 30,
                            /* inflation factor */ 1.4,
                  /* finger comp/stretch thresh */ 3.0,
                      /* finger smooth strength */ 1.0,
                         /* finger smooth iters */ fingerSmoothingIters);
        }                                             
    
        if (stop_at > 2) {
                /* very inflated */
                fprintf(stderr, "%20s\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b",
                        "Very inflating...   ");
                inflate_surface_and_smooth_fingers(polygons,
                  /*                     cycles */ 4,
                  /* regular smoothing strength */ 1.0,
                  /*    regular smoothing iters */ 30,
                  /*           inflation factor */ 1.1,
                  /* finger comp/stretch thresh */ 3.0,
                  /*     finger smooth strength */ 1.0,
                  /*        finger smooth iters */ 0);
        }                                             
    
        if (stop_at > 3) {
                /* high smooth */
                fingerSmoothingIters = 0;
                if (enableFingerSmoothing)
                        fingerSmoothingIters = 60;
	        fprintf(stderr, "%20s\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b",
                        "High smoothing...   ");
                inflate_surface_and_smooth_fingers(polygons,
                  /*                     cycles */ 6,
                  /* regular smoothing strength */ 1.0,
                  /*    regular smoothing iters */ 60,
                  /*           inflation factor */ 1.6,
                  /* finger comp/stretch thresh */ 3.0,
                  /*     finger smooth strength */ 1.0,
                  /*        finger smooth iters */ fingerSmoothingIters);
        }

        if (stop_at > 4) {
                /* ellipsoid */
                fprintf(stderr, "%20s\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b",
                        "Ellipsoid...        ");
                inflate_surface_and_smooth_fingers(polygons,
                  /*                     cycles */ 6,
                  /* regular smoothing strength */ 1.0,
                  /*    regular smoothing iters */ 50,
                  /*           inflation factor */ 1.4,
                  /* finger comp/stretch thresh */ 4.0,
                  /*     finger smooth strength */ 1.0,
                  /*        finger smooth iters */ fingerSmoothingIters);
                convert_ellipsoid_to_sphere_with_surface_area(polygons,
                                                              surfarea);
        }
    
        fprintf(stderr, "Done                \n");

        compute_polygon_normals(polygons);
        output_graphics_any_format(output_file, format, 1, object_list);
        return(0);
}
