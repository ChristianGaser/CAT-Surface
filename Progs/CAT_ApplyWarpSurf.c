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

#include "CAT_Map.h"
#include "CAT_SurfaceIO.h"
#include "CAT_Surf.h"
#include "dartel/dartel.h"

#define  BINTREE_FACTOR   0.5

int deform = 0;

static ArgvInfo argTable[] = {
  {"-deform", ARGV_CONSTANT, (char *) TRUE, (char *) &deform,
     "Expect deformation field instead of flow field."},
   {NULL, ARGV_END, NULL, NULL, NULL}
};

void
usage(char *executable)
{
        char *usage_str = "\n\
Usage: %s  sphere_file FlowField_file output_sphere_file\n\n";

        fprintf(stderr, usage_str, executable);
}

int 
main(int argc, char *argv[])
{
        char              *surface_file, *output_surface_file, *flow_file;
        FILE              *infp;
        File_formats      format;
        int               n_objects, i, xy_size;
        polygons_struct   *sphere;
        object_struct     **objects;
        double            *inflow, *flow, *flow1;
        int               size_map[2], shift[2];

        /* get the arguments from the command line */
        if (ParseArgv(&argc, argv, argTable, 0) || argc < 4) {
                usage(argv[0]);
                fprintf(stderr, "     %s -help\n\n", argv[0]);
                exit(EXIT_FAILURE);
        }

        initialize_argument_processing(argc, argv);

        if (!get_string_argument(NULL, &surface_file) ||
            !get_string_argument(NULL, &flow_file) ||
            !get_string_argument(NULL, &output_surface_file)) {
                usage(argv[0]);
                exit(EXIT_FAILURE);
        }

        if (input_graphics_any_format(surface_file, &format,
                                      &n_objects, &objects ) != OK ||
            n_objects != 1 || get_object_type(objects[0]) != POLYGONS) {
                fprintf(stderr,
                        "Surface file must contain one polygons struct\n");
                exit(EXIT_FAILURE);
        }

        sphere = get_polygons_ptr(objects[0]);

        if ((infp = fopen(flow_file, "rb")) == NULL) {
                fprintf(stderr, "Error: Couldn't read file %s.\n", flow_file);
                exit(EXIT_FAILURE);
        }

        fread(&size_map, 2, sizeof(int), infp);
        fread(&shift, 2, sizeof(int), infp);
        xy_size = size_map[0] * size_map[1];

        flow  = (double *)malloc(sizeof(double) * xy_size * 2);
        flow1 = (double *)malloc(sizeof(double) * xy_size * 2);
        inflow  = (double *)malloc(sizeof(double) * xy_size * 2);
        fread(inflow, xy_size*2, sizeof(double), infp);
        fclose(infp);

        if (deform) {
                apply_warp(sphere, sphere, inflow, size_map, 0);  
        } else {
                expdef(size_map, 10, inflow, flow, flow1, (double *)0, (double *)0);  
                apply_warp(sphere, sphere, flow, size_map, 0);  
        }

        free(flow);
        free(flow1);
        free(inflow);

        if (output_graphics_any_format(output_surface_file, format, n_objects,
                                       objects, NULL) != OK)
                exit(EXIT_FAILURE);

        return(EXIT_SUCCESS);
}
