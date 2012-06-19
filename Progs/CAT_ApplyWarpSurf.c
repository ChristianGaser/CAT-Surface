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

#define  BINTREE_FACTOR   0.5

/* defaults */
char *sphere_file = NULL;

static ArgvInfo argTable[] = {
  {"-sphere", ARGV_STRING, (char *) 1, (char *) &sphere_file, 
     "Sphere for input surface."},
   {NULL, ARGV_END, NULL, NULL, NULL}
};

void
usage(char *executable)
{
        char *usage_str = "\n\
Usage: %s  surface.obj FlowField output.obj\n\n";

        fprintf(stderr, usage_str, executable);
}

int 
main(int argc, char *argv[])
{
        char              *surface_file, *output_file, *flow_file;
        FILE              *infp;
        File_formats      format;
        int               n_objects, i, xy_size, shift[2];
        polygons_struct   *polygons, *sphere;
        object_struct     **objects;
        double            *inflow, *flow, *flow1;
        int               size_map[2];

        initialize_argument_processing(argc, argv);

        if (!get_string_argument(NULL, &surface_file) ||
            !get_string_argument(NULL, &flow_file) ||
            !get_string_argument(NULL, &output_file)) {
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

        polygons = get_polygons_ptr(objects[0]);

        /* read sphere for input surface */
        if (sphere_file != NULL) {
                if (input_graphics_any_format(sphere_file, &format,
                                      &n_objects, &objects) != OK)
                        exit(EXIT_FAILURE);
                sphere = get_polygons_ptr(objects[0]);
        } else 
                sphere = (polygons_struct *) 0;

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

        expdef(size_map, 10, inflow, flow, flow1, (double *)0, (double *)0);  
        free(flow1);
        free(inflow);

        apply_warp(polygons, sphere, flow, size_map, 0);  

        if (output_graphics_any_format(output_file, format, n_objects,
                                       objects) != OK)
                exit(EXIT_FAILURE);

        free(flow);

        return(EXIT_SUCCESS);
}