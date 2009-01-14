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

#define  BINTREE_FACTOR   0.5

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
        polygons_struct   *polygons;
        object_struct     **objects;
        double            *inflow, *flow, *flow1;
        int               size_map[2];

        initialize_argument_processing(argc, argv);

        if (!get_string_argument(NULL, &surface_file) ||
            !get_string_argument(NULL, &flow_file) ||
            !get_string_argument(NULL, &output_file)) {
                usage(argv[0]);
                return(1);
        }

        if (input_graphics_any_format(surface_file, &format,
                                      &n_objects, &objects ) != OK ||
            n_objects != 1 || get_object_type(objects[0]) != POLYGONS) {
                fprintf(stderr,
                        "Surface file must contain one polygons struct\n");
                return(1);
        }

        polygons = get_polygons_ptr(objects[0]);

        if ((infp = fopen(flow_file, "rb")) == NULL) {
                fprintf(stderr, "Error: Couldn't read file %s.\n", flow_file);
                return(1);
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

        apply_warp(polygons, flow, size_map, shift);  

        if (output_graphics_any_format(output_file, format, n_objects,
                                       objects) != OK)
                return(1);

        free(flow);

        return(0);
}
