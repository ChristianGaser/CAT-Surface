/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 *
 */

/*
 * Heat kernel smoothing is based on matlab code from Moo K. Chung:
 * Chung, M.K., Robbins,S., Dalton, K.M., Davidson, R.J., Evans, A.C. (2005) 
 * Cortical thickness analysis in autism via heat kernel smoothing. NeuroImage. 
 * http://www.stat.wisc.edu/~mchung/papers/ni_heatkernel.pdf
 */

#include <bicpl.h>

#include "CAT_Smooth.h"
#include "CAT_SurfaceIO.h"
#include "CAT_Curvature.h"

void
usage(char *executable)
{
        char *usage_str = "\n\
Usage: %s object_file output_file fwhm [gausscurv_threshold]\n\n\
     Diffusion smoothing of surface points w.r.t. gaussian curvature using\n\
     heat kernel.\n\n";

        fprintf(stderr, usage_str, executable);
}

int
main(int argc, char *argv[])
{
        char             *input_file, *output_file;
        int              n_objects;
        int              *n_neighbours, **neighbours;
        File_formats     format;
        object_struct    **object_list;
        polygons_struct  *polygons;
        double             fwhm;
        double           *gc_strength, gc_threshold;

        initialize_argument_processing(argc, argv);

        if (!get_string_argument(NULL, &input_file) ||
            !get_string_argument(NULL, &output_file)) {
                usage(argv[0]);
                exit(EXIT_FAILURE);
        }

        get_real_argument(25.0, &fwhm);
        get_real_argument(0.01, &gc_threshold);

        if (input_graphics_any_format(input_file, &format, &n_objects,
                                      &object_list) != OK ||
            n_objects != 1 || get_object_type(object_list[0]) != POLYGONS) {
                fprintf(stderr, "Error reading %s.\n", input_file);
                exit(EXIT_FAILURE);
        }

        polygons = get_polygons_ptr(object_list[0]);

        gc_strength = (double *)malloc(sizeof(double)*polygons->n_points);
            
        get_all_polygon_point_neighbours(polygons, &n_neighbours, &neighbours);

        get_polygon_vertex_curvatures_cg(polygons, n_neighbours, neighbours,
                                         0.0, 1, gc_strength);

        smooth_heatkernel(polygons, NULL, fwhm);
        
        compute_polygon_normals(polygons);

        if(output_graphics_any_format(output_file, format, 1, 
                        object_list) != OK)
                exit(EXIT_FAILURE);
        
        free(gc_strength);

        return(EXIT_SUCCESS);
}
