/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 */

#include <volume_io/internal_volume_io.h>
#include <bicpl.h>

#include "CAT_Map2d.h"
#include "CAT_Surf.h"

#define BINTREE_FACTOR   0.5
        
void
usage(char *executable)
{
        char *usage_str = "\n\
Usage: %s input.txt ShiftField output.txt\n\n\
     Applies deformations of warping to surface values.\n\n";

        fprintf(stderr, usage_str, executable);
}

int
main(int argc, char *argv[])
{
        char                 *output_file, *vector_file, *values_file;
        FILE                 *infp;
        File_formats         format;
        polygons_struct      unit_sphere;
        int                  i, j, x, y;
        int                  degree, poly, size, ind, n_done;
        int                  n_values;
        double               *inflow, *flow, *flow1;
        BOOLEAN              use_volume;
        Point                unit_point, on_sphere_point, centre;
        Point                poly_points[1000], centroid;
        double               u, v, *values, *input_values, **sheet;
        double               weights[1000], value, indx, indy;
        object_struct        *object;    
        Vector               normal;
        progress_struct      progress;
        double               inflow_x, inflow_y, ux, vy;
        int                  size_map[2], shift[2];

        initialize_argument_processing(argc, argv);

        if (!get_string_argument(NULL, &values_file) ||
            !get_string_argument(NULL, &vector_file) ||
            !get_string_argument(NULL, &output_file)) {
                usage(argv[0]);
                return(1);
        }


        if (input_values_any_format(values_file, &n_values, &input_values) != OK)
                return(1);

        if ((infp = fopen(vector_file, "rb")) == NULL) {
                fprintf(stderr, "Error: Couldn't read file %s.\n", vector_file);
                return(1);
        }

        fread(&size_map, 2, sizeof(int), infp);
        fread(&shift, 2, sizeof(int), infp);

        inflow  = (double *)malloc(sizeof(double)*size_map[0]*size_map[1]*2);
        flow    = (double *)malloc(sizeof(double)*size_map[0]*size_map[1]*2);
        flow1   = (double *)malloc(sizeof(double)*size_map[0]*size_map[1]*2);
        fread(inflow, size_map[0]*size_map[1]*2, sizeof(double), infp);
        fclose(infp);

        expdef(size_map, 10, inflow, flow, flow1, (double *)0, (double *)0);  
        free(flow1);
        free(inflow);

        /* create unit sphere with same number of triangles as skin surface */
        fill_Point(centre, 0.0, 0.0, 0.0);

        create_tetrahedral_sphere(&centre, 1.0, 1.0, 1.0,
                                  2*(n_values-2), &unit_sphere);

        create_polygons_bintree(&unit_sphere,
                                round((double) unit_sphere.n_items *
                                      BINTREE_FACTOR));

        values = (double *)malloc(sizeof(double)*unit_sphere.n_points);
    
        initialize_progress_report(&progress, FALSE, size_map[0],
                                   "Mapping to sheet");

        ALLOC2D(sheet, size_map[0],size_map[1]);

        inflow_x = (double)size_map[0] - 1.0;
        inflow_y = (double)size_map[1] - 1.0;

        initialize_progress_report(&progress, FALSE, unit_sphere.n_points,
                                   "Mapping to sphere");
	
        /* remap to sphere */
        for (i = 0; i < unit_sphere.n_points; i++ ) {
                point_to_uv(&unit_sphere.points[i], &u, &v);
        
                indx = u*inflow_x;
    	        indy = v*inflow_y;    
    	        ind  = (int)round(indx) + size_map[0]*(int)round(indy);

    	        ux = (flow[ind] - 1.0 - indx + shift[0])/inflow_x;
    	        vy = (flow[ind + size_map[0]*size_map[1]] - 1.0 - indy +
                     shift[1])/inflow_y;
    
    	        u += ux;
    	        v += vy;

    	        // wrap borders
                while (u < 0.0)  u += 1.0;
               	while (u >= 1.0) u -= 1.0;
                if (v < 0.0)     v = 0.0;
                if (v > 1.0)     v = 1.0;
	
                uv_to_point(u, v, &unit_point);

                poly = find_closest_polygon_point(&unit_point, &unit_sphere,
                                                  &on_sphere_point);
      
                size = get_polygon_points(&unit_sphere, poly, poly_points);

                get_polygon_interpolation_weights(&on_sphere_point, size,
                                                  poly_points, weights);

                value = 0.0;
                for (j = 0; j < size; j++) {
                        ind = unit_sphere.indices[
                              POINT_INDEX(unit_sphere.end_indices,poly,j)];
                        value += weights[j] * input_values[ind];
                }
                values[i] = value;
 
                update_progress_report(&progress, i + 1);
        }

        output_values_any_format(output_file, unit_sphere.n_points, values);

        terminate_progress_report(&progress);

        delete_polygons(&unit_sphere);
        free(flow);
        free(values);
 
        return(0);
}
