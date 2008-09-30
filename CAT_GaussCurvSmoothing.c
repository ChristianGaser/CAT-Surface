/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 */

#include  <volume_io/internal_volume_io.h>
#include  <bicpl.h>

#include "CAT_Curvature.h"
#include "CAT_Blur2d.h"
#include "CAT_Surf.h"
#include "CAT_SurfaceIO.h"
    
void
gaussian_curv_smoothing(polygons_struct *polygons, double *strength, int iters)
{
        int     i, j, k, l;
        int     n1, n2, next;
        int     *n_neighbours, **neighbours;
        double  *area_values;
        Point   pts[1000];
        double tileAreas[32], tileCenters[32*3];
        double xyz[3], pt1[3], pt2[3], pt3[3];
        double totalArea, weight, invstr;

        create_polygon_point_neighbours(polygons, TRUE, &n_neighbours,
                                        &neighbours, NULL, NULL);
    
        for (k = 1; k < iters; k++) {

                for (i = 0; i < polygons->n_points; i++) {
        
                        if (strength[i] == 0)
                                continue; /* skip this point */

                        totalArea = 0.0;    
                    
                        /* Get 2 consecutive neighbors of this node */
                        for (j = 0; j < n_neighbours[i]; j++) {
                                n1 = neighbours[i][j];
                                next = j + 1;
                                if (next >= n_neighbours[i])
                                        next = 0;
                                n2 = neighbours[i][next];
            
                                /* Area of the triangle */
                                pts[0] = polygons->points[i];
                                pts[1] = polygons->points[n1];
                                pts[2] = polygons->points[n2];
                                tileAreas[j] = get_polygon_surface_area(3, pts);
                                totalArea += tileAreas[j];
                
                                /* Save center of this tile */
                                to_array(&polygons->points[i], pt1);
                                to_array(&polygons->points[n1], pt2);
                                to_array(&polygons->points[n2], pt3);
                                for (l = 0; l < 3; l++) {
                                        tileCenters[j*3+l] = (pt1[l] +
                                                              pt2[l] +
                                                              pt3[l]) / 3.0;
                                }
                        }

                        /* Compute the influence of the neighboring nodes */
                        for (j = 0; j < 3; j++)
                                xyz[j] = 0.0;
                        for (j = 0; j < n_neighbours[i]; j++) {
                                if (tileAreas[j] > 0.0) {
                                        weight = tileAreas[j] / totalArea;
                                        for (l = 0; l < 3; l++) {
                                                xyz[l] += weight *
                                                          tileCenters[j*3+l];
                                        }
                                }
                        }
                        /* Update the nodes position */
                        to_array(&polygons->points[i], pt1);
                        invstr = 1.0 - strength[i];
                        for (l = 0; l < 3; l++) {
                                pt1[l] = (pt1[l] * invstr) +
                                         (xyz[l] * strength[i]);
                        }
                        from_array(pt1, &polygons->points[i]);
                }
        
        }
}


int
main(int argc, char *argv[])
{
        STRING               object_filename, output_filename;
        FILE                 *file;
        File_formats         format;
        int                  i, poly, n_objects, n_iters;
        int                  *n_neighbours, **neighbours;
        object_struct        **objects;
        polygons_struct      *polygons;
        double               *gc_strength, gc_threshold;

        initialize_argument_processing(argc, argv);

        if(!get_string_argument(NULL, &object_filename) || (!get_string_argument(NULL, &output_filename))) {
                printf("Usage: %s  input.obj output.obj [n_iters] [gausscurv_threshold]\n", argv[0]);
                return(1);
        }

        get_int_argument(250, &n_iters);
        get_real_argument(0.1, &gc_threshold);
    
        if(input_graphics_any_format(object_filename, &format, &n_objects, &objects) != OK)
                return(1);

        if(n_objects != 1 || get_object_type(objects[0]) != POLYGONS) {
                printf("File must contain 1 polygons object.\n");
                return(1);
        }

        polygons = get_polygons_ptr(objects[0]);

        gc_strength      = (double *)malloc(sizeof(double)*polygons->n_points);
            
        get_all_polygon_point_neighbours(polygons, &n_neighbours, &neighbours);

        get_polygon_vertex_curvatures_cg(polygons, n_neighbours, neighbours,
                                         0.0, 1, gc_strength);
                                         
        for (i=0; i<polygons->n_points; i++) {
                /* use absolute gaussian gc_strength if values are above threshold otherwise don't smooth */
                gc_strength[i] = (fabs(gc_strength[i]) > gc_threshold) ? fabs(gc_strength[i]) : 0.0;       
                /* get sure that values are <= 1.0 */
                gc_strength[i] = (gc_strength[i] > 1.0) ? 1.0 : gc_strength[i];       
        }
        
        /* smooth according to strength */
        gaussian_curv_smoothing(polygons, gc_strength, n_iters);
                        
        compute_polygon_normals(polygons);
        output_graphics_any_format(output_filename, format, 1, objects);

        delete_object_list(n_objects, objects);
        
        free(gc_strength);

        return(0);
}
