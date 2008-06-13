/* Rachel Yotter - rachel.yotter@uni-jena.de                                 */
/* Department of Psychiatry                                                  */
/* University of Jena                                                        */
/*                                                                           */
/* Copyright Christian Gaser, University of Jena.                            */

/* Program to make a conformal map more isometric.                           */

#include <volume_io/internal_volume_io.h>
#include <bicpl.h>
#include <CAT_Curvature.h>
#include <CAT_Blur2d.h>
#include <CAT_Surf.h>
#include <CAT_SurfaceIO.h>



void usage(char *executable) {
    fprintf(stderr, 
      "\nUsage: %s [options] infile.obj conformalmap.obj outfile.obj\n\n\
         Adjust a conformal map to be more isometric (area preserving).\n\
         The original mesh is infile.obj and the conformal map of the mesh\n\
         is conformalmap.obj.  Results are saved in outfile.obj.\n\n",
                     executable);
}


int main(int argc, char** argv) {
    STRING ifname, cmfname, ofname;
    object_struct **objects;
    polygons_struct *ipolygons, *polygons;
    int     *n_neighbours, **neighbours;
    int n_objects;
    File_formats format;
    int i, n, it, p;
    Point pts[3];
    
    /* Get arguments */

    initialize_argument_processing(argc, argv);
    if (!get_string_argument(NULL, &ifname) ||
            !get_string_argument(NULL, &cmfname) ||
            !get_string_argument(NULL, &ofname)) {
        usage(argv[0]);
        return(1);
    }

    if (input_graphics_any_format(ifname, &format, &n_objects, &objects) != OK) {
        print("Error reading input file\n");
        return(1);
    }

    if (n_objects != 1 || get_object_type(objects[0]) != POLYGONS) {
        print("Input file must contain one polygon object.\n");
        return(1);
    }

    ipolygons = get_polygons_ptr(objects[0]);
    compute_polygon_normals(ipolygons);
    
    if (input_graphics_any_format(cmfname, &format, &n_objects, &objects) != OK) {
        print("Error reading conformal map file\n");
        return(1);
    }

    if (n_objects != 1 || get_object_type(objects[0]) != POLYGONS) {
        print("Conformal map file must contain one polygon object.\n");
        return(1);
    }

    polygons = get_polygons_ptr(objects[0]);
    compute_polygon_normals(polygons);

    if (polygons->n_points != ipolygons->n_points
            || polygons->n_items != ipolygons->n_items) {
        print("Input mesh and conformal map mesh do not match.\n");
        return(1);
    }
    
    float sphereRadius = get_largest_dist(polygons);
    
    create_polygon_point_neighbours(polygons, TRUE, &n_neighbours,
                                     &neighbours, NULL, NULL);
    
    float beta = 0.5; /* 1 = all area, 0 = all angle */

    for (it = 0; it < 100; it++) {
        float area_distortion = 0.0;

        for (p = 0; p < polygons->n_points; p++) {
            if (n_neighbours[p] > 1) { 
                
                float tileAreas[n_neighbours[p]];
                float itileAreas[n_neighbours[p]];
                float tileAngles[n_neighbours[p]];
                float itileAngles[n_neighbours[p]];
                float tileCenters[n_neighbours[p]*3];
                float totalArea = 0.0;
                float itotalArea = 0.0;
                    
                // Get 2 consecutive neighbors of this node
                for (n = 0; n < n_neighbours[p]; n++) {        
                    int n1 = neighbours[p][n];
                    int next = (n + 1) % n_neighbours[p];
                    int n2 = neighbours[p][next];
                    
                    // area and angle of the triangle
                    pts[0] = polygons->points[p];
                    pts[1] = polygons->points[n1];
                    pts[2] = polygons->points[n2];
                    tileAreas[n] = get_polygon_surface_area(3, pts);
                    totalArea += tileAreas[n];
                    tileAngles[n] = get_angle_between_points(&pts[1], &pts[0], &pts[2]);

                    // area and angles of the original triangle
                    pts[0] = ipolygons->points[p];
                    pts[1] = ipolygons->points[n1];
                    pts[2] = ipolygons->points[n2];
                    itileAreas[n] = get_polygon_surface_area(3, pts);
                    itotalArea += itileAreas[n];
                    itileAngles[n] = get_angle_between_points(&pts[1], &pts[0], &pts[2]);


                    area_distortion += fabs(itotalArea - totalArea)/n_neighbours[p];

                    // Save center of this tile
                    for (i = 0; i < 3; i++) {
                        tileCenters[n*3+i] = (Point_coord(polygons->points[p],i) + 
                                              Point_coord(polygons->points[n1],i) +
                                              Point_coord(polygons->points[n2],i)) / 3.0;
                    }
                }
                    
                // Compute the influence of the neighboring nodes
                float xyz[3] = {0.0, 0.0, 0.0};
                for (n = 0; n <  n_neighbours[p]; n++) {
                    if (itileAreas[n] > 0.0) {
                        float weight = (tileAreas[n] / totalArea);
                        /* weight +=  itileAngles[n]/(2*PI);
                        if (itotalArea != totalArea) {
                            weight *= 0.8;
                            weight += 0.2*(itileAreas[n] - tileAreas[n])/(itotalArea - totalArea);
                        } */
                        for (i = 0; i < 3; i++) {
                            xyz[i] += weight * tileCenters[n*3+i];
                        }
                    }
                }
                // Update the nodes position
                for (i = 0; i < 3; i++) {
                    Point_coord(polygons->points[p],i) = xyz[i];
                }
            }
        }
    }

    convert_ellipsoid_to_sphere_with_surface_area(polygons, 100.0);

    compute_polygon_normals(polygons);
    output_graphics_any_format(ofname, format, 1, objects);

    delete_object_list(n_objects, objects);

    return(0);

}


