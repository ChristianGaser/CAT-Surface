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



/* calc_polygon_areas(): calculate the polygon areas, store in areas */
void calc_polygon_areas(polygons_struct *polygons, double *areas) {
    Point tpoints[3];
    int p;

    for (p = 0; p < polygons->n_items; p++) {
        tpoints[0] = polygons->points[polygons->indices[
                            POINT_INDEX(polygons->end_indices,p,0)]];
        tpoints[1] = polygons->points[polygons->indices[
                            POINT_INDEX(polygons->end_indices,p,1)]];
        tpoints[2] = polygons->points[polygons->indices[
                            POINT_INDEX(polygons->end_indices,p,2)]];
        areas[p] = get_polygon_surface_area(3, tpoints);
    }
}

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
    int n_objects;
    File_formats format;
    int i, it, p;
    
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
    
    double iareas[ipolygons->n_items];
    calc_polygon_areas(ipolygons, iareas);

    double alpha = 2; /* how fast to change the point positions */
    float sphereRadius = get_largest_dist(polygons);
    float surfaceArea = get_polygons_surface_area(polygons);

    for (i = 0; i < 100; i++) { /* 5 iterations */
        double area_distortion = 0;
        for (it = 0; it < polygons->n_items; it++) {
            Point tpoints[3];

            for (p = 0; p < 3; p++) {
                tpoints[p] = polygons->points[polygons->indices[
                                    POINT_INDEX(polygons->end_indices,it,p)]];
            }

            double area = get_polygon_surface_area(3, tpoints);
            if (iareas[it] > 0 && area > 0) {
                area_distortion += fabs(log10(area/iareas[it]));
            }

            Point center;
            for (p = 0; p < 3; p++) {
                center.coords[p] = (tpoints[0].coords[p]
                                 + tpoints[1].coords[p]
                                 + tpoints[2].coords[p])/3.0;
            }

            double weight = alpha * (iareas[it] - area) / (area + iareas[it]); 
            for (p = 0; p < 3; p++) {
                Vector dir;
                Point np;

                SUB_POINTS(dir, tpoints[p], center);
                SCALE_VECTOR(dir, dir, 1 + weight / MAGNITUDE(dir));
                ADD_POINT_VECTOR(tpoints[p], center, dir);
                set_vector_length(&tpoints[p], sphereRadius);

                polygons->points[polygons->indices[
                                 POINT_INDEX(polygons->end_indices,it,p)]]
                                      = tpoints[p];
            }

        }
        alpha *= 0.95;
        convert_ellipsoid_to_sphere_with_surface_area(polygons, surfaceArea);
        printf("ad = %f\n", area_distortion);
    }

    compute_polygon_normals(polygons);

    output_graphics_any_format(ofname, format, 1, objects);
}


