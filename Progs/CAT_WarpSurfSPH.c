/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 *
 */

#include    <bicpl.h>
#include    <float.h>

#include    "CAT_Surf.h"
#include    "CAT_SPH.h"
    
void
usage(char *executable)
{
        static  char  *usage_str = "\n\
Usage: %s  source.obj source_sphere.obj target.obj target_sphere.obj warped.obj [cutoff]\n\
        Use spherical harmonic coefficients to warp a surface to a given template.\n\
\n\n";

        fprintf(stderr, usage_str, executable);
}
int
main(int argc, char *argv[])
{
        char             *source_file, *source_sphere_file, *warped_sphere_file;
        char             *target_file, *target_sphere_file, *warped_file;
        File_formats     format;
        polygons_struct  *source, *source_sphere, *warped_sphere;
        polygons_struct  *target, *target_sphere, *warped;
        int              x, y, i, j, l, m;
        int              xy_size, n_objects, bandwidth, bandwidth2, cutoff, n_points;
        int              size_map[2], n_triangles, n_polygons, poly;
        double           *srcoeffsx, *sicoeffsx, *trcoeffsx, *ticoeffsx;
        double           *srcoeffsy, *sicoeffsy, *trcoeffsy, *ticoeffsy;
        double           *srcoeffsz, *sicoeffsz, *trcoeffsz, *ticoeffsz;
        double           *srx, *sry, *srz;
        double           *trx, *try, *trz;
        double           H00, H01, H10, H11, valuex, valuey, valuez;
        double           u, v, u2, v2, xp, yp, xm, ym;
        object_struct    **objects, *object, *object2;
        int              dataformat;
        Point            sphere_point, new_point, closest;    

        initialize_argument_processing(argc, argv);

        if (!get_string_argument(NULL, &source_file) ||
            !get_string_argument(NULL, &source_sphere_file) ||
            !get_string_argument(NULL, &target_file) ||
            !get_string_argument(NULL, &target_sphere_file) ||    
            !get_string_argument(NULL, &warped_file) ||    
            !get_string_argument(NULL, &warped_sphere_file)) {
                usage(argv[0]);
                exit(EXIT_FAILURE);
        }

        get_int_argument(32, &cutoff);
        bandwidth = 256;
        
        bandwidth2  = bandwidth*2;
        n_triangles = 81920;

        if (input_graphics_any_format(target_file, &format,
                                      &n_objects, &objects) != OK)
                exit(EXIT_FAILURE);

        /* check that the surface file contains a polyhedron */
        if (n_objects != 1 || get_object_type(objects[0]) != POLYGONS) {
                print("Surface file must contain 1 polygons object.\n");
                exit(EXIT_FAILURE);
        }

        /* get a pointer to the surface */
        target = get_polygons_ptr(objects[0]);

        if (input_graphics_any_format(target_sphere_file, &format,
                                      &n_objects, &objects) != OK)
                exit(EXIT_FAILURE);
        
        /* check that the surface file contains a polyhedron */
        if (n_objects != 1 || get_object_type(objects[0]) != POLYGONS) {
                print("Surface file must contain 1 polygons object.\n");
                exit(EXIT_FAILURE);
        }

        /* get a pointer to the surface */
        target_sphere = get_polygons_ptr(objects[0]);

        if (input_graphics_any_format(source_file, &format,
                                      &n_objects, &objects) != OK)
                exit(EXIT_FAILURE);
                
        /* check that the surface file contains a polyhedron */
        if (n_objects != 1 || get_object_type(objects[0]) != POLYGONS) {
                print("Surface file must contain 1 polygons object.\n");
                exit(EXIT_FAILURE);
        }
        
        /* get a pointer to the surface */
        source = get_polygons_ptr(objects[0]);
    
        if (input_graphics_any_format(source_sphere_file, &format,
                                      &n_objects, &objects) != OK)
                exit(EXIT_FAILURE);
        
        /* check that the surface file contains a polyhedron */
        if (n_objects != 1 || get_object_type(objects[0]) != POLYGONS) {
                print("Surface file must contain 1 polygons object.\n");
                exit(EXIT_FAILURE);
        }
        
        /* get a pointer to the surface */
        source_sphere = get_polygons_ptr(objects[0]);

        /* allocate memory for source and target data */
        srx       = (double *) malloc(sizeof(double) * bandwidth2*bandwidth2);
        sry       = (double *) malloc(sizeof(double) * bandwidth2*bandwidth2);
        srz       = (double *) malloc(sizeof(double) * bandwidth2*bandwidth2);
        trx       = (double *) malloc(sizeof(double) * bandwidth2*bandwidth2);
        try       = (double *) malloc(sizeof(double) * bandwidth2*bandwidth2);
        trz       = (double *) malloc(sizeof(double) * bandwidth2*bandwidth2);
        
        /* allocate memory for source and target SPH coefficients */
        srcoeffsx = (double *) malloc(sizeof(double) * bandwidth*bandwidth);
        srcoeffsy = (double *) malloc(sizeof(double) * bandwidth*bandwidth);
        srcoeffsz = (double *) malloc(sizeof(double) * bandwidth*bandwidth);
        sicoeffsx = (double *) malloc(sizeof(double) * bandwidth*bandwidth);
        sicoeffsy = (double *) malloc(sizeof(double) * bandwidth*bandwidth);
        sicoeffsz = (double *) malloc(sizeof(double) * bandwidth*bandwidth);
        trcoeffsx = (double *) malloc(sizeof(double) * bandwidth*bandwidth);
        trcoeffsy = (double *) malloc(sizeof(double) * bandwidth*bandwidth);
        trcoeffsz = (double *) malloc(sizeof(double) * bandwidth*bandwidth);
        ticoeffsx = (double *) malloc(sizeof(double) * bandwidth*bandwidth);
        ticoeffsy = (double *) malloc(sizeof(double) * bandwidth*bandwidth);
        ticoeffsz = (double *) malloc(sizeof(double) * bandwidth*bandwidth);
        		
        /* dataformat indicates real data */
        dataformat = 1;

        /* get source SPH coefficients */
        get_equally_sampled_coords_of_polygon(source, source_sphere, bandwidth, srx, sry, srz);
        get_sph_coeffs_of_realdata(srx, bandwidth, dataformat, srcoeffsx, sicoeffsx);
        get_sph_coeffs_of_realdata(sry, bandwidth, dataformat, srcoeffsy, sicoeffsy);
        get_sph_coeffs_of_realdata(srz, bandwidth, dataformat, srcoeffsz, sicoeffsz);
  
        /* get target SPH coefficients */
        get_equally_sampled_coords_of_polygon(target, target_sphere, bandwidth, trx, try, trz);
        get_sph_coeffs_of_realdata(trx, bandwidth, dataformat, trcoeffsx, ticoeffsx);
        get_sph_coeffs_of_realdata(try, bandwidth, dataformat, trcoeffsy, ticoeffsy);
        get_sph_coeffs_of_realdata(trz, bandwidth, dataformat, trcoeffsz, ticoeffsz);

        /* filter source and target SPH coefficients */
        if (cutoff > 0 && cutoff < bandwidth) {
                butterworth_filter(bandwidth, cutoff, srcoeffsx, srcoeffsx);
                butterworth_filter(bandwidth, cutoff, srcoeffsy, srcoeffsy);
                butterworth_filter(bandwidth, cutoff, srcoeffsz, srcoeffsz);
                butterworth_filter(bandwidth, cutoff, sicoeffsx, sicoeffsx);
                butterworth_filter(bandwidth, cutoff, sicoeffsy, sicoeffsy);
                butterworth_filter(bandwidth, cutoff, sicoeffsz, sicoeffsz);
                butterworth_filter(bandwidth, cutoff, trcoeffsx, trcoeffsx);
                butterworth_filter(bandwidth, cutoff, trcoeffsy, trcoeffsy);
                butterworth_filter(bandwidth, cutoff, trcoeffsz, trcoeffsz);
                butterworth_filter(bandwidth, cutoff, ticoeffsx, ticoeffsx);
                butterworth_filter(bandwidth, cutoff, ticoeffsy, ticoeffsy);
                butterworth_filter(bandwidth, cutoff, ticoeffsz, ticoeffsz);
        }

        /* calculate difference of SPH coefficients */
        for (l = 0; l < bandwidth; l++) {
                for (m = -l; m < l+1; m++) {
                        i = seanindex(m, l, bandwidth);
                        srcoeffsx[i] = trcoeffsx[i] - srcoeffsx[i];
                        srcoeffsy[i] = trcoeffsy[i] - srcoeffsy[i];
                        srcoeffsz[i] = trcoeffsz[i] - srcoeffsz[i];
                        sicoeffsx[i] = ticoeffsx[i] - sicoeffsx[i];
                        sicoeffsy[i] = ticoeffsy[i] - sicoeffsy[i];
                        sicoeffsz[i] = ticoeffsz[i] - sicoeffsz[i];
                }
        }

        /* get real data from difference SPH coefficients */
        get_realdata_from_sph_coeffs(trx, bandwidth, dataformat, srcoeffsx, sicoeffsx);
        get_realdata_from_sph_coeffs(try, bandwidth, dataformat, srcoeffsy, sicoeffsy);
        get_realdata_from_sph_coeffs(trz, bandwidth, dataformat, srcoeffsz, sicoeffsz);
        

        /* add difference to source data */
        for (i=0; i<bandwidth2*bandwidth2; i++) {
                srx[i] += trx[i];
                sry[i] += try[i];
                srz[i] += trz[i];        
        }

        /* create sphere */
        object = create_object(POLYGONS);
        warped = get_polygons_ptr(object);
        copy_polygons(source_sphere,warped);
        
        object2 = create_object(POLYGONS);
        warped_sphere = get_polygons_ptr(object2);
        copy_polygons(source_sphere,warped_sphere);

        /* set radius to 1 */
        for (i = 0; i < source_sphere->n_points; i++) 
                set_vector_length(&source_sphere->points[i], 1.0);
        for (i = 0; i < target_sphere->n_points; i++) 
                set_vector_length(&target_sphere->points[i], 1.0);

        create_polygons_bintree(warped, round((double) warped->n_items *
                                              BINTREE_FACTOR));
        create_polygons_bintree(target, round((double) target->n_items *
                                              BINTREE_FACTOR));

        size_map[0] = bandwidth2;
        size_map[1] = bandwidth2;

        for (i = 0; i < source->n_points; i++) {

                fill_Point(sphere_point, Point_x(source_sphere->points[i]), Point_y(source_sphere->points[i]), Point_z(source_sphere->points[i]));
                point_to_uv(&sphere_point, &u, &v);

                /* interpolate points */
                xp = u*((double)bandwidth2) - 0.5;
                yp = v*((double)bandwidth2) - 0.5;

                x = (int) floor(xp); xp -= x; xm = 1.0 - xp;
                y = (int) floor(yp); yp -= y; ym = 1.0 - yp;
        
                H00 = srx[bound(x,   y,   size_map)];
                H01 = srx[bound(x,   y+1, size_map)];
                H10 = srx[bound(x+1, y,   size_map)];
                H11 = srx[bound(x+1, y+1, size_map)];

                valuex = ym * (xm * H00 + xp * H10) + 
                         yp * (xm * H01 + xp * H11);
 
                H00 = sry[bound(x,   y,   size_map)];
                H01 = sry[bound(x,   y+1, size_map)];
                H10 = sry[bound(x+1, y,   size_map)];
                H11 = sry[bound(x+1, y+1, size_map)];

                valuey = ym * (xm * H00 + xp * H10) + 
                         yp * (xm * H01 + xp * H11);

                H00 = srz[bound(x,   y,   size_map)];
                H01 = srz[bound(x,   y+1, size_map)];
                H10 = srz[bound(x+1, y,   size_map)];
                H11 = srz[bound(x+1, y+1, size_map)];

                valuez = ym * (xm * H00 + xp * H10) + 
                         yp * (xm * H01 + xp * H11);

                fill_Point(new_point, valuex, valuey, valuez);
                warped->points[i] = new_point;
}
        warped = get_polygons_ptr(object);
        copy_polygons(target_sphere,warped);
        warped_sphere = get_polygons_ptr(object2);
        copy_polygons(target_sphere,warped_sphere);

                poly = find_closest_polygon_point(&warped->points[i], target, &closest);
//                n_points = get_polygon_points(target_sphere, poly, &sphere_point);
                sphere_point = target_sphere->points[target_sphere->indices[POINT_INDEX(target_sphere->end_indices,poly,0)]];

                point_to_uv(&sphere_point, &u2, &v2);
                
//                fprintf(stderr,"%3.2f %3.2f %g\n",u,u2,u-u2);
                u = u2; v = v2;
                /* wrap borders */
                if (v < 0.0) {
                        v = -v;
                        u += 0.5;
                }
                if (v > 1.0) {
                        v = 2 - v;
                        u += 0.5;
                }
                while (u < 0.0)  u += 1.0;
                while (u >= 1.0) u -= 1.0;
                
                uv_to_point(u, v, &new_point);
                warped_sphere->points[i] = new_point;
                
        }


        if (output_graphics_any_format(warped_file, ASCII_FORMAT,
                                       1, &object) != OK)
                exit(EXIT_FAILURE);

        if (output_graphics_any_format(warped_sphere_file, ASCII_FORMAT,
                                       1, &object2) != OK)
                exit(EXIT_FAILURE);

        delete_object_list(n_objects, objects);

        free(srx); free(sry); free(srz);
        free(trx); free(try); free(trz);
        free(srcoeffsx); free(srcoeffsy); free(srcoeffsz);
        free(sicoeffsx); free(sicoeffsy); free(sicoeffsz);
        free(trcoeffsx); free(trcoeffsy); free(trcoeffsz);
        free(ticoeffsx); free(ticoeffsy); free(ticoeffsz);
               
        return(EXIT_SUCCESS);
}
