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
        char             *source_file, *source_sphere_file, *warped_sphere;
        char             *target_file, *target_sphere_file, *warped_file;
        File_formats     format;
        polygons_struct  *source, *source_sphere;
        polygons_struct  *target, *target_sphere, *warped;
        int              x, y, i, j, l, m;
        int              xy_size, n_objects, bandwidth, bandwidth2, cutoff;
        int              size_map[2], n_triangles, n_polygons;
        double           *srcoeffsx, *sicoeffsx, *trcoeffsx, *ticoeffsx;
        double           *srcoeffsy, *sicoeffsy, *trcoeffsy, *ticoeffsy;
        double           *srcoeffsz, *sicoeffsz, *trcoeffsz, *ticoeffsz;
        double           *rdatax0, *rdatay0, *rdataz0;
        double           *rdatax, *rdatay, *rdataz;
        double           H00, H01, H10, H11, valuex, valuey, valuez;
        double           u, v, x0, y0, xp, yp, xm, ym;
        object_struct    **objects, *object;
        int              dataformat;
        Point            centre, new_point;    

        initialize_argument_processing(argc, argv);

        if (!get_string_argument(NULL, &source_file) ||
            !get_string_argument(NULL, &source_sphere_file) ||
            !get_string_argument(NULL, &target_file) ||
            !get_string_argument(NULL, &target_sphere_file) ||    
            !get_string_argument(NULL, &warped_file) ||    
            !get_string_argument(NULL, &warped_sphere)) {
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

        rdatax0   = (double *) malloc(sizeof(double) * bandwidth2*bandwidth2);
        rdatay0   = (double *) malloc(sizeof(double) * bandwidth2*bandwidth2);
        rdataz0   = (double *) malloc(sizeof(double) * bandwidth2*bandwidth2);
        rdatax    = (double *) malloc(sizeof(double) * bandwidth2*bandwidth2);
        rdatay    = (double *) malloc(sizeof(double) * bandwidth2*bandwidth2);
        rdataz    = (double *) malloc(sizeof(double) * bandwidth2*bandwidth2);
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

        get_equally_sampled_coords_of_polygon(source, source_sphere, bandwidth,
                                              rdatax0, rdatay0, rdataz0);
        get_sph_coeffs_of_realdata(rdatax0, bandwidth, dataformat,
                                   srcoeffsx, sicoeffsx);
        get_sph_coeffs_of_realdata(rdatay0, bandwidth, dataformat,
                                   srcoeffsy, sicoeffsy);
        get_sph_coeffs_of_realdata(rdataz0, bandwidth, dataformat,
                                   srcoeffsz, sicoeffsz);
  
        get_equally_sampled_coords_of_polygon(target, target_sphere, bandwidth,
                                              rdatax, rdatay, rdataz);
        get_sph_coeffs_of_realdata(rdatax, bandwidth, dataformat,
                                   trcoeffsx, ticoeffsx);
        get_sph_coeffs_of_realdata(rdatay, bandwidth, dataformat,
                                   trcoeffsy, ticoeffsy);
        get_sph_coeffs_of_realdata(rdataz, bandwidth, dataformat,
                                   trcoeffsz, ticoeffsz);

        if (cutoff > 0 && cutoff < bandwidth) {
                butterworth_filter(bandwidth, cutoff, srcoeffsx, srcoeffsx);
                butterworth_filter(bandwidth, cutoff, srcoeffsy, srcoeffsy);
                butterworth_filter(bandwidth, cutoff, srcoeffsz, srcoeffsz);
                butterworth_filter(bandwidth, cutoff, sicoeffsx, sicoeffsx);
                butterworth_filter(bandwidth, cutoff, sicoeffsy, sicoeffsy);
                butterworth_filter(bandwidth, cutoff, sicoeffsz, sicoeffsz);
        }

        if (cutoff > 0 && cutoff < bandwidth) {
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

        get_realdata_from_sph_coeffs(rdatax, bandwidth, 1,
                                     srcoeffsx, sicoeffsx);
        get_realdata_from_sph_coeffs(rdatay, bandwidth, 1,
                                     srcoeffsy, sicoeffsy);
        get_realdata_from_sph_coeffs(rdataz, bandwidth, 1,
                                     srcoeffsz, sicoeffsz);
        
        object = create_object(POLYGONS);
        warped = get_polygons_ptr(object);
        
        fill_Point(centre, 0.0, 0.0, 0.0);
        create_tetrahedral_sphere(&centre, 1.0, 1.0, 1.0, n_triangles, warped);

        create_polygons_bintree(warped, round((double) warped->n_items *
                                              BINTREE_FACTOR));

        size_map[0] = bandwidth2;
        size_map[1] = bandwidth2;

        for (i = 0; i < warped->n_points; i++) {
                point_to_uv(&warped->points[i], &u, &v);

                /* interpolate points */
                x0 = (double) (bandwidth2 - 1) * u;
                y0 = (double) (bandwidth2 - 1) * v;
                x = (int) x0;
                y = (int) y0;
                xp = x0 - x;
                yp = y0 - y;
                xm = 1.0 - xp;
                ym = 1.0 - yp;
        
                H00 = rdatax[bound(x,   y,   size_map)];
                H01 = rdatax[bound(x,   y+1, size_map)];
                H10 = rdatax[bound(x+1, y,   size_map)];
                H11 = rdatax[bound(x+1, y+1, size_map)];

                valuex = ym * (xm * H00 + xp * H10) + 
                         yp * (xm * H01 + xp * H11);
 
                H00 = rdatay[bound(x,   y,   size_map)];
                H01 = rdatay[bound(x,   y+1, size_map)];
                H10 = rdatay[bound(x+1, y,   size_map)];
                H11 = rdatay[bound(x+1, y+1, size_map)];

                valuey = ym * (xm * H00 + xp * H10) + 
                         yp * (xm * H01 + xp * H11);

                H00 = rdataz[bound(x,   y,   size_map)];
                H01 = rdataz[bound(x,   y+1, size_map)];
                H10 = rdataz[bound(x+1, y,   size_map)];
                H11 = rdataz[bound(x+1, y+1, size_map)];

                valuez = ym * (xm * H00 + xp * H10) + 
                         yp * (xm * H01 + xp * H11);

                fill_Point(new_point, valuex, valuey, valuez);
                warped->points[i] = new_point;
        }

        compute_polygon_normals(warped);

        if (output_graphics_any_format(warped_sphere, ASCII_FORMAT,
                                       1, &object) != OK)
                exit(EXIT_FAILURE);

        for (i=0; i<bandwidth2*bandwidth2; i++) {
                rdatax[i] += rdatax0[i];
                rdatay[i] += rdatay0[i];
                rdataz[i] += rdataz0[i];        
        }

        create_tetrahedral_sphere(&centre, 1.0, 1.0, 1.0, n_triangles, warped);

        for (i = 0; i < warped->n_points; i++) {
                point_to_uv(&warped->points[i], &u, &v);

                /* interpolate points */
                x0 = (double) (bandwidth2 - 1) * u;
                y0 = (double) (bandwidth2 - 1) * v;
                x = (int) x0;
                y = (int) y0;
                xp = x0 - x;
                yp = y0 - y;
                xm = 1.0 - xp;
                ym = 1.0 - yp;
        
                H00 = rdatax[bound(x,   y,   size_map)];
                H01 = rdatax[bound(x,   y+1, size_map)];
                H10 = rdatax[bound(x+1, y,   size_map)];
                H11 = rdatax[bound(x+1, y+1, size_map)];

                valuex = ym * (xm * H00 + xp * H10) + 
                         yp * (xm * H01 + xp * H11);
 
                H00 = rdatay[bound(x,   y,   size_map)];
                H01 = rdatay[bound(x,   y+1, size_map)];
                H10 = rdatay[bound(x+1, y,   size_map)];
                H11 = rdatay[bound(x+1, y+1, size_map)];

                valuey = ym * (xm * H00 + xp * H10) + 
                         yp * (xm * H01 + xp * H11);

                H00 = rdataz[bound(x,   y,   size_map)];
                H01 = rdataz[bound(x,   y+1, size_map)];
                H10 = rdataz[bound(x+1, y,   size_map)];
                H11 = rdataz[bound(x+1, y+1, size_map)];

                valuez = ym * (xm * H00 + xp * H10) + 
                         yp * (xm * H01 + xp * H11);

                fill_Point(new_point, valuex, valuey, valuez);
                warped->points[i] = new_point;
        }

        compute_polygon_normals(warped);

        if (output_graphics_any_format(warped_file, ASCII_FORMAT,
                                       1, &object) != OK)
                exit(EXIT_FAILURE);


        delete_object_list(n_objects, objects);

        free(rdatax0); free(rdatay0); free(rdataz0);
        free(rdatax); free(rdatay); free(rdataz);
        free(srcoeffsx); free(srcoeffsy); free(srcoeffsz);
        free(sicoeffsx); free(sicoeffsy); free(sicoeffsz);
        free(trcoeffsx); free(trcoeffsy); free(trcoeffsz);
        free(ticoeffsx); free(ticoeffsy); free(ticoeffsz);
                
        return(EXIT_SUCCESS);
}
