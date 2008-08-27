/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 */

#include    <volume_io/internal_volume_io.h>
#include    <bicpl.h>
#include    <float.h>

#include    "CAT_Surf.h"
#include    "CAT_SPH.h"
    
void
usage(char *executable)
{
        static  char  *usage_str = "\n\
Usage: %s  source.obj source_sphere.obj target.obj target_sphere.obj warped.obj\n\
        Use spherical harmonic coefficients to warp a surface to a given template.\n\
\n\n";

        fprintf(stderr, usage_str, executable);
}
int
main(int argc, char *argv[])
{
        STRING               source_filename, source_sphere_filename, target_filename, target_sphere_filename, warped_filename;
        File_formats         format;
        polygons_struct      *polygons_source, *polygons_target, *polygons_source_sphere;
        polygons_struct      *polygons_target_sphere, *polygons_warped;
        int                  x, y, i, j, xy_size, n_objects, bandwidth, bandwidth2;
        int                  size_map[2], n_triangles, n_polygons;
        double               *srcoeffsx, *sicoeffsx, *srcoeffsy, *sicoeffsy, *srcoeffsz, *sicoeffsz;
        double               *trcoeffsx, *ticoeffsx, *trcoeffsy, *ticoeffsy, *trcoeffsz, *ticoeffsz;
        double               *rdatax, *rdatay, *rdataz, u, v;
        double               H00, H01, H10, H11, valuex, valuey, valuez;
        double               x0, y0, xp, yp, xm, ym;
        object_struct        **objects, *object;
        int                  dataformat, cutoff;
        Point                centre, new_point;    

        bandwidth   = 256;
        bandwidth2  = bandwidth*2;
        n_triangles = 81920;

        initialize_argument_processing(argc, argv);

        if (!get_string_argument(NULL, &source_filename) ||
            !get_string_argument(NULL, &source_sphere_filename) ||
            !get_string_argument(NULL, &target_filename) ||
            !get_string_argument(NULL, &target_sphere_filename) ||    
            !get_string_argument(NULL, &warped_filename)) {
                usage(argv[0]);
                return(1);
        }
        
        if (input_graphics_any_format(target_filename, &format, &n_objects, &objects) != OK)
                return(1);
        /* check that the surface file contains a polyhedron */
        if (n_objects != 1 || get_object_type(objects[0]) != POLYGONS) {
                print("Surface file must contain 1 polygons object.\n");
                return(1);
        }

        /* get a pointer to the surface */
        polygons_target = get_polygons_ptr(objects[0]);

        if (input_graphics_any_format(target_sphere_filename, &format, &n_objects, &objects) != OK)
                return(1);
        
        /* check that the surface file contains a polyhedron */
        if (n_objects != 1 || get_object_type(objects[0]) != POLYGONS) {
                print("Surface file must contain 1 polygons object.\n");
                return(1);
        }

        /* get a pointer to the surface */
        polygons_target_sphere = get_polygons_ptr(objects[0]);

        if (input_graphics_any_format(source_filename, &format, &n_objects, &objects) != OK)
                return(1);
                
        /* check that the surface file contains a polyhedron */
        if (n_objects != 1 || get_object_type(objects[0]) != POLYGONS) {
                print("Surface file must contain 1 polygons object.\n");
                return(1);
        }
        
        /* get a pointer to the surface */
        polygons_source = get_polygons_ptr(objects[0]);
    
        if (input_graphics_any_format(source_sphere_filename, &format, &n_objects, &objects) != OK)
                return(1);
        
        /* check that the surface file contains a polyhedron */
        if (n_objects != 1 || get_object_type(objects[0]) != POLYGONS) {
                print("Surface file must contain 1 polygons object.\n");
                return(1);
        }
        
        /* get a pointer to the surface */
        polygons_source_sphere = get_polygons_ptr(objects[0]);

        rdatax    = (double *) malloc(sizeof(double)*bandwidth2*bandwidth2);
        rdatay    = (double *) malloc(sizeof(double)*bandwidth2*bandwidth2);
        rdataz    = (double *) malloc(sizeof(double)*bandwidth2*bandwidth2);
        srcoeffsx = (double *) malloc(sizeof(double)*bandwidth*bandwidth);
        srcoeffsy = (double *) malloc(sizeof(double)*bandwidth*bandwidth);
        srcoeffsz = (double *) malloc(sizeof(double)*bandwidth*bandwidth);
        sicoeffsx = (double *) malloc(sizeof(double)*bandwidth*bandwidth);
        sicoeffsy = (double *) malloc(sizeof(double)*bandwidth*bandwidth);
        sicoeffsz = (double *) malloc(sizeof(double)*bandwidth*bandwidth);
        trcoeffsx = (double *) malloc(sizeof(double)*bandwidth*bandwidth);
        trcoeffsy = (double *) malloc(sizeof(double)*bandwidth*bandwidth);
        trcoeffsz = (double *) malloc(sizeof(double)*bandwidth*bandwidth);
        ticoeffsx = (double *) malloc(sizeof(double)*bandwidth*bandwidth);
        ticoeffsy = (double *) malloc(sizeof(double)*bandwidth*bandwidth);
        ticoeffsz = (double *) malloc(sizeof(double)*bandwidth*bandwidth);
        		
        get_equally_sampled_coords_of_polygon(polygons_source, polygons_source_sphere, bandwidth, rdatax, rdatay, rdataz);

        // dataformat indicates real data
        dataformat = 1;
        get_sph_coeffs_of_realdata(rdatax, bandwidth, dataformat, srcoeffsx, sicoeffsx);
        get_sph_coeffs_of_realdata(rdatay, bandwidth, dataformat, srcoeffsy, sicoeffsy);
        get_sph_coeffs_of_realdata(rdataz, bandwidth, dataformat, srcoeffsz, sicoeffsz);
  
        get_equally_sampled_coords_of_polygon(polygons_target, polygons_target_sphere, bandwidth, rdatax, rdatay, rdataz);
        get_sph_coeffs_of_realdata(rdatax, bandwidth, dataformat, trcoeffsx, ticoeffsx);
        get_sph_coeffs_of_realdata(rdatay, bandwidth, dataformat, trcoeffsy, ticoeffsy);
        get_sph_coeffs_of_realdata(rdataz, bandwidth, dataformat, trcoeffsz, ticoeffsz);

        cutoff = 4;
        replaceSPH(bandwidth, cutoff, srcoeffsx, trcoeffsx);
        replaceSPH(bandwidth, cutoff, srcoeffsy, trcoeffsy);
        replaceSPH(bandwidth, cutoff, srcoeffsz, trcoeffsz);
        replaceSPH(bandwidth, cutoff, sicoeffsx, ticoeffsx);
        replaceSPH(bandwidth, cutoff, sicoeffsy, ticoeffsy);
        replaceSPH(bandwidth, cutoff, sicoeffsz, ticoeffsz);

        get_realdata_from_sph_coeffs(rdatax, bandwidth, 1, srcoeffsx, sicoeffsx);
        get_realdata_from_sph_coeffs(rdatay, bandwidth, 1, srcoeffsy, sicoeffsy);
        get_realdata_from_sph_coeffs(rdataz, bandwidth, 1, srcoeffsz, sicoeffsz);
        
        object = create_object(POLYGONS);
        polygons_warped = get_polygons_ptr(object);
        
        fill_Point(centre, 0.0, 0.0, 0.0);
        create_tetrahedral_sphere(&centre, 1.0, 1.0, 1.0,
                n_triangles, polygons_warped);

        create_polygons_bintree(polygons_warped,
                round((double) polygons_warped->n_items * BINTREE_FACTOR));

        size_map[0] = bandwidth2;
        size_map[1] = bandwidth2;

        for (i=0; i<polygons_warped->n_points; i++) {

                point_to_uv(&polygons_warped->points[i], &u, &v);

                /* interpolate points */
                x0 = (double)(bandwidth2 - 1) * u;
                y0 = (double)(bandwidth2 - 1) * v;
                x = (int) x0;
                y = (int) y0;
                xp = x0 - x;
                yp = y0 - y;
                xm = 1.0 - xp;
                ym = 1.0 - yp;
        
                H00 = rdatax[bound(x,  y,  size_map)];
                H01 = rdatax[bound(x,  y+1,size_map)];
                H10 = rdatax[bound(x+1,y,  size_map)];
                H11 = rdatax[bound(x+1,y+1,size_map)];

                valuex = (ym * (xm * H00 + xp * H10) + 
                        yp * (xm * H01 + xp * H11));
 
                H00 = rdatay[bound(x,  y,  size_map)];
                H01 = rdatay[bound(x,  y+1,size_map)];
                H10 = rdatay[bound(x+1,y,  size_map)];
                H11 = rdatay[bound(x+1,y+1,size_map)];

                valuey = (ym * (xm * H00 + xp * H10) + 
                        yp * (xm * H01 + xp * H11));

                H00 = rdataz[bound(x,  y,  size_map)];
                H01 = rdataz[bound(x,  y+1,size_map)];
                H10 = rdataz[bound(x+1,y,  size_map)];
                H11 = rdataz[bound(x+1,y+1,size_map)];

                valuez = (ym * (xm * H00 + xp * H10) + 
                        yp * (xm * H01 + xp * H11));

                fill_Point(new_point,valuex,valuey,valuez);
                polygons_warped->points[i] = new_point;
        }

        compute_polygon_normals(polygons_warped);

        if (output_graphics_any_format(warped_filename, ASCII_FORMAT, 1, &object) != OK)
                return(1);


        delete_object_list(n_objects, objects);
        free(rdatax);
        free(rdatay);
        free(rdataz);
        free(srcoeffsx);
        free(srcoeffsy);
        free(srcoeffsz);
        free(sicoeffsx);
        free(sicoeffsy);
        free(sicoeffsz);
        free(trcoeffsx);
        free(trcoeffsy);
        free(trcoeffsz);
        free(ticoeffsx);
        free(ticoeffsy);
        free(ticoeffsz);
                
        return(0);
}
