/* Rachel Yotter - rachel.yotter@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id: CAT_FixTopology.c 89 2009-01-27 14:43:59Z raytrace $
 */

#include <volume_io/internal_volume_io.h>
#include <bicpl.h>
#include <ParseArgv.h>

#include "CAT_SheetIO.h"
#include "CAT_Map2d.h"
#include "CAT_Blur2d.h"
#include "CAT_SPH.h"

#define BW1 128
#define BW2 1024
#define N_TRIANGLES 327680
#define DATAFORMAT 1 /* 1 = real data, 0 = complex data */

/* argument defaults */
double t1_threshold = 0.0;
char *t1_file   = NULL;

/* the argument table */
ArgvInfo argTable[] = {
  {"-t1", ARGV_STRING, (char *) 1, 
    (char *) &t1_file,
    "Optional T1-image to weight topology correction. If values in T1-image are above the threshold defined with -th, \n\
           then the radius of the sphere is decreased by 1% to resample the most outside points. This indicates a hole. \n\
           For T1-values below the threshold this area is assumed to be a handle and the sphere radius is increased by 1% \n\
           to resample the most inside points."},
  { "-th", ARGV_FLOAT, (char *) 1, 
    (char *) &t1_threshold,
    "Threshold between GM and CSF." },
  { NULL, ARGV_END, NULL, NULL, NULL }
};

void
usage(char *executable)
{
        static char *usage_str = "\n\
Usage: %s surface.obj sphere.obj output.obj\n\
Correct the topology of a brain surface.\n\\n\n";

       fprintf(stderr, usage_str, executable);
}

object_struct *
get_sph_object(polygons_struct *polygons, polygons_struct *sphere,
               int bw, int bw_lim)
{
        double *rcx, *icx, *rcy, *icy, *rcz, *icz;
        double *rdatax, *rdatay, *rdataz;
        polygons_struct *sph_polygons;
        object_struct *objects;
        int bw2, i;

        bw2 = bw * bw;

        rdatax = (double *) malloc(sizeof(double) * 4 * bw2);
        rdatay = (double *) malloc(sizeof(double) * 4 * bw2);
        rdataz = (double *) malloc(sizeof(double) * 4 * bw2);
        rcx    = (double *) malloc(sizeof(double) * bw2);
        rcy    = (double *) malloc(sizeof(double) * bw2);
        rcz    = (double *) malloc(sizeof(double) * bw2);
        icx    = (double *) malloc(sizeof(double) * bw2);
        icy    = (double *) malloc(sizeof(double) * bw2);
        icz    = (double *) malloc(sizeof(double) * bw2);

        get_equally_sampled_coords_of_polygon(polygons, sphere, bw,
                                              rdatax, rdatay, rdataz);

        get_sph_coeffs_of_realdata(rdatax, bw, DATAFORMAT, rcx, icx);
        get_sph_coeffs_of_realdata(rdatay, bw, DATAFORMAT, rcy, icy);
        get_sph_coeffs_of_realdata(rdataz, bw, DATAFORMAT, rcz, icz);

        if (bw_lim > 0 && bw_lim < bw) {
                fprintf(stderr,"Limit BW..");
                butterworth_filter(bw, bw_lim, rcx, rcx);
                butterworth_filter(bw, bw_lim, rcy, rcy);
                butterworth_filter(bw, bw_lim, rcz, rcz);
                butterworth_filter(bw, bw_lim, icx, icx);
                butterworth_filter(bw, bw_lim, icy, icy);
                butterworth_filter(bw, bw_lim, icz, icz);
        }
    
        get_realdata_from_sph_coeffs(rdatax, bw, DATAFORMAT, rcx, icx);
        get_realdata_from_sph_coeffs(rdatay, bw, DATAFORMAT, rcy, icy);
        get_realdata_from_sph_coeffs(rdataz, bw, DATAFORMAT, rcz, icz);

        objects = create_object(POLYGONS);
        sph_polygons = get_polygons_ptr(objects);

        sample_sphere_from_sph(rdatax, rdatay, rdataz,
                               sph_polygons, N_TRIANGLES, bw);

        free(rcx);
        free(rcy);
        free(rcz);
        free(icx);
        free(icy);
        free(icz);
        free(rdatax);
        free(rdatay);
        free(rdataz);

        return objects;
}

static int
compare(const void *a, const void *b)
{
        double aa = * (const double *) a;
        double bb = * (const double *) b;
        if (aa < bb)
                return -1;
        else if (aa > bb)
                return 1;
        else
                return 0;
}

void
set_weights_of_neighbours(int **neighbours, int *n_neighbours, int p,
                          double *sharpness, double *weights)
{
        int n, idx;
        double w;

        for (n = 0; n < n_neighbours[p]; n++) {
                idx = neighbours[p][n];
                if (weights[idx] == 1 || sharpness[idx] < 45)
                        continue;

                weights[idx] = 1;
                set_weights_of_neighbours(neighbours, n_neighbours, idx,
                                          sharpness, weights);
        }
}

object_struct *
fix_topology_sph(polygons_struct *polygons, polygons_struct *sphere)
{
        object_struct *bw1_objects, *bw2_objects;
        polygons_struct *bw1, *bw2;
        double *sharpness, *weights;
        int *n_neighbours, **neighbours;
        int i, p, count;

        bw1_objects = get_sph_object(polygons, sphere, BW1, BW1);
        bw1 = get_polygons_ptr(bw1_objects);
        bw2_objects = get_sph_object(polygons, sphere, BW2, BW2);
        bw2 = get_polygons_ptr(bw2_objects);

        create_polygon_point_neighbours(bw2, TRUE, &n_neighbours,
                                        &neighbours, NULL, NULL);

        sharpness = (double *) malloc(sizeof(double) * bw1->n_points);
        weights = (double *) malloc(sizeof(double) * bw1->n_points);
        memset(weights, 0.0, sizeof(double) * bw1->n_points);

        compute_local_sharpness(bw2, n_neighbours, neighbours, sharpness);

        for (p = 0; p < bw2->n_points; p++) {
                if (sharpness[p] < 70.0)
                        continue;

                weights[p] = 1;
                set_weights_of_neighbours(neighbours, n_neighbours, p,
                                          sharpness, weights);
        }

        /* combine the surfaces based on the weights */
        count = 0;
        for (p = 0; p < bw2->n_points; p++) {
                for (i = 0; i < 3; i++) {
                        Point_coord(bw2->points[p], i) =
                            ((1 - weights[p]) * Point_coord(bw2->points[p],i)) +
                            ((    weights[p]) * Point_coord(bw1->points[p],i));
                }
        }

        return bw2_objects;
}

int
main(int argc, char *argv[])
{
        char                 *surface_file, *sphere_file, *output_file;
        File_formats         format;
        int                  i, j, n_objects;
        double               r, value, avgt1, *t1value;
        Volume               volume;
        polygons_struct      *polygons, *sphere;
        object_struct        **poly_objects, **objects;

        /* Call ParseArgv */
        if (ParseArgv(&argc, argv, argTable, 0) || argc != 4) {
                usage(argv[0]);
                fprintf(stderr, "       %s -help\n\n", argv[0]);
                return(1);
        }

        initialize_argument_processing(argc, argv);

        if (!get_string_argument(NULL, &surface_file) ||
            !get_string_argument(NULL, &sphere_file) ||
            !get_string_argument(NULL, &output_file)) {
                usage(argv[0]);
                return(1);
        }
     
        if (input_graphics_any_format(surface_file, &format,
                                      &n_objects, &poly_objects) != OK)
                return(1);

        /* check that the surface file contains a polyhedron */
        if (n_objects != 1 || get_object_type(poly_objects[0]) != POLYGONS) {
                printf("Surface file must contain 1 polygons object.\n");
                return(1);
        }
        /* get a pointer to the surface */
        polygons = get_polygons_ptr(poly_objects[0]);
    
        if (input_graphics_any_format(sphere_file, &format,
                                      &n_objects, &objects) != OK)
                return(1);

        /* check that the surface file contains a polyhedron */
        if (n_objects != 1 || get_object_type(objects[0]) != POLYGONS) {
                printf("Surface file must contain 1 polygons object.\n");
                return(1);
        }
        /* get a pointer to the surface */
        sphere = get_polygons_ptr(objects[0]);

        /* check that surface and sphere are same size */
        if (polygons->n_items != sphere->n_items) {
                printf("Surface and sphere must have same size.\n");
                return(1);
        }
    
        /* An optional T1 image can be used for changing the radius of the
         * sphere to influence the resampling. We assume that values in the
         * T1 image below a given threshold (=CSF) indicate a handle. The
         * radius in these areas is increased by 1% to get sure that the
         * bottom is resampled which equals a cutting, becasue the most
         * inside points are resampled.  In contrast T1-values above the
         * threshold (=GM/WM) point to a hole and the radius of the sphere
         * is decrease by 1%. Thus, the resampling is done form the most
         * outside points, which equals a filling.
         */
        if (t1_file != NULL) {
                if (input_volume(t1_file, 3, File_order_dimension_names,
                                 NC_UNSPECIFIED, FALSE, 0.0, 0.0,
                                 TRUE, &volume, NULL) != OK)
                        return(1);
                t1value = (double *) malloc(sizeof(double) * sphere->n_points);
                
                avgt1 = 0.0;
                for (i = 0; i < sphere->n_points; i++) {
                        evaluate_volume_in_world(volume,
                                                 RPoint_x(sphere->points[i]),
                                                 RPoint_y(sphere->points[i]),
                                                 RPoint_z(sphere->points[i]),
                                                 0, FALSE, 0.0, &value, NULL,
                                                 NULL, NULL, NULL, NULL,
                                                 NULL, NULL, NULL, NULL);
                        t1value[i] = value;
                        avgt1 += value;
                }
                if (t1_threshold == 0.0)
                        avgt1 /= sphere->n_points;
                else
                        avgt1 = t1_threshold;

                for (i = 0; i < sphere->n_points; i++) {
                        if (t1value[i] < avgt1)
                                r = 1.01;
                        else
                                r = 0.99;
                        for (j = 0; j < 3; j++) 
                                Point_coord(sphere->points[i], j) *= r;                
                }
                free(t1value);
        }

        *objects = fix_topology_sph(polygons, sphere);

        if (output_graphics_any_format(output_file, ASCII_FORMAT, 1,
                                       objects) != OK)
                return(1);

        /* clean up */
    
        return(0);    
}
