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

#define BW 1024
#define LIM 64
#define N_TRIANGLES 327680
#define DATAFORMAT 1 /* 1 = real data, 0 = complex data */

/* argument defaults */
double t1_threshold = 0.0;
char *t1_file = NULL;
BOOLEAN gauss_smooth = 0; /* gaussian curvature weighted smoothing on/off */

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
  { "-smooth", ARGV_CONSTANT, (char *) TRUE, 
    (char *) &gauss_smooth,
    "Turn on Gaussian curvature weighted smoothing (recommended for CS)." },
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
                          double *sharpness, double *weights, int level)
{
        int n, idx;
        double w;

        if (level == 0)
                return;

        for (n = 0; n < n_neighbours[p]; n++) {
                idx = neighbours[p][n];

                if (weights[idx] == 1)
                        continue;

                weights[idx] = 1;
                set_weights_of_neighbours(neighbours, n_neighbours, idx,
                                          sharpness, weights, level - 1);
        }
}

void
smooth_topology_sph(polygons_struct *polygons, polygons_struct *smoothed_polys)
{
        double *sharpness, *weights;
        Point *smooth_pts;
        Real value;
        int *n_neighbours, **neighbours;
        int i, p;

        create_polygon_point_neighbours(polygons, TRUE, &n_neighbours,
                                        &neighbours, NULL, NULL);

        smooth_pts = (Point *) malloc(sizeof(Point) * polygons->n_points);
        sharpness = (double *) malloc(sizeof(double) * polygons->n_points);
        weights = (double *) malloc(sizeof(double) * polygons->n_points);
        memset(weights, 0.0, sizeof(double) * polygons->n_points);

        compute_local_sharpness(polygons, n_neighbours, neighbours, sharpness);

        for (p = 0; p < polygons->n_points; p++) {
                if (sharpness[p] > 110.0) {
                        weights[p] = 1;
                        set_weights_of_neighbours(neighbours, n_neighbours, p,
                                                  sharpness, weights, 2);
                }
        }

        /* combine the surfaces based on the weights */
        for (p = 0; p < polygons->n_points; p++) {
                for (i = 0; i < 3; i++) {
                        Point_coord(polygons->points[p], i) =
                            ((1 - weights[p]) * Point_coord(polygons->points[p],i)) +
                            ((    weights[p]) * Point_coord(smoothed_polys->points[p],i));
                }
        }

        /* remove any additional peaks caused by the above step */
        compute_local_sharpness(polygons, n_neighbours, neighbours, sharpness);

        for (p = 0; p < polygons->n_points; p++) {
                if (sharpness[p] > 110.0 && weights[p] == 0) {
                        for (i = 0; i < 3; i++) {
                                Point_coord(polygons->points[p], i) =
                                      Point_coord(smoothed_polys->points[p],i);
                        }
                }
        }

        /* Gaussian curvature smoothing */
        if (gauss_smooth) {
                get_polygon_vertex_curvatures_cg(polygons, n_neighbours, neighbours,
                                                 0.0, 1, weights);
                for (p = 0; p < polygons->n_points; p++)
                        weights[p] = (fabs(weights[p]) > 0.01)
                                     ? 5.0*fabs(weights[p]) : 0.0;

                for (p = 0; p < polygons->n_points; p++) {
                        /* smooth only if gaussian curv strength is > 0 */
                        if (weights[p] > 0.0) {
                                heatkernel_blur_points(polygons->n_points,
                                                       polygons->points, NULL,
                                                       n_neighbours[p],
                                                       neighbours[p],
                                                       p, weights[p],
                                                       &smooth_pts[p],
                                                       &value);
                        } else smooth_pts[p] = polygons->points[p];
                }
                for (p = 0; p < polygons->n_points; p++)
                        polygons->points[p] = smooth_pts[p];
        }

        free(smooth_pts);
        free(sharpness);
        free(weights);
}

object_struct *
fix_topology_sph(polygons_struct *polygons, polygons_struct *sphere)
{
        object_struct *objects, *lim_objects;
        polygons_struct *sph_polys, *lim_polys;
        double *rcx, *icx, *rcy, *icy, *rcz, *icz;
        double *lrcx, *licx, *lrcy, *licy, *lrcz, *licz;
        double *rdatax, *rdatay, *rdataz;
        int bw2;

        bw2 = BW * BW;

        rdatax = (double *) malloc(sizeof(double) * 4 * bw2);
        rdatay = (double *) malloc(sizeof(double) * 4 * bw2);
        rdataz = (double *) malloc(sizeof(double) * 4 * bw2);
        rcx    = (double *) malloc(sizeof(double) * bw2);
        rcy    = (double *) malloc(sizeof(double) * bw2);
        rcz    = (double *) malloc(sizeof(double) * bw2);
        icx    = (double *) malloc(sizeof(double) * bw2);
        icy    = (double *) malloc(sizeof(double) * bw2);
        icz    = (double *) malloc(sizeof(double) * bw2);
        lrcx   = (double *) malloc(sizeof(double) * bw2);
        lrcy   = (double *) malloc(sizeof(double) * bw2);
        lrcz   = (double *) malloc(sizeof(double) * bw2);
        licx   = (double *) malloc(sizeof(double) * bw2);
        licy   = (double *) malloc(sizeof(double) * bw2);
        licz   = (double *) malloc(sizeof(double) * bw2);

        get_equally_sampled_coords_of_polygon(polygons, sphere, BW,
                                              rdatax, rdatay, rdataz);

        get_sph_coeffs_of_realdata(rdatax, BW, DATAFORMAT, rcx, icx);
        get_sph_coeffs_of_realdata(rdatay, BW, DATAFORMAT, rcy, icy);
        get_sph_coeffs_of_realdata(rdataz, BW, DATAFORMAT, rcz, icz);

        get_realdata_from_sph_coeffs(rdatax, BW, DATAFORMAT, rcx, icx);
        get_realdata_from_sph_coeffs(rdatay, BW, DATAFORMAT, rcy, icy);
        get_realdata_from_sph_coeffs(rdataz, BW, DATAFORMAT, rcz, icz);

        objects = create_object(POLYGONS);
        sph_polys = get_polygons_ptr(objects);
        sample_sphere_from_sph(rdatax, rdatay, rdataz,
                               sph_polys, N_TRIANGLES, BW);

        butterworth_filter(BW, LIM, rcx, lrcx);
        butterworth_filter(BW, LIM, rcy, lrcy);
        butterworth_filter(BW, LIM, rcz, lrcz);
        butterworth_filter(BW, LIM, icx, licx);
        butterworth_filter(BW, LIM, icy, licy);
        butterworth_filter(BW, LIM, icz, licz);
    
        get_realdata_from_sph_coeffs(rdatax, BW, DATAFORMAT, lrcx, licx);
        get_realdata_from_sph_coeffs(rdatay, BW, DATAFORMAT, lrcy, licy);
        get_realdata_from_sph_coeffs(rdataz, BW, DATAFORMAT, lrcz, licz);

        lim_objects = create_object(POLYGONS);
        lim_polys = get_polygons_ptr(lim_objects);
        sample_sphere_from_sph(rdatax, rdatay, rdataz,
                               lim_polys, N_TRIANGLES, BW);

        free(rcx); free(rcy); free(rcz);
        free(icx); free(icy); free(icz);
        free(lrcx); free(lrcy); free(lrcz);
        free(licx); free(licy); free(licz);
        free(rdatax); free(rdatay); free(rdataz);

        smooth_topology_sph(sph_polys, lim_polys);
        return objects;
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
