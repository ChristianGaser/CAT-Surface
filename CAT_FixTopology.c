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

#define DATAFORMAT 1 /* 1 = real data, 0 = complex data */
#define DEBUG 1

/* argument defaults */
int bw = 1024;
int lim = 64;
int n_triangles = 327680;
BOOLEAN gauss_smooth = 0; /* gaussian curvature weighted smoothing on/off */

#define PINF  1.7976931348623157e+308 /* for doubles */
#define NINF -1.7976931348623157e+308 /* for doubles */

/* the argument table */
ArgvInfo argTable[] = {
  { "-bw", ARGV_INT, (char *) 1, 
    (char *) &bw,
    "Bandwidth of coefficients for spherical harmonic expansion." },
  { "-lim", ARGV_INT, (char *) 1, 
    (char *) &lim,
    "Limit bandwidth of spherical harmonic expansion." },
  { "-n", ARGV_INT, (char *) 1, 
    (char *) &n_triangles,
    "Number of triangles for sampled surface." },
  { "-smooth", ARGV_CONSTANT, (char *) TRUE, 
    (char *) &gauss_smooth,
    "Turn on Gaussian curvature weighted smoothing (recommended for central surface)." },
  { NULL, ARGV_END, NULL, NULL, NULL }
};

void
usage(char *executable)
{
        static char *usage_str = "\n\
Usage: %s surface.obj sphere.obj output.obj\n\
Correct the topology of a brain surface.\n\n\n";

       fprintf(stderr, usage_str, executable);
}

void
set_weights_of_neighbours(polygons_struct *surface, polygons_struct *lbw,
                          int **neighbours, int *n_neighbours, int p,
                          double *weights, int level)
{
        int n, idx;
        double dist;

        if (level == 0) {
                for (n = 0; n < n_neighbours[p]; n++) {
                        idx = neighbours[p][n];

                        if (weights[idx] == 1)
                                continue;

                        dist = distance_between_points(&surface->points[idx],
                                                       &lbw->points[idx]);
                        if (dist > 2.0) { /* prevent discontinuities */
                                weights[idx] = 1;
                                set_weights_of_neighbours(surface, lbw,
                                                          neighbours,
                                                          n_neighbours, idx,
                                                          weights, 0);
                        }
                }
                return;
        }

        for (n = 0; n < n_neighbours[p]; n++) {
                idx = neighbours[p][n];

                if (weights[idx] == 1)
                        continue;

                weights[idx] = 1;
                set_weights_of_neighbours(surface, lbw,
                                          neighbours, n_neighbours, idx,
                                          weights, level - 1);
        }
}

void
resample_defects_sph(polygons_struct *sphere, int *defects, int *remap_defects,
                     int n_items)
{
        object_struct **objects;
        polygons_struct *remap;
        Point centre;

        objects = (object_struct **) malloc(sizeof(object_struct *));
        *objects = create_object(POLYGONS);
        remap = get_polygons_ptr(*objects);
        initialize_polygons(remap, WHITE, NULL);

        fill_Point(centre, 0.0, 0.0, 0.0);
        create_tetrahedral_sphere(&centre, 1.0, 1.0, 1.0, n_items, remap);

        memset(remap_defects, 0, sizeof(int) * remap->n_points);

        remap_defect_points(sphere, defects, remap, remap_defects);

        delete_object_list(1, objects);
}


void
smooth_topology_sph(polygons_struct *surface, polygons_struct *sphere,
                    polygons_struct *hbw, polygons_struct *lbw)
{
        double *sharpness, *weights;
        int *n_neighbours, **neighbours;
        int i, p;
        FILE *fp;
        int n_defects, *defects, *hbw_defects;

        create_polygon_point_neighbours(sphere, TRUE, &n_neighbours,
                                        &neighbours, NULL, NULL);

        defects = (int *) malloc(sizeof(int) * sphere->n_points);
        hbw_defects = (int *) malloc(sizeof(int) * hbw->n_points);

        if (DEBUG) fprintf(stderr,"find_topological_defects...\n");
        n_defects = find_topological_defects(sphere, defects,
                                             n_neighbours, neighbours);
        if (n_defects == 0) return; /* nothing to be done! */

        delete_polygon_point_neighbours(sphere, n_neighbours,
                                        neighbours, NULL, NULL);

        create_polygon_point_neighbours(hbw, TRUE, &n_neighbours,
                                        &neighbours, NULL, NULL);

        if (DEBUG) fprintf(stderr,"resample_defects_sph...\n");
        resample_defects_sph(sphere, defects, hbw_defects, hbw->n_items);
        expand_defects(hbw, hbw_defects, 0, 1, n_neighbours, neighbours);

        if (open_file("defects.txt", WRITE_FILE, ASCII_FORMAT, &fp) != OK) {
                exit(0);
        }

        for (p = 0; p < hbw->n_points; p++)
                fprintf(fp, " %d.0\n", hbw_defects[p]);
        fclose(fp);

        sharpness = (double *) malloc(sizeof(double) * hbw->n_points);
        weights = (double *) malloc(sizeof(double) * hbw->n_points);
        memset(weights, 0.0, sizeof(double) * hbw->n_points);

        if (DEBUG) fprintf(stderr,"compute_local_sharpness...\n");
        compute_local_sharpness(hbw, n_neighbours, neighbours, sharpness);

        if (DEBUG) fprintf(stderr,"patch_points...\n");
        for (p = 0; p < hbw->n_points; p++) {
                if (sharpness[p] <= 60)
                        continue; /* skip, this point doesn't need patching */

                for (i = 0; i < n_defects; i++) {
                        if (hbw_defects[p] != 0) {
                                weights[p] = 1; /* patch this one */
                                set_weights_of_neighbours(hbw, lbw,
                                                          neighbours,
                                                          n_neighbours,
                                                          p, weights, 1);
                                break;
                        }
                }
        }

        /* combine the surfaces based on the weights */
        for (p = 0; p < hbw->n_points; p++) {
                for (i = 0; i < 3; i++) {
                        Point_coord(hbw->points[p], i) =
                            ((1 - weights[p]) * Point_coord(hbw->points[p],i)) +
                            ((    weights[p]) * Point_coord(lbw->points[p],i));
                }
        }
        fclose(fp);

        /* remove any self-intersections caused by the above step */
        if (DEBUG) fprintf(stderr,"repair_selfintersections...\n");
        repair_selfintersections(hbw, n_neighbours, neighbours);

        /* clean up */
        delete_polygon_point_neighbours(hbw, n_neighbours,
                                        neighbours, NULL, NULL);

        free(sharpness);
        free(weights);
        free(defects);
        free(hbw_defects);
}

object_struct **
fix_topology_sph(polygons_struct *surface, polygons_struct *sphere)
{
        object_struct **hbw_objects, **lbw_objects;
        polygons_struct *hbw, *lbw;
        double *rcx, *icx, *rcy, *icy, *rcz, *icz;
        double *lrcx, *licx, *lrcy, *licy, *lrcz, *licz;
        double *rdatax, *rdatay, *rdataz;
        int bw2;

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
        lrcx   = (double *) malloc(sizeof(double) * bw2);
        lrcy   = (double *) malloc(sizeof(double) * bw2);
        lrcz   = (double *) malloc(sizeof(double) * bw2);
        licx   = (double *) malloc(sizeof(double) * bw2);
        licy   = (double *) malloc(sizeof(double) * bw2);
        licz   = (double *) malloc(sizeof(double) * bw2);

        if (DEBUG) fprintf(stderr,"get_equally_sampled_coords_of_polygon...\n");
        get_equally_sampled_coords_of_polygon(surface, sphere, bw,
                                              rdatax, rdatay, rdataz);

        if (DEBUG) fprintf(stderr,"get_sph_coeffs_of_realdata (hbw)...\n");
        get_sph_coeffs_of_realdata(rdatax, bw, DATAFORMAT, rcx, icx);
        get_sph_coeffs_of_realdata(rdatay, bw, DATAFORMAT, rcy, icy);
        get_sph_coeffs_of_realdata(rdataz, bw, DATAFORMAT, rcz, icz);

        if (DEBUG) fprintf(stderr,"get_realdata_from_sph_coeffs (hbw)...\n");
        get_realdata_from_sph_coeffs(rdatax, bw, DATAFORMAT, rcx, icx);
        get_realdata_from_sph_coeffs(rdatay, bw, DATAFORMAT, rcy, icy);
        get_realdata_from_sph_coeffs(rdataz, bw, DATAFORMAT, rcz, icz);

        hbw_objects = (object_struct **) malloc(sizeof(object_struct *));
        *hbw_objects = create_object(POLYGONS);
        hbw = get_polygons_ptr(*hbw_objects);
        if (DEBUG) fprintf(stderr,"sample_sphere_from_sph (hbw)...\n");
        sample_sphere_from_sph(rdatax, rdatay, rdataz, hbw, n_triangles, bw);

        if (DEBUG) fprintf(stderr,"butterworth_filter...\n");
        butterworth_filter(bw, lim, rcx, lrcx);
        butterworth_filter(bw, lim, rcy, lrcy);
        butterworth_filter(bw, lim, rcz, lrcz);
        butterworth_filter(bw, lim, icx, licx);
        butterworth_filter(bw, lim, icy, licy);
        butterworth_filter(bw, lim, icz, licz);
    
        if (DEBUG) fprintf(stderr,"get_realdata_from_sph_coeffs (lbw)...\n");
        get_realdata_from_sph_coeffs(rdatax, bw, DATAFORMAT, lrcx, licx);
        get_realdata_from_sph_coeffs(rdatay, bw, DATAFORMAT, lrcy, licy);
        get_realdata_from_sph_coeffs(rdataz, bw, DATAFORMAT, lrcz, licz);

        lbw_objects = (object_struct **) malloc(sizeof(object_struct *));
        *lbw_objects = create_object(POLYGONS);
        lbw = get_polygons_ptr(*lbw_objects);
        if (DEBUG) fprintf(stderr,"sample_sphere_from_sph (lbw)...\n");
        sample_sphere_from_sph(rdatax, rdatay, rdataz,
                               lbw, n_triangles, bw);

        free(rcx); free(rcy); free(rcz);
        free(icx); free(icy); free(icz);
        free(lrcx); free(lrcy); free(lrcz);
        free(licx); free(licy); free(licz);
        free(rdatax); free(rdatay); free(rdataz);

        if (DEBUG) fprintf(stderr,"smooth_topology_sph...\n");
        smooth_topology_sph(surface, sphere, hbw, lbw);

        delete_object_list(1, lbw_objects);

        return hbw_objects;
}

int
main(int argc, char *argv[])
{
        char                 *surface_file, *sphere_file, *output_file;
        File_formats         format;
        int                  n_objects;
        Volume               volume;
        polygons_struct      *surface, *sphere;
        object_struct        **surf_objects, **sphere_objects, **objects;

        /* Call ParseArgv */
        if (ParseArgv(&argc, argv, argTable, 0) || argc != 4) {
                usage(argv[0]);
                fprintf(stderr, "       %s -help\n\n", argv[0]);
                exit(EXIT_FAILURE);
        }

        initialize_argument_processing(argc, argv);

        if (!get_string_argument(NULL, &surface_file) ||
            !get_string_argument(NULL, &sphere_file) ||
            !get_string_argument(NULL, &output_file)) {
                usage(argv[0]);
                exit(EXIT_FAILURE);
        }
     
        if (input_graphics_any_format(surface_file, &format,
                                      &n_objects, &surf_objects) != OK)
                exit(EXIT_FAILURE);

        /* check that the surface file contains a polyhedron */
        if (n_objects != 1 || get_object_type(surf_objects[0]) != POLYGONS) {
                fprintf(stderr,"Surface file must contain 1 polygons object.\n");
                exit(EXIT_FAILURE);
        }
        /* get a pointer to the surface */
        surface = get_polygons_ptr(surf_objects[0]);
    
        if (input_graphics_any_format(sphere_file, &format,
                                      &n_objects, &sphere_objects) != OK)
                exit(EXIT_FAILURE);

        /* check that the surface file contains a polyhedron */
        if (n_objects != 1 || get_object_type(sphere_objects[0]) != POLYGONS) {
                fprintf(stderr,"Surface file must contain 1 polygons object.\n");
                exit(EXIT_FAILURE);
        }
        /* get a pointer to the surface */
        sphere = get_polygons_ptr(sphere_objects[0]);

        /* check that surface and sphere are same size */
        if (surface->n_items != sphere->n_items) {
                fprintf(stderr,"Surface and sphere must have same size.\n");
                exit(EXIT_FAILURE);
        }
    
        objects = fix_topology_sph(surface, sphere);

        if (output_graphics_any_format(output_file, ASCII_FORMAT, 1,
                                       objects) != OK)
                exit(EXIT_FAILURE);

        /* clean up */

        delete_object_list(1, surf_objects);
        delete_object_list(1, sphere_objects);
        delete_object_list(1, objects);
    
        return(EXIT_SUCCESS);    
}
