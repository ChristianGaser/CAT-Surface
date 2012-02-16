/* Rachel Yotter - rachel.yotter@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 */

#include <volume_io/internal_volume_io.h>
#include <bicpl.h>
#include <bicpl/deform.h>
#include <ParseArgv.h>

#include "CAT_SheetIO.h"
#include "CAT_Map2d.h"
#include "CAT_Blur2d.h"
#include "CAT_SPH.h"
#include "CAT_Intersect.h"
#include "CAT_SurfaceIO.h"
#include "CAT_Defect.h"

#define DATAFORMAT 1 /* 1 = real data, 0 = complex data */
#define DEBUG 1
#define DUMP_FILES 0

#define FLAG_MODIFY 0
#define FLAG_PRESERVE 1

/* argument defaults */
int bw = 1024;
int lim = 64;
int n_triangles = 327680;
char *t1_file = NULL;
char *reparam_file = NULL;

Volume volume;

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
  {"-t1", ARGV_STRING, (char *) 1,
    (char *) &t1_file,
    "Optional T1-image for post-harmonic topology correction."},
  {"-sphere", ARGV_STRING, (char *) 1, (char *) &reparam_file,
     "Sphere object for reparameterization."},
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
add_neighbours(polygons_struct *surface, polygons_struct *lbw,
               int **neighbours, int *n_neighbours, int p, int *flag, int level)
{
        int n, idx;
        double dist;

        if (level == 0) {
                for (n = 0; n < n_neighbours[p]; n++) {
                        idx = neighbours[p][n];

                        if (flag[idx] == 1)
                                continue;

                        dist = distance_between_points(&surface->points[idx],
                                                       &lbw->points[idx]);
                        if (dist > 2.0) { /* prevent discontinuities */
                                flag[idx] = 1;
                                add_neighbours(surface, lbw, neighbours,
                                               n_neighbours, idx, flag, 0);
                        }
                }
                return;
        }

        for (n = 0; n < n_neighbours[p]; n++) {
                idx = neighbours[p][n];

                if (flag[idx] == 1)
                        continue;

                flag[idx] = 1;
                add_neighbours(surface, lbw, neighbours, n_neighbours, idx,
                               flag, level - 1);
        }
}

void
resample_defects_sph(polygons_struct *sphere, int *defects, int *polydefects,
                     int *remap_defects, int *remap_polydefects, int n_items)
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

        remap_defect(sphere, defects, polydefects, remap, remap_defects,
                     remap_polydefects);

        delete_object_list(1, objects);
}

void
sph_postcorrect(polygons_struct *surface, polygons_struct *sphere, int *defects,
                int *polydefects, int n_defects, int *holes,
                double t1_threshold, polygons_struct *hbw, polygons_struct *lbw)
{
        double *sharpness;
        int *n_neighbours, **neighbours;
        int p, d, val, *flag;
        int *hbw_defects, *hbw_polydefects, *hbw_holes;
        Point *pts;
        deform_struct deform;

        /* remap the defects onto the spherical harmonic reconstruction */
        create_polygon_point_neighbours(hbw, TRUE, &n_neighbours,
                                        &neighbours, NULL, NULL);

        if (DEBUG) fprintf(stderr,"resample_defects_sph...\n");
        hbw_defects = (int *) malloc(sizeof(int) * hbw->n_points);
        hbw_polydefects = (int *) malloc(sizeof(int) * hbw->n_items);
        resample_defects_sph(sphere, defects, polydefects,
                             hbw_defects, hbw_polydefects, hbw->n_items);
        expand_defects(hbw, hbw_defects, hbw_polydefects, 0, 2,
                       n_neighbours, neighbours);

        /* remap the hole flags */
        hbw_holes = (int *) malloc(sizeof(int) * hbw->n_points);
        for (d = 1; d <= n_defects; d++) {
                for (p = 0; p < sphere->n_points; p++) {
                        if (defects[p] == d && holes[p] != 0) {
                                val = holes[p];
                                break;
                        }
                }
 
                for (p = 0; p < hbw->n_points; p++) {
                        if (hbw_defects[p] == d)
                                hbw_holes[p] = val;
                }
        }

        /* optionally output the defects list */
        if (DUMP_FILES) {
                output_values_any_format("hbw_defects.txt", hbw->n_points,
                                         hbw_defects, TYPE_INTEGER);
        }

        /* preserve original points */
        pts = (Point *) malloc(sizeof(Point) * hbw->n_points);
        for (p = 0; p < hbw->n_points; p++)
                pts[p] = hbw->points[p];

        if (DEBUG) fprintf(stderr,"compute_local_sharpness...\n");
        sharpness = (double *) malloc(sizeof(double) * hbw->n_points);
        compute_local_sharpness(hbw, n_neighbours, neighbours, sharpness);

        if (DEBUG) fprintf(stderr,"patch_points...\n");
        flag = (int *) malloc(sizeof(int) * hbw->n_points);
        memset(flag, 0, sizeof(int) * hbw->n_points);
        for (p = 0; p < hbw->n_points; p++) {
                if (hbw_holes[p] == VENTRICLE || hbw_holes[p] == LARGE_DEFECT) {
                        flag[p] = 1;
                } else if (sharpness[p] > 60 && hbw_defects[p] != 0) {
                        flag[p] = 1; /* patch this one */
                        add_neighbours(hbw, lbw, neighbours, n_neighbours, p,
                                       flag, 1);
                }
        }

        /* combine the surfaces based on the modify flag */
        for (p = 0; p < hbw->n_points; p++) {
                if (flag[p] == 1)
                        hbw->points[p] = lbw->points[p];
        }

        /* find remaining self-intersections */
        fprintf(stderr,"Skip errornous find_selfintersections function\n");
//        n_defects = find_selfintersections(hbw, hbw_defects, hbw_polydefects);
        n_defects = join_intersections(hbw, hbw_defects, hbw_polydefects,
                                       n_neighbours, neighbours);

        fprintf(stderr,"%d self intersection(s) to repair\n", n_defects);

        /* patch self-intersections first */
        if (DEBUG) fprintf(stderr, "patch_selfintersections...\n");
        n_defects = patch_selfintersections(hbw, lbw, hbw_defects,
                                            hbw_polydefects, n_defects,
                                            n_neighbours, neighbours);
        fprintf(stderr,"Post-patch: %d self intersection(s) remaining\n", n_defects);

        /* smooth out remaining self-intersections */
        if (DEBUG) fprintf(stderr,"smooth_selfintersections...\n");
        n_defects = smooth_selfintersections(hbw, hbw_defects, hbw_polydefects,
                                             n_defects, n_neighbours,
                                             neighbours, 200);

        /* correct modified points using original T1 image */
        if (t1_file != NULL) {
                initialize_deformation_parameters(&deform);

                deform.fractional_step = 0.1;
                deform.max_step = 0.1;
                deform.max_search_distance = 15;
                deform.degrees_continuity = 0;
                deform.max_iterations = 40;
                deform.movement_threshold = 0.01;
                deform.stop_threshold = 0.0;

                deform.deform_data.type = VOLUME_DATA;

                deform.deform_data.volume = volume;
                deform.deform_data.label_volume = (Volume) NULL;

                if (add_deformation_model(&deform.deformation_model,
                                          -1, 0.5, "avg", -0.1, 0.1) != OK)
                        exit(EXIT_FAILURE);

                set_boundary_definition(&deform.boundary_definition,
                                        t1_threshold, t1_threshold,
                                        0, 0, 'n', 0);

                memset(flag, FLAG_MODIFY, sizeof(int) * hbw->n_points);
                for (p = 0; p < hbw->n_points; p++) {
                        if (hbw_holes[p] == VENTRICLE ||
                            (hbw_holes[p] == 0 &&
                             EQUAL_POINTS(hbw->points[p], pts[p])))
                                flag[p] = FLAG_PRESERVE;
                }

                if (DEBUG) fprintf(stderr,"deform_polygons_points...\n");
                deform_polygons_points(hbw, &deform, flag);

                n_defects = find_selfintersections(hbw, hbw_defects,
                                                   hbw_polydefects);

                fprintf(stderr,"%d intersection(s) to repair\n", n_defects);
                n_defects = join_intersections(hbw, hbw_defects,
                                               hbw_polydefects,
                                               n_neighbours, neighbours);
                fprintf(stderr,"%d intersection(s) to repair\n", n_defects);

                /* smooth out remaining self-intersections */
                if (DEBUG) fprintf(stderr,"smooth_selfintersections...\n");
                n_defects = smooth_selfintersections(hbw, hbw_defects,
                                                     hbw_polydefects,
                                                     n_defects, n_neighbours,
                                                     neighbours, 200);

        }

        if (DUMP_FILES) { 
                output_values_any_format("hbw_modpts.txt", hbw->n_points, flag,
                                 TYPE_INTEGER);
        }

        if (DEBUG) fprintf(stderr,"compute_polygon_normals...\n");
        compute_polygon_normals(hbw);

        /* clean up */ /*
        if (DEBUG) fprintf(stderr,"delete_polygon_point_neighbours...\n");
        delete_polygon_point_neighbours(hbw, n_neighbours,
                                        neighbours, NULL, NULL);
        */

        free(sharpness);
        free(flag);
        free(hbw_defects);
        free(hbw_polydefects);
        free(hbw_holes);
        free(pts);
}

object_struct **
fix_topology_sph(polygons_struct *surface, polygons_struct *sphere)
{
        object_struct **hbw_objects, **lbw_objects, **reparam_objects;
        polygons_struct *hbw, *lbw, *reparam;
        double *rcx, *icx, *rcy, *icy, *rcz, *icz;
        double *lrcx, *licx, *lrcy, *licy, *lrcz, *licz;
        double *rdatax, *rdatay, *rdataz;
        int bw2;
        int *defects, *polydefects, *holes, n_defects, n_objects, p;
        int *n_neighbours, **neighbours;
        double t1_threshold;
        File_formats format;

        /* find defects in original uncorrected surface */
        create_polygon_point_neighbours(sphere, TRUE, &n_neighbours,
                                        &neighbours, NULL, NULL);

        defects = (int *) malloc(sizeof(int) * sphere->n_points);
        polydefects = (int *) malloc(sizeof(int) * sphere->n_items);

        if (DEBUG) fprintf(stderr,"(orig) find_topological_defects...\n");
        n_defects = find_topological_defects(surface, sphere, defects,
                                             n_neighbours, neighbours);

        if (DEBUG) fprintf(stderr,"%d topological defects\n", n_defects);

        /* label defects as holes or handles */
        holes = (int *) malloc(sizeof(int) * sphere->n_points);
        if (t1_file != NULL) {
                t1_threshold = get_holes_handles(surface, sphere, defects,
                                                 n_defects, holes, volume,
                                                 n_neighbours, neighbours);

                if (DEBUG) fprintf(stderr,"T1 threshold = %f\n", t1_threshold);
                if (DUMP_FILES) { 
                        output_values_any_format("orig_holes.txt",
                                                 sphere->n_points, holes,
                                                 TYPE_INTEGER);
                }
        } else {
                for (p = 0; p < sphere->n_points; p++)
                        holes[p] = 2; /* always cut */
                t1_threshold = -1;
        }

        delete_polygon_point_neighbours(sphere, n_neighbours,
                                        neighbours, NULL, NULL);

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

        if (DEBUG) fprintf(stderr,"get_equally_sampled_coords_holes...\n");
        get_equally_sampled_coords_holes(surface, sphere, defects, n_defects,
                                         holes, bw, rdatax, rdatay, rdataz);

        if (DEBUG) fprintf(stderr,"get_sph_coeffs_of_realdata (hbw)...\n");
        get_sph_coeffs_of_realdata(rdatax, bw, DATAFORMAT, rcx, icx);
        get_sph_coeffs_of_realdata(rdatay, bw, DATAFORMAT, rcy, icy);
        get_sph_coeffs_of_realdata(rdataz, bw, DATAFORMAT, rcz, icz);

        hbw_objects = (object_struct **) malloc(sizeof(object_struct *));
        *hbw_objects = create_object(POLYGONS);
        hbw = get_polygons_ptr(*hbw_objects);

        if (reparam_file != NULL) {
                if (input_graphics_any_format(reparam_file, &format, &n_objects,
                                              &reparam_objects) != OK)
                        exit(EXIT_FAILURE);

                /* check that the surface file contains a polyhedron */
                if (n_objects != 1 ||
                    get_object_type(reparam_objects[0]) != POLYGONS) {
                        fprintf(stderr,"Reparam sphere file must contain 1 polygons object.\n");
                        exit(EXIT_FAILURE);
                }
                reparam = get_polygons_ptr(reparam_objects[0]);
                for (p = 0; p < reparam->n_points; p++)
                        set_vector_length(&reparam->points[p], 1.0);
                n_triangles = reparam->n_items;
        } else {
                reparam = NULL;
        }
        if (DEBUG) fprintf(stderr,"sample_sphere_from_sph (hbw)...\n");
        sample_sphere_from_sph(rdatax, rdatay, rdataz, hbw,
                               n_triangles, reparam, bw);

        if (DUMP_FILES) {
                output_graphics_any_format("hbw.obj", ASCII_FORMAT, 1,
                                           hbw_objects);
        }

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
                               lbw, n_triangles, reparam, bw);

        if (DUMP_FILES) {
                output_graphics_any_format("lbw.obj", ASCII_FORMAT, 1,
                                           lbw_objects);
        }

        free(rcx); free(rcy); free(rcz);
        free(icx); free(icy); free(icz);
        free(lrcx); free(lrcy); free(lrcz);
        free(licx); free(licy); free(licz);
        free(rdatax); free(rdatay); free(rdataz);

        sph_postcorrect(surface, sphere, defects, polydefects, n_defects, holes,
                        t1_threshold, hbw, lbw);

        delete_object_list(1, lbw_objects);
        free(defects);
        free(polydefects);
        free(holes);

        return hbw_objects;
}

int
main(int argc, char *argv[])
{
        char                 *surface_file, *sphere_file, *output_file;
        File_formats         format;
        int                  n_objects;
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

        if (t1_file != NULL) {
                if (input_volume_all(t1_file, 3, File_order_dimension_names,
                                     NC_UNSPECIFIED, FALSE, 0.0, 0.0, TRUE,
                                     &volume, NULL) == ERROR) {
                        fprintf(stderr, "Error opening T1 file: %s\n", t1_file);
                        exit(EXIT_FAILURE);
                }
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

        if (t1_file != NULL)
                delete_volume(volume);
    
        return(EXIT_SUCCESS);    
}
