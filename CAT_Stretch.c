/* Rachel Yotter - rachel.yotter@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 */

/* Program to make a conformal map more isometric. */

#include <volume_io/internal_volume_io.h>
#include <volume_io/geometry.h>
#include <bicpl.h>
#include <ParseArgv.h>

#include "CAT_Curvature.h"
#include "CAT_Blur2d.h"
#include "CAT_Surf.h"
#include "CAT_SurfaceIO.h"

#define PINF  1.7976931348623157e+308 /* for doubles */
#define NINF -1.7976931348623157e+308 /* for doubles */

int Iter = 100; /* number of iterations */
BOOLEAN Select = 0; /* by default, apply to all triangles */
double Threshold = 0; /* by default, apply to all triangles */
BOOLEAN LargeOnly = 0; /* only correct larger-than-orig triangles */
BOOLEAN Quiet = 0; /* turn progress reports on and off */

struct ptinfo {
        double *lengths;
        Vector *norm;
        double *areas;
};

static ArgvInfo argTable[] = {
  {"-iter", ARGV_INT, (char *) 1, (char *) &Iter,
       "number of iterations."},
  {"-select", ARGV_CONSTANT, (char *) TRUE, (char *) &Select,
    "apply stretch algorithm only if it lowers local area distortion" },
  {"-large", ARGV_CONSTANT, (char *) TRUE, (char *) &LargeOnly,
    "apply stretch algorithm only to triangles with positive area distortion" },
  {"-quiet", ARGV_CONSTANT, (char *) TRUE, (char *) &Quiet,
    "turn off progress reports" },
  {"-thresh", ARGV_FLOAT, (char *) 1, (char *) &Threshold,
       "apply stretch algorithm only if local area distortion is above threshold."},
   {NULL, ARGV_END, NULL, NULL, NULL}
};

void
usage(char *executable)
{
        char *usage_str =
"\nUsage: %s [options] infile.obj conformalmap.obj outfile.obj\n\n\
    Adjust a conformal map to be more isometric (area preserving) using\n\
    the stretch correction algorithm.  The original mesh is infile.obj and\n\
    the conformal map of the mesh is conformalmap.obj.  Results are saved\n\
    in outfile.obj.\n\n";

        fprintf(stderr, usage_str, executable);
}

int
main(int argc, char** argv)
{
        char               *input_file, *cmap_file, *output_file;
        object_struct      **objects;
        polygons_struct    *ipolygons, *polygons;
        int                *n_neighbours, **neighbours, n_objects;
        File_formats       format;
        int                i, n, it, p;
        Point              pts[3], npt, dir, newcenter;
        Vector             norm;
        double             dist, radius;
        double             adistort, newdistort, iarea, localarea[128];
        int                n1, n2, nidx, count, ok;
        struct ptinfo      **ipolyinfo;
        progress_struct    progress;

        if (ParseArgv(&argc, argv, argTable, 0) || (argc < 3)) {
                usage(argv[0]);
                fprintf(stderr, "       %s -help\n\n", argv[0]);
                return(1);
        }

        initialize_argument_processing(argc, argv);
        if (!get_string_argument(NULL, &input_file) ||
            !get_string_argument(NULL, &cmap_file) ||
            !get_string_argument(NULL, &output_file)) {
                fprintf(stderr, "\nUsage: %s [options] infile.obj conformalmap.obj outfile.obj\n", argv[0]);
                return(1);
        }

        if (input_graphics_any_format(input_file, &format, &n_objects,
                                      &objects) != OK) {
                printf("Error reading input file\n");
                return(1);
        }

        if (n_objects != 1 || get_object_type(objects[0]) != POLYGONS) {
                printf("Input file must contain one polygon object.\n");
                return(1);
        }

        ipolygons = get_polygons_ptr(objects[0]);
        compute_polygon_normals(ipolygons);

        if (input_graphics_any_format(cmap_file, &format, &n_objects,
                                      &objects) != OK) {
                printf("Error reading conformal map file\n");
                return(1);
        }

        if (n_objects != 1 || get_object_type(objects[0]) != POLYGONS) {
                printf("Conformal map file must contain one polygon object.\n");
                return(1);
        }

        polygons = get_polygons_ptr(objects[0]);
        compute_polygon_normals(polygons);

        if (polygons->n_points != ipolygons->n_points
            || polygons->n_items != ipolygons->n_items) {
                printf("Input mesh and conformal map mesh do not match.\n");
                return(1);
        }
    
        create_polygon_point_neighbours(polygons, TRUE, &n_neighbours,
                                         &neighbours, NULL, NULL);

        /* map the spherical map to a sphere of twice the surface area */
        radius = sqrt(get_polygons_surface_area(ipolygons) / (2.0 * PI));
        for (p = 0; p < polygons->n_points; p++)
                set_vector_length(&polygons->points[p], radius);


        /* first get some info about the original mesh */
        ipolyinfo = (struct ptinfo **)
                    malloc(sizeof(struct ptinfo *) * ipolygons->n_points);

        for (p = 0; p < ipolygons->n_points; p++) {
                if (n_neighbours[p] <= 1)
                        continue; /* skip this point */

                ipolyinfo[p] = (struct ptinfo *) malloc(sizeof(struct ptinfo));

                ipolyinfo[p]->areas = (double *) malloc(sizeof(double) *
                                                        n_neighbours[p]);

                ipolyinfo[p]->lengths = (double *) malloc(sizeof(double) *
                                                          n_neighbours[p]);

                ipolyinfo[p]->norm = (Vector *) malloc(sizeof(Vector) *
                                                       n_neighbours[p]);


                /* Get 2 consecutive neighbors of this node */
                for (n = 0; n < n_neighbours[p]; n++) {
                        n1 = neighbours[p][n];
                        n2 = neighbours[p][(n + 1) % n_neighbours[p]];

                        /* area of the triangle */
                        pts[0] = ipolygons->points[p];
                        pts[1] = ipolygons->points[n1];
                        pts[2] = ipolygons->points[n2];
                        ipolyinfo[p]->areas[n] = get_polygon_surface_area(3,
                                                                       pts) * 2;

                        /* normal of the triangle */
                        pts[0] = polygons->points[p];
                        pts[1] = polygons->points[n1];
                        pts[2] = polygons->points[n2];
                        find_polygon_normal(3, pts, &ipolyinfo[p]->norm[n]);

                        /* lengths to neighbouring nodes */
                        SUB_POINTS(dir, ipolygons->points[n1],
                                        ipolygons->points[p]);
                        ipolyinfo[p]->lengths[n] = MAGNITUDE(dir);
                }
        }
        if (Quiet == 0)
                initialize_progress_report(&progress, FALSE, Iter, "Stretch");

for (it = 1; it <= Iter; it++) {
        count = 0;

        for (p = 0; p < polygons->n_points; p++) {
                if (n_neighbours[p] <= 1)
                        continue; /* skip this point */

                if (Threshold > 0 || Select == 1 || LargeOnly == 1) {
                        adistort = 0.0;

                        pts[0] = polygons->points[p];
                        for (n = 0; n < n_neighbours[p]; n++) {
                                n1 = neighbours[p][n];
                                n2 = neighbours[p][(n + 1) % n_neighbours[p]];

                                /* area of the triangle */
                                pts[1] = polygons->points[n1];
                                pts[2] = polygons->points[n2];
                                localarea[n] = get_polygon_surface_area(3, pts);

                                iarea = ipolyinfo[p]->areas[n];
                                if (iarea > 0)
                                        adistort += log10(localarea[n]/iarea);
                        }
                        if (LargeOnly == 1) {
                                if (adistort + 1 < Threshold)
                                        continue; /* skip this point */
                                adistort = fabs(adistort / n_neighbours[p]);
                        } else {
                                adistort = fabs(adistort / n_neighbours[p]);
                                if (adistort < Threshold)
                                        continue; /* skip this point */
                       }
                }

                /* Algorithm - Stretch Optimization */
                fill_Point(newcenter, 0.0, 0.0, 0.0);

                for (n = 0; n <  n_neighbours[p]; n++) {
                        nidx = neighbours[p][n];
                        npt = polygons->points[nidx];

                        SUB_POINTS(dir, npt, polygons->points[p]);

                        dist = MAGNITUDE(dir) - ipolyinfo[p]->lengths[n];
                        dist /= MAGNITUDE(dir);
                        SCALE_VECTOR(dir, dir, dist);
                        ADD_POINT_VECTOR(newcenter, newcenter, dir);
                }

                Point_x(newcenter) /= n_neighbours[p];
                Point_y(newcenter) /= n_neighbours[p];
                Point_z(newcenter) /= n_neighbours[p];
                ADD_POINT_VECTOR(newcenter, polygons->points[p], newcenter);
                set_vector_length(&newcenter, radius);

                /* check for flips */
                pts[0] = newcenter;
                ok = 1;
                for (n = 0; n < n_neighbours[p]; n++) {
                        n1 = neighbours[p][n];
                        n2 = neighbours[p][(n + 1) % n_neighbours[p]];

                        /* normal of the triangle */
                        pts[1] = polygons->points[n1];
                        pts[2] = polygons->points[n2];
                        find_polygon_normal(3, pts, &norm);
                        if (acos(DOT_VECTORS(norm, ipolyinfo[p]->norm[n])) > 0.01) {
                                ok = 0;
                                break;
                        }
                }
                if (ok == 0)
                        continue; /* flipped, so skip it */

                if (Select == 1) {
                        /* see if the new center reduces area distortion! */
                        newdistort = 0; 
                        for (n = 0; n < n_neighbours[p]; n++) {
                                n1 = neighbours[p][n];
                                n2 = neighbours[p][(n + 1) % n_neighbours[p]];
                        
                                /* area of the triangle */
                                pts[0] = newcenter;
                                pts[1] = polygons->points[n1];
                                pts[2] = polygons->points[n2];
                                localarea[n] = get_polygon_surface_area(3, pts);

                                iarea = ipolyinfo[p]->areas[n];
                                if (iarea > 0)
                                        newdistort += log10(localarea[n]/iarea);
                        }
                        newdistort = fabs(newdistort / n_neighbours[p]);
                        if (newdistort >= adistort) {
                                continue; /* skip this point */
                        }
                }

                count++;
                fill_Point(polygons->points[p], Point_x(newcenter),
                           Point_y(newcenter), Point_z(newcenter));
                pts[0] = newcenter;
                for (n = 0; n < n_neighbours[p]; n++) {
                        n1 = neighbours[p][n];
                        n2 = neighbours[p][(n + 1) % n_neighbours[p]];

                        pts[1] = polygons->points[n1];
                        pts[2] = polygons->points[n2];
                        find_polygon_normal(3, pts, &ipolyinfo[p]->norm[n]);
                }
        }
        if (Quiet == 0)
                update_progress_report(&progress, it);

        if (count == 0) {
                printf("Iterations = %d\n", it);
                break; /* done! */
        }
}
        if (Quiet == 0)
                terminate_progress_report(&progress);

        radius = sqrt(get_polygons_surface_area(ipolygons) / (4.0 * PI));
        for (p = 0; p < polygons->n_points; p++)
                set_vector_length(&polygons->points[p], radius);

        compute_polygon_normals(polygons);
        output_graphics_any_format(output_file, format, 1, objects);

        delete_object_list(n_objects, objects);

        return(0);
}
