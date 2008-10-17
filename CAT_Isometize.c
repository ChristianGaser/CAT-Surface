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

#define SMOOTHING 1
#undef AREADISTORTION /* don't run this algorithm by default */
#define STRETCH 1

#define PINF  1.7976931348623157e+308 /* for doubles */
#define NINF -1.7976931348623157e+308 /* for doubles */

int iter = 100; /* number of iterations */

struct ptinfo {
        double *areas;
        double *lengths;
};

static ArgvInfo argTable[] = {
  {"-iter", ARGV_INT, (char *) 1, (char *) &iter,
       "number of iterations."},
   {NULL, ARGV_END, NULL, NULL, NULL}
};

void
usage(char *executable)
{
        char *usage_str =
"\nUsage: %s [options] infile.obj conformalmap.obj outfile.obj\n\n\
    Adjust a conformal map to be more isometric (area preserving).\n\
    The original mesh is infile.obj and the conformal map of the mesh\n\
    is conformalmap.obj.  Results are saved in outfile.obj.\n\n";

        fprintf(stderr, usage_str, executable);
}

int
main(int argc, char** argv)
{
        char               *input_file, *cmap_file, *output_file;
        object_struct      **objects;
        polygons_struct    *ipolygons, *polygons;
        int                *n_neighbours, **neighbours;
        int                n_objects;
        File_formats       format;
        int                i, n, it, p;
        Point              pts[3];
        Point              npt, dir, newcenter;
        double             areas[128];
        double             iarea, centers[384], ratio, totalArea;
        double             adistort, newdistort, totaldistort, weight;
        double             dist, radius, maxdist = 0.05;
        double             bounds[6];
        double             xyz1[3], xyz2[3];
        int                n1, n2, nidx;
        int                Acount, Dcount, Scount;
        struct ptinfo      **ipolyinfo;

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
        radius = sqrt(get_polygons_surface_area(polygons) / (4.0 * PI));

        if (polygons->n_points != ipolygons->n_points
            || polygons->n_items != ipolygons->n_items) {
                printf("Input mesh and conformal map mesh do not match.\n");
                return(1);
        }
    
        create_polygon_point_neighbours(polygons, TRUE, &n_neighbours,
                                         &neighbours, NULL, NULL);

        ratio = get_polygons_surface_area(ipolygons) /
                get_polygons_surface_area(polygons);

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

                /* Get 2 consecutive neighbors of this node */
                for (n = 0; n < n_neighbours[p]; n++) {
                        n1 = neighbours[p][n];
                        n2 = neighbours[p][(n + 1) % n_neighbours[p]];

                        /* area of the triangle */
                        pts[0] = ipolygons->points[p];
                        pts[1] = ipolygons->points[n1];
                        pts[2] = ipolygons->points[n2];
                        ipolyinfo[p]->areas[n] = get_polygon_surface_area(3,
                                                                          pts);

                        /* lengths to neighbouring nodes */
                        SUB_POINTS(dir, ipolygons->points[n1],
                                        ipolygons->points[p]);
                        ipolyinfo[p]->lengths[n] = MAGNITUDE(dir);
                }
        }

for (it = 1; it <= iter; it++) {
        Acount = Dcount = Scount = 0;
        for (p = 0; p < polygons->n_points; p++) {
                if (n_neighbours[p] <= 1)
                        continue; /* skip this point */

                adistort = 0.0;
                totaldistort = 0.0;
                totalArea = 0.0;

                /* Get 2 consecutive neighbors of this node */
                pts[0] = polygons->points[p];
                for (n = 0; n < n_neighbours[p]; n++) {
                        n1 = neighbours[p][n];
                        n2 = neighbours[p][(n + 1) % n_neighbours[p]];

                        /* area of the triangle */
                        pts[1] = polygons->points[n1];
                        pts[2] = polygons->points[n2];
                        areas[n] = get_polygon_surface_area(3, pts);

                        iarea = ipolyinfo[p]->areas[n];

                        totalArea += areas[n];
                        if (iarea > 0) {
                                adistort += log10(ratio*areas[n]/iarea);
                                totaldistort += areas[n]/iarea;
                        }

                        /* Save center of this tile */
                        centers[n*3    ] = (Point_x(polygons->points[p]) +
                                            Point_x(polygons->points[n1]) +
                                            Point_x(polygons->points[n2]))/3.0;
                        centers[n*3 + 1] = (Point_y(polygons->points[p]) +
                                            Point_y(polygons->points[n1]) +
                                            Point_y(polygons->points[n2]))/3.0;
                        centers[n*3 + 2] = (Point_z(polygons->points[p]) +
                                            Point_z(polygons->points[n1]) +
                                            Point_z(polygons->points[n2]))/3.0;
                }
                adistort = fabs(adistort / n_neighbours[p]);

#ifdef SMOOTHING
                /* Algorithm #1 - Area Smoothing */
                for (i = 0; i < 3; i++)
                        xyz1[i] = 0.0;
                for (n = 0; n <  n_neighbours[p]; n++) {
                        iarea = ipolyinfo[p]->areas[n];
                        if (totalArea > 0) {
                                weight = areas[n]/totalArea;
                                for (i = 0; i < 3; i++)
                                        xyz1[i] += weight * centers[n*3+i];
                        }
                }

                /* see if the new centers reduces area distortion! */
                fill_Point(newcenter, xyz1[0], xyz1[1], xyz1[2]);
                set_vector_length(&newcenter, radius);
                newdistort = 0;
                for (n = 0; n < n_neighbours[p]; n++) {
                        n1 = neighbours[p][n];
                        n2 = neighbours[p][(n + 1) % n_neighbours[p]];

                        /* area of the triangle */
                        pts[0] = newcenter;
                        pts[1] = polygons->points[n1];
                        pts[2] = polygons->points[n2];
                        areas[n] = get_polygon_surface_area(3, pts);

                        iarea = ipolyinfo[p]->areas[n];

                        if (iarea > 0)
                                newdistort += log10(ratio*areas[n]/iarea);
                }

                newdistort = fabs(newdistort / n_neighbours[p]);
                if (newdistort < adistort) {
                        Acount++;
                        fill_Point(polygons->points[p],
                                   xyz1[0], xyz1[1], xyz1[2]);
                        continue;
                }
#endif

#ifdef AREADISTORTION
                /* Algorithm #2 - Area Distortion */
                for (i = 0; i < 3; i++)
                        xyz2[i] = 0.0;
                for (n = 0; n <  n_neighbours[p]; n++) {
                        iarea = ipolyinfo[p]->areas[n];
                        if (iarea > 0) {
                                weight = (areas[n]/iarea) / totaldistort;
                                for (i = 0; i < 3; i++)
                                        xyz2[i] += weight * centers[n*3+i];
                        }
                }
                fill_Point(newcenter, xyz2[0], xyz2[1], xyz2[2]);
                set_vector_length(&newcenter, radius);
                newdistort = 0;
                for (n = 0; n < n_neighbours[p]; n++) {
                        n1 = neighbours[p][n];
                        n2 = neighbours[p][(n + 1) % n_neighbours[p]];

                        /* area of the triangle */
                        pts[0] = newcenter;
                        pts[1] = polygons->points[n1];
                        pts[2] = polygons->points[n2];
                        areas[n] = get_polygon_surface_area(3, pts);

                        iarea = ipolyinfo[p]->areas[n];

                        if (iarea > 0)
                                newdistort += log10(ratio*areas[n]/iarea);
                }

                newdistort = fabs(newdistort / n_neighbours[p]);
                if (newdistort < adistort) {
                        Dcount++;
                        fill_Point(polygons->points[p],
                                   xyz2[0], xyz2[1], xyz2[2]);
                        continue;
                }
#endif

#ifdef STRETCH
                /* Algorithm #3 - Stretch Optimization */
                fill_Point(newcenter, 0.0, 0.0, 0.0);

                bounds[0] = bounds[1] = Point_x(polygons->points[p]);
                bounds[2] = bounds[3] = Point_y(polygons->points[p]);
                bounds[4] = bounds[5] = Point_z(polygons->points[p]);

                for (n = 0; n <  n_neighbours[p]; n++) {
                        nidx = neighbours[p][n];
                        npt = polygons->points[nidx];

                        SUB_POINTS(dir, npt, polygons->points[p]);
                        bounds[0] = (bounds[0] < Point_x(npt)) ?
                                    bounds[0] : Point_x(npt);
                        bounds[1] = (bounds[1] > Point_x(npt)) ?
                                    bounds[1] : Point_x(npt);
                        bounds[2] = (bounds[2] < Point_y(npt)) ?
                                    bounds[2] : Point_y(npt);
                        bounds[3] = (bounds[3] > Point_y(npt)) ?
                                    bounds[3] : Point_y(npt);
                        bounds[4] = (bounds[4] < Point_z(npt)) ?
                                    bounds[4] : Point_z(npt);
                        bounds[5] = (bounds[5] > Point_z(npt)) ?
                                    bounds[5] : Point_z(npt);

                        dist = MAGNITUDE(dir) - ipolyinfo[p]->lengths[n];
                        if (dist > maxdist)
                                dist = maxdist;
                        dist /= MAGNITUDE(dir);
                        SCALE_VECTOR(dir, dir, dist);
                        ADD_POINT_VECTOR(newcenter, newcenter, dir);
                }
                for (i = 0; i < 3; i++)
                        Point_coord(newcenter, i) = Point_coord(newcenter, i) /
                                                    n_neighbours[p];

                set_vector_length(&newcenter, radius);

                bounds[0] -= Point_x(polygons->points[p]);
                bounds[1] -= Point_x(polygons->points[p]);
                bounds[2] -= Point_y(polygons->points[p]);
                bounds[3] -= Point_y(polygons->points[p]);
                bounds[4] -= Point_z(polygons->points[p]);
                bounds[5] -= Point_z(polygons->points[p]);

                if (Point_x(newcenter) < bounds[0]/10 ||
                    Point_x(newcenter) > bounds[1]/10 ||
                    Point_y(newcenter) < bounds[2]/10 ||
                    Point_y(newcenter) > bounds[3]/10 ||
                    Point_z(newcenter) < bounds[4]/10 ||
                    Point_z(newcenter) > bounds[5]/10) {
                        continue; /* algorithm gave an aggressive answer */
                }

                ADD_POINT_VECTOR(newcenter, polygons->points[p], newcenter);
                /* see if the new center reduces area distortion! */
                newdistort = 0;
                for (n = 0; n < n_neighbours[p]; n++) {
                        n1 = neighbours[p][n];
                        n2 = neighbours[p][(n + 1) % n_neighbours[p]];

                        /* area of the triangle */
                        pts[0] = newcenter;
                        pts[1] = polygons->points[n1];
                        pts[2] = polygons->points[n2];
                        areas[n] = get_polygon_surface_area(3, pts);

                        iarea = ipolyinfo[p]->areas[n];

                        if (iarea > 0)
                                newdistort += log10(ratio*areas[n]/iarea);
                }
                newdistort = fabs(newdistort / n_neighbours[p]);
                if (newdistort < adistort) {
                        Scount++;
                        fill_Point(polygons->points[p], Point_x(newcenter),
                                   Point_y(newcenter), Point_z(newcenter));
                }
#endif

        }
        /* printf("A = %d, D = %d, S = %d\n", Acount, Dcount, Scount); */
        if (Acount == 0 && Dcount == 0 && Scount == 0) break; /* done! */
}

        convert_ellipsoid_to_sphere_with_surface_area(polygons, 100.0);

        compute_polygon_normals(polygons);
        output_graphics_any_format(output_file, format, 1, objects);

        delete_object_list(n_objects, objects);

        return(0);
}
