/* Rachel Yotter - rachel.yotter@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 */

/* Program to make a conformal map more isometric using a smoothing
 * algorithm. */

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
 
int Iter = 10000; /* number of iterations */
BOOLEAN Select = 0; /* by default, apply to all triangles */
double Threshold = 0; /* by default, apply to all triangles */
BOOLEAN Quiet = 0; /* turn progress reports on and off */

static ArgvInfo argTable[] = {
  {"-iter", ARGV_INT, (char *) 1, (char *) &Iter,
       "number of iterations."},
  {"-select", ARGV_CONSTANT, (char *) TRUE, (char *) &Select,
    "apply smoothing only if it lowers local area distortion" },
  {"-quiet", ARGV_CONSTANT, (char *) TRUE, (char *) &Quiet,
    "turn off progress reports" },
  {"-thresh", ARGV_FLOAT, (char *) 1, (char *) &Threshold,
       "apply smoothing only if local area distortion is above threshold."},
   {NULL, ARGV_END, NULL, NULL, NULL}
};

void
usage(char *executable)
{
        char *usage_str =
"\nUsage: %s [options] infile.obj conformalmap.obj outfile.obj\n\n\
    Adjust a conformal map to be more isometric (area preserving), using the\n\
    smoothing algorithm.  The original original mesh is infile.obj and the\n\
    conformal map of the mesh is conformalmap.obj.  Results are saved in\n\
    outfile.obj.\n\n";

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
        Point              pts[3], newcenter;
        double             **areas, localarea[128], centers[128], totalArea;
        double             adistort, newdistort, weight, radius, xyz[3];
        int                n1, n2, count;
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

        /* map the spherical map to a sphere of the same surface area */
        radius = sqrt(get_polygons_surface_area(ipolygons) / (4.0 * PI));
        for (p = 0; p < polygons->n_points; p++)
                set_vector_length(&polygons->points[p], radius);

        /* first get some info about the original mesh */
        areas = (double **) malloc(sizeof(double *) * ipolygons->n_points);

        for (p = 0; p < ipolygons->n_points; p++) {
                if (n_neighbours[p] <= 1)
                        continue; /* skip this point */

                areas[p] = (double *) malloc(sizeof(double) * n_neighbours[p]);

                /* Get 2 consecutive neighbors of this node */
                for (n = 0; n < n_neighbours[p]; n++) {
                        n1 = neighbours[p][n];
                        n2 = neighbours[p][(n + 1) % n_neighbours[p]];

                        /* area of the triangle */
                        pts[0] = ipolygons->points[p];
                        pts[1] = ipolygons->points[n1];
                        pts[2] = ipolygons->points[n2];
                        areas[p][n] = get_polygon_surface_area(3, pts);
                }
        }

        if (Quiet == 0)
                initialize_progress_report(&progress, FALSE, Iter, "Smooth");

for (it = 1; it <= Iter; it++) {
        count = 0;
        for (p = 0; p < polygons->n_points; p++) {
                if (n_neighbours[p] <= 1)
                        continue; /* skip this point */

                adistort = 0.0;
                totalArea = 0.0;

                /* Get 2 consecutive neighbors of this node */
                pts[0] = polygons->points[p];
                for (n = 0; n < n_neighbours[p]; n++) {
                        n1 = neighbours[p][n];
                        n2 = neighbours[p][(n + 1) % n_neighbours[p]];

                        /* area of the triangle */
                        pts[1] = polygons->points[n1];
                        pts[2] = polygons->points[n2];
                        localarea[n] = get_polygon_surface_area(3, pts);

                        totalArea += localarea[n];
                        if (areas[p][n] > 0)
                                adistort += log10(localarea[n]/areas[p][n]);

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
                if (adistort < Threshold)
                        continue; /* skip this point */

                /* Algorithm #1 - Area Smoothing */
                for (i = 0; i < 3; i++)
                        xyz[i] = 0.0;
                for (n = 0; n <  n_neighbours[p]; n++) {
                        if (totalArea > 0) {
                                weight = localarea[n]/totalArea;
                                for (i = 0; i < 3; i++)
                                        xyz[i] += weight * centers[n*3+i];
                        }
                }

                if (Select == 0) {
                        count++;
                        fill_Point(polygons->points[p], xyz[0], xyz[1], xyz[2]);
                        continue;
                }

                /* see if the new center reduces area distortion! */
                fill_Point(newcenter, xyz[0], xyz[1], xyz[2]);
                newdistort = 0;
                for (n = 0; n < n_neighbours[p]; n++) {
                        n1 = neighbours[p][n];
                        n2 = neighbours[p][(n + 1) % n_neighbours[p]];

                        /* area of the triangle */
                        pts[0] = newcenter;
                        pts[1] = polygons->points[n1];
                        pts[2] = polygons->points[n2];
                        localarea[n] = get_polygon_surface_area(3, pts);

                        if (areas[p][n] > 0)
                                newdistort += log10(localarea[n]/areas[p][n]);
                }

                newdistort = fabs(newdistort / n_neighbours[p]);
                if (newdistort < adistort) {
                        count++;
                        fill_Point(polygons->points[p], xyz[0], xyz[1], xyz[2]);
                }

        }

        for (p = 0; p < polygons->n_points; p++)
                set_vector_length(&polygons->points[p], radius);

        if (Quiet == 0)
                update_progress_report(&progress, it);

        if (count == 0) {
                printf("Iterations = %d\n", it);
                break; /* done! */
        }
}
        if (Quiet == 0)
                terminate_progress_report(&progress);

        compute_polygon_normals(polygons);
        output_graphics_any_format(output_file, format, 1, objects);

        delete_object_list(n_objects, objects);

        return(0);
}
