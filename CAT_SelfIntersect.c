/* Rachel Yotter - rachel.yotter@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id: CAT_SelfIntersections.c 121 2009-07-27 14:42:19Z raytrace $
 *
 */

#include <volume_io/internal_volume_io.h>
#include <volume_io/geometry.h>
#include <bicpl.h>
#include <ParseArgv.h>

#include "CAT_Curvature.h"
#include "CAT_Blur2d.h"
#include "CAT_Octree.h"
#include "CAT_SurfaceIO.h"

#define PINF  1.7976931348623157e+308 /* for doubles */
#define NINF -1.7976931348623157e+308 /* for doubles */


void
usage(char *executable)
{
        char *usage_str =
"\nUsage: %s [options] object_file output_file\n\n\
    Find the number of self-intersections and mark the intersecting\n\
    triangles in the output_file.\n\n";

        fprintf(stderr, usage_str, executable);
}

int
intersect_triangle_triangle(int poly0, int poly1, polygons_struct *surface)
{
        int size, i, p;
        Point pts[3];
        int result;

        size = GET_OBJECT_SIZE(*surface, poly0);
        for (i = 0; i < size; i++) {
                pts[i] = surface->points[surface->indices[
                                  POINT_INDEX(surface->end_indices, poly0, i)]];
        }

        for (i = 1; i < size; i++) {
                result = intersect_segment_triangle(pts[i-1], pts[i], poly1,
                                                    surface);
                if (result > 0)
                        return 1;
        }

        return 0;
}


// intersect_RayTriangle(): intersect a ray with a 3D triangle
//    Input:  a ray R, and a triangle T
//    Output: *I = intersection point (when it exists)
//    Return: -1 = triangle is degenerate (a segment or point)
//             0 = disjoint (no intersect)
//             1 = intersect in unique point I1
//             2 = are in the same plane
int
intersect_segment_triangle(Point p0, Point p1, int poly,
                           polygons_struct *surface)
{
        Vector   u, v, n;             // triangle vectors
        Vector   dir, w0, w;          // ray vectors
        Vector   zero;
        float    r, a, b;             // params to calc ray-plane intersect
        float    uu, uv, vv, wu, wv, D;
        float    s, t;
        Point    pts[3], I;
        int      size, i;
        double   area;

        size = GET_OBJECT_SIZE(*surface, poly);
        if (size != 3) return 0;
        for (i = 0; i < size; i++) {
                pts[i] = surface->points[surface->indices[
                                  POINT_INDEX(surface->end_indices, poly, i)]];
                if (EQUAL_POINTS(pts[i], p0)) return 0;
                if (EQUAL_POINTS(pts[i], p1)) return 0;
        }

        area = get_polygon_surface_area(size, pts);
        if (area == 0) return 0;

        /* get triangle edge vectors and plane normal */
        SUB_POINTS(u, pts[1], pts[0]);
        SUB_POINTS(v, pts[2], pts[0]);

        CROSS_VECTORS(n, u, v);

        fill_Vector(zero, 0.0, 0.0, 0.0);
        if (EQUAL_VECTORS(n, zero))
                return -1; /* triangle is degenerate */

        SUB_POINTS(dir, p1, p0); /* ray direction vector */
        SUB_POINTS(w0, p0, pts[0]);

        a = -DOT_VECTORS(n, w0);
        b = DOT_VECTORS(n, dir);

        if (fabs(b) < 1e-6) { /* ray is parallel to triangle plane */
                if (a == 0) { /* ray lies in triangle plane */
                        return 2;
                } else return 0;             /* ray disjoint from plane */
        }

        /* get intersect point of ray with triangle plane */
        r = a / b;
        if (r < 0.0 || r > 1.0)  /* no intersect */
                return 0;

        
        SCALE_VECTOR(dir, dir, r);
        ADD_POINT_VECTOR(I, p0, dir); /* intersect point of ray and plane */

        /* is I inside T? */
        uu = DOT_VECTORS(u, u);
        uv = DOT_VECTORS(u,v);
        vv = DOT_VECTORS(v,v);
        SUB_POINTS(w, I, pts[0]);
        wu = DOT_VECTORS(w,u);
        wv = DOT_VECTORS(w,v);
        D = uv * uv - uu * vv;

        // get and test parametric coords
        s = (uv * wv - vv * wu) / D;
        if (s < 0.0 || s > 1.0) /* I is outside T */
                return 0;
        t = (uv * wu - uu * wv) / D;
        if (t < 0.0 || (s + t) > 1.0)  /* I is outside T */
                return 0;

        /*printf("(%f, %f, %f) => (%f, %f, %f)\n", Point_x(p0), Point_y(p0), Point_z(p0), Point_x(p1), Point_y(p1), Point_z(p1));
        printf("(%f, %f, %f)\n", Point_x(pts[0]), Point_y(pts[0]), Point_z(pts[0]));
        printf("(%f, %f, %f)\n", Point_x(pts[1]), Point_y(pts[1]), Point_z(pts[1]));
        printf("(%f, %f, %f)\n", Point_x(pts[2]), Point_y(pts[2]), Point_z(pts[2]));
        */
        return 1; /* I is in T */
}

int
main(int argc, char** argv)
{
        char               *in_file, *out_file;
        object_struct      **objects;
        polygons_struct    *polygons;
        double             *intersectflag, **bounds, x, y, z;
        int                n_objects;
        File_formats       format;
        int                p, p2, n_intersects, size, i, pt;
        progress_struct    progress;


        initialize_argument_processing(argc, argv);
        if (!get_string_argument(NULL, &in_file) ||
            !get_string_argument(NULL, &out_file)) {
                fprintf(stderr, "\nUsage: %s object_file output_file\n",
                                argv[0]);
                return(1);
        }

        if (input_graphics_any_format(in_file, &format, &n_objects,
                                      &objects) != OK) {
                printf("Error reading input file %s\n", in_file);
                return(1);
        }

        if (n_objects != 1 || get_object_type(objects[0]) != POLYGONS) {
                printf("Input file must contain one polygon object.\n");
                return(1);
        }

        polygons = get_polygons_ptr(objects[0]);

        intersectflag = (double *) malloc(sizeof(double) * polygons->n_points);
        memset(intersectflag, 0, sizeof(double) * polygons->n_points);

        bounds = (double **) malloc(sizeof(double *) * polygons->n_items);
        for (p = 0; p < polygons->n_items; p++)
                bounds[p] = (double *) malloc(sizeof(double) * 6);

        /* get bounding boxes for all triangles to make it run FASTER */
        for (p = 0; p < polygons->n_items; p++) {
                size = GET_OBJECT_SIZE(*polygons, p);
                pt = polygons->indices[POINT_INDEX(polygons->end_indices,p,0)];

                bounds[p][0] = Point_x(polygons->points[pt]);
                bounds[p][1] = Point_x(polygons->points[pt]);
                bounds[p][2] = Point_y(polygons->points[pt]);
                bounds[p][3] = Point_y(polygons->points[pt]);
                bounds[p][4] = Point_z(polygons->points[pt]);
                bounds[p][5] = Point_z(polygons->points[pt]);

                for (i = 1; i < size; i++) {
                        pt = polygons->indices[
                                     POINT_INDEX(polygons->end_indices, p, i)];
                        x = Point_x(polygons->points[pt]);
                        y = Point_y(polygons->points[pt]);
                        z = Point_z(polygons->points[pt]);

                        bounds[p][0] = bounds[p][0] < x ? bounds[p][0] : x;
                        bounds[p][1] = bounds[p][1] > x ? bounds[p][1] : x;
                        bounds[p][2] = bounds[p][2] < y ? bounds[p][2] : y;
                        bounds[p][3] = bounds[p][3] > y ? bounds[p][3] : y;
                        bounds[p][4] = bounds[p][4] < z ? bounds[p][4] : z;
                        bounds[p][5] = bounds[p][5] > z ? bounds[p][5] : z;
                }
        }

        initialize_progress_report(&progress, FALSE, polygons->n_items,
                                   "Self-Intersection Test");

        n_intersects = 0;
        for (p = 0; p < polygons->n_items; p++) {
                for (p2 = p+1; p2 < polygons->n_items; p2++) {
                        if (xintersect(bounds[p], bounds[p2]) == 0 ||
                            yintersect(bounds[p], bounds[p2]) == 0 ||
                            zintersect(bounds[p], bounds[p2]) == 0)
                                continue;
                        if (intersect_triangle_triangle(p, p2, polygons)) {
                                n_intersects++;

                                size = GET_OBJECT_SIZE(*polygons, p);

                                for (i = 0; i < size; i++) {
                                        pt = polygons->indices[
                                             POINT_INDEX(polygons->end_indices,
                                                         p, i)];
                                        intersectflag[pt] = 1;
                                }
                                size = GET_OBJECT_SIZE(*polygons, p2);

                                for (i = 0; i < size; i++) {
                                        pt = polygons->indices[
                                             POINT_INDEX(polygons->end_indices,
                                                         p2, i)];
                                        intersectflag[pt] = 1;
                                }
                        }

                }
                update_progress_report(&progress, p);
        }

        terminate_progress_report(&progress);

        output_values_any_format(out_file, polygons->n_points, intersectflag);

        printf("Self-Intersections:  %d\n", n_intersects);

        delete_object_list(n_objects, objects);

        free(intersectflag);

        return(0);
}
