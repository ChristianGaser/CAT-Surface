/* Rachel Yotter - rachel.yotter@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id: CAT_FractalSurface.c 140 2009-11-09 12:54:40Z raytrace $
 *
 */

/* Generate a von Koch surface.  Code is not optimized and is fairly
 * ugly, for internal use only!  Please do not distribute.
 */

#include <volume_io/internal_volume_io.h>
#include <volume_io/geometry.h>
#include <bicpl.h>
#include <ParseArgv.h>

#include "Cat_Surf.h"
#include "CAT_SurfaceIO.h"

int level = 0; /* default: generate a normal tetrahedron */
double Scale = 2; /* default: projections are half of "real" von Koch surface */
BOOLEAN squareflag = 0; /* default: tetrahedral output */
BOOLEAN comboflag = 0; /* default: regular surface */
BOOLEAN roiflag = 0; /* default: no ROI output */
BOOLEAN flipflag = 0; /* fractal structures point outwards or alternates */

static ArgvInfo argTable[] = {
  {"-level", ARGV_INT, (char *) 0, (char *) &level,
    "Number of times to subdivide tetrahedral surfaces" },
  { "-square", ARGV_CONSTANT, (char *) TRUE, (char *) &squareflag,
    "Cubic topology rather than tetrahedral." },
  { "-combo", ARGV_CONSTANT, (char *) TRUE, (char *) &comboflag,
    "Different fractal topology for each surface." },
  { "-roi", ARGV_CONSTANT, (char *) TRUE, (char *) &roiflag,
    "Output ROIs for each face into roi.txt file." },
  { "-flip", ARGV_CONSTANT, (char *) TRUE, (char *) &flipflag,
    "Alternate inward/outward projections of fractal structures." },
  {"-scale", ARGV_FLOAT, (char *) 0, (char *) &Scale,
    "Divide the projections by this amount... 1 gives a true von Koch surface" },
   {NULL, ARGV_END, NULL, NULL, NULL}
};

/* Create a fractal surface using triangles and tetrahedron base */
object_struct **
create_von_koch_surface(int iterations)
{
        polygons_struct *polygons;
        object_struct **object;
        int iter, n_points, n_items, n_p, n_i, *indices, idx;
        int n, n2, i, a, b, c, p, p2, count, offset, *pflag, *flip;
        Point *points, tri[3], median, dir;
        Vector *normals, normal;
        double len, dist, scale;
        int *polynum, *pointnum;
        FILE *fp;

        /* overestimate the # of points and triangles initially */
        //n_points = 4 * pow(6, iterations);
        //n_items = 4 * pow(6, iterations);
        n_points = 4 * pow(6, 7);
        n_items = 4 * pow(6, 7);

        points = (Point *) malloc(sizeof(Point) * n_points);
        pflag = (int *) malloc(sizeof(int) * n_points);
        memset(pflag, 0, sizeof(int) * n_points);
        indices = (int *) malloc(sizeof(int) * n_items * 3);
        normals = (Vector *) malloc(sizeof(Vector) * n_items);

        if (flipflag) {
                flip = (int *) malloc(sizeof(int) * n_items);
                memset(flip, 0, sizeof(int) * n_items);
        }

        if (roiflag) {
                polynum = (int *) malloc(sizeof(int) * n_items);
                pointnum = (int *) malloc(sizeof(int) * n_points);
        }

        /* create the base tetrahedron */
        n_points = 4; n_items = 4;
        fill_Point(points[0], 1.0, 1.0, 1.0);
        fill_Point(points[1], -1.0, -1.0, 1.0);
        fill_Point(points[2], 1.0, -1.0, -1.0);
        fill_Point(points[3], -1.0, 1.0, -1.0);
        indices[0] = 0; indices[1] = 1; indices[2] = 2;
        indices[3] = 0; indices[4] = 2; indices[5] = 3;
        indices[6] = 0; indices[7] = 3; indices[8] = 1;
        indices[9] = 1; indices[10] = 3; indices[11] = 2;
        if (roiflag) {
                polynum[0] = 1; polynum[1] = 2; polynum[2] = 4; polynum[3] = 8;
        }

        SUB_POINTS(tri[0], points[1], points[0]);
        len = sqrt(DOT_POINTS(tri[0], tri[0]));

        /* cut each triangle into 6 triangles */
        for (iter = 0; iter < iterations; iter++) {
                len /= 2;

                for (n = 0; n < n_items; n++) {
                        tri[0] = points[indices[n*3]];
                        tri[1] = points[indices[n*3 + 1]];
                        tri[2] = points[indices[n*3 + 2]];
                        find_polygon_normal(3, tri, &normals[n]);
                }

                n_i = n_items;
                for (n = 0; n < n_i; n++) {
                        a = indices[n*3];
                        b = indices[n*3 + 1];
                        c = indices[n*3 + 2];

                        INTERPOLATE_POINTS(points[n_points],
                                           points[a], points[b], 0.5);
                        pflag[n_points] = 1;
                        INTERPOLATE_POINTS(points[n_points+1],
                                           points[b], points[c], 0.5);
                        pflag[n_points + 1] = 1;
                        INTERPOLATE_POINTS(points[n_points+2],
                                           points[a], points[c], 0.5);
                        pflag[n_points + 2] = 1;
                        n_points += 3;

                        INTERPOLATE_POINTS(median, points[n_points-3],
                                           points[n_points-2], 0.5);
                        SUB_POINTS(median, median, points[n_points-1]);
                        SCALE_POINT(median, median, 2.0f/3.0f);

                        /* find the peak point */
                        points[n_points] = points[n_points-1];
                        ADD_POINTS(points[n_points], points[n_points], median);
                        /* a real von koch surface */
                        //SCALE_POINT(normal, normals[n], sqrt(6.0f)*len/3.0f);
                        /* our modified one to avoid intersections */
                        if (!comboflag) {
                                SCALE_POINT(normal, normals[n], sqrt(6.0f)*len /
                                                                (3.0f * Scale));
                        } else { /* make combination surface */
                                if (n == 0 || (n >= 4 && n <= 8) ||
                                    (n >= 24 && n <= 28) ||
                                    (n >= 44 && n <= 68)) { /* 1.5 */
                                        SCALE_POINT(normal, normals[n],
                                                    sqrt(6.0f)*len / 4.5f);
                                } else if (n == 1 || (n >= 9 && n <= 13) ||
                                           (n >= 29 && n <= 33) ||
                                           (n >= 69 && n <= 93)) { /* 2.0 */
                                        SCALE_POINT(normal, normals[n],
                                                    sqrt(6.0f)*len / 6.0f);
                                } else if (n == 2 || (n >= 14 && n <= 18) ||
                                           (n >= 34 && n <= 38) ||
                                           (n >= 94 && n <= 118)) { /* 3.0 */
                                        SCALE_POINT(normal, normals[n],
                                                    sqrt(6.0f)*len / 9.0f);
                                } else if (n == 3 || (n >= 19 && n <= 23) ||
                                           (n >= 39 && n <= 43) ||
                                           (n >= 119 && n <= 143)) { /* 4.0 */
                                        SCALE_POINT(normal, normals[n],
                                                    sqrt(6.0f)*len / 12.0f);
                                }
                        }

                        if (flipflag && flip[n])
                                SCALE_POINT(normal, normal, -1.0f);
                        ADD_POINTS(points[n_points], points[n_points], normal);
                        n_points++;

                        /* add the triangle indices */
                        indices[n*3 + 1] = n_points - 4;
                        indices[n*3 + 2] = n_points - 2;

                        indices[n_items*3] = n_points - 4;
                        indices[n_items*3 + 1] = b;
                        indices[n_items*3 + 2] = n_points - 3;
                        if (flipflag) flip[n_items] = (flip[n] == 0);
                        if (roiflag) polynum[n_items] = polynum[n];
                        n_items++;
                        indices[n_items*3] = n_points - 3;
                        indices[n_items*3 + 1] = c;
                        indices[n_items*3 + 2] = n_points - 2;
                        if (flipflag) flip[n_items] = (flip[n] == 0);
                        if (roiflag) polynum[n_items] = polynum[n];
                        n_items++;
                        indices[n_items*3] = n_points - 4;
                        indices[n_items*3 + 1] = n_points - 1;
                        indices[n_items*3 + 2] = n_points - 2;
                        if (flipflag) flip[n_items] = (flip[n] == 0);
                        if (roiflag) polynum[n_items] = polynum[n];
                        n_items++;
                        indices[n_items*3] = n_points - 4;
                        indices[n_items*3 + 1] = n_points - 3;
                        indices[n_items*3 + 2] = n_points - 1;
                        if (flipflag) flip[n_items] = (flip[n] == 0);
                        if (roiflag) polynum[n_items] = polynum[n];
                        n_items++;
                        indices[n_items*3] = n_points - 1;
                        indices[n_items*3 + 1] = n_points - 3;
                        indices[n_items*3 + 2] = n_points - 2;
                        if (flipflag) flip[n_items] = (flip[n] == 0);
                        if (roiflag) polynum[n_items] = polynum[n];
                        n_items++;

                        if (flipflag) flip[n] = (flip[n] == 0);
                }
                /* remove duplicate points */
                for (p = 0; p < n_points; p++) {
                        if (pflag[p] == 0) continue; /* skip */
                        for (p2 = p+1; p2 < n_points; p2++) {
                                dist = distance_between_points(&points[p],
                                                               &points[p2]);
                                if (dist < 0.01*len) { /* delete it */
                                        n_points--;
                                        for (n = 0; n < n_items * 3; n++) {
                                                if (indices[n] == p2)
                                                        indices[n] = p;
                                                else if (indices[n] == n_points)
                                                        indices[n] = p2;
                                        }
                                        points[p2] = points[n_points];
                                        p2--;
                                }
                        }
                }
        }

        /* clean up any rounding ickiness */
        dist = 0.01*len;
        for (p = 0; p < n_points; p++) {
                if (fabs(Point_x(points[p])) < dist)
                        Point_x(points[p]) = 0;
                if (fabs(Point_y(points[p])) < dist)
                        Point_y(points[p]) = 0;
                if (fabs(Point_z(points[p])) < dist)
                        Point_z(points[p])  = 0;
        }

        /* cut each triangle into 4 triangles */
        for (; iter < 7; iter++) {
                len /= 2;

                n_i = n_items;
                for (n = 0; n < n_i; n++) {
                        a = indices[n*3];
                        b = indices[n*3 + 1];
                        c = indices[n*3 + 2];

                        INTERPOLATE_POINTS(points[n_points],
                                           points[a], points[b], 0.5);
                        INTERPOLATE_POINTS(points[n_points+1],
                                           points[b], points[c], 0.5);
                        INTERPOLATE_POINTS(points[n_points+2],
                                           points[a], points[c], 0.5);
                        n_points += 3;

                        /* add the triangle indices */
                        indices[n*3 + 1] = n_points - 3;
                        indices[n*3 + 2] = n_points - 1;

                        indices[n_items*3] = n_points - 3;
                        indices[n_items*3 + 1] = b;
                        indices[n_items*3 + 2] = n_points - 2;
                        if (roiflag) polynum[n_items] = polynum[n];
                        n_items++;
                        indices[n_items*3] = n_points - 2;
                        indices[n_items*3 + 1] = c;
                        indices[n_items*3 + 2] = n_points - 1;
                        if (roiflag) polynum[n_items] = polynum[n];
                        n_items++;
                        indices[n_items*3 + 2] = n_points - 3;
                        indices[n_items*3] = n_points - 2;
                        indices[n_items*3 + 1] = n_points - 1;
                        if (roiflag) polynum[n_items] = polynum[n];
                        n_items++;
                }
                /* remove duplicate points */
                for (p = 0; p < n_points; p++) {
                        for (p2 = p+1; p2 < n_points; p2++) {
                                dist = distance_between_points(&points[p],
                                                               &points[p2]);
                                if (dist < 0.01*len) { /* delete it */
                                        n_points--;
                                        for (n = 0; n < n_items * 3; n++) {
                                                if (indices[n] == p2)
                                                        indices[n] = p;
                                                else if (indices[n] == n_points)
                                                        indices[n] = p2;
                                        }
                                        points[p2] = points[n_points];
                                        p2--;
                                }
                        }
                }
        }

        /* build surface */
        object = (object_struct **) malloc(sizeof(object_struct *));
        *object = create_object(POLYGONS);
        polygons = get_polygons_ptr(*object);
        initialize_polygons(polygons, WHITE, NULL);

        Surfprop_a(polygons->surfprop) = 0.3;
        Surfprop_d(polygons->surfprop) = 0.6;
        Surfprop_s(polygons->surfprop) = 0.6;
        Surfprop_se(polygons->surfprop) = 60;
        Surfprop_t(polygons->surfprop) = 1.0;

        polygons->colour_flag = 0;
        polygons->line_thickness = 1.0f;

        polygons->points = (Point *) malloc(sizeof(Point) * n_points);
        polygons->indices = (int *) malloc(sizeof(int) * 3 * n_items);
        polygons->end_indices = (int *) malloc(sizeof(int) * n_items);
        polygons->normals = (Vector *) malloc(sizeof(Vector) * n_points);
        polygons->n_points = n_points;
        polygons->n_items = n_items;

        for (i = 0; i < polygons->n_points; i++)
                polygons->points[i] = points[i];

        for (i = 0; i < n_items; i++) {
                polygons->end_indices[i] = 3 * (i + 1);
                polygons->indices[i] = indices[i];
        }

        for (i = n_items; i < n_items * 3; i++) 
                polygons->indices[i] = indices[i];

       if (roiflag) {
               for (i = 0; i < n_items; i++) {
                        for (n = 0; n < 3; n++) {
                                idx = polygons->indices[i*3 + n];
                                pointnum[idx] = pointnum[idx] | polynum[i];
                        }
                }
                fp = fopen("roi.txt", "w");
                for (i = 0; i < n_points; i++)
                        fprintf(fp, " %d\n", pointnum[i]);
                fclose(fp);
        }

        for (i = 0; i < n_items; i++) {
                polygons->end_indices[i] = 3 * (i + 1);
                polygons->indices[i] = indices[i];
        }

        for (i = n_items; i < n_items * 3; i++) 
                polygons->indices[i] = indices[i];

        compute_polygon_normals(polygons);

        free(points);
        free(pflag);
        free(indices);
        free(normals);

        if (flipflag) free(flip);

        if (roiflag) {
                free(polynum);
                free(pointnum);
        }

        return(object);
}

/* Create a square von koch surface */
object_struct **
create_square_von_koch_surface(int iterations)
{
        polygons_struct *polygons;
        object_struct **object;
        int iter, n_points, n_items, n_p, n_squares, *indices;
        int n, n2, i, a, b, c, d, p, p2, count, offset, *flip;
        Point *points, tri[3], a_pt, b_pt, c_pt, d_pt, median, dir;
        Vector *normals, normal;
        double len, dist, scale1 = 0.3f, scale2 = 0.6f;
        int *squarenum, *pointnum;
        FILE *fp;

        /* overestimate the # of points and triangles initially */
        n_points = 2 + 6 * pow(20, 4);
        n_items = 6 + 6 * pow(26, 4);

        points = (Point *) malloc(sizeof(Point) * n_points);
        indices = (int *) malloc(sizeof(int) * n_items * 3);
        normals = (Vector *) malloc(sizeof(Vector) * n_items);

        if (flipflag) {
                flip = (int *) malloc(sizeof(int) * n_items / 2);
                memset(flip, 0, sizeof(int) * n_items / 2);
        }

        squarenum = (int *) malloc(sizeof(int) * n_items / 2);
        pointnum = (int *) malloc(sizeof(int) * n_points);

        /* create the base cube */
        n_points = 8; n_items = 12;
        fill_Point(points[0], -1.0, -1.0, -1.0);
        fill_Point(points[1], -1.0, -1.0, 1.0);
        fill_Point(points[2], -1.0, 1.0, -1.0);
        fill_Point(points[3], -1.0, 1.0, 1.0);
        fill_Point(points[4], 1.0, -1.0, -1.0);
        fill_Point(points[5], 1.0, -1.0, 1.0);
        fill_Point(points[6], 1.0, 1.0, -1.0);
        fill_Point(points[7], 1.0, 1.0, 1.0);
        indices[0] = 0; indices[1] = 4; indices[2] = 1;
        indices[3] = 1; indices[4] = 4; indices[5] = 5;
        indices[6] = 4; indices[7] = 6; indices[8] = 5;
        indices[9] = 5; indices[10] = 6; indices[11] = 7;
        indices[12] = 6; indices[13] = 2; indices[14] = 7;
        indices[15] = 7; indices[16] = 2; indices[17] = 3;
        indices[18] = 2; indices[19] = 0; indices[20] = 3;
        indices[21] = 3; indices[22] = 0; indices[23] = 1;
        indices[24] = 1; indices[25] = 5; indices[26] = 3;
        indices[27] = 3; indices[28] = 5; indices[29] = 7;
        indices[30] = 2; indices[31] = 6; indices[32] = 0;
        indices[33] = 0; indices[34] = 6; indices[35] = 4;

        squarenum[0] = 1; squarenum[1] = 2; squarenum[2] = 4;
        squarenum[3] = 8; squarenum[4] = 16; squarenum[5] = 32;

        SUB_POINTS(tri[0], points[1], points[0]);
        len = sqrt(DOT_POINTS(tri[0], tri[0]));

        /* cut each square into thirds */
        for (iter = 0; iter < iterations; iter++) {
                len /= 3;
                n_squares = n_items / 2;
                for (n = 0; n < n_squares; n++) {
                        tri[0] = points[indices[n*6]];
                        tri[1] = points[indices[n*6 + 1]];
                        tri[2] = points[indices[n*6 + 2]];
                        find_polygon_normal(3, tri, &normals[n]);
                }

                for (n = 0; n < n_squares; n++) {
                        b = indices[n*6];
                        c = indices[n*6 + 1];
                        a = indices[n*6 + 2];
                        d = indices[n*6 + 5];
                        a_pt = points[a];
                        b_pt = points[b];
                        c_pt = points[c];
                        d_pt = points[d];

                        /* make perimeter points */
                        INTERPOLATE_POINTS(points[n_points], a_pt, d_pt,
                                           1.0/3.0);
                        INTERPOLATE_POINTS(points[n_points+1], a_pt, d_pt,
                                           2.0/3.0);
                        INTERPOLATE_POINTS(points[n_points+5], d_pt, c_pt,
                                           1.0/3.0);
                        INTERPOLATE_POINTS(points[n_points+9], d_pt, c_pt,
                                           2.0/3.0);
                        INTERPOLATE_POINTS(points[n_points+2], a_pt, b_pt,
                                           1.0/3.0);
                        INTERPOLATE_POINTS(points[n_points+6], a_pt, b_pt,
                                           2.0/3.0);
                        INTERPOLATE_POINTS(points[n_points+10], b_pt, c_pt,
                                           1.0/3.0);
                        INTERPOLATE_POINTS(points[n_points+11], b_pt, c_pt,
                                           2.0/3.0);

                        /* make the four inner points */
                        INTERPOLATE_POINTS(points[n_points+3],
                                           points[n_points+2],
                                           points[n_points+5], 1.0/3.0);
                        INTERPOLATE_POINTS(points[n_points+4],
                                           points[n_points+2],
                                           points[n_points+5], 2.0/3.0);
                        INTERPOLATE_POINTS(points[n_points+7],
                                           points[n_points+6],
                                           points[n_points+9], 1.0/3.0);
                        INTERPOLATE_POINTS(points[n_points+8],
                                           points[n_points+6],
                                           points[n_points+9], 2.0/3.0);

                        /* a real von koch surface */
                        //SCALE_VECTOR(normal, normals[n], len);
                        /* our modified one to avoid intersections */
                        if (!comboflag) {
                                SCALE_VECTOR(normal, normals[n], len / Scale);
                        } else { /* make combination surface */
                                if (squarenum[n] == 1) { /* 1.5 */
                                        SCALE_VECTOR(normal, normals[n],
                                                    len / 1.5f);
                                } else if (squarenum[n] == 2) { /* 2.0 */
                                        SCALE_VECTOR(normal, normals[n],
                                                    len / 2.0f);
                                } else if (squarenum[n] == 4) { /* 2.5 */
                                        SCALE_VECTOR(normal, normals[n],
                                                    len / 2.5f);
                                } else if (squarenum[n] == 8) { /* 3.0 */
                                        SCALE_VECTOR(normal, normals[n],
                                                    len / 3.0f);
                                } else if (squarenum[n] == 16) { /* 3.5 */
                                        SCALE_VECTOR(normal, normals[n],
                                                    len / 3.5f);
                                } else if (squarenum[n] == 32) { /* 4.0 */
                                        SCALE_VECTOR(normal, normals[n],
                                                    len / 4.0f);
                                } else {
                                        printf("ERROR! %d\n", n);
                                }
                        }
                        if (flipflag && flip[n])
                                SCALE_VECTOR(normal, normal, -1.0);

                        /* make the other four points */
                        ADD_POINTS(points[n_points+12], points[n_points+3],
                                   normal);
                        ADD_POINTS(points[n_points+13], points[n_points+4],
                                   normal);
                        ADD_POINTS(points[n_points+14], points[n_points+7],
                                   normal);
                        ADD_POINTS(points[n_points+15], points[n_points+8],
                                   normal);

                        /* add the triangle indices */
                        indices[n*6] = n_points + 2;
                        indices[n*6 + 1] = indices[n*6 + 4] = n_points + 3;
                        indices[n*6 + 5] = n_points;

                        indices[n_items*3] = n_points + 3;
                        indices[n_items*3 + 1] = n_points + 4;
                        indices[n_items*3 + 2] = n_points;
                        if (flipflag) flip[n_items/2] = (flip[n] == 0);
                        squarenum[n_items/2] = squarenum[n];
                        n_items++;
                        indices[n_items*3] = n_points;
                        indices[n_items*3 + 1] = n_points + 4;
                        indices[n_items*3 + 2] = n_points + 1;
                        n_items++;
                        indices[n_items*3] = n_points + 4;
                        indices[n_items*3 + 1] = n_points + 5;
                        indices[n_items*3 + 2] = n_points + 1;
                        if (flipflag) flip[n_items/2] = (flip[n] == 0);
                        squarenum[n_items/2] = squarenum[n];
                        n_items++;
                        indices[n_items*3] = n_points + 1;
                        indices[n_items*3 + 1] = n_points + 5;
                        indices[n_items*3 + 2] = d;
                        n_items++;
                        indices[n_items*3] = n_points + 6;
                        indices[n_items*3 + 1] = n_points + 7;
                        indices[n_items*3 + 2] = n_points + 2;
                        if (flipflag) flip[n_items/2] = (flip[n] == 0);
                        squarenum[n_items/2] = squarenum[n];
                        n_items++;
                        indices[n_items*3] = n_points + 2;
                        indices[n_items*3 + 1] = n_points + 7;
                        indices[n_items*3 + 2] = n_points + 3;
                        n_items++;

                        /* center cube */
                        indices[n_items*3] = n_points + 3;
                        indices[n_items*3 + 1] = n_points + 7;
                        indices[n_items*3 + 2] = n_points + 12;
                        if (flipflag) flip[n_items/2] = flip[n];
                        squarenum[n_items/2] = squarenum[n];
                        n_items++;
                        indices[n_items*3] = n_points + 12;
                        indices[n_items*3 + 1] = n_points + 7;
                        indices[n_items*3 + 2] = n_points + 14;
                        n_items++;
                        indices[n_items*3] = n_points + 7;
                        indices[n_items*3 + 1] = n_points + 8;
                        indices[n_items*3 + 2] = n_points + 14;
                        if (flipflag) flip[n_items/2] = flip[n];
                        squarenum[n_items/2] = squarenum[n];
                        n_items++;
                        indices[n_items*3] = n_points + 14;
                        indices[n_items*3 + 1] = n_points + 8;
                        indices[n_items*3 + 2] = n_points + 15;
                        n_items++;
                        indices[n_items*3] = n_points + 8;
                        indices[n_items*3 + 1] = n_points + 4;
                        indices[n_items*3 + 2] = n_points + 15;
                        if (flipflag) flip[n_items/2] = flip[n];
                        squarenum[n_items/2] = squarenum[n];
                        n_items++;
                        indices[n_items*3] = n_points + 15;
                        indices[n_items*3 + 1] = n_points + 4;
                        indices[n_items*3 + 2] = n_points + 13;
                        n_items++;
                        indices[n_items*3] = n_points + 4;
                        indices[n_items*3 + 1] = n_points + 3;
                        indices[n_items*3 + 2] = n_points + 13;
                        if (flipflag) flip[n_items/2] = flip[n];
                        squarenum[n_items/2] = squarenum[n];
                        n_items++;
                        indices[n_items*3] = n_points + 13;
                        indices[n_items*3 + 1] = n_points + 3;
                        indices[n_items*3 + 2] = n_points + 12;
                        n_items++;
                        indices[n_items*3] = n_points + 14;
                        indices[n_items*3 + 1] = n_points + 15;
                        indices[n_items*3 + 2] = n_points + 12;
                        if (flipflag) flip[n_items/2] = flip[n];
                        squarenum[n_items/2] = squarenum[n];
                        n_items++;
                        indices[n_items*3] = n_points + 12;
                        indices[n_items*3 + 1] = n_points + 15;
                        indices[n_items*3 + 2] = n_points + 13;
                        n_items++;

                        indices[n_items*3] = n_points + 8;
                        indices[n_items*3 + 1] = n_points + 9;
                        indices[n_items*3 + 2] = n_points + 4;
                        if (flipflag) flip[n_items/2] = (flip[n] == 0);
                        squarenum[n_items/2] = squarenum[n];
                        n_items++;
                        indices[n_items*3] = n_points + 4;
                        indices[n_items*3 + 1] = n_points + 9;
                        indices[n_items*3 + 2] = n_points + 5;
                        n_items++;
                        indices[n_items*3] = b;
                        indices[n_items*3 + 1] = n_points + 10;
                        indices[n_items*3 + 2] = n_points + 6;
                        if (flipflag) flip[n_items/2] = (flip[n] == 0);
                        squarenum[n_items/2] = squarenum[n];
                        n_items++;
                        indices[n_items*3] = n_points + 6;
                        indices[n_items*3 + 1] = n_points + 10;
                        indices[n_items*3 + 2] = n_points + 7;
                        n_items++;
                        indices[n_items*3] = n_points + 10;
                        indices[n_items*3 + 1] = n_points + 11;
                        indices[n_items*3 + 2] = n_points + 7;
                        if (flipflag) flip[n_items/2] = (flip[n] == 0);
                        squarenum[n_items/2] = squarenum[n];
                        n_items++;
                        indices[n_items*3] = n_points + 7;
                        indices[n_items*3 + 1] = n_points + 11;
                        indices[n_items*3 + 2] = n_points + 8;
                        n_items++;
                        indices[n_items*3] = n_points + 11;
                        indices[n_items*3 + 1] = c;
                        indices[n_items*3 + 2] = n_points + 8;
                        if (flipflag) flip[n_items/2] = (flip[n] == 0);
                        squarenum[n_items/2] = squarenum[n];
                        n_items++;
                        indices[n_items*3] = n_points + 8;
                        indices[n_items*3 + 1] = c;
                        indices[n_items*3 + 2] = n_points + 9;
                        n_items++;

                        n_points += 16;
                        if (flipflag) flip[n] = (flip[n] == 0);

                }

                /* remove duplicate points */
                for (p = 0; p < n_points; p++) {
                        for (p2 = p+1; p2 < n_points; p2++) {
                                dist = distance_between_points(&points[p],
                                                               &points[p2]);
                                if (dist < 0.01*len) { /* delete it */
                                        n_points--;
                                        for (n = 0; n < n_items * 3; n++) {
                                                if (indices[n] == p2)
                                                        indices[n] = p;
                                                else if (indices[n] == n_points)
                                                        indices[n] = p2;
                                        }
                                        points[p2] = points[n_points];
                                        p2--;
                                }
                        }
                }
        }

        /* clean up any rounding ickiness */
        dist = 0.01*len;
        for (p = 0; p < n_points; p++) {
                if (fabs(Point_x(points[p])) < dist)
                        Point_x(points[p]) = 0;
                if (fabs(Point_y(points[p])) < dist)
                        Point_y(points[p]) = 0;
                if (fabs(Point_z(points[p])) < dist)
                        Point_z(points[p])  = 0;
        }

        /* increase the resolution */
        for (; iter < 4; iter++) {
                len /= 3;
                n_squares = n_items / 2;
                for (n = 0; n < n_squares; n++) {
                        tri[0] = points[indices[n*6]];
                        tri[1] = points[indices[n*6 + 1]];
                        tri[2] = points[indices[n*6 + 2]];
                        find_polygon_normal(3, tri, &normals[n]);
                }

                for (n = 0; n < n_squares; n++) {
                        b = indices[n*6];
                        c = indices[n*6 + 1];
                        a = indices[n*6 + 2];
                        d = indices[n*6 + 5];
                        a_pt = points[a];
                        b_pt = points[b];
                        c_pt = points[c];
                        d_pt = points[d];

                        /* make perimeter points */
                        INTERPOLATE_POINTS(points[n_points], a_pt, d_pt,
                                           1.0/3.0);
                        INTERPOLATE_POINTS(points[n_points+1], a_pt, d_pt,
                                           2.0/3.0);
                        INTERPOLATE_POINTS(points[n_points+5], d_pt, c_pt,
                                           1.0/3.0);
                        INTERPOLATE_POINTS(points[n_points+9], d_pt, c_pt,
                                           2.0/3.0);
                        INTERPOLATE_POINTS(points[n_points+2], a_pt, b_pt,
                                           1.0/3.0);
                        INTERPOLATE_POINTS(points[n_points+6], a_pt, b_pt,
                                           2.0/3.0);
                        INTERPOLATE_POINTS(points[n_points+10], b_pt, c_pt,
                                           1.0/3.0);
                        INTERPOLATE_POINTS(points[n_points+11], b_pt, c_pt,
                                           2.0/3.0);

                        /* make the four inner points */
                        INTERPOLATE_POINTS(points[n_points+3],
                                           points[n_points+2],
                                           points[n_points+5], 1.0/3.0);
                        INTERPOLATE_POINTS(points[n_points+4],
                                           points[n_points+2],
                                           points[n_points+5], 2.0/3.0);
                        INTERPOLATE_POINTS(points[n_points+7],
                                           points[n_points+6],
                                           points[n_points+9], 1.0/3.0);
                        INTERPOLATE_POINTS(points[n_points+8],
                                           points[n_points+6],
                                           points[n_points+9], 2.0/3.0);

                        /* add the triangle indices */
                        indices[n*6] = n_points + 2;
                        indices[n*6 + 1] = indices[n*6 + 4] = n_points + 3;
                        indices[n*6 + 5] = n_points;

                        indices[n_items*3] = n_points + 3;
                        indices[n_items*3 + 1] = n_points + 4;
                        indices[n_items*3 + 2] = n_points;
                        n_items++;
                        indices[n_items*3] = n_points;
                        indices[n_items*3 + 1] = n_points + 4;
                        indices[n_items*3 + 2] = n_points + 1;
                        squarenum[n_items/2] = squarenum[n];
                        n_items++;
                        indices[n_items*3] = n_points + 4;
                        indices[n_items*3 + 1] = n_points + 5;
                        indices[n_items*3 + 2] = n_points + 1;
                        n_items++;
                        indices[n_items*3] = n_points + 1;
                        indices[n_items*3 + 1] = n_points + 5;
                        indices[n_items*3 + 2] = d;
                        squarenum[n_items/2] = squarenum[n];
                        n_items++;
                        indices[n_items*3] = n_points + 6;
                        indices[n_items*3 + 1] = n_points + 7;
                        indices[n_items*3 + 2] = n_points + 2;
                        n_items++;
                        indices[n_items*3] = n_points + 2;
                        indices[n_items*3 + 1] = n_points + 7;
                        indices[n_items*3 + 2] = n_points + 3;
                        squarenum[n_items/2] = squarenum[n];
                        n_items++;

                        /* center square */
                        indices[n_items*3] = n_points + 7;
                        indices[n_items*3 + 1] = n_points + 8;
                        indices[n_items*3 + 2] = n_points + 3;
                        n_items++;
                        indices[n_items*3] = n_points + 3;
                        indices[n_items*3 + 1] = n_points + 8;
                        indices[n_items*3 + 2] = n_points + 4;
                        squarenum[n_items/2] = squarenum[n];
                        n_items++;

                        indices[n_items*3] = n_points + 8;
                        indices[n_items*3 + 1] = n_points + 9;
                        indices[n_items*3 + 2] = n_points + 4;
                        n_items++;
                        indices[n_items*3] = n_points + 4;
                        indices[n_items*3 + 1] = n_points + 9;
                        indices[n_items*3 + 2] = n_points + 5;
                        squarenum[n_items/2] = squarenum[n];
                        n_items++;
                        indices[n_items*3] = b;
                        indices[n_items*3 + 1] = n_points + 10;
                        indices[n_items*3 + 2] = n_points + 6;
                        n_items++;
                        indices[n_items*3] = n_points + 6;
                        indices[n_items*3 + 1] = n_points + 10;
                        indices[n_items*3 + 2] = n_points + 7;
                        squarenum[n_items/2] = squarenum[n];
                        n_items++;
                        indices[n_items*3] = n_points + 10;
                        indices[n_items*3 + 1] = n_points + 11;
                        indices[n_items*3 + 2] = n_points + 7;
                        n_items++;
                        indices[n_items*3] = n_points + 7;
                        indices[n_items*3 + 1] = n_points + 11;
                        indices[n_items*3 + 2] = n_points + 8;
                        squarenum[n_items/2] = squarenum[n];
                        n_items++;
                        indices[n_items*3] = n_points + 11;
                        indices[n_items*3 + 1] = c;
                        indices[n_items*3 + 2] = n_points + 8;
                        n_items++;
                        indices[n_items*3] = n_points + 8;
                        indices[n_items*3 + 1] = c;
                        indices[n_items*3 + 2] = n_points + 9;
                        squarenum[n_items/2] = squarenum[n];
                        n_items++;

                        n_points += 12;

                }

                /* remove duplicate points */
                for (p = 0; p < n_points; p++) {
                        for (p2 = p+1; p2 < n_points; p2++) {
                                dist = distance_between_points(&points[p],
                                                               &points[p2]);
                                if (dist < 0.01*len) { /* delete it */
                                        n_points--;
                                        for (n = 0; n < n_items * 3; n++) {
                                                if (indices[n] == p2)
                                                        indices[n] = p;
                                                else if (indices[n] == n_points)
                                                        indices[n] = p2;
                                        }
                                        points[p2] = points[n_points];
                                        p2--;
                                }
                        }
                }
        }

        /* build surface */
        object = (object_struct **) malloc(sizeof(object_struct *));
        *object = create_object(POLYGONS);
        polygons = get_polygons_ptr(*object);
        initialize_polygons(polygons, WHITE, NULL);

        Surfprop_a(polygons->surfprop) = 0.3;
        Surfprop_d(polygons->surfprop) = 0.6;
        Surfprop_s(polygons->surfprop) = 0.6;
        Surfprop_se(polygons->surfprop) = 60;
        Surfprop_t(polygons->surfprop) = 1.0;

        polygons->colour_flag = 0;
        polygons->line_thickness = 1.0f;

        polygons->points = (Point *) malloc(sizeof(Point) * n_points);
        polygons->indices = (int *) malloc(sizeof(int) * 3 * n_items);
        polygons->end_indices = (int *) malloc(sizeof(int) * n_items);
        polygons->normals = (Vector *) malloc(sizeof(Vector) * n_points);
        polygons->n_points = n_points;
        polygons->n_items = n_items;

        for (i = 0; i < polygons->n_points; i++)
                polygons->points[i] = points[i];

        for (i = 0; i < n_items; i++) {
                polygons->end_indices[i] = 3 * (i + 1);
                polygons->indices[i] = indices[i];
        }

        for (i = n_items; i < n_items * 3; i++) 
                polygons->indices[i] = indices[i];

        compute_polygon_normals(polygons);

        if (roiflag) {
               for (i = 0; i < n_items; i++) {
                        for (n = 0; n < 3; n++) {
                                p = polygons->indices[i*3 + n];
                                pointnum[p] = pointnum[p] | squarenum[i/2];
                        }
                }
                fp = fopen("roi.txt", "w");
                for (i = 0; i < n_points; i++)
                        fprintf(fp, " %d\n", pointnum[i]);
                fclose(fp);
        }

        free(points);
        free(indices);
        free(normals);

        if (flipflag) free(flip);

        free(squarenum);
        free(pointnum);

        return(object);
}

void
usage(char *executable)
{
        char *usage_str =
"\nUsage: %s [options] outfile.obj\n\n\
    Output a fractal surface (von Koch surface).  Specify the number of\n\
    subdivisions with the -level argument.\n\n";

        fprintf(stderr, usage_str, executable);
}

int
main(int argc, char** argv)
{
        char               *output_file;
        object_struct      **objects;
        polygons_struct    *polygons;
        File_formats       format;

        if (ParseArgv(&argc, argv, argTable, 0) || (argc < 2)) {
                usage(argv[0]);
                fprintf(stderr, "       %s -help\n\n", argv[0]);
                exit(EXIT_FAILURE);
        }

        initialize_argument_processing(argc, argv);
        if (!get_string_argument(NULL, &output_file)) {
                fprintf(stderr, "\nUsage: %s [options] outfile.obj\n", argv[0]);
                exit(EXIT_FAILURE);
        }

        if (squareflag) {
                objects = create_square_von_koch_surface(level);
        } else {
                objects = create_von_koch_surface(level);
        }

        if(output_graphics_any_format(output_file, ASCII_FORMAT,
                                      1, objects) != OK)
                exit(EXIT_FAILURE);

        delete_object_list(1, objects);

        return(EXIT_SUCCESS);
}

