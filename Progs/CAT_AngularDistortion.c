/* Rachel Yotter - rachel.yotter@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Rachel Yotter, University of Jena.
 * $Id$
 *
 */

#include <bicpl.h>
#include <ParseArgv.h>

#include "CAT_Surf.h"
#include "CAT_SurfaceIO.h"

BOOLEAN PerPoly = 0; /* by default, generate angle distortion per point */

/* the argument table */
ArgvInfo argTable[] = {
  { "-polygon", ARGV_CONSTANT, (char *) 1,
  (char *) &PerPoly,
  "Calculate angular distortion on a per-polygon basis." },

  { NULL, ARGV_END, NULL, NULL, NULL }
};


/*
 * Calculation of angular distortion based on Ju et al 2005.
 * All calculations are in radians.
 */
double
angular_distortion(float *a1, float *b1, float *c1,
           float *a2, float *b2, float *c2)
{
    double ab1, ac1, bc1, ab2, ac2, bc2;
    double A1, B1, C1, A2, B2, C2;

    ab1 = sqrt(pow(b1[0] - a1[0], 2) +
           pow(b1[1] - a1[1], 2) +
           pow(b1[2] - a1[2], 2));
    ac1 = sqrt(pow(c1[0] - a1[0], 2) +
           pow(c1[1] - a1[1], 2) +
           pow(c1[2] - a1[2], 2));
    bc1 = sqrt(pow(c1[0] - b1[0], 2) +
           pow(c1[1] - b1[1], 2) +
           pow(c1[2] - b1[2], 2));
    ab2 = sqrt(pow(b2[0] - a2[0], 2) +
           pow(b2[1] - a2[1], 2) +
           pow(b2[2] - a2[2], 2));
    ac2 = sqrt(pow(c2[0] - a2[0], 2) +
           pow(c2[1] - a2[1], 2) +
           pow(c2[2] - a2[2], 2));
    bc2 = sqrt(pow(c2[0] - b2[0], 2) +
           pow(c2[1] - b2[1], 2) +
           pow(c2[2] - b2[2], 2));

    if (ab1 > ac1 && ab1 > bc1) {
        C1 = acos( (bc1*bc1 + ac1*ac1 - ab1*ab1) / (2*bc1*ac1) );
        A1 = asin(bc1 * sin(C1)/ab1);
        B1 = PI - A1 - C1;
    } else if (ac1 > ab1 && ac1 > bc1) {
        B1 = acos( (bc1*bc1 + ab1*ab1 - ac1*ac1) / (2*bc1*ab1) );
        A1 = asin(bc1 * sin(B1)/ac1);
        C1 = PI - A1 - B1;
    } else {
        A1 = acos( (ab1*ab1 + ac1*ac1 - bc1*bc1) / (2*ab1*ac1) );
        B1 = asin(ac1 * sin(A1)/bc1);
        C1 = PI - A1 - B1;
    }

    if (ab2 > ac2 && ab2 > bc2) {
        C2 = acos( (bc2*bc2 + ac2*ac2 - ab2*ab2) / (2*bc2*ac2) );
        A2 = asin(bc2 * sin(C2)/ab2);
        B2 = PI - A2 - C2;
    } else if (ac2 > ab2 && ac2 > bc2) {
        B2 = acos( (bc2*bc2 + ab2*ab2 - ac2*ac2) / (2*bc2*ab2) );
        A2 = asin(bc2 * sin(B2)/ac2);
        C2 = PI - A2 - B2;
    } else {
        A2 = acos( (ab2*ab2 + ac2*ac2 - bc2*bc2) / (2*ab2*ac2) );
        B2 = asin(ac2 * sin(A2)/bc2);
        C2 = PI - A2 - B2;
    }
    
    A1 *= RAD_TO_DEG;
    A2 *= RAD_TO_DEG;
    B1 *= RAD_TO_DEG;
    B2 *= RAD_TO_DEG;
    C1 *= RAD_TO_DEG;
    C2 *= RAD_TO_DEG;

    return(fabs(A2 - A1) + fabs(B2 - B1) + fabs(C2 - C1));
}


int
main(int argc, char *argv[])
{
    char         *object_file, *object2_file, *output_surface_file;
    FILE         *fp;
    File_formats     format;
    int          n_objects, n_obj;
    int          i, tp[3], p, size, poly;
    object_struct    **objects, **objects2;
    polygons_struct    *polygons, *polygons2;
    double         ad, total_distortion = 0;
    double         *ad_values;
    int          *n_polys;

    /* Call ParseArgv */
    if (ParseArgv(&argc, argv, argTable, 0) || (argc < 3)) {
        fprintf(stderr,"\nUsage: %s [options] surface_file surface_file2 output_values_file\n", argv[0]);
        fprintf( stderr,"\nCalculate angular distortion between two surfaces.\n");
        fprintf(stderr, "   %s -help\n\n", argv[0]);
        exit(EXIT_FAILURE);
    }

    initialize_argument_processing(argc, argv);

    if (!get_string_argument(NULL, &object_file) ||
      !get_string_argument(NULL, &object2_file) ||
      !get_string_argument(NULL, &output_surface_file)) {
        fprintf(stderr,
            "Usage: %s  object_file object_file2 output_surface_file\n",
            argv[0]);
        exit(EXIT_FAILURE);
    }

    if (input_graphics_any_format(object_file, &format,
                    &n_objects, &objects) != OK) {
        exit(EXIT_FAILURE);
    }

    if (n_objects != 1 || get_object_type(objects[0]) != POLYGONS) {
        printf("File must contain 1 polygons object.\n");
        exit(EXIT_FAILURE);
    }

    if (input_graphics_any_format(object2_file, &format,
                    &n_objects, &objects2) != OK) {
        exit(EXIT_FAILURE);
    }

    if (n_objects != 1 || get_object_type(objects2[0]) != POLYGONS) {
        printf("File must contain 1 polygons object.\n");
        exit(EXIT_FAILURE);
    }

    fp = fopen(output_surface_file, "w");
    if (!fp) exit(EXIT_FAILURE);

    polygons = get_polygons_ptr(objects[0]);
    polygons2 = get_polygons_ptr(objects2[0]);

    if (polygons->n_items != polygons2->n_items ||
      polygons->n_points != polygons2->n_points) {
        fprintf(stderr, "Input polygons don't match. Exiting.\n");
        exit(EXIT_FAILURE);
    }

    if (PerPoly) {
        ALLOC(ad_values, polygons->n_items);
        n_obj = polygons->n_items;
    } else {    
        ALLOC(ad_values, polygons->n_points);
        ALLOC(n_polys, polygons->n_points);
        for (i = 0; i < polygons->n_points; i++) {
            n_polys[i] = 0;
            ad_values[i] = 0;
        }
        n_obj = polygons->n_points;
    }

    /* walk through facets */
    for (poly = 0; poly < polygons->n_items; poly++) {
        size = GET_OBJECT_SIZE(*polygons, poly);
        if (size != 3) {
            printf("Mesh must only contain triangles. Exiting..\n");
            exit(EXIT_FAILURE);
        }

        /* walk through facet points */
        for (i = 0; i < size; i++)
            tp[i] = polygons->indices[
                   POINT_INDEX(polygons->end_indices, poly, i)];

        ad = angular_distortion(polygons->points[tp[0]].coords,
                    polygons->points[tp[1]].coords,
                    polygons->points[tp[2]].coords,
                    polygons2->points[tp[0]].coords,
                    polygons2->points[tp[1]].coords,
                    polygons2->points[tp[2]].coords);

        if (!isnan(ad))
            total_distortion += ad;

        if (PerPoly) {
            ad_values[poly] = ad;
        } else {
            for (i = 0; i < size; i++) {
                n_polys[tp[i]]++;
                ad_values[tp[i]] += ad;
            }
        }
    }

    for (i = 0; i < n_obj; i++) {
        if (!PerPoly && n_polys[i] > 0) /* average distortion */
            ad_values[i] /= n_polys[i];

        if ((fprintf(fp, " %g", ad_values[i]) <= 0 ) ||
          (fprintf(fp, "\n" ) <= 0))
            exit(EXIT_FAILURE);
    }

    printf("Angular distortion: %f\n", total_distortion /
                       (3 * polygons->n_items));
    
    fclose(fp);

    delete_object_list(n_objects, objects);
    delete_object_list(n_objects, objects2);

    FREE(ad_values);
    if (!PerPoly)
        FREE(n_polys);

    return(EXIT_SUCCESS);
}
