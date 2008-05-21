/* Rachel Yotter - rachel.yotter@uni-jena.de                                 */
/* Department of Psychiatry                                                  */
/* University of Jena                                                        */
/*                                                                           */
/* Copyright Rachel Yotter, University of Jena.                              */

#include  <volume_io/internal_volume_io.h>
#include  <bicpl.h>
#include  <CAT_Surf.h>

/**
 * Calculation of angular distortion based on Ju et al 2005.
 * All calculations are in radians.
 */
double angular_distortion(float *a1, float *b1, float *c1,
                          float *a2, float *b2, float *c2) {
    double ab1, ac1, bc1, ab2, ac2, bc2;
    double A1, B1, C1, A2, B2, C2;

    ab1 = sqrt(pow(b1[0] - a1[0], 2)
             + pow(b1[1] - a1[1], 2)
             + pow(b1[2] - a1[2], 2));
    ac1 = sqrt(pow(c1[0] - a1[0], 2)
             + pow(c1[1] - a1[1], 2)
             + pow(c1[2] - a1[2], 2));
    bc1 = sqrt(pow(c1[0] - b1[0], 2)
             + pow(c1[1] - b1[1], 2)
             + pow(c1[2] - b1[2], 2));
    ab2 = sqrt(pow(b2[0] - a2[0], 2)
             + pow(b2[1] - a2[1], 2)
             + pow(b2[2] - a2[2], 2));
    ac2 = sqrt(pow(c2[0] - a2[0], 2)
             + pow(c2[1] - a2[1], 2)
             + pow(c2[2] - a2[2], 2));
    bc2 = sqrt(pow(c2[0] - b2[0], 2)
             + pow(c2[1] - b2[1], 2)
             + pow(c2[2] - b2[2], 2));

    if (ab1 > ac1 && ab1 > bc1) {
        C1 = acos( (bc1*bc1 + ac1*ac1 - ab1*ab1) / (2*bc1*ac1) );
        A1 = asin(bc1*sin(C1)/ab1);
        B1 = PI - A1 - C1;
    } else if (ac1 > ab1 && ac1 > bc1) {
        B1 = acos( (bc1*bc1 + ab1*ab1 - ac1*ac1) / (2*bc1*ab1) );
        A1 = asin(bc1*sin(B1)/ac1);
        C1 = PI - A1 - B1;
    } else {
        A1 = acos( (ab1*ab1 + ac1*ac1 - bc1*bc1) / (2*ab1*ac1) );
        B1 = asin(ac1*sin(A1)/bc1);
        C1 = PI - A1 - B1;
    }

    if (ab2 > ac2 && ab2 > bc2) {
        C2 = acos( (bc2*bc2 + ac2*ac2 - ab2*ab2) / (2*bc2*ac2) );
        A2 = asin(bc2*sin(C2)/ab2);
        B2 = PI - A2 - C2;
    } else if (ac2 > ab2 && ac2 > bc2) {
        B2 = acos( (bc2*bc2 + ab2*ab2 - ac2*ac2) / (2*bc2*ab2) );
        A2 = asin(bc2*sin(B2)/ac2);
        C2 = PI - A2 - B2;
    } else {
        A2 = acos( (ab2*ab2 + ac2*ac2 - bc2*bc2) / (2*ab2*ac2) );
        B2 = asin(ac2*sin(A2)/bc2);
        C2 = PI - A2 - B2;
    }

    return(fabs(A2 - A1) + fabs(B2 - B1) + fabs(C2 - C1));
}


int main(int argc, char *argv[]) {
    STRING               iob_fname1, iob_fname2, out_fname;
    FILE                 *file;
    File_formats         format;
    int                  n_objects1, n_objects2;
    int                  pnt, tp[3], p, size, fac;
    object_struct        **objects1, **objects2;
    polygons_struct      *polygons1, *polygons2;
    double               ang_dist, total_ang_dist = 0;

    initialize_argument_processing(argc, argv);

    if (!get_string_argument(NULL, &iob_fname1) ||
            !get_string_argument(NULL, &iob_fname2) ||
            !get_string_argument(NULL, &out_fname)) {
        print_error(
              "Usage: %s  object_file object_file2 output_file\n",
              argv[0] );
        return(1);
    }

    if (input_graphics_any_format(iob_fname1, &format, &n_objects1, &objects1) != OK) {
        return(1);
    }

    if (n_objects1 != 1 || get_object_type(objects1[0]) != POLYGONS) {
        print("File must contain 1 polygons object.\n");
        return(1);
    }

    if (input_graphics_any_format(iob_fname2, &format, &n_objects2, &objects2) != OK) {
        return(1);
    }

    if (n_objects2 != 1 || get_object_type(objects2[0]) != POLYGONS) {
        print("File must contain 1 polygons object.\n");
        return(1);
    }

    if (open_file(out_fname, WRITE_FILE, ASCII_FORMAT, &file) != OK) {
        return(1);
    }

    polygons1 = get_polygons_ptr(objects1[0]);
    polygons2 = get_polygons_ptr(objects2[0]);

    if (polygons1->n_items != polygons2->n_items) {
        print("The number of facets must be the same for both objects.\n");
        return(1);
    }

    for (fac = 0; fac < polygons1->n_items; fac++) { /* walk through facets */
        size = GET_OBJECT_SIZE(*polygons1, fac);
        if (size != GET_OBJECT_SIZE(*polygons2, fac)) {
            print("Objects don't match. Exiting...\n");
            return(1);
        }
        if (size != 3) {
            print("Mesh must only contain triangles. Exiting...\n");
            return(1);
        }

        for (pnt = 0; pnt < size; pnt++) { /* walk through facet points */
            int p1 = polygons1->indices[POINT_INDEX(polygons1->end_indices, fac, pnt)];
            int p2 = polygons2->indices[POINT_INDEX(polygons2->end_indices, fac, pnt)];
            if (p1 != p2) {
                print("Objects don't match. Exiting...\n");
                return(1);
            }
            tp[pnt] = p1; /* save facet points */
        }

        ang_dist = angular_distortion(polygons1->points[tp[0]].coords,
                                      polygons1->points[tp[1]].coords,
                                      polygons1->points[tp[2]].coords,
                                      polygons2->points[tp[0]].coords,
                                      polygons2->points[tp[1]].coords,
                                      polygons2->points[tp[2]].coords);
        total_ang_dist += ang_dist;
        if (output_double(file, ang_dist) != OK || output_newline(file) != OK) {
            return(1);
        }
    }

    print("Angular distortion: %f\n", total_ang_dist/(3*polygons1->n_items));
      
    close_file(file);

    delete_object_list(n_objects1, objects1);
    delete_object_list(n_objects2, objects2);

    return(0);
}

