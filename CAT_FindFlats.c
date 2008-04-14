/* Rachel Yotter - rachel.yotter@uni-jena.de                                 */
/* Department of Psychiatry                                                  */
/* University of Jena                                                        */
/*                                                                           */
/* Copyright Rachel Yotter, University of Jena.                              */


#include  <volume_io/internal_volume_io.h>
#include  <bicpl.h>
#include  <CAT_Blur2d.h>
#include  <CAT_Curvature.h>

#define DISTANCE 3.0
#define PINF 1.7976931348623157e+308 /* for doubles */
#define NINF -1.7976931348623157e+308 /* for doubles */


struct list {
    void *val;
    struct list *next;
};

private int maxi(double a[], int size) {
    int i, maxi;

    maxi = 0;
    for (i = 1; i < size; i++) {
        if (a[i] > a[maxi]) maxi = i;
    }
    return(maxi);
}

    

private void usage (STRING executable) {
    STRING   usage_str = "\n\
Usage: %s  object_file output_file [n_flats] [area] [curvtype]\n\n\
     Using mean curvature, CAT_FindFlats finds a few regions of points where the\n\
     curvature is close to zero. Prints the point number, the area of minimal\n\
     flatness, and the average curvature values.  Generates a colormap file for\n\
     viewing the flats.  Can specify the number of flats to look at [default: 2,\n\
     range is 1-10], the approximate area of the flat in mm2, and the type of\n\
     curvature calculation [default: 0].\n\n\
     curvtype:  0 - mean curvature (averaged over 3mm, in degrees)\n\
                1 - gaussian curvature\n\
                2 - curvedness\n\
                3 - shape index\n\
                4 - mean curvature (in radians)\n\n";

    print_error( usage_str, executable );
}

int main (int argc, char *argv[]) {
    STRING               ifname, ofname;
    FILE                 *file;
    File_formats         format;
    int                  i, vertex, flat, curvidx;
    int                  n_objects, n_flats, total_flats, curv_type;
    int                  *n_neighbours, **neighbours, *flats, *flatmap;
    object_struct        **objects;
    polygons_struct      *polygons;
    double               *curvatures, *curvscore;
    progress_struct      progress;

    initialize_argument_processing(argc, argv);

    if (!get_string_argument(NULL, &ifname) ||
            !get_string_argument(NULL, &ofname) ) {
        usage(argv[0]);
        return(1);
    }

    get_int_argument(2, &total_flats);
    get_int_argument(0, &curv_type);

    if (total_flats < 1 || total_flats > 10) {
        print("Out of bounds value for the number of flats.\n");
        return(1);
    }

    if (input_graphics_file(ifname, &format, &n_objects, &objects) != OK) {
        return(1);
    }

    if (n_objects != 1 || get_object_type(objects[0]) != POLYGONS) {
        print("File must contain 1 polygons object.\n");
        return(1);
    }

    polygons = get_polygons_ptr(objects[0]);
    compute_polygon_normals(polygons);
    get_all_polygon_point_neighbours(polygons, &n_neighbours, &neighbours);

    ALLOC(curvatures, polygons->n_points);

    initialize_progress_report(&progress, FALSE, polygons->n_items,
                                "Finding Regions of Low Curvature" );
    
    get_polygon_vertex_curvatures_cg(polygons, n_neighbours, neighbours,
                                   DISTANCE, curv_type, curvatures);
    /* use absolute value */
    for (i=0; i < polygons->n_points; i++) curvatures[i] = fabs(curvatures[i]);

    n_flats = 100;
    if (n_flats > polygons->n_points) n_flats = polygons->n_points;

    ALLOC(flats,n_flats);
    ALLOC(curvscore,n_flats);
    ALLOC(flatmap,polygons->n_points);

    curvidx = 0;
    for (i=0; i < n_flats; i++) curvscore[i] = PINF;
    

    /* get the points with the smallest curvature values */
    for (vertex=0; vertex < polygons->n_points; vertex++) {
        if (curvatures[vertex] < curvscore[curvidx]) {
            curvscore[curvidx] = curvatures[vertex];
            flats[curvidx] = vertex;
            curvidx = maxi(curvscore, n_flats);
        }
    } /* for vertex */

    /* calculate their "curvature scores" based on nearest neighbors */
    for (vertex=0; vertex < n_flats; vertex++) {
        for (i=0; i < n_neighbours[vertex]; i++) {
            curvscore[i] += curvatures[neighbours[vertex][i]];
        }
    }

    for (i=0; i < polygons->n_points; i++) flatmap[i] = 0;

    /* find the N flattest spots, color them */
    for (flat=0; flat < total_flats; flat++) {
        curvidx = 0;
        for (i=1; i < n_flats; i++) {
            if (flatmap[flats[i]] == 0 && curvscore[i] < curvscore[curvidx]) {
                curvidx = i;
            }
        }
        print("Vertex %d, curvscore %f\n",flats[curvidx], curvscore[curvidx]);
        flatmap[flats[curvidx]] = 1;
        for (i=0; i < n_neighbours[flats[curvidx]]; i++) {
            flatmap[neighbours[flats[curvidx]][i]] = 1;
        }
    }

    if (open_file(ofname, WRITE_FILE, ASCII_FORMAT, &file) != OK) {
        return(1);
    }

    for (i=0; i < polygons->n_points; i++) {
        if (output_int(file, flatmap[i]) != OK ||
            output_newline(file) != OK) {
            break;
        }
    }

    close_file(file);

    delete_object_list(n_objects, objects);
    FREE(curvatures);
    
    return(0);
}
