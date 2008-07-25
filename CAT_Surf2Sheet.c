/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 */

#include <volume_io/internal_volume_io.h>
#include <bicpl.h>
#include <float.h>
#include <ParseArgv.h>

#include "CAT_SheetIO.h"
#include "CAT_Map2d.h"

#define  BINTREE_FACTOR   0.5

double fwhm         = 0.0;
double curv_dist    = 3.0;
char * values_file  = NULL;
char * surface_file = NULL;
int    sz_map[2]    = {512, 256};

static ArgvInfo argTable[] = {
  {"-values", ARGV_STRING, (char *) 1, (char *) &values_file, 
     "Optional file with values for mapping."},
  {"-surf", ARGV_STRING, (char *) 1, (char *) &surface_file, 
     "Surface for mapping."},
  {"-distance", ARGV_FLOAT, (char *) 1, (char *) &curv_dist,
     "Average curvature around distance."},
  {"-fwhm", ARGV_FLOAT, (char *) 1, (char *) &fwhm,
     "Filter size for curvature map in FWHM."},
  {"-sz", ARGV_INT, (char *) 2, (char *) &sz_map,
     "Size of PGM-map."},
   {NULL, ARGV_END, NULL, NULL, NULL}
};

void
usage(char *executable)
{
        char *usage_str = "\n\
Usage: %s -surf surface.obj|-values values_file output.pgm\n\n\
     Maps a surface to a flat sheet image. If values_file is not specified\n\
     then mean curvature is used as color. In this case you have to define\n\
     surface_file.  Sheet image will be saved as binary PGM-file (default\n\
     size 512 x 256).\n\n";

        fprintf(stderr, usage_str, executable);
}

int
main(int argc, char *argv[])
{
        char                 *output_file;
        File_formats         format;
        polygons_struct      *polygons;
        int                  i, n_objects, x, y;
        int                  n_values;
        int                  *n_neighbours, **neighbours;
        double               *values, *data, mn, mx;
        Point                centre;
        object_struct        **objects, *object;
        BOOLEAN              values_specified;

        initialize_argument_processing(argc, argv);

        if (ParseArgv(&argc, argv, argTable, 0) ||
            (values_file != NULL && surface_file != NULL) ||
            (values_file == NULL && surface_file == NULL)) {
                fprintf(stderr, "Either surface or value file must be defined.\n");
                usage(argv[0]);
                fprintf(stderr, "     %s -help\n\n", argv[0]);
                exit(EXIT_FAILURE);
        }

        if (!get_string_argument(NULL, &output_file)) {
                usage(argv[0]);
                return(1);
        }

        if (values_file == NULL) {
                values_specified = FALSE;
        } else {
                if (input_texture_values(values_file, &n_values, &values) != OK)
                        return(1);
                values_specified = TRUE;
        }

        if (surface_file == NULL) {
                /* if no surface_file is given then create tetra */
                object = create_object(POLYGONS);
                polygons = get_polygons_ptr(object);
                fill_Point(centre, 0.0, 0.0, 0.0);
                create_tetrahedral_sphere(&centre, 1, 1, 1, (2 * (n_values-2)),
                                          polygons);
        } else {
                if (values_specified) {
                        if (input_texture_values(values_file, &n_values,
                                                 &values) != OK)
                                return(1);
                }
                if (input_graphics_any_format(surface_file, &format,
                                              &n_objects, &objects) != OK)
                        return(1);

                /* check that the surface file contains a polyhedron */
                if (n_objects != 1 || get_object_type(objects[0]) != POLYGONS) {
                        printf("Surface file must contain 1 polygon object.\n");
                        return(1);
                }
                /* get a pointer to the surface */
                polygons = get_polygons_ptr(objects[0]);
        }

        if (!values_specified) {
                create_polygon_point_neighbours(polygons, TRUE, &n_neighbours,
                                                &neighbours, FALSE, NULL);
                values = (double *) malloc(sizeof(double) * polygons->n_points);
                get_polygon_vertex_curvatures(polygons,
                                              n_neighbours, neighbours,
                                              curv_dist, 0.0, values);
        }

        data = (double *) malloc(sizeof(double) * sz_map[0] * sz_map[1]);
        map_smoothed_curvature_to_sphere(polygons, values, data, fwhm, sz_map);

        /* scale data to uint8 range */
        mn = FLT_MAX; mx = -FLT_MAX;
        for (i = 0; i < sz_map[0]*sz_map[1]; i++) {
                if (data[i] > mx) mx = data[i];
                if (data[i] < mn) mn = data[i];
        }

        for (i = 0; i < sz_map[0]*sz_map[1]; i++) 
                data[i] = 255.0 * (data[i] - mn) / (mx - mn);

        if (write_pgm(output_file, data, sz_map[0], sz_map[1]) != 0)
                return(1);

        delete_object_list(n_objects, objects);
        free(data);

        if (!values_specified)
                free(values);

        return(0);
}
