/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 *
 */

#include <volume_io/internal_volume_io.h>
#include <bicpl.h>
#include <float.h>
#include <ParseArgv.h>

#include "CAT_Map2d.h"
#include "CAT_Curvature.h"
#include "CAT_SurfaceIO.h"

#define  BINTREE_FACTOR   0.5

double fwhm         = 0.0;
char * values_file  = NULL;
char * surface_file = NULL;
char * sphere_file  = NULL;
int    sz_map[2]    = {512, 256};
int    curvtype     = 0;

static ArgvInfo argTable[] = {
  {"-values", ARGV_STRING, (char *) 1, (char *) &values_file, 
     "Optional file with values for mapping."},
  {"-surf", ARGV_STRING, (char *) 1, (char *) &surface_file, 
     "Surface for mapping."},
  {"-sphere", ARGV_STRING, (char *) 1, (char *) &sphere_file, 
     "Sphere for surface."},
  {"-type", ARGV_INT, (char *) 1, (char *) &curvtype,
     "Curvature type\n\t0 - mean curvature (averaged over 3mm, in degrees)\n\t1 - gaussian curvature\n\t2 - curvedness\n\t3 - shape index\n\t4 - mean curvature (in radians)."},
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
Usage: %s -surf surface.obj|-values values_file [-sphere sphere_file] output.pgm\n\n\
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
        polygons_struct      *polygons, *sphere;
        int                  i, n_objects, x, y;
        int                  n_values;
        int                  *n_neighbours, **neighbours;
        double               *values, *data, mn, mx, distance;
        Point                centre;
        object_struct        **objects, *object;

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
                exit(EXIT_FAILURE);
        }

        if (sphere_file != NULL) {
                if (input_graphics_any_format(sphere_file, &format, &n_objects, &objects) != OK)
                        exit(EXIT_FAILURE);
                /* get a pointer to the sphere */
                sphere = get_polygons_ptr(objects[0]);
        } 

        if (surface_file != NULL) {
                if (input_graphics_any_format(surface_file, &format,
                                              &n_objects, &objects) != OK)
                        exit(EXIT_FAILURE);

                /* check that the surface file contains a polyhedron */
                if (n_objects != 1 || get_object_type(objects[0]) != POLYGONS) {
                        printf("Surface file must contain 1 polygon object.\n");
                        exit(EXIT_FAILURE);
                }
                /* get a pointer to the surface */
                polygons = get_polygons_ptr(objects[0]);
        }
        
        if (values_file == NULL) {
                create_polygon_point_neighbours(polygons, TRUE, &n_neighbours,
                                                &neighbours, FALSE, NULL);
                values = (double *) malloc(sizeof(double) * polygons->n_points);
                if (curvtype == 0)
                        distance = 3.0;
                else distance = 0.0;
                get_polygon_vertex_curvatures_cg(polygons,
                                              n_neighbours, neighbours,
                                              distance, curvtype, values);
        } else {
                if (input_values_any_format(values_file, &n_values, &values) != OK)
                        exit(EXIT_FAILURE);
                if (sphere_file == NULL) {
                        /* if no surface_file and sphere_file is given then create tetra */
                        object = create_object(POLYGONS);
                        polygons = get_polygons_ptr(object);
                        fill_Point(centre, 0.0, 0.0, 0.0);
                        create_tetrahedral_sphere(&centre, 1.0, 1.0, 1.0, (2 * (n_values-2)),
                                  sphere);   
                } 
        }
        data = (double *) malloc(sizeof(double) * sz_map[0] * sz_map[1]);
        
        if (values_file == NULL)
                map_smoothed_curvature_to_sphere(polygons, sphere, values, data, fwhm, sz_map, curvtype);
        else
                map_smoothed_curvature_to_sphere(sphere, sphere, values, data, fwhm, sz_map, curvtype);

        /* scale data to uint8 range */
        mn = FLT_MAX; mx = -FLT_MAX;
        for (i = 0; i < sz_map[0]*sz_map[1]; i++) {
                if (data[i] > mx) mx = data[i];
                if (data[i] < mn) mn = data[i];
        }

        for (i = 0; i < sz_map[0]*sz_map[1]; i++) 
                data[i] = 255.0 * (data[i] - mn) / (mx - mn);

        if (write_pgm(output_file, data, sz_map[0], sz_map[1]) != 0)
                exit(EXIT_FAILURE);

        delete_object_list(n_objects, objects);
        free(data);

        if (values_file == NULL)
                free(values);

        return(EXIT_SUCCESS);
}
