/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 *
 */

#include <bicpl.h>
#include <ParseArgv.h>

#include "CAT_Surf.h"
#include "CAT_SurfaceIO.h"

enum { Tfs, Tlink};

/* argument defaults */
int  dist_func        = Tfs;     /* default distance function: Freesurfer method */
char *thickness_file  = NULL;    /* thickness file for estimating inner and outer surface from central surface */
char *position_file  = NULL;
int   check_intersect = 0;

/* the argument table */
static ArgvInfo argTable[] = {
  {"-thickness", ARGV_STRING, (char *) 1, (char *) &thickness_file, 
     "Additional thickness file for internally estimating inner and outer surface based on central surface and thickness."},
  {"-position", ARGV_STRING, (char *) 1, (char *) &position_file, 
     "Additional position file for estimating whether inner or outer border is reached."},
  {"-mean", ARGV_CONSTANT, (char *) Tfs, (char *) &dist_func,
     "Calculate mean of closest distance between surface 1 and 2 and vice versa (Tfs, default)." },
  {"-link", ARGV_CONSTANT, (char *) Tlink, (char *) &dist_func,
     "Calculate the linked (exact) distance between both surfaces (Tlink)." },
  {"-check_intersect", ARGV_CONSTANT, (char *) TRUE, (char *) &check_intersect,
     "Correct self intersections if you use thickness file for internally estimating inner and outer surface."},
  { NULL, ARGV_END, NULL, NULL, NULL }
};

void
usage(char *executable)
{
        char *usage_str = "\n\
Usage: CAT_SurfDistance [options] surface_file surface_file2 output_values_file\n\n\
   or: CAT_SurfDistance [options] -thickness thickness_file surface_file output_values_file\n\
     Calculate linked or Freesurfer distance between two surfaces with the same number of \n\
     vertices (i.e. inner and outer surface) on a point-by-point basis.\n\
     Instead of defining two surfaces you can also use the thickness flag and the central \n\
     surface to internally create inner and outer surface.\n\
     This function is primarily thought to estimate the Tfs-distance from Freesurfer for \n\
     already existing Tpbt-distances from CAT12.\n";

        fprintf(stderr, usage_str, executable);
}

int
main(int argc, char *argv[])
{
    int n_objects, n_values;
    int i;
    float *positions;
    double shift, max_distance = 0;
    double *extents, *distance, *thickness_values;
    char *object_file, *object2_file, *output_surface_file;
    FILE *fp;
    File_formats format;
    object_struct **objects, **objects2;
    polygons_struct *polygons, *polygons2;
    nifti_image *nii_ptr;

    /* Call ParseArgv */
    if (ParseArgv(&argc, argv, argTable, 0) ||
      ((argc < 4) && (thickness_file == NULL)) ||
      ((argc < 3) && (thickness_file != NULL))) {
        usage(argv[0]);
        fprintf(stderr, "       %s -help\n\n", argv[0]);
        exit(EXIT_FAILURE);
    }

    initialize_argument_processing(argc, argv);

    get_string_argument(NULL, &object_file);
    if (thickness_file == NULL) {
        get_string_argument(NULL, &object2_file);
        if (position_file != NULL) {
            fprintf(stderr, "Position file can only be used with thickness option.\n");
            return (EXIT_FAILURE);
        }
    } else {
        if (input_values_any_format(thickness_file, &n_values, &thickness_values) != OK)
            exit(EXIT_FAILURE);
    }
    get_string_argument(NULL, &output_surface_file);

    if (input_graphics_any_format(object_file, &format,
                                  &n_objects, &objects) != OK) {
        exit(EXIT_FAILURE);
    }

    if (n_objects != 1 || get_object_type(objects[0]) != POLYGONS) {
        printf("File must contain 1 polygons object.\n");
        exit(EXIT_FAILURE);
    }

    polygons  = get_polygons_ptr(objects[0]);

    if (thickness_file == NULL) { /* two surfaces defined */ 
        if (input_graphics_any_format(object2_file, &format,
                                      &n_objects, &objects2) != OK) {
            exit(EXIT_FAILURE);
        }
    } else { /* thickness flag and one surface defined */
        extents = (double *) malloc(sizeof(double) * polygons->n_points);

        if (position_file != NULL) {
            nii_ptr = read_nifti_float(position_file, &positions, 0);
            if (!nii_ptr) {
                fprintf(stderr, "Error reading %s.\n", position_file);
                return (EXIT_FAILURE);
            }
        } else {
            nii_ptr = NULL;
            positions = NULL;
        }
        
        /* If position file is defined we go a bit further than half thickness */
        shift = (position_file != NULL) ? 0.55 : 0.5;

        /* obtain pial surface */
        for (i = 0; i < polygons->n_points; i++) extents[i] = shift;
        objects2 = central_to_new_pial(polygons, thickness_values, extents, 
            positions, nii_ptr, check_intersect);
        
        /* obtain white surface */
        for (i = 0; i < polygons->n_points; i++) extents[i] = -shift;
        objects = central_to_new_pial(polygons, thickness_values, extents, 
            NULL, NULL, check_intersect);
        polygons = get_polygons_ptr(objects[0]);

        free(extents);
    }

    polygons2 = get_polygons_ptr(objects2[0]);
    
    if (n_objects != 1 || get_object_type(objects2[0]) != POLYGONS) {
        printf("File must contain 1 polygons object.\n");
        exit(EXIT_FAILURE);
    }

    if ((polygons->n_items != polygons2->n_items ||
                  polygons->n_points != polygons2->n_points)) {
        fprintf(stderr, "Number of vertices between surfaces don't match. Exiting.\n");
        exit(EXIT_FAILURE);
    }

    distance = (double *) malloc(sizeof(double) * polygons->n_points);

    if (dist_func == Tfs) /* mean of both distances */
        max_distance = compute_point_distance_mean(polygons, polygons2, distance, 0);
    else                  /* linked distance */
        max_distance = compute_point_distance(polygons, polygons2, distance, 0);

    if (output_values_any_format(output_surface_file, polygons->n_points,
                                 distance, TYPE_DOUBLE) != OK) {
        exit(EXIT_FAILURE);
    }

    delete_object_list(n_objects, objects);
    delete_object_list(n_objects, objects2);

    free(distance);
    return(EXIT_SUCCESS);
}
