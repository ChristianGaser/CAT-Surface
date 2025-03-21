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
#include "CAT_Vol.h"
#include "CAT_SurfaceIO.h"
#include "CAT_Smooth.h"
#include "CAT_Deform.h"

/* argument defaults */
int   check_intersect = 0;
int  verbose = 0; 

/* the argument table */
static ArgvInfo argTable[] = {
  {"-check_intersect", ARGV_CONSTANT, (char *) TRUE, (char *) &check_intersect,
   "Correct self intersections"},
  {"-verbose", ARGV_CONSTANT, (char *) TRUE, (char *) &verbose,
     "Be verbose."},
   {NULL, ARGV_END, NULL, NULL, NULL}
};

void
usage(char *executable)
{
    static char *usage_str = "\n\
Usage: %s  surface_file thickness_file output_surface_file [extent]\n\
Estimate pial or white surface from central surface using cortical thickness values. In order to estimate the pial surface an extent of 0.5 (default) should be used, while an extent of -0.5 results in the estimation of the white matter surface.\n\
The equi-volume model optionally allows to correct the position of the surface around gyri and sulci. The area of the inner (white) and outer (pial) surface is used for this correction.\n\
Furthermore, you can weight the extent of equi-volume correction which is helpful to correct the initial central surface in CAT12 in highly folded areas with high mean curvature.\n\n";

    fprintf(stderr, usage_str, executable);
}

int
main(int argc, char *argv[])
{
    int p, n_objects, n_values;
    double *thickness_values, *extents;
    double value, surface_area, pos;
    float *labels;
    Status status;
    char *src_file, *pial_file, *white_file, *values_file, *label_file;
    File_formats format;
    nifti_image *nii_ptr;
    polygons_struct *polygons;
    object_struct **object_list, **objects_out;

    /* get the arguments from the command line */
    if (ParseArgv(&argc, argv, argTable, 0)) {
        usage(argv[0]);
        fprintf(stderr, "   %s -help\n\n", argv[0]);
        exit(EXIT_FAILURE);
    }

    initialize_argument_processing(argc, argv);
    
    if (!get_string_argument( NULL, &src_file) ||
      !get_string_argument( NULL, &values_file) ||
      !get_string_argument( NULL, &label_file) ||
      !get_string_argument( NULL, &pial_file) ||
      !get_string_argument( NULL, &white_file)) {
        usage(argv[0]);
        fprintf(stderr, "   %s -help\n\n", argv[0]);
        exit(EXIT_FAILURE);
    }

    if (input_graphics_any_format(src_file, &format, &n_objects,
                    &object_list) != OK)
        exit(EXIT_FAILURE);

    if (input_values_any_format(values_file, &n_values, &thickness_values) != OK)
        exit(EXIT_FAILURE);
    
    if (n_objects > 1) {
        fprintf(stderr,"Only one object allowed.\n");
        exit(EXIT_FAILURE);
    }

    polygons = get_polygons_ptr(object_list[0]);

    if (polygons->n_points != n_values) {
        fprintf(stderr,"Number of points differs from number of values.\n");
        exit(EXIT_FAILURE);
    }

    nii_ptr = read_nifti_float(label_file, &labels, 0);
    if (!nii_ptr) {
        fprintf(stderr, "Error reading %s.\n", label_file);
        return (EXIT_FAILURE);
    }

    polygons_struct *polygons_pial, *polygons_white;
    object_struct *object_pial = create_object(POLYGONS);
    object_struct *object_white = create_object(POLYGONS);

    extents = (double *) malloc(sizeof(double) * polygons->n_points);

    for (p = 0; p < polygons->n_points; p++) extents[p] = 0.35;
    objects_out = central_to_new_pial(polygons, thickness_values, extents, NULL, NULL, 0, verbose);
    object_pial = objects_out[0];

    for (p = 0; p < polygons->n_points; p++) extents[p] = -0.4;
    objects_out = central_to_new_pial(polygons, thickness_values, extents, NULL, NULL, 0, verbose);
    object_white = objects_out[0];

    polygons_pial = get_polygons_ptr(object_pial);
    polygons_white = get_polygons_ptr(object_white);

    int it = 200;
    double sigma = 0.4;
    double w[4] = {0.2, 0.1, 0.25, 0.05};
    surf_deform_dual(polygons_pial, polygons_white, labels, nii_ptr, 
                      w, sigma, 0.001, 0.999, thickness_values, it, verbose);

    if(output_graphics_any_format(pial_file, format, 1, &object_pial, NULL) != OK)
        exit(EXIT_FAILURE);

    if(output_graphics_any_format(white_file, format, 1, &object_white, NULL) != OK)
        exit(EXIT_FAILURE);
            
    free(extents);
    free(object_pial);
    free(object_white);
    return(status != OK);
}
