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

/* argument defaults */
int   check_intersect = 0;
int   equivol = 0;
double weight = 1.0;
char *label_file  = NULL;

/* the argument table */
static ArgvInfo argTable[] = {
  {"-label", ARGV_STRING, (char *) 1, (char *) &label_file, 
     "Additional label file for estimating  whether inner or outer border is reached."},
  {"-equivolume", ARGV_CONSTANT, (char *) TRUE, (char *) &equivol,
   "Use equi-volume model by Bok (1929) to correct distances/layers. The correction is based on Waehnert et al. (2014)."},
  {"-weight", ARGV_FLOAT, (char *) TRUE, (char *) &weight,
   "Weight between equi-volume model (1.0) and original data (0.0)."},
  {"-check_intersect", ARGV_CONSTANT, (char *) TRUE, (char *) &check_intersect,
   "Correct self intersections"},
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
    double *thickness_values, *area_inner, *area_outer, *extents;
    double value, surface_area, pos;
    float *labels, extent;
    Status status;
    char *src_file, *out_file, *values_file;
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
      !get_string_argument( NULL, &out_file)) {
        usage(argv[0]);
        fprintf(stderr, "   %s -help\n\n", argv[0]);
        exit(EXIT_FAILURE);
    }

    get_real_argument(0.5, &extent);
    
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

    extents = (double *) malloc(sizeof(double) * polygons->n_points);
        
    if ((equivol) && (extent > -0.5) && (extent < 0.5)) {
        area_inner = (double *) malloc(sizeof(double) * polygons->n_points);
        area_outer = (double *) malloc(sizeof(double) * polygons->n_points);
        
        /* pial (outer) surface */
        surface_area = get_area_of_points_central_to_pial(polygons, area_outer, thickness_values, 0.5);

        /* white (inner) surface */
        surface_area = get_area_of_points_central_to_pial(polygons, area_inner, thickness_values, -0.5);
        
        /* get relative position inside cortical band */
        pos = extent + 0.5;
        
        for (p = 0; p < polygons->n_points; p++) {
            /* no change needed if inner and outer area are equal */
            if (area_outer[p] == area_inner[p])
                extents[p] = extent;
            else {
                /* eq. 10 from Waehnert et al. 2014 */
                extents[p] = (1.0/(area_outer[p] - area_inner[p]))*
                       (sqrt((pos*area_outer[p]*area_outer[p]) + ((1.0 - pos)*area_inner[p]*area_inner[p])) - area_inner[p]);
                extents[p] -= 0.5; /* subtract offset of 0.5 that was added to pos */      
                extents[p] = weight*extents[p] + (1.0 - weight)*extent;        
            }        
        }
        free(area_inner);
        free(area_outer);
    } else {
        for (p = 0; p < polygons->n_points; p++)
            extents[p] = extent;
    }
    
    if (label_file != NULL) {
        nii_ptr = read_nifti_float(label_file, &labels, 0);
        if (!nii_ptr) {
            fprintf(stderr, "Error reading %s.\n", label_file);
            return (EXIT_FAILURE);
        }
    } else {
        nii_ptr = NULL;
        labels = NULL;
    }
    
    objects_out = central_to_new_pial(polygons, thickness_values, extents, labels, nii_ptr, check_intersect);

    if(output_graphics_any_format(out_file, format, 1, objects_out, NULL) != OK)
        exit(EXIT_FAILURE);
        
    free(extents);
    
    return(status != OK);
}
