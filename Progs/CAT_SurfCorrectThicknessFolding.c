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

#include "CAT_SurfaceIO.h"
#include "CAT_Math.h"
#include "CAT_CorrectThicknessFolding.h"

double max_dist = 6.0;           /* maximal thickness */
double slope = 0.0;              /* thickness-dependent correction slope */

/* the argument table */
static ArgvInfo argTable[] = {
  {"-max", ARGV_FLOAT, (char *) TRUE, (char *) &max_dist, "Define maximum thickness, where all values exceeding that will be cut."},
    {"-slope", ARGV_FLOAT, (char *) TRUE, (char *) &slope, "Linear weighting slope: smaller thickness values are corrected less strongly than larger thickness values. 0 disables weighting."},
  { NULL, ARGV_END, NULL, NULL, NULL }
};

void
usage(char *executable)
{
    char *usage_str = "\n\
Usage: %s   [options] surface_file thickness_file output_thickness_file\n\n\
   This program corrects cortical thickness that is influenced by folding using\n\
   the ideas from that paper:\n\
   Nagehan Demirci, Timothy S. Coalson, Maria A. Holland, David C. Van Essen, Matthew F. Glasser\n\
   Compensating Cortical Thickness for Cortical Folding-Related Variation\n\
   https://doi.org/10.1101/2025.05.03.651968.\n";

    fprintf(stderr, usage_str, executable);
}
  
int
main(int argc, char *argv[])
{
    double *thickness;
    int n_objects, n_vals;
    char *surface_file, *output_thickness_file, *thickness_file;
    File_formats format;
    object_struct **objects;
    polygons_struct *polygons;

    /* Call ParseArgv */
    if (ParseArgv(&argc, argv, argTable, 0)) {
        usage(argv[0]);
        fprintf(stderr, "       %s -help\n\n", argv[0]);
        exit(EXIT_FAILURE);
    }

    initialize_argument_processing(argc, argv);

    if (!get_string_argument( NULL, &surface_file) ||
        !get_string_argument( NULL, &thickness_file) ||
        !get_string_argument( NULL, &output_thickness_file)) {
          usage(argv[0]);
          exit(EXIT_FAILURE);
    }

    if (input_graphics_any_format(surface_file, &format, &n_objects,
                    &objects) != OK)
        exit(EXIT_FAILURE);

    if (input_values_any_format(thickness_file, &n_vals, &thickness) != OK)
        exit(EXIT_FAILURE);

    polygons = get_polygons_ptr(objects[0]);

    if (CAT_CorrectThicknessFoldingWeighted(polygons, n_vals, thickness, slope) != OK)
        exit(EXIT_FAILURE);

    clip_data(thickness, n_vals, 1E-10, max_dist, DT_FLOAT64); 
    
    output_values_any_format(output_thickness_file, n_vals,
                 thickness, TYPE_DOUBLE);

    delete_object_list(n_objects, objects);
    free(thickness);

    return(OK);

}
