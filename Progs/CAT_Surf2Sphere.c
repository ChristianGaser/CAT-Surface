/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * most of the code is modified from
 * caret/caret_brain_set/BrainModelSurface.cxx.
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 *
 */

#include <bicpl.h>
#include <float.h>
#include <ParseArgv.h>

#include "CAT_Surf.h"
#include "CAT_SurfaceIO.h"
#include "CAT_SphereRelax.h"
#include "CAT_SphereUntangle.h"
#include "CAT_ArgvCheck.h"

int relax = 0;
int keep_folds = 0;

static ArgvInfo argTable[] = {
    {"-relax", ARGV_CONSTANT, (char *)TRUE, (char *)&relax,
     "Reduce the metric distortion of the sphere afterwards, as FreeSurfer's\n\t"
     "mris_sphere does: distances between neighbours (2 rings) are pulled\n\t"
     "towards their distances on the input surface, with the gradient\n\t"
     "averaged over 1024 down to 0 neighbour iterations. The inflation leaves\n\t"
     "the poles compressed and the medial wall expanded, which spherical\n\t"
     "registration inherits; use it for spheres that are registered. With it\n\t"
     "stop_at 6 is enough (10 gives the same). Needs stop_at 5 or more."},
    {"-keep-folds", ARGV_CONSTANT, (char *)TRUE, (char *)&keep_folds,
     "Do not remove the folded (negative-area) triangles of the sphere. By\n\t"
     "default every cluster of folded triangles is re-embedded on a disk of\n\t"
     "the mesh around it, as the last step once the surface is a sphere\n\t"
     "(stop_at 5 or more, after -relax); a sphere without folds is not\n\t"
     "changed."},
    {NULL, ARGV_END, NULL, NULL, NULL}};

void
usage(char *executable)
{
    char *usage_str = "\n\
Usage: %s [-relax] [-keep-folds] surface_file output_surface_file [stop_at] [verbose]\n\n\
   Maps a surface to a sphere using the caret inflating approach.\n\n\
   The inflating can be limited using stop_at (default 5), where\n\
     1  - Low smooth\n\
     2  - Inflating\n\
     3  - Very inflating\n\
     4  - High smoothing\n\
     5  - Ellipsoid\n\n\
     >5 - Additional areal smoothing to obtain more equally sized meshes\n\
      1000 iterations are used for each stop_at value > 5, thus \n\
      stop_at of 6 means 1000 iterations, stop_at of 7 means 2000 iterations, etc.\n\
      The recommended stop_at value is 10 that corresponds to 5000 iterations.\n\n\
   -relax reduces the metric distortion of the sphere afterwards (see -help).\n\
   Folded triangles of the sphere are removed at the end unless -keep-folds\n\
   is given.\n";

    fprintf(stderr, usage_str, executable);
}

int
main(int argc, char *argv[])
{
    char *input_file, *output_surface_file;
    int n_objects, stop_at, verbose, unfold;
    File_formats format;
    object_struct **object_list;
    polygons_struct *polygons, original;

    if (ParseArgv(&argc, argv, argTable, 0))
    {
        usage(argv[0]);
        exit(EXIT_FAILURE);
    }

    /* a leftover that looks like an option is a typo or a removed flag;
       ParseArgv() keeps it, so it would be read as a file name below */
    if (cat_check_unknown_options(argc, argv, argv[0]))
    {
        usage(argv[0]);
        exit(EXIT_FAILURE);
    }

    initialize_argument_processing(argc, argv);
    if (!get_string_argument(NULL, &input_file) ||
      !get_string_argument(NULL, &output_surface_file)) {
        usage(argv[0]);
        exit(EXIT_FAILURE);
    }

    get_int_argument(5, &stop_at);
    get_int_argument(0, &verbose);

    if (relax && stop_at < 5)
    {
        fprintf(stderr, "-relax needs a sphere: stop_at must be 5 or more.\n");
        exit(EXIT_FAILURE);
    }

    if (input_graphics_any_format(input_file, &format, &n_objects,
                    &object_list) != OK || n_objects != 1 ||
      get_object_type(object_list[0]) != POLYGONS) {
        fprintf(stderr, "Error reading %s.\n", input_file);
        exit(EXIT_FAILURE);
    }

    polygons = get_polygons_ptr(object_list[0]);

    /* the relaxation and the fold removal compare the sphere with the
       surface it maps */
    unfold = !keep_folds && stop_at >= 5;
    if (relax || unfold)
        copy_polygons(polygons, &original);

    surf_to_sphere(polygons, stop_at, verbose);

    if (relax)
    {
        CAT_SphereRelaxOpts opts;

        CAT_SphereRelaxOptionsInit(&opts);
        opts.verbose = verbose;
        if (CAT_SphereRelaxMetric(&original, polygons, &opts) != 0)
            exit(EXIT_FAILURE);
    }

    if (unfold)
    {
        int left = CAT_SphereRemoveFolds(&original, polygons, verbose);

        if (left < 0)
            exit(EXIT_FAILURE);
        if (left > 0)
            fprintf(stderr, "Warning: %d folded or degenerate triangles are left.\n",
                    left);
    }
    if (relax || unfold)
        delete_polygons(&original);

    if(output_graphics_any_format(output_surface_file, format, 1,
            object_list, NULL) != OK)
        exit(EXIT_FAILURE);

    delete_object_list(n_objects, object_list);

    return(EXIT_SUCCESS);
}
