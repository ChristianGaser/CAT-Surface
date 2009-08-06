/* Rachel Yotter - rachel.yotter@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id: CAT_DumpTopoErrors.c 89 2009-01-27 14:43:59Z raytrace $
 */

#include <volume_io/internal_volume_io.h>
#include <bicpl.h>
#include <ParseArgv.h>

#include "CAT_SheetIO.h"
#include "CAT_Map2d.h"
#include "CAT_Blur2d.h"
#include "CAT_SPH.h"
#include "CAT_Octree.h"

#define PINF  1.7976931348623157e+308 /* for doubles */
#define NINF -1.7976931348623157e+308 /* for doubles */


void
usage(char *executable)
{
        char *usage_str =
"\nUsage: %s [options] object_file output_file\n\n\
    Locate and mark topological errors\n\n";

        fprintf(stderr, usage_str, executable);
}

int
main(int argc, char *argv[])
{
        char                 *surface_file, *out_file;
        File_formats         format;
        int                  p, ndefects, *defects, n_objects;
        polygons_struct      *polygons;
        object_struct        **poly_objects, **objects;
        FILE                 *fp;

        initialize_argument_processing(argc, argv);

        if (!get_string_argument(NULL, &surface_file) ||
            !get_string_argument(NULL, &out_file)) {
                usage(argv[0]);
                return(1);
        }
     
        if (input_graphics_any_format(surface_file, &format,
                                      &n_objects, &poly_objects) != OK)
                return(1);

        /* check that the surface file contains a polyhedron */
        if (n_objects != 1 || get_object_type(poly_objects[0]) != POLYGONS) {
                printf("Surface file must contain 1 polygons object.\n");
                return(1);
        }
        /* get a pointer to the surface */
        polygons = get_polygons_ptr(poly_objects[0]);
    
        defects = (int *) malloc(sizeof(int) * polygons->n_points);
        ndefects = find_topological_defects(polygons, defects);
        expand_defects(polygons, defects, 0, 2);

        printf("%d errors found\n", ndefects);

        if (open_file(out_file, WRITE_FILE, ASCII_FORMAT, &fp) != OK)
                exit(0);

        for (p = 0; p < polygons->n_points; p++)
                fprintf(fp, " %d.0\n", defects[p]);
        fclose(fp);

        /* clean up */
        free(defects);
    
        return(0);
}
