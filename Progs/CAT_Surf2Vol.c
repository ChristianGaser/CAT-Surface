/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * most of the code is modified from
 * caret/caret_brain_set/BrainModelSurface.cxx.
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id: CAT_Surf2Vol.c 330 2014-11-06 13:13:05Z gaser $
 *
 */

#include <bicpl.h>
#include <float.h>

#include "CAT_Surf.h"
#include "CAT_SurfaceIO.h"

void
usage(char *executable)
{
        char *usage_str = "\n\
Usage: %s surface_file output_volume_file\n\n\
     Maps a surface to a volume and fill voxel inside with ones.\n";

        fprintf(stderr, usage_str, executable);
}

int
main(int argc, char *argv[])
{
        char             *input_file, *output_volume_file;
        int              n_objects, i;
        File_formats     format;
        object_struct    **object_list;
        polygons_struct  *polygons;
        Volume           volume, label_volume;        
        int              x, y, z, range_changed[2][3], value[3], sizes[3];
        double           separations[3], voxel[3], world[3];
        double           bounds[6], label_val;
    
        initialize_argument_processing(argc, argv);

        if (!get_string_argument(NULL, &input_file) ||
            !get_string_argument(NULL, &output_volume_file)) {
                usage(argv[0]);
                exit(EXIT_FAILURE);
        }
    
        if (input_graphics_any_format(input_file, &format, &n_objects,
                                      &object_list) != OK || n_objects != 1 ||
            get_object_type(object_list[0]) != POLYGONS) {
                fprintf(stderr, "Error reading %s.\n", input_file);
                exit(EXIT_FAILURE);
        }

        polygons = get_polygons_ptr(object_list[0]);
        
        /* create volume inside surface bounds */
        get_bounds(polygons, bounds);
        volume = create_volume( 3, XYZ_dimension_names, NC_BYTE, FALSE,
                            0.0, 255.0 );
        
        /* prepare volume parameters */
        for (i=0; i<3; i++) {
                separations[i] = 0.75;
                voxel[i] = -0.5;
                world[i] = bounds[2*i] - separations[i];
                value[i] = 2.0;
                sizes[i] = ROUND((bounds[2*i+1] - bounds[2*i]) / separations[i]) + 2;
        }

        set_volume_separations( volume, separations );    
        set_volume_sizes( volume, sizes );
        set_volume_voxel_range( volume, 0.0, 255.0 );
        set_volume_real_range( volume, 0, 255.0 );
        set_volume_translation( volume, voxel, world );

        alloc_volume_data( volume );
        
        /* label volume according to surface */
        label_val = 1.0;
        label_volume = create_label_volume( volume, NC_BYTE );
        scan_object_to_volume( object_list[0], volume, label_volume, (int) label_val, 0.0 );
        
        /* fill inside volume */
        fill_connected_voxels( volume, label_volume, EIGHT_NEIGHBOURS,
                           value, 0, 0, 1, 0.0, -1.0, range_changed );
                           
        /* label is inverse, correct thgis */
        for (x = 0; x < sizes[0]; x++)
          for (y = 0; y < sizes[1]; y++)
            for (z = 0; z < sizes[2]; z++)
            set_volume_real_value( label_volume, x, y, z, 0, 0, (1 - get_volume_real_value( label_volume, x, y, z, 0, 0 )));
                

        output_volume(output_volume_file, NC_BYTE, 0, 0.0, 0.0, label_volume, "Surf2Vol\n", NULL); 

        delete_volume( volume );
        delete_volume( label_volume );

        return(EXIT_SUCCESS);
}
