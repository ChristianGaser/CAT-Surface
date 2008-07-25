#include "topologysimplification/utils.h"
#include "topologysimplification/SurfaceObj.h"
#include "topologysimplification/FastLoop.h"
#include "topologysimplification/PatchDisk.h"

extern "C"
{
        #include  <bicpl.h>

        #include  "CAT_SurfaceIO.h"
}

void
usage(char* executable)
{
        char*   usage_str = "\n\
Usage: %s  input.obj output.obj\n\
    Correct topology of a surface to a sphere (genus zero) by using the\n\
    genetic approach of Segonne et al.\n\
\n\n";

        fprintf(stderr, usage_str, executable);
}

int
main(int argc, char *argv[])
{
        char                 *input_file, *output_file;
        File_formats         format;
        int                  status, n_objects;
        object_struct        **objects;

        initialize_argument_processing(argc, argv);

        if (!get_string_argument(NULL, &input_file) ||
            !get_string_argument(NULL, &output_file)) {
                usage(argv[0]);
                return(1);
        }

        if (input_graphics_any_format(input_file, &format, &n_objects,
                                      &objects) != OK) {
                fprintf(stderr, "Error while reading %s.\n", input_file);
                return(1);
        }
    
        // create a surfaace
        Surface surface=Surface();
        
        surface.ObjToSurface(objects[0]);
        
        int euler = surface.GetEuler();
        
        cout << "Euler number is " << surface.euler << " ( = " <<
                surface.nvertices << " - " << surface.nedges <<
                " + " << surface.nfaces << " ) " << endl;        

        // the patching disks
        surface.disk = new PatchDisk[4];
        for (int n = 0; n < 4; n++)
                surface.disk[n].Create(n);

        surface.InitSurface();
        surface.CorrectTopology();

        objects[0] = surface.SurfaceToObj();

        if (output_graphics_any_format(output_file, format, 1, objects) != OK)
                return(1);
    
        delete_object_list(n_objects, objects);

        return 0;
}
