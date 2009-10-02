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

#include "CAT_SurfaceIO.h"

#define GET_grid_POINT(result, grid_origin, normal, length) \
{ \
        (result)[X] = RPoint_x(grid_origin) + (length) * Vector_x(normal); \
        (result)[Y] = RPoint_y(grid_origin) + (length) * Vector_y(normal); \
        (result)[Z] = RPoint_z(grid_origin) + (length) * Vector_z(normal); \
}

#define F_AVERAGE 0
#define F_RANGE   1
#define F_MAX     2
#define F_MIN     3
#define F_EXP     4
#define F_SUM     5

#define LOG05     -0.69314718
#define MAX_N_ARRAY 250

/* argument defaults */
int  degrees_continuity = -1; /* interpolation - default: nearest neighbour */
int  map_func = F_MAX;        /* default mapping function: maximum value */
Real grid_res = 1.0;          /* resolution of grid along normals in mm */
Real grid_length = 10.0;      /* length of grid along normals */
Real grid_origin = 0.0;       /* origin of grid along normals */
Real range[2] = {FLT_MAX, FLT_MAX};
Real exp_half = FLT_MAX;

/* the argument table */
ArgvInfo argTable[] = {
   {"-res", ARGV_FLOAT, (char *) 1, (char *) &grid_res,
       "Resolution of grid along normals [mm]."},
   {"-origin", ARGV_FLOAT, (char *) 1, (char *) &grid_origin,
       "Origin (start point) of grid along normals [mm]. Give negative values for \n\t\torigin outside the surface."},
   {"-length", ARGV_FLOAT, (char *) 1, (char *) &grid_length,
       "Lenght of grid along normals [mm]."},
   {NULL, ARGV_HELP, (char *) NULL, (char *) NULL, 
       "Interpolation options:"},
  { "-linear", ARGV_CONSTANT, (char *) 0, 
    (char *) &degrees_continuity,
    "Use linear interpolation." },
  { "-nearest_neighbour", ARGV_CONSTANT, (char *) -1, 
    (char *) &degrees_continuity,
    "Use nearest neighbour interpolation (Default)." },
  { "-cubic", ARGV_CONSTANT, (char *) 2,
    (char *) &degrees_continuity,
        "Use cubic interpolation." },
   {NULL, ARGV_HELP, (char *) NULL, (char *) NULL, 
       "Mapping function options:"},
  { "-average", ARGV_CONSTANT, (char *) 0, 
    (char *) &map_func,
    "Use average for mapping along normals." },
  { "-range", ARGV_FLOAT, (char *) 2, 
    (char *) range,
    "Count number of values in range for mapping along normals. If any value is out of range \n\t\tvalues will be counted only until this point" },
  { "-max", ARGV_CONSTANT, (char *) 2, 
    (char *) &map_func,
    "Use maximum value for mapping along normals (Default). Optionally a 2nd volume can be defined to output its value at the maximum value of the 1st volume." },
  { "-min", ARGV_CONSTANT, (char *) 3, 
    (char *) &map_func,
    "Use minimum value for mapping along normals. Optionally a 2nd volume can be defined to output its value at the minimum value of the 1st volume." },
  { "-exp", ARGV_FLOAT, (char *) 1, 
    (char *) &exp_half,
    "Use exponential average of values for mapping along normals. The argument defines the \n\t\tdistance in mm where values are decayed to 50% (recommended value is 10mm)." },
  { "-sum", ARGV_CONSTANT, (char *) 5, 
    (char *) &map_func,
    "Use sum of values for mapping along normals." },
  { NULL, ARGV_END, NULL, NULL, NULL }
};


Real
evaluate_function(Real val_array[], int n_val, int map_func, Real exp_array[], int index[])
{
        int   i, in_range;
        Real  result;
        
        index[0] = 0;
    
        switch (map_func) {
        case F_AVERAGE:
                result = 0.0;
                for (i = 0; i < n_val; i++)
                        result += val_array[i]; 
                result /= (Real) n_val;
                break;
        case F_RANGE:
                /*
                 * count only if values are in range
                 * until any value is out of range
                 */
                in_range = 0;
                result = 0.0;
                for (i = 0; i < n_val; i++) {
                        /* stop counting if values are leaving range */
                        if (in_range == 1 &&
                            (val_array[i] < range[0] ||
                             val_array[i] > range[1]))
                                break;
                        /* are values for the first time in range? */
                        if (val_array[i] >= range[0] &&
                            val_array[i] <= range[1])
                                in_range = 1;
                        /* count values in range */
                        if (in_range)
                                result += grid_res;
                }
                break;
        case F_MAX:
                result = -FLT_MAX;
                for (i = 0; i < n_val; i++) {
                        if (val_array[i] > result) {
                                result = val_array[i];
                                index[0] = i;
                        }
                } 
                break;
        case F_MIN:
                result = FLT_MAX;
                for (i = 0; i < n_val; i++) {
                        if (val_array[i] < result) {
                                result = val_array[i];
                                index[0] = i;
                        }
                }
                break;
        case F_EXP:
                /* exponential average */
                result = 0.0;
                for (i = 0; i < n_val; i++)
                        result += val_array[i]*exp_array[i];
                break;
        case F_SUM:
                result = 0.0;
                for (i = 0; i < n_val; i++)
                        result += val_array[i]; 
                break;
        }
        return(result);
}

int
main(int argc, char *argv[])
{
        char                 *input_volume_file, *input_volume_file2, *object_file;
        char                 *output_file, *output_file2;
        File_formats         format;
        Volume               volume, volume2;
        int                  i, j, index, n_values, n_objects;
        object_struct        **objects;
        polygons_struct      *polygons;
        Real                 value, value2, *values, *values2, voxel[N_DIMENSIONS];
        Real                 val_array[MAX_N_ARRAY], length_array[MAX_N_ARRAY];
        Real                 exp_sum, exp_array[MAX_N_ARRAY];
        Vector               normal;

        /* Call ParseArgv */
        if (ParseArgv(&argc, argv, argTable, 0) || ((argc != 4) && (argc != 6))) {
                fprintf(stderr, "\nUsage2: %s [options] <object.obj> <volume.mnc> <output.txt> [<volume2.mnc> <output2.txt>]\n\n", argv[0]);
                fprintf(stderr, "Map data from a minc volume to a surface.\n");
                exit(EXIT_FAILURE);
        }

        initialize_argument_processing(argc, argv);

        if (!get_string_argument(NULL, &object_file) ||
            !get_string_argument(NULL, &input_volume_file) ||
            !get_string_argument(NULL, &output_file)) {
                fprintf(stderr, "\nUsage: %s [options] <object.obj> <volume.mnc> <output.txt> [<volume2.mnc> <output2.txt>]\n\n", argv[0]);
                exit(EXIT_FAILURE);
        }
    
        /* get optional arguments for 2nd volume and output */ 
        get_string_argument(NULL, &input_volume_file2);
        get_string_argument(NULL, &output_file2);
        
        /* calculate number of values needed to subdivide grid along normals */ 
        n_values = (grid_length/grid_res) + 1;
    
        /* check minimum number of values */
        if (n_values < 2) {
                fprintf(stderr, "Resolution of grid is too low.\n");
                exit(EXIT_FAILURE);
        }
    
        /* .. and also check maximum number of values */
        if (n_values > MAX_N_ARRAY) {
                fprintf(stderr, "Resolution of grid is too high.\n");
                exit(EXIT_FAILURE);
        }

        /* if exp is given use exponential mapping function */
        if (exp_half != FLT_MAX)
                map_func = F_EXP;

        /* if range is given use range mapping function */
        if (range[0] != FLT_MAX && range[1] != FLT_MAX)
                map_func = F_RANGE;
        
        /* check range values */
        if (range[0] > range[1]) {
                fprintf(stderr, "First range value is larger than second.\n");
                exit(EXIT_FAILURE);
        }

        /* initialize values for lengths (starting with origin) */
        fprintf(stderr, "Calculate values along:\n");
        for (j = 0; j < n_values; j++) {
                length_array[j] = grid_origin +
                                  ((Real)j / (Real)(n_values-1) * grid_length);
                fprintf(stderr,"%3.2f ",length_array[j]);
        }
        fprintf(stderr, "\n");

        /* calculate exponential decay if exp function is defined */
        if (exp_half != FLT_MAX) {
                exp_sum = 0.0;
                for (j = 0; j < n_values; j++) {
                        exp_array[j] = exp(LOG05 / exp_half * length_array[j]);
                        exp_sum += exp_array[j];
                }
                /* scale sum of exponential function to 1 */
                for (j = 0; j < n_values; j++)
                        exp_array[j] /= exp_sum;
        }
    
        if (input_volume(input_volume_file, 3, File_order_dimension_names,
                         NC_UNSPECIFIED, FALSE, 0.0, 0.0,
                         TRUE, &volume, NULL) != OK)
                exit(EXIT_FAILURE);

        if (input_graphics_any_format(object_file, &format,
                                      &n_objects, &objects) != OK)
                exit(EXIT_FAILURE);

        if (output_file2 != NULL) {  
                /* check that optional 2nd volume was used either with min or max mapping function */
                if ((map_func != F_MAX) && (map_func != F_MIN)) {
                        fprintf(stderr, "For 2nd volume only min/max is allowed as mapping function.\n");
                        exit(EXIT_FAILURE);
                }
                
                if (input_volume(input_volume_file2, 3, File_order_dimension_names,
                         NC_UNSPECIFIED, FALSE, 0.0, 0.0,
                         TRUE, &volume2, NULL) != OK)
                exit(EXIT_FAILURE);
                ALLOC(values2, polygons->n_points);
        }

        if (n_objects != 1)
                printf("Warning, more than one object in file: %s\n",
                      object_file);

        polygons = get_polygons_ptr(objects[0]);

        compute_polygon_normals(polygons);

        ALLOC(values, polygons->n_points);

        for (i = 0; i < polygons->n_points; i++) {
                /* look only for inward normals */
                SCALE_VECTOR(normal, polygons->normals[i], -1.0);
                for (j = 0; j < n_values; j++) {
                        /* get point from origin in normal direction */
                        GET_grid_POINT(voxel, polygons->points[i],
                                       normal, length_array[j]);
                        evaluate_volume_in_world(volume, voxel[X], voxel[Y],
                                                 voxel[Z], degrees_continuity, 
                                                 FALSE, 0.0, &value, NULL,
                                                 NULL, NULL, NULL, NULL, NULL,
                                                 NULL, NULL, NULL);
                        val_array[j] = value;
                }
                /* evaluate function */
                value = evaluate_function(val_array, n_values,
                                          map_func, exp_array, &index);
                values[i] = value;
                
                /* get optional values for 2nd volume according to index of 1st volume */
                if (output_file2 != NULL) {
                        GET_grid_POINT(voxel, polygons->points[i],
                                       normal, length_array[index]);
                        evaluate_volume_in_world(volume2, voxel[X], voxel[Y],
                                                 voxel[Z], degrees_continuity, 
                                                 FALSE, 0.0, &value2, NULL,
                                                 NULL, NULL, NULL, NULL, NULL,
                                                 NULL, NULL, NULL);
                        values2[i] = value2;
                }
        }

        output_values_any_format(output_file, polygons->n_points,
                                 values, TYPE_REAL);

        if (output_file2 != NULL) {  
                output_values_any_format(output_file2, polygons->n_points,
                                         values2, TYPE_REAL);
                FREE(values2);
        }

        FREE(values);
        delete_object_list(n_objects, objects);
        return(EXIT_SUCCESS);
}
