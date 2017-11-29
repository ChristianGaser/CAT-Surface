/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id: CAT_3dVol2Surf.c 378 2017-08-08 08:23:22Z gaser $
 *
 */

#include <bicpl.h>
#include <float.h>
#include <ParseArgv.h>

#include "CAT_SurfaceIO.h"
#include "CAT_NiftiIO.h"

#define GET_grid_POINT(result, grid_start, normal, length) \
{ \
        (result)[X] = RPoint_x(grid_start) + length * Vector_x(normal); \
        (result)[Y] = RPoint_y(grid_start) + length * Vector_y(normal); \
        (result)[Z] = RPoint_z(grid_start) + length * Vector_z(normal); \
}

#define F_AVERAGE  0
#define F_RANGE    1
#define F_MAXABS   2
#define F_MAX      3
#define F_MIN      4
#define F_EXP      5
#define F_SUM      6
#define F_WAVERAGE 7

#define LOG05       -0.69314718
#define MAX_N_ARRAY 250

/* argument defaults */
int  degrees_continuity = 0;          /* interpolation - default: linear */
int  grid_steps         = 10;         /* number of grid steps */
double grid_start       = 0.0;        /* start point (origin) of grid along normals */
double grid_end         = 5.0;        /* end point of grid along normals */
double offset_value     = 0.0;        /* offset according to thickness that is given with offset option */
int  map_func           = F_MAXABS;   /* default mapping function: (absolute) maximum value */
double frange[2]        = {FLT_MAX, FLT_MAX};
double exp_half         = FLT_MAX;
char *thickness_file    = NULL;       /* thickness file for restricting mapping inside defined thickness */
char *offset_file       = NULL;       /* thickness file for defining offset to (central) surface */
char *pre_outname       = "vol2surf"; /* prepend to output file name */
char *post_outname      = "";         /* append to output file name */

/* the argument table */
ArgvInfo argTable[] = {
  {"-prepend", ARGV_STRING, (char *) 1, (char *) &pre_outname, 
     "Prepend to output filename."},
  {"-append", ARGV_STRING, (char *) 1, (char *) &post_outname, 
     "Append to output filename."},
  {"-start", ARGV_FLOAT, (char *) 1, (char *) &grid_start,
       "Start point (origin) of grid along normals [mm]. Give negative values for a start point\n\t\t     outside of the surface (outwards).\n\t\t     If thickness is used to define mapping the grid resolution is considered\n\t\t     as normalized value according to the cortical thickness."},
  {"-steps", ARGV_INT, (char *) 1, (char *) &grid_steps,
       "Number of grid steps."},
  {"-end", ARGV_FLOAT, (char *) 1, (char *) &grid_end,
       "End point of the grid along the surface normals (pointing inwards) in mm."},
  {"-thickness", ARGV_STRING, (char *) 1, (char *) &thickness_file, 
     "Additional thickness file for mapping inside defined normalized cortical thickness of GM.\n\t\t     If this option is used then -start and -end will be handled as normalized (relative) values:\n\t\t     e.g. start=-0.5, steps=11 and end=0.5 for a central surface will map all values inside the GM-band (-0.5:0.1:0.5)\n\t\t     that is defined using the normalized cortical thickness."},
  {"-offset", ARGV_STRING, (char *) 1, (char *) &offset_file, 
     "Additional thickness file defining an offset according to the given surface.\n\t\t     If this option is used then also use the option -offset_value to define the offset (default 0)."},
  {"-offset_value", ARGV_FLOAT, (char *) 1, (char *) &offset_value,
       "Offset to the surface according to a thickness file. A value of 0.5 means that the \n\t\t     WM surface will be used if a central surface is used as input (adding half of the thickness).\n\t\t     A negative value of -0.5 can be used to define the pial surface."},
  {NULL, ARGV_HELP, (char *) NULL, (char *) NULL, 
       "Interpolation options:"},
  { "-linear", ARGV_CONSTANT, (char *) 0, 
    (char *) &degrees_continuity,
    "Use linear interpolation (Default)." },
  { "-nearest_neighbour", ARGV_CONSTANT, (char *) -1, 
    (char *) &degrees_continuity,
    "Use nearest neighbour interpolation." },
  { "-cubic", ARGV_CONSTANT, (char *) 2,
    (char *) &degrees_continuity,
        "Use cubic interpolation." },
   {NULL, ARGV_HELP, (char *) NULL, (char *) NULL, 
       "Mapping function options:"},
  { "-avg", ARGV_CONSTANT, (char *) F_AVERAGE, 
    (char *) &map_func,
    "Use average for mapping along normals." },
  { "-weighted_avg", ARGV_CONSTANT, (char *) F_WAVERAGE, 
    (char *) &map_func,
    "Use weighted average with gaussian kernel for mapping along normals.\n\t\t     The kernel is so defined that values at the boundary are weighted with 50% while the center is weighted with 100%" },
  { "-range", ARGV_FLOAT, (char *) F_RANGE, 
    (char *) frange,
    "Count number of values in range for mapping along normals. If any value is out of range \n\t\t     values will be counted only until this point" },
  { "-maxabs", ARGV_CONSTANT, (char *) F_MAXABS, 
    (char *) &map_func,
    "Use absolute maximum value for mapping along normals (Default)." },
  { "-max", ARGV_CONSTANT, (char *) F_MAX, 
    (char *) &map_func,
    "Use maximum value for mapping along normals." },
  { "-min", ARGV_CONSTANT, (char *) F_MIN, 
    (char *) &map_func,
    "Use minimum value for mapping along normals." },
  { "-exp", ARGV_FLOAT, (char *) F_EXP, 
    (char *) &exp_half,
    "Use exponential average of values for mapping along normals. The argument defines the \n\t\t     distance in mm where values are decayed to 50% (recommended value is 10mm)." },
  { "-sum", ARGV_CONSTANT, (char *) F_SUM, 
    (char *) &map_func,
    "Use sum of values for mapping along normals." },
  { NULL, ARGV_END, NULL, NULL, NULL }
};


Real
evaluate_function(double val_array[], int n_val, int map_func, double kernel[], int index[])
{
        int   i, in_range;
        double  result;
        
        index[0] = 0;
    
        switch (map_func) {
        case F_AVERAGE:
                result = 0.0;
                for (i = 0; i < n_val; i++)
                        result += val_array[i]; 
                result /= (Real) n_val;
                break;
        case F_WAVERAGE:
                result = 0.0;
                for (i = 0; i < n_val; i++)
                        result += val_array[i]*kernel[i];
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
                            (val_array[i] < frange[0] ||
                             val_array[i] > frange[1]))
                                break;
                        /* are values for the first time in range? */
                        if (val_array[i] >= frange[0] &&
                            val_array[i] <= frange[1])
                                in_range = 1;
                        /* count values in range */
                        if (in_range)
                                result++;
                }
                break;
        case F_MAXABS:
                result = 0;
                for (i = 0; i < n_val; i++) {
                        if (fabs(val_array[i]) > fabs(result)) {
                                result = val_array[i];
                                index[0] = i;
                        }
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
                        result += val_array[i]*kernel[i];
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
        char                 **volume_files, *surface_file, *output_values_file, buffer[1024];
        File_formats         format;
        Volume               volume;
        int                  i, j, k, index, n_thickness_values, n_objects;
        object_struct        **objects;
        polygons_struct      *polygons;
        double               value, *values, *thickness, voxel[N_DIMENSIONS];
        double               val_array[MAX_N_ARRAY], length_array[MAX_N_ARRAY];
        double               sum, kernel[MAX_N_ARRAY], fwhm, x;
        Vector               normal;

        /* Call ParseArgv */
        if (ParseArgv(&argc, argv, argTable, 0)) {
                fprintf(stderr, "\nUsage: %s [options] surface_file volume_file(s)\n\n", argv[0]);
                fprintf(stderr, "Map data from a volume to a surface.\n");
                exit(EXIT_FAILURE);
        }

        initialize_argument_processing(argc, argv);

        if (!get_string_argument(NULL, &surface_file) || argc < 3) {
                fprintf(stderr, "\nUsage: %s [options] surface_file volume_file(s)\n\n", argv[0]);
                fprintf(stderr, "Map data from a volume to a surface.\n");
                exit(EXIT_FAILURE);
        }

        volume_files = &argv[1];
        
        /* check maximum number of values */
        if (grid_steps > MAX_N_ARRAY) {
                fprintf(stderr, "Resolution of grid is too high.\n");
                exit(EXIT_FAILURE);
        }

        /* if exp is given use exponential mapping function */
        if (exp_half != FLT_MAX)
                map_func = F_EXP;

        /* if range is given use range mapping function */
        if (frange[0] != FLT_MAX && frange[1] != FLT_MAX)
                map_func = F_RANGE;
        
        /* check range values */
        if (frange[0] > frange[1]) {
                fprintf(stderr, "First range value is larger than second.\n");
                exit(EXIT_FAILURE);
        }

        /* initialize values for lengths (starting with origin) */
        if (thickness_file == NULL)
                fprintf(stdout, "Calculate values along absolute positions [mm]:\n");
        else {
                /* set offset_value to 0 if thickness flag is defined too */
                if (offset_value != 0.0) {
                        offset_value = 0.0;
                        fprintf(stdout, "Offset value can only be defined together with offset flag.\n");
                }
                if (input_values_any_format(thickness_file, &n_thickness_values, &thickness) != OK)
                        exit(EXIT_FAILURE);

                fprintf(stdout, "Calculate values along relative position using thickness:\n");
        }

        /* use offset file to get other surfaces */
        if (offset_file != NULL) {
                /* check that not both thickness and offset are defined */
                if (thickness_file != NULL) {
                        fprintf(stderr, "Please only define either thickness or offset.\n");
                        exit(EXIT_FAILURE);
                }
                if (input_values_any_format(offset_file, &n_thickness_values, &thickness) != OK)
                        exit(EXIT_FAILURE);

                fprintf(stdout, "Use relative offset of %g of thickness to the surface:\n", offset_value);
        }

        for (j = 0; j < grid_steps; j++) {
                length_array[j] = grid_start;

                /* only use grid calculation if more than 1 value is given */
                if (grid_steps > 1) length_array[j] += ((Real)j / (Real)(grid_steps-1) * (grid_end - grid_start));
                
                fprintf(stdout,"%3.2f ",length_array[j]);
        }
        fprintf(stdout, "\n");

        /* calculate exponential decay if exp function is defined */
        if (exp_half != FLT_MAX) {
                sum = 0.0;
                for (j = 0; j < grid_steps; j++) {
                        kernel[j] = exp(LOG05 / exp_half * length_array[j]);
                        sum += kernel[j];
                }
                /* scale sum of exponential function to 1 */
                for (j = 0; j < grid_steps; j++)
                        kernel[j] /= sum;
        }
    
        /* calculate gaussian kernel if weighted average function is defined */
        if (map_func == F_WAVERAGE) {
                sum = 0.0;
                fwhm = sqrt((double)grid_steps/2.5); /* fwhm is approximated that extreme values at the borders are weighted with 50% and 
                                center with 100% */
                for (i = 0; i < grid_steps; i++) {
                        x = ((double)i+1) - ((double)grid_steps + 1.0)/2.0;
                        kernel[i] = (1.0/sqrt(6.28*fwhm))*exp(-(x*x)/(2.0*fwhm));
                        sum += kernel[i];
                }
                
                /* scale sum of gaussian kernel to 1 */
                for (i = 0; i < grid_steps; i++)
                        kernel[i] /= sum;
        }

        if (input_graphics_any_format(surface_file, &format,
                                      &n_objects, &objects) != OK)
                exit(EXIT_FAILURE);

        if (n_objects != 1)
                printf("Warning, more than one object in file: %s\n",
                      surface_file);

        polygons = get_polygons_ptr(objects[0]);

        /* check whether # of thickness values fits to surface */
        if ((thickness_file != NULL) || (offset_file != NULL)) {
                if (polygons->n_points != n_thickness_values) {
                        fprintf(stderr,"Number of surface points differs from number of thickness values.\n");
                        exit(EXIT_FAILURE);
                }
        }

        compute_polygon_normals(polygons);

        ALLOC(values, polygons->n_points);

        for (k = 0; k < argc-2; k++) {

                if (input_volume_all(volume_files[k+1], 3, File_order_dimension_names,
                                     NC_UNSPECIFIED, FALSE, 0.0, 0.0,
                                     TRUE, &volume, NULL) != OK)
                        exit(EXIT_FAILURE);
        
                for (i = 0; i < polygons->n_points; i++) {
                        /* look only for inward normals */
                        SCALE_VECTOR(normal, polygons->normals[i], -1.0);
                        

                        for (j = 0; j < grid_steps; j++) {
                        
                                /* get point from origin in normal direction */
                                if (thickness_file == NULL) {
                                        if (offset_file != NULL) {
                                                GET_grid_POINT(voxel, polygons->points[i],
                                                        normal, length_array[j] + (offset_value*thickness[i])); 
                                        } else {
                                                GET_grid_POINT(voxel, polygons->points[i],
                                                        normal, length_array[j]); 
                                        }
                                } else { /* relate grid position to thickness values */
                                        GET_grid_POINT(voxel, polygons->points[i],
                                               normal, length_array[j]*thickness[i]);
                                }
                                
                                evaluate_volume_in_world(volume, voxel[X], voxel[Y],
                                                         voxel[Z], degrees_continuity, 
                                                         FALSE, 0.0, &value, NULL,
                                                         NULL, NULL, NULL, NULL, NULL,
                                                         NULL, NULL, NULL);
                                val_array[j] = value;
                        }
                        
                        /* evaluate function */
                        value = evaluate_function(val_array, grid_steps,
                                                  map_func, kernel, &index);
                        values[i] = value;
                        
                }
        
                /* remove file extension and add output name */
                strcpy(buffer,volume_files[k+1]);
                buffer[strlen(buffer)-4] = 0;
                output_values_file = create_string(pre_outname);
                concat_to_string(&output_values_file, "_");
                concat_to_string(&output_values_file, buffer);
                concat_to_string(&output_values_file, post_outname);
                
                fprintf(stdout, "Saving %s\n", output_values_file);
                output_values_any_format(output_values_file, polygons->n_points,
                                         values, TYPE_DOUBLE);

        }
                
        FREE(values);
        delete_object_list(n_objects, objects);
        return(EXIT_SUCCESS);
}
