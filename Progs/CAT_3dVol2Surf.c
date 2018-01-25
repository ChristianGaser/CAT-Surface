/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 *
 */

#include <bicpl.h>
#include <float.h>
#include <ParseArgv.h>

#include "CAT_SurfaceIO.h"
#include "CAT_NiftiIO.h"
#include "CAT_Surf.h"

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
#define PI2         6.28319
#define MAX_N_ARRAY 250

//#define DEBUG 1

/* argument defaults */
int  degrees_continuity = 0;        /* interpolation - default: linear */
int  grid_steps         = 7;        /* number of grid steps */
int   equivol           = 0;        /* sse equi-volume model by Bok */
double grid_start       = -0.5;     /* start point (origin) of grid along normals */
double grid_end         = 0.5;      /* end point of grid along normals */
double offset_value     = 0.0;      /* offset according to thickness that is given with offset option */
int  map_func           = F_MAXABS; /* default mapping function: (absolute) maximum value */
double frange[2]        = {FLT_MAX, FLT_MAX};
double exp_half         = FLT_MAX;
char *thickness_file    = NULL;     /* thickness file for restricting mapping inside defined thickness */
char *offset_file       = NULL;     /* thickness file for defining offset to (central) surface */

/* the argument table */
ArgvInfo argTable[] = {
  {"-start", ARGV_FLOAT, (char *) 1, (char *) &grid_start,
       "Start point (origin) of grid along normals [mm]. Give negative values for a start point\n\t\t     outside of the surface (outwards).\n\t\t     If thickness is used to define mapping the grid resolution is considered\n\t\t     as normalized value according to the cortical thickness."},
  {"-steps", ARGV_INT, (char *) 1, (char *) &grid_steps,
       "Number of grid steps."},
  {"-end", ARGV_FLOAT, (char *) 1, (char *) &grid_end,
       "End point of the grid along the surface normals (pointing inwards) in mm."},
  {"-thickness", ARGV_STRING, (char *) 1, (char *) &thickness_file, 
     "Additional thickness file for mapping inside defined normalized cortical thickness of GM.\n\t\t     If this option is used then -start and -end will be handled as normalized (relative) values:\n\t\t     e.g. start=-0.5, steps=7 and end=0.5 for a central surface will map all values inside the GM-band (-0.5:1/6:0.5)\n\t\t     that is defined using the normalized cortical thickness."},
  {"-offset", ARGV_STRING, (char *) 1, (char *) &offset_file, 
     "Additional thickness file defining an offset according to the given surface.\n\t\t     If this option is used then also use the option -offset_value to define the offset (default 0)."},
  {"-offset_value", ARGV_FLOAT, (char *) 1, (char *) &offset_value,
       "Offset to the surface according to a thickness file. A value of 0.5 means that the \n\t\t     WM surface will be used if a central surface is used as input (adding half of the thickness).\n\t\t     A negative value of -0.5 can be used to define the pial surface."},
  {"-equivolume", ARGV_CONSTANT, (char *) TRUE, (char *) &equivol,
       "Use equi-volume model by Bok (1929) to correct distances/layers. The correction is based on Waehnert et al. (2014).\n\t\t     Using this option the mappings for each defined step are saved separately and\n\t\t     no special mapping functions will be used. \n\t\t     This option can only be used if a thickness file is defined."},
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
    "Use absolute maximum value for mapping along normals (Default). Optionally a 2nd volume can be defined\n\t\t     to output its value at the maximum value of the 1st volume." },
  { "-max", ARGV_CONSTANT, (char *) F_MAX, 
    (char *) &map_func,
    "Use maximum value for mapping along normals. Optionally a 2nd volume can be defined\n\t\t     to output its value at the maximum value of the 1st volume." },
  { "-min", ARGV_CONSTANT, (char *) F_MIN, 
    (char *) &map_func,
    "Use minimum value for mapping along normals. Optionally a 2nd volume can be defined to\n\t\t     output its value at the minimum value of the 1st volume." },
  { "-exp", ARGV_FLOAT, (char *) F_EXP, 
    (char *) &exp_half,
    "Use exponential average of values for mapping along normals. The argument defines the \n\t\t     distance in mm where values are decayed to 50% (recommended value is 10mm)." },
  { "-sum", ARGV_CONSTANT, (char *) F_SUM, 
    (char *) &map_func,
    "Use sum of values for mapping along normals." },
  { NULL, ARGV_END, NULL, NULL, NULL }
};


double
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
                result /= (double) n_val;
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
        char                 *volume_file, *object_file;
        char                 *output_values_file;
        char                 *tmp_string, ext[5];
        File_formats         format;
        Volume               volume, volume2;
        int                  i, j, index, n_thickness_values, n_objects, grid_steps1, grid_increase;
        object_struct        **objects;
        polygons_struct      *polygons;
        double               value, value2, voxel[N_DIMENSIONS];
        double               *area_inner, *area_outer, *values, **values2, *thickness;
        double               val_array[MAX_N_ARRAY], length_array[MAX_N_ARRAY];
        double               sum, x, sigma, kernel[MAX_N_ARRAY];
        double               grid_start1, grid_end1, step_size, pos;
        Vector               normal;

        /* Call ParseArgv */
        if (ParseArgv(&argc, argv, argTable, 0)) {
                fprintf(stdout, "\nUsage: %s [options] surface_file volume_file output_values_file\n\n", argv[0]);
                fprintf(stdout, "Map data from a volume to a surface.\n");
                exit(EXIT_FAILURE);
        }

        initialize_argument_processing(argc, argv);

        if (!get_string_argument(NULL, &object_file) || !get_string_argument(NULL, &volume_file) 
            || !get_string_argument(NULL, &output_values_file)) {
                fprintf(stdout, "\nUsage: %s [options] surface_file volume_file output_values_file\n\n", argv[0]);
                fprintf(stdout, "Map data from a volume to a surface.\n");
                exit(EXIT_FAILURE);
        }
                
        /* we need larger values because gaussian kernel exceeds defined grid-values */
        if (map_func == F_WAVERAGE) {
                grid_increase = round(2.0*(grid_steps - 1.0)/3.0);
                step_size = (grid_end - grid_start)/((double)grid_steps - 1.0);
                
                /* force even numbers */
                if (grid_increase % 2) grid_increase++;
                
                /* extend grid values to cover almost a whole gaussian curve */
                grid_steps1 = grid_steps + grid_increase;
                grid_start1 = grid_start - (double)grid_increase * step_size / 2.0;
                grid_end1   = grid_end   + (double)grid_increase * step_size / 2.0;
#ifdef DEBUG
                printf("grid_increase: %d\tgrid_start: %3.3f\tgrid_end: %3.3f\n",grid_increase,grid_start1,grid_end1);
#endif
        } else {
                grid_steps1 = grid_steps;
                grid_start1 = grid_start;
                grid_end1   = grid_end;
        }

        /* check that the option for equivolume is used together with the thickness option */
        if ((equivol)  && (thickness_file == NULL)) {
                fprintf(stderr, "You have to define a thickness file for the equivolume model.\n");
                exit(EXIT_FAILURE);
        }

        /* check maximum number of values */
        if (grid_steps1 > MAX_N_ARRAY) {
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
                        fprintf(stdout, "Offset value can only be defined together with offset flag and cannot be combined with thickness flag.\n");
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
        
        for (j = 0; j < grid_steps1; j++) {
                length_array[j] = grid_start1;

                /* only use grid calculation if more than 1 value is given */
                if (grid_steps1 > 1) length_array[j] += ((double)j / (double)(grid_steps1-1) * (grid_end1 - grid_start1));
                
                if ((length_array[j] >= grid_start) && ((length_array[j] - grid_end) < 1e-10))
                        fprintf(stdout,"%3.2f ",length_array[j]);
        }
        fprintf(stdout, "\n");

        /* calculate exponential decay if exp function is defined */
        if (exp_half != FLT_MAX) {
                sum = 0.0;
                for (j = 0; j < grid_steps1; j++) {
                        kernel[j] = exp(LOG05 / exp_half * length_array[j]);
                        sum += kernel[j];
                }
                /* scale sum of exponential function to 1 */
                for (j = 0; j < grid_steps1; j++) {
                        kernel[j] /= sum;
#ifdef DEBUG
                        printf("%g ",kernel[j]);
#endif
                }
#ifdef DEBUG
                printf("\n");
#endif
        }
    
        /* calculate gaussian kernel if weighted average function is defined */
        if (map_func == F_WAVERAGE) {
                sum = 0.0;
                sigma = -1.0/(2.0 * LOG05); 
                for (i = 0; i < grid_steps1; i++) {
                        x = ((2.0*(double)i) / ((double)grid_steps - 1.0)) - ((double)grid_steps1 - 1.0)/((double)grid_steps - 1.0);
                        kernel[i] = (1.0/sqrt(PI2*sigma))*exp(-(x*x)/(2.0*sigma));
                        sum += kernel[i];
                }
                
                /* scale sum of gaussian kernel to 1 */
                for (j = 0; j < grid_steps1; j++) {
                        kernel[j] /= sum;
#ifdef DEBUG
                        printf("%g ",kernel[j]);
#endif
                }
#ifdef DEBUG
                printf("\n");
#endif
        }

        if (input_volume_all(volume_file, 3, File_order_dimension_names,
                             NC_UNSPECIFIED, FALSE, 0.0, 0.0,
                             TRUE, &volume, NULL) != OK)
                exit(EXIT_FAILURE);

        if (input_graphics_any_format(object_file, &format,
                                      &n_objects, &objects) != OK)
                exit(EXIT_FAILURE);

        if (n_objects != 1)
                printf("Warning, more than one object in file: %s\n",
                      object_file);

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
        
        if (equivol) {
                ALLOC2D(values2, polygons->n_points, grid_steps1);

                /* get point area of pial (outer) surface */
                area_outer = (double *) malloc(sizeof(double) * polygons->n_points);
                get_area_of_points_central_to_pial(polygons, area_outer, thickness, 0.5);

                /* get point area of white (inner) surface */
                area_inner = (double *) malloc(sizeof(double) * polygons->n_points);
                get_area_of_points_central_to_pial(polygons, area_inner, thickness, -0.5);
        }

        for (i = 0; i < polygons->n_points; i++) {
        
                /* look only for inward normals */
                SCALE_VECTOR(normal, polygons->normals[i], -1.0);
                
                for (j = 0; j < grid_steps1; j++) {

                        /* get point from origin in normal direction */
                        if (thickness_file == NULL) {
                                if (offset_file != NULL) {
                                        pos = length_array[j] + (offset_value*thickness[i]); 
                                } else {
                                        pos = length_array[j]; 
                                }
                        } else { /* relate grid position to thickness values */
                                
                                if (equivol) {
                                        /* get relative position inside cortical band */
                                        pos = length_array[j] + 0.5;
                                        /* eq. 10 from Waehnert et al. 2014 */
                                        pos = (1.0/(area_outer[i]-area_inner[i]))*
                                                     (sqrt((pos*area_outer[i]*area_outer[i]) + ((1.0-pos)*area_inner[i]*area_inner[i]))-area_inner[i]);
                                        pos = (pos - 0.5)*thickness[i]; /* subtract offset of 0.5 that was added to pos */            
                                } else pos = length_array[j]*thickness[i];
                        }
                        
                        GET_grid_POINT(voxel, polygons->points[i], normal, pos);

                        evaluate_volume_in_world(volume, voxel[X], voxel[Y],
                                                 voxel[Z], degrees_continuity, 
                                                 FALSE, 0.0, &value, NULL,
                                                 NULL, NULL, NULL, NULL, NULL,
                                                 NULL, NULL, NULL);
                                                 
                        if (isnan(value)) value = 0.0;
                        if (equivol) values2[i][j] = value;
                        
                        val_array[j] = value;
                }

                if (equivol==0)
                        /* evaluate function */
                        values[i] = evaluate_function(val_array, grid_steps1,
                                          map_func, kernel, &index);
        }

        if (equivol) {
                ALLOC(tmp_string, string_length(output_values_file)+3);
                
                /* remove potential extension for output name */
                strcpy(tmp_string,output_values_file);
                if (filename_extension_matches(output_values_file,"txt")) {
                        output_values_file[string_length(output_values_file)-4] = '\0';
                        strcpy(ext,".txt");
                } else  strcpy(ext,"");

                /* prepare numbered output name and write values */
                for (j = 0; j < grid_steps1; j++) {
                        (void) sprintf(tmp_string,"%s_%d%s",output_values_file,j+1,ext);
                        for (i = 0; i < polygons->n_points; i++) values[i] = values2[i][j];

                        output_values_any_format(tmp_string, polygons->n_points,
                                 values, TYPE_DOUBLE);
                }
                free(values2);
                free(area_inner);
                free(area_outer);
                        
        } else  output_values_any_format(output_values_file, polygons->n_points,
                                 values, TYPE_DOUBLE);

        free(values);
        
        delete_object_list(n_objects, objects);
        return(EXIT_SUCCESS);
}

