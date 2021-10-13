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
#include "CAT_Resample.h"

#define GET_grid_POINT(result, grid_start, normal, length) \
{ \
        (result)[X] = RPoint_x(grid_start) + length * Vector_x(normal); \
        (result)[Y] = RPoint_y(grid_start) + length * Vector_y(normal); \
        (result)[Z] = RPoint_z(grid_start) + length * Vector_z(normal); \
}

enum { F_AVERAGE, F_MEDIAN, F_RANGE, F_COUNT, F_MAXABS, F_MAX, F_MIN, F_EXP, F_SUM, F_WAVERAGE, F_MULTI };

#define LOG05       -0.69314718
#define PI2         6.28319
#define MAX_N_ARRAY 250

//#define DEBUG 1

/* argument defaults */
int  degrees_continuity = 0;        /* interpolation - default: linear */
int  grid_steps         = 7;        /* number of grid steps */
int   equivol           = 0;        /* use equi-volume approach by Bok, otherwise an equi-distance approach is used */
double grid_start       = -0.5;     /* start point (origin) of grid along normals */
double grid_end         = 0.5;      /* end point of grid along normals */
double offset_value     = 0.0;      /* offset according to thickness that is given with offset option */
int  map_func           = F_MAXABS; /* default mapping function: (absolute) maximum value */
double frange[2]        = {-FLT_MAX, FLT_MAX};
double frange_count[2]  = {-FLT_MAX, FLT_MAX};
double exp_half         = FLT_MAX;
char *thickness_file    = NULL;     /* thickness file for restricting mapping inside defined thickness */
char *offset_file       = NULL;     /* thickness file for defining offset to (central) surface */
char *sphere_src_file   = NULL;     /* source sphere for resampling */
char *sphere_trg_file   = NULL;     /* target source for resampling */
char *annot_file        = NULL;     /* annotation atlas file */

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
  {"-annot", ARGV_STRING, (char *) 1, (char *) &annot_file, 
     "Annotation atlas file for ROI partitioning."},
  {"-sphere_src", ARGV_STRING, (char *) 1, (char *) &sphere_src_file, 
     "Source sphere file for resampling of annotation file. This is usually the sphere of the input surface file."},
  {"-sphere_trg", ARGV_STRING, (char *) 1, (char *) &sphere_trg_file, 
     "Target sphere file for resampling of annotation file. This is usually the sphere of the fsaverage file."},
  {"-equivolume", ARGV_CONSTANT, (char *) TRUE, (char *) &equivol,
       "Use equi-volume approach by Bok (1929) to correct distances/layers. The correction is based on Waehnert et al. (2014).\n\t\t     This option can only be used in conjuntion with a thickness file."},
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
  { "-median", ARGV_CONSTANT, (char *) F_MEDIAN, 
    (char *) &map_func,
    "Use median for mapping along normals." },
  { "-weighted_avg", ARGV_CONSTANT, (char *) F_WAVERAGE, 
    (char *) &map_func,
    "Use weighted average with gaussian kernel for mapping along normals.\n\t\t     The kernel is so defined that values at the boundary are weighted with 50% while the center is weighted with 100%" },
  { "-range-count", ARGV_FLOAT, (char *) 2, 
    (char *) frange_count,
    "Assign a value of 1 if at least one value is in the range for assignment along the normal, 0 otherwise." },
  { "-range", ARGV_FLOAT, (char *) 2, 
    (char *) frange,
    "Assign a value of 1 if at least one value is in the range for assignment along the normal, 0 otherwise." },
  { "-maxabs", ARGV_CONSTANT, (char *) F_MAXABS, 
    (char *) &map_func,
    "Use absolute maximum value for mapping along normals (Default)." },
  { "-max", ARGV_CONSTANT, (char *) F_MAX, 
    (char *) &map_func,
    "Use maximum value for mapping along normals." },
  { "-min", ARGV_CONSTANT, (char *) F_MIN, 
    (char *) &map_func,
    "Use minimum value for mapping along normals." },
  { "-exp", ARGV_FLOAT, (char *) 1, 
    (char *) &exp_half,
    "Use exponential average of values for mapping along normals. The argument defines the \n\t\t     distance in mm where values are decayed to 50% (recommended value is 10mm)." },
  { "-sum", ARGV_CONSTANT, (char *) F_SUM, 
    (char *) &map_func,
    "Use sum of values for mapping along normals." },
  { "-multi", ARGV_CONSTANT, (char *) F_MULTI, 
    (char *) &map_func,
    "Map data for each grid step separately and save file with indicated grid value. Please note that this option is intended for high-resolution (f)MRI data only." },
  { NULL, ARGV_END, NULL, NULL, NULL }
};


/* qicksort */
void swap(double *a, double *b)
{
  double t=*a; *a=*b; *b=t;
}

void sort(double arr[], int beg, int end)
{
  if (end > beg + 1)
  {
    double piv = arr[beg];
    int l = beg + 1, r = end;
    while (l < r)
    {
      if (arr[l] <= piv)
        l++;
      else
        swap(&arr[l], &arr[--r]);
    }
    swap(&arr[--l], &arr[beg]);
    sort(arr, beg, l);
    sort(arr, r, end);
  }
}

double
evaluate_function(double val_array[], int n_val, int map_func, double kernel[], int index[])
{
        int   i, in_range;
        double  result;
        double *data_sort;
        
        index[0] = 0;
    
        switch (map_func) {
        case F_AVERAGE:
                result = 0.0;
                for (i = 0; i < n_val; i++)
                        result += val_array[i]; 
                result /= (double) n_val;
                break;
        case F_MEDIAN:
                data_sort = (double *) malloc(sizeof(double) * n_val);
                for (i = 0; i < n_val; i++) 
                       data_sort[i] = val_array[i];
         
                sort(data_sort, 0, n_val);
                result = data_sort[(int)(n_val/2)];
                free(data_sort);
                break;
        case F_WAVERAGE:
                result = 0.0;
                for (i = 0; i < n_val; i++)
                        result += val_array[i]*kernel[i];
                break;
        case F_RANGE:
                /*
                 * set to 1 if at least one value is in range, 0 otherwise
                 */
                result = 0.0;
                for (i = 0; i < n_val; i++) {
                        /* are values in range? */
                        if (val_array[i] > frange[0] && val_array[i] < frange[1])
                                result = 1.0;
                }
                break;
        case F_COUNT:
                /*
                 * count only if values are in range
                 * until any value is out of range
                 */
                in_range = 0;
                result = 0.0;
                for (i = 0; i < n_val; i++) {
                        /* stop counting if values are leaving range */
                        if (in_range == 1 &&
                            (val_array[i] < frange_count[0] ||
                             val_array[i] > frange_count[1]))
                                break;
                        /* are values for the first time in range? */
                        if (val_array[i] > frange_count[0] &&
                            val_array[i] < frange_count[1])
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
        char                 **volume_file, *object_file;
        char                 *output_values_file;
        char                 *tmp_string, ext[5];
        File_formats         format;
        Volume               volume;
        int                  i, j, k, index, n_thickness_values, grid_steps1, grid_increase;
        int                  n_objects, n_arrays, n_labels, *in_annot, n_values;
        int                  *out_annot;
        object_struct        **objects, **objects_src_sphere, **objects_trg_sphere;
        polygons_struct      *polygons, *src_sphere, *trg_sphere;
        double               value, voxel[N_DIMENSIONS], *values_atlas;
        double               *area_inner, *area_outer, *values, **values2d, *thickness;
        double               val_array[MAX_N_ARRAY], length_array[MAX_N_ARRAY];
        double               sum, x, sigma, kernel[MAX_N_ARRAY];
        double               grid_start1, grid_end1, step_size, pos;
        double               *input_values, *resampled_values, *roi_values;
        Vector               normal;
        ATABLE               *atable;
        FILE                 *fp;

        /* Call ParseArgv */
        if (ParseArgv(&argc, argv, argTable, 0)) {
                fprintf(stdout, "\nUsage: %s [options] surface_file volume_file(s) output_values_file\n\n", argv[0]);
                fprintf(stdout, "Map data from a volume to a surface.\n");
                exit(EXIT_FAILURE);
        }

        initialize_argument_processing(argc, argv);

        if (!get_string_argument(NULL, &object_file) || argc < 4) {
                fprintf(stdout, "\nUsage: %s [options] surface_file volume_file(s) output_values_file\n\n", argv[0]);
                fprintf(stdout, "Map data from a volume to a surface.\n");
                exit(EXIT_FAILURE);
        }
                
        volume_file = &argv[1];
        output_values_file = volume_file[argc-2];

        if (map_func == F_MULTI && argc > 4) {
                fprintf(stdout, "Multiple volumes cannot be used with multi-grid option.\n");
                exit(EXIT_FAILURE);
        }
        
        if (sphere_src_file != NULL  && sphere_trg_file != NULL && annot_file != NULL) {
                if ((filename_extension_matches(output_values_file,"txt") != 1) &&
                    (filename_extension_matches(output_values_file,"csv") != 1)) {
                        fprintf(stdout, "Extension of output files for ROI partitioning has to be csv or txt.\n");
                        exit(EXIT_FAILURE);
                }
        } else if (argc > 4) {
                fprintf(stdout, "Multiple volumes can only be used in conjunction with atlas annotations.\n");
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
                fprintf(stderr, "You have to define a thickness file for the equi-volume approach.\n");
                exit(EXIT_FAILURE);
        }

        if ((sphere_src_file != NULL  && (sphere_trg_file == NULL || annot_file == NULL)) ||
            (sphere_trg_file != NULL  && (sphere_src_file == NULL || annot_file == NULL)) ||
            (annot_file != NULL  && (sphere_trg_file == NULL || sphere_src_file == NULL)) ) {
                fprintf(stderr, "You have to define all three files for using annotation files:\n\
                        -sphere_src, -sphere_trg and -annot.\n");
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
        if (frange[0] != -FLT_MAX || frange[1] != FLT_MAX) {
                map_func = F_RANGE;
                /* check range values */
                if (frange[0] > frange[1]) {
                        fprintf(stderr, "First range value is larger than second.\n");
                        exit(EXIT_FAILURE);
                }
        }
                
        /* if range is given use range mapping function */
        if (frange_count[0] != -FLT_MAX || frange_count[1] != FLT_MAX) {
                map_func = F_COUNT;
                /* check range values */
                if (frange_count[0] > frange_count[1]) {
                        fprintf(stderr, "First range value is larger than second.\n");
                        exit(EXIT_FAILURE);
                }
        }
        /* initialize values for lengths (starting with origin) */
        if (thickness_file == NULL) {
                if (map_func == F_MULTI)
                        fprintf(stdout, "Save values for absolute positions [mm]:\n");
                else    fprintf(stdout, "Calculate values along absolute positions [mm]:\n");
        } else {
                /* set offset_value to 0 if thickness flag is defined too */
                if (offset_value != 0.0) {
                        offset_value = 0.0;
                        fprintf(stdout, "Offset value can only be defined together with offset flag and cannot be combined with thickness flag.\n");
                }
                if (input_values_any_format(thickness_file, &n_thickness_values, &thickness) != OK)
                        exit(EXIT_FAILURE);

                if (map_func == F_MULTI)
                        fprintf(stdout, "Save values for relative position using thickness:\n");
                else    fprintf(stdout, "Calculate values along relative position using thickness:\n");
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
        
        if (map_func == F_MULTI)
                ALLOC2D(values2d, polygons->n_points, grid_steps1);

        if (equivol) {
                /* get point area of pial (outer) surface */
                area_outer = (double *) malloc(sizeof(double) * polygons->n_points);
                get_area_of_points_central_to_pial(polygons, area_outer, thickness, 0.5);

                /* get point area of white (inner) surface */
                area_inner = (double *) malloc(sizeof(double) * polygons->n_points);
                get_area_of_points_central_to_pial(polygons, area_inner, thickness, -0.5);
        }

        /* read source and target sphere and annotation file if defined */
        if (sphere_src_file != NULL  && sphere_trg_file != NULL && annot_file != NULL) {

                if (map_func == F_MULTI) {
                        fprintf(stderr, "No multiple grid mapping possible for atlas annotation files.\n");
                        exit(EXIT_FAILURE);
                }

                if (input_graphics_any_format(sphere_trg_file, &format, &n_objects,
                    &objects_trg_sphere) != OK || n_objects != 1 ||
                    get_object_type(objects_trg_sphere[0]) != POLYGONS ) {
                        fprintf(stderr, "File %s must contain 1 polygons object.\n",
                        sphere_trg_file);
                        exit(EXIT_FAILURE);
                }


                if (filename_extension_matches(annot_file, "annot")) {
                        read_annotation_table(annot_file, &n_arrays, &in_annot, &n_labels, &atable);

                        if ((fp = fopen(output_values_file, "w")) == 0) {
                                fprintf(stderr, "output_values_file: Couldn't open file %s.\n",
                                        output_values_file);
                                return(EXIT_FAILURE);
                        }

                        trg_sphere = get_polygons_ptr(objects_trg_sphere[0]);
                        resampled_values = (double *) malloc(sizeof(double) * polygons->n_points);
                        values_atlas     = (double *) malloc(sizeof(double) * n_labels);
                        input_values     = (double *) malloc(sizeof(double) * trg_sphere->n_points);
                        roi_values       = (double *) malloc(sizeof(double) * n_labels);
                        
                        fprintf(fp,"File");
                        for (i = 0 ; i < n_labels ; i++) 
                                fprintf(fp,",%s",atable[i].name);
                        fprintf(fp,"\n");

                        
                        for (i = 0 ; i < trg_sphere->n_points ; i++) 
                                input_values[i] = (double)in_annot[i];
                                
                } else {
                        fprintf(stderr, "Only annotation files accepted.\n");
                        exit(EXIT_FAILURE);
                }

                if (input_graphics_any_format(sphere_src_file, &format, &n_objects,
                    &objects_src_sphere) != OK || n_objects != 1 ||
                    get_object_type(objects_src_sphere[0]) != POLYGONS ) {
                        fprintf(stderr, "File %s must contain 1 polygons object.\n",
                        sphere_src_file);
                        exit(EXIT_FAILURE);
                }
                src_sphere = get_polygons_ptr(objects_src_sphere[0]);
        
                objects_trg_sphere = resample_surface_to_target_sphere(trg_sphere, trg_sphere, 
                        src_sphere, input_values, resampled_values, 1);
                        
        } 

        for (k = 0; k < argc-3; k++) {
                if (input_volume_all(volume_file[k+1], 3, File_order_dimension_names,
                                     NC_UNSPECIFIED, FALSE, 0.0, 0.0,
                                     TRUE, &volume, NULL) != OK)
                        exit(EXIT_FAILURE);

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
                                                /* check that inner and outer surface area are not equal */
                                                if ((area_outer[i]-area_inner[i]) != 0) {
                                                        /* get relative position inside cortical band */
                                                        pos = length_array[j] + 0.5;
                
                                                        /* eq. 10 from Waehnert et al. 2014 */
                                                        pos = (1.0/(area_outer[i]-area_inner[i]))*
                                                                     (sqrt((pos*area_outer[i]*area_outer[i]) + ((1.0-pos)*area_inner[i]*area_inner[i]))-area_inner[i]);
                                                                     
                                                        /* subtract offset of 0.5 that was added to pos */
                                                        pos = (pos - 0.5)*thickness[i];     
                                                } else  pos = length_array[j]*thickness[i];        
                                        } else pos = length_array[j]*thickness[i];
                                }                                
                                GET_grid_POINT(voxel, polygons->points[i], normal, pos);
        
                                evaluate_volume_in_world(volume, voxel[X], voxel[Y],
                                                         voxel[Z], degrees_continuity, 
                                                         FALSE, 0.0, &value, NULL,
                                                         NULL, NULL, NULL, NULL, NULL,
                                                         NULL, NULL, NULL);
                                                         
                                if (isnan(value)) value = 0.0;
                                if (map_func == F_MULTI) values2d[i][j] = value;
                                
                                val_array[j] = value;
                        }
                                
                        if (map_func != F_MULTI)
                                /* evaluate function */
                                values[i] = evaluate_function(val_array, grid_steps1,
                                                  map_func, kernel, &index);
                }
                
                if (annot_file != NULL) {
                        for (j = 0; j < n_labels; j++) {
                                roi_values[j] = 0.0;
                                n_values = 0;
                                for (i = 0; i < polygons->n_points; i++) {
                                        if (round(resampled_values[i]) == atable[j].annotation) {
                                                roi_values[j] += values[i]; 
                                                n_values++;
                                        }
                                }
                                roi_values[j] /= (double) n_values;
                        }
                        fprintf(fp,"%s",volume_file[k+1]);
                        for (j = 0 ; j < n_labels ; j++) 
                                fprintf(fp,",%g",roi_values[j]);
                        fprintf(fp,"\n");
                }

        
                /* read source and target sphere and annotation file if defined */
                if (sphere_src_file != NULL  && sphere_trg_file != NULL && annot_file != NULL) {
                
                        output_values_any_format(output_values_file, n_labels,
                                         roi_values, TYPE_DOUBLE);
                                         
                } else if (map_func == F_MULTI) {
                        ALLOC(tmp_string, string_length(output_values_file)+10);
                        
                        /* remove potential extension for output name */
                        strcpy(tmp_string,output_values_file);
                        if (filename_extension_matches(output_values_file,"txt")) {
                                output_values_file[string_length(output_values_file)-4] = '\0';
                                strcpy(ext,".txt");
                        } else  strcpy(ext,"");
        
                        /* prepare numbered output name and write values */
                        for (j = 0; j < grid_steps1; j++) {
                                if (equivol)
                                        (void) sprintf(tmp_string,"%s_ev%3.2f%s",output_values_file,length_array[j],ext);
                                else    (void) sprintf(tmp_string,"%s_ed%3.2f%s",output_values_file,length_array[j],ext);
                                for (i = 0; i < polygons->n_points; i++) values[i] = values2d[i][j];
        
                                output_values_any_format(tmp_string, polygons->n_points,
                                         values, TYPE_DOUBLE);
                        }
                        free(values2d);
                        free(tmp_string);
                                
                } else  {
                        if (filename_extension_matches(output_values_file, "annot")) {
                                n_arrays = polygons->n_points;
                                out_annot  = (int *) malloc(n_arrays * sizeof(int));
                                
                                n_labels = 0;
                                for (i = 0 ; i < n_arrays ; i++) {
                                        out_annot[i] = (int)round(values[i]);
                                        n_labels = MAX(n_labels, out_annot[i]);
                                }

                                atable = (ATABLE*) malloc(n_labels * sizeof(ATABLE));
                                for (i = 0 ; i < n_labels ; i++) {
                                        sprintf(atable[i].name,"Label %d",i+1);
                                        atable[i].annotation = i+1;
                                        /* g/b colors should be set to 0 because annotation label is built on 
                                           atable[i].annotation = atable[i].r+atable[i].g*256+atable[i].b*65536 */
                                        atable[i].r = i+1;
                                        atable[i].g = 0;
                                        atable[i].b = 0;
                                }

                                write_annotation_table(output_values_file, n_arrays, out_annot, n_labels, atable);
                                
                                free(atable);
                                free(out_annot);
                        } else  output_values_any_format(output_values_file, polygons->n_points,
                                                 values, TYPE_DOUBLE);
                }
        }
        
        if (sphere_src_file != NULL  && sphere_trg_file != NULL && annot_file != NULL) {
                free(input_values);
                free(resampled_values);
                free(values_atlas);
                free(roi_values);
                fclose(fp);
        }
        
        free(values);
        if (equivol) {
                free(area_inner);
                free(area_outer);
        }
        delete_object_list(n_objects, objects);

        return(EXIT_SUCCESS);
}

