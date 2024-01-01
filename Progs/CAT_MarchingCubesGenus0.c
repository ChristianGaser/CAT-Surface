/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 *
 */

#include <bicpl/marching.h>
#include <ParseArgv.h>
#include "genus0.h"
#include "CAT_Separate.h"
#include "CAT_NiftiIO.h"
#include "CAT_SurfaceIO.h"
#include "CAT_Surf.h"
#include "CAT_Smooth.h"
#include "CAT_Curvature.h"
#include "CAT_Vol.h"

#define   CHUNK_SIZE    1000000

private void extract_isosurface(
    Volume           volume,
    double           min_label,
    double           max_label,
    int              spatial_axes[],
    General_transform    *voxel_to_world_transform,
    Marching_cubes_methods method,
    BOOLEAN          binary_flag,
    double           min_threshold,
    double           max_threshold,
    double           valid_low,
    double           valid_high,
    polygons_struct  *polygons);
    
private void extract_surface(
    Marching_cubes_methods method,
    BOOLEAN          binary_flag,
    double           min_threshold,
    double           max_threshold,
    double           valid_low,
    double           valid_high,
    int              x_size,
    int              y_size,
    double           ***slices,
    double           min_label,
    double           max_label,
    double           ***label_slices,
    int              slice_index,
    BOOLEAN          right_handed,
    int              spatial_axes[],
    General_transform    *voxel_to_world_transform,
    int              ***point_ids[],
    polygons_struct  *polygons);

static char *dimension_names_3D[] = { MIzspace, MIyspace, MIxspace };
static char *dimension_names[] = { MIyspace, MIxspace };

/* argument defaults */
double min_threshold = 0.5;
double fwhm = 3.0;
double pre_fwhm = 2.0;
double scl_open = -1;
int median_correction = 1;
int use_distopen = 1;
int any_genus = 0;
int verbose = 0;
int use_thickness = 1;

/* the argument table */
static ArgvInfo argTable[] = {
  {"-thresh", ARGV_FLOAT, (char *) TRUE, (char *) &min_threshold,
    "Define the volume threshold, also known as the isovalue.\n\
     This value is crucial for initial image thresholding."},
  
  {"-pre-fwhm", ARGV_FLOAT, (char *) TRUE, (char *) &pre_fwhm,
    "Specify the Full Width Half Maximum (FWHM) for the preprocessing\n\
     smoothing filter. This helps in preserving gyri and sulci by\n\
     creating a weighted average between original and smoothed\n\
     images based on their distance to the threshold (isovalue)."},
  
  {"-fwhm", ARGV_FLOAT, (char *) TRUE, (char *) &fwhm,
    "Set FWHM for surface smoothing. This aids in correcting the mesh\n\
     in folded areas like gyri and sulci. Note: Do not use smoothing\n\
     sizes > 3 mm for reliable compensation in these areas."},
  
  {"-scl-opening", ARGV_FLOAT, (char *) TRUE, (char *) &scl_open,
    "Manually set the scaling factor for morphological opening. This\n\
     affects the isovalue used only for the opening process. Use -1\n\
     for automatic estimation."},
  
  {"-no-median", ARGV_CONSTANT, (char *) FALSE, (char *) &median_correction,
    "Disable the median filter typically used outside sulcal areas."},
  
  {"-no-distopen", ARGV_CONSTANT, (char *) FALSE, (char *) &use_distopen,
    "Turn off the additional morphological opening feature."},
  
  {"-no-genus0", ARGV_CONSTANT, (char *) TRUE, (char *) &any_genus,
    "Disable the genus0 functionality. This option skips topology\n\
     correction steps."},
  
  {"-no-thickness", ARGV_CONSTANT, (char *) TRUE, (char *) &any_genus,
    "Avoid using cortical thickness for local correction in\n\
     additional morphological opening."},
  
  {"-verbose", ARGV_CONSTANT, (char *) TRUE, (char *) &verbose,
    "Enable verbose mode for detailed output during processing."},

    {NULL, ARGV_END, NULL, NULL, NULL}
};


private void
usage(
    char *executable)
{
    char *usage_str = "\n\
Usage: CAT_MarchingCubesGenus0 input.nii output_surface_file\n\
\n\
    This method generates a mesh with an Euler number of 2 (genus 0) from the\n\
    thresholded volume. The process involves:\n\
    \n\
    1. **Preprocessing with Smoothing Filter:**\n\
       - Apply a smoothing filter to the input image to remove outliers.\n\
       - Use a weighted average of the original and smoothed images to\n\
         preserve gyri and sulci.\n\
       - Weighting is based on the distance to the isovalue.\n\
       - Weights range from 0 (for intensities at 0 or 1) to 1 (intensity\n\
         equals the isovalue), with intermediate values scaled linearly.\n\
       - Weighting effect is enhanced by squaring the value.\n\
    \n\
    2. **Morphological Opening:**\n\
       - Apply additional morphological opening, scaled by `scl_open`,\n\
         to prevent gyri fusion and minimize local artifacts.\n\
       - Opening strength is determined by analyzing the impact of\n\
         different `scl_open` values and tracking RMSE changes.\n\
    \n\
    3. **Extraction of the Largest Component:**\n\
       - Extract the largest component for further processing.\n\
    \n\
    4. **Mesh Smoothing:**\n\
       - Smooth the extracted mesh.\n\
    \n\
    5. **Mesh Correction in Folded Areas:**\n\
       - Correct the mesh in areas with folds, particularly in gyri and\n\
         sulci, to counterbalance the averaging effect from smoothing.\n\
       - Use mean curvature average as a folding measure to estimate\n\
         necessary compensation.\n\
       - Compensation degree is auto-calculated based on deviation\n\
         from the defined isovalue.\n\n";

    print_error(usage_str, executable);
}

/* convert subscripts to linear index */
static int 
sub2ind(
    int horiz, 
    int vert, 
    int depth, 
    int img_horiz, 
    int img_vert)
{  
    return(horiz+(vert+depth*img_vert)*img_horiz);
}


int   
main(
    int   argc,
    char  *argv[])
{
    char            *input_filename, *output_filename;
    Volume          volume;
    double          min_label, max_label, start_scl_open, dist;
    double          valid_low, valid_high, val, RMSE, sum_RMSE;
    double          voxelsize[N_DIMENSIONS];
    int             i, j, k, c, spatial_axes[N_DIMENSIONS];
    int             n_out, EC, sizes[MAX_DIMENSIONS];
    int             nvol, ind, count, stop_distopen, replace = 0;
    Marching_cubes_methods    method;
    object_struct       *object, **object2, *object3;
    General_transform   voxel_to_world_transform;
    polygons_struct     *polygons;
    unsigned short      *input;
    unsigned char       *input_uint8, *ref_uint8;
    float               *input_float, *input_filtered, *dist_CSF, *dist_WM, *GMT;

    /* get the arguments from the command line */
    if (ParseArgv(&argc, argv, argTable, 0)) {
            usage(argv[0]);
            fprintf(stderr, "     %s -help\n\n", argv[0]);
            exit(EXIT_FAILURE);
    }


    initialize_argument_processing(argc, argv);

    if (!get_string_argument(NULL, &input_filename) ||
        !get_string_argument(NULL, &output_filename))
    {
        usage(argv[0]);
        fprintf(stderr, "     %s -help\n\n", argv[0]);
        return(1);
    }
    
    /* marching cubes without holes */
    method = (Marching_cubes_methods) 1;

    valid_low  =  0.0;
    valid_high = -1.0;
    min_label  =  0.0;
    max_label  = -1.0;

    if (input_volume_all(input_filename, 3, dimension_names_3D,
              NC_UNSPECIFIED, FALSE, 0.0, 0.0,
              TRUE, &volume, NULL) != OK)
        return(1);

    copy_general_transform(get_voxel_to_world_transform(volume), &voxel_to_world_transform);

    for (c = 0; c < N_DIMENSIONS; c++)
        spatial_axes[c] = volume->spatial_axes[c];

    /* It is really weird, but only this combination worked */
    spatial_axes[0] = 0;
    spatial_axes[1] = 2;
    spatial_axes[2] = 1;

    /* get voxel size for optional morphological opening */
    get_volume_separations(volume, voxelsize);

    object  = create_object(POLYGONS);
    object3 = create_object(POLYGONS);
        
    get_volume_sizes(volume, sizes);
    nvol = sizes[0]*sizes[1]*sizes[2];

    input_float = (float *)malloc(nvol*sizeof(float));
    input       = (unsigned short *) malloc(nvol*sizeof(unsigned short));  
    input_uint8 = (unsigned char  *) malloc(nvol*sizeof(unsigned char));  
    ref_uint8   = (unsigned char  *) malloc(nvol*sizeof(unsigned char));  

    for (i = 0; i < sizes[0]; i++)
    for (j = 0; j < sizes[1]; j++)
    for (k = 0; k < sizes[2]; k++)
    {
        ind = sub2ind(i,j,k,sizes[0],sizes[1]);
        input_float[ind]  = get_volume_real_value(volume, i, j, k, 0, 0);
    }

    /* Preprocessing Step: Smoothing Filter
       - Purpose: To remove outliers in the input image.
       - Method: A weighted average is calculated between the original image and the 
         smoothed image to protect the structural integrity of gyri and sulci.
       - Weight Estimation: Based on the proximity to the isovalue.
         * When the intensity equals the isovalue, the weight is set to 1.
         * For intensities at 0 or 1, the weight is 0.
         * All other intensities are scaled linearly between 0 and 1.
       - Enhanced Weighting: Weights are squared to emphasize larger weightings,
         providing a more robust distinction between regions of interest and outliers.
    */
    if (pre_fwhm > 0.0) {
        input_filtered = (float *)malloc(nvol*sizeof(float));
    
        for (i = 0; i < nvol; i++)
            input_filtered[i] = input_float[i];
        double s[] = {pre_fwhm, pre_fwhm, pre_fwhm};
        smooth_float(input_filtered, sizes, voxelsize, s, 0);
        
        /* Protect values in sulci and gyri and weight areas with filtered values 
          depending on distance isovalue (threshold) */
        float weight;
        for (i = 0; i < nvol; i++) {
            /* estimate weight using distance to isovalue: weight will be 1 of intensity
               equals isovalue and is 0 for intensities of 0 or 1. Everything inbetween
               will have a linear scaling between 0..1 */
            weight = 2.0*(0.5 - fabs(input_float[i] - min_threshold));
            
            /* emphasize large weightings by using the squared value */
            weight *= weight;
            
            /* weighted average of filtered and original input */
            input_float[i] = weight*input_float[i] + (1.0 - weight)*input_filtered[i];
        }
        free(input_filtered);
    }

    /* Analyzing the Impact of scl_open Values
       - Range: scl_open values are analyzed between 1.5 and 0.5 using distopen.
       - Purpose: To track RMSE (Root Mean Square Error) changes across these values.
       - Observations:
         * Large scl_open values (closer to 1.5): These significantly and stably 
           affect the entire image, leading to smaller gyri and wider sulci.
         * Small scl_open values (closer to 0.5): These cause only localized changes,
           impacting primarily areas with artifacts (e.g., vessels, poor skull-stripping),
           resulting in smaller alterations.
       - Optimal Value Determination: The ideal scl_open value is identified at the 
         point where a significant decrease in RMSE is observed. This indicates effective 
         artifact removal while preserving the overall structure of gyri and sulci.
    */

    /* default scl_open is negative that indicates automatically search for 
       optimal scl_open value
    */
    if (scl_open < 0)
    {
        stop_distopen = 0;
        start_scl_open = 1.5;
    } else {
        stop_distopen = 1;
        start_scl_open = scl_open;
    }
    
    sum_RMSE = 0.0;
    count = 0;

        
    /* estimate cortical thickness for local correction of intensities for morphological opening */
    if (use_distopen && use_thickness && 0) {
        dist_CSF = (float *)malloc(sizeof(float)*nvol);
        dist_WM  = (float *)malloc(sizeof(float)*nvol);
        GMT      = (float *)malloc(sizeof(float)*nvol);
        
        /* check for memory faults */
        if ((dist_CSF == NULL) || (dist_WM == NULL) || (GMT == NULL)) {
            fprintf(stderr,"Memory allocation error\n");
            exit(EXIT_FAILURE);
        }
        
        /* initialize distances */
        for (i = 0; i < nvol; i++) {
            dist_CSF[i] = 0.0;
            dist_WM[i]  = 0.0;
        }

        /* prepare map outside CSF and mask to obtain distance map for CSF */
        for (i = 0; i < nvol; i++) {
            GMT[i] = (input_float[i] < 0.001) ? 1.0f : 0.0f;
            input_uint8[i]  = (input_float[i] < 1.0) ? 1 : 0;
        }
        
        for (i = 0; i < sizes[0]; i++)
        for (j = 0; j < sizes[1]; j++)
        for (k = 0; k < sizes[2]; k++)
        {
            ind = sub2ind(i,j,k,sizes[0],sizes[1]);
            set_volume_real_value(volume, i, j, k, 0, 0, GMT[ind]);
        }
        output_volume_all("test.nii", NC_FLOAT, 0, 0.0, 100000.0, volume, "test\n", NULL);

        /* obtain CSF distance map */
        vbdist(GMT, input_uint8, sizes, voxelsize, replace);
        for (i = 0; i < nvol; i++)
            dist_CSF[i] = GMT[i];

        for (i = 0; i < sizes[0]; i++)
        for (j = 0; j < sizes[1]; j++)
        for (k = 0; k < sizes[2]; k++)
        {
            ind = sub2ind(i,j,k,sizes[0],sizes[1]);
            set_volume_real_value(volume, i, j, k, 0, 0, GMT[ind]);
        }
        output_volume_all("test2.nii", NC_FLOAT, 0, 0.0, 100000.0, volume, "test\n", NULL);

        /* prepare map outside WM and mask to obtain distance map for WM */
        for (i = 0; i < nvol; i++) {
            GMT[i] = (input_float[i] > 0.999) ? 1.0f : 0.0f;
            input_uint8[i]  = (input_float[i] > 0.0) ? 1 : 0;
        }

        /* obtain WM distance map */
        vbdist(GMT, input_uint8, sizes, voxelsize, replace);
        for (i = 0; i < nvol; i++)
            dist_WM[i] = GMT[i];

        projection_based_thickness(input_float, dist_WM, dist_CSF, GMT, sizes, voxelsize); 

        for (i = 0; i < sizes[0]; i++)
        for (j = 0; j < sizes[1]; j++)
        for (k = 0; k < sizes[2]; k++)
        {
            ind = sub2ind(i,j,k,sizes[0],sizes[1]);
            set_volume_real_value(volume, i, j, k, 0, 0, GMT[ind]);
        }
        
        output_volume_all("test.nii", NC_FLOAT, 0, 0.0, 100000.0, volume, "test\n", NULL);
    }
                       
    for (scl_open = start_scl_open; scl_open > 0.4; scl_open -= 0.1)
    {
      
        /* Skip morphological opening if distopen is disabled */
        if (!use_distopen) scl_open = 1.0;
      
        /* We first apply a slightly different threshold for initial mask 
           to allow to control amount of morphological opening */
        for (i = 0; i < nvol; i++)
            input_uint8[i] = (double)input_float[i] >= scl_open*min_threshold ? 1 : 0;
    
        /* Interrupt here if distopen is disabled and use default scl_open value */
        if (!use_distopen) break;

        /* Optional morphological opening with distance criteria (distopen)
           Additionaly adapt dist parameter w.r.t. scl_open  */
        dist = 0.75/(scl_open*min_threshold);
        distopen(input_uint8, sizes, voxelsize, dist, 0.0, DT_UINT8);

        /* Apply threshold to original input image, but keep any changes from 
           the above distopen. This ensures correct position using the original
           isoval, but opens (glued) areas using a different threshold from above */
        for (i = 0; i < nvol; i++)
            input_uint8[i] = (double)input_float[i] >= min_threshold ? ((input_uint8[i] == 0) ? 0 : 1) : 0;
        
        /* Stop after one additional iteration */
        if (stop_distopen > 0) break;
        
        /* Calulate RMSE between actual and previous distopen and stop if
           changes in RMSE are getting much smaller to obtain the optimal 
           scl_open parameter */
        if (count)
        {
            RMSE = 0.0;
            for (i = 0; i < nvol; i++)
            {
                val = (double)input_uint8[i] - (double)ref_uint8[i];  
                RMSE += val*val;
            } 
            RMSE = sqrt(RMSE/(double)nvol);
            sum_RMSE += RMSE;
            
            if (verbose) fprintf(stderr,"%5.4f\t%5.4f\t%5.4f\n",scl_open,sum_RMSE/RMSE/(double)count,RMSE);
            
            /* Indicate stop if changes are getting smaller by a factor of 1.5 */
            if (sum_RMSE/RMSE/(double)count > 1.5) {
                if (!verbose) fprintf(stderr,"%5.4f\t%5.4f\t%5.4f\n",scl_open,sum_RMSE/RMSE/(double)count,RMSE);
                break;    
            }
        }
        
        /* save previous image after distopen */ 
        for (i = 0; i < nvol; i++)
            ref_uint8[i] = input_uint8[i];   

        count++;        

    }

    if (use_distopen && use_thickness && 0) {
        free(GMT);
        free(dist_CSF);
        free(dist_WM);
    }
    
    genus0parameters g0[1]; /* need an instance of genus0 parameters */

    genus0init(g0); /* initialize the instance, set default parameters */
    
    /* we need uint16 for genus0 approach */
    for (i = 0; i < nvol; i++)
        input[i] = (unsigned short)input_uint8[i];
    free(ref_uint8);

    /* set some parameters/options for the firt iteration */
    for(j = 0; j <N_DIMENSIONS; j++) g0->dims[j] = sizes[j];
    g0->connected_component = 1;
    g0->value = 1;
    g0->contour_value = 1;
    g0->any_genus = any_genus;
    g0->biggest_component = 1;
    g0->pad[0] = g0->pad[1] = g0->pad[2] = 2;
    g0->ijk2ras = NULL;
    g0->verbose = 0;
    g0->return_surface = 0;
    g0->extraijkscale[2] = 1;
        
    /* don't call loop if genus0 is not forced */
    if (any_genus) count = 10; else count = 0;
    
    EC = -1;
    
    /* repeat until EC is 2 or max. count is reached */
    while ((EC != 2) && (count < 10)) {        
        /* call genus0 for the 1st time */
        g0->input = input;
        g0->cut_loops = 0;
        g0->connectivity = 6;
        g0->alt_value = 1;
        g0->alt_contour_value = 1;
    
        /* call the function! */
        if (genus0(g0)) return(1); /* check for error */
    
        /* save results as next input */
        for (i = 0; i < nvol; i++)
            input[i] = g0->output[i];
            
        /* call genus0 a 2nd time with other parameters */
        g0->input = input;
        g0->cut_loops = 1;
        g0->connectivity = 18;
        g0->alt_value = 0;
        g0->alt_contour_value = 0;
    
        if (genus0(g0)) return(1); 
    
        /* apply median-correction after 2nd iteration */
        if ((median_correction) && (count > 1)) {
            /* use previous output for filtering */
            for (i = 0; i < nvol; i++)
                    input_uint8[i] = (unsigned char)g0->output[i];
            
            median3(input_uint8, sizes, DT_UINT8);

            /* replace with its median filtered version outside sulcal areas */
            for (i = 0; i < nvol; i++)
                g0->output[i] = (unsigned short)(input_float[i] >= min_threshold ? input_uint8[i] : 0);
        }
        
        for (i = 0; i < sizes[0]; i++)
        for (j = 0; j < sizes[1]; j++)
        for (k = 0; k < sizes[2]; k++)
        {
            ind = sub2ind(i,j,k,sizes[0],sizes[1]);
            set_volume_real_value(volume, i, j, k, 0, 0, g0->output[ind]);                    
        }

        /* extract surface to check euler number */
        extract_isosurface(volume,
                  min_label, max_label,
                  spatial_axes,
                  &voxel_to_world_transform,
                  method, FALSE,
                  0.5, 0.5,
                  valid_low, valid_high, get_polygons_ptr(object));
    
        check_polygons_neighbours_computed(get_polygons_ptr(object));
        n_out = separate_polygons(get_polygons_ptr(object), -1, &object2);
          
        triangulate_polygons(get_polygons_ptr(object2[0]), get_polygons_ptr(object3));
        polygons = get_polygons_ptr(object3);
        EC = euler_characteristic(polygons);
        count++;

        /* save results as next input */
        for (i = 0; i < nvol; i++)
            input[i] = g0->output[i];
    
    }

    free(input_uint8);
    free(input_float);
    
    if (n_out > 2) fprintf(stderr,"Extract largest of %d components.\n",n_out);
    fprintf(stderr,"Euler characteristics after %d iterations is %d.\n", count, EC);

    /* Mesh Correction in Folded Areas
       - Objective: To compensate for the averaging effect observed in gyri and sulci.
       - Method: Utilization of a folding measure, specifically the mean curvature average, 
         to estimate the necessary compensation in these areas.
       - Compensation Estimation: Automatically calculated based on the difference 
         between the actual mesh curvature and the predefined isovalue. This approach 
         ensures accurate correction in folded regions, maintaining the integrity of 
         gyri and sulci structures.
    */
    if (fwhm > 0.0) {
        smooth_heatkernel(polygons, NULL, fwhm);
        correct_mesh_folding(polygons, NULL, volume, min_threshold);
    }

    (void) output_graphics_any_format(output_filename, ASCII_FORMAT, 1, &object3, NULL);

    delete_volume(volume);
    delete_marching_cubes_table();
    delete_general_transform(&voxel_to_world_transform);

    free(input);

    return(0);
}


private void 
clear_slice(
    Volume volume,
    double **slice)
{
    int    x, y, sizes[MAX_DIMENSIONS];

    get_volume_sizes(volume, sizes);

    for (x = 0; x < sizes[0]; x++)
    for (y = 0; y < sizes[1]; y++)
    {
        slice[x][y] = 0.0;
    }
}

private void 
input_slice(
    Volume volume,
    double **slice,
    int    z)
{
    int    x, y, sizes[MAX_DIMENSIONS];

    get_volume_sizes(volume, sizes);

    for (x = 0; x < sizes[0]; x++)
    for (y = 0; y < sizes[1]; y++)
    {
        slice[x][y] = get_volume_real_value(volume, x, y, z, 0, 0);
    }
}

private double
get_slice_value(
    double ***slices,
    int    x_size,
    int    y_size,
    int    z,
    int    x,
    int    y)
{
    if (x < 0 || x >= x_size || y < 0 || y >= y_size)
        return(0.0);
    else
        return(slices[z][x][y]);
}

private void
clear_points(
    int x_size,
    int y_size,
    int max_edges,
    int ***point_ids)
{
    int x, y, edge;

    for (x = 0; x < x_size+2; x++)
    for (y = 0; y < y_size+2; y++)
    for (edge = 0; edge < max_edges; edge++)
    {
        point_ids[x][y][edge] = -1;
    }
}

private void
get_world_point(
    double slice,
    double x,
    double y,
    int    spatial_axes[],
    General_transform *voxel_to_world_transform,
    Point  *point)
{
    int    c;
    double xw, yw, zw;
    double real_voxel[N_DIMENSIONS], voxel_pos[N_DIMENSIONS];

    real_voxel[0] = slice;
    real_voxel[1] = x;
    real_voxel[2] = y;

    for (c = 0; c < N_DIMENSIONS; c++)
        voxel_pos[c] = ((spatial_axes[c] >= 0) ? real_voxel[spatial_axes[c]] : 0.0);

    general_transform_point(voxel_to_world_transform,
                              voxel_pos[X], voxel_pos[Y], voxel_pos[Z],
                              &xw, &yw, &zw);

    fill_Point(*point, xw, yw, zw);
}

private void
extract_isosurface(
    Volume  volume,
    double  min_label,
    double  max_label,
    int   spatial_axes[],
    General_transform *voxel_to_world_transform,
    Marching_cubes_methods    method,
    BOOLEAN binary_flag,
    double  min_threshold,
    double  max_threshold,
    double  valid_low,
    double  valid_high,
    polygons_struct *polygons)
{
    int   n_slices, sizes[MAX_DIMENSIONS], x_size, y_size, slice;
    int   ***point_ids[2], ***tmp_point_ids;
    int   max_edges;
    double **slices[2], **tmp_slices;
    double **label_slices[2];
    progress_struct progress;
    Surfprop spr;
    Point  point000, point100, point010, point001;
    Vector   v100, v010, v001, perp;
    BOOLEAN  right_handed;

    get_world_point(0.0, 0.0, 0.0, spatial_axes, voxel_to_world_transform, &point000);
    get_world_point(1.0, 0.0, 0.0, spatial_axes, voxel_to_world_transform, &point100);
    get_world_point(0.0, 1.0, 0.0, spatial_axes, voxel_to_world_transform, &point010);
    get_world_point(0.0, 0.0, 1.0, spatial_axes, voxel_to_world_transform, &point001);

    SUB_POINTS(v100, point100, point000);
    SUB_POINTS(v010, point010, point000);
    SUB_POINTS(v001, point001, point000);
    CROSS_VECTORS(perp, v100, v010);

    right_handed = DOT_VECTORS(perp, v001) >= 0.0;

    get_volume_sizes(volume, sizes);
    x_size = sizes[X];
    y_size = sizes[Y];
    n_slices = sizes[Z];

    ALLOC2D(slices[0], x_size, y_size);
    ALLOC2D(slices[1], x_size, y_size);

    max_edges = get_max_marching_edges(method);

    ALLOC3D(point_ids[0], x_size+2, y_size+2, max_edges);
    ALLOC3D(point_ids[1], x_size+2, y_size+2, max_edges);

    clear_slice(volume, slices[1]);

    clear_points(x_size, y_size, max_edges, point_ids[0]);
    clear_points(x_size, y_size, max_edges, point_ids[1]);

    Surfprop_a(spr) = 0.3f;
    Surfprop_d(spr) = 0.6f;
    Surfprop_s(spr) = 0.6f;
    Surfprop_se(spr)= 30.0f;
    Surfprop_t(spr) = 1.0f;
    initialize_polygons(polygons, WHITE, &spr);

    initialize_progress_report(&progress, FALSE, n_slices+1, "Extracting Surface");

    for (slice = 0; slice < n_slices; slice++)
    {
        tmp_slices = slices[0];
        slices[0] = slices[1];
        slices[1] = tmp_slices;
        if (slice < n_slices - 1)
            input_slice(volume, slices[1], slice);
        else
            clear_slice(volume, slices[1]);

        tmp_point_ids = point_ids[0];
        point_ids[0] = point_ids[1];
        point_ids[1] = tmp_point_ids;
        clear_points(x_size, y_size, max_edges, point_ids[1]);

        extract_surface(method, binary_flag, min_threshold, max_threshold,
            valid_low, valid_high,
            x_size, y_size, slices,
            min_label, max_label, label_slices, slice - 1,
            right_handed, spatial_axes, voxel_to_world_transform,
            point_ids, polygons);

        update_progress_report(&progress, slice+2);
    }

    terminate_progress_report(&progress);

    if (polygons->n_points > 0)
    {
        ALLOC(polygons->normals, polygons->n_points);
        compute_polygon_normals(polygons);
    }

    FREE2D(slices[0]);
    FREE2D(slices[1]);

    FREE3D(point_ids[0]);
    FREE3D(point_ids[1]);
}

private int
get_point_index(
    int x,
    int y,
    int slice_index,
    int x_size,
    int y_size,
    voxel_point_type *point,
    double corners[2][2][2],
    int    spatial_axes[],
    General_transform   *voxel_to_world_transform,
    BOOLEAN  binary_flag,
    double   min_threshold,
    double   max_threshold,
    int    ***point_ids[],
    polygons_struct *polygons)
{
    int    voxel[N_DIMENSIONS], edge, point_index;
    int    edge_voxel[N_DIMENSIONS];
    double v[N_DIMENSIONS];
    Point  world_point;
    Point_classes point_class;

    voxel[X] = x + point->coord[X];
    voxel[Y] = y + point->coord[Y];
    voxel[Z] = point->coord[Z];
    edge = point->edge_intersected;

    point_index = point_ids[voxel[Z]][voxel[X]+1][voxel[Y]+1][edge];
    if (point_index < 0)
    {
        edge_voxel[X] = point->coord[X];
        edge_voxel[Y] = point->coord[Y];
        edge_voxel[Z] = point->coord[Z];
        point_class = get_isosurface_point(corners, edge_voxel, edge,
                        binary_flag,
                        min_threshold, max_threshold, v);

        get_world_point(v[Z] + (Real) slice_index,
                        v[X] + (Real) x, v[Y] + (Real) y,
                        spatial_axes, voxel_to_world_transform, &world_point);

        point_index = polygons->n_points;
        ADD_ELEMENT_TO_ARRAY(polygons->points, polygons->n_points,
                        world_point, CHUNK_SIZE);

        point_ids[voxel[Z]][voxel[X]+1][voxel[Y]+1][edge] = point_index;
    }

    return(point_index);
}

private void
extract_surface(
    Marching_cubes_methods    method,
    BOOLEAN           binary_flag,
    double            min_threshold,
    double            max_threshold,
    double            valid_low,
    double            valid_high,
    int             x_size,
    int             y_size,
    double            ***slices,
    double            min_label,
    double            max_label,
    double            ***label_slices,
    int             slice_index,
    BOOLEAN           right_handed,
    int             spatial_axes[],
    General_transform     *voxel_to_world_transform,
    int             ***point_ids[],
    polygons_struct       *polygons)
{
    int         x, y, *sizes, tx, ty, tz, n_polys, ind;
    int         p, point_index, poly, size, start_points, dir;
    voxel_point_type  *points;
    double        corners[2][2][2], label;
    BOOLEAN       valid;

    for (x = -1; x < x_size; x++)
    for (y = -1; y < y_size; y++) {
        valid = TRUE;
        for (tx = 0; tx < 2; tx++)
        for (ty = 0; ty < 2; ty++)
        for (tz = 0; tz < 2; tz++)
        {
            corners[tx][ty][tz] = get_slice_value(slices, x_size, y_size,
                tz, x + tx, y + ty);
            if (valid_low <= valid_high &&
                   (corners[tx][ty][tz] < min_threshold ||
                  corners[tx][ty][tz] > max_threshold) &&
                   (corners[tx][ty][tz] < valid_low ||
                  corners[tx][ty][tz] > valid_high))
                valid = FALSE;

            if (min_label <= max_label)
            {
                label = get_slice_value(label_slices, x_size, y_size,
                          tz, x + tx, y + ty);
                if (label < min_label || label > max_label)
                    corners[tx][ty][tz] = 0.0;
            }
        }

        if (!valid)
            continue;

        n_polys = compute_isosurface_in_voxel(method, x, y, slice_index,
                          corners, binary_flag, min_threshold,
                          max_threshold, &sizes, &points);

        if (n_polys == 0)
            continue;

        if (right_handed)
        {
            start_points = 0;
            dir = 1;
        }
        else
        {
            start_points = sizes[0]-1;
            dir = -1;
        }

        for (poly = 0; poly < n_polys; poly++)
        {
            size = sizes[poly];

            start_new_polygon(polygons);

            /*--- orient polygons properly */

            for (p = 0; p < size; p++)
            {
                ind = start_points + p * dir;
                point_index = get_point_index(x, y, slice_index,
                          x_size, y_size, &points[ind], corners,
                          spatial_axes, voxel_to_world_transform,
                          binary_flag, min_threshold, max_threshold,
                          point_ids, polygons);

                ADD_ELEMENT_TO_ARRAY(polygons->indices,
                          polygons->end_indices[polygons->n_items-1],
                          point_index, CHUNK_SIZE);
            }

            if (right_handed)
                start_points += size;
            else if (poly < n_polys-1)
                start_points += sizes[poly+1];
        }
    }
}

