/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 *
 */

#include <ParseArgv.h>
#include "CAT_Separate.h"
#include "CAT_NiftiLib.h"
#include "CAT_SurfaceIO.h"
#include "CAT_Surf.h"
#include "CAT_Smooth.h"
#include "CAT_Curvature.h"
#include "CAT_Vol.h"
#include "CAT_MarchingCubes.h"
#include "genus0.h"

/* argument defaults */
double min_threshold = 0.5;
double post_fwhm = 2.0;
double pre_fwhm = -1.0;
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
     images based on their distance to the threshold (isovalue).\n\
     A negative value will force masked smoothing, which may\n\
     preserves gyri and sulci even better."},

  {"-post-fwhm", ARGV_FLOAT, (char *) TRUE, (char *) &post_fwhm,
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

    fprintf(stderr,"%s\n %s\n",usage_str, executable);
}


int   
main(
    int   argc,
    char  *argv[])
{
    char                *input_filename, *output_filename;
    double              min_label, max_label, start_scl_open, dist;
    double              valid_low, valid_high, val, RMSE, sum_RMSE;
    double              voxelsize[N_DIMENSIONS];
    int                 i, j, k, c;
    int                 n_out, EC, sizes[MAX_DIMENSIONS];
    int                 nvol, ind, count, stop_distopen, replace = 0;
    Marching_cubes_methods    method;
    object_struct       *object, **object2, *object3;
    polygons_struct     *polygons;
    unsigned short      *input;
    unsigned char       *input_uint8, *vol_uint8;
    float               *input_float, *vol_float, *dist_CSF, *dist_WM, *GMT;
    nifti_image         *nii_ptr;
    mat44               nii_mat;

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

    nii_ptr = read_nifti_float(input_filename, &input_float, 0);
    if (!nii_ptr) {
            fprintf(stderr,"Error reading %s.\n", input_filename);
            return(EXIT_FAILURE);
    }

    nii_mat = nii_ptr->sto_xyz; /* 4x4 transformation matrix */
    
    sizes[0] = nii_ptr->nx;
    sizes[1] = nii_ptr->ny;
    sizes[2] = nii_ptr->nz;
    voxelsize[0] = nii_ptr->dx;
    voxelsize[1] = nii_ptr->dy;
    voxelsize[2] = nii_ptr->dz;

    object  = create_object(POLYGONS);
    object3 = create_object(POLYGONS);
    polygons = get_polygons_ptr(object);
        
    nvol = sizes[0]*sizes[1]*sizes[2];

    vol_float   = (float *)malloc(nvol*sizeof(float));
    input       = (unsigned short *) malloc(nvol*sizeof(unsigned short));  
    input_uint8 = (unsigned char  *) malloc(nvol*sizeof(unsigned char));  
    vol_uint8   = (unsigned char  *) malloc(nvol*sizeof(unsigned char));

    if (!vol_float || !input || !input_uint8 || !vol_uint8) {
        fprintf(stderr,"Memory allocation error\n");
        exit(EXIT_FAILURE);
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
    if (pre_fwhm != 0.0) {
    
        for (i = 0; i < nvol; i++)
            vol_float[i] = input_float[i];
        double s[] = {fabs(pre_fwhm), fabs(pre_fwhm), fabs(pre_fwhm)};
        smooth_float(vol_float, sizes, voxelsize, s, (pre_fwhm < 0.0));
        
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
            input_float[i] = weight*input_float[i] + (1.0 - weight)*vol_float[i];
        }
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
        if (!dist_CSF || !dist_WM || !GMT) {
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
        
        /* obtain CSF distance map */
        vbdist(GMT, input_uint8, sizes, voxelsize, replace);
        for (i = 0; i < nvol; i++)
            dist_CSF[i] = GMT[i];

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

    }

    for (scl_open = start_scl_open; scl_open > 0.4; scl_open -= 0.1)
    {
        /* Skip morphological opening if distopen is disabled */
        if (!use_distopen) scl_open = 1.0;
      
        /* We first apply a slightly different threshold for initial mask 
           to allow to control amount of morphological opening */
        for (i = 0; i < nvol; i++)
            input_uint8[i] = (double)input_float[i] >= (scl_open*min_threshold) ? 1 : 0;
    
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
                val = (double)input_uint8[i] - (double)vol_uint8[i];  
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
            vol_uint8[i] = input_uint8[i];   

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
    free(vol_uint8);

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
        
        /* keep largest cluster after thresholding */
        keep_largest_cluster(g0->output, min_threshold, sizes, DT_UINT16, 0, 1);

        for (i = 0; i < nvol; i++)
            vol_float[i] = (float)g0->output[i];

        /* extract surface to check euler number */
        extract_isosurface(vol_float, sizes,
                  min_label, max_label,
                  nii_mat,
                  method, FALSE,
                  0.5, 0.5,
                  valid_low, valid_high, get_polygons_ptr(object));
    
        compute_polygon_normals(polygons);

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
    
    if (n_out > 1) fprintf(stderr,"Extract largest of %d components.\n",n_out);
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
    if (post_fwhm > 0.0) {
        smooth_heatkernel(polygons, NULL, post_fwhm);
        correct_mesh_folding(polygons, NULL, input_float, nii_ptr, min_threshold);
    }

    output_graphics_any_format(output_filename, ASCII_FORMAT, 1, &object3, NULL);

    free(input);
    free(input_float);
    free(vol_float);
    delete_marching_cubes_table();

    return(0);
}
