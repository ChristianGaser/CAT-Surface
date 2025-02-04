/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
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
double local_smoothing = 10.0;
double min_threshold = 0.5;
double post_fwhm = 2.0;
double pre_fwhm = 2.0;
double scl_open = 0.9;
int median_correction = 2;
int use_distopen = 1;
int any_genus = 0;
int verbose = 0;
int iter = 10;

/* the argument table */
static ArgvInfo argTable[] = {
  {"-thresh", ARGV_FLOAT, (char *) TRUE, (char *) &min_threshold,
    "Define the volume threshold, also known as the isovalue.\n\
     This value is crucial for initial image thresholding."},
  
  {"-pre-fwhm", ARGV_FLOAT, (char *) TRUE, (char *) &pre_fwhm,
    "Specify the Full Width Half Maximum (FWHM) for the preprocessing\n\
     smoothing filter. This helps in preserving gyri and sulci by\n\
     creating a weighted average between original and smoothed\n\
     images based on the gradient of the input image. Areas with \n\
     topology artefacts are often characterized by large gradients,\n\
     thus smoothing in these areas tries to prevent these artefacts."},

  {"-post-fwhm", ARGV_FLOAT, (char *) TRUE, (char *) &post_fwhm,
    "Set FWHM for surface smoothing. This aids in correcting the mesh\n\
     in folded areas like gyri and sulci. Note: Do not use smoothing\n\
     sizes > 3 mm for reliable compensation in these areas."},
  
  {"-median-filter", ARGV_INT, (char *) TRUE, (char *) &median_correction,
    "Specify the number of iterations to apply a median filter to areas\n\
     where the gradient of the thresholded image indicates larger clusters.\n\
     These clusters may point to potential topology artifacts and regions\n\
     with high local variations. This process helps to smooth these areas, \n\
     improving the quality of the surface reconstruction in subsequent steps."},
  
  {"-scl-opening", ARGV_FLOAT, (char *) TRUE, (char *) &scl_open,
    "Manually set the scaling factor for morphological opening. This\n\
     affects the isovalue for the opening process.\n\
     Use -1 for automatic estimation."},
  
  {"-iter", ARGV_INT, (char *) TRUE, (char *) &iter,
    "Number of iterations."},
  
  {"-local-smoothing", ARGV_FLOAT, (char *) TRUE, (char *) &local_smoothing,
    "Apply local surface smoothing to resulting surface in areas where the distance\n\
     between the surface and a shifted surface is below the expected distance,\n\
     which often happens due to self intersections of the surface."},
  
  {"-no-distopen", ARGV_CONSTANT, (char *) FALSE, (char *) &use_distopen,
    "Turn off the additional morphological opening feature."},
  
  {"-no-genus0", ARGV_CONSTANT, (char *) TRUE, (char *) &any_genus,
    "Disable the genus0 functionality. This option skips topology\n\
     correction steps."},
  
  {"-verbose", ARGV_CONSTANT, (char *) TRUE, (char *) &verbose,
    "Enable verbose mode for detailed output during processing."},

    {NULL, ARGV_END, NULL, NULL, NULL}
};


private void
usage(
    char *executable)
{
    char *usage_str = "\n\
Usage: CAT_VolMarchingCubes input.nii output_surface_file\n\
\n\
    This method generates a mesh with an Euler number of 2 (genus 0) from the\n\
    thresholded volume. The process involves:\n\
    \n\
    1. **Preprocessing with Smoothing Filter:**\n\
       - Apply a smoothing filter to the input image to remove outliers.\n\
       - Use a weighted average of the original and smoothed images to\n\
         preserve gyri and sulci.\n\
       - Weighting is based on the gradient of the input image.\n\
       - Weights range from 0 (areas with low gradient) to 1 (areas\n\
         with large gradient), with intermediate values scaled linearly.\n\
       - Weighting effect is enhanced by squaring the value.\n\
    \n\
    2. **Preprocessing with Median Filter:**\n\
       - Apply an iterative median filter to areas where the gradient of \n\
         the thresholded image indicates larger clusters.\n\
       - Use a weighted average of the original and median filterd images.\n\
       - Weighting is estimated using gradient of the input image and.\n\
         morphological operations to find larger clusters\n\
    \n\
    3. **Morphological Opening:**\n\
       - Apply additional morphological opening, scaled by `scl_open`,\n\
         to prevent gyri fusion and minimize local artifacts.\n\
       - Opening strength is determined by analyzing the impact of\n\
         different `scl_open` values and tracking RMSE changes.\n\
    \n\
    4. **Extraction of the Largest Component:**\n\
       - Extract the largest component for further processing.\n\
    \n\
    5. **Mesh Smoothing:**\n\
       - Smooth the extracted mesh.\n\
    \n\
    6. **Mesh Correction in Folded Areas:**\n\
       - Correct the mesh in areas with folds, particularly in gyri and\n\
         sulci, to counterbalance the averaging effect from smoothing.\n\
       - Use mean curvature average as a folding measure to estimate\n\
         necessary compensation.\n\
       - Compensation degree is auto-calculated based on deviation\n\
         from the defined isovalue.\n\
    7. **Mesh Correction in Areas with Self Intersections:**\n\
       - Apply local surface smoothing to resulting surface in areas where\n\
         the distance between the surface and a shifted surface is below \n\
         the expected distance, which often happens due to self intersections\n\
         of the surface.\n\n";

    fprintf(stderr,"%s\n %s\n",usage_str, executable);
}

#define IDX(x, y, z, nx, ny) ((z) * (nx) * (ny) + (y) * (nx) + (x))


void correct_topology(float *volume, float thresh, int dims[3], int conn_arr[2], int n_loops) {
    float *vol_euler;
    unsigned short *vol_bin;
    int loop, i, x, y, z, nx, ny, nz, iter, val_corr;
    int V, E, F, conn, chi_local, n_errors, cube[8];
    int edges18[18][2] = {
        {0,1}, {1,3}, {3,2}, {2,0}, // Bottom face edges
        {4,5}, {5,7}, {7,6}, {6,4}, // Top face edges
        {0,4}, {1,5}, {2,6}, {3,7}, // Vertical edges
        {0,5}, {1,4}, {2,7}, {3,6}, // Diagonal edges across faces
        {0,7}, {3,4}  // Extended diagonal connections
    };
    int edges26[26][2] = {
        {0,1}, {1,3}, {3,2}, {2,0}, // Bottom face edges
        {4,5}, {5,7}, {7,6}, {6,4}, // Top face edges
        {0,4}, {1,5}, {2,6}, {3,7}, // Vertical edges
        {0,5}, {1,4}, {2,7}, {3,6}, // Diagonal edges across faces
        {0,7}, {3,4}, {2,5}, {1,6}, // Extended diagonal connections
        {0,6}, {1,7}, {2,4}, {3,5}  // Corner-to-corner vertex edges
    };
    // Count faces (6 faces in a cube)
    int faces[6][4] = {
        {0,1,3,2}, {4,5,7,6}, // XY faces
        {0,2,6,4}, {1,3,7,5}, // YZ faces
        {0,1,5,4}, {2,3,7,6}  // XZ faces
    };
    
    vol_euler = (float *)malloc(sizeof(float)*dims[0]*dims[1]*dims[2]);
    vol_bin = (unsigned short *)malloc(sizeof(unsigned short)*dims[0]*dims[1]*dims[2]);

    nx = dims[0], ny = dims[1], nz = dims[2];

    // Initialize Euler map
    for (i = 0; i < nx * ny * nz; i++)
        vol_euler[i] = (volume[i] >= thresh) ? 1.0 : 0.0;

    for (iter = 0; iter < 50; iter++) {
        n_errors = 0;
        for (loop = 0; loop < n_loops; loop++) {
    
            // Threshold volume
            for (i = 0; i < nx * ny * nz; i++)
                vol_bin[i] = (volume[i] > thresh) ? 1.0 - (float)loop : 0.0 + (float)loop;
                
            conn = conn_arr[((loop + iter) % 2)];
            val_corr = (loop == 0) ? 0.0 : 1.0;

            for (z = 0; z < nz - 1; z++) {
                for (y = 0; y < ny - 1; y++) {
                    for (x = 0; x < nx - 1; x++) {
                            
                        // Extract 2x2x2 voxel cube
                        cube[0] = vol_bin[IDX(x, y, z, nx, ny)];
                        cube[1] = vol_bin[IDX(x+1, y, z, nx, ny)];
                        cube[2] = vol_bin[IDX(x, y+1, z, nx, ny)];
                        cube[3] = vol_bin[IDX(x+1, y+1, z, nx, ny)];
                        cube[4] = vol_bin[IDX(x, y, z+1, nx, ny)];
                        cube[5] = vol_bin[IDX(x+1, y, z+1, nx, ny)];
                        cube[6] = vol_bin[IDX(x, y+1, z+1, nx, ny)];
                        cube[7] = vol_bin[IDX(x+1, y+1, z+1, nx, ny)];
                        
                        // Count vertices, edges, faces
                        V = 0; E = 0; F = 0;
        
                        // Count vertices (fully occupied corners)
                        for (i = 0; i < 8; i++) V += cube[i];
                        
                        if (conn == 18) {
                            for (i = 0; i < 18; i++) 
                                if (cube[edges18[i][0]] && cube[edges18[i][1]]) E++;
                        } else if (conn == 26) {
                            for (i = 0; i < 26; i++) 
                                if (cube[edges26[i][0]] && cube[edges26[i][1]]) E++;
                        }
        
                        for (i = 0; i < 6; i++)
                            if (cube[faces[i][0]] && cube[faces[i][1]] && cube[faces[i][2]] && cube[faces[i][3]]) F++;
        
                        // Compute local Euler number
                        chi_local = V - E + F;
        
                        if (vol_bin[IDX(x, y, z, nx, ny)] > 0) {
                            if (((conn == 18) &&  (chi_local ==  2)) || 
                                ((conn == 26) &&  (chi_local == -6))) {
                                vol_euler[IDX(x,   y,   z,   nx, ny)] = val_corr;
                                vol_euler[IDX(x+1, y,   z,   nx, ny)] = val_corr;
                                vol_euler[IDX(x,   y+1, z,   nx, ny)] = val_corr;
                                vol_euler[IDX(x+1, y+1, z,   nx, ny)] = val_corr;
                                vol_euler[IDX(x,   y,   z+1, nx, ny)] = val_corr;
                                vol_euler[IDX(x+1, y,   z+1, nx, ny)] = val_corr;
                                vol_euler[IDX(x,   y+1, z+1, nx, ny)] = val_corr;
                                vol_euler[IDX(x+1, y+1, z+1, nx, ny)] = val_corr;
                                n_errors++;
                            }
                        }
                    }
                }
            }
        }
        //printf("%d\n", n_errors);
        for (i = 0; i < dims[0]*dims[1]*dims[2]; i++)
            volume[i] = vol_euler[i];
            
        if (n_errors == 0) break;
    }
    
    free(vol_euler);
    free(vol_bin);

}

int   
main(
    int   argc,
    char  *argv[])
{
    char                *input_filename, *output_filename;
    char                out_diff[1024], out_bin[1024];;
    double              min_label, max_label, start_scl_open, dist;
    double              valid_low, valid_high, val, RMSE, sum_RMSE;
    double              *values, *extents, voxelsize[N_DIMENSIONS];
    double              x, y, z, min_value;
    double              max_value, mean_grad, max_grad;
    float               *input_float, *vol_float;
    float               *grad, weight;
    int                 i, j, k, c, n_values, conn_arr[2], n_loops;
    int                 n_out, EC, sizes[MAX_DIMENSIONS];
    int                 nvol, ind, count, stop_distopen, replace = 0;
    unsigned short      *input_uint16;
    unsigned char       *input_uint8, *vol_uint8;
    Marching_cubes_methods    method;
    object_struct       *object, **object2, *object3, **object4;
    polygons_struct     *polygons, *smooth_polygons;
    nifti_image         *nii_ptr, *out_ptr;
    mat44               nii_mat;

    initialize_argument_processing(argc, argv);

    /* get the arguments from the command line */
    if (ParseArgv(&argc, argv, argTable, 0)) {
        usage(argv[0]);
        fprintf(stderr, "     %s -help\n\n", argv[0]);
        exit(EXIT_FAILURE);
    }

    if (!get_string_argument(NULL, &input_filename) ||
        !get_string_argument(NULL, &output_filename))
    {
        usage(argv[0]);
        fprintf(stderr, "     %s -help\n\n", argv[0]);
        return(1);
    }
    
    if(argc > 3)
        (void) sprintf(out_diff, "%s", argv[3]); 

    if(argc > 4)
        (void) sprintf(out_bin, "%s", argv[4]); 

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

    grad         = (float *)malloc(nvol*sizeof(float));
    vol_float    = (float *)malloc(nvol*sizeof(float));
    input_uint16 = (unsigned short *) malloc(nvol*sizeof(unsigned short));  
    input_uint8  = (unsigned char  *) malloc(nvol*sizeof(unsigned char));  
    vol_uint8    = (unsigned char  *) malloc(nvol*sizeof(unsigned char));

    if (!vol_float || !input_uint16 || !input_uint8 || !vol_uint8) {
        fprintf(stderr,"Memory allocation error\n");
        exit(EXIT_FAILURE);
    }

    /* We have to scale image to a maximum of 1.0 if isoval > 1.0 */
    if (min_threshold > 1.0) {
        max_value = get_max_float(input_float, nvol, 0);
        for (i = 0; i < nvol; i++) input_float[i] /= (float)max_value;
        min_threshold = min_threshold/max_value;
    }
    
    /* Preprocessing Step: Smoothing Filter
       - Purpose: To remove outliers in the input image.
       - Method: A weighted average is calculated between the original image and the 
         smoothed image to protect the structural integrity of gyri and sulci.
       - Weight Estimation: Based on the gradient of the input image.
         * For values below mean gradient, the weight is 0.
       - Enhanced Weighting: Weights are squared to emphasize larger weightings,
         providing a more robust distinction between regions of interest and outliers.
    */
    if (pre_fwhm != 0.0) {
    
        for (i = 0; i < nvol; i++)
            vol_float[i] = input_float[i];
        double s[] = {pre_fwhm, pre_fwhm, pre_fwhm};
        smooth3(vol_float, sizes, voxelsize, s, 0, DT_FLOAT32);
        
        /* estimate magnitude of gradient for weighting the smoothing */
        gradient3D_magnitude(input_float, grad, sizes);

        /* calculate mean of gradient (where gradient > 0) */
        mean_grad = get_mean_float(grad, nvol, 1);
        max_grad  = get_max_float(grad, nvol, 1);
        
        /* Protect values in sulci and gyri and weight areas with filtered values 
          depending on gradient of input image */
        for (i = 0; i < nvol; i++) {
            /* estimate weight using gradient of input image and limit range to 0..1 */
            weight = (grad[i] - mean_grad)/(max_grad - mean_grad);
            weight = (weight < 0.0) ? 0.0f : weight;
            weight = (weight > 1.0) ? 1.0f : weight;
            
            /* emphasize large weightings by using the squared value */
            weight *= weight;
            
            /* weighted average of filtered and original input */
            input_float[i] = (1.0 - weight)*input_float[i] + weight*vol_float[i];
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

    /* apply cluster function the 1st time and keep largest cluster after thresholding */
    keep_largest_cluster(input_float, min_threshold, sizes, DT_FLOAT32, 0, 1, 18);
    fill_holes(input_float, min_threshold, sizes, DT_FLOAT32);

    conn_arr[0] = 18; conn_arr[1] = 26;
    n_loops = 2;
    correct_topology(input_float, min_threshold, sizes, conn_arr, n_loops);

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
            
            if (verbose) fprintf(stderr,"%5.2f\t%5.4f\t%5.4f\n",scl_open,sum_RMSE/RMSE/(double)count,RMSE);
            
            /* Indicate stop if changes are getting smaller by a factor of 1.5 */
            if (sum_RMSE/RMSE/(double)count > 1.5) {
                if (!verbose) fprintf(stderr,"Final threshold for distopen: %5.2f\n",scl_open);
                break;    
            }
        }
        
        /* save previous image after distopen */ 
        for (i = 0; i < nvol; i++)
            vol_uint8[i] = input_uint8[i];   

        count++;        
    }
    
    for (i = 0; i < nvol; i++)
        grad[i] = 0.0;

    genus0parameters g0[1]; /* need an instance of genus0 parameters */

    genus0init(g0); /* initialize the instance, set default parameters */
    
    /* we need uint16 for genus0 approach */
    for (i = 0; i < nvol; i++)
        input_uint16[i] = (unsigned short)input_uint8[i];

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
    
    count = 0;
    EC = -1;
    
    /* repeat until EC is 2 or max. count is reached */
    while ((EC != 2) && (count < iter)) {        
        /* call genus0 for the 1st time */
        g0->input = input_uint16;
        g0->cut_loops = 0;
        g0->connectivity = 6;
        g0->alt_value = 1;
        g0->alt_contour_value = 1;

        /* call the function */
        if (~any_genus) {
            if (genus0(g0)) return(1); /* check for error */
        
            /* save results as next input */
            for (i = 0; i < nvol; i++)
                input_uint16[i] = g0->output[i];
        }
        
        for (i = 0; i < nvol; i++)
            grad[i] += (float)(count + 1)*((float)g0->output[i] - (float)g0->input[i]);   

        /* call genus0 a 2nd time with other parameters */
        g0->input = input_uint16;
        g0->cut_loops = 1;
        g0->connectivity = 18;
        g0->alt_value = 0;
        g0->alt_contour_value = 0;
    
        if (any_genus) {
            for (i = 0; i < nvol; i++)
                g0->output[i] = input_uint16[i];
            count = 10;
        } else {
            if (genus0(g0)) return(1); 
        }

        for (i = 0; i < nvol; i++)
            grad[i] += (float)(count + 1)*((float)g0->output[i] - (float)g0->input[i]);   
 
        /* apply median-correction and only consider the dilated areas */
        if (median_correction) {
            /* find areas that were corrected for topology artefacts and dilate them */
            for (i = 0; i < nvol; i++)
                vol_uint8[i] = (unsigned char)(input_uint16[i] != g0->output[i]);
                
            morph_dilate(vol_uint8, sizes, 4, 0.5, DT_UINT8);

            /* use previous output for filtering */
            for (i = 0; i < nvol; i++)
                input_uint16[i] = g0->output[i];
            
            /* apply iterative median filter */
            for (i = 0; i < median_correction; i++) median3(input_uint16, NULL, sizes, DT_UINT16);

            /* replace genus0 output with its median filtered version in (dilated)
               areas with topology artefacts */
            for (i = 0; i < nvol; i++)
                g0->output[i] = (unsigned short)(((vol_uint8[i] > 0)) ? input_uint16[i] : g0->output[i]);
        }
        
        /* apply cluster function a 2nd time and keep largest cluster after thresholding */
        keep_largest_cluster(g0->output, min_threshold, sizes, DT_UINT16, 0, 1, 18);
        printf(".");
        fill_holes(g0->output, min_threshold, sizes, DT_UINT16);

        for (i = 0; i < nvol; i++)
            vol_float[i] = (float)g0->output[i];
            
        /* extract surface to check euler number */
        extract_isosurface(vol_float, sizes,
                  min_label, max_label,
                  nii_mat,
                  method, FALSE,
                  min_threshold, min_threshold,
                  valid_low, valid_high, polygons, verbose);
        compute_polygon_normals(polygons);
        check_polygons_neighbours_computed(polygons);
        n_out = separate_polygons(polygons, -1, &object2);          
        triangulate_polygons(get_polygons_ptr(object2[0]), get_polygons_ptr(object3));
        polygons = get_polygons_ptr(object3);

        EC = euler_characteristic(polygons);
        count++;
        if (verbose) fprintf(stderr,"Euler characteristics after %d iterations is %d.\n", count, EC);

        /* save results as next input */
        for (i = 0; i < nvol; i++)
            input_uint16[i] = g0->output[i];
    
    }

    free(input_uint8);

    if (verbose && (n_out > 1)) fprintf(stderr,"Extract largest of %d components.\n",n_out);

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

    /* Mesh Correction in Areas with Self Intersections
       - Objective: To correct areas with self intersections
       - Method: Apply local surface smoothing to resulting surface in areas 
         where the distance between the surface and a shifted surface is below 
         the expected distance, which often happens due to self intersections of 
         the surface.
    */
    if (local_smoothing > 0.0) {
        values  = (double *) malloc(sizeof(double) * polygons->n_points);
        extents = (double *) malloc(sizeof(double) * polygons->n_points);

        smooth_polygons = (polygons_struct *) malloc(sizeof(polygons_struct));
        copy_polygons(polygons, smooth_polygons);

        for (i = 0; i < polygons->n_points; i++) {
            extents[i] = 0.1;
            values[i]  = 3.0;
        }

        object4 = central_to_new_pial(polygons, values, extents, NULL, NULL, 0);
        compute_exact_hausdorff(polygons, get_polygons_ptr(object4[0]), values, 0);
        smooth_heatkernel(polygons, values, 5.0);
        
        min_value = 0.25;
        max_value = 0.3;
        for (i = 0; i < polygons->n_points; i++) {
            /* scale values between 0.25..0.3 and force min=0 and max=1 */
            values[i] = (values[i] - min_value)/(max_value - min_value);
            values[i] = MIN(1.0, values[i]);
            values[i] = MAX(0.0, values[i]);
            values[i] = values[i]*values[i];
        }
        
        /* smooth values to obtain a smooth border */
        smooth_heatkernel(polygons, values, 5.0);
        
        /* obtain smoothed surface */
        smooth_heatkernel(smooth_polygons, NULL, local_smoothing);

        /* use smoothed or original surface w.r.t. weighting */ 
        for (i = 0; i < polygons->n_points; i++)
            Point_x(polygons->points[i]) = values[i]*Point_x(polygons->points[i]) + (1.0 - values[i])*Point_x(smooth_polygons->points[i]);

        compute_polygon_normals(smooth_polygons);
        
        free(smooth_polygons);
        free(values);
        free(extents);
    }

    if(argc > 3) {
        out_ptr = nifti_copy_nim_info(nii_ptr);
        if (!write_nifti_float(out_diff, grad, DT_FLOAT32, 1.0, sizes, voxelsize, out_ptr)) 
            exit(EXIT_FAILURE);
    }

    if(argc > 4) {
        out_ptr = nifti_copy_nim_info(nii_ptr);
        if (!write_nifti_float(out_bin, vol_float, DT_FLOAT32, 1.0, sizes, voxelsize, out_ptr)) 
            exit(EXIT_FAILURE);
    }

    compute_polygon_normals(get_polygons_ptr(object3));
    output_graphics_any_format(output_filename, ASCII_FORMAT, 1, &object3, NULL);

    free(grad);
    free(vol_uint8);
    free(input_uint16);
    free(input_float);
    free(vol_float);
    delete_marching_cubes_table();

    return(0);
}
