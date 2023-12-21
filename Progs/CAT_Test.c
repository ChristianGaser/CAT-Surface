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
#include "CAT_NiftiLib.h"
#include "CAT_SurfaceIO.h"
#include "CAT_Surf.h"
#include "CAT_Smooth.h"
#include "CAT_Curvature.h"
#include "CAT_Vol.h"
#include "MarchingCubes.h"

/* argument defaults */
double min_threshold = 0.5;
double fwhm   = 3.0;
double scl_open = -1;
int median_correction  = 1;
int distopen  = 1;
int any_genus = 0;
int verbose = 0;
int use_thickness = 1;

/* the argument table */
static ArgvInfo argTable[] = {
  {"-thresh", ARGV_FLOAT, (char *) TRUE, (char *) &min_threshold,
      "Volume threshold (i.e. isovalue)."},
  {"-fwhm", ARGV_FLOAT, (char *) TRUE, (char *) &fwhm,
      "Create a slightly smoothed surface that is corrected in folded areas to compensate for the averaging effect in gyri and sulci.\n\
       Do not use smoothing sizes > 3 mm because that compensation only works reliably for smaller smoothing sizes."},
  {"-scl-opening", ARGV_FLOAT, (char *) TRUE, (char *) &scl_open,
      "Manually define scaling factor for morphological opening that is used to change isovalue for opening only (default -1 for automatic estimation)."},
  {"-no-median", ARGV_CONSTANT, (char *) FALSE, (char *) &median_correction,
      "Disable median filter that is used outside sulcal areas."},
  {"-no-distopen", ARGV_CONSTANT, (char *) FALSE, (char *) &distopen,
      "Disable additional morphological opening."},
  {"-no-genus0", ARGV_CONSTANT, (char *) TRUE, (char *) &any_genus,
      "Allows to switch off genus0 functionality."},
  {"-no-thickness", ARGV_CONSTANT, (char *) TRUE, (char *) &any_genus,
      "Do not use cortical thickness for local correction of additional morphological opening."},
  {"-verbose", ARGV_CONSTANT, (char *) TRUE, (char *) &verbose,
      "Be verbose."},
    {NULL, ARGV_END, NULL, NULL, NULL}
};


private void
usage(
    char *executable)
{
    char *usage_str = "\n\
Usage: CAT_MarchingCubesGenus0 input.nii output_surface_file\n\
\n\
    This method creates a mesh with Euler number 2 (genus 0) of the\n\
    thresholded volume. Here are the steps involved:\n\
    1. Apply additional morphological opening (using scl_open as scaling factor)\n\
       to prevent glued gyri and minimizes local artifacts caused by small vessels\n\
       or mis-segmentations. The strength of the opening is automatically\n\
       estimated by analyzing the impact of different scl_open values\n\
       and tracking the changes in RMSE between these values.\n\
    2. Extract largest component.\n\
    3. Smooth the resulting mesh.\n\
    4. Correct mesh in folded areas to compensate for the averaging effect\n\
       in gyri and sulci. We use a folding measure (i.e. mean curvature\n\
       averaged) to estimate the compensation. The amount of compensation\n\
       is automatically estimated using the difference to the defined\n\
       isovalue.\n\n";

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
    double          min_label, max_label, start_scl_open, dist, xw, yw, zw;
    double          valid_low, valid_high, val, RMSE, sum_RMSE;
    double          voxelsize[N_DIMENSIONS];
    int             i, j, k, c, spatial_axes[N_DIMENSIONS];
    int             n_out, EC, sizes[MAX_DIMENSIONS];
    int             nvol, ind, count, stop_distopen;
    object_struct       *object, **object2, *object3;
    General_transform   voxel_to_world_transform;
    polygons_struct     *polygons;
    unsigned short      *input;
    unsigned char       *input_uint8, *ref_uint8;
    unsigned int        *input_uint16;
    float               *input_float, *dist_CSF, *dist_WM, *GMT;
    Point               point;
    MCB                 *mcb ;
    nifti_image         *nii_ptr;
    double              real_voxel[N_DIMENSIONS], voxel_pos[N_DIMENSIONS];

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
    
    valid_low  =  0.0;
    valid_high = -1.0;
    min_label  =  0.0;
    max_label  = -1.0;

    nii_ptr = read_nifti_float(input_filename, &input_float, 0);
    if(nii_ptr == NULL) {
            fprintf(stderr,"Error reading %s.\n", input_filename);
            return(EXIT_FAILURE);
    }
    
    sizes[0] = nii_ptr->nx;
    sizes[1] = nii_ptr->ny;
    sizes[2] = nii_ptr->nz;
    voxelsize[0] = nii_ptr->dx;
    voxelsize[1] = nii_ptr->dy;
    voxelsize[2] = nii_ptr->dz;

    /* It is really weird, but only this combination worked */
    spatial_axes[0] = 2;
    spatial_axes[1] = 0;
    spatial_axes[2] = 1;

    object  = create_object(POLYGONS);
    object3 = create_object(POLYGONS);
    polygons = get_polygons_ptr(object);
        
    nvol = sizes[0]*sizes[1]*sizes[2];

    input       = (unsigned short *) malloc(nvol*sizeof(unsigned short));  
    input_uint8 = (unsigned char  *) malloc(nvol*sizeof(unsigned char));  
    ref_uint8   = (unsigned char  *) malloc(nvol*sizeof(unsigned char));  

    mcb = MarchingCubes(-1, -1, -1);
    set_resolution(mcb, sizes[0], sizes[1], sizes[2]);
    init_all(mcb);

    /* We analyze the impact of different scl_open values ranging from 1.5 to 0.5 
       using distopen. By doing so, we can track the changes in RMSE (Root Mean 
       Square Error) between these values. When using large scl_open values, the 
       changes in RMSE are relatively significant and stable because they affect 
       the whole image resulting in smaller gyri and wider sulci. In contrast, 
       using smaller scl_open values results in only local changes that are much 
       smaller since only regions with artifacts such as vessels or poor skull-
       stripping are affected. The optimal scl_open value is determined by identifying 
       the point where the RMSE decreases significantly, indicating successful 
       artifact removal while maintaining global gyri and sulci characteristics. 
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

    if (!distopen) fprintf(stderr,"%5s\t%5s\t%5s\n","Dist","avgRMSE","RMSE");
    
    /* estimate cortical thickness for local correction of intensities for morphological opening */
    if (distopen && use_thickness && 0) {
        dist_CSF = (float *)malloc(sizeof(float)*nvol);
        dist_WM  = (float *)malloc(sizeof(float)*nvol);
        GMT      = (float *)malloc(sizeof(float)*nvol);
        input_uint16 = (unsigned int *)malloc(sizeof(unsigned int)*nvol);
        
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
            input_uint16[i]  = (input_float[i] < 1.0) ? 1 : 0;
        }
        
        /* obtain CSF distance map */
        vbdist(GMT, input_uint16, sizes, voxelsize);
        for (i = 0; i < nvol; i++)
            dist_CSF[i] = GMT[i];

        /* prepare map outside WM and mask to obtain distance map for WM */
        for (i = 0; i < nvol; i++) {
            GMT[i] = (input_float[i] > 0.999) ? 1.0f : 0.0f;
            input_uint16[i]  = (input_float[i] > 0.0) ? 1 : 0;
        }

        /* obtain WM distance map */
        vbdist(GMT, input_uint16, sizes, voxelsize);
        for (i = 0; i < nvol; i++)
            dist_WM[i] = GMT[i];

        projection_based_thickness(input_float, dist_WM, dist_CSF, GMT, sizes, voxelsize); 

    }
                       
    for (scl_open = start_scl_open; scl_open > 0.4; scl_open -= 0.1)
    {
      
        /* Skip morphological opening if distopen is disabled */
        if (!distopen) scl_open = 1.0;
      
        /* We first apply a slightly different threshold for initial mask 
           to allow to control amount of morphological opening */
        for (i = 0; i < nvol; i++)
            input_uint8[i] = (double)input_float[i] >= scl_open*min_threshold ? 1 : 0;
    
        /* Interrupt here if distopen is disabled and use default scl_open value */
        if (!distopen) break;

        /* Optional morphological opening with distance criteria (distopen)
           Additionaly adapt dist parameter w.r.t. scl_open  */
        dist = 0.75/(scl_open*min_threshold);
        distopen_uint8(input_uint8, sizes, voxelsize, dist, 0.0);

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

    if (distopen && use_thickness && 0) {
        free(GMT);
        free(dist_CSF);
        free(dist_WM);
        free(input_uint16);
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
            
            median3_uint8(input_uint8, sizes);

            /* replace with its median filtered version outside sulcal areas */
            for (i = 0; i < nvol; i++)
                g0->output[i] = (unsigned short)(input_float[i] >= min_threshold / 1.0 ? input_uint8[i] : 0);
        }
        
        for (i = 0; i < nvol; i++)
            mcb->data[i]  = (double)(g0->output[i] - 0.5);

        /* extract surface to check euler number */
        run(mcb) ;
    
    /* convert mcb structure to BIC polygon data */
    polygons->n_items = mcb->ntrigs;
    polygons->n_points = mcb->nverts;
    ALLOC(polygons->points, polygons->n_points);
    ALLOC(polygons->normals, polygons->n_points);
    ALLOC(polygons->end_indices, polygons->n_items);
    polygons->bintree = (bintree_struct_ptr) NULL;
    
    for (i = 0; i < polygons->n_items; i++)
        polygons->end_indices[i] = (i + 1) * 3;
        
    ALLOC(polygons->indices,
        polygons->end_indices[polygons->n_items-1]);
        
    for (i = 0; i < polygons->n_points; i++) {
        real_voxel[0] = mcb->vertices[i].x;
        real_voxel[1] = mcb->vertices[i].y;
        real_voxel[2] = mcb->vertices[i].z;

        for(j= 0; j <N_DIMENSIONS; j++ )
        {
            if( spatial_axes[j] >= 0 )
                voxel_pos[j] = real_voxel[spatial_axes[j]];
            else
                voxel_pos[j] = 0.0;
        }
/*        general_transform_point( &voxel_to_world_transform,
                             voxel_pos[X], voxel_pos[Y], voxel_pos[Z],
                             &xw, &yw, &zw );

        Point_x(point) = xw;
        Point_y(point) = yw;
        Point_z(point) = zw;
*/
        Point_x(point) = voxel_pos[X];
        Point_y(point) = -voxel_pos[Y];
        Point_z(point) = voxel_pos[Z];
        polygons->points[i] = point;
    }
    
    for (i = 0; i < polygons->n_items; i++) {
        polygons->indices[POINT_INDEX(polygons->end_indices, i, 0)] = mcb->triangles[i].v3;
        polygons->indices[POINT_INDEX(polygons->end_indices, i, 1)] = mcb->triangles[i].v2;
        polygons->indices[POINT_INDEX(polygons->end_indices, i, 2)] = mcb->triangles[i].v1;
    }

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
    free(input_float);
    
    if (n_out > 2) fprintf(stderr,"Extract largest of %d components.\n",n_out);
    fprintf(stderr,"Euler characteristics after %d iterations is %d.\n", count, EC);

    /* Correct mesh in folded areas to compensate for the averaging effect in gyri and sulci.
       We use a folding measure (i.e. mean curvature averaged) to estimate the compensation. 
       The amount of compensation is automatically estimated using the difference to the defined 
       isovalue. */
    if (fwhm > 0.0) {
        smooth_heatkernel(polygons, NULL, fwhm);
//        correct_mesh_folding(polygons, NULL, volume, min_threshold);
    }

    (void) output_graphics_any_format(output_filename, ASCII_FORMAT, 1, &object3, NULL);

    free(input);
    clean_temps(mcb);

    return(0);
}



