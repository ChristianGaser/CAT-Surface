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
double fwhm   = 3.0;
double dist   = -1;
int median_correction  = 1;
int distopen  = 1;
int any_genus = 0;

/* the argument table */
static ArgvInfo argTable[] = {
  {"-thresh", ARGV_FLOAT, (char *) TRUE, (char *) &min_threshold,
      "Volume threshold (i.e. isovalue)."},
  {"-fwhm", ARGV_FLOAT, (char *) TRUE, (char *) &fwhm,
      "Create a slightly smoothed surface that is corrected in folded areas to compensate for the averaging effect in gyri and sulci.\n\
       Do not use smoothing sizes > 3 mm because that compensation only works reliably for smaller smoothing sizes."},
  {"-dist", ARGV_FLOAT, (char *) TRUE, (char *) &dist,
      "Manually define threshold for morphological opening (default -1 for automatic estimation)."},
  {"-no-median", ARGV_CONSTANT, (char *) FALSE, (char *) &median_correction,
      "Disable median filter that is used outside sulcal areas."},
  {"-no-distopen", ARGV_CONSTANT, (char *) FALSE, (char *) &distopen,
      "Disable additional morphological opening."},
  {"-no-genus0", ARGV_CONSTANT, (char *) TRUE, (char *) &any_genus,
      "Allows to switch off genus0 functionality."},
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
      1. Apply additional morphological opening (dist) to prevent glued gyri\n\
         and minimizes locl artifacts caused by small vessels or\n\
         mis-segmentations. The strength of the opening is automatically\n\
         estimated by analyzing the impact of different dist values\n\
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
    Volume          volume;
    double          min_label, max_label, start_dist;
    double          valid_low, valid_high, val, RMSE, sum_RMSE;
    double          voxelsize[N_DIMENSIONS];
    int             i, j, k, c, spatial_axes[N_DIMENSIONS];
    int             n_out, EC, sizes[MAX_DIMENSIONS];
    int             nvol, ind, count, stop_distopen;
    Marching_cubes_methods    method;
    object_struct       *object, **object2, *object3;
    General_transform   voxel_to_world_transform;
    polygons_struct     *polygons;
    unsigned short      *input;
    unsigned char       *input_uint8, *ref_uint8;
    float               *input_float;

    /* get the arguments from the command line */
    if (ParseArgv(&argc, argv, argTable, 0)) {
            usage(argv[0]);
            fprintf(stderr, "     %s -help\n\n", argv[0]);
            exit(EXIT_FAILURE);
    }

    initialize_argument_processing(argc, argv);

    if(!get_string_argument(NULL, &input_filename) ||
        !get_string_argument(NULL, &output_filename))
    {
        usage(argv[0]);
        fprintf(stderr, "     %s -help\n\n", argv[0]);
        return(1);
    }
    
    /* marching cubes without holes */
    method = (Marching_cubes_methods) 1;

    valid_low =  0.0;
    valid_high  = -1.0;
    min_label =  0.0;
    max_label = -1.0;

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

    input_float  = (float *)calloc(nvol,sizeof(float));
    for (i = 0; i < sizes[0]; i++) {
        for (j = 0; j < sizes[1]; j++) {
            for (k = 0; k < sizes[2]; k++) {
                ind = sub2ind(i,j,k,sizes[0],sizes[1]);
                input_float[ind]  = get_volume_real_value(volume, i, j, k, 0, 0);
            }
        }
    }
    
    input       = (unsigned short *) calloc(nvol,sizeof(unsigned short));  
    input_uint8 = (unsigned char  *) calloc(nvol,sizeof(unsigned char));  
    ref_uint8   = (unsigned char  *) calloc(nvol,sizeof(unsigned char));  

    /* We analyze the impact of different dist values ranging from 1.8 to 0.5 
       using distopen. By doing so, we can track the changes in RMSE (Root Mean 
       Square Error) between these values. When using large dist values, the 
       changes in RMSE are relatively significant and stable because they affect 
       the whole image resulting in smaller gyri and wider sulci. In contrast, 
       using smaller dist values results in only local changes that are much 
       smaller since only regions with artifacts such as vessels or poor skull-
       stripping are affected. The optimal dist value is determined by identifying 
       the point where the RMSE decreases significantly, indicating successful 
       artifact removal while maintaining global gyri and sulci characteristics. 
    */

    /* default dist is negative that indicates automatically search for 
       optimal dist value
    */
    if (dist < 0) {
        stop_distopen = 0;
        start_dist = 1.8;
    } else {
        stop_distopen = 1;
        start_dist = dist;
    }
    
    sum_RMSE = 0.0;
    count  = 0;

    if (!distopen) fprintf(stderr,"%5s\t%5s\t%5s\n","Dist","avgRMSE","RMSE");                       
    for (dist = start_dist; dist > 0.4; dist -= 0.1) {
      
        /* skip morphological opening if distopen is disabled */
        if (!distopen) dist = 1.0;
      
        /* we first apply a slightly lower threshold for initial mask 
           to allow to control amount of morphological opening */
        for (i = 0; i < nvol; i++) {            
            if ((double)input_float[i] >= dist*min_threshold)
                input_uint8[i] = 1;
            else  input_uint8[i] = 0;
        }
    
        /* interrupt here if distopen is disabled */
        if (!distopen) break;

        /* optional morphological opening with distance criteria (distopen)
           use fixed dist of 1.5 because we control strength with the 
           previous step above. */
        distopen_uint8(input_uint8, sizes, voxelsize, 1.5, 0.0);

        for (i = 0; i < nvol; i++) {            
            if ((double)input_float[i] >= min_threshold) {
                if (input_uint8[i] == 0)
                    input_uint8[i] = 0;
                else  input_uint8[i] = 1;
            } else input_uint8[i] = 0;
        }
        
        /* stop after one additional iteration */
        if (stop_distopen > 0) break;
        
        /* Calulate RMSE between actual and previous distopen and stop if
           changes in RMSE are getting much smaller to obtain the optimal 
           dist parameter */
        if (count) {
            RMSE = 0.0;
            for (i = 0; i < nvol; i++) {
                val = (double)input_uint8[i] - (double)ref_uint8[i];  
                RMSE += val*val;
            } 
            RMSE = sqrt(RMSE/(double)nvol);
            sum_RMSE += RMSE;
            
            /* indicate stop if changes are getting smaller by a factor of 1.5 */
            if (sum_RMSE/RMSE/(double)count > 1.5) {
                fprintf(stderr,"%5.4f\t%5.4f\t%5.4f\n",dist,sum_RMSE/RMSE/(double)count,RMSE);    
                break;    
            }
        }
        
        /* save previous image after distopen */ 
        for (i = 0; i < nvol; i++) ref_uint8[i] = input_uint8[i];   

        count++;        

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
            for (i = 0; i < nvol; i++) {
                if (input_float[i] > min_threshold/1.0)
                    g0->output[i] = (unsigned short)input_uint8[i];
            }
        }
        
        for (i = 0; i < sizes[0]; i++) {
            for (j = 0; j < sizes[1]; j++) {
                for (k = 0; k < sizes[2]; k++) {
                    ind = sub2ind(i,j,k,sizes[0],sizes[1]);
                    set_volume_real_value(volume, i, j, k, 0, 0, g0->output[ind]);                    
                }
            }
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
    
    if(n_out > 2) fprintf(stderr,"Extract largest of %d components.\n",n_out);
    fprintf(stderr,"Euler characteristics after %d iterations is %d.\n", count, EC);

    /* Correct mesh in folded areas to compensate for the averaging effect in gyri and sulci.
       We use a folding measure (i.e. mean curvature averaged) to estimate the compensation. 
       The amount of compensation is automatically estimated using the difference to the defined 
       isovalue. */
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
        slice[x][y] = 0.0;
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
        slice[x][y] = get_volume_real_value(volume, x, y, z, 0, 0);
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
    if(x < 0 || x >= x_size || y < 0 || y >= y_size)
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
    {
        for (y = 0; y < y_size+2; y++)
        {
            for (edge = 0; edge < max_edges; edge++)
                point_ids[x][y][edge] = -1;
        }
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
    {
        if(spatial_axes[c] >= 0)
            voxel_pos[c] = real_voxel[spatial_axes[c]];
        else
            voxel_pos[c] = 0.0;
    }

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
        if(slice < n_slices - 1)
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

    if(polygons->n_points > 0)
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
    if(point_index < 0)
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
    {
        for (y = -1; y < y_size; y++)
        {
            valid = TRUE;
            for (tx = 0; tx < 2; tx++)
            for (ty = 0; ty < 2; ty++)
            for (tz = 0; tz < 2; tz++)
            {
                corners[tx][ty][tz] = get_slice_value(slices, x_size, y_size,
                    tz, x + tx, y + ty);
                if(valid_low <= valid_high &&
                       (corners[tx][ty][tz] < min_threshold ||
                      corners[tx][ty][tz] > max_threshold) &&
                       (corners[tx][ty][tz] < valid_low ||
                      corners[tx][ty][tz] > valid_high))
                    valid = FALSE;

                if(min_label <= max_label)
                {
                    label = get_slice_value(label_slices, x_size, y_size,
                              tz, x + tx, y + ty);
                    if(label < min_label || label > max_label)
                        corners[tx][ty][tz] = 0.0;
                }
            }

            if(!valid)
                continue;

            n_polys = compute_isosurface_in_voxel(method, x, y, slice_index,
                              corners, binary_flag, min_threshold,
                              max_threshold, &sizes, &points);

            if(n_polys == 0)
                continue;

            if(right_handed)
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

                if(right_handed)
                    start_points += size;
                else if(poly < n_polys-1)
                    start_points += sizes[poly+1];
            }
        }
    }
}

