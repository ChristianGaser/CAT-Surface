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

#define   CHUNK_SIZE      1000000

private void extract_isosurface(
        Volume                 volume,
        double                 min_label,
        double                 max_label,
        int                    spatial_axes[],
        General_transform      *voxel_to_world_transform,
        Marching_cubes_methods method,
        BOOLEAN                binary_flag,
        double                 min_threshold,
        double                 max_threshold,
        double                 valid_low,
        double                 valid_high,
        polygons_struct        *polygons);
        
private void extract_surface(
        Marching_cubes_methods method,
        BOOLEAN                binary_flag,
        double                 min_threshold,
        double                 max_threshold,
        double                 valid_low,
        double                 valid_high,
        int                    x_size,
        int                    y_size,
        double                 ***slices,
        double                 min_label,
        double                 max_label,
        double                 ***label_slices,
        int                    slice_index,
        BOOLEAN                right_handed,
        int                    spatial_axes[],
        General_transform      *voxel_to_world_transform,
        int                    ***point_ids[],
        polygons_struct        *polygons);

void distopen_ushort(unsigned short *vol, int dims[3], double voxelsize[3], int niter, double th);

static char *dimension_names_3D[] = { MIzspace, MIyspace, MIxspace };
static char *dimension_names[] = { MIyspace, MIxspace };

/* argument defaults */
double min_threshold = 0.5;
double fwhm = 3.0;
int    int_method = 1;
int    force_genus0 = 1;
int    n_opening = 0;

/* the argument table */
static ArgvInfo argTable[] = {
    {"-thresh", ARGV_FLOAT, (char *) TRUE, (char *) &min_threshold,
          "Marching cubes method: \n\t\t0 - marching cubes\n\t\t1 - marching cubes without holes\n\t\t2 - marching tetra"},
    {"-method", ARGV_INT, (char *) TRUE, (char *) &int_method,
          "Marching cubes method: \n\t\t0 - marching cubes\n\t\t1 - marching cubes without holes\n\t\t2 - marching tetra"},
    {"-force_genus0", ARGV_INT, (char *) TRUE, (char *) &force_genus0,
          "Allows to switch off genus0 functionality by setting to 0."},
    {"-fwhm", ARGV_FLOAT, (char *) TRUE, (char *) &fwhm,
          "Create a slightly smoothed surface that is corrected folded areas to compensate for the averaging effect in gyri and sulci.\n\
          Do not use smoothing sizes > 3 mm because that compensation only works reliably for smaller smoothing."},
    {"-n_opening", ARGV_INT, (char *) TRUE, (char *) &n_opening,
          "Additional number of opening iterations to prevent glued gyri."},
      {NULL, ARGV_END, NULL, NULL, NULL}
};


private void
usage(
        char *executable)
{
        char *usage_str = "\n\
Usage: CAT_MarchingCubesGenus0 input.nii output_surface_file threshold\n\
\n\
          Creates a polygonal surface of either the thresholded volume\n\
          and extracts the largest component.\n\n";

        print_error(usage_str, executable);
}

int   
main(
        int     argc,
        char    *argv[])
{
        char                      *input_filename, *output_filename;
        Volume                    volume;
        double                    min_label, max_label;
        double                    valid_low, valid_high, val;
        int                       i, j, k, c, spatial_axes[N_DIMENSIONS];
        int                       n_out, sizes[MAX_DIMENSIONS];
        Marching_cubes_methods    method;
        object_struct             *object, **object2, *object3;
        General_transform         voxel_to_world_transform;
        polygons_struct           *polygons;
        unsigned short            *input;

        /* get the arguments from the command line */
        if (ParseArgv(&argc, argv, argTable, 0)) {
                        usage(argv[0]);
                        fprintf(stderr, "         %s -help\n\n", argv[0]);
                        exit(EXIT_FAILURE);
        }

        initialize_argument_processing(argc, argv);

        if(!get_string_argument(NULL, &input_filename) ||
                !get_string_argument(NULL, &output_filename))
        {
                usage(argv[0]);
                fprintf(stderr, "         %s -help\n\n", argv[0]);
                return(1);
        }
        
        method = (Marching_cubes_methods) int_method;

        valid_low   =  0.0;
        valid_high =  -1.0;
        min_label   =  0.0;
        max_label   = -1.0;

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

        object  = create_object(POLYGONS);
        object3 = create_object(POLYGONS);

        if (force_genus0) {
                get_volume_sizes(volume, sizes);
                input = (unsigned short *)calloc(sizes[0]*sizes[1]*sizes[2],sizeof(unsigned short));    
        
                for (i = 0; i < sizes[0]; i++) {
                        for (j = 0; j < sizes[1]; j++) {
                                for (k = 0; k < sizes[2]; k++) {
                                        val = get_volume_real_value(volume, i, j, k, 0, 0);
                                        if (val > min_threshold)
                                                input[i + j*sizes[0] + k*sizes[0]*sizes[1]] = 1;
                                        else input[i + j*sizes[0] + k*sizes[0]*sizes[1]] = 0;
                                }
                        }
                }
        
                double voxelsize[3] = {1.0,1.0,1.0};
                
                /* We have to use a threshold of 0 because of data tpye (short) */
                if(n_opening > 0) fprintf(stderr,"Number of morphological openings: %d\n",n_opening);
                distopen_ushort(input, sizes, voxelsize, n_opening, 0);

                genus0parameters g0[1];   /* need an instance of genus0parameters */
        
                genus0init(g0);   /* initialize the instance, set default parameters */
        
                /* set some parameters/options */
                for(j= 0; j <N_DIMENSIONS; j++) g0->dims[j] = sizes[j];
                g0->connected_component = 1;
                g0->input = input;
                g0->cut_loops = 0;
                g0->connectivity = 6;
                g0->value = 1;
                g0->alt_value=1;
                g0->contour_value=1;
                g0->alt_contour_value=1;
                g0->any_genus = 0;
                g0->biggest_component = 1;
                g0->pad[0] = g0->pad[1] = g0->pad[2] = 2;
                g0->ijk2ras = NULL;
                g0->verbose = 0;
                g0->return_surface = 0;
                g0->extraijkscale[2] = 1;
        
                /* call the function! */
                if (genus0(g0)) return(1); /* check for error */
        
                /* save results as next input */
                for (i = 0; i < sizes[0]*sizes[1]*sizes[2]; i++)
                        input[i] = (unsigned int)g0->output[i];
        
                /* call genus0 a 2nd time with other parameters */
                g0->cut_loops = 1;
                g0->connectivity = 18;
                g0->value = 1;
                g0->alt_value=0;
                g0->contour_value=1;
                g0->alt_contour_value=0;
        
                if (genus0(g0)) return(1); 
        
                free(input);
        
                for (i = 0; i < sizes[0]; i++) {
                        for (j = 0; j < sizes[1]; j++) {
                                for (k = 0; k < sizes[2]; k++) {
                                        set_volume_real_value(volume, i, j, k, 0, 0, g0->output[i + j*sizes[0] + k*sizes[0]*sizes[1]]);
                                }
                        }
                }
        }
        extract_isosurface(volume,
                            min_label, max_label,
                            spatial_axes,
                            &voxel_to_world_transform,
                            method, FALSE,
                            0.5, 0.5,
                            valid_low, valid_high, get_polygons_ptr(object));

        check_polygons_neighbours_computed(get_polygons_ptr(object));
        n_out = separate_polygons(get_polygons_ptr(object), -1, &object2);

        if(n_out > 2) fprintf(stderr,"Extract largest of %d components.\n",n_out);
    
        triangulate_polygons(get_polygons_ptr(object2[0]), get_polygons_ptr(object3));
        polygons = get_polygons_ptr(object3);
        
        fprintf(stderr,"Euler characteristics is %d...\n", euler_characteristic(polygons));

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

        return(0);
}

/* estimate minimum of A and its index in A */
void
pmin(float A[], int sA, float *minimum, int *index)
{
        int i; 
        *minimum = FLT_MAX; *index = 0;
        for (i=0; i<sA; i++) {
                if ((A[i]>0.0) && (*minimum>A[i])) { 
                        *minimum = A[i]; 
                        *index      = i;
                }
        }
}

/* estimate x,y,z position of index i in an array size sx,sxy=sx*sy... */
void
ind2sub(int i,int *x,int *y, int *z, int sxy, int sy) {
        *z = (int)floor((double)i / (double)sxy) +1; 
        i = i % (sxy);
        *y = (int)floor((double)i / (double)sy) +1;             
        *x = i % sy + 1;
}

void
vbdist(float *V, unsigned int *IO, int *dims, double *voxelsize) 
{
    
        /* main information about input data (size, dimensions, ...) */
        const int         nvol = dims[0]*dims[1]*dims[2];
        const int         x   = dims[0];
        const int         y   = dims[1];
        const int         xy = x*y;
    
        float s1 = fabs((float)voxelsize[0]);
        float s2 = fabs((float)voxelsize[1]);
        float s3 = fabs((float)voxelsize[2]);
        const float     s12   = (float) sqrt((double)    s1*s1   + s2*s2); /* xy - voxel size */
        const float     s13   = (float) sqrt((double)    s1*s1   + s3*s3); /* xz - voxel size */
        const float     s23   = (float) sqrt((double)    s2*s2   + s3*s3); /* yz - voxel size */
        const float     s123 = (float) sqrt((double) s12*s12 + s3*s3); /* xyz - voxel size */
        
        /* indices of the neighbor Ni (index distance) and euclidean distance NW */
        const int     NI[] = { 0, -1,-x+1, -x,-x-1, -xy+1,-xy,-xy-1, -xy+x+1,-xy+x,-xy+x-1, -xy-x+1,-xy-x,-xy-x-1}; 
        const float   ND[] = {0.0, s1, s12, s2, s12, s13, s3,   s13, s123, s23, s123, s123, s23, s123};
        const int     sN = sizeof(NI)/4; /* division by 4 to get from the number of bytes to the number of elements */ 
        float DN[sN];
        float DNm = FLT_MAX;
        int   i, n, ni, DNi;
    
        /* data */
        float        *D;
        unsigned int *I;
        
        D = (float *)malloc(sizeof(float)*nvol);
        I = (unsigned int *)malloc(sizeof(unsigned int)*nvol);
        
        /* initialisation */
        for (i=0; i<nvol; i++) {
                if ((V[i]==0.0) || isnan(V[i])) D[i]=FLT_MAX; else D[i]=0.0; 
                I[i]=(unsigned int)i;
        }
    
        int u,v,w,nu,nv,nw; 
        /* forward direction that consider all points smaller than i */
        for (i=0; i<nvol; i++) {
            if (D[i]>0) {
                    ind2sub(i,&u,&v,&w,xy,x);
                    
                    /* read neighbor values */
                    for (n=0;n<sN;n++){
                            ni = i + NI[n];
                            ind2sub(ni,&nu,&nv,&nw,xy,x);
                            if ((ni<0) || (ni>=nvol) || (abs(nu-u)>1) || (abs(nv-v)>1) || (abs(nw-w)>1)) ni=i;
                            DN[n] = D[ni] + ND[n];
                    }
        
                    /* find minimum distance within the neighborhood */
                    pmin(DN,sN,&DNm,&DNi);
        
                    /* update values */
                    if (DNi>0) {
                            I[i] = (unsigned int)   I[i+NI[DNi]];
                            D[i] = DNm; 
                            ind2sub((int)I[i],&nu,&nv,&nw,xy,x); 
                            D[i] = sqrt(pow((float)(u-nu)*s1,2) + pow((float)(v-nv)*s2,2) + pow((float)(w-nw)*s3,2));
                    }
              }
        }
        
        /* backward direction that consider all points larger than i */
        for (i=nvol-1;i>=0;i--) {
                if (D[i]>0) {
                        ind2sub(i,&u,&v,&w,xy,x);
            
                        /* read neighbour values */
                        for (n=0;n<sN;n++) {
                                ni = i - NI[n];
                                ind2sub(ni,&nu,&nv,&nw,xy,x);
                                if ((ni<0) || (ni>=nvol) || (abs(nu-u)>1) || (abs(nv-v)>1) || (abs(nw-w)>1)) ni=i;
                                DN[n] = D[ni] + ND[n];
                        }
            
                        /* find minimum distance within the neighborhood */
                        pmin(DN,sN,&DNm,&DNi);
            
                        /* update values */
                        if (DNi>0) {
                                I[i] = (unsigned int)   I[i-NI[DNi]];
                                D[i] = DNm; 
                                ind2sub((int)I[i],&nu,&nv,&nw,xy,x); 
                                D[i] = sqrt(pow((float)(u-nu)*s1,2) + pow((float)(v-nv)*s2,2) + pow((float)(w-nw)*s3,2));
                        }
                }
        }
    
        for (i=0; i<nvol; i++) {
                V[i] = D[i];
                if (IO!=NULL) IO[i] = I[i] + 0;
        }
            
        free(D);
        free(I);
}

void
distopen_ushort(unsigned short *vol, int dims[3], double voxelsize[3], int niter, double th)
{
        float *buffer;
        int i,j;
        unsigned short max_vol;
        int nvol = dims[0]*dims[1]*dims[2];
        
        if (niter < 1) return;
    
        for (i=0; i<nvol; i++) max_vol = MAX(max_vol,vol[i]);
        th *= (double)max_vol;
        
        buffer = (float *)malloc(sizeof(float)*nvol);
    
        if(buffer == NULL) {
                fprintf(stderr,"Memory allocation error\n");
                exit(EXIT_FAILURE);
        }
        
        /* threshold input */
        for (i=0; i<nvol; i++)
                buffer[i] = 1.0 - (float)((double)vol[i]>th);
    
                    
        vbdist(buffer, NULL, dims, voxelsize);
    
        for (i=0; i<nvol; i++)
                buffer[i] = buffer[i] > (float)niter;
    
    
        vbdist(buffer, NULL, dims, voxelsize);
        for (i=0; i<nvol; i++)
                buffer[i] = buffer[i] < (float)niter;
    
    
        /* return image */
        for (i=0; i<nvol; i++)
                vol[i] = (unsigned short)buffer[i];
    
        free(buffer);
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
        int     spatial_axes[],
        General_transform *voxel_to_world_transform,
        Marching_cubes_methods    method,
        BOOLEAN binary_flag,
        double  min_threshold,
        double  max_threshold,
        double  valid_low,
        double  valid_high,
        polygons_struct *polygons)
{
        int     n_slices, sizes[MAX_DIMENSIONS], x_size, y_size, slice;
        int     ***point_ids[2], ***tmp_point_ids;
        int     max_edges;
        double **slices[2], **tmp_slices;
        double **label_slices[2];
        progress_struct progress;
        Surfprop spr;
        Point    point000, point100, point010, point001;
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
        General_transform     *voxel_to_world_transform,
        BOOLEAN  binary_flag,
        double   min_threshold,
        double   max_threshold,
        int      ***point_ids[],
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
        BOOLEAN                   binary_flag,
        double                    min_threshold,
        double                    max_threshold,
        double                    valid_low,
        double                    valid_high,
        int                       x_size,
        int                       y_size,
        double                    ***slices,
        double                    min_label,
        double                    max_label,
        double                    ***label_slices,
        int                       slice_index,
        BOOLEAN                   right_handed,
        int                       spatial_axes[],
        General_transform         *voxel_to_world_transform,
        int                       ***point_ids[],
        polygons_struct           *polygons)
{
        int               x, y, *sizes, tx, ty, tz, n_polys, ind;
        int               p, point_index, poly, size, start_points, dir;
        voxel_point_type  *points;
        double            corners[2][2][2], label;
        BOOLEAN           valid;

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

