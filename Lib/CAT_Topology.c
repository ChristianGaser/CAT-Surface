/* Rachel Yotter - rachel.yotter@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 */

#include <bicpl.h>
#include <bicpl/deform.h>

#include "CAT_SPH.h"
#include "CAT_Surf.h"
#include "CAT_Intersect.h"
#include "CAT_SurfaceIO.h"
#include "CAT_Curvature.h"
#include "CAT_Defect.h"
#include "CAT_DeformPolygons.h"
#include "CAT_Refine.h"

#define DATAFORMAT 1 /* 1 = real data, 0 = complex data */
#define DEBUG 0
#define DUMP_FILES 0

#define index(A,B,C,DIM) ((C)*DIM[0]*DIM[1] + (B)*DIM[0] + (A))

#ifndef MIN
#define MIN(A,B) ((A) > (B) ? (B) : (A))
#endif

#ifndef MAX
#define MAX(A,B) ((A) > (B) ? (A) : (B))
#endif


/* based on caret code for TFCE calculation in BrainModelVolumeTFCE.cxx */
void get_cluster_size(unsigned int *bw, int dims[3])
{
        unsigned int valToAdd, *bw2;
        int i, j, k, ti, tj, tk, maxi, maxj, maxk, mini, minj, mink, ind, ind1, growingInd, growingCur;
        int numVoxels = dims[0] * dims[1] * dims[2];
        char *flagUsed;
        short *growing;
         
        bw2 = (unsigned int*) malloc(numVoxels*sizeof(unsigned int));
        flagUsed = (char*) malloc(numVoxels*sizeof(char));
        growing  = (short*)malloc(numVoxels*3*sizeof(short));

        for (i = 0; i < numVoxels; ++i) flagUsed[i] = 0;
                
        for (k = 0; k < dims[2]; ++k) for (j = 0; j < dims[1]; ++j) for (i = 0; i < dims[0]; ++i)
        {
                ind = k*(dims[0]*dims[1]) + (j*dims[0]) + i;
                                                
                if (!flagUsed[ind] && bw[ind] > 0)
                {
                        flagUsed[ind] = 1;
                        growingInd = 3;
                        growingCur = 0;
                        growing[0] = i;
                        growing[1] = j;
                        growing[2] = k;
                        
                        while (growingCur < growingInd)
                        {
                                maxi = MIN(dims[0], growing[growingCur          ] + 2);
                                maxj = MIN(dims[1], growing[growingCur + 1] + 2);
                                maxk = MIN(dims[2], growing[growingCur + 2] + 2);
                                
                                mini = MAX(0, growing[growingCur                ] - 1);
                                minj = MAX(0, growing[growingCur + 1] - 1);
                                mink = MAX(0, growing[growingCur + 2] - 1);
                                
                                for (tk = mink; tk < maxk; ++tk) for (tj = minj; tj < maxj; ++tj) for (ti = mini; ti < maxi; ++ti)
                                {
                                        ind1 = tk*(dims[0]*dims[1]) + (tj*dims[0]) + ti;
                                        
                                        if (!flagUsed[ind1] && bw[ind1] > 0)
                                        {
                                                flagUsed[ind1] = 1;
                                                growing[growingInd              ] = ti;
                                                growing[growingInd + 1] = tj;
                                                growing[growingInd + 2] = tk;
                                                growingInd += 3;
                                        }
                                }
                                growingCur += 3;
                        }
                        
                        growingCur = 0;
                        valToAdd = (unsigned int)(growingInd / 3.0);

                        while (growingCur < growingInd)
                        {
                                bw2[growing[growingCur + 2]*(dims[0]*dims[1])+(growing[growingCur + 1]*dims[0])+growing[growingCur]] = valToAdd;
                                growingCur += 3;
                        }
                }                
        }
         
        /* copy cluster size values to input */
        for (i = 0; i < numVoxels; ++i) bw[i] = bw2[i];

        free(flagUsed);
        free(growing);
        free(bw2);
}

void
add_neighbours(polygons_struct *surface, polygons_struct *lbw,
               int **neighbours, int *n_neighbours, int p, int *flag, int level)
{
        int n, idx;
        double dist;

        if (level == 0) {
                for (n = 0; n < n_neighbours[p]; n++) {
                        idx = neighbours[p][n];

                        if (flag[idx] == 1)
                                continue;

                        dist = distance_between_points(&surface->points[idx],
                                                       &lbw->points[idx]);
                        if (dist > 2.0) { /* prevent discontinuities */
                                flag[idx] = 1;
                                add_neighbours(surface, lbw, neighbours,
                                               n_neighbours, idx, flag, 0);
                        }
                }
                return;
        }

        for (n = 0; n < n_neighbours[p]; n++) {
                idx = neighbours[p][n];

                if (flag[idx] == 1)
                        continue;

                flag[idx] = 1;
                add_neighbours(surface, lbw, neighbours, n_neighbours, idx,
                               flag, level - 1);
        }
}

void
resample_defects_sph(polygons_struct *sphere, int *defects, int *polydefects,
                     int *remap_defects, int *remap_polydefects, int n_items)
{
        object_struct **objects;
        polygons_struct *remap;
        Point centre;

        objects = (object_struct **) malloc(sizeof(object_struct *));
        *objects = create_object(POLYGONS);
        remap = get_polygons_ptr(*objects);
        initialize_polygons(remap, WHITE, NULL);

        fill_Point(centre, 0.0, 0.0, 0.0);
        create_tetrahedral_sphere(&centre, 1.0, 1.0, 1.0, n_items, remap);

        memset(remap_defects, 0, sizeof(int) * remap->n_points);

        remap_defect(sphere, defects, polydefects, remap, remap_defects,
                     remap_polydefects);

        delete_object_list(1, objects);
}

void
changed_voxels_between_surfaces(object_struct *object1, object_struct *object2, double *changes)
{
        Volume    label_volume1, volume1, label_volume2;        
        int       i, label, sizes[3], x, y, z;
        int       range_changed[2][3], value[3];
        int       *n_neighbours, **neighbours, *flag;
        double    separations[3], voxel[3], world[3];
        double    bounds[6], value1, value2;
        unsigned int *cluster;
        
        /* create volume inside surface bounds */
        get_bounds(get_polygons_ptr(object1), bounds);
        volume1  = create_volume( 3, XYZ_dimension_names, NC_BYTE, FALSE,
                            0.0, 255.0 );
        
        /* prepare volume parameters */
        for (i=0; i<3; i++) {
                separations[i] = 0.75;
                voxel[i] = -0.5;
                world[i] = bounds[2*i] - separations[i];
                value[i] = 2.0;
                sizes[i] = ROUND((bounds[2*i+1] - bounds[2*i]) / separations[i]) + 2;
        }

        cluster = (unsigned int *) malloc(sizeof(unsigned int) * sizes[0]*sizes[1]*sizes[2]);

        set_volume_separations( volume1, separations );    
        set_volume_sizes( volume1, sizes );
        set_volume_voxel_range( volume1, 0.0, 255.0 );
        set_volume_real_range( volume1, 0, 255.0 );
        set_volume_translation( volume1, voxel, world );

        alloc_volume_data( volume1 );
        
        /* label volume according to surface */
        label_volume1 = create_label_volume( volume1, NC_SHORT );
        label_volume2 = create_label_volume( volume1, NC_SHORT );
        scan_object_to_volume( object1, volume1, label_volume1, 1, 0.0 );
        fill_connected_voxels( volume1, label_volume1, EIGHT_NEIGHBOURS,
                           value, 0, 0, 1, 0.0, -1.0, range_changed );
        scan_object_to_volume( object2, volume1, label_volume2, 1, 0.0 );
        fill_connected_voxels( volume1, label_volume2, EIGHT_NEIGHBOURS,
                           value, 0, 0, 1, 0.0, -1.0, range_changed );

        if (DUMP_FILES) { 
                (void) output_volume( "label_volume1.mnc", NC_SHORT, TRUE,
                          0.0, 0.0,
                          label_volume1, "label1\n", NULL );

                (void) output_volume( "label_volume2.mnc", NC_SHORT, TRUE,
                          0.0, 0.0,
                          label_volume2, "label2\n", NULL );
        }
        
        for (x = 0; x < sizes[0]; x++)
          for (y = 0; y < sizes[1]; y++)
            for (z = 0; z < sizes[2]; z++) {
        
                value1 = get_volume_real_value( label_volume1, x, y, z, 0, 0 );
                value2 = get_volume_real_value( label_volume2, x, y, z, 0, 0 );

                if (value1 != value2) 
                        cluster[index(x,y,z,sizes)] = 255;
                else    cluster[index(x,y,z,sizes)] = 0;
        }

        get_cluster_size(cluster, sizes);

        for (x = 0; x < sizes[0]; x++)
          for (y = 0; y < sizes[1]; y++)
            for (z = 0; z < sizes[2]; z++) {
                        set_volume_real_value( label_volume1, x, y, z, 0, 0, (double)cluster[index(x,y,z,sizes)]);
        }

        for (i = 0; i < get_polygons_ptr(object1)->n_points; i++) {
                evaluate_volume_in_world(label_volume1,
                                         RPoint_x(get_polygons_ptr(object1)->points[i]),
                                         RPoint_y(get_polygons_ptr(object1)->points[i]),
                                         RPoint_z(get_polygons_ptr(object1)->points[i]),
                                         0, FALSE, 0.0, &value1,
                                         NULL, NULL, NULL, NULL, NULL,
                                         NULL, NULL, NULL, NULL);
                changes[i] = value1;
        }

        if (DUMP_FILES) 
                (void) output_volume( "cluster.mnc", NC_SHORT, TRUE,
                          0.0, 0.0,
                          label_volume1, "diff\n", NULL );

        
        delete_volume( volume1 );
        delete_volume( label_volume1 );
        delete_volume( label_volume2 );
        free(cluster);
}

void
surface_deform(object_struct *object, polygons_struct *hbw, int *hbw_defects)
{
        Volume    volume, label_volume, tmp;        
        int       i, label, sizes[3];
        int       range_changed[2][3], value[3];
        int       *n_neighbours, **neighbours, *flag;
        double    separations[3], voxel[3], world[3];
        double    bounds[6], label_val, threshold;
        deform_struct deform;
        polygons_struct *hbw_restore;
        
        if (hbw_defects != NULL) {
                hbw_restore = (polygons_struct *) malloc(sizeof(polygons_struct));
                copy_polygons(hbw, hbw_restore);
        }
        
        /* create volume inside surface bounds */
        get_bounds(get_polygons_ptr(object), bounds);
        volume = create_volume( 3, XYZ_dimension_names, NC_BYTE, FALSE,
                            0.0, 255.0 );
        
        /* prepare volume parameters */
        for (i=0; i<3; i++) {
                separations[i] = 0.75;
                voxel[i] = -0.5;
                world[i] = bounds[2*i] - separations[i];
                value[i] = 2.0;
                sizes[i] = ROUND((bounds[2*i+1] - bounds[2*i]) / separations[i]) + 2;
        }

        set_volume_separations( volume, separations );    
        set_volume_sizes( volume, sizes );
        set_volume_voxel_range( volume, 0.0, 255.0 );
        set_volume_real_range( volume, 0, 255.0 );
        set_volume_translation( volume, voxel, world );

        alloc_volume_data( volume );
        
        /* label volume according to surface */
        label_val = 127.0;
        label_volume = create_label_volume( volume, NC_BYTE );
        scan_object_to_volume( object, volume, label_volume, (int) label_val, 0.0 );
        
        /* fill inside volume */
        fill_connected_voxels( volume, label_volume, EIGHT_NEIGHBOURS,
                           value, 0, 0, label_val*2.0, 0.0, -1.0, range_changed );
                           
        tmp = create_box_filtered_volume(label_volume, NC_BYTE, FALSE,
                                                 0.0, 0.0, 2, 2, 2);

        delete_volume(label_volume);
        label_volume = tmp;

        initialize_deformation_parameters(&deform);
        deform.fractional_step = 0.1;
        deform.max_step = 0.1;
        deform.max_search_distance = 5;
        deform.degrees_continuity = 0;
        deform.max_iterations = 100;
        deform.movement_threshold = 0.01;
        deform.stop_threshold = 0.0;
        deform.deform_data.type = VOLUME_DATA;
        deform.deform_data.volume = label_volume;
        deform.deform_data.label_volume = (Volume) NULL;

        if (add_deformation_model(&deform.deformation_model, -1, 0.5, "avg", -0.1, 0.1) != OK)
                exit(EXIT_FAILURE);

        threshold = 0.75*label_val;
        set_boundary_definition(&deform.boundary_definition, threshold, threshold, 0, 0, 'n', 0);

        if (DEBUG) fprintf(stderr,"deform_polygons_points...\n");
        deform_polygons(hbw, &deform);

        if (hbw_defects != NULL) {
                flag = (int *) malloc(sizeof(int) * hbw->n_points);
                memset(flag, 0, sizeof(int) * hbw->n_points);
                create_polygon_point_neighbours(hbw, TRUE, &n_neighbours,
                                        &neighbours, NULL, NULL);
                for (i = 0; i < hbw->n_points; i++) {
                        if  (hbw_defects[i] != 0) {
                                flag[i] = 1; /* patch this one */
                                add_neighbours(hbw, hbw_restore, neighbours, n_neighbours, i,
                                       flag, 1);
                        }
                }

                /* combine the surfaces based on the modify flag */
                for (i = 0; i < hbw->n_points; i++) {
                        if (flag[i] == 0)
                                hbw->points[i] = hbw_restore->points[i]; /* restore */
                }
                free(flag);
                free(hbw_restore);
        }
        
        delete_volume( volume );
        delete_volume( label_volume );
}

void
sph_postcorrect(polygons_struct *surface, polygons_struct *sphere, int *defects, int *polydefects, 
                int n_defects, int *holes, polygons_struct *hbw, polygons_struct *lbw, int do_surface_deform)
{
        object_struct *surface_object;
        double *sharpness;
        int *n_neighbours, **neighbours;
        int p, d, val, poly, *flag, smooth_iters;
        int *hbw_defects, *hbw_polydefects, *hbw_holes;

        /* remap the defects onto the spherical harmonic reconstruction */
        create_polygon_point_neighbours(hbw, TRUE, &n_neighbours,
                                        &neighbours, NULL, NULL);

        if (DEBUG) fprintf(stderr,"resample_defects_sph...\n");
        hbw_defects =     (int *) malloc(sizeof(int) * hbw->n_points);
        hbw_polydefects = (int *) malloc(sizeof(int) * hbw->n_items);
        resample_defects_sph(sphere, defects, polydefects,
                             hbw_defects, hbw_polydefects, hbw->n_items);
        expand_defects(hbw, hbw_defects, hbw_polydefects, 0, 2,
                       n_neighbours, neighbours);

        /* remap the hole flags */
        hbw_holes = (int *) malloc(sizeof(int) * hbw->n_points);
        memset(hbw_holes, 0, sizeof(int) * hbw->n_points);
        for (d = 1; d <= n_defects; d++) {
                for (p = 0; p < sphere->n_points; p++) {
                        if (defects[p] == d && holes[p] != 0) {
                                val = holes[p];
                                break;
                        }
                }
 
                for (p = 0; p < hbw->n_points; p++) {
                        if (hbw_defects[p] == d)
                                hbw_holes[p] = val;
                }
        }

        /* optionally output the defects list */
        if (DUMP_FILES) {
                output_values_any_format("hbw_defects.txt", hbw->n_points,
                                         hbw_defects, TYPE_INTEGER);
        }

        if (DEBUG) fprintf(stderr,"compute_local_sharpness...\n");
        sharpness = (double *) malloc(sizeof(double) * hbw->n_points);
        compute_local_sharpness(hbw, n_neighbours, neighbours, sharpness);

        if (DEBUG) fprintf(stderr,"patch_points...\n");
        flag = (int *) malloc(sizeof(int) * hbw->n_points);
        memset(flag, 0, sizeof(int) * hbw->n_points);
        for (p = 0; p < hbw->n_points; p++) {
                if (hbw_holes[p] == VENTRICLE || hbw_holes[p] == LARGE_DEFECT) {
                        flag[p] = 1;
                } else if (sharpness[p] > 60 && hbw_defects[p] != 0) {
                        flag[p] = 1; /* patch this one */
                        add_neighbours(hbw, lbw, neighbours, n_neighbours, p,
                                       flag, 1);
                }
        }

        /* combine the surfaces based on the modify flag */
        for (p = 0; p < hbw->n_points; p++) {
                if (flag[p] == 1)
                        hbw->points[p] = lbw->points[p];
        }

        /* find remaining self-intersections */
        n_defects = find_selfintersections(hbw, hbw_defects, hbw_polydefects);
        n_defects = join_intersections(hbw, hbw_defects, hbw_polydefects,
                                       n_neighbours, neighbours);

        fprintf(stderr,"%d self intersection(s) to repair\n", n_defects);

        /* patch self-intersections first */
        if (DEBUG) fprintf(stderr, "patch_selfintersections...\n");
        n_defects = patch_selfintersections(hbw, lbw, hbw_defects,
                                            hbw_polydefects, n_defects,
                                            n_neighbours, neighbours);
        fprintf(stderr,"Post-patch: %d self intersection(s) remaining\n", n_defects);

        /* smooth out remaining self-intersections */
        smooth_iters = 200;
        if (DEBUG) fprintf(stderr,"smooth_selfintersections with %d iterations...\n", smooth_iters);
        n_defects = smooth_selfintersections(hbw, hbw_defects, hbw_polydefects,
                                             n_defects, n_neighbours,
                                             neighbours, smooth_iters);


        if (DUMP_FILES) { 
                output_values_any_format("hbw_modpts.txt", hbw->n_points, flag,
                                 TYPE_INTEGER);
        }

        if ( do_surface_deform ) {
                surface_object = create_object(POLYGONS);
                copy_polygons(surface, get_polygons_ptr(surface_object));
                surface_deform(surface_object, hbw, hbw_holes);
        }

        if (DEBUG) fprintf(stderr,"compute_polygon_normals...\n");
        compute_polygon_normals(hbw);

        free(sharpness);
        free(flag);
        free(hbw_defects);
        free(hbw_polydefects);
        free(hbw_holes);
}

object_struct **
fix_topology_sph(polygons_struct *surface, polygons_struct *sphere, int n_triangles, int bw, int lim, 
        char *reparam_file, double max_refine_length, int do_surface_deform, int force)
{
        object_struct **hbw_objects, **lbw_objects, **reparam_objects, *surface_object;
        polygons_struct *hbw, *lbw, *reparam, refined;
        Point *length_points;
        double *rcx, *icx, *rcy, *icy, *rcz, *icz;
        double *lrcx, *licx, *lrcy, *licy, *lrcz, *licz;
        double *rdatax, *rdatay, *rdataz;
        int bw2, i, p, d, n_done;
        int *defects, *polydefects, *holes, n_defects, n_objects;
        int *n_neighbours, **neighbours;
        double *HD_hole, *HD_handle, *curv;
        double sum_HD_handle, sum_HD_hole;
        File_formats format;

        defects     = (int *) malloc(sizeof(int) * surface->n_points);
        holes       = (int *) malloc(sizeof(int) * surface->n_points);
        polydefects = (int *) malloc(sizeof(int) * surface->n_items);
        HD_handle   = (double *) malloc(sizeof(double) * surface->n_points);
        HD_hole     = (double *) malloc(sizeof(double) * surface->n_points);
        curv        = (double *) malloc(sizeof(double) * surface->n_points);
        hbw_objects = (object_struct **) malloc(sizeof(object_struct *));

        /* find defects in original uncorrected surface */
        create_polygon_point_neighbours(sphere, TRUE, &n_neighbours,
                                        &neighbours, NULL, NULL);

        if (DEBUG) fprintf(stderr,"find_topological_defects...\n");
        n_defects = find_topological_defects(surface, sphere, defects,
                                             n_neighbours, neighbours);
        fprintf(stderr,"%d topological defects\n", n_defects);

        update_polydefects(surface, defects, polydefects);

        if (reparam_file != NULL) {
                if (input_graphics_any_format(reparam_file, &format, &n_objects,
                                              &reparam_objects) != OK)
                        exit(EXIT_FAILURE);

                /* check that the surface file contains a polyhedron */
                if (n_objects != 1 ||
                    get_object_type(reparam_objects[0]) != POLYGONS) {
                        fprintf(stderr,"Reparam sphere file must contain 1 polygons object.\n");
                        exit(EXIT_FAILURE);
                }
                reparam = get_polygons_ptr(reparam_objects[0]);
                for (p = 0; p < reparam->n_points; p++)
                        set_vector_length(&reparam->points[p], 1.0);
                n_triangles = reparam->n_items;
        } else reparam = NULL;
        
        bw2 = bw * bw;

        rdatax = (double *) malloc(sizeof(double) * 4 * bw2);
        rdatay = (double *) malloc(sizeof(double) * 4 * bw2);
        rdataz = (double *) malloc(sizeof(double) * 4 * bw2);
        rcx    = (double *) malloc(sizeof(double) * bw2);
        rcy    = (double *) malloc(sizeof(double) * bw2);
        rcz    = (double *) malloc(sizeof(double) * bw2);
        icx    = (double *) malloc(sizeof(double) * bw2);
        icy    = (double *) malloc(sizeof(double) * bw2);
        icz    = (double *) malloc(sizeof(double) * bw2);
        lrcx   = (double *) malloc(sizeof(double) * bw2);
        lrcy   = (double *) malloc(sizeof(double) * bw2);
        lrcz   = (double *) malloc(sizeof(double) * bw2);
        licx   = (double *) malloc(sizeof(double) * bw2);
        licy   = (double *) malloc(sizeof(double) * bw2);
        licz   = (double *) malloc(sizeof(double) * bw2);

        if (!force) {
            /* label defects first always as handles and in the 2nd step as holes */
            for (p = 0; p < surface->n_points; p++) {
                    if (defects[p] > 0) 
                            holes[p] = HANDLE; 
                    else    holes[p] = 0;
            }
            
            if (DEBUG) fprintf(stderr,"get_equally_sampled_coords_holes (handles)...\n");
            get_equally_sampled_coords_holes(surface, sphere, defects, n_defects,
                                             holes, bw, rdatax, rdatay, rdataz);
    
            if (DEBUG) fprintf(stderr,"get_sph_coeffs_of_realdata (handles)...\n");
            get_sph_coeffs_of_realdata(rdatax, bw, DATAFORMAT, rcx, icx);
            get_sph_coeffs_of_realdata(rdatay, bw, DATAFORMAT, rcy, icy);
            get_sph_coeffs_of_realdata(rdataz, bw, DATAFORMAT, rcz, icz);
    
            if (DEBUG) fprintf(stderr,"sample_sphere_from_sph (handles)...\n");
            *hbw_objects = create_object(POLYGONS);
            hbw = get_polygons_ptr(*hbw_objects);
            sample_sphere_from_sph(rdatax, rdatay, rdataz, hbw,
                                   n_triangles, reparam, bw);
    
            surface_object = create_object(POLYGONS);
            copy_polygons(surface, get_polygons_ptr(surface_object));
    
            get_polygon_vertex_curvatures_cg(surface, n_neighbours, neighbours,
                                             3.0, 0, curv);
    
            /* calculate Hausdorff distance between SPH-reparameterized and original surface */ 
            compute_point_hausdorff(surface, hbw, HD_handle, 0);
            
            /* multiply Hausdorff distance and curvature to differentiate between gyri and sulci */
            for (p = 0; p < surface->n_points; p++) 
                    HD_handle[p] *= curv[p];        
    
            if (DUMP_FILES) {
                    output_graphics_any_format("hbw_handle.obj", ASCII_FORMAT, 1,
                                               hbw_objects, NULL);
                    output_values_any_format("HD_handle.txt", surface->n_points,
                                             HD_handle, TYPE_DOUBLE);
            }
    
            /* label all defects as holes */
            for (p = 0; p < surface->n_points; p++) {
                    if (defects[p] > 0) 
                            holes[p] = HOLE; 
                    else    holes[p] = 0;
            }
        
            if (DEBUG) fprintf(stderr,"get_equally_sampled_coords_holes (holes)...\n");
            get_equally_sampled_coords_holes(surface, sphere, defects, n_defects,
                                             holes, bw, rdatax, rdatay, rdataz);
    
            if (DEBUG) fprintf(stderr,"get_sph_coeffs_of_realdata  (holes)...\n");
            get_sph_coeffs_of_realdata(rdatax, bw, DATAFORMAT, rcx, icx);
            get_sph_coeffs_of_realdata(rdatay, bw, DATAFORMAT, rcy, icy);
            get_sph_coeffs_of_realdata(rdataz, bw, DATAFORMAT, rcz, icz);
    
            if (DEBUG) fprintf(stderr,"sample_sphere_from_sph  (holes)...\n");
            *hbw_objects = create_object(POLYGONS);
            hbw = get_polygons_ptr(*hbw_objects);
            sample_sphere_from_sph(rdatax, rdatay, rdataz, hbw,
                                   n_triangles, reparam, bw);
    
            /* calculate Hausdorff distance between SPH-reparameterized and original surface */ 
            compute_point_hausdorff(surface, hbw, HD_hole, 0);
    
            /* multiply Hausdorff distance and curvature to differentiate between gyri and sulci */
            for (p = 0; p < surface->n_points; p++) 
                    HD_hole[p] *= curv[p];
    
            if (DUMP_FILES) {
                    output_graphics_any_format("hbw_hole.obj", ASCII_FORMAT, 1,
                                               hbw_objects, NULL);
                    output_values_any_format("HD_hole.txt", surface->n_points,
                                             HD_hole, TYPE_DOUBLE);
            }
            
            /* label defects as handles or holes depending on their hausdorff distance multiplied by curvature inside the defect */ 
            for (d = 1; d < n_defects+1; d++) {
                    sum_HD_handle = 0.0;
                    sum_HD_hole   = 0.0;
                    int defect_size = 0.0;
                    for (p = 0; p < surface->n_points; p++) {
                            if (defects[p] == d) {
                                    sum_HD_handle += HD_handle[p];
                                    sum_HD_hole   += HD_hole[p];
                                    defect_size++;
                            }
                    } 
                    if (DEBUG) fprintf(stderr,"%d %d %g %g\n",d,defect_size,sum_HD_handle,sum_HD_hole);
                    for (p = 0; p < surface->n_points; p++) {
                            if (defects[p] == d) {
                                    /* use handles or holes depending on hausdorff 
                                       distance multiplied by curvature between surfaces */
                                    if (sum_HD_handle >= sum_HD_hole)
                                            holes[p] = HANDLE;
                                    else    holes[p] = HOLE;
                            }
                    }
            }
            
        } else {
            /* label defects according to force parameter to either use 
                holes or handles by default */
            for (p = 0; p < surface->n_points; p++) {
                    if (defects[p] > 0) 
                            holes[p] = force; 
                    else    holes[p] = 0;
            }
        }

        if (DUMP_FILES) {
                output_values_any_format("defects.txt", surface->n_points,
                                         holes, TYPE_INTEGER);
        }

        if (DEBUG) fprintf(stderr,"get_equally_sampled_coords_holes (final)...\n");
        get_equally_sampled_coords_holes(surface, sphere, defects, n_defects,
                                         holes, bw, rdatax, rdatay, rdataz);

        if (DEBUG) fprintf(stderr,"get_sph_coeffs_of_realdata (final)...\n");
        get_sph_coeffs_of_realdata(rdatax, bw, DATAFORMAT, rcx, icx);
        get_sph_coeffs_of_realdata(rdatay, bw, DATAFORMAT, rcy, icy);
        get_sph_coeffs_of_realdata(rdataz, bw, DATAFORMAT, rcz, icz);
        
        if (DEBUG) fprintf(stderr,"sample_sphere_from_sph (final)...\n");
        *hbw_objects = create_object(POLYGONS);
        hbw = get_polygons_ptr(*hbw_objects);
        sample_sphere_from_sph(rdatax, rdatay, rdataz, hbw,
                               n_triangles, reparam, bw);

        if (DUMP_FILES) 
                output_graphics_any_format("hbw.obj", ASCII_FORMAT, 1,
                                           hbw_objects, NULL);

        if (DEBUG) fprintf(stderr,"butterworth_filter...\n");
        butterworth_filter(bw, lim, rcx, lrcx);
        butterworth_filter(bw, lim, rcy, lrcy);
        butterworth_filter(bw, lim, rcz, lrcz);
        butterworth_filter(bw, lim, icx, licx);
        butterworth_filter(bw, lim, icy, licy);
        butterworth_filter(bw, lim, icz, licz);
    
        if (DEBUG) fprintf(stderr,"get_realdata_from_sph_coeffs (lbw)...\n");
        get_realdata_from_sph_coeffs(rdatax, bw, DATAFORMAT, lrcx, licx);
        get_realdata_from_sph_coeffs(rdatay, bw, DATAFORMAT, lrcy, licy);
        get_realdata_from_sph_coeffs(rdataz, bw, DATAFORMAT, lrcz, licz);

        lbw_objects = (object_struct **) malloc(sizeof(object_struct *));
        *lbw_objects = create_object(POLYGONS);
        lbw = get_polygons_ptr(*lbw_objects);
        if (DEBUG) fprintf(stderr,"sample_sphere_from_sph (lbw)...\n");
        sample_sphere_from_sph(rdatax, rdatay, rdataz,
                               lbw, n_triangles, reparam, bw);

        if (DUMP_FILES) 
                output_graphics_any_format("lbw.obj", ASCII_FORMAT, 1,
                                           lbw_objects, NULL);

        free(rcx);  free(rcy);  free(rcz);
        free(icx);  free(icy);  free(icz);
        free(lrcx); free(lrcy); free(lrcz);
        free(licx); free(licy); free(licz);
        free(rdatax); free(rdatay); free(rdataz);

        /* make post correction */
        sph_postcorrect(surface, sphere, defects, polydefects, n_defects, holes,
                        hbw, lbw, do_surface_deform);

        /* make refinement to guarantee small sized vertices */ 
        if (max_refine_length > 0.0) {
                SET_ARRAY_SIZE( length_points, 0, hbw->n_points, DEFAULT_CHUNK_SIZE );
                for_less( i, 0, hbw->n_points )
                        length_points[i] = hbw->points[i];

                do {
                        n_done = refine_mesh( &length_points, hbw, max_refine_length,
                              &refined );

                        delete_polygons( hbw );
                        *hbw = refined;
                } while( n_done > 0 );

                print( "Resampled into %d polygons.\n", hbw->n_items );
        
        }
        
        delete_object_list(1, lbw_objects);
        free(defects);
        free(polydefects);
        free(holes);
        free(HD_handle);
        free(HD_hole);
        free(curv);

        return hbw_objects;
}
