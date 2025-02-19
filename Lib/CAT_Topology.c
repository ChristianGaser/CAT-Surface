/* Rachel Yotter - rachel.yotter@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id$
 */

#include <bicpl.h>

#include "CAT_SPH.h"
#include "CAT_Surf.h"
#include "CAT_Intersect.h"
#include "CAT_SurfaceIO.h"
#include "CAT_Curvature.h"
#include "CAT_Defect.h"
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
                maxi = MIN(dims[0], growing[growingCur     ] + 2);
                maxj = MIN(dims[1], growing[growingCur + 1] + 2);
                maxk = MIN(dims[2], growing[growingCur + 2] + 2);
                
                mini = MAX(0, growing[growingCur        ] - 1);
                minj = MAX(0, growing[growingCur + 1] - 1);
                mink = MAX(0, growing[growingCur + 2] - 1);
                
                for (tk = mink; tk < maxk; ++tk) for (tj = minj; tj < maxj; ++tj) for (ti = mini; ti < maxi; ++ti)
                {
                    ind1 = tk*(dims[0]*dims[1]) + (tj*dims[0]) + ti;
                    
                    if (!flagUsed[ind1] && bw[ind1] > 0)
                    {
                        flagUsed[ind1] = 1;
                        growing[growingInd       ] = ti;
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

    if (DEBUG) printf("remap_defect...\n");
    remap_defect(sphere, defects, polydefects, remap, remap_defects,
           remap_polydefects);

    delete_object_list(1, objects);
}

void
sph_postcorrect(polygons_struct *surface, polygons_struct *sphere, int *defects, int *polydefects, 
        int n_defects, int *holes, polygons_struct *hbw, polygons_struct *lbw)
{
    object_struct *surface_object;
    double *sharpness;
    int *n_neighbours, **neighbours;
    int p, d, val, poly, *flag, smooth_iters;
    int *hbw_defects, *hbw_polydefects, *hbw_holes;

    /* remap the defects onto the spherical harmonic reconstruction */
    if (DEBUG) printf("create_polygon_point_neighbours...\n");
    create_polygon_point_neighbours(hbw, TRUE, &n_neighbours,
                    &neighbours, NULL, NULL);

    if (DEBUG) printf("resample_defects_sph...\n");
    hbw_defects =   (int *) malloc(sizeof(int) * hbw->n_points);
    hbw_polydefects = (int *) malloc(sizeof(int) * hbw->n_items);
    resample_defects_sph(sphere, defects, polydefects,
               hbw_defects, hbw_polydefects, hbw->n_items);
    if (DEBUG) printf("expand_defects...\n");
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

    if (DEBUG) printf("compute_local_sharpness...\n");
    sharpness = (double *) malloc(sizeof(double) * hbw->n_points);
    compute_local_sharpness(hbw, n_neighbours, neighbours, sharpness);

    if (DEBUG) printf("patch_points...\n");
    flag = (int *) malloc(sizeof(int) * hbw->n_points);
    memset(flag, 0, sizeof(int) * hbw->n_points);
    for (p = 0; p < hbw->n_points; p++) {
        if (sharpness[p] > 60 && hbw_defects[p] != 0) {
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
    if (DEBUG) printf("find_selfintersections...\n");
    n_defects = find_selfintersections(hbw, hbw_defects, hbw_polydefects);
    n_defects = join_intersections(hbw, hbw_defects, hbw_polydefects,
                    n_neighbours, neighbours);

    printf("%d self intersection(s) to repair\n", n_defects);

    /* patch self-intersections first */
    if (DEBUG) printf("patch_selfintersections...\n");
    n_defects = patch_selfintersections(hbw, lbw, hbw_defects,
                      hbw_polydefects, n_defects,
                      n_neighbours, neighbours);
    printf("Post-patch: %d self intersection(s) remaining\n", n_defects);

    /* smooth out remaining self-intersections */
    smooth_iters = 50;
    if (DEBUG) printf("smooth_selfintersections with %d iterations...\n", smooth_iters);
    n_defects = smooth_selfintersections(hbw, hbw_defects, hbw_polydefects,
                       n_defects, n_neighbours,
                       neighbours, smooth_iters);


    if (DUMP_FILES) { 
        output_values_any_format("hbw_modpts.txt", hbw->n_points, flag,
                 TYPE_INTEGER);
    }

    if (DEBUG) printf("compute_polygon_normals...\n");
    compute_polygon_normals(hbw);

    free(sharpness);
    free(flag);
    free(hbw_defects);
    free(hbw_polydefects);
    free(hbw_holes);
}

object_struct **
fix_topology_sph(polygons_struct *surface, polygons_struct *sphere, int n_triangles, int bw, int lim, 
    char *reparam_file, double max_refine_length, int force, double laplace_thresh)
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
    double *HD_hole, *HD_handle, *curv, *tmp_vol;
    double sum_HD_handle, sum_HD_hole;
    unsigned char *mask;
    File_formats format;

    defects   = (int *) malloc(sizeof(int) * surface->n_points);
    holes    = (int *) malloc(sizeof(int) * surface->n_points);
    polydefects = (int *) malloc(sizeof(int) * surface->n_items);
    HD_handle  = (double *) malloc(sizeof(double) * surface->n_points);
    HD_hole   = (double *) malloc(sizeof(double) * surface->n_points);
    curv    = (double *) malloc(sizeof(double) * surface->n_points);
    hbw_objects = (object_struct **) malloc(sizeof(object_struct *));

    /* find defects in original uncorrected surface */
    create_polygon_point_neighbours(sphere, TRUE, &n_neighbours,
                    &neighbours, NULL, NULL);

    if (DEBUG) printf("find_topological_defects...\n");
    n_defects = find_topological_defects(surface, sphere, defects,
                       n_neighbours, neighbours);
    printf("%d topological defect(s)\n", n_defects);

    if (DEBUG) printf("update_polydefects...\n");
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
    rcx   = (double *) malloc(sizeof(double) * bw2);
    rcy   = (double *) malloc(sizeof(double) * bw2);
    rcz   = (double *) malloc(sizeof(double) * bw2);
    icx   = (double *) malloc(sizeof(double) * bw2);
    icy   = (double *) malloc(sizeof(double) * bw2);
    icz   = (double *) malloc(sizeof(double) * bw2);
    lrcx  = (double *) malloc(sizeof(double) * bw2);
    lrcy  = (double *) malloc(sizeof(double) * bw2);
    lrcz  = (double *) malloc(sizeof(double) * bw2);
    licx  = (double *) malloc(sizeof(double) * bw2);
    licy  = (double *) malloc(sizeof(double) * bw2);
    licz  = (double *) malloc(sizeof(double) * bw2);

    memset(holes, 0, sizeof(int) * surface->n_points);
    if (force ) {
      /* label defects according to force parameter to either use 
        holes or handles by default */
      for (p = 0; p < surface->n_points; p++) {
          if (defects[p] > 0) 
              holes[p] = force; 
          else  holes[p] = 0;
      }
    }

    if (DUMP_FILES) {
        output_values_any_format("defects.txt", surface->n_points,
                     holes, TYPE_INTEGER);
    }

    if (DEBUG) printf("get_equally_sampled_coords_holes...\n");
    get_equally_sampled_coords_holes(surface, sphere, defects, n_defects,
                     holes, bw, rdatax, rdatay, rdataz, force);

    if (DEBUG) printf("laplace2d...\n");
    if (laplace_thresh > 0.000001 && laplace_thresh < 0.5) {
    
        /* obtain a filtered map */
        tmp_vol = laplace2d(rdataz, NULL, 2*bw, 2*bw, 0.49);
        
        /* estimate mask by using the difference between orignal and filtered map */
        for (i = 0; i < 4*bw2; i++) tmp_vol[i] = (rdataz[i] - tmp_vol[i]);
        mask = threshold_image(tmp_vol, 2*bw, 2*bw, 1);
        
        /* filter inside mask which now defines areas with local deviations that are not continuous */
        rdataz = laplace2d(rdataz, mask, 2*bw, 2*bw, laplace_thresh);

        /* do the same for y-coordinates */
        tmp_vol = laplace2d(rdatay, NULL, 2*bw, 2*bw, 0.49);
        for (i = 0; i < 4*bw2; i++) tmp_vol[i] -= rdatay[i];
        mask = threshold_image(tmp_vol, 2*bw, 2*bw, 1);
        rdatay = laplace2d(rdatay, mask, 2*bw, 2*bw, laplace_thresh);

        /* and finally run filtering twice for x-coordinates because we have to consider both
          hemisphere where the x-coordinates are flipped */
        tmp_vol = laplace2d(rdatax, NULL, 2*bw, 2*bw, 0.49);
        for (i = 0; i < 4*bw2; i++) tmp_vol[i] -= rdatax[i];
        mask = threshold_image(tmp_vol, 2*bw, 2*bw, 1);
        rdatax = laplace2d(rdatax, mask, 2*bw, 2*bw, laplace_thresh);

        tmp_vol = laplace2d(rdatax, NULL, 2*bw, 2*bw, 0.49);
        for (i = 0; i < 4*bw2; i++) tmp_vol[i] = (rdatax[i] - tmp_vol[i]);
        mask = threshold_image(tmp_vol, 2*bw, 2*bw, 1);
        rdatax = laplace2d(rdatax, mask, 2*bw, 2*bw, laplace_thresh);

        free(tmp_vol);
        free(mask);

    }

    if (DEBUG) printf("get_sph_coeffs_of_realdata...\n");
    get_sph_coeffs_of_realdata(rdatax, bw, DATAFORMAT, rcx, icx);
    get_sph_coeffs_of_realdata(rdatay, bw, DATAFORMAT, rcy, icy);
    get_sph_coeffs_of_realdata(rdataz, bw, DATAFORMAT, rcz, icz);
    
    if (DEBUG) printf("sample_sphere_from_sph...\n");
    *hbw_objects = create_object(POLYGONS);
    hbw = get_polygons_ptr(*hbw_objects);
    sample_sphere_from_sph(rdatax, rdatay, rdataz, hbw,
                n_triangles, reparam, bw);

    if (DUMP_FILES) 
        output_graphics_any_format("hbw.obj", ASCII_FORMAT, 1,
                      hbw_objects, NULL);

    if (DEBUG) printf("butterworth_filter...\n");
    butterworth_filter(bw, lim, rcx, lrcx);
    butterworth_filter(bw, lim, rcy, lrcy);
    butterworth_filter(bw, lim, rcz, lrcz);
    butterworth_filter(bw, lim, icx, licx);
    butterworth_filter(bw, lim, icy, licy);
    butterworth_filter(bw, lim, icz, licz);
  
    if (DEBUG) printf("get_realdata_from_sph_coeffs (lbw)...\n");
    get_realdata_from_sph_coeffs(rdatax, bw, DATAFORMAT, lrcx, licx);
    get_realdata_from_sph_coeffs(rdatay, bw, DATAFORMAT, lrcy, licy);
    get_realdata_from_sph_coeffs(rdataz, bw, DATAFORMAT, lrcz, licz);

    lbw_objects = (object_struct **) malloc(sizeof(object_struct *));
    *lbw_objects = create_object(POLYGONS);
    lbw = get_polygons_ptr(*lbw_objects);
    
    if (DEBUG) printf("sample_sphere_from_sph (lbw)...\n");
    sample_sphere_from_sph(rdatax, rdatay, rdataz,
                lbw, n_triangles, reparam, bw);

    if (DUMP_FILES) 
        output_graphics_any_format("lbw.obj", ASCII_FORMAT, 1,
                      lbw_objects, NULL);

    free(rcx); free(rcy); free(rcz);
    free(icx); free(icy); free(icz);
    free(lrcx); free(lrcy); free(lrcz);
    free(licx); free(licy); free(licz);
    free(rdatax); free(rdatay); free(rdataz);

    /* make post correction */
    if (DEBUG) printf("sph_postcorrect...\n");
    sph_postcorrect(surface, sphere, defects, polydefects, n_defects, holes,
            hbw, lbw);

    /* make refinement to guarantee small sized vertices */ 
    if (DEBUG) printf("refine_mesh...\n");
    if (max_refine_length > 0.0) {
        SET_ARRAY_SIZE( length_points, 0, hbw->n_points, DEFAULT_CHUNK_SIZE );
        for_less( i, 0, hbw->n_points )
            length_points[i] = hbw->points[i];

        do {
            n_done = refine_mesh( &length_points, hbw, max_refine_length,
               &refined, 0.0 );

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