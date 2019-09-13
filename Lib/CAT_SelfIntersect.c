/*
   Copyright Alan C. Evans
   Professor of Neurology
   McGill University
*/

#include  <fit_3d.h>

#define  MAX_THRESHOLD_N_POLYGONS   200
#define  THRESHOLD_N_POLYGONS   60

#ifdef PRINT_DIST
#define PRINT_DIST
static Real sum_dist;
#endif

#define  A0  (1 << 0)
#define  A1  (1 << 1)
#define  A2  (1 << 2)
#define  B0  (1 << 3)
#define  B1  (1 << 4)
#define  B2  (1 << 5)

private  void    recursive_find_close_pairs( 
    Real               min_distance,
    Real               search_distance,
    Real               parameters[],
    poly_info_struct   poly_info[],
    int                **list_of_polys,
    unsigned int       **poly_classes,
    int                *n_alloced,
    int                offset_index,
    int                n_polygons,
    int                next_offset_index,
    Real               line_dir[],
    int                *n_pairs_ptr,
    int                *n_pairs_alloc,
    int                *p1s_ptr[],
    int                *p2s_ptr[],
    unsigned char      *cases_ptr[],
    float              *min_line_dists_ptr[],
    Smallest_int       used_flags[],
    Real               *closest_distance )
{
    int               parms1[3], parms2[3];
    Real              mid_plane, min_distance_sq, search_distance_sq;
    Real              movement, movement_sq, dx, dy, dz, d, t_min, dist;
    poly_info_struct  *info_ptr;
    poly_info_struct  *static_compare[MAX_THRESHOLD_N_POLYGONS];
    poly_info_struct  **compare;
    int               static_compare_index[MAX_THRESHOLD_N_POLYGONS];
    int               *compare_index;
    int               n1, n2;
    float             mid_plane_plus, mid_plane_minus;
    Real              min_range[N_DIMENSIONS], max_range[N_DIMENSIONS];
    float             x1_min, x1_max, y1_min, y1_max, z1_min, z1_max;
    float             fl_search_distance;
    float             x1_low, x1_high, y1_low, y1_high, z1_low, z1_high;
    int               split_axis, dim, *list, poly1, p1, p2, poly;
    unsigned int      *classes, *src_classes, class_bit;
    unsigned int      *new_classes1, *new_classes2;
    int               i, i1, i2, diff;
    int               *src_list, *new_list1, *new_list2;
    int               new_n_polygons1, new_n_polygons2;
    int               n11, n12, n21, n22;
    int               parm_p1, parm_p2, parm_n11, parm_n12, parm_n21, parm_n22;
    BOOLEAN           on_left, on_right;
    unsigned int      cl;
    unsigned char     which_case;
    int               count[17], c1, c2, start;
    int               static_list[MAX_THRESHOLD_N_POLYGONS];
    int               *sorted_list;
    int               n_to_compare;
    static int        threshold_n_polygons = THRESHOLD_N_POLYGONS;
    static BOOLEAN    first = TRUE;

    if( first ) {
        first = FALSE;
        if( getenv( "THRESHOLD_N_POLYGONS1" ) == NULL ||
            sscanf( getenv("THRESHOLD_N_POLYGONS1" ), "%d",
                    &threshold_n_polygons ) != 1  ||
            threshold_n_polygons > MAX_THRESHOLD_N_POLYGONS )
            threshold_n_polygons = THRESHOLD_N_POLYGONS;
    }

    if( n_polygons == 0 )
        return;

    if( n_polygons > threshold_n_polygons ) {

      // Compute the true min/max bounding box of this list of polygons.
      list = &(*list_of_polys)[offset_index];
      for( dim = 0; dim < N_DIMENSIONS; dim++ ) {
        min_range[dim] = poly_info[list[0]].low_limits[dim];
        max_range[dim] = poly_info[list[0]].high_limits[dim];
      }

      for_less( i, 1, n_polygons ) {
        poly = list[i];
        for( dim = 0; dim < N_DIMENSIONS; dim++ ) {
          if( poly_info[poly].low_limits[dim] < min_range[dim] ) {
            min_range[dim] = poly_info[poly].low_limits[dim];
          } else {
            if( poly_info[poly].high_limits[dim] > max_range[dim] ) {
              max_range[dim] = poly_info[poly].high_limits[dim];
            }
          }
        }
      }
   
      split_axis = 0;

      for_less( dim, 1, N_DIMENSIONS ) {
        if( max_range[dim] - min_range[dim] >
            max_range[split_axis] - min_range[split_axis] ) {
          split_axis = dim;
        }
      }

      if( max_range[split_axis] - min_range[split_axis] >= 2.0*search_distance ) {
        if( next_offset_index + 2 * n_polygons > *n_alloced ) {
            *n_alloced = next_offset_index + 2 * n_polygons;
            REALLOC( *list_of_polys, *n_alloced );
            REALLOC( *poly_classes, *n_alloced );
        }

        // Compute the mid_plane by averaging the centroids
        // of the bounding boxes of the triangles. This gives 
        // a much better balance for the number of polygons
        // per branch.
        mid_plane = 0.0;
        for_less( i, 0, n_polygons ) {
          poly = list[i];
          mid_plane += poly_info[poly].low_limits[split_axis] +
                       poly_info[poly].high_limits[split_axis];
        }
        mid_plane = 0.5 * mid_plane / ( (float)n_polygons );
        mid_plane_plus = (float) (mid_plane + search_distance / 2.0);
        mid_plane_minus = (float) (mid_plane - search_distance / 2.0);

        /*--- do everything to left */

        src_list = &(*list_of_polys)[offset_index];
        src_classes = &(*poly_classes)[offset_index];

        new_list1 = &(*list_of_polys)[next_offset_index];
        new_classes1 = &(*poly_classes)[next_offset_index];
        new_n_polygons1 = 0;

        new_list2 = &(*list_of_polys)[next_offset_index+n_polygons];
        new_classes2 = &(*poly_classes)[next_offset_index+n_polygons];
        new_n_polygons2 = 0;

        class_bit = (1u << split_axis);

        for_less( i, 0, n_polygons ) {
            poly = src_list[i];
            on_left = ( poly_info[poly].low_limits[split_axis] <= mid_plane_plus );
            on_right = ( poly_info[poly].high_limits[split_axis] >= mid_plane_minus );
            if( on_left ) {
                new_list1[new_n_polygons1] = poly;
                new_classes1[new_n_polygons1] = src_classes[i];
                ++new_n_polygons1;
            }
            if( on_right ) {
                new_list2[new_n_polygons2] = poly;
                if( on_left )
                    new_classes2[new_n_polygons2] = (src_classes[i] |class_bit);
                else
                    new_classes2[new_n_polygons2] = src_classes[i];
                ++new_n_polygons2;
            }
        }

        diff = n_polygons - new_n_polygons1;
        for_less( i, 0, new_n_polygons2 ) {
            new_list2[i-diff] = new_list2[i];
            new_classes2[i-diff] = new_classes2[i];
        }

        // Avoid creating a branch with too few triangles or with
        // nearly the same number of triangles as its parent. We
        // split if in between 15% and 85% of n_polygons (that's
        // the corresponding value for 0.1275). Work in float to
        // avoid int overflow on large meshes.
        float r1 = (float)new_n_polygons1 / (float)n_polygons;
        float r2 = (float)new_n_polygons2 / (float)n_polygons;

        if( r1*(1.0-r1) > 0.1275 && r2*(1.0-r2) > 0.1275 ) {

            recursive_find_close_pairs( min_distance, search_distance,
                        parameters, poly_info, list_of_polys, poly_classes,
                        n_alloced, next_offset_index, new_n_polygons1,
                        next_offset_index + new_n_polygons1 + new_n_polygons2,
                        line_dir, n_pairs_ptr, n_pairs_alloc, p1s_ptr, p2s_ptr,
                        cases_ptr, min_line_dists_ptr, used_flags,
                        closest_distance );

            recursive_find_close_pairs( min_distance, search_distance,
                        parameters, poly_info, list_of_polys, poly_classes,
                        n_alloced, next_offset_index + new_n_polygons1,
                        new_n_polygons2,
                        next_offset_index + new_n_polygons1 + new_n_polygons2,
                        line_dir, n_pairs_ptr, n_pairs_alloc, p1s_ptr, p2s_ptr, 
                        cases_ptr, min_line_dists_ptr, used_flags,
                        closest_distance );

            return;
        }
      }
    }

    min_distance_sq = min_distance * min_distance;
    search_distance_sq = search_distance * search_distance;
    fl_search_distance = (float) search_distance;
    list = &(*list_of_polys)[offset_index];
    classes = &(*poly_classes)[offset_index];

    if( n_polygons > MAX_THRESHOLD_N_POLYGONS ) {
        ALLOC( sorted_list, n_polygons );
        ALLOC( compare, n_polygons );
        ALLOC( compare_index, n_polygons );
    } else {
        sorted_list = static_list;
        compare = static_compare;
        compare_index = static_compare_index;
    }
        
    for_less( i1, 0, 16 )
        count[i1] = 0;

    for_less( i1, 0, n_polygons )
        ++count[classes[i1]];

    for_less( i1, 1, 16 )
        count[i1] += count[i1-1];

    for_less( i1, 0, n_polygons ) {
        cl = classes[i1];
        --count[cl];
        sorted_list[count[cl]] = list[i1];
    }

    count[16] = n_polygons;

    for_less( c1, 0, 16 ) {
        if( count[c1] == count[c1+1] )
            continue;

        n_to_compare = 0;

        for_less( c2, c1, 16 ) {
            if( (c1 & c2) != 0 || count[c2] == count[c2+1] )
                continue;

            for_less( i2, count[c2], count[c2+1] ) {
                compare_index[n_to_compare] = sorted_list[i2];
                compare[n_to_compare] = &poly_info[sorted_list[i2]];
                ++n_to_compare;
            }
        }

        for_less( i1, count[c1], count[c1+1] ) {
            poly1 = sorted_list[i1];
            p1 = poly_info[poly1].p1;
            n11 = poly_info[poly1].n11;
            n12 = poly_info[poly1].n12;

            used_flags[p1] = TRUE;
            used_flags[n11] = TRUE;
            used_flags[n12] = TRUE;

            parm_p1 = IJ( p1, 0, 3 );
            parm_n11 = IJ( n11, 0, 3 );
            parm_n12 = IJ( n12, 0, 3 );

            x1_min = poly_info[poly1].low_limits[X];
            x1_max = poly_info[poly1].high_limits[X];
            y1_min = poly_info[poly1].low_limits[Y];
            y1_max = poly_info[poly1].high_limits[Y];
            z1_min = poly_info[poly1].low_limits[Z];
            z1_max = poly_info[poly1].high_limits[Z];

            x1_low = x1_min - fl_search_distance;
            x1_high = x1_max + fl_search_distance;
            y1_low = y1_min - fl_search_distance;
            y1_high = y1_max + fl_search_distance;
            z1_low = z1_min - fl_search_distance;
            z1_high = z1_max + fl_search_distance;

            if( c1 == 0 )
                start = i1 + 1; 
            else
                start = 0; 

            for_less( i2, start, n_to_compare ) {
                info_ptr = compare[i2];

                p2 = info_ptr->p1;
                n21 = info_ptr->n11;
                n22 = info_ptr->n12;
                if( info_ptr->low_limits[X] > x1_high ||
                    info_ptr->high_limits[X] < x1_low ||
                    info_ptr->low_limits[Y] > y1_high ||
                    info_ptr->high_limits[Y] < y1_low ||
                    info_ptr->low_limits[Z] > z1_high ||
                    info_ptr->high_limits[Z] < z1_low ||
                    used_flags[p2] || used_flags[n21] || used_flags[n22] )
                {
                    continue;
                }

                parm_p2 = 3 * p2;
                parm_n21 = 3 * n21;
                parm_n22 = 3 * n22;
                if( sq_triangle_triangle_dist_estimate(
                                                &parameters[parm_p1],
                                                &parameters[parm_n11],
                                                &parameters[parm_n12],
                                                &parameters[parm_p2],
                                                &parameters[parm_n21],
                                                &parameters[parm_n22],
                                                search_distance_sq ) ) {
                    continue;
                }

                which_case = 0;
                dist = sq_triangle_triangle_dist( &parameters[parm_p1],
                                                  &parameters[parm_n11],
                                                  &parameters[parm_n12],
                                                  &parameters[parm_p2],
                                                  &parameters[parm_n21],
                                                  &parameters[parm_n22],
                                                  &which_case );

                if( dist >= search_distance_sq ) {
                    continue;
                }

#ifdef PRINT_DIST
sum_dist += dist;
#endif

                if( *closest_distance < 0.0 || dist < *closest_distance )
                    *closest_distance = dist;

                if( (*n_pairs_ptr) >= (*n_pairs_alloc) ) {
                  // Re-allocate more space and copy. Add 10% more.
                  *n_pairs_alloc = (int)( 1.1 * (float)(*n_pairs_alloc) );
                  REALLOC( *p1s_ptr, *n_pairs_alloc );
                  REALLOC( *p2s_ptr, *n_pairs_alloc );
                  REALLOC( *cases_ptr, *n_pairs_alloc );
                  REALLOC( *min_line_dists_ptr, *n_pairs_alloc );
                }

                (*min_line_dists_ptr)[*n_pairs_ptr] = (float)dist;
                (*p1s_ptr)[*n_pairs_ptr] = poly1;
                (*p2s_ptr)[*n_pairs_ptr] = compare_index[i2];
                (*cases_ptr)[*n_pairs_ptr] = which_case;
                ++(*n_pairs_ptr);
            }

            used_flags[p1] = FALSE;
            used_flags[n11] = FALSE;
            used_flags[n12] = FALSE;
        }
    }

    if( n_polygons > MAX_THRESHOLD_N_POLYGONS ) {
        FREE( sorted_list );
        FREE( compare );
        FREE( compare_index );
    }
}


// Find pairs of triangular faces that are within the distance
// min_distance from one another.
// In:
//   min_distance : distance threshold
//   n_points     : number of nodes
//   n_triangles  : number of triangles
//   triangles    : connectivity of the triangles
//   parameters   : (x,y,z) coords of the nodes
//   max_movement :
//   line_dir     : derivative (gradient) at the nodes
// Out:
//   n_pairs_ptr : number of pairs
//   p1s_ptr     : first triangle of each pair
//   p2s_ptr     : second triangle of each pair
//   cases_ptr   : some flag for each pair (intersection pattern)
//   dist_ptr    : true min distance between two triangles
//   min_line_dists_ptr : some distance between the two triangles
//                        of each pair
//   closest_distance : 
// Return:
//   none
//
private  void   find_self_intersect_candidates(
    Real               min_distance,
    int                *n_pairs_ptr,
    int                *n_pairs_alloc,
    int                *p1s_ptr[],
    int                *p2s_ptr[],
    unsigned char      *cases_ptr[],
    float              *min_line_dists_ptr[],
    int                n_points,
    int                n_triangles,
    int                triangles[],
    int              * n_neighbours,
    int              * neighbours[],
    Real               parameters[],
    Real               max_movement,
    Real               line_dir[],
    Real               *closest_distance )
{
    poly_info_struct  *poly_info;
    int               poly1, n_polygons;
    int               p1, dim, p, n;
    int               n11, n12;
    int               i, j;
    int               *list_of_polys, n_alloced;
    unsigned int      *poly_classes;
    Real              search_distance;
    Smallest_int      *used_flags;

#ifdef PRINT_DIST
sum_dist = 0.0;
#endif

    search_distance = min_distance + max_movement;

    // Number of triangles.
    n_polygons = n_triangles;

    // Store the triangles in poly_info in the same order
    // as they were read. Also, find the min/max bounding
    // box for each triangle.
    ALLOC( poly_info, n_polygons );
    for( p = 0; p < n_polygons; p++ ) {
      // n = 0;
      // if( triangles[3*p+1] < triangles[3*p+n] ) n = 1;
      // if( triangles[3*p+2] < triangles[3*p+n] ) n = 2;
      // p1 = triangles[3*p+n];
      // n11 = triangles[3*p+(n+1)%3];
      // n12 = triangles[3*p+(n+2)%3];
      // poly_info[p].p1 = p1;           // want p1 < n11 && p1 < n12
      //  poly_info[p].n11 = n11;
      // poly_info[p].n12 = n12;
      //
      // make sure that triangles are ordered in the same way as poly_info
      // (this is such that which_case can be applied in the evaluation using 
      // the same ordering with triangles as it is defined with poly_info).
      poly_info[p].p1 = triangles[3*p];
      poly_info[p].n11 = triangles[3*p+1];
      poly_info[p].n12 = triangles[3*p+2];

      p1 = poly_info[p].p1;
      n11 = poly_info[p].n11;
      n12 = poly_info[p].n12;

      for_less( dim, 0, N_DIMENSIONS ) {
        poly_info[p].low_limits[dim] = (float) MIN3( parameters[IJ(p1,dim,3)],
                                                     parameters[IJ(n11,dim,3)],
                                                     parameters[IJ(n12,dim,3)] );
        poly_info[p].high_limits[dim] = (float) MAX3( parameters[IJ(p1,dim,3)],
                                                      parameters[IJ(n11,dim,3)],
                                                      parameters[IJ(n12,dim,3)] );
      }
    }

    // Set-up arrays for the partitionning of the domain.
    // It seems that 4*n_polygons is big enough to hold 
    // the temporary lists of polygons as the tree is
    // being constructed.

    n_alloced = 4 * n_polygons;
    ALLOC( list_of_polys, n_alloced );
    ALLOC( poly_classes, n_alloced );

    for_less( poly1, 0, n_polygons ) {
        list_of_polys[poly1] = poly1;
        p1 = poly_info[poly1].p1;
        n11 = poly_info[poly1].n11;
        n12 = poly_info[poly1].n12;

        if( line_dir[IJ(p1,0,3)] == 0.0 &&
            line_dir[IJ(p1,1,3)] == 0.0 &&
            line_dir[IJ(p1,2,3)] == 0.0 &&
            line_dir[IJ(n11,0,3)] == 0.0 &&
            line_dir[IJ(n11,1,3)] == 0.0 &&
            line_dir[IJ(n11,2,3)] == 0.0 &&
            line_dir[IJ(n12,0,3)] == 0.0 &&
            line_dir[IJ(n12,1,3)] == 0.0 &&
            line_dir[IJ(n12,2,3)] == 0.0 ) {
            poly_classes[poly1] = 8;   // derivative = 0
        } else {
            poly_classes[poly1] = 0;   // non-zero derivative/gradient at node
        }
    }

    ALLOC( used_flags, n_points );
    for_less( p1, 0, n_points )
        used_flags[p1] = FALSE;

    *closest_distance = -1.0;
    *n_pairs_ptr = 0;

    recursive_find_close_pairs( min_distance, search_distance,
                        parameters, poly_info, &list_of_polys, &poly_classes,
                        &n_alloced, 0, n_polygons, n_polygons,
                        line_dir, n_pairs_ptr, n_pairs_alloc,
                        p1s_ptr, p2s_ptr, cases_ptr,
                        min_line_dists_ptr, used_flags,
                        closest_distance );
    FREE( used_flags );

    if( n_alloced > 0 ) {
        FREE( list_of_polys );
        FREE( poly_classes );
    }

    if( *closest_distance > 0.0 )
        *closest_distance = sqrt( *closest_distance );

#if 0

    // Note by Claude: I have no clue why I once added this
    // piece of code below, since it ignores nearby points 
    // and leads to self-intersections in CLASP. Maybe the
    // code below can give a better count??? Can't remember.

    // Ignore cases for which the vertices involved are 
    // directly connected via another vertex.

    int new_n_pairs_ptr = 0;
    for( n = 0; n < *n_pairs_ptr; n++ ) {
      int poly1 = (*p1s_ptr)[n];
      int poly2 = (*p2s_ptr)[n];
      int which_case = (*cases_ptr)[n];
      
      int a0 = poly_info[poly1].p1;
      int a1 = poly_info[poly1].n11;
      int a2 = poly_info[poly1].n12;
      int b0 = poly_info[poly2].p1;
      int b1 = poly_info[poly2].n11;
      int b2 = poly_info[poly2].n12;

      int nn[6];
      nn[0] = ( which_case & A0 ) ? ( a0 ) : ( -1 );
      nn[1] = ( which_case & A1 ) ? ( a1 ) : ( -1 );
      nn[2] = ( which_case & A2 ) ? ( a2 ) : ( -1 );
      nn[3] = ( which_case & B0 ) ? ( b0 ) : ( -1 );
      nn[4] = ( which_case & B1 ) ? ( b1 ) : ( -1 );
      nn[5] = ( which_case & B2 ) ? ( b2 ) : ( -1 );
      
      int count = 6;
      for( i = 0; i < 6; i++ ) {
        if( nn[i] == -1 ) count--;
      }

      // check if any of the a's shares a vertex neighbour with
      // any of the b's. Note that the a's are connected to each
      // other and so are the b's.

      int common = 0, ki, kj;
      if( count == 3 || count == 4 ) {
        for( i = 0; i < 3; i++ ) {
          if( nn[i] == -1 ) continue;
          for( j = 3; j < 6; j++ ) {
            if( nn[j] == -1 ) continue;
            for( kj = 0; kj < n_neighbours[nn[j]]; kj++ ) {
              for( ki = 0; ki < n_neighbours[nn[i]]; ki++ ) {
                if( neighbours[nn[i]][ki] == neighbours[nn[j]][kj] ) common = 1;
              }
            }
          } 
          break;
        }
      }
      if( common == 0 ) {
        (*p1s_ptr)[new_n_pairs_ptr] = (*p1s_ptr)[n];
        (*p2s_ptr)[new_n_pairs_ptr] = (*p2s_ptr)[n];
        (*cases_ptr)[new_n_pairs_ptr] = (*cases_ptr)[n];
        (*min_line_dists_ptr)[new_n_pairs_ptr] = (*min_line_dists_ptr)[n];
        new_n_pairs_ptr++;
      }
    }
    *n_pairs_ptr = new_n_pairs_ptr;
#endif

    FREE( poly_info );

#ifdef PRINT_DIST
print( "Sum dist: %.15g\n", sum_dist );
#endif
}

// 
// Loop over all surfaces of the dataset to find, for each
// surface, the pairs of triangular faces that are within the 
// distance max_movement from one another.
// 
// In:
//   deform       : 
//   start_parameter : 
//   parameters   : (x,y,z) coords of the nodes
//   max_movement :
//   line_dir     : derivative (gradient) at the nodes
//   si_lookup    : table of surface intersection data
// Out:
//   si_lookup    : surface intersection data found
// Return:
//   closest_distance from all surfaces
//

public  Real  recompute_self_intersects(
    Deform_struct                 *deform,
    int                           start_parameter[],
    Real                          parameters[],
    Real                          max_movement,
    Real                          line_dir[],
    self_intersect_lookup_struct  **si_lookup ) {

    int                    i, w, surface;
    Real                   *this_parms, *this_line, closest_distance;
    self_intersect_struct  *self;
    Real                   min_distance, close;

    closest_distance = -1.0;

    for_less( surface, 0, deform->n_surfaces ) {
        this_parms = &parameters[start_parameter[surface]];
        if( line_dir == NULL ) {
            this_line = NULL;
        } else {
            this_line = &line_dir[start_parameter[surface]];
        }

        for_less( i, 0, deform->surfaces[surface].n_self_intersects ) {
            self = &deform->surfaces[surface].self_intersects[i];

            min_distance = 0.0;
            for_less( w, 0, self->n_weights ) {
                min_distance = MAX( min_distance, self->min_distances[w] );
            }

            find_self_intersect_candidates(
                       min_distance,
                       &si_lookup[surface][i].n_pairs,
                       &si_lookup[surface][i].n_pairs_alloc,
                       &si_lookup[surface][i].p1s,
                       &si_lookup[surface][i].p2s,
                       &si_lookup[surface][i].cases,
                       &si_lookup[surface][i].min_line_dists,
                       deform->surfaces[surface].surface.n_points,
                       deform->surfaces[surface].surface.n_polygons,
                       deform->surfaces[surface].surface.triangles,
                       deform->surfaces[surface].surface.n_neighbours,
                       deform->surfaces[surface].surface.neighbours,
                       this_parms, max_movement, this_line,
                       &close );

            if( close >= 0.0 &&
                (closest_distance < 0.0 || close < closest_distance) )
                closest_distance = close;
        }
    }

    return( closest_distance );
}


public  void  delete_self_intersect_lookup(
    self_intersect_lookup_struct  *si_lookup )
{
    if( si_lookup->n_pairs > 0 ) {
        FREE( si_lookup->p1s );
        FREE( si_lookup->p2s );
        FREE( si_lookup->cases );
        FREE( si_lookup->min_line_dists );
        si_lookup->n_pairs = 0;
        si_lookup->n_pairs_alloc = 0;
    }
}

public  int  get_n_self_intersect_candidate(
    self_intersect_lookup_struct  *si_lookup )
{
    return( si_lookup->n_pairs );
}

// ---------------------------------------------------------------------
// Find the minimal distance between between the two triangles
// (p1,n11,n12) and (p2,n21,n22), if these two triangles are
// close enough to intersect given the current distance of
// movement.
//
public  BOOLEAN   test_self_intersect_candidate(
    int                           p1,
    int                           n11,
    int                           n12,
    int                           p2,
    int                           n21,
    int                           n22,
    unsigned char *               which_case,
    Real                          min_line_dists,
    Real                          dist_from_computed_self_intersect,
    Real                          parameters[],
    Smallest_int                  active_flags[],
    Real                          *dist_sq ) {

    if( dist_from_computed_self_intersect > 0.0 &&
        min_line_dists > dist_from_computed_self_intersect ) {
        return( FALSE );
    }

    if( active_flags ) {
      if( !active_flags[p1] && !active_flags[n11] && !active_flags[n12] &&
          !active_flags[p2] && !active_flags[n21] && !active_flags[n22] ) {
        return( FALSE );
      }
    }

    *dist_sq = sq_triangle_triangle_dist( &parameters[p1*3],
                                          &parameters[n11*3],
                                          &parameters[n12*3],
                                          &parameters[p2*3],
                                          &parameters[n21*3],
                                          &parameters[n22*3], which_case );

    return( TRUE );
}

public  Real   get_self_intersect_deriv_factor(
    BOOLEAN                       use_square_flag,
    Real                          min_distance,
    Real                          weight,
    Real                          dist ) {

    Real  diff, factor;

    if( dist == 0.0 ) {
      factor = 0.0;
    } else {
      diff = min_distance - dist;
#ifdef USE_CUBE
#define USE_CUBE
      factor = -1.5 * diff * diff * weight / ( dist * min_distance );
#else
      // This 1/dist term probably comes from taking the derivative
      // of the distance squared. Claude
      //   d/dx(D^2) = 2*D*dD/dx = ...
      // so dD/dx = 1/(2*D)*(...)

      if( use_square_flag ) {
        factor = -diff * weight / dist;
      } else {
        if( diff < 0.0 ) {
            factor = 0.50 / dist;
        } else {
            factor = -0.50 / dist;
        }
      }
#endif
    }
    return( factor );
}

public  void   get_self_intersect_deriv(
    int                           p1,
    int                           n11,
    int                           n12,
    int                           p2,
    int                           n21,
    int                           n22,
    Real                          factor,
    Real                          deriv[],
    Real                          tri_tri_deriv[6][3] ) {

    if( factor == 0.0 ) return;

    deriv[IJ(p1,0,3)] += factor * tri_tri_deriv[0][0];
    deriv[IJ(p1,1,3)] += factor * tri_tri_deriv[0][1];
    deriv[IJ(p1,2,3)] += factor * tri_tri_deriv[0][2];

    deriv[IJ(n11,0,3)] += factor * tri_tri_deriv[1][0];
    deriv[IJ(n11,1,3)] += factor * tri_tri_deriv[1][1];
    deriv[IJ(n11,2,3)] += factor * tri_tri_deriv[1][2];

    deriv[IJ(n12,0,3)] += factor * tri_tri_deriv[2][0];
    deriv[IJ(n12,1,3)] += factor * tri_tri_deriv[2][1];
    deriv[IJ(n12,2,3)] += factor * tri_tri_deriv[2][2];

    deriv[IJ(p2,0,3)] += factor * tri_tri_deriv[3][0];
    deriv[IJ(p2,1,3)] += factor * tri_tri_deriv[3][1];
    deriv[IJ(p2,2,3)] += factor * tri_tri_deriv[3][2];

    deriv[IJ(n21,0,3)] += factor * tri_tri_deriv[4][0];
    deriv[IJ(n21,1,3)] += factor * tri_tri_deriv[4][1];
    deriv[IJ(n21,2,3)] += factor * tri_tri_deriv[4][2];

    deriv[IJ(n22,0,3)] += factor * tri_tri_deriv[5][0];
    deriv[IJ(n22,1,3)] += factor * tri_tri_deriv[5][1];
    deriv[IJ(n22,2,3)] += factor * tri_tri_deriv[5][2];
}
