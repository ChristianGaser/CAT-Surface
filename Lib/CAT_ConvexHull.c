/* Christian Gaser - christian.gaser@uni-jena.de
 * Department of Psychiatry
 * University of Jena
 *
 * Copyright Christian Gaser, University of Jena.
 * $Id: CAT_ConvexHull.c 267 2012-06-21 14:57:59Z gaser $
 *
 */

#include <bicpl.h>

#include "Cat_Surf.h"
#include "CAT_SurfaceIO.h"
#include "CAT_ConvexHull.h"
#include "CAT_Resample.h"

// -------------------------------------------------------------------
// Construct the edges around each node. The edges are sorted to
// make an ordered closed loop.
//
private int get_surface_neighbours( polygons_struct * surface,
                                    int * n_neighbours_return[],
                                    int ** neighbours_return[] ) {

  int    i, j, k, jj;
  int  * tri;
  int  * n_ngh;
  int ** ngh;
  int  * ngh_array;

  // Check if all polygons are triangles.

  if( 3 * surface->n_items != surface->end_indices[surface->n_items-1] ) {
    printf( "Surface must contain only triangular polygons.\n" );
    return ERROR;
  }

  // Check if the node numbering starts at 0 or 1.

  int min_idx, max_idx;

  min_idx = 100*surface->n_points;  // anything big
  max_idx = 0;                      // anything small

  for( i = 0; i < 3*surface->n_items; i++ ) {
    if( surface->indices[i] < min_idx ) min_idx = surface->indices[i];
    if( surface->indices[i] > max_idx ) max_idx = surface->indices[i];
  }

  // Shift numbering to start at zero, for array indexing. Note
  // that we don't care if surface->indices array is modified.

  if( min_idx != 0 ) {
    for( i = 0; i < 3*surface->n_items; i++ ) {
      surface->indices[i] -= min_idx;
    }
  }

  // Count number of triangles attached to each node.

  ALLOC( n_ngh, surface->n_points );
  ALLOC( ngh, surface->n_points );
  ALLOC( ngh_array, 3*surface->n_items );

  for( i = 0; i < surface->n_points; i++ ) {
    n_ngh[i] = 0;
  }

  for( i = 0; i < 3*surface->n_items; i++ ) {
    n_ngh[surface->indices[i]]++;
    ngh_array[i] = -1;
  }

  int max_ngh = 0;
  int sum_ngh = 0;
  for( i = 0; i < surface->n_points; i++ ) {
    ngh[i] = &(ngh_array[sum_ngh]);
    sum_ngh += n_ngh[i];
    max_ngh = MAX( max_ngh, n_ngh[i] );
  }

  // At first, store the indices of the triangles in the neighbours.
  for( i = 0; i < surface->n_items; i++ ) {
    for( j = 0; j < 3; j++ ) {
      jj = surface->indices[3*i+j];
      for( k = 0; k < n_ngh[jj]; k++ ) {
        if( ngh[jj][k] == -1 ) {
          ngh[jj][k] = i;
          break;
        }
      }
    }
  }

  // Now create a sort closed loop of the node neighbours.
  // This is needed by the parametric=0 FEM algorithm.
  //
  //         1 ----- 2
  //          /\   /\
  //         /  \ /  \
  //       0 ----P---- 3
  //         \  / \  /
  //          \/   \/
  //         5 ----- 4
  //

  int * tmp;
  ALLOC( tmp, 2*max_ngh );

  for( i = 0; i < surface->n_points; i++ ) {
    for( k = 0; k < n_ngh[i]; k++ ) {
      tri = &(surface->indices[3*ngh[i][k]]);
      for( j = 0; j < 3; j++ ) {
        if( tri[j] == i ) break;
      }
      tmp[2*k+0] = tri[(j+1)%3];
      tmp[2*k+1] = tri[(j+2)%3];
    }

    ngh[i][0] = tmp[0];
    ngh[i][1] = tmp[1];
    for( k = 2; k < n_ngh[i]; k++ ) {
      for( j = 1; j < n_ngh[i]; j++ ) {
        if( tmp[2*j] == ngh[i][k-1] || tmp[2*j+1] == ngh[i][k-1] ) {
          if( tmp[2*j] == ngh[i][k-1] ) {
            ngh[i][k] = tmp[2*j+1];
          } else {
            ngh[i][k] = tmp[2*j];
          }
          tmp[2*j] = -1;
          tmp[2*j+1] = -1;
          break;
        }
      }
    }
  }

  *n_neighbours_return = n_ngh;
  *neighbours_return = ngh;

  FREE( tmp );

  return OK;

}

private  int  get_points_of_region(
    polygons_struct *surface,
    Point   points[] ) {

    int     n_points;

    Point  * coords;             // coordinates
    Vector * normals;            // normal vectors
    int    * n_ngh = NULL;       // node neighbours (inverse connectivity)
    int   ** ngh = NULL;
    int      i, j;
    int      n_convex = 0;


    // Find the convex points for a surface.
    // Make a copy of the coordinates and the normals, since
    // delete_object_list will destroy them.

fprintf(stderr,".");
    n_points = surface->n_points;
    ALLOC( points, surface->n_points );
    ALLOC( normals, surface->n_points );
fprintf(stderr,".");
    for( i = 0; i < n_points; i++ ) {
      (points)[i].coords[0] = surface->points[i].coords[0];
      (points)[i].coords[1] = surface->points[i].coords[1];
      (points)[i].coords[2] = surface->points[i].coords[2];
      (normals)[i].coords[0] = surface->normals[i].coords[0];
      (normals)[i].coords[1] = surface->normals[i].coords[1];
      (normals)[i].coords[2] = surface->normals[i].coords[2];
    }

fprintf(stderr,".");
    get_surface_neighbours( surface, &n_ngh, &ngh );

fprintf(stderr,".");
    for( i = 0; i < n_points; i++ ) {
fprintf(stderr,":");
      double plane_const = coords[i].coords[0] * normals[i].coords[0] +
                           coords[i].coords[1] * normals[i].coords[1] +
                           coords[i].coords[2] * normals[i].coords[2];
      for( j = 0; j < n_ngh[i]; j++ ) {
fprintf(stderr,"%d/%d ",i,j);
        double check = coords[ngh[i][j]].coords[0] * normals[i].coords[0] +
                       coords[ngh[i][j]].coords[1] * normals[i].coords[1] +
                       coords[ngh[i][j]].coords[2] * normals[i].coords[2];
        if( check > plane_const ) break;
      }
      if( j == n_ngh[i] ) {
        ADD_ELEMENT_TO_ARRAY( points, n_convex, coords[i], DEFAULT_CHUNK_SIZE);
      }
    }
fprintf(stderr,".");
    if( coords ) FREE( coords );
    if( normals ) FREE( normals );
    if( n_ngh ) FREE( n_ngh );
    if( ngh ) {
      FREE( ngh[0] );   // this is ngh_array
      FREE( ngh );
    }
fprintf(stderr,".");
    n_points = n_convex;
 
    return( n_points );

}

private  Real  compute_clockwise_degrees( Real x, Real y )
{
    Real   degrees;

    if( x >= -TOLERANCE_DISTANCE && x <= TOLERANCE_DISTANCE )
    {
        if( y < -TOLERANCE_DISTANCE )
            return( 90.0 );
        else if( y > TOLERANCE_DISTANCE )
            return( 270.0 );
        else
            return( 0.0 );
    }
    else if( y >= -TOLERANCE_DISTANCE && y <= TOLERANCE_DISTANCE)
    {
        if( x > 0.0 )
            return( 0.0 );
        else
            return( 180.0 );
    }
    else
    {
        degrees = - RAD_TO_DEG * (Real) atan2( (double) y, (double) x );

        if( degrees < 0.0 )
            degrees += 360.0;

        return( degrees );
    }
}

private  int  find_limit_plane(
    int              n_points,
    Point            points[],
    Smallest_int     point_flags[],
    Point            *centre,
    Vector           *hinge,
    Vector           *normal )
{
    int      i, best_ind;
    Vector   horizontal, vertical, offset;
    Real     angle, best_angle, x, y;
    BOOLEAN  first;

    best_angle = 0.0;
    best_ind = -1;

    NORMALIZE_VECTOR( horizontal, *normal );
    CROSS_VECTORS( vertical, *normal, *hinge );
    NORMALIZE_VECTOR( vertical, vertical );

    first = TRUE;
    double plane_constant = -distance_from_plane( centre, &vertical, 0.0 );

    /* these are the total number of inputs points */
    for_less( i, 0, n_points ) {
        if( point_flags[i] & POINT_DISCARDED )
            continue;

        SUB_VECTORS( offset, points[i], *centre );
        x = -DOT_VECTORS( horizontal, offset );
        y = DOT_VECTORS( vertical, offset );

        if( x >= -TOLERANCE_2D && x <= TOLERANCE_2D &&
            y >= -TOLERANCE_2D && y <= TOLERANCE_2D )
            continue;

        angle = compute_clockwise_degrees( x, y ) - 180.0;
        if( angle < 0.0 )
            angle += 360.0;

        if( first || angle < best_angle )
        {
            if( angle < 90.0 - 0.1 || angle > 270.0 + 0.1 )
            {
                handle_internal_error( "find_limit_plane angle" );
                exit(1);
            }
            else
            {        
                best_angle = angle;
                best_ind = i;
                first = FALSE;
            }
        }
        else if( angle == best_angle )
        {
            if( distance_between_points( centre, &points[i] ) <
                distance_between_points( centre, &points[best_ind] ) )
            {
                best_ind = i;
            }
        }
    }

    if( best_ind < 0 ) {
        handle_internal_error( "find_limit_plane" );
        exit(1);
    }

    return( best_ind );
}

private  int  get_polygon_point_index(
    polygons_struct  *polygons,
    Point            points[],
    int              new_indices[],
    int              v )
{
    if( new_indices[v] < 0 )
    {
        new_indices[v] = polygons->n_points;
        ADD_ELEMENT_TO_ARRAY( polygons->points, polygons->n_points,
                              points[v], DEFAULT_CHUNK_SIZE );
    }

    return( new_indices[v] );
}

private  int  add_polygon(
    polygons_struct  *polygons,
    int              n_vertices,
    int              vertices[] )
{
    int   i, n_indices;

    n_indices = NUMBER_INDICES( *polygons );

    ADD_ELEMENT_TO_ARRAY( polygons->end_indices, polygons->n_items,
                          n_indices + n_vertices, DEFAULT_CHUNK_SIZE );

    SET_ARRAY_SIZE( polygons->indices, n_indices, n_indices + n_vertices,
                    DEFAULT_CHUNK_SIZE );

    for_less( i, 0, n_vertices ) {
        polygons->indices[n_indices+i] = vertices[i];

    }

    return( polygons->n_items - 1 );
}

private  int get_edge_key(
    polygons_struct              *polygons,
    int                          poly,
    int                          edge ) {

    int          p0, p1, size;

    size = GET_OBJECT_SIZE( *polygons, poly );

    p0 = polygons->indices[POINT_INDEX(polygons->end_indices,poly,edge)];
    p1 = polygons->indices[POINT_INDEX(polygons->end_indices,poly,
                           (edge+1)%size)];

    return( IJ(MIN( p0, p1 ),MAX( p0, p1 ),KeyFactor) );
}


private  void  add_edge_to_list(
    queue_struct                 *queue,
    hash_table_struct            *edge_table,
    polygons_struct              *polygons,
    int                          poly,
    int                          edge )
{
    int          key;
    edge_struct  edge_ptr;
    queue_entry  entry;

    key = get_edge_key( polygons, poly, edge );

    if( lookup_in_hash_table( edge_table, key, NULL ) ) {
      /* This is quite convoluted, but it works. CL */
      remove_from_hash_table( edge_table, key, (void *) &edge_ptr );
if( edge_ptr.ref_count == 2 ) {
  printf( " bad count for edge %d %d\n", key/KeyFactor, key%KeyFactor );
  // exit(1);
}
      edge_ptr.ref_count++;
      insert_in_hash_table( edge_table, key, (void *) &edge_ptr );
    } else {
      edge_ptr.ref_count = 1;
      insert_in_hash_table( edge_table, key, (void*)&edge_ptr );

      entry.poly = poly;
      entry.edge = edge;
      INSERT_IN_QUEUE( *queue, entry );
    }

}

private  int  get_convex_hull_2d(
    int              n_points,
    Real             x[],
    Real             y[],
    int              hull_indices[],
    int              e0,
    int              e1 ) {

    int      i, j, min_ind, n_in_hull, current_ind, best_ind;
    Real     dx, dy, best_len, cross;

    n_in_hull = 0;
    if( e0 < 0 && e1 < 0 ) {
      min_ind = 0;
      for_less( i, 1, n_points ) {
          if( x[i] < x[min_ind] )
              min_ind = i;
          else if( x[i] == x[min_ind] && y[i] < y[min_ind] )
              min_ind = i;
      }
      current_ind = min_ind;
    } else {
      min_ind = e1;
      hull_indices[n_in_hull] = e1;
      n_in_hull++;
      current_ind = e0;
    }

    do {
      if( n_in_hull >= n_points ) {
            handle_internal_error( "get_convex_hull_2d" );
            printf( "\n" );
            for( i = 0; i < n_in_hull; i++ ) {
              printf( " %g %g\n", x[hull_indices[i]], y[hull_indices[i]] );
            }
            exit(1);
        }
        hull_indices[n_in_hull] = current_ind;
        ++n_in_hull;
        best_len = 0.0;
        best_ind = -1;

        for_less( i, 0, n_points ) {

          if( i == current_ind ) continue;

          dx = x[i] - x[current_ind];
          dy = y[i] - y[current_ind];

          for( j = 0; j < n_points; j++ ) {
            if( j == current_ind || j == i ) continue;
            cross = dx * ( y[j] - y[current_ind] ) - dy * ( x[j] - x[current_ind] );
            if( cross < -TOLERANCE_DISTANCE ) break;
          }
          if( j == n_points ) {
            if( dx * dx + dy * dy > best_len ) {
              best_len = dx * dx + dy * dy;
              best_ind = i;
            }
          }
        }
        if( best_ind >= 0 ) {
          current_ind = best_ind;
        } else {
          printf( "could not find best index\n" );
          exit(1);
        }
    } while( current_ind != min_ind );

    return( n_in_hull );
}

private  int   get_plane_polygon_vertices(
    int              n_points,
    Point            points[],
    Smallest_int     point_flags[],
    int              p0,
    int              p1,
    int              p2,
    int              e0,
    int              e1,
    int              vertices[] )
{
    int      i, n_in_hull, n_in_plane, *plane_points, *hull_points;
    Real     plane_constant, horiz_constant, *x, *y, dist;
    Vector   v01, v02, normal, offset, horizontal, vertical;

    SUB_POINTS( v01, points[p1], points[p0] );
    SUB_POINTS( v02, points[p2], points[p0] );
    CROSS_VECTORS( normal, v01, v02 );
    NORMALIZE_VECTOR( normal, normal );

    if( e0 < 0 && e1 < 0 ) {
      NORMALIZE_VECTOR( horizontal, v01 );
      CROSS_VECTORS( vertical, normal, horizontal );
      NORMALIZE_VECTOR( vertical, vertical );
    } else {
      SUB_POINTS( v01, points[e1], points[e0] );
      NORMALIZE_VECTOR( vertical, v01 );
      CROSS_VECTORS( horizontal, vertical, normal );
      NORMALIZE_VECTOR( horizontal, horizontal );
    }

    plane_constant = -distance_from_plane( &points[p0], &normal, 0.0 );
    horiz_constant = -distance_from_plane( &points[p0], &horizontal, 0.0 );

    n_in_plane = 0;
    plane_points = NULL;

    for_less( i, 0, n_points ) {
        if( point_flags[i] & POINT_DISCARDED )
            continue;
        if( i == p0 || i == p1 || i == p2 ) {
          ADD_ELEMENT_TO_ARRAY( plane_points, n_in_plane, i, 10 );
        } else {
          dist = distance_from_plane( &points[i], &normal, plane_constant );
          if( dist >= -TOLERANCE_DISTANCE ) {
            dist = distance_from_plane( &points[i], &horizontal, horiz_constant );
            if( dist >= -TOLERANCE_DISTANCE ) {
              ADD_ELEMENT_TO_ARRAY( plane_points, n_in_plane, i, 10 );
            }
          }
        }
    }

    if( n_in_plane < 3 )
        handle_internal_error( "get_plane_polygon_vertices" );

    ALLOC( x, n_in_plane );
    ALLOC( y, n_in_plane );
    ALLOC( hull_points, n_in_plane );

    int i0 = -1, i1 = -1;
    for_less( i, 0, n_in_plane ) {
        SUB_POINTS( offset, points[plane_points[i]], points[p0] );
        x[i] = DOT_VECTORS( offset, horizontal );
        y[i] = DOT_VECTORS( offset, vertical );
        if( ! (point_flags[plane_points[i]] & POINT_USED_IN_CONVEX_HULL) )
            point_flags[plane_points[i]] |= POINT_DISCARDED;
        if( plane_points[i] == e0 ) i0 = i;
        if( plane_points[i] == e1 ) i1 = i;
    }

    n_in_hull = get_convex_hull_2d( n_in_plane, x, y, hull_points, i0, i1 );

    for_less( i, 0, n_in_hull ) {
        // Correction is not so effective in single precision.
        int ii = plane_points[hull_points[i]];
        dist = distance_from_plane( &points[ii], &normal, plane_constant );
        points[ii].coords[0] -= dist * normal.coords[0];
        points[ii].coords[1] -= dist * normal.coords[1];
        points[ii].coords[2] -= dist * normal.coords[2];
        point_flags[ii] = POINT_USED_IN_CONVEX_HULL;
        vertices[i] = ii;
    }

    FREE( hull_points );
    FREE( x );
    FREE( y );
    FREE( plane_points );

    return( n_in_hull );
}

private  void  get_convex_hull(
    int              n_points,
    Point            points[],
    polygons_struct  *polygons )
{

    int                          i, min_ind, ind, second_ind, size;
    int                          n_bad_ref_count, n_edges;
    Vector                       hinge, new_hinge, normal, new_normal;
    int                          *new_indices, other_index, key;
    int                          poly, edge, new_poly;
    int                          n_vertices, *vertices;
    int                          *poly_vertices;
    Point                        *poly_points;
    Smallest_int                 *point_flags;
    queue_entry                  entry;
    queue_struct                 queue;
    edge_struct                  edge_ptr;
    hash_table_struct            edge_table;
    hash_table_pointer           hash_ptr;

    initialize_polygons( polygons, WHITE, NULL );

    if( n_points == 0 )
        return;

    min_ind = 0;
    for_less( i, 0, n_points ) {
        if( i == 0 || Point_x(points[i]) < Point_x(points[min_ind]) )
            min_ind = i;
        else if( Point_x(points[i]) == Point_x(points[min_ind]) &&
                 Point_y(points[i]) <  Point_y(points[min_ind]) )
            min_ind = i;
        else if( Point_x(points[i]) == Point_x(points[min_ind]) &&
                 Point_y(points[i]) == Point_y(points[min_ind]) &&
                 Point_z(points[i]) <  Point_z(points[min_ind]) )
            min_ind = i;
    }

    fill_Vector( hinge, 0.0, 0.0, 1.0 );
    fill_Vector( normal, -1.0, 0.0, 0.0 );

    ALLOC( point_flags, n_points );

    for_less( i, 0, n_points )
        point_flags[i] = FALSE;

    ind = find_limit_plane( n_points, points, point_flags,
                            &points[min_ind], &hinge, &normal );

    SUB_POINTS( new_hinge, points[ind], points[min_ind] );
    CROSS_VECTORS( new_normal, hinge, new_hinge );

    second_ind = find_limit_plane( n_points, points, point_flags,
                                   &points[ind], &new_hinge, &new_normal );

    if( min_ind == ind || min_ind == second_ind || ind == second_ind )
        handle_internal_error( "get_convex_hull" );

    ALLOC( vertices, n_points );
    ALLOC( poly_vertices, n_points );
    ALLOC( poly_points, n_points );

    n_vertices = get_plane_polygon_vertices( n_points, points, point_flags, 
                                             min_ind, ind, second_ind,
                                             -1, -1, vertices );
   
    poly = add_polygon( polygons, n_vertices, vertices );

    INITIALIZE_QUEUE( queue );

    initialize_hash_table( &edge_table, 1000, sizeof(edge_struct),
                           ENLARGE_THRESHOLD, NEW_DENSITY );

    for_less( i, 0, n_vertices ) {
        add_edge_to_list( &queue, &edge_table, polygons, poly, i );
    }

    while( !IS_QUEUE_EMPTY( queue ) ) {
        REMOVE_FROM_QUEUE( queue, entry );

        poly = entry.poly;
        edge = entry.edge;

        key = get_edge_key( polygons, poly, edge );

        if( !lookup_in_hash_table( &edge_table, key, (void *) &edge_ptr ) ) {
            handle_internal_error( "Convex hull" );
        }

        if( edge_ptr.ref_count >= 2 )
            continue;

        size = GET_OBJECT_SIZE( *polygons, poly );

        for_less( i, 0, size ) {
            poly_vertices[i] = polygons->indices[
                            POINT_INDEX(polygons->end_indices,poly,i)];
            poly_points[i].coords[0] = points[poly_vertices[i]].coords[0];
            poly_points[i].coords[1] = points[poly_vertices[i]].coords[1];
            poly_points[i].coords[2] = points[poly_vertices[i]].coords[2];
        }

        SUB_POINTS( hinge, points[poly_vertices[edge]],
                           points[poly_vertices[(edge+1)%size]] );

        find_polygon_normal( size, poly_points, &normal );

        ind = find_limit_plane( n_points, points, point_flags,
                                &points[poly_vertices[edge]],
                                &hinge, &normal );

        other_index = ind;

        n_vertices = get_plane_polygon_vertices( n_points, points, point_flags,
                 poly_vertices[edge], other_index, poly_vertices[(edge+1)%size], 
                 poly_vertices[edge], poly_vertices[(edge+1)%size],
                 vertices );

        new_poly = add_polygon( polygons, n_vertices, vertices );

        for( i = 0; i < n_vertices; i++ ) {
            add_edge_to_list( &queue, &edge_table, polygons, new_poly, i );
        }
    }

    DELETE_QUEUE( queue );
    if( vertices ) FREE( vertices );
    if( poly_vertices ) FREE( poly_vertices );
    if( poly_points ) FREE( poly_points );
    if( point_flags ) FREE( point_flags );

    initialize_hash_pointer( &hash_ptr );

    n_bad_ref_count = 0;
    n_edges = 0;
    while( get_next_hash_entry( &edge_table, &hash_ptr, (void *) &edge_ptr ) ) {
        if( edge_ptr.ref_count != 2 ) {
            printf( "bad ref_count is %d for edge %d\n", edge_ptr.ref_count, n_edges );
            ++n_bad_ref_count;
        }
        ++n_edges;
    }

    delete_hash_table( &edge_table );

    if( n_bad_ref_count > 0 )
        print( "N ref counts != 2: %d/%d\n", n_bad_ref_count, n_edges );

    // Renumber the vertices of the convex hull locally.

    ALLOC( new_indices, n_points );
    for_less( i, 0, n_points ) {
        new_indices[i] = -1;
    }
    for( i = 0; i < polygons->end_indices[polygons->n_items-1]; i++ ) {
      polygons->indices[i] = get_polygon_point_index( polygons, points, new_indices, 
                                                      polygons->indices[i] );
    }

    if( new_indices ) FREE( new_indices );

    if( polygons->n_points > 0 ) {
        ALLOC( polygons->normals, polygons->n_points );
        compute_polygon_normals( polygons );
    }
}

int
calculate_convex_hull(
    polygons_struct *surface)
{

    int            n_points;
    Point          *points;
    object_struct  *object;
 
    n_points = get_points_of_region( surface, points );
    KeyFactor = n_points;

    object = create_object( POLYGONS );

    get_convex_hull( n_points, points, get_polygons_ptr(object) );

    check_polygons_neighbours_computed( get_polygons_ptr(object) );

    (void) output_graphics_file( "test.obj", ASCII_FORMAT, 1, &object );

    return( 0 );
}

double
gyrification_index_sph(polygons_struct *surface, polygons_struct *sphere,
                      char *file, int n_triangles, polygons_struct *reparam,
                      int smoothflag)
{
        polygons_struct *polygons;
        object_struct **object;
        double *rcx, *icx, *rcy, *icy, *rcz, *icz;
        double *lrcx, *licx, *lrcy, *licy, *lrcz, *licz;
        double *rdatax, *rdatay, *rdataz, *spectral_power;
        double *areas, *orig_areas, gi, orig_area, **sph_areas, *local_gi;
        int it, bw2, l, m, i;
        char str[80];

        bw2 = BW * BW;

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

        get_equally_sampled_coords_of_polygon(surface, sphere, BW,
                                              rdatax, rdatay, rdataz);

        get_sph_coeffs_of_realdata(rdatax, BW, DATAFORMAT, rcx, icx);
        get_sph_coeffs_of_realdata(rdatay, BW, DATAFORMAT, rcy, icy);
        get_sph_coeffs_of_realdata(rdataz, BW, DATAFORMAT, rcz, icz);

        get_realdata_from_sph_coeffs(rdatax, BW, DATAFORMAT, rcx, icx);
        get_realdata_from_sph_coeffs(rdatay, BW, DATAFORMAT, rcy, icy);
        get_realdata_from_sph_coeffs(rdataz, BW, DATAFORMAT, rcz, icz);

        object = (object_struct **) malloc(sizeof(object_struct *));
        *object = create_object(POLYGONS);
        polygons = get_polygons_ptr(*object);

        sample_sphere_from_sph(rdatax, rdatay, rdataz, polygons,
                               n_triangles, reparam, BW);

        orig_areas = (double *) malloc(sizeof(double) * polygons->n_items);
        orig_area = get_area_of_polygons(polygons, orig_areas);

        delete_polygons(polygons);

        areas = (double *) malloc(sizeof(double) * 10);
        sph_areas = (double **) malloc(sizeof(double *) * 10);
        for (i = 0; i < 10; i++) {
                sph_areas[i] = (double *) malloc(sizeof(double) *
                                                 polygons->n_items);
        }

        for (it = 0; it < 10; it++) {
/*                limit_bandwidth(BW, (int) bws[it], rcx, lrcx);
                limit_bandwidth(BW, (int) bws[it], rcy, lrcy);
                limit_bandwidth(BW, (int) bws[it], rcz, lrcz);
                limit_bandwidth(BW, (int) bws[it], icx, licx);
                limit_bandwidth(BW, (int) bws[it], icy, licy);
                limit_bandwidth(BW, (int) bws[it], icz, licz);
*/
                get_realdata_from_sph_coeffs(rdatax, BW, DATAFORMAT,
                                             lrcx, licx);
                get_realdata_from_sph_coeffs(rdatay, BW, DATAFORMAT,
                                             lrcy, licy);
                get_realdata_from_sph_coeffs(rdataz, BW, DATAFORMAT,
                                             lrcz, licz);

                sample_sphere_from_sph(rdatax, rdatay, rdataz,
                                       polygons, n_triangles, reparam, BW);

                areas[it] = get_area_of_polygons(polygons, sph_areas[it]) /
                            orig_area;
                for (i = 0; i < polygons->n_items; i++) {
                        if (orig_areas[i] != 0)
                                sph_areas[it][i] /= orig_areas[i];
                }

                delete_polygons(polygons);
        }

        local_gi = (double *) malloc(sizeof(double) * reparam->n_points);
//        get_localgi(reparam, bws, sph_areas, 10, local_gi, smoothflag);
//        output_values_any_format(file, reparam->n_points, local_gi,
//                                 TYPE_DOUBLE);
 
//        gi = get_globalgi(bws, areas, 10);

        free(rcx); free(rcy); free(rcz); 
        free(icx); free(icy); free(icz); 
        free(lrcx); free(lrcy); free(lrcz);
        free(licx); free(licy); free(licz);
        free(rdatax); free(rdatay); free(rdataz);
        
        return gi;
}
