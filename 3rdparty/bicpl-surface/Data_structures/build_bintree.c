/* ----------------------------------------------------------------------------
@COPYRIGHT  :
              Copyright 1993,1994,1995 David MacDonald,
              McConnell Brain Imaging Centre,
              Montreal Neurological Institute, McGill University.
              Permission to use, copy, modify, and distribute this
              software and its documentation for any purpose and without
              fee is hereby granted, provided that the above copyright
              notice appear in all copies.  The author and McGill University
              make no representations about the suitability of this
              software for any purpose.  It is provided "as is" without
              express or implied warranty.
---------------------------------------------------------------------------- */

#include "bicpl_internal.h"

#ifndef lint
static char rcsid[] = "$Header: /private-cvsroot/libraries/bicpl/Data_structures/build_bintree.c,v 1.13 2005/08/17 22:31:12 bert Exp $";
#endif

#define  NODE_VISIT_COST        0.02
#define  NET_CHANGE_THRESHOLD   0.0

typedef  struct
{
    bintree_node_struct    **ptr_to_node;
    range_struct           limits;
} leaf_queue_type;

typedef     PRIORITY_QUEUE_STRUCT( leaf_queue_type )   leaf_queue_struct;

/*---------------- private functions ------------------------- */

static  void  subdivide_bintree(
    bintree_struct_ptr   bintree,
    int                  max_nodes,
    int                  n_objects,
    range_struct         bound_vols[] );
static  void  split_node(
    range_struct          bound_vols[],
    bintree_node_struct   **ptr_to_node,
    range_struct          *limits,
    int                   *n_nodes,
    range_struct          *left_limits,
    Real                  *left_cost,
    range_struct          *right_limits,
    Real                  *right_cost );
static  void  create_leaf_queue(
    leaf_queue_struct   *leaf_queue );
static  void  delete_leaf_queue(
    leaf_queue_struct   *leaf_queue );
static  BOOLEAN  leaf_queue_empty(
    leaf_queue_struct   *leaf_queue );
static  void  insert_in_leaf_queue(
    leaf_queue_struct     *leaf_queue,
    bintree_node_struct   **ptr_to_node,
    Real                  node_cost,
    range_struct          *limits );
static  void  remove_from_leaf_queue(
    leaf_queue_struct     *leaf_queue,
    bintree_node_struct   ***ptr_to_node,
    range_struct          *limits );
static  BOOLEAN  node_tightly_bounds_object(
    range_struct  *bound_vol,
    range_struct  *limits );
static  void  recursive_efficiency_count(
    bintree_node_struct   *node,
    range_struct          *limits,
    Real                  *avg_nodes_visited,
    Real                  *avg_objects_visited );
static  Real  node_visit_estimation(
    range_struct  *limits );

/* ---------------------------------------------------------- */

#define  FACTOR 1.0e-4

/* ----------------------------- MNI Header -----------------------------------
@NAME       : create_object_bintree
@INPUT      : n_objects
              bound_vols
              max_nodes
@OUTPUT     : bintree
@RETURNS    : 
@DESCRIPTION: Creates a bintree of objects specified by their bounding volumes.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  void  create_object_bintree(
    int                  n_objects,
    range_struct         bound_vols[],
    bintree_struct_ptr   bintree,
    int                  max_nodes )
{
    int       i, c;
    Real      avg_nodes, avg_objects, limits[N_DIMENSIONS][2], size;
#ifdef  DEBUG
    Real      best_objects;
#endif

    for_less( i, 0, n_objects )
    {
        for_less( c, 0, N_DIMENSIONS )
        {
            size = (Real) bound_vols[i].limits[c][1] -
                   (Real) bound_vols[i].limits[c][0];
            bound_vols[i].limits[c][0] -= (float) (size * FACTOR);
            bound_vols[i].limits[c][1] += (float) (size * FACTOR);
        }

        if( i == 0 )
        {
            limits[X][0] = (Real) bound_vols[i].limits[X][0];
            limits[Y][0] = (Real) bound_vols[i].limits[Y][0];
            limits[Z][0] = (Real) bound_vols[i].limits[Z][0];
            limits[X][1] = (Real) bound_vols[i].limits[X][1];
            limits[Y][1] = (Real) bound_vols[i].limits[Y][1];
            limits[Z][1] = (Real) bound_vols[i].limits[Z][1];
        }
        else
        {
            for_less( c, 0, N_DIMENSIONS )
            {
                if( (Real) bound_vols[i].limits[c][0] < limits[c][0] )
                    limits[c][0] = (Real) bound_vols[i].limits[c][0];
                if( (Real) bound_vols[i].limits[c][1] > limits[c][1] )
                    limits[c][1] = (Real) bound_vols[i].limits[c][1];
            }
        }
    }

    initialize_bintree( limits[X][0], limits[X][1],
                        limits[Y][0], limits[Y][1],
                        limits[Z][0], limits[Z][1], bintree );

    subdivide_bintree( bintree, max_nodes, n_objects, bound_vols );

    evaluate_bintree_efficiency( bintree, &avg_nodes, &avg_objects );

#ifdef  DEBUG
    best_objects = 0.0;
    for_less( i, 0, n_objects )
        best_objects += node_visit_estimation( &bound_vols[i] );
    best_objects /= node_visit_estimation( &bintree->range );

    print( "Est Nodes Visit: %g    Est Objects Visit %g  (Best possible: %g\n",
           avg_nodes, avg_objects, best_objects );

    if( n_objects > 100 && avg_objects > (Real) n_objects * 0.1 )
    {
        print( "Warning, bintree not efficient: n_objects = %d, ",
               n_objects );
        print( "avg_nodes %g, avg_objects %g\n", avg_nodes, avg_objects );
    }
#endif
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : subdivide_bintree
@INPUT      : bintree
              max_nodes
              n_objects
              bound_vols
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: Creates a subtree by prioritized subdivision of leaves until
              the max_nodes is reached.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

static  void  subdivide_bintree(
    bintree_struct_ptr   bintree,
    int                  max_nodes,
    int                  n_objects,
    range_struct         bound_vols[] )
{
    bintree_node_struct     **ptr_to_node;
    bintree_node_struct     **ptr_to_left_child, **ptr_to_right_child;
    Real                    left_cost, right_cost;
    range_struct            left_limits, right_limits;
    leaf_queue_struct       leaf_queue;
    int                     n_nodes;
    int                     i, *object_list;
    range_struct            limits;
    progress_struct         progress;

    ALLOC( object_list, n_objects );
    for_less( i, 0, n_objects )
        object_list[i] = i;

    bintree->root = create_bintree_leaf( 0.0, n_objects, object_list );

    FREE( object_list );

    create_leaf_queue( &leaf_queue );

    get_bintree_limits( bintree, &limits );

    insert_in_leaf_queue( &leaf_queue, &bintree->root, 1.0, &limits );

    initialize_progress_report( &progress, FALSE, max_nodes+2,
                                "Creating bintree" );

    n_nodes = 1;

    while( n_nodes < max_nodes && !leaf_queue_empty( &leaf_queue ) )
    {
        remove_from_leaf_queue( &leaf_queue, &ptr_to_node, &limits );

        split_node( bound_vols, ptr_to_node, &limits, &n_nodes,
                    &left_limits, &left_cost, &right_limits, &right_cost );

        if( get_bintree_left_child_ptr( *ptr_to_node, &ptr_to_left_child ) )
        {
            insert_in_leaf_queue( &leaf_queue, ptr_to_left_child, left_cost,
                                  &left_limits );
        }

        if( get_bintree_right_child_ptr( *ptr_to_node, &ptr_to_right_child ) )
        {
            insert_in_leaf_queue( &leaf_queue, ptr_to_right_child, right_cost,
                                  &right_limits );
        }

        update_progress_report( &progress, n_nodes );
    }

    terminate_progress_report( &progress );

    delete_leaf_queue( &leaf_queue );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : get_planes
@INPUT      : axis_index
              n_objects
              list
              bound_vols
              right_plane
@OUTPUT     : left_plane
              n_left
              n_overlap
              n_right
@RETURNS    : 
@DESCRIPTION: Computes the plane position for a node.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

static  void   get_planes(
    int             axis_index,
    int             n_objects,
    int             list[],
    range_struct    bound_vols[],
    Real            right_plane,
    Real            *left_plane,
    int             *n_left,
    int             *n_overlap,
    int             *n_right )
{
    int        i;
    Real       left_side, right_side;
    BOOLEAN    found;

    found = FALSE;
    *left_plane = (Real) bound_vols[list[0]].limits[axis_index][0];

    for_less( i, 0, n_objects )
    {
        left_side = (Real) bound_vols[list[i]].limits[axis_index][0];
        right_side = (Real) bound_vols[list[i]].limits[axis_index][1];
        if( left_side < right_plane &&
            (!found || right_side > *left_plane ) )
        {
            *left_plane = right_side;
            found = TRUE;
        }
    }

    *n_left = 0;
    *n_overlap = 0;
    *n_right = 0;

    for_less( i, 0, n_objects )
    {
        left_side = (Real) bound_vols[list[i]].limits[axis_index][0];
        right_side = (Real) bound_vols[list[i]].limits[axis_index][1];
        if( left_side < right_plane )
            ++(*n_left);
        else if( !found || right_side > *left_plane )
            ++(*n_right);
        else
            ++(*n_overlap);
    }
}

#define  THRESHOLD  30

/* ----------------------------- MNI Header -----------------------------------
@NAME       : find_best_split_node_for_axis
@INPUT      : axis_index
              n_objects
              object_list
              bound_vols
              limits
@OUTPUT     : best_left_plane
              best_right_plane
@RETURNS    : 
@DESCRIPTION: Finds a good split plane for the bintree.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

static  Real  find_best_split_node_for_axis(
    int                   axis_index,
    int                   n_objects,
    int                   object_list[],
    range_struct          bound_vols[],
    range_struct          *limits,
    Real                  *best_left_plane,
    Real                  *best_right_plane )
{
    int                  i, n_left, n_right, n_overlap, step;
    Real                 left_plane, right_plane;
    float                save_plane;
    Real                 cost, best_cost, size_left, size_right;

    best_cost = 0.0;

    if( n_objects > THRESHOLD )
        step = n_objects+1;
    else
        step = 1;

    for( i = 0;   i <= n_objects;  i += step )
    {
        if( n_objects > THRESHOLD )
            right_plane = ((Real) limits->limits[axis_index][0] +
                           (Real) limits->limits[axis_index][1]) / 2.0;
        else if( i < n_objects )
        {
            right_plane = (Real) bound_vols[object_list[i]].
                                          limits[axis_index][0];
        }
        else
            right_plane = (Real) limits->limits[axis_index][1];

        get_planes( axis_index, n_objects, object_list, bound_vols,
                    right_plane, &left_plane,
                    &n_left, &n_overlap, &n_right );

        save_plane = limits->limits[axis_index][1];
        limits->limits[axis_index][1] = (float) left_plane;
        size_left = node_visit_estimation( limits );
        limits->limits[axis_index][1] = save_plane;

        save_plane = limits->limits[axis_index][0];
        limits->limits[axis_index][0] = (float) right_plane;
        size_right = node_visit_estimation( limits );
        limits->limits[axis_index][0] = save_plane;

        if( size_left < size_right )
            n_left += n_overlap;
        else
            n_right += n_overlap;

        cost = size_left * (Real) n_left + size_right * (Real) n_right;

        if( i == 0 || cost < best_cost )
        {
            best_cost = cost;
            *best_right_plane = right_plane;
            *best_left_plane = left_plane;
        }
    }

    return( best_cost );
}

#define  SWAP( a, b, type ) \
    { \
        type  _tmp; \
        _tmp = (a); \
        (a) = (b); \
        (b) = _tmp; \
    }

#ifdef  DEBUG
static  void  check_objects(
    range_struct          *limits,
    int                   n_objects,
    int                   object_list[],
    range_struct          bound_vols[] )
{
    int   i, c;
    for_less( i, 0, n_objects )
        for_less( c, 0, N_DIMENSIONS )
            if( bound_vols[object_list[i]].limits[c][0] < limits->limits[c][0]||
                bound_vols[object_list[i]].limits[c][1] > limits->limits[c][1] )
                handle_internal_error( "check_objects" );
}
#endif

/* ----------------------------- MNI Header -----------------------------------
@NAME       : split_node
@INPUT      : 
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: 
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

static  void  split_node(
    range_struct          bound_vols[],
    bintree_node_struct   **ptr_to_node,
    range_struct          *limits,
    int                   *n_nodes,
    range_struct          *left_limits,
    Real                  *left_cost,
    range_struct          *right_limits,
    Real                  *right_cost )
{
    int                  i, c, n_objects;
    int                  *object_list;
    int                  axis_index, top, bottom, start_dim, end_dim;
    BOOLEAN              overlap_in_left, tightly_bounds;
    Real                 left_plane, right_plane, split_position;
    Real                 min_plane, max_plane, best, cost;
    Real                 test_left, test_right;
    Real                 node_cost, previous_cost, net_change;
    bintree_node_struct  *left_child, *right_child;

    n_objects = get_bintree_leaf_objects( *ptr_to_node, &object_list );

    if( n_objects == 0 )
        return;

    tightly_bounds = TRUE;

    for_less( i, 0, n_objects )
    {
        if( !node_tightly_bounds_object( &bound_vols[object_list[i]], limits ) )
        {
            tightly_bounds = FALSE;
            break;
        }
    }

    if( tightly_bounds )
        return;

    best = 0.0;

    if( n_objects > THRESHOLD )
    {
        start_dim = 0;
        for_less( c, 0, N_DIMENSIONS )
        {
            if( limits->limits[c][1] - limits->limits[c][0] >
                limits->limits[start_dim][1] - limits->limits[start_dim][0] )
            {
                start_dim = c;
            }
        }
        end_dim = start_dim + 1;
    }
    else
    {
        start_dim = 0;
        end_dim = N_DIMENSIONS;
    }

    for_less( c, start_dim, end_dim )
    {
        cost = find_best_split_node_for_axis( c, n_objects, object_list,
                                              bound_vols, limits,
                                              &test_left, &test_right );

        if( c == start_dim || cost < best )
        {
            best = cost;
            axis_index = c;
            left_plane = test_left;
            right_plane = test_right;
        }
    }

    previous_cost = (Real) n_objects * node_visit_estimation( limits );
    node_cost = NODE_VISIT_COST * node_visit_estimation( limits );

    net_change = node_cost + cost - previous_cost;

    if( net_change >= NET_CHANGE_THRESHOLD )
        return;

    min_plane = (Real) limits->limits[axis_index][0];
    max_plane = (Real) limits->limits[axis_index][1];

    if( left_plane - min_plane < max_plane - right_plane )
        overlap_in_left = TRUE;
    else
        overlap_in_left = FALSE;

    bottom = 0;
    top = n_objects - 1;

    while( bottom <= top )
    {
        while( bottom < n_objects &&
               (Real) bound_vols[object_list[bottom]].limits[axis_index][1] <=
                                                       left_plane &&
               (overlap_in_left ||
                (Real) bound_vols[object_list[bottom]].limits[axis_index][0] <
                                                       right_plane) )
            ++bottom;
        while( top >= 0 &&
               (Real) bound_vols[object_list[top]].limits[axis_index][0] >=
                                                     right_plane &&
               (!overlap_in_left ||
               (Real) bound_vols[object_list[top]].limits[axis_index][1] >
                                                     left_plane) )
            --top;
        if( bottom < top )
        {
            SWAP( object_list[bottom], object_list[top], int )
        }
    }

    if( (bottom == 0 && right_plane == (Real) limits->limits[axis_index][0]) ||
        (bottom == n_objects &&
            left_plane == (Real) limits->limits[axis_index][1]) )
    {
        return;
    }

    if( bottom > 0 )
    {
        *left_limits = *limits;
        left_limits->limits[axis_index][1] = (float) left_plane;

        *left_cost = (Real) bottom * node_visit_estimation( left_limits );

        left_child = create_bintree_leaf( left_plane, bottom, object_list );

#ifdef  DEBUG
        check_objects( left_limits, bottom, object_list, bound_vols );
#endif

        ++(*n_nodes);
    }
    else
        left_child = (bintree_node_struct *) NULL;

    if( bottom < n_objects )
    {
        *right_limits = *limits;
        right_limits->limits[axis_index][0] = (float) right_plane;

        *right_cost = (Real) (n_objects - bottom) *
                           node_visit_estimation( right_limits );

        right_child = create_bintree_leaf( right_plane, n_objects - bottom,
                                           &object_list[bottom] );
#ifdef  DEBUG
        check_objects( right_limits, n_objects - bottom, &object_list[bottom],
                       bound_vols );
#endif
        ++(*n_nodes);
    }
    else
        right_child = (bintree_node_struct *) NULL;

    if( left_child == (bintree_node_struct *) NULL &&
        right_child == (bintree_node_struct *) NULL )
        handle_internal_error( "Split bintree node" );

    /* --- replace the leaf with an internal node */

    split_position = get_node_split_position( *ptr_to_node );
    delete_bintree_node( *ptr_to_node );
    *ptr_to_node = create_bintree_internal_node( axis_index,
                                                 split_position,
                                                 left_child,
                                                 right_child );
}

static  void  create_leaf_queue(
    leaf_queue_struct   *leaf_queue )
{
    INITIALIZE_PRIORITY_QUEUE( *leaf_queue );
}

static  void  delete_leaf_queue(
    leaf_queue_struct   *leaf_queue )
{
    DELETE_PRIORITY_QUEUE( *leaf_queue );
}

static  BOOLEAN  leaf_queue_empty(
    leaf_queue_struct   *leaf_queue )
{
    return( IS_PRIORITY_QUEUE_EMPTY( *leaf_queue ) );
}

static  void  insert_in_leaf_queue(
    leaf_queue_struct     *leaf_queue,
    bintree_node_struct   **ptr_to_node,
    Real                  node_cost,
    range_struct          *limits )
{
    leaf_queue_type   entry;

    entry.ptr_to_node = ptr_to_node;
    entry.limits = *limits;

    INSERT_IN_PRIORITY_QUEUE( *leaf_queue, entry, node_cost );
}

static  void  remove_from_leaf_queue(
    leaf_queue_struct     *leaf_queue,
    bintree_node_struct   ***ptr_to_node,
    range_struct          *limits )
{
    leaf_queue_type   entry;
    Real              node_cost;

    REMOVE_FROM_PRIORITY_QUEUE( *leaf_queue, entry, node_cost );

    if( node_cost < 0.0 )
        handle_internal_error( "remove from leaf queue" );

    *ptr_to_node = entry.ptr_to_node;
    *limits = entry.limits;
}

static  BOOLEAN  node_tightly_bounds_object(
    range_struct  *bound_vol,
    range_struct  *limits )
{
    int      l, c;
    
    for_less( c, 0, N_DIMENSIONS )
    {
        for_less( l, 0, 2 )
        {
            if( bound_vol->limits[c][l] != limits->limits[c][l] )
            {
                return( FALSE );
            }
        }
    }

    return( TRUE );
}

BICAPI  void  evaluate_bintree_efficiency(
    bintree_struct_ptr   bintree,
    Real                 *avg_nodes_visited,
    Real                 *avg_objects_visited )
{
    Real           n_visits_top_level;
    range_struct   limits;

    *avg_nodes_visited = 0.0;
    *avg_objects_visited = 0.0;

    limits = bintree->range;

    recursive_efficiency_count( bintree->root, &limits,
                                avg_nodes_visited, avg_objects_visited );

    n_visits_top_level = node_visit_estimation( &bintree->range );

    *avg_nodes_visited /= n_visits_top_level;
    *avg_objects_visited /= n_visits_top_level;
}

static  void  recursive_efficiency_count(
    bintree_node_struct   *node,
    range_struct          *limits,
    Real                  *avg_nodes_visited,
    Real                  *avg_objects_visited )
{
    float                 save_limit;
    Real                  node_visit_estimate;
    bintree_node_struct   *left_child, *right_child;
    Real                  left_limit, right_limit;
    int                   *object_list, axis_index;

    node_visit_estimate = node_visit_estimation( limits );

    *avg_nodes_visited += node_visit_estimate;

    if( bintree_node_is_leaf( node ) )
    {
        *avg_objects_visited += node_visit_estimate *
                     (Real) get_bintree_leaf_objects( node, &object_list );
    }
    else
    {
        axis_index = get_node_split_axis( node );

        if( get_bintree_left_child( node, &left_child ) )
        {
            left_limit = get_node_split_position( left_child );

            save_limit = limits->limits[axis_index][1];
            limits->limits[axis_index][1] = (float) left_limit;

            recursive_efficiency_count( left_child, limits,
                                        avg_nodes_visited,
                                        avg_objects_visited );
            limits->limits[axis_index][1] = save_limit;
        }

        if( get_bintree_right_child( node, &right_child ) )
        {
            right_limit = get_node_split_position( right_child );

            save_limit = limits->limits[axis_index][0];
            limits->limits[axis_index][0] = (float) right_limit;

            recursive_efficiency_count( right_child, limits,
                                        avg_nodes_visited,
                                        avg_objects_visited );
            limits->limits[axis_index][0] = save_limit;
        }
    }
}

static  Real  node_visit_estimation(
    range_struct  *limits )
{
#ifdef  VOLUME_EST
    return( range_volume(limits) );
#else
    return( range_surface_area(limits) );
#endif
}
