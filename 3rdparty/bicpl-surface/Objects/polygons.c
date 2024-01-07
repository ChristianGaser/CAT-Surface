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
static char rcsid[] = "$Header: /private-cvsroot/libraries/bicpl/Objects/polygons.c,v 1.44 2005/08/17 22:28:27 bert Exp $";
#endif


/*! \defgroup grp_polygons Polygons
 *  \ingroup grp_bicobj
 *
 * A polygons structure contains a list of vertices and
 * facets constructed from those vertices.
 *
 * There are \em n_points vertices and \em n_items
 * facets stored in the structure.
 *
 * @{
 */


static  void  reverse_polygon_order(
    polygons_struct   *polygons,
    int               poly );
static  Real  estimate_polygon_curvature(
    Point   *point,
    int     n_neighs,
    Point   neighs[],
    Point   *centroid,
    Vector  *normal,
    Real    *base_length );

/* ----------------------------- MNI Header -----------------------------------
@NAME       : initialize_polygons
@INPUT      : col
              spr   : optional
@OUTPUT     : polygons
@RETURNS    : 
@DESCRIPTION: Initializes the polygons to an empty set.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */
/*!
 * \brief generate empty set of polygons
 *
 * Generates a valid polygons instance, with the given colour
 * and surface properties.  
 *
 * @param polygons pointer to allocated structure
 * @param col colour for polygons
 * @param spr surface properties (or NULL for defaults)
 */
BICAPI  void  initialize_polygons(
    polygons_struct   *polygons,
    Colour            col,
    Surfprop          *spr )
{
    ALLOC( polygons->colours, 1 );

    polygons->colour_flag = ONE_COLOUR;

    polygons->colours[0] = col;

    if( spr != (Surfprop *) 0 )
        polygons->surfprop = *spr;
    else
        get_default_surfprop( &polygons->surfprop );

    polygons->n_points = 0;
    polygons->points = (Point *) 0;
    polygons->normals = (Vector *) 0;

    polygons->line_thickness = 1.0f;

    polygons->n_items = 0;
    polygons->end_indices = (int *) 0;
    polygons->indices = (int *) 0;

    polygons->visibilities = (Smallest_int *) 0;
    polygons->neighbours = (int *) 0;
    polygons->bintree = (bintree_struct_ptr) NULL;
}


/*!
 * \brief generate polygons with fixed-size facets
 *
 * Partially initialize a polygons instance.  To finish
 * the description, the points, normals, and indices arrays
 * must be filled in.  The rest of the fields (including
 * end_indices) are filled in by this routine.
 *
 * @param polygons pointer to allocated structure
 * @param col colour for polygons
 * @param spr surface properties (or NULL)
 * @param n_points number of vertices
 * @param n_polygons number of facets
 * @param size number of vertices per facet
 */
BICAPI  void  initialize_polygons_with_size(
    polygons_struct   *polygons,
    Colour            col,
    Surfprop          *spr,
    int               n_points,
    int               n_polygons,
    int               size )
{
    int   p;

    initialize_polygons( polygons, col, spr );

    polygons->n_points = n_points;
    ALLOC( polygons->points, polygons->n_points );
    ALLOC( polygons->normals, polygons->n_points );

    polygons->n_items = n_polygons;

    ALLOC( polygons->end_indices, polygons->n_items );

    for_less( p, 0, polygons->n_items )
        polygons->end_indices[p] = (p+1) * size;

    ALLOC( polygons->indices, polygons->end_indices[polygons->n_items-1] );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : free_polygon_neighbours
@INPUT      : polygons
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: Deletes the polygon neighbours.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  void  free_polygon_neighbours(
    polygons_struct   *polygons )
{
    if( polygons->neighbours != (int *) NULL )
    {
        FREE( polygons->neighbours );
        polygons->neighbours = (int *) NULL;
    }
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : delete_polygons
@INPUT      : polygons
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: Deletes the polygons.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  void  delete_polygons( polygons_struct *polygons )
{
    free_colours( polygons->colour_flag, polygons->colours, polygons->n_points,
                  polygons->n_items );

    if( polygons->n_points > 0 )
        FREE( polygons->points );

    if( polygons->n_points > 0 && polygons->normals != (Vector *) 0 )
        FREE( polygons->normals );

    if( polygons->n_items > 0 )
        FREE( polygons->indices );

    if( polygons->n_items > 0 )
        FREE( polygons->end_indices );

    if( polygons->visibilities != (Smallest_int *) 0 )
        FREE( polygons->visibilities );

    free_polygon_neighbours( polygons );

    delete_bintree_if_any( &polygons->bintree );

    polygons->visibilities = (Smallest_int *) 0;
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : copy_polygons
@INPUT      : src
@OUTPUT     : dest
@RETURNS    : 
@DESCRIPTION: Copies the polygons.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  void  copy_polygons(
    polygons_struct   *src,
    polygons_struct   *dest )
{
    int   i, n_indices, n_colours;

    *dest = *src;

    n_colours = get_n_colours( src->colour_flag, src->n_points, src->n_items );

    ALLOC( dest->colours, n_colours );
    for_less( i, 0, n_colours )
        dest->colours[i] = src->colours[i];

    ALLOC( dest->points, src->n_points );
    if( src->normals != (Vector *) NULL )
        ALLOC( dest->normals, src->n_points );

    for_less( i, 0, src->n_points )
    {
        dest->points[i] = src->points[i];
        if( src->normals != (Vector *) NULL )
            dest->normals[i] = src->normals[i];
        
    }

    ALLOC( dest->end_indices, src->n_items );
    for_less( i, 0, src->n_items )
        dest->end_indices[i] = src->end_indices[i];

    n_indices = NUMBER_INDICES( *src );
    ALLOC( dest->indices, n_indices );

    for_less( i, 0, n_indices )
        dest->indices[i] = src->indices[i];

    dest->visibilities = (Smallest_int *) 0;
    dest->neighbours = (int *) 0;
    dest->bintree = (bintree_struct_ptr) NULL;
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : set_polygon_per_item_colours
@INPUT      : polygons
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: Sets the polygons to per polygon colours.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  void  set_polygon_per_item_colours(
    polygons_struct   *polygons )
{
    Colour   colour;
    int      i;

    if( polygons->colour_flag != PER_ITEM_COLOURS )
    {
        colour = polygons->colours[0];
        polygons->colour_flag = PER_ITEM_COLOURS;
        REALLOC( polygons->colours, polygons->n_items );

        for_less( i, 0, polygons->n_items )
            polygons->colours[i] = colour;

    }
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : create_polygons_visibilities
@INPUT      : polygons
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: Allocates a visibility flag for each polygon.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  void  create_polygons_visibilities(
    polygons_struct   *polygons )
{
    int     i;

    if( polygons->visibilities == (Smallest_int *) 0 &&
        polygons->n_items > 0 )
    {
        ALLOC( polygons->visibilities, polygons->n_items );

        for_less( i, 0, polygons->n_items )
            polygons->visibilities[i] = TRUE;
    }
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : set_polygons_visibilities
@INPUT      : polygons
              state
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: Sets visibilities of all polygons to the state.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  void  set_polygons_visibilities(
    polygons_struct   *polygons,
    BOOLEAN           state )
{
    int   i;

    if( polygons->visibilities != (Smallest_int *) 0 )
    {
        for_less( i, 0, polygons->n_items )
            polygons->visibilities[i] = (Smallest_int) state;
    }
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : start_new_polygon
@INPUT      : polygons
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: Starts a new polygon in the structure, with size zero.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  void  start_new_polygon(
    polygons_struct *polygons )
{
    int      n_indices;

    n_indices = NUMBER_INDICES( *polygons );

    ADD_ELEMENT_TO_ARRAY( polygons->end_indices, polygons->n_items,
                          n_indices, DEFAULT_CHUNK_SIZE );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : add_point_to_polygon
@INPUT      : polygons
              point
              normal
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: Adds the given point and normal to the last polygon.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  void  add_point_to_polygon(
    polygons_struct   *polygons,
    Point             *point,
    Vector            *normal )
{
    int      n_points;

    if( polygons->n_items == 0 )
        start_new_polygon( polygons );

    if( polygons->n_points > 1 )
    {
        if( (normal != (Vector *) 0 && polygons->normals == (Vector *) 0) ||
            (normal == (Vector *) 0 && polygons->normals != (Vector *) 0) )
        {
            print_error(
              "Error: be consistent with normals in add_point_to_polygon.\n" );
        }
    }

    ADD_ELEMENT_TO_ARRAY( polygons->indices,
                          polygons->end_indices[polygons->n_items-1],
                          polygons->n_points, DEFAULT_CHUNK_SIZE );

    if( normal != (Vector *) 0 )
    {
        n_points = polygons->n_points;
        ADD_ELEMENT_TO_ARRAY( polygons->normals, n_points,
                              *normal, DEFAULT_CHUNK_SIZE );
    }

    ADD_ELEMENT_TO_ARRAY( polygons->points, polygons->n_points,
                          *point, DEFAULT_CHUNK_SIZE );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : get_polygon_points
@INPUT      : polygons
              poly
@OUTPUT     : points
@RETURNS    : number points
@DESCRIPTION: Passes back the points of the specified poly.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  int  get_polygon_points(
    polygons_struct   *polygons,
    int               poly,
    Point             points[] )
{
    int      size, p;

    size = GET_OBJECT_SIZE( *polygons, poly );

    for_less( p, 0, size )
    {
        points[p] = polygons->points[
                 polygons->indices[POINT_INDEX(polygons->end_indices,poly,p)]];
    }

    return( size );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : get_polygon_centroid
@INPUT      : polygons
              poly
@OUTPUT     : centroid
@RETURNS    : 
@DESCRIPTION: Passes back the centroid of the specified poly.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  void  get_polygon_centroid(
    polygons_struct   *polygons,
    int               poly,
    Point             *centroid )
{
    int      size, p;
    Point    point;

    fill_Point( *centroid, 0.0, 0.0, 0.0 );

    size = GET_OBJECT_SIZE( *polygons, poly );

    for_less( p, 0, size )
    {
        point = polygons->points[
                  polygons->indices[POINT_INDEX(polygons->end_indices,poly,p)]];
        ADD_POINTS( *centroid, *centroid, point );
    }

    if( size > 0 )
        SCALE_POINT( *centroid, *centroid, 1.0 / (Real) size );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : find_vertex_index
@INPUT      : polygons
              poly
              point_index
@OUTPUT     : 
@RETURNS    : index in poly of the point
@DESCRIPTION: Finds which vertex index of the poly corresponds to the
              point_index.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  int  find_vertex_index(
    polygons_struct  *polygons,
    int              poly,
    int              point_index )
{
    int   e, ind, size, p;

    ind = -1;

    size = GET_OBJECT_SIZE( *polygons, poly );

    for_less( e, 0, size )
    {
        p = polygons->indices[POINT_INDEX(polygons->end_indices,poly,e)];

        if( p == point_index )
        {
            ind = e;
            break;
        }
    }

    return( ind );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : find_edge_index
@INPUT      : polygons
              poly 
              point_index1 
              point_index2
@OUTPUT     : 
@RETURNS    : edge index
@DESCRIPTION: Finds the index of the edge in the poly between the two points.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  int  find_edge_index(
    polygons_struct  *polygons,
    int              poly,
    int              point_index1,
    int              point_index2 )
{
    int   e, ind, size, p1, p2;

    ind = -1;

    size = GET_OBJECT_SIZE( *polygons, poly );

    p2 = polygons->indices[POINT_INDEX(polygons->end_indices,poly,0)];
    for_less( e, 0, size )
    {
        p1 = p2;
        p2 = polygons->indices[POINT_INDEX(polygons->end_indices,poly,
                               (e+1)%size)];

        if( (p1 == point_index1 && p2 == point_index2) ||
            (p1 == point_index2 && p2 == point_index1) )
        {
            ind = e;
            break;
        }
    }

    return( ind );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : find_polygon_with_edge
@INPUT      : polygons
              point_index1
              point_index2
@OUTPUT     : poly_containing_edge
              edge_index
@RETURNS    : TRUE if found
@DESCRIPTION: Finds a polygon containing the given edge.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  BOOLEAN  find_polygon_with_edge(
    polygons_struct  *polygons,
    int              point_index1,
    int              point_index2,
    int              *poly_containing_edge,
    int              *edge_index )
{
    int   poly;

    for_less( poly, 0, polygons->n_items )
    {
        *edge_index = find_edge_index( polygons, poly,
                                       point_index1, point_index2 );

        if( *edge_index >= 0 )
        {
            *poly_containing_edge = poly;
            break;
        }
    }

    return( poly < polygons->n_items );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : find_polygon_with_vertex
@INPUT      : polygons
              point_index
@OUTPUT     : poly_index
              vertex_index
@RETURNS    : TRUE if found
@DESCRIPTION: Searches for a polygon containing the point.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  BOOLEAN  find_polygon_with_vertex(
    polygons_struct   *polygons,
    int               point_index,
    int               *poly_index,
    int               *vertex_index )
{
    int      poly, size, i;
    BOOLEAN  found;

    found = FALSE;

    for_less( poly, 0, polygons->n_items )
    {
        size = GET_OBJECT_SIZE( *polygons, poly );
        for_less( i, 0, size )
        {
            if( polygons->indices[POINT_INDEX(polygons->end_indices,poly,i)] ==
                point_index )
            {
                *poly_index = poly;
                *vertex_index = i;
                found = TRUE;
                break;
            }
        }
    }

    return( found );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : lookup_polygon_vertex
@INPUT      : polygons
              point
@OUTPUT     : point_index
@RETURNS    : TRUE if successful
@DESCRIPTION: Finds the point index of a given point.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  BOOLEAN  lookup_polygon_vertex(
    polygons_struct   *polygons,
    Point             *point,
    int               *point_index )
{
    int      i;

    for_less( i, 0, polygons->n_points )
    {
        if( EQUAL_POINTS( polygons->points[i], *point ) )
        {
            *point_index = i;
            break;
        }
    }

    return( i < polygons->n_points );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : find_next_edge_around_point
@INPUT      : polygons
              poly
              index_1
              index_2
@OUTPUT     : next_poly
              next_index_1
              next_index_2
@RETURNS    : TRUE if found
@DESCRIPTION: Given a polygon and two vertex indices within the poly,
              corresponding to an edge, finds the neighbouring poly, and
              the two vertex indices corresponding to the next 
              edge around the point corresponding to poly,index_1.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  BOOLEAN  find_next_edge_around_point(
    polygons_struct   *polygons,
    int               poly,
    int               index_1,
    int               index_2,
    int               *next_poly,
    int               *next_index_1,
    int               *next_index_2 )
{
    int                    size, edge, point_index, neighbour_point_index;
    int                    next_neigh_index;

    point_index = polygons->indices[
                       POINT_INDEX( polygons->end_indices, poly, index_1 )];

    neighbour_point_index = polygons->indices[
                       POINT_INDEX( polygons->end_indices, poly, index_2 )];

    size = GET_OBJECT_SIZE( *polygons, poly );

    if( index_2 == (index_1+1) % size )
        edge = index_1;
    else
        edge = index_2;
 
    *next_poly = polygons->neighbours[
                      POINT_INDEX(polygons->end_indices,poly,edge)];

    if( *next_poly >= 0 )
    {
        size = GET_OBJECT_SIZE(*polygons,*next_poly);
        *next_index_1 = find_vertex_index( polygons, *next_poly, point_index );

        *next_index_2 = (*next_index_1 + 1) % size;
        next_neigh_index = polygons->indices[
                           POINT_INDEX( polygons->end_indices, *next_poly,
                                        *next_index_2 )];
        if( next_neigh_index == neighbour_point_index )
        {
            *next_index_2 = (*next_index_1 + size- 1) % size;
        }
    }

    return( *next_poly >= 0 );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : get_neighbours_of_point
@INPUT      : polygons
              poly
              vertex_index
              max_neighbours  - max storage 
@OUTPUT     : neighbours
              interior_point
@RETURNS    : number of neighbours
@DESCRIPTION: Gets the points which are the neighbours of the specified point.
              interior_point is set to true if the point is totally
              surrounded by polygons.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : Sept. 12, 1996  D. MacDonald  : fixed it for the case where not
                                              a closed surface.
---------------------------------------------------------------------------- */

BICAPI  int  get_neighbours_of_point(
    polygons_struct   *polygons,
    int               poly,
    int               vertex_index,
    int               neighbours[],
    int               max_neighbours,
    BOOLEAN           *interior_point )
{
    int      p, n_neighbours, current_poly, current_index_within_poly;
    int      size, neighbour_index_within_poly;
    BOOLEAN  found;

    size = GET_OBJECT_SIZE( *polygons, poly );

    current_poly = poly;
    current_index_within_poly = vertex_index;
    neighbour_index_within_poly = (vertex_index-1+size)%size;

    if( max_neighbours > 0 )
    {
        neighbours[0] = polygons->indices[POINT_INDEX(polygons->end_indices,
                                          poly,neighbour_index_within_poly)];
    }

    n_neighbours = 1;

    do
    {
        found = find_next_edge_around_point( polygons,
                        current_poly, current_index_within_poly,
                        neighbour_index_within_poly,
                        &current_poly, &current_index_within_poly,
                        &neighbour_index_within_poly );

        if( found && current_poly != poly )
        {
            if( n_neighbours < max_neighbours )
            {
                neighbours[n_neighbours] = polygons->indices[
                       POINT_INDEX(polygons->end_indices,current_poly,
                                   neighbour_index_within_poly)];
            }

            ++n_neighbours;
        }
    }
    while( found && current_poly != poly );

    if( !found )
    {
        current_poly = poly;
        current_index_within_poly = vertex_index;
        neighbour_index_within_poly = (vertex_index+1+size)%size;

        if( n_neighbours < max_neighbours )
        {
            for_down( p, n_neighbours, 1 )
                neighbours[p] = neighbours[p-1];
            neighbours[0] = polygons->indices[
                            POINT_INDEX(polygons->end_indices,
                                        poly,neighbour_index_within_poly)];
        }
        ++n_neighbours;

        do
        {
            found = find_next_edge_around_point( polygons,
                            current_poly, current_index_within_poly,
                            neighbour_index_within_poly,
                            &current_poly, &current_index_within_poly,
                            &neighbour_index_within_poly );

            if( found && current_poly != poly )
            {
                if( n_neighbours < max_neighbours )
                {
                    for_down( p, n_neighbours, 1 )
                        neighbours[p] = neighbours[p-1];

                    neighbours[0] = polygons->indices[
                           POINT_INDEX(polygons->end_indices,current_poly,
                                       neighbour_index_within_poly)];
                }

                ++n_neighbours;
            }
        }
        while( found && current_poly != poly );

        if( current_poly == poly )
            print_error( "get_neighbours_of_point: topology_error" );
    }

    *interior_point = found;

    return( n_neighbours );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : get_polygons_around_vertex
@INPUT      : polygons
              poly
              vertex_index
              n_polys_alloced - max space for poly_indices
@OUTPUT     : poly_indices
              closed_flag
@RETURNS    : number of polygons
@DESCRIPTION: Returns a list of polygons, in order, around the vertex.  If
              the vertex is surrounded by polygons, closed_flag will be set
              to true.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  int  get_polygons_around_vertex(
    polygons_struct   *polygons,
    int               poly,
    int               vertex_index,
    int               poly_indices[],
    int               n_polys_alloced,
    BOOLEAN           *closed_flag )
{
    int      current_poly, current_index_within_poly;
    int      size, neighbour_index_within_poly;
    int      n_polys, dir;
    BOOLEAN  found;

    size = GET_OBJECT_SIZE( *polygons, poly );

    poly_indices[0] = poly;
    n_polys = 1;

    for( dir = -1;  dir <= 1;  dir +=2 )
    {
        current_poly = poly;
        current_index_within_poly = vertex_index;
        neighbour_index_within_poly = (vertex_index + size + dir) % size;

        do
        {
            found = find_next_edge_around_point( polygons,
                            current_poly, current_index_within_poly,
                            neighbour_index_within_poly,
                            &current_poly, &current_index_within_poly,
                            &neighbour_index_within_poly );

            if( found && current_poly != poly )
            {
                if( n_polys < n_polys_alloced )
                {
                    poly_indices[n_polys] = current_poly;
                    ++n_polys;
                }
            }
        }
        while( found && current_poly != poly );

        if( found ) break;
    }

    if( dir == -1 )
        *closed_flag = TRUE;
    else
        *closed_flag = FALSE;

    return( n_polys );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : compute_polygon_normal
@INPUT      : polygons
              poly
@OUTPUT     : normal
@RETURNS    : 
@DESCRIPTION: Computes the polygon normal.
              Vertices after the 1000th are silently ignored.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  void  compute_polygon_normal(
    polygons_struct  *polygons,
    int              poly,
    Vector           *normal )
{
#define  MAX_TEMP_STORAGE  1000
    int                e, size, point_index;
    Point              polygon[MAX_TEMP_STORAGE];

    size = GET_OBJECT_SIZE( *polygons, poly );
    if( size > MAX_TEMP_STORAGE )
        size = MAX_TEMP_STORAGE;

    for_less( e, 0, size )
    {
        point_index =
              polygons->indices[POINT_INDEX(polygons->end_indices,poly,e)];

        polygon[e] = polygons->points[point_index];
    }

    find_polygon_normal( size, polygon, normal );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : compute_polygon_normals
@INPUT      : polygons
@OUTPUT     : polygons
@RETURNS    : 
@DESCRIPTION: Computes the normals at the vertices of the polygons.  Each is
              the average of the normals of the polygonal faces
              touching the vertex, weighted by the interior angle of the
	      polygon at the vertex.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  void  compute_polygon_normals(
    polygons_struct  *polygons )
{
    int                e, poly, size, point_index, prev_index, next_index;
    Real               scale;
    Vector             normal, normal_scaled;
    progress_struct    progress;

    for_less( point_index, 0, polygons->n_points )
        fill_Vector( polygons->normals[point_index], 0.0, 0.0, 0.0 );

    initialize_progress_report( &progress, FALSE, polygons->n_items,
                                "Computing Normals" );

    for_less( poly, 0, polygons->n_items )
    {
        compute_polygon_normal( polygons, poly, &normal );

        NORMALIZE_VECTOR( normal, normal );

        size = GET_OBJECT_SIZE( *polygons, poly );

        point_index =
             polygons->indices[POINT_INDEX(polygons->end_indices,poly,size-1)];
        next_index =
             polygons->indices[POINT_INDEX(polygons->end_indices,poly,0)];

        for_less( e, 0, size )
        {
            prev_index = point_index;
            point_index = next_index;
            next_index =
                  polygons->indices[POINT_INDEX(polygons->end_indices,poly,
                                                (e+1)%size)];

            /*--- weight this contribution by the size of the angle around
                  the vertex */

            scale = get_angle_between_points( &polygons->points[prev_index],
                                              &polygons->points[point_index],
                                              &polygons->points[next_index] );

            scale = FABS( scale );
            if( scale > PI )
                scale = scale - PI;

            SCALE_VECTOR( normal_scaled, normal, scale );

            ADD_VECTORS( polygons->normals[point_index],
                         polygons->normals[point_index],
                         normal_scaled );
        }

        update_progress_report( &progress, poly + 1 );
    }

    terminate_progress_report( &progress );

    for_less( point_index, 0, polygons->n_points )
    {
        NORMALIZE_VECTOR( polygons->normals[point_index],
                          polygons->normals[point_index] );
    }
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : average_polygon_normals
@INPUT      : polygons
              n_iters
              neighbour_weight
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: Averages the normals of the polygons by iteratively setting each
              normal to some interpolation between its current value and the
              average of its neighbours.  The interpolation is controlled by
              neighbour_weight which should be between 0 and 1.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  void  average_polygon_normals(
    polygons_struct  *polygons,
    int              n_iters,
    Real             neighbour_weight )
{
    Real               avg_dot_prod;
    int                e, poly, size, point_index, neigh_index, i;
    Vector             *neigh_normal_sum, *new_normals, new_normal;
    progress_struct    progress;

    if( polygons->n_points <= 0 || polygons->n_items <= 0 )
        return;

    compute_polygon_normals( polygons );

    if( n_iters <= 0 )
        return;

    ALLOC( new_normals, polygons->n_points );
    for_less( point_index, 0, polygons->n_points )
         new_normals[point_index] = polygons->normals[point_index];

    ALLOC( neigh_normal_sum, polygons->n_points );

    for_less( i, 0, n_iters )
    {
        for_less( point_index, 0, polygons->n_points )
            fill_Vector( neigh_normal_sum[point_index], 0.0, 0.0, 0.0 );

        initialize_progress_report( &progress, FALSE, polygons->n_items,
                                    "Averaging Normals" );

        for_less( poly, 0, polygons->n_items )
        {
            size = GET_OBJECT_SIZE( *polygons, poly );

            for_less( e, 0, size )
            {
                point_index = polygons->indices[
                           POINT_INDEX(polygons->end_indices,poly,e)];
                neigh_index = polygons->indices[
                           POINT_INDEX(polygons->end_indices,poly,(e+1)%size)];

                ADD_VECTORS( neigh_normal_sum[point_index],
                             neigh_normal_sum[point_index],
                             new_normals[neigh_index] );
                ADD_VECTORS( neigh_normal_sum[neigh_index],
                             neigh_normal_sum[neigh_index],
                             new_normals[point_index] );
            }

            update_progress_report( &progress, poly + 1 );
        }

        terminate_progress_report( &progress );

        avg_dot_prod = 0.0;

        for_less( point_index, 0, polygons->n_points )
        {
            NORMALIZE_VECTOR( neigh_normal_sum[point_index],
                              neigh_normal_sum[point_index] );

            INTERPOLATE_VECTORS( new_normal,
                                 polygons->normals[point_index],
                                 neigh_normal_sum[point_index],
                                 neighbour_weight );

            NORMALIZE_VECTOR( new_normal, new_normal );

            avg_dot_prod += DOT_VECTORS( new_normal, new_normals[point_index] );
            new_normals[point_index] = new_normal;
        }

        avg_dot_prod /= (Real) polygons->n_points;

        print( "Average dot product: %g\n", avg_dot_prod );
    }

    for_less( point_index, 0, polygons->n_points )
         polygons->normals[point_index] = new_normals[point_index];

    FREE( neigh_normal_sum );
    FREE( new_normals );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : get_plane_polygon_intersection
@INPUT      : normal
              d
              polygons
              poly
@OUTPUT     : intersection_points
@RETURNS    : number of intersection points
@DESCRIPTION: Computes the intersection of a plane with a polygon, returning
              0 or 2 points, (or more if non-convex polygon).
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  BOOLEAN  get_plane_polygon_intersection(
    Vector           *normal,
    Real             d,
    polygons_struct  *polygons,
    int              poly,
    Point            intersection_points[] )
{
    int       i1, i2, edge, size, n_intersections;

    n_intersections = 0;

    size = GET_OBJECT_SIZE( *polygons, poly );

    for_less( edge, 0, size )
    {
        i1 = polygons->indices[POINT_INDEX(polygons->end_indices,poly,edge)];
        i2 = polygons->indices[
                  POINT_INDEX(polygons->end_indices,poly,(edge+1)%size)];

        if( get_plane_segment_intersection( normal, d, &polygons->points[i1],
                                        &polygons->points[i2],
                                        &intersection_points[n_intersections] ))
        {
            ++n_intersections;

            if( n_intersections == 2 )
                break;
        }
    }

    return( n_intersections == 2 );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : get_plane_segment_intersection
@INPUT      : normal
              d
              p1
              p2
@OUTPUT     : intersection_point
@RETURNS    : TRUE if intersects
@DESCRIPTION: Tests for a plane and segment intersection.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  BOOLEAN  get_plane_segment_intersection(
    Vector           *normal,
    Real             d,
    Point            *p1,
    Point            *p2,
    Point            *intersection_point )
{
    Real     dist1, dist2, t;
    BOOLEAN  intersects;

    dist1 = DIST_FROM_PLANE( *normal, d, *p1 );
    dist2 = DIST_FROM_PLANE( *normal, d, *p2 );

    if( dist1 == 0.0 )
        t = 0.0;
    else if( dist2 == 0.0 )
        t = 1.0;
    else if( dist1 == dist2 )
        t = -1.0;
    else
        t = dist1 / (dist1 - dist2);

    if( t >= 0.0 && t <= 1.0 )
    {
        INTERPOLATE_POINTS( *intersection_point, *p1, *p2, t );

        intersects = TRUE;
    }
    else
        intersects = FALSE;

    return( intersects );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : reverse_polygons_vertices
@INPUT      : polygons
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: Reverses the order of the polygon vertices.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  void   reverse_polygons_vertices(
    polygons_struct  *polygons )
{
    int         poly;

    if( polygons->neighbours != NULL )
        FREE( polygons->neighbours );

    for_less( poly, 0, polygons->n_items )
        reverse_polygon_order( polygons, poly );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : make_polygons_front_facing
@INPUT      : polygons
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: Makes all polygons front facing, based on the normals at the
              vertices.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  void   make_polygons_front_facing(
    polygons_struct  *polygons )
{
    int         poly;

    if( polygons->neighbours != (int *) 0 )
        FREE( polygons->neighbours );

    for_less( poly, 0, polygons->n_items )
    {
        if( polygon_is_back_facing( polygons, poly ) )
            reverse_polygon_order( polygons, poly );
    }
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : polygon_is_back_facing
@INPUT      : polygons
              poly
@OUTPUT     : 
@RETURNS    : TRUE if polygon is back facing
@DESCRIPTION: Uses the average normal of the polygon vertices to see if the
              polygon is back facing.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  BOOLEAN  polygon_is_back_facing(
    polygons_struct   *polygons,
    int               poly )
{
    int      i, size, point_index;
    Vector   avg_vertex_normal, polygon_normal;

    compute_polygon_normal( polygons, poly, &polygon_normal );

    size = GET_OBJECT_SIZE( *polygons, poly );

    fill_Vector( avg_vertex_normal, 0.0, 0.0, 0.0 );

    for_less( i, 0, size )
    {
        point_index = polygons->indices[
                          POINT_INDEX(polygons->end_indices,poly,i)];
        ADD_VECTORS( avg_vertex_normal, avg_vertex_normal,
                     polygons->normals[point_index] );
    }

    return( DOT_VECTORS( avg_vertex_normal, polygon_normal ) > 0.0 );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : reverse_polygon_order
@INPUT      : polygons
              poly
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: Reverses the order of the specified poly.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

static  void  reverse_polygon_order(
    polygons_struct   *polygons,
    int               poly )
{
    int   i, size, start, tmp_swap;

    size = GET_OBJECT_SIZE( *polygons, poly );

    start = POINT_INDEX( polygons->end_indices, poly, 0 );

    for_less( i, 0, size / 2 )
    {
        tmp_swap = polygons->indices[start+(size-1 - i)];
        polygons->indices[start+(size-1 - i)] = polygons->indices[start+i];
        polygons->indices[start+i] = tmp_swap;
    }
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : compute_points_centroid_and_normal
@INPUT      : polygons
              point_index
              n_neighbours
              neighbours
@OUTPUT     : centroid
              normal
              base_length
              curvature
@RETURNS    : 
@DESCRIPTION: Computes the centroid and normal of the neighbours of a vertex,
              as well as a measure of the size of the polygon defined by the
              neighbours, and the relative curvature of the surface at the
              vertex.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  void  compute_points_centroid_and_normal(
    polygons_struct  *polygons,
    int              point_index,
    int              n_neighbours,
    int              neighbours[],
    Point            *centroid,
    Vector           *normal,
    Real             *base_length,
    Real             *curvature )
{
#define  MAX_NEIGHBOURS   1000
    int              i;
    Point            neigh_points[MAX_NEIGHBOURS];

    if( n_neighbours > 2 )
    {
        for_less( i, 0, n_neighbours )
            neigh_points[i] = polygons->points[neighbours[i]];

        get_points_centroid( n_neighbours, neigh_points, centroid );

        find_polygon_normal( n_neighbours, neigh_points, normal );

        *curvature = estimate_polygon_curvature( &polygons->points[point_index],
                                       n_neighbours, neigh_points,
                                       centroid, normal, base_length );
    }
    else
    {
        *centroid = polygons->points[point_index];
        fill_Vector( *normal, 0.0, 0.0, 0.0 );
        *base_length = 1.0;
        *curvature = 0.0;
    }
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : compute_points_centroid_and_normal
@INPUT      : polygons
              poly
              vertex_index
              point_index
@OUTPUT     : centroid
              normal
              base_length
              curvature
@RETURNS    : 
@DESCRIPTION: Computes the centroid and normal of the neighbours of a vertex,
              as well as a measure of the size of the polygon defined by the
              neighbours, and the relative curvature of the surface at the
              vertex.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  void  compute_polygon_point_centroid(
    polygons_struct  *polygons,
    int              poly,
    int              vertex_index,
    int              point_index,
    Point            *centroid,
    Vector           *normal,
    Real             *base_length,
    Real             *curvature )
{
    int              n_neighbours;
    int              neighbours[MAX_NEIGHBOURS];
    BOOLEAN          interior_point;

    n_neighbours = get_neighbours_of_point( polygons, poly, vertex_index,
                                            neighbours, MAX_NEIGHBOURS,
                                            &interior_point );

    compute_points_centroid_and_normal( polygons, point_index,
                                        n_neighbours, neighbours,
                                        centroid, normal, base_length,
                                        curvature );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : estimate_polygon_curvature
@INPUT      : point
              n_neighs
              neighs
              centroid
              normal
@OUTPUT     : base_length
@RETURNS    : curvature
@DESCRIPTION: Computes the curvature at a vertex of a polygon.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

static  Real  estimate_polygon_curvature(
    Point   *point,
    int     n_neighs,
    Point   neighs[],
    Point   *centroid,
    Vector  *normal,
    Real    *base_length )
{
    int    i;
    Real   curvature, len;
    Vector to_point;

    len = 0.0;
    for_less( i, 0, n_neighs )
        len += distance_between_points( &neighs[i], centroid );

    if( n_neighs > 0 )
        len = 2.0 * len / (Real) n_neighs;

    if( len == 0.0 )
        len = 1.0;

    *base_length = len;

    SUB_POINTS( to_point, *point, *centroid );

    curvature = MAGNITUDE( to_point ) / len;

    if( DOT_VECTORS( to_point, *normal ) < 0.0 )
        curvature = -curvature;

    return( curvature );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : compute_polygon_vertex_curvature
@INPUT      : polygons
              point_index
@OUTPUT     : 
@RETURNS    : curvature
@DESCRIPTION: Computes the curvature of the surface at a specified vertex.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  Real  compute_polygon_vertex_curvature(
    polygons_struct  *polygons,
    int              point_index )
{
    Real      curvature, base_length;
    int       poly, vertex;
    Point     centroid;
    Vector    normal;

    if( !find_polygon_with_vertex( polygons, point_index, &poly, &vertex ) )
    {
        handle_internal_error( "compute_polygon_vertex_curvature" );
        return( 0.0 );
    }

    compute_polygon_point_centroid( polygons, poly, vertex, point_index,
                                    &centroid, &normal, &base_length,
                                    &curvature );

    return( curvature );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : get_opposite_point
@INPUT      : polygons
              poly
              edge
@OUTPUT     : point
@RETURNS    : 
@DESCRIPTION: Gets the point opposite the given edge.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

static  void  get_opposite_point(
    polygons_struct  *polygons,
    int              poly,
    int              edge,
    Point            *point )
{
    int  v, size;

    size = GET_OBJECT_SIZE( *polygons, poly );

    if( size == 3 )
        v = (edge + 2) % size;
    else
        v = (edge + size/2) % size;

    *point = polygons->points[polygons->indices[
                         POINT_INDEX( polygons->end_indices, poly, v)]];
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : get_polygon_edge_angle
@INPUT      : polygons
              poly
              edge
@OUTPUT     : 
@RETURNS    : angle between -PI and PI
@DESCRIPTION: Computes the angle across the edge.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  Real  get_polygon_edge_angle(
    polygons_struct  *polygons,
    int              poly,
    int              edge )
{
    int     size, i, point_index1, point_index2, neighbour_poly;
    Real    angle, edge_len_squared, scale, x, y;
    Point   p1, p2, poly1_point, poly2_point;
    Vector  v1, v2, normal, diff, edge_vec;

    neighbour_poly = polygons->neighbours[POINT_INDEX( polygons->end_indices,
                                                       poly, edge )];

    if( neighbour_poly < 0 )
        return( PI );

    size = GET_OBJECT_SIZE( *polygons, poly );

    point_index1 = polygons->indices[ POINT_INDEX( polygons->end_indices,
                                                   poly, edge )];
    point_index2 = polygons->indices[ POINT_INDEX( polygons->end_indices,
                                                   poly, (edge+1)%size )];
    p1 = polygons->points[point_index1];
    p2 = polygons->points[point_index2];

    get_opposite_point( polygons, poly, edge, &poly1_point );

    i = find_edge_index( polygons, neighbour_poly, point_index1, point_index2 );
    get_opposite_point( polygons, neighbour_poly, i, &poly2_point );

    SUB_POINTS( edge_vec, p2, p1 );
    SUB_POINTS( v1, poly1_point, p1 );
    SUB_POINTS( v2, poly2_point, p1 );

    edge_len_squared = DOT_VECTORS( edge_vec, edge_vec );
    if( edge_len_squared == 0.0 )
        edge_len_squared = 1.0;

    scale = DOT_VECTORS( v1, edge_vec ) / edge_len_squared;
    SCALE_VECTOR( diff, edge_vec, scale );
    SUB_VECTORS( v1, v1, diff );
    NORMALIZE_VECTOR( v1, v1 );

    CROSS_VECTORS( normal, edge_vec, v1 );
    NORMALIZE_VECTOR( normal, normal );

    scale = DOT_VECTORS( v2, edge_vec ) / edge_len_squared;
    SCALE_VECTOR( diff, edge_vec, scale );
    SUB_VECTORS( v2, v2, diff );

    x = DOT_VECTORS( v2, v1 );
    y = -DOT_VECTORS( v2, normal );

    angle = compute_clockwise_rotation( x, y );

    if( angle < 0.0 )
        angle += 2.0 * PI;

    return( angle );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : polygons_are_same_topology
@INPUT      : p1
              p2
@OUTPUT     : 
@RETURNS    : TRUE if same
@DESCRIPTION: Tests if two polygons are the same topology.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  BOOLEAN  polygons_are_same_topology(
    polygons_struct  *p1,
    polygons_struct  *p2 )
{
    return( objects_are_same_topology( p1->n_points, p1->n_items,
                                       p1->end_indices, p1->indices,
                                       p2->n_points, p2->n_items,
                                       p2->end_indices, p2->indices ) );
}


/*!
 * @}
 * End of grp_polygons.
 */
