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
static char rcsid[] = "$Header: /private-cvsroot/libraries/bicpl/Objects/lines.c,v 1.20 2005/08/17 22:28:26 bert Exp $";
#endif

/* ----------------------------- MNI Header -----------------------------------
@NAME       : initialize_lines
@INPUT      : col
@OUTPUT     : lines
@RETURNS    : 
@DESCRIPTION: Initializes the lines to empty.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  void  initialize_lines(
    lines_struct    *lines,
    Colour          col )
{
    ALLOC( lines->colours, 1 );

    lines->colour_flag = ONE_COLOUR;
    lines->colours[0] = col;

    lines->line_thickness = 1.0f;
    lines->n_points = 0;
    lines->n_items = 0;

    lines->bintree = (bintree_struct_ptr) NULL;
}

BICAPI  void  initialize_lines_with_size(
    lines_struct    *lines,
    Colour          col,
    int             size,
    BOOLEAN         closed )
{
    int  i, n_indices;

    initialize_lines( lines, col );

    if( closed )
        n_indices = size + 1;
    else
        n_indices = size;

    lines->n_points = size;
    ALLOC( lines->points, size );

    lines->n_items = 1;
    ALLOC( lines->end_indices, 1 );
    lines->end_indices[0] = n_indices;

    ALLOC( lines->indices, n_indices );

    for_less( i, 0, n_indices )
        lines->indices[i] = i % size;
}


/* ----------------------------- MNI Header -----------------------------------
@NAME       : delete_lines
@INPUT      : lines
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: Frees up memory associated with the lines.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  void  delete_lines( lines_struct *lines )
{
    free_colours( lines->colour_flag, lines->colours, lines->n_points,
                  lines->n_items );

    if( lines->n_points > 0 )
        FREE( lines->points );

    if( lines->n_items > 0 )
        FREE( lines->end_indices );

    if( lines->n_items > 0 )
        FREE( lines->indices );

    delete_bintree_if_any( &lines->bintree );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : start_new_line
@INPUT      : lines
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: Creates a new line in the structure of size 0, in preparation
              for adding points to the line.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  void  start_new_line( lines_struct *lines )
{
    int      n_indices;

    n_indices = NUMBER_INDICES( *lines );

    ADD_ELEMENT_TO_ARRAY( lines->end_indices, lines->n_items,
                          n_indices, DEFAULT_CHUNK_SIZE );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : add_point_to_line
@INPUT      : lines
              point
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: Adds a point to the line.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  void  add_point_to_line(
    lines_struct   *lines,
    Point          *point )
{
    if( lines->n_items == 0 )
        start_new_line( lines );

    ADD_ELEMENT_TO_ARRAY( lines->indices, lines->end_indices[lines->n_items-1],
                          lines->n_points, DEFAULT_CHUNK_SIZE );

    ADD_ELEMENT_TO_ARRAY( lines->points, lines->n_points,
                          *point, DEFAULT_CHUNK_SIZE );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : get_line_segment_index
@INPUT      : lines
              obj_index
@OUTPUT     : line
              seg
@RETURNS    : 
@DESCRIPTION: Converts an obj_index to a line index, and a line segment
              index.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  void  get_line_segment_index(
    lines_struct  *lines,
    int           obj_index,
    int           *line,
    int           *seg )
{
    int    first, last, mid, n_objs;

    first = 0;
    last = lines->n_items-1;

    while( first < last )
    {
        mid = (first + last) / 2;

        n_objs = lines->end_indices[mid] - mid - 1;
        if( obj_index < n_objs )
            last = mid;
        else
            first = mid+1;
    }

    *line = first;

    if( first == 0 )
        *seg = obj_index;
    else
        *seg = obj_index - (lines->end_indices[first-1] - first);
}

static  void  (*bintree_delete_function) ( bintree_struct_ptr* ) = NULL;

/* ----------------------------- MNI Header -----------------------------------
@NAME       : set_bintree_delete_function
@INPUT      : func
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: Sets the bintree delete function for objects.  This is done
              so that if the application does not use any bintree code, the
              bintree library does not get linked into the executable, simply
              because the delete objects code might have to delete a bintree.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  void  set_bintree_delete_function(
    void  (*func)( bintree_struct_ptr * ) )
{
    bintree_delete_function = func;
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : delete_bintree_if_any
@INPUT      : bintree
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: Deletes the bintree, if non-null.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  void  delete_bintree_if_any(
    bintree_struct_ptr   *bintree )
{
    if( bintree_delete_function != NULL )
        (*bintree_delete_function) ( bintree );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : get_lines_length
@INPUT      : lines
@OUTPUT     : 
@RETURNS    : length
@DESCRIPTION: Gets the total length of the lines.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  Real  get_lines_length(
    lines_struct  *lines )
{
    Real  len;
    int   line, i, p0, p1, size;

    len = 0.0;

    for_less( line, 0, lines->n_items )
    {
        size = GET_OBJECT_SIZE( *lines, line );

        for_less( i, 0, size-1 )
        {
            p0 = lines->indices[POINT_INDEX(lines->end_indices,line,i)];
            p1 = lines->indices[POINT_INDEX(lines->end_indices,line,i+1)];
            len += distance_between_points( &lines->points[p0],
                                            &lines->points[p1] );
        }
    }

    return( len );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : 
@INPUT      : 
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: 
@METHOD     : 
@GLOBALS    : 
@CALLS      :  
@CREATED    : Jul. 16, 1997    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  void  get_lines_arc_point(
    lines_struct  *lines,
    Real          arc_length,
    Point         *point )
{
    Real     len, segment_length, ratio;
    int      line, i, p0, p1, size;
    BOOLEAN  found;

    if( arc_length < 0.0 )
    {
        print_error( "get_lines_arc_point: arc_length < 0.0, using 0.0\n" );
        arc_length = 0.0;
    }

    len = 0.0;
    found = FALSE;

    for_less( line, 0, lines->n_items )
    {
        size = GET_OBJECT_SIZE( *lines, line );

        for_less( i, 0, size-1 )
        {
            p0 = lines->indices[POINT_INDEX(lines->end_indices,line,i)];
            p1 = lines->indices[POINT_INDEX(lines->end_indices,line,i+1)];
            segment_length = distance_between_points( &lines->points[p0],
                                                      &lines->points[p1] );
            len += segment_length;

            if( len > arc_length )
            {
                found = TRUE;
                break;
            }
        }

        if( found )
            break;
    }

    if( !found )
    {
        print_error( "get_lines_arc_point: arc_length too large, using end\n" );
        *point = lines->points[p1];
        return;
    }

    ratio = (len - arc_length) / segment_length;
    INTERPOLATE_POINTS( *point, lines->points[p1], lines->points[p0], ratio );
}
