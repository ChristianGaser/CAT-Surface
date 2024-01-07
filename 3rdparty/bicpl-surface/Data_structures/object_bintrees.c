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
static char rcsid[] = "$Header: /private-cvsroot/libraries/bicpl/Data_structures/object_bintrees.c,v 1.13 2005/08/17 22:31:12 bert Exp $";
#endif

/* ----------------------------- MNI Header -----------------------------------
@NAME       : delete_the_bintree
@INPUT      : bintree
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: Deletes the bintree.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : 1993            David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  void  delete_the_bintree(
    bintree_struct_ptr  *bintree )
{
    if( *bintree != (bintree_struct_ptr) NULL )
    {
        delete_bintree( *bintree );

        FREE( *bintree );
        *bintree = (bintree_struct_ptr) NULL;
    }
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : check_install_bintree_delete_function
@INPUT      : 
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: Checks whether the bintree delete function has been installed.
              In order to avoid linking in all the bintree stuff whenever
              object_struct's are used, the delete function is called as
              a pointer.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : 1993            David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

static  void  check_install_bintree_delete_function( void )
{
    static  BOOLEAN  first = TRUE;

    if( first )
    {
        first = FALSE;
        set_bintree_delete_function( delete_the_bintree );
    }
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : allocate_bintree
@INPUT      : 
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: Allocates a bintree structure.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : 1993            David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  bintree_struct_ptr allocate_bintree( void )
{
    bintree_struct_ptr   bintree;

    ALLOC( bintree, 1 );

    return( bintree );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : create_lines_bintree
@INPUT      : lines
              max_nodes
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: Creates a bintree for the lines, storing it in lines->bintree.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : 1993            David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  void  create_lines_bintree(
    lines_struct   *lines,
    int            max_nodes )
{
    Real             radius;
    int              line, size, n_segments, seg, object_id;
    range_struct     *bound_vols;
    Point            min_range, max_range;
    Point            points[2];

    check_install_bintree_delete_function();

    lines->bintree = allocate_bintree();

    n_segments = 0;
    for_less( line, 0, lines->n_items )
        n_segments += GET_OBJECT_SIZE( *lines, line ) - 1;

    ALLOC( bound_vols, n_segments );

    radius = (Real) lines->line_thickness;

    object_id = 0;
    for_less( line, 0, lines->n_items )
    {
        size = GET_OBJECT_SIZE( *lines, line );

        for_less( seg, 0, size - 1 )
        {
            points[0] = lines->points[lines->indices[
                          POINT_INDEX(lines->end_indices,line,seg)]];
            points[1] = lines->points[lines->indices[
                          POINT_INDEX(lines->end_indices,line,seg+1)]];

            get_range_points( 2, points, &min_range, &max_range );
            bound_vols[object_id].limits[X][0] =
                                 (float) ((Real) Point_x(min_range) - radius);
            bound_vols[object_id].limits[Y][0] =
                                 (float) ((Real) Point_y(min_range) - radius);
            bound_vols[object_id].limits[Z][0] =
                                 (float) ((Real) Point_z(min_range) - radius);
            bound_vols[object_id].limits[X][1] =
                                 (float) ((Real) Point_x(max_range) + radius);
            bound_vols[object_id].limits[Y][1] =
                                 (float) ((Real) Point_y(max_range) + radius);
            bound_vols[object_id].limits[Z][1] =
                                 (float) ((Real) Point_z(max_range) + radius);
            ++object_id;
        }
    }

    create_object_bintree( n_segments, bound_vols,
                           lines->bintree, max_nodes );

    FREE( bound_vols );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : create_polygons_bintree
@INPUT      : polygons
              max_nodes
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: Creates a bintree for the polygons, storing it in
              polygons->bintree.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : 1993            David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  void  create_polygons_bintree(
    polygons_struct   *polygons,
    int               max_nodes )
{
    int              poly, size;
    range_struct     *bound_vols;
    Point            min_range, max_range;
    Point            points[MAX_POINTS_PER_POLYGON];

    check_install_bintree_delete_function();

    polygons->bintree = allocate_bintree();

    ALLOC( bound_vols, polygons->n_items );

    for_less( poly, 0, polygons->n_items )
    {
        size = get_polygon_points( polygons, poly, points );

        get_range_points( size, points, &min_range, &max_range );
        bound_vols[poly].limits[X][0] = Point_x(min_range);
        bound_vols[poly].limits[Y][0] = Point_y(min_range);
        bound_vols[poly].limits[Z][0] = Point_z(min_range);
        bound_vols[poly].limits[X][1] = Point_x(max_range);
        bound_vols[poly].limits[Y][1] = Point_y(max_range);
        bound_vols[poly].limits[Z][1] = Point_z(max_range);
    }

    create_object_bintree( polygons->n_items, bound_vols,
                           polygons->bintree, max_nodes );

    FREE( bound_vols );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : create_quadmesh_bintree
@INPUT      : quadmesh
              max_nodes
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: Creates a bintree for the quadmesh, storing it in
              quadmesh->bintree.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : 1993            David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  void  create_quadmesh_bintree(
    quadmesh_struct   *quadmesh,
    int               max_nodes )
{
    int              i, j, m, n, obj_index;
    range_struct     *bound_vols;
    Point            min_range, max_range;
    Point            points[4];

    check_install_bintree_delete_function();

    quadmesh->bintree = allocate_bintree();

    get_quadmesh_n_objects( quadmesh, &m, &n );

    ALLOC( bound_vols, m * n );

    for_less( i, 0, m )
    {
        for_less( j, 0, n )
        {
            obj_index = IJ( i, j, n );
            get_quadmesh_patch( quadmesh, i, j, points );

            get_range_points( 4, points, &min_range, &max_range );

            bound_vols[obj_index].limits[X][0] = Point_x(min_range);
            bound_vols[obj_index].limits[Y][0] = Point_y(min_range);
            bound_vols[obj_index].limits[Z][0] = Point_z(min_range);
            bound_vols[obj_index].limits[X][1] = Point_x(max_range);
            bound_vols[obj_index].limits[Y][1] = Point_y(max_range);
            bound_vols[obj_index].limits[Z][1] = Point_z(max_range);
        }
    }

    create_object_bintree( m * n, bound_vols,
                           quadmesh->bintree, max_nodes );

    FREE( bound_vols );
}
