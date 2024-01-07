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
static char rcsid[] = "$Header: /private-cvsroot/libraries/bicpl/Geometry/tubes.c,v 1.7 2005/08/17 22:30:25 bert Exp $";
#endif

static  void  generate_tube(
    int    n_points,
    Point  points[],
    int    n_around,
    Real   radius,
    Point  tube_points[],
    Vector tube_normals[] );

/* ----------------------------- MNI Header -----------------------------------
@NAME       : create_tube
@INPUT      : n_points
              points
              n_around
              radius
              colour
@OUTPUT     : quadmesh
@RETURNS    : 
@DESCRIPTION: Creates a quadmesh cylindrical tube corresponding to the points.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

static  void  create_tube(
    int               n_points,
    Point             points[],
    int               n_around,
    Real              radius,
    Colour            colour,
    quadmesh_struct   *quadmesh )
{
    Surfprop    spr;

    get_default_surfprop( &spr );

    initialize_quadmesh( quadmesh, colour, &spr, n_points, n_around );
    quadmesh->m_closed = FALSE;
    quadmesh->n_closed = TRUE;

    generate_tube( n_points, points, n_around, radius,
                   quadmesh->points, quadmesh->normals );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : convert_lines_to_tubes
@INPUT      : lines
              n_around
              radius
@OUTPUT     : quadmeshes
@RETURNS    : number of quadmeshes
@DESCRIPTION: Creates a set of quadmeshes representing tubes corresponding
              to the lines structure, returning the number of quadmeshes
              (the number of separate lines).
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  int  convert_lines_to_tubes(
    lines_struct     *lines,
    int              n_around,
    Real             radius,
    quadmesh_struct  *quadmeshes[] )
{
    int              i, l, line_size;
    quadmesh_struct  *quadmesh;
    Point            *points;

    if( lines->n_items == 0 )
        return( 0 );

    ALLOC( *quadmeshes, lines->n_items );

    for_less( l, 0, lines->n_items )
    {
        line_size = GET_OBJECT_SIZE( *lines, l );

        ALLOC( points, line_size );

        for_less( i, 0, line_size )
        {
            points[i] = lines->points[lines->indices[
                           POINT_INDEX(lines->end_indices,l,i)]];
        }

        quadmesh = &(*quadmeshes)[l];

        create_tube( line_size, points, n_around, radius, lines->colours[0],
                     quadmesh );

        FREE( points );
    }

    return( lines->n_items );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : get_direction
@INPUT      : n_points
              points
              i
              wrap_around
@OUTPUT     : dir
@RETURNS    : 
@DESCRIPTION: Gets the direction of the line segment from i to i+1 of the
              points, with or without wraparound.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

static  void  get_direction(
    int      n_points,
    Point    points[],
    int      i,
    BOOLEAN  wrap_around,
    Vector   *dir )
{
    int      end;

    end = i;

    do
    {
        ++end;
        if( end == n_points && wrap_around )
            end = 0;
    }
    while( end < n_points && end != i && EQUAL_POINTS(points[i],points[end]) );

    if( end >= n_points || end == i )
    {
        fill_Vector( *dir, 0.0, 0.0, 0.0 );
    }
    else
    {
        SUB_POINTS( *dir, points[end], points[i] );

        NORMALIZE_VECTOR( *dir, *dir );
    }
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : fill_in_ellipse_points
@INPUT      : centre
              n_around
              x
              y
              hor
              vert
@OUTPUT     : tube_points
              tube_normals
@RETURNS    : 
@DESCRIPTION: Given x,y positions in 2D space and a horizontal and vertical
              vector, fills in the 3D ellipse points into the tube points
              and normals
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

static  void   fill_in_ellipse_points(
    Point   tube_points[],
    Vector  tube_normals[],
    Point   *centre,
    int     n_around,
    Real    x[],
    Real    y[],
    Vector  *hor,
    Vector  *vert )
{
    int      i;
    Vector   h, v, offset;

    for_less( i, 0, n_around )
    {
        SCALE_VECTOR( h, *hor, x[i] );
        SCALE_VECTOR( v, *vert, y[i] );

        ADD_VECTORS( offset, h, v );

        ADD_POINT_VECTOR( tube_points[i], *centre, offset );
        NORMALIZE_VECTOR( tube_normals[i], offset );
    }
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : project_vector_to_plane
@INPUT      : v
              direction
              normal
@OUTPUT     : projected
@RETURNS    : 
@DESCRIPTION: Projects a vector onto a plane, along the direction.  Finds the
              point on the line through v in the direction direction, which
              is on the plane perpendicular to normal.  Solves the equation 
              (V + t * Direction) dot Normal = 0
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

static  void  project_vector_to_plane(
    Vector   *v,
    Vector   *direction,
    Vector   *normal,
    Vector   *projected )
{
    Vector   offset;
    Real     t, n_dot_d, n_dot_v;

    n_dot_d = DOT_VECTORS( *normal, *direction );

    if( n_dot_d == 0.0 )
    {
        print_error( "Error in project_vector_to_plane\n" );
    }
    else
    {
        n_dot_v = DOT_VECTORS( *normal, *v );

        t =  - n_dot_v / n_dot_d;

        SCALE_VECTOR( offset, *direction, t );
        ADD_VECTORS( *projected, *v, offset );
    }
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : generate_tube
@INPUT      : n_points
              points
              n_around
              radius
@OUTPUT     : tube_points
              tube_normals
@RETURNS    : 
@DESCRIPTION: Generates a tube quadmesh from a set of points.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

static  void  generate_tube(
    int    n_points,
    Point  points[],
    int    n_around,
    Real   radius,
    Point  tube_points[],
    Vector tube_normals[] )
{
    int      i;
    Real     *x, *y, angle;
    Vector   hor, vert, dir, prev_dir, normal;
    BOOLEAN  wrap_around;

    ALLOC( x, n_around );
    ALLOC( y, n_around );

    for_less( i, 0, n_around )
    {
        angle = (Real) (n_around - 1 - i) / (Real) n_around * 2.0 * PI;
        x[i] = radius * cos( angle );
        y[i] = radius * sin( angle );
    }

    wrap_around = EQUAL_POINTS( points[0], points[n_points-1] );

    get_direction( n_points, points, 0, wrap_around, &dir );
    create_orthogonal_vector( &dir, &hor );

    if( wrap_around )
        get_direction( n_points, points, n_points-2, wrap_around, &prev_dir );
    else
        prev_dir = dir;

    for_less( i, 0, n_points )
    {
        CROSS_VECTORS( vert, prev_dir, hor );
        NORMALIZE_VECTOR( vert, vert );

        CROSS_VECTORS( hor, vert, prev_dir );
        NORMALIZE_VECTOR( hor, hor );

        ADD_VECTORS( normal, prev_dir, dir );
        NORMALIZE_VECTOR( normal, normal );

        project_vector_to_plane( &hor, &prev_dir, &normal, &hor );
        project_vector_to_plane( &vert, &prev_dir, &normal, &vert );

        fill_in_ellipse_points( &tube_points[IJ(i,0,n_around)],
                                &tube_normals[IJ(i,0,n_around)],
                                &points[i], n_around, x, y, &hor, &vert );

        if( i < n_points-1 )
        {
            prev_dir = dir;

            if( wrap_around || i < n_points-2 )
                get_direction( n_points, points, i+1, wrap_around, &dir );
        }
    }

    FREE( x );
    FREE( y );
}
