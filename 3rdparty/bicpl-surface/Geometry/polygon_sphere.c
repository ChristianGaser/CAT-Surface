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
static char rcsid[] = "$Header: /private-cvsroot/libraries/bicpl/Geometry/polygon_sphere.c,v 1.17 2005/08/17 22:30:25 bert Exp $";
#endif

static  int  get_n_sphere_points(
    int   n_up,
    int   n_around );
static  void  get_sphere_point(
    Real     up,
    Real     around,
    Point    *centre,
    Real     x_size,
    Real     y_size,
    Real     z_size,
    Point    *point );
static  void  get_subdivided_point(
    int      up,
    int      around,
    Point    input_points[],
    int      input_n_up,
    Point    *point );

/* ----------------------------- MNI Header -----------------------------------
@NAME       : create_polygons_sphere
@INPUT      : centre
              x_size
              y_size
              z_size
              n_up
              n_around
              subdividing_flag
@OUTPUT     : polygons
@RETURNS    : 
@DESCRIPTION: Creates a tessellated sphere, using lines of latitude and
              longitude to create triangles around the poles, and quadrilaterals
              everywhere else.  If subdividing_flag is true, then we are
              subdividing the polygons and are maintaining its sphere topology.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  void  create_polygons_sphere(
    Point            *centre,
    Real             x_size,
    Real             y_size,
    Real             z_size,
    int              n_up,
    int              n_around,
    BOOLEAN          subdividing_flag,
    polygons_struct  *polygons )
{
    int      point_index, top_point_index, bottom_point_index;
    int      point_index1, point_index2, point_index3, point_index4;
    int      up, around, n_circum, next_around, n_indices, end, start, a;
    int      n_around_top;
    Point    *input_points;
    Colour   save_colour;
    int      input_n_up;
    Real     up_pos, around_pos;

    if( subdividing_flag )
    {
        if( !get_tessellation_of_polygons_sphere( polygons, &input_n_up ))
        {
            print_error( "Not a sphere topology polygon.\n" );
            return;
        }

        n_up = 2 * input_n_up;
        n_around = 2 * n_up;

        save_colour = polygons->colours[0];

        input_points = polygons->points;

        /* so delete polygons frees a valid pointer */

        ALLOC( polygons->points, 1 );

        delete_polygons( polygons );
    }

    initialize_polygons( polygons, WHITE, (Surfprop *) NULL );

    if( subdividing_flag )
        polygons->colours[0] = save_colour;

    polygons->n_points = get_n_sphere_points( n_up, n_around );

    ALLOC( polygons->points, polygons->n_points );
    ALLOC( polygons->normals, polygons->n_points );

    /* ------ build points ------ */

    for_less( up, 0, n_up+1 )
    {
        up_pos = (Real) up / (Real) n_up;

        if( up == 0 || up == n_up )
            n_circum = 1;
        else
            n_circum = n_around;

        for_less( around, 0, n_circum )
        {
            around_pos = (Real) around / (Real) n_around;

            point_index = get_sphere_point_index( up, around, n_up, n_around );

            if( subdividing_flag )
            {
                get_subdivided_point( up, around, input_points, input_n_up,
                                      &polygons->points[point_index] );
            }
            else
            {
                get_sphere_point( up_pos, around_pos, centre, x_size, y_size,
                                  z_size, &polygons->points[point_index] );
            }
        }
    }

    if( subdividing_flag )
        FREE( input_points );

    n_indices = 0;

    /* ------ build indices for top ------ */

    top_point_index = get_sphere_point_index( 0, 0, n_up, n_around );

#ifdef  N_AROUND_TOP
    n_around_top = N_AROUND_TOP;
#else
    n_around_top = n_around;
#endif

    for_less( a, 0, n_around_top )
    {
        ADD_ELEMENT_TO_ARRAY( polygons->indices, n_indices,
                              top_point_index, DEFAULT_CHUNK_SIZE );

        start = n_around * a / n_around_top;
        end = n_around * (a+1) / n_around_top + 1;
        if( end > n_around + 1 )
            end = n_around + 1;
        for_less( around, start, end )
        {
            point_index = get_sphere_point_index( 1, around % n_around,
                                                  n_up, n_around );

            ADD_ELEMENT_TO_ARRAY( polygons->indices, n_indices,
                                  point_index, DEFAULT_CHUNK_SIZE );
        }

        ADD_ELEMENT_TO_ARRAY( polygons->end_indices, polygons->n_items,
                              n_indices, DEFAULT_CHUNK_SIZE );
    }

    /* ------ build indices for main part ------ */

    for_less( up, 1, n_up-1 )
    {
        for_less( around, 0, n_around )
        {
            next_around = (around + 1) % n_around;
            point_index1 = get_sphere_point_index( up, around, n_up, n_around );
            point_index2 = get_sphere_point_index( up+1, around, n_up, n_around );
            point_index3 = get_sphere_point_index( up+1, next_around, n_up, n_around );
            point_index4 = get_sphere_point_index( up, next_around, n_up, n_around );

            ADD_ELEMENT_TO_ARRAY( polygons->indices, n_indices,
                                  point_index1, DEFAULT_CHUNK_SIZE );
            ADD_ELEMENT_TO_ARRAY( polygons->indices, n_indices,
                                  point_index2, DEFAULT_CHUNK_SIZE );
            ADD_ELEMENT_TO_ARRAY( polygons->indices, n_indices,
                                  point_index3, DEFAULT_CHUNK_SIZE );
            ADD_ELEMENT_TO_ARRAY( polygons->indices, n_indices,
                                  point_index4, DEFAULT_CHUNK_SIZE );

            ADD_ELEMENT_TO_ARRAY( polygons->end_indices, polygons->n_items,
                                  n_indices, DEFAULT_CHUNK_SIZE );
        }
    }

    /* ------ build indices for bottom ------ */

    bottom_point_index = get_sphere_point_index( n_up, 0, n_up, n_around );

    for_less( a, 0, n_around_top )
    {
        ADD_ELEMENT_TO_ARRAY( polygons->indices, n_indices,
                              bottom_point_index, DEFAULT_CHUNK_SIZE );

        end = n_around * a / n_around_top;
        start = n_around * (a+1) / n_around_top;
        if( start > n_around )
            start = n_around;
        for( around = start;  around >= end;  --around )
        {
            point_index = get_sphere_point_index( n_up-1, around % n_around,
                                                  n_up, n_around );

            ADD_ELEMENT_TO_ARRAY( polygons->indices, n_indices,
                                  point_index, DEFAULT_CHUNK_SIZE );
        }

        ADD_ELEMENT_TO_ARRAY( polygons->end_indices, polygons->n_items,
                              n_indices, DEFAULT_CHUNK_SIZE );
    }
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : get_sphere_point_index
@INPUT      : up
              around
              n_up
              n_around
@OUTPUT     : 
@RETURNS    : index
@DESCRIPTION: Computes the point index in the sphere tessellation, where
              up is in range 0 to n_up, and around is in the range of 0
              to n_around-1.  up = 0 and up == n_up correspond to the 2 poles.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  int  get_sphere_point_index(
    int   up,
    int   around,
    int   n_up,
    int   n_around )
{
    int   point_index;

    if( up < 0 || up > n_up || around < 0 || around >= n_around )
    {
        print_error( "up %d/%d     around %d/%d\n", up, n_up,
                      around, n_around );
        handle_internal_error( "get_sphere_point_index" );
    }

    if( up == 0 )
    {
        point_index = 0;
    }
    else if( up == n_up )
    {
        point_index = 1;
    }
    else
    {
        point_index = 2 + IJ( up-1, around, n_around );
    }

    return( point_index );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : get_n_sphere_points
@INPUT      : n_up
              n_around
@OUTPUT     : 
@RETURNS    : number of points
@DESCRIPTION:
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

static  int  get_n_sphere_points(
    int   n_up,
    int   n_around )
{
    return( 2 + (n_up-1) * n_around );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : get_sphere_point
@INPUT      : up        : between 0 and 1
              around    : between 0 and 1
              centre
              x_size
              y_size
              z_size
@OUTPUT     : point
@RETURNS    : 
@DESCRIPTION: Computes a point on a sphere.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

static  void  get_sphere_point(
    Real     up,
    Real     around,
    Point    *centre,
    Real     x_size,
    Real     y_size,
    Real     z_size,
    Point    *point )
{
    Real      dx, dy, dz;

    dx = x_size * cos( (double) up * PI );
    dy = y_size * sin( (double) up * PI ) * cos( (double) around * 2.0 *PI);
    dz = z_size * sin( (double) up * PI ) * sin( (double) around * 2.0 *PI);

    fill_Point( *point, (Real) Point_x(*centre) + dx,
                        (Real) Point_y(*centre) + dy,
                        (Real) Point_z(*centre) + dz );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : get_subdivided_point
@INPUT      : up
              around
              input_points
              input_n_up
@OUTPUT     : point
@RETURNS    : 
@DESCRIPTION: Finds the point in the subdivided sphere topology.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

static  void  get_subdivided_point(
    int      up,
    int      around,
    Point    input_points[],
    int      input_n_up,
    Point    *point )
{
    int       input_up, input_around;
    BOOLEAN   up_midpoint_flag, around_midpoint_flag;
    Point     corner_below, corner_across, corner_below_across;

    input_up = up / 2;
    up_midpoint_flag = (up & 1) == 1;
    input_around = around / 2;
    around_midpoint_flag = (around & 1) == 1;

    *point = input_points[get_sphere_point_index( input_up,
                             input_around, input_n_up, 2 * input_n_up )];

    if( up_midpoint_flag )
    {
        corner_below = input_points[get_sphere_point_index( input_up+1,
                         input_around, input_n_up, 2 * input_n_up )];
        INTERPOLATE_POINTS( *point, *point, corner_below, 0.5 );
    }

    if( around_midpoint_flag )
    {
        corner_across = input_points[get_sphere_point_index( input_up,
           (input_around + 1) % (2*input_n_up), input_n_up, 2*input_n_up )];

        if( up_midpoint_flag )
        {
            corner_below_across = input_points[get_sphere_point_index(
               input_up+1, (input_around + 1) % (2*input_n_up),
               input_n_up, 2*input_n_up )];
            INTERPOLATE_POINTS( corner_across, corner_across,
                                corner_below_across, 0.5 );
        }

        INTERPOLATE_POINTS( *point, *point, corner_across, 0.5 );
    }
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : is_this_sphere_topology
@INPUT      : polygons
@OUTPUT     : 
@RETURNS    : TRUE or FALSE
@DESCRIPTION: Returns true if the polygons are of sphere topology.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  BOOLEAN  is_this_sphere_topology(
    polygons_struct  *polygons )
{
    int     tess;

    return( get_tessellation_of_polygons_sphere( polygons, &tess ) );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : get_tessellation_of_polygons_sphere
@INPUT      : polygons
@OUTPUT     : tess
@RETURNS    : TRUE or FALSE
@DESCRIPTION: Determines if the polygons is of sphere topology and if so,
              what tessellation.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  BOOLEAN  get_tessellation_of_polygons_sphere(
    polygons_struct  *polygons,
    int              *tess )
{
    Real    size; 
    int     int_size, n_around_top, n_triangles, i;
    BOOLEAN is_sphere;

    is_sphere = FALSE;

    if( polygons->n_items > 10 )
    {
#ifdef  N_AROUND_TOP
        size = 1 + sqrt( (double) polygons->n_items / 2.0 + 1.0 -
                        (double) N_AROUND_TOP );
        n_around_top = N_AROUND_TOP;
#else
        size = sqrt( (double) polygons->n_items / 2.0 );
        n_around_top = 2 * ROUND( size );
#endif

        if( IS_INT(size) )
        {
            *tess = (int) size;

            int_size = (int) size;
            while( (int_size & 1) == 0 )
            {
                int_size >>= 1;
            }

            is_sphere = (int_size == 1);
        }

        if( is_sphere )
        {
            n_triangles = 0;
            for_less( i, 0, polygons->n_items )
            {
                if( GET_OBJECT_SIZE( *polygons, i ) == 3 )
                    ++n_triangles;
            }

            if( n_triangles != 2 * n_around_top )
                is_sphere = FALSE;
        }
    }

    return( is_sphere );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : get_tessellation_with_n_points
@INPUT      : n_points
@OUTPUT     : 
@RETURNS    : tess
@DESCRIPTION: Given the number of points in the polygons, finds a sphere
              tessellation n_up that corresponds to this.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  int  get_tessellation_with_n_points(
    int   n_points )
{
    Real  a, b, c, s1, s2;
    int   n_sol;

#ifdef  N_AROUND_TOP
    a = 2.0;
    b = -6.0;
    c = 2.0 * (Real) N_AROUND_TOP + 2.0 - (Real) n_points;
#else
    a = 2.0;
    b = -2.0;
    c = 2.0 - (Real) n_points;
#endif

    n_sol = solve_quadratic( a, b, c, &s1, &s2 );

    if( n_sol == 1 || (n_sol == 2 && s2 <= 0.0) )
    {
        if( s1 > 0.0 )
            return( ROUND( s1 ) );
    }
    else if( n_sol == 2 )
        return( ROUND( s2 ) );

    return( 0 );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : half_sample_sphere_tessellation
@INPUT      : polygons
@OUTPUT     : half
@RETURNS    : 
@DESCRIPTION: Half samples the sphere tessellation polygon, maintaining a
              sphere topology.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :         1993    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

BICAPI  void  half_sample_sphere_tessellation(
    polygons_struct   *polygons,
    polygons_struct   *half )
{
    static Point   centre = { 0.0f, 0.0f, 0.0f };
    int            n_up, n_around, half_n_up, half_n_around;
    int            up, around, n_circum, point_index, half_point_index;

    if( get_tessellation_of_polygons_sphere( polygons, &n_up ) )
    {
        n_around = 2 * n_up;
        half_n_up = n_up / 2;
        half_n_around = 2 * half_n_up;
        create_polygons_sphere( &centre, 1.0, 1.0, 1.0, half_n_up,
                                half_n_around, FALSE, half );

        half->surfprop = polygons->surfprop;

        if( polygons->colour_flag == ONE_COLOUR ||
            polygons->colour_flag == PER_ITEM_COLOURS )
        {
            half->colour_flag = ONE_COLOUR;
            half->colours[0] = polygons->colours[0];
        }
        else
        {
            REALLOC( half->colours, half->n_points );
        }

        for_less( up, 0, half_n_up+1 )
        {
            if( up == 0 || up == half_n_up )
                n_circum = 1;
            else
                n_circum = half_n_around;

            for_less( around, 0, n_circum )
            {
                half_point_index = get_sphere_point_index( up, around,
                                                  half_n_up, half_n_around );
                point_index = get_sphere_point_index( 2 * up, 2 * around,
                                                  n_up, n_around );

                half->points[half_point_index] = polygons->points[point_index];

                if( half->colour_flag == PER_VERTEX_COLOURS )
                    half->colours[half_point_index] =
                                                polygons->colours[point_index];
            }
        }
    }
}
