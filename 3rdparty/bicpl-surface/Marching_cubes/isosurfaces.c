#include "bicpl_internal.h"
#include "bicpl/marching.h"

static  int  get_polygons(
    Marching_cubes_methods  method,
    int                     x,
    int                     y,
    int                     z,
    Real                    corners[2][2][2],
    Real                    isovalue,
    int                     *sizes[],
    voxel_point_type        *points[] )
{
    int           n_polygons;
    static   int  static_sizes[4] = { 3, 3, 3, 3 };

    switch( method )
    {
    case MARCHING_CUBES:
        n_polygons = compute_isotriangles_in_voxel( corners, isovalue, points );
        *sizes = static_sizes;
        break;

    case MARCHING_NO_HOLES:
        n_polygons = get_holeless_isosurface_polygons( corners, isovalue, sizes,
                                                    points );
        break;

    case MARCHING_TETRA:
        n_polygons = get_tetra_isosurface_polygons( x, y, z, corners,
                                                    isovalue, sizes, points );
        break;
    }

    return( n_polygons );
}

static  BOOLEAN  is_binary_inside(
    Real   value,
    Real   min_value,
    Real   max_value )
{
    return( min_value <= value && value <= max_value );
}

BICAPI  int  compute_isosurface_in_voxel(
    Marching_cubes_methods  method,
    int                     x,
    int                     y,
    int                     z,
    Real                    corners[2][2][2],
    BOOLEAN                 binary_flag,
    Real                    min_value,
    Real                    max_value,
    int                     *sizes[],
    voxel_point_type        *points[] )
{
    int   i, j, k;
    Real  binary_corners[2][2][2];

    if( binary_flag )
    {
        for_less( i, 0, 2 )
            for_less( j, 0, 2 )
                for_less( k, 0, 2 )
                {
                    if( is_binary_inside( corners[i][j][k],
                                          min_value, max_value ) )
                    {
                        binary_corners[i][j][k] = 1.0;
                    }
                    else
                        binary_corners[i][j][k] = 0.0;
                }

        return( get_polygons( method, x, y, z, binary_corners,
                              0.5, sizes, points ) );
    }
    else
        return( get_polygons( method, x, y, z, corners, min_value,
                              sizes, points ) );
}

BICAPI  Point_classes  get_isosurface_point(
    Real              corners[2][2][2],
    int               voxel[],
    int               edge_intersected,
    BOOLEAN           binary_flag,
    Real              min_value,
    Real              max_value,
    Real              point[] )
{
    int             dim;
    int             v1[N_DIMENSIONS], v2[N_DIMENSIONS];
    int             offset[N_DIMENSIONS];
    Real            alpha, val1, val2;
    Point_classes   point_class;

    v1[0] = voxel[0];
    v1[1] = voxel[1];
    v1[2] = voxel[2];

    translate_from_edge_index( edge_intersected, offset );

    for_less( dim, 0, N_DIMENSIONS )
        v2[dim] = v1[dim] + offset[dim];

    val1 = corners[v1[0]][v1[1]][v1[2]];
    val2 = corners[v2[0]][v2[1]][v2[2]];

    if( binary_flag )
    {
        if( is_binary_inside( val1, min_value, max_value ) !=
            is_binary_inside( val2, min_value, max_value ) )
        {
            point_class = ON_EDGE;
            for_less( dim, 0, N_DIMENSIONS )
                point[dim] = (Real) (v1[dim] + v2[dim]) / 2.0;
        }
        else
            point_class = (Point_classes) -1;
    }
    else
    {
        if( val1 == min_value )
        {
            point_class = ON_FIRST_CORNER;
            alpha = 0.0;
        }
        else if( val2 == min_value )
        {
            point_class = ON_SECOND_CORNER;
            alpha = 1.0;
        }
        else if( (val1 < min_value && val2 > min_value) ||
                 (val1 > min_value && val2 < min_value) )
        {
            point_class = ON_EDGE;
            alpha = (min_value - val1) / (val2 - val1);
        }
        else
            point_class = (Point_classes) -1;

        if( point_class >= (Point_classes) 0 )
        {
            for_less( dim, 0, N_DIMENSIONS )
            {
                point[dim] = INTERPOLATE( alpha,
                                          (Real) v1[dim], (Real) v2[dim] );
            }
        }
    }

    return( point_class );
}

BICAPI  int  get_max_marching_edges(
    Marching_cubes_methods  method )
{
    switch( method )
    {
    case MARCHING_CUBES:
    case MARCHING_NO_HOLES:
        return( 3 );

    case MARCHING_TETRA:
        return( N_MARCHING_TETRA_EDGES );

    default:
        print_error( "Invalid parameter in get_max_marching_edges.\n" );
        return( 0 );
    }
}

BICAPI  int  get_max_marching_polygons_per_voxel(
    Marching_cubes_methods  method )
{
    switch( method )
    {
    case MARCHING_CUBES:
        return( 4 );

    case MARCHING_NO_HOLES:
        return( 4 );

    case MARCHING_TETRA:
        return( 5 );

    default:
        print_error(
               "Invalid parameter in get_max_marching_polygons_per_voxel.\n" );
        return( 0 );
    }
}

