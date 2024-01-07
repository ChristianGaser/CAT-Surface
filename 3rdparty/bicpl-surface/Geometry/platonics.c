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
static char rcsid[] = "$Header: /private-cvsroot/libraries/bicpl/Geometry/platonics.c,v 1.5 2005/08/17 22:30:25 bert Exp $";
#endif

static  void  fill_in_indices(
    polygons_struct  *polygons,
    int              n_indices,
    int              indices[] )
{
    int   ind;

    if( polygons->end_indices[polygons->n_items-1] != n_indices ) 
    {
        handle_internal_error( "create_unit_(platonic)" );
        return;
    }

    for_less( ind, 0, n_indices )
        polygons->indices[ind] = indices[ind];

    compute_polygon_normals( polygons );
}

BICAPI  void  create_unit_tetrahedron(
    polygons_struct  *polygons )
{
    Real         x, z, c;
    static  int  indices[] = {
                                 0, 1, 2,
                                 0, 2, 3,
                                 0, 3, 1,
                                 1, 3, 2
                             };

    initialize_polygons_with_size( polygons, WHITE, NULL, 4, 4, 3 );

    x = sqrt( 8.0 ) / 3.0;
    z = - 1.0 / 3.0;
    c = sqrt( 6.0 ) / 3.0;

    fill_Point( polygons->points[0],      0.0, 0.0, 1.0 );
    fill_Point( polygons->points[1],        x, 0.0,   z );
    fill_Point( polygons->points[2], -0.5 * x,   c,   z );
    fill_Point( polygons->points[3], -0.5 * x,  -c,   z );

    fill_in_indices( polygons, SIZEOF_STATIC_ARRAY(indices), indices );
}

BICAPI  void  create_unit_cube(
    polygons_struct  *polygons )
{
    Real         d;
    static  int  indices[] = {
                                 0, 1, 3, 2,
                                 4, 6, 7, 5,
                                 0, 4, 5, 1,
                                 2, 3, 7, 6,
                                 0, 2, 6, 4,
                                 1, 5, 7, 3
                             };

    initialize_polygons_with_size( polygons, WHITE, NULL, 8, 6, 4 );

    d = 1.0 / sqrt( 3.0 );

    fill_Point( polygons->points[0], -d, -d, -d );
    fill_Point( polygons->points[1], -d, -d,  d );
    fill_Point( polygons->points[2], -d,  d, -d );
    fill_Point( polygons->points[3], -d,  d,  d );
    fill_Point( polygons->points[4],  d, -d, -d );
    fill_Point( polygons->points[5],  d, -d,  d );
    fill_Point( polygons->points[6],  d,  d, -d );
    fill_Point( polygons->points[7],  d,  d,  d );

    fill_in_indices( polygons, SIZEOF_STATIC_ARRAY(indices), indices );
}

BICAPI  void  create_unit_octohedron(
    polygons_struct  *polygons )
{
    static  int  indices[] = {
                                 0, 3, 4,
                                 0, 4, 2,
                                 0, 2, 5,
                                 0, 5, 3,
                                 1, 3, 5,
                                 1, 5, 2,
                                 1, 2, 4,
                                 1, 4, 3
                             };

    initialize_polygons_with_size( polygons, WHITE, NULL, 6, 8, 3 );

    fill_Point( polygons->points[0], -1.0,  0.0,  0.0 );
    fill_Point( polygons->points[1],  1.0,  0.0,  0.0 );
    fill_Point( polygons->points[2],  0.0, -1.0,  0.0 );
    fill_Point( polygons->points[3],  0.0,  1.0,  0.0 );
    fill_Point( polygons->points[4],  0.0,  0.0, -1.0 );
    fill_Point( polygons->points[5],  0.0,  0.0,  1.0 );

    fill_in_indices( polygons, SIZEOF_STATIC_ARRAY(indices), indices );
}

BICAPI  void  create_unit_icosahedron(
    polygons_struct  *polygons )
{
    Real         a, b;
    static  int  indices[] = {
                                2, 11, 3,
                                2, 5, 11,
                                11, 5, 10,
                                5, 0, 10,
                                0, 5, 4,
                                4, 5, 2,
                                4, 2, 9,
                                9, 2, 3,
                                9, 3, 6,
                                9, 6, 8,
                                4, 9, 8,
                                4, 8, 0,
                                0, 8, 1,
                                1, 8, 6,
                                1, 6, 7,
                                7, 6, 3,
                                11, 7, 3,
                                11, 10, 7,
                                10, 1, 7,
                                10, 0, 1
                             };

    initialize_polygons_with_size( polygons, WHITE, NULL, 12, 20, 3 );

    a = sqrt( (5.0+sqrt(5.0))/10.0 );
    b = sqrt( (5.0-sqrt(5.0))/10.0 );

    fill_Point( polygons->points[0],  a,  b, 0.0 );
    fill_Point( polygons->points[1],  a, -b, 0.0 );
    fill_Point( polygons->points[2], -a,  b, 0.0 );
    fill_Point( polygons->points[3], -a, -b, 0.0 );
    fill_Point( polygons->points[4], 0.0,  a,  b );
    fill_Point( polygons->points[5], 0.0,  a, -b );
    fill_Point( polygons->points[6], 0.0, -a,  b );
    fill_Point( polygons->points[7], 0.0, -a, -b );
    fill_Point( polygons->points[8],  b, 0.0,  a );
    fill_Point( polygons->points[9], -b, 0.0,  a );
    fill_Point( polygons->points[10],  b, 0.0, -a );
    fill_Point( polygons->points[11], -b, 0.0, -a );

    fill_in_indices( polygons, SIZEOF_STATIC_ARRAY(indices), indices );
}
